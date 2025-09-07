import pandas as pd
import os
import numpy as np
import re
import cobra
import io


class model:
    """ a COMETS metabolic model to use in a layout
    
        A model contains information including the stoichiometric matrix,
        reaction boundaries and Michaelis-Menten parameters, founder populations,
        and various attributes related to motion in space.
        
        A model is usually initiated by supplying a cobrapy model, and to be used
        in a COMETS simulation must minimally be given an initial population:
            
        Parameters
        ----------
        
        model : cobra.Model or str, optional
            Either a cobrapy Model, a path to a .cmd COMETS model, a path to a
            cobra sbml model which could be loaded with cobra.io.read_sbml_model()
            or None.
            
        
        Attributes
        ----------
            
        initial_pop : list
            list of lists of spatial founder pops. e.g.[[x, y, grams], [x,y,grams]]
        id : str
            the name of the model. should be unique. 
        reactions : pandas.DataFrame
            DataFrame containing reaction information including boundaries. 
        smat : pandas.DataFrame
            DataFrame containing the stoichiometric matrix
        metabolites : pandas.DataFrame
            single-column DataFrame containing the names of the metabolites
        signals : pandas.Dataframe
            DataFrame on reactions whose bounds depend on metabolite concentrations
        default_vmax : float
            default uptake vmax in mmol / gDW /hr for Michaelis-Menten kinetics
        default_km : float
            default half-max metabolite concentration (M) for M-M kinetics
        default_hill : float
            default hill coefficient for hill-like uptake kinetics
        default_bounds : list
            list of two floats that are the default lower and upper reaction bounds
        obj_style : str
            one of MAXIMIZE_OBJECTIVE_FLUX (fba) or MAX_OBJECTIVE_MIN_TOTAL (pfba)
        optimizer : str
            one of "GLOP", "GUROBI" or "GLPK". not all functionality works with GLPK
            
        Examples
        --------
        
            >>> import cobra.test
            >>> import cometspy as c
            >>> ecoli = cobra.test.create_test_model("ecoli")
            >>> model = c.model(ecoli)
            >>> model.initial_pop = [0, 0, 1.e-12] # puts 1.e-12 gDW of biomass at 0,0
            >>> print(model.id)
            iJO1366

    """
    def __init__(self, model : cobra.Model = None, randomtag : bool = False):
        self.initial_pop = [[0, 0, 0.0]]

        if randomtag:
            self.id = '_' + hex(id(self))
        else:
            self.id = ""

        self.reactions = pd.DataFrame(columns=['REACTION_NAMES', 'ID',
                                               'LB', 'UB', 'EXCH',
                                               'EXCH_IND', 'V_MAX',
                                               'KM', 'HILL'])
        self.smat = pd.DataFrame(columns=['metabolite',
                                          'rxn',
                                          's_coef'])
        self.metabolites = pd.DataFrame(columns=['METABOLITE_NAMES'])
        self.signals = pd.DataFrame(columns=['REACTION_NUMBER',
                                             'EXCH_IND',
                                             'BOUND',
                                             'FUNCTION',
                                             'PARAMETERS',
                                             'REACTION_NAMES', 'EXCH'],
                                    dtype=object)
        # multitoxins is a special case of signaling
        self.multitoxins = pd.DataFrame(columns = ['REACTION_NUMBER',
                                                   'EXCH_INDS',
                                                   'BOUND',
                                                   'KMS',
                                                   'HILLS',
                                                   'VMAX',
                                                   'REACTION_NAME',
                                                   'EXCH_NAMES'],
                                        dtype = object)
        self.light = []

        self.vmax_flag = False
        self.km_flag = False
        self.hill_flag = False
        self.convection_flag = False
        self.light_flag = False

        self.nonlinear_diffusion_flag = False
        self.neutral_drift_flag = False
        self.noise_variance_flag = False
        self.default_vmax = 10
        self.default_km = 1
        self.default_hill = 1
        self.default_bounds = [0, 1000]
        self.objectives = {} #Empty dictionary
        self.objective = None
        self.optimizer = 'GUROBI'
        self.obj_style = 'MAXIMIZE_OBJECTIVE_FLUX'

        self.multispecies_convmodel_flag = False



        if model is not None:
            if isinstance(model, cobra.Model):
                self.load_cobra_model(model, randomtag)
            else:  # assume it is a path
                if model[-3:] == "cmd":
                    self.read_comets_model(model, randomtag)
                else:
                    self.read_cobra_model(model, randomtag)

    def get_reaction_names(self) -> list:
        """ returns a list of reaction names"""
        return(list(self.reactions['REACTION_NAMES']))
    
    def add_multitoxin(self, rxn_num, exch_ids, bound,
                       vmax, kms, hills):
        """
        Adds a signaling relationship where >= 1 signal close a bound 
        
        This incorporates a hill-function reduction of (usually) an upper bound
        It is useful for incorporating multiple, additively-functioning signals.
        
        bound = vmax * MULT((km_i^h_i) / (km_i^h_i + M_i^h_i)) 
        where MULT is a multiplicative function, M_i is the molarity of the 
        metabolite in the given exch_id, and km_i and h_i are the km and h for
        for that signal. 

        Parameters
        ----------
        rxn_num : int
            number of the affected reaction
        exch_ids : list[int]
            exchange numbers of the signals
        bound : str
            "ub" or "lb"
        vmax : float
            maximum bound in absence of signals
        kms : list[float]
            kms for the signals
        hills : list[float]
            hill coefficients for the signals

        Returns
        -------
        None.

        """
        rxn_name = self.reactions.loc[self.reactions.ID == rxn_num,
                              'REACTION_NAMES'].item()
        rxn_num = str(rxn_num)
        exch_names = [list(self.get_exchange_metabolites())[exch_id-1] for
                      exch_id in exch_ids]
        new_row = pd.DataFrame({'REACTION_NUMBER': rxn_num,
                                'EXCH_INDS': 1,
                                'BOUND': bound,
                                'KMS': 1,
                                'HILLS': 1,
                                'VMAX' :vmax,
                                'REACTION_NAME': rxn_name,
                                'EXCH_NAMES': ""},
                               index=[0],
                               dtype=object)
        new_row.loc[0, 'KMS'] = kms
        new_row.loc[0, 'HILLS'] = hills
        new_row.loc[0, "EXCH_INDS"] = exch_ids
        new_row.loc[0, 'EXCH_NAMES'] = exch_names
        self.multitoxins = self.multitoxins.append(new_row, ignore_index=True)

    def add_signal(self, rxn_num, exch_id, bound,
                   function, parms):
        """adds a signal to the reaction rxn_num that changes bounds
        
            Parameters
            ----------
            rxn_num : int or 'death'
                reaction affected (int) or whether the signal causes death 'death'
            exch_id : int
                the number of the exchange reaction. see model.metabolites
            bound : str  'ub' or 'lb'
                specifies whether the upper or lower bound is affected
            function : str 'linear' or 'bounded_linear' or 'generalized_logistic'
                specifies the function relating the metabolite conc. to the bound
            parms : list (float)
                a list of floats for the function 
                
        """

        if str(rxn_num).lower().strip() == 'death':
            rxn_name = 'death'
            rxn_num = 'death'
        else:
            rxn_name = self.reactions.loc[self.reactions.ID == rxn_num,
                                          'REACTION_NAMES'].item()
            rxn_num = str(rxn_num)
	# note: rxn_num matches the pandas DataFrame in this object.
	#       however, COMETS actually wants that rxn_num - 1, which
	#       is why it does that during the saving process, and why
	#	you will see it that way in saved model files.
        exch_name = list(self.get_exchange_metabolites())[exch_id-1]
        new_row = pd.DataFrame({'REACTION_NUMBER': rxn_num,
                                'EXCH_IND': exch_id,
                                'BOUND': bound,
                                'FUNCTION': function,
                                'PARAMETERS': 1,
                                'REACTION_NAMES': rxn_name,
                                'EXCH': exch_name},
                               index=[0],
                               dtype=object)
        new_row.loc[0, 'PARAMETERS'] = parms
        self.signals = self.signals.append(new_row, ignore_index=True)

    def add_neutral_drift_parameter(self, neutralDriftSigma):
        """ sets the neutralDriftSigma parameter """
        if not isinstance(neutralDriftSigma, float):
            raise ValueError("neutralDriftSigma must be a float")
        #self.neutral_drift_flag = True
        self.neutralDriftSigma = neutralDriftSigma

    def add_neutral_drift(self, neutralDrift):
        """ sets the neutralDriftSigma parameter """
        if not isinstance(neutralDrift, bool):
            raise ValueError("neutralDriftSigma must be a boolean")
        self.neutral_drift_flag = neutralDrift

    def add_nonlinear_diffusion_parameters(self,
                                           D0 : float=1.,
                                           Dk : float =1.,
                                           exponent : float=1.,
                                           hilln : float=10.,
                                           hillk : float=0.9):
        """ sets the model to use non-linear diffusion biomass spread.
            
            This also requires one set the biomassMotionStyle to 
            'ConvNonlin Diffusion 2D' in the comets.params object (see Example)
            
            Parameters
            ----------
            D0 : float, optional
            Dk : float, optional
            exponent : float, optional
            hilln : float, optional
            hillk : float, optional
            
            Example
            --------
            
            >>> import cobra.test
            >>> import cometspy as c
            >>> model = c.model(cobra.test.create_test_model("ecoli"))
            >>> model.add_nonlinear_diffusion_parameters(1., 1., 2., 5., 0.5)
            >>> params = c.params()
            >>> params.set_param("biomassMotionStyle", "ConvNonlin Diffusion 2D")
            
        """
        
        for parm in [D0, Dk,
                     exponent, hilln,
                     hillk]:
            if not isinstance(parm, float):
                raise ValueError('all nonlinear diffusion terms must be float')
        self.nonlinear_diffusion_flag = True
        self.nonlinear_diffusion_parameters = {'convNonLinDiffZero': D0,
                                               'convNonlinDiffN': Dk,
                                               'convNonlinDiffExponent': exponent,
                                               'convNonlinDiffHillN': hilln,
                                               'convNonlinDiffHillK': hillk}
    #H. Shi Jan 2023
    def add_nonlinear_diffusion_chemotaxis_parameters(self,
                                           hilln : float=1.0,
                                           hillk : float=0.0):
        """ sets the model to use non-linear diffusion chemotaxis biomass spread.
            
            This also requires one set the biomassMotionStyle to 
            'Nonlin Diff Chemotaxis 2D' in the comets.params object (see Example)
            
            Parameters
            ----------
            hilln : float, optional
            hillk : float, optional
            
            Example
            --------
            
            >>> import cobra.test
            >>> import cometspy as c
            >>> model = c.model(cobra.test.create_test_model("ecoli"))
            >>> model.add_nonlinear_diffusion_parameters(5., 0.5)
            >>> params = c.params()
            >>> params.set_param("biomassMotionStyle", "Nonlin Diff Chemotaxis 2D")
            
        """
        
        for parm in [hilln,
                     hillk]:
            if not isinstance(parm, float):
                raise ValueError('all nonlinear diffusion chemotaxis terms must be float')
        self.nonlinear_diffusion_ctx_flag = True
        self.nonlinear_diffusion_ctx_parameters = {
                                               'NonlinDiffCtxHillN': hilln,
                                               'NonlinDiffCtxHillK': hillk}
    def add_light(self, reaction : str, 
                  abs_coefficient : float, 
                  abs_base : float):
        """Causes a reaction to function in response to light
            
            Parameters
            ----------
            
            reaction : str
                the name of the reaction affected
            abs_coefficient : float
                the absorption coefficient
            abs_base : float
                absorption baseline
                
        """
        if (reaction not in self.reactions['REACTION_NAMES'].values):
            raise ValueError('the reaction is not present in the model')
        self.light.append([reaction, abs_coefficient, abs_base])
        self.light_flag = True

    def add_convection_parameters(self, packedDensity : float =1.,
                                  elasticModulus : float =1.,
                                  frictionConstant : float =1.,
                                  convDiffConstant : float =1.):
        """ sets the parameters for biomass spread to use convection
        
            This also requires one set the biomassMotionStyle to 
            'Convection 2D' in the comets.params object (see Example)
            
            Parameters
            ----------
            
            packedDensity : float, optional
            elasticModulus : float, optional
            frictionConstant : float, optional
            convDiffConstant : float, optional
            
            Example
            -------
            
            import cobra.test
            import cometspy as c
            model = c.model(cobra.test.create_test_model("ecoli"))
            model.add_convection_parameters(1., 0.5, 2., 10.e-5)
            params = c.params()
            params.set_param("biomassMotionStyle", "Convection 2D") 
            
        """ 
        if not isinstance(packedDensity, float):
            raise ValueError('packed_density must be a float')
        if not isinstance(elasticModulus, float):
            raise ValueError('elasticModulus must be a float')
        if not isinstance(frictionConstant, float):
            raise ValueError('frictionConstant must be a float')
        if not isinstance(convDiffConstant, float):
            raise ValueError('convDiffConstant must be a float')
        self.convection_flag = True
        self.convection_parameters = {'packedDensity': packedDensity,
                                      'elasticModulus': elasticModulus,
                                      'frictionConstant': frictionConstant,
                                      'convDiffConstant': convDiffConstant}
        
    def add_multspecies_convmodel_parameters(self, pressureKappa : float, pressureExponent : float, packBiomass : float, maxPressure : float):
        if not isinstance(pressureKappa, float):
            raise ValueError('pressureKappa must be a float')
        if not isinstance(pressureExponent, float):
            raise ValueError('pressureExponent must be a float')
        if not isinstance(packBiomass, float):
            raise ValueError('packBiomass must be a float')
        if not isinstance(maxPressure, float):
            raise ValueError('maxPressure must be a float')
        self.multispecies_convmodel_flag = True
        self.multimodel_parameters = {'pressureKappa': pressureKappa,
                                      'pressureExponent': pressureExponent,
                                      'packBiomass': packBiomass,
                                      'maxPressure': maxPressure}

    def add_noise_variance_parameter(self, noiseVariance : float):
        """ sets the noise variance parameter 
        
            Parameters
            ----------
            noiseVariance : float
            
        """
        if not isinstance(noiseVariance, float):
            raise ValueError('noiseVariance must be a float')
        self.noise_variance_flag = True
        self.noise_variance = noiseVariance
        
    def ensure_sinks_are_not_exchanges(self, suffix : str = '_c'):
        """ set exchange reactions ending in suffix (default = '_c') to sink
        
            Parameters
            ----------
            suffix : str, optional
                the suffix to look for in exchange reactions. default = '_c'

        """
        rxn_names = [rxn for rxn in self.reactions.REACTION_NAMES.to_list() if rxn[-2:] == suffix]
        for rxn in rxn_names:
            self.reactions.loc[self.reactions.REACTION_NAMES == rxn, 'EXCH'] = False
            self.reactions.loc[self.reactions.REACTION_NAMES == rxn, 'EXCH_IND'] = 0

    def open_exchanges(self, lower_bound : float = -1000., upper_bound : float = 1000.):
        """ sets all exchange boundaries to the specifed values 
            
            This method is useful when first making a COMETS model so that medium
            definitions are not carried over from cobra into COMETS. 
            
            Parameters
            ----------
            
            lower_bound : float, optional
                default = -1000.   in units of mmol / gDW / hr
            upper_bound : float, optional
                default = 1000.    in units of mmol / gDW / hr
                
        """
        
        self.reactions.loc[self.reactions.EXCH,'LB'] = lower_bound
        self.reactions.loc[self.reactions.EXCH,'UB'] = upper_bound

    def get_exchange_metabolites(self) -> list:
        """ returns a list of the names of the exchange metabolites """
        exchmets = pd.merge(self.reactions.loc[self.reactions['EXCH'], 'ID'],
                            self.smat,
                            left_on='ID', right_on='rxn',
                            how='inner')['metabolite']
        exchmets = self.metabolites.iloc[exchmets-1]
        return(exchmets.METABOLITE_NAMES)

    def change_bounds(self, reaction : str, lower_bound : float, upper_bound : float):
        """ changes the bounds of the specified reaction
        
            Parameters
            ----------
            
            reaction : str
                the name of the reaction
            lower_bound : float
                the new value of the reaction's lower bound in mmol / gDW / hr
            upper_bound : float
                the new value of the reaction's upper bound in mmol / gDW / hr
                
        """
        if reaction not in self.reactions['REACTION_NAMES'].values:
            print('reaction couldnt be found')
            return
        self.reactions.loc[self.reactions['REACTION_NAMES'] == reaction,
                           'LB'] = lower_bound
        self.reactions.loc[self.reactions['REACTION_NAMES'] == reaction,
                           'UB'] = upper_bound

    def get_bounds(self, reaction : str) -> tuple:
        """returns a tuple (lb, ub) for the specified reaction
        
            Parameters
            ----------
            
            reaction : str
                the name of the reaction
        
        """
        if reaction not in self.reactions['REACTION_NAMES'].values:
            print('reaction couldnt be found')
            return
        lb = float(self.reactions.loc[self.reactions[
            'REACTION_NAMES'] == reaction, 'LB'])
        ub = float(self.reactions.loc[self.reactions[
            'REACTION_NAMES'] == reaction, 'UB'])
        return((lb, ub))

    def change_vmax(self, reaction : str, vmax : float):
        """ changes the vmax of the specified reaction. 
        
            Parameters
            ----------
            
            reaction : str
                the name of the reaction
            vmax : float
                the new vmax value for the reaction, for Monod (Michaelis-Menten)
                kinetics
        
        """
        if reaction not in self.reactions['REACTION_NAMES'].values:
            print('reaction couldnt be found')
            return
        self.vmax_flag = True
        self.reactions.loc[self.reactions[
            'REACTION_NAMES'] == reaction, 'V_MAX'] = vmax

    def change_km(self, reaction, km):
        """ changes the km of the specified reaction. 
        
            Parameters
            ----------
            
            reaction : str
                the name of the reaction
            km : float
                the new km value for the reaction, for Monod (Michaelis-Menten)
                kinetics
        
        """
        if reaction not in self.reactions['REACTION_NAMES'].values:
            print('reaction couldnt be found')
            return
        self.km_flag = True
        self.reactions.loc[self.reactions[
            'REACTION_NAMES'] == reaction, 'KM'] = km

    def change_hill(self, reaction, hill):
        """ changes the hill coefficient of the specified reaction. 
        
            Parameters
            ----------
            
            reaction : str
                the name of the reaction
            hill : float
                the new hill coef. for the reaction, for Hill-like kinetics
        
        """
        if reaction not in self.reactions['REACTION_NAMES'].values:
            print('reaction couldnt be found')
            return
        self.hill_flag = True
        self.reactions.loc[self.reactions[
            'REACTION_NAMES'] == reaction, 'HILL'] = hill

    def change_optimizer(self, optimizer):
        """ changes the optimizer to a specified one. 
        
            Parameters
            ----------
            
            optimizer : str
                the name of the optimizer
        
        """
        self.optimizer = optimizer

    def change_objective(self, reaction, weight):
        """ changes the list of objective reactions

            Parameters
            ----------

            reaction : str
                the name of the reaction
            weight : float
                the weight for the reaction (0 to remove it as an objective)

        """

        rxnIdx = self.reactions.loc[self.reactions['REACTION_NAMES'] == reaction]['ID'].iloc[0]

        if (weight == 0):
            if (rxnIdx in self.objectives):
                del self.objectives[rxnIdx]
            return

        self.objectives[rxnIdx] = weight

    def change_biomass(self, reaction):
        """ changes the biomass that's logged

            Parameters
            ----------

            reaction : str
                the name of the reaction
        """

        self.biomass = self.reactions.loc[self.reactions['REACTION_NAMES'] == reaction]['ID'].iloc[0]

    def change_objective_style(self, style):
        """ changes the objective style

            Parameters
            ----------

            reaction : str
                the new objective style
        """

        self.obj_style = style

    def change_maintenance(self, reaction, flux):
        """ changes the maintenance reaction and its flux

            Parameters
            ----------

            reaction : str
                the name of the reaction
            flux : float
                the minimum flux for the cell to survive

        """

        rxnIdx = self.reactions.loc[self.reactions['REACTION_NAMES'] == reaction]['ID'].iloc[0]

        self.maintenanceRxn = rxnIdx
        self.maintenanceFlux = flux

    def read_cobra_model(self, path : str, randomtag : bool = False):
        """ reads a cobra model from a file and loads it into this model 
            
            This is an alternative way to initialize a COMETS model, by supplying
            it with the path to a cobra model that can be read using
            cobra.io.read_sbml_model.
            
            Parameters
            ----------
            
            path : str
                path to a cobra model in sbml format
            
            Example
            -------
            
            >>> import cometspy as c
            >>> path_to_file = "./my_new_model.xml"
            >>> model = c.model()
            >>> model.read_cobra_model(path_to_file)
        
        """
        curr_m = cobra.io.read_sbml_model(path)
        self.load_cobra_model(curr_m, randomtag)

    def load_cobra_model(self, curr_m : cobra.Model, randomtag : bool = False):
        """ creates the COMETS model from the supplied cobra model
        
            This is usually used internally when creating a COMETS model. 
            
            Parameters
            ----------
            
            curr_m : cobra.Model
                the cobra model object to be converted to the COMETS model object
                
            Example
            -------
            
            >>> import cometspy as c
            >>> import cobra.test
            >>> ecoli = cobra.test.create_test_model('ecoli')
            >>> model = c.model()
            >>> model.load_cobra_model(ecoli)
            
        """ 
        self.id = curr_m.id
        if randomtag:
            self.id = self.id + '_' + hex(id(self))

        # reactions and their features
        reaction_list = curr_m.reactions
        self.reactions['REACTION_NAMES'] = [str(x).split(':')[0] for
                                            x in reaction_list]
        self.reactions['ID'] = [k for k in
                                range(1, len(reaction_list)+1)]
        self.reactions['LB'] = [x.lower_bound for x in reaction_list]
        self.reactions['UB'] = [x.upper_bound for x in reaction_list]

        self.reactions['EXCH'] = [True if (len(k.metabolites) == 1) &
                                  (list(k.metabolites.
                                        values())[0] == (-1)) &
                                  ('DM_' not in k.id)
                                  else False for k in reaction_list]

        exch = self.reactions.loc[self.reactions['EXCH'], 'ID'].tolist()
        self.reactions['EXCH_IND'] = [exch.index(x)+1
                                      if x in exch else 0
                                      for x in self.reactions['ID']]

        self.reactions['V_MAX'] = [k.Vmax
                                   if hasattr(k, 'Vmax')
                                   else float('NaN')
                                   for k in reaction_list]

        if not self.reactions.V_MAX.isnull().all():
            self.vmax_flag = True

        self.reactions['KM'] = [k.Km
                                if hasattr(k, 'Km')
                                else float('NaN')
                                for k in reaction_list]

        if not self.reactions.KM.isnull().all():
            self.km_flag = True

        self.reactions['HILL'] = [k.Hill
                                  if hasattr(k, 'Hill')
                                  else float('NaN')
                                  for k in reaction_list]

        if not self.reactions.HILL.isnull().all():
            self.hill_flag = True

        if self.vmax_flag:
            if hasattr(curr_m, 'default_vmax'):
                self.default_vmax = curr_m.default_vmax

        if self.km_flag:
            if hasattr(curr_m, 'default_km'):
                self.default_km = curr_m.default_km

        if self.hill_flag:
            if hasattr(curr_m, 'default_hill'):
                self.default_hill = curr_m.default_hill

        # Metabolites
        metabolite_list = curr_m.metabolites
        self.metabolites['METABOLITE_NAMES'] = [str(x) for
                                                x in metabolite_list]

        # S matrix
        for index, row in self.reactions.iterrows():
            rxn = curr_m.reactions.get_by_id(
                row['REACTION_NAMES'])
            rxn_num = row['ID']
            rxn_mets = [1+list(self.metabolites[
                'METABOLITE_NAMES']).index(
                x.id) for x in rxn.metabolites]
            met_s_coefs = list(rxn.metabolites.values())

            cdf = pd.DataFrame({'metabolite': rxn_mets,
                                'rxn': [rxn_num]*len(rxn_mets),
                                's_coef': met_s_coefs})
            cdf = cdf.sort_values('metabolite')
            # skip concatenation on first row (prevents pandas FutureWarning)
            if index == 0:
                self.smat = cdf
            else:
                self.smat = pd.concat([self.smat, cdf])
                

        self.smat = self.smat.sort_values(by=['metabolite', 'rxn'])

        # The rest of stuff
        if hasattr(curr_m, 'default_bounds'):
            self.default_bounds = curr_m.default_bounds

        # set objective(s) using COBRA reactions objective coefficient information
        obj = {str(x).split(':')[0]:x.objective_coefficient 
               for x in reaction_list 
               if x.objective_coefficient != 0}
        obj = {rx: -1 if coef < 0 else 1 for rx, coef in obj.items()}

        self.objective = [int(self.reactions[self.reactions.REACTION_NAMES == rx]['ID'].iloc[0]) * coef for rx, coef in obj.items()]

        if hasattr(curr_m, 'comets_optimizer'):
            self.optimizer = curr_m.comets_optimizer

        if hasattr(curr_m, 'comets_obj_style'):
            self.obj_style = curr_m.comets_obj_style

    def read_comets_model(self, path : str, randomtag : bool = False):
        """ an alternative way to create a COMETS model by reading from a file
        
            This allows a user to load previously-saved COMETS models. The 
            contents populate this object's attributes.
            
            Parameters
            ----------
            
            path : str
                the path to the COMETS model file
                
            Example
            -------
            
            >>> import cometspy as c
            >>> path_to_model = "./iJO1366.cmd" # this must actually exist
            >>> model = c.model()
            >>> model.read_comets_model(path_to_model)
            
        """        
        self.id = os.path.splitext(os.path.basename(path))[0] + self.id
        if randomtag:
            self.id = self.id + '_' + hex(id(self))

        # in this way, its robust to empty lines:
        m_f_lines = [s for s in _read_file(path).splitlines() if s]
        m_filedata_string = os.linesep.join(m_f_lines)
        ends = []
        for k in range(0, len(m_f_lines)):
            if '//' in m_f_lines[k]:
                ends.append(k)

        # '''----------- S MATRIX ------------------------------'''
        lin_smat = re.split('SMATRIX',
                            m_filedata_string)[0].count('\n')
        lin_smat_end = next(x for x in ends if x > lin_smat)

        self.smat = pd.read_csv(io.StringIO('\n'.join(m_f_lines[
            lin_smat:lin_smat_end])),
                           delimiter=r'\s+',
                           skipinitialspace=True)
        self.smat.columns = ['metabolite', 'rxn', 's_coef']

        # '''----------- REACTIONS AND BOUNDS-------------------'''
        lin_rxns = re.split('REACTION_NAMES',
                            m_filedata_string)[0].count('\n')
        lin_rxns_end = next(x for x in
                            ends if x > lin_rxns)

        rxn = pd.read_csv(io.StringIO('\n'.join(m_f_lines[
            lin_rxns:lin_rxns_end])),
                          delimiter=r'\s+',
                          skipinitialspace=True)

        rxn['ID'] = range(1, len(rxn)+1)

        lin_bnds = re.split('BOUNDS',
                            m_filedata_string)[0].count('\n')
        lin_bnds_end = next(x for x in ends if x > lin_bnds)

        bnds = pd.read_csv(io.StringIO('\n'.join(m_f_lines[
            lin_bnds:lin_bnds_end])),
                           delimiter=r'\s+',
                           skipinitialspace=True)

        default_bounds = [float(bnds.columns[1]),
                          float(bnds.columns[2])]

        bnds.columns = ['ID', 'LB', 'UB']
        reactions = pd.merge(rxn, bnds,
                             left_on='ID', right_on='ID',
                             how='left')
        reactions.LB.fillna(default_bounds[0], inplace=True)
        reactions.UB.fillna(default_bounds[1], inplace=True)

        # '''----------- METABOLITES ---------------------------'''
        lin_mets = re.split('METABOLITE_NAMES',
                            m_filedata_string)[0].count('\n')
        lin_mets_end = next(x for x in ends if x > lin_mets)

        metabolites = pd.read_csv(io.StringIO('\n'.join(m_f_lines[
            lin_mets:lin_mets_end])),
                                  delimiter=r'\s+',
                                  skipinitialspace=True)

        # '''----------- EXCHANGE RXNS -------------------------'''
        lin_exch = re.split('EXCHANGE_REACTIONS',
                            m_filedata_string)[0].count('\n')+1
        exch = [int(k) for k in re.findall(r'\S+',
                                           m_f_lines[lin_exch].
                                           strip())]

        reactions['EXCH'] = [True if x in exch else False
                             for x in reactions['ID']]
        reactions['EXCH_IND'] = [exch.index(x)+1
                                 if x in exch else 0
                                 for x in reactions['ID']]

        # '''----------- VMAX VALUES --------------------------'''
        if 'VMAX_VALUES' in m_filedata_string:
            self.vmax_flag = True
            lin_vmax = re.split('VMAX_VALUES',
                                m_filedata_string)[0].count('\n')
            lin_vmax_end = next(x for x in ends if x > lin_vmax)

            Vmax = pd.read_csv(io.StringIO('\n'.join(m_f_lines[
                lin_vmax:lin_vmax_end])),
                               delimiter=r'\s+',
                               skipinitialspace=True)

            Vmax.columns = ['EXCH_IND', 'V_MAX']

            reactions = pd.merge(reactions, Vmax,
                                 left_on='EXCH_IND',
                                 right_on='EXCH_IND',
                                 how='left')
            self.default_vmax = float(m_f_lines[lin_vmax-1].split()[1])
        else:
            reactions['V_MAX'] = np.NaN

        # '''----------- VMAX VALUES --------------------------'''
        if 'KM_VALUES' in m_filedata_string:
            self.km_flag = True
            lin_km = re.split('KM_VALUES',
                              m_filedata_string)[0].count('\n')
            lin_km_end = next(x for x in ends if x > lin_km)

            Km = pd.read_csv(io.StringIO('\n'.join(m_f_lines[
                lin_km:lin_km_end])),
                             delimiter=r'\s+',
                             skipinitialspace=True)
            Km.columns = ['EXCH_IND', 'KM']

            reactions = pd.merge(reactions, Km,
                                 left_on='EXCH_IND',
                                 right_on='EXCH_IND',
                                 how='left')
            self.default_km = float(m_f_lines[lin_km-1].split()[1])
        else:
            reactions['KM'] = np.NaN

        # '''----------- VMAX VALUES --------------------------'''
        if 'HILL_COEFFICIENTS' in m_filedata_string:
            self.hill_flag = True
            lin_hill = re.split('HILL_COEFFICIENTS',
                                m_filedata_string)[0].count('\n')
            lin_hill_end = next(x for x in ends if x > lin_hill)

            Hill = pd.read_csv(io.StringIO('\n'.join(m_f_lines[
                lin_hill:lin_hill_end])),
                               delimiter=r'\s+',
                               skipinitialspace=True)
            Hill.columns = ['EXCH_IND', 'HILL']

            reactions = pd.merge(reactions, Hill,
                                 left_on='EXCH_IND',
                                 right_on='EXCH_IND',
                                 how='left')
            self.default_hill = float(m_f_lines[lin_hill-1].split()[1])
        else:
            reactions['HILL'] = np.NaN

        # '''----------- OBJECTIVE -----------------------------'''
        lin_obj = re.split('OBJECTIVE',
                           m_filedata_string)[0].count('\n')+1
        self.objective = int(m_f_lines[lin_obj].strip())

        # '''----------- OBJECTIVE STYLE -----------------------'''
        if 'OBJECTIVE_STYLE' in m_filedata_string:
            lin_obj_st = re.split('OBJECTIVE_STYLE',
                                  m_filedata_string)[0].count(
                                      '\n')+1
            self.obj_style = m_f_lines[lin_obj_st].strip()

        # '''----------- OPTIMIZER -----------------------------'''
        if 'OPTIMIZER' in m_filedata_string:
            lin_opt = re.split('OPTIMIZER',
                               m_filedata_string)[0].count('\n')
            self.optimizer = m_f_lines[lin_opt].split()[1]

        # '''--------------neutral drift------------------------'''
        if "neutralDrift" in m_filedata_string:
            lin_obj_st = re.split('neutralDrift',
                                  m_filedata_string)[0].count(
                                      '\n')
            if "TRUE" == m_f_lines[lin_obj_st].strip().split()[1].upper():
                self.neutral_drift_flag = True
                self.neutralDriftSigma = 0.
        if "neutralDriftsigma" in m_filedata_string:
            lin_opt = re.split('neutralDriftsigma',
                               m_filedata_string)[0].count('\n')
            self.neutralDriftSigma = float(m_f_lines[lin_opt].split()[1])

        # '''--------------convection---------------------------'''
        for parm in ['packedDensity', 'elasticModulus',
                     'frictionConstant', 'convDiffConstant']:
            if parm in m_filedata_string:
                lin_obj_st = re.split(parm,
                                      m_filedata_string)[0].count(
                                          '\n')
                parm_value = float(m_f_lines[lin_obj_st].strip().split()[1])
                try:
                    self.convection_parameters[parm] = parm_value
                except:  # TODO change bare except statements to ifelse
                    self.convection_flag = True
                    self.convection_parameters = {'packedDensity': 1.,
                                                  'elasticModulus': 1.,
                                                  'frictionConstant': 1.,
                                                  'convDiffConstant': 1.}
                    self.convection_parameters[parm] = parm_value

        # '''--------------non-linear diffusion---------------------------'''
        for parm in ['convNonLinDiffZero', 'convNonlinDiffN', 'convNonlinDiffExponent',
                     'convNonlinDiffHillN', 'convNonlinDiffHillK']:
            if parm in m_filedata_string:
                lin_obj_st = re.split(parm,
                                      m_filedata_string)[0].count(
                                          '\n')
                parm_value = float(m_f_lines[lin_obj_st].strip().split()[1])
                try:
                    self.nonlinear_diffusion_parameters[parm] = parm_value
                except:  # TODO change bare except statements to ifelse
                    self.nonlinear_diffusion_flag = True
                    self.nonlinear_diffusion_parameters = {'convNonLinDiffZero': 1.,
                                                           'convNonlinDiffN': 1.,
                                                           'convNonlinDiffExponent': 1.,
                                                           'convNonlinDiffHillN': 10.,
                                                           'convNonlinDiffHillK': .9}
                    self.nonlinear_diffusion_parameters[parm] = parm_value

        # '''--------------non-linear diffusion chemotaxis---------------------------'''
        #H. Shi Jan 2023
        for parm in ['NonlinDiffCtxHillN', 'NonlinDiffCtxHillK']:
            if parm in m_filedata_string:
                lin_obj_st = re.split(parm,
                                      m_filedata_string)[0].count(
                                          '\n')
                parm_value = float(m_f_lines[lin_obj_st].strip().split()[1])
                try:
                    self.nonlinear_diffusion_parameters[parm] = parm_value
                except:  # TODO change bare except statements to ifelse
                    self.nonlinear_diffusion_ctx_flag = True
                    self.nonlinear_diffusion_ctx_parameters = {'NonlinDiffCtxHillN': 1.0,
                                               'NonlinDiffCtxHillK': 0.0}
                    self.nonlinear_diffusion_ctx_parameters[parm] = parm_value  

        # '''-----------noise variance-----------------'''
        if 'noiseVariance' in m_filedata_string:
            lin_obj_st = re.split('noiseVariance',
                                  m_filedata_string)[0].count(
                                      '\n')
            noiseVariance = float(m_f_lines[lin_obj_st].strip().split()[1])

            self.noise_variance_flag = True
            self.noise_variance = noiseVariance
        # assign the dataframes we just built
        self.reactions = reactions
        self.metabolites = metabolites

    def delete_comets_model(self, working_dir = None):
        """ deletes a file version of this model if it exists.
        
            Parameters
            ----------
            
                working_dir : str, optional
                    the directory where the comets model file exists
                
        """
        path_to_delete = ""
        if working_dir is not None:
            path_to_delete = working_dir
        path_to_delete = path_to_delete + self.id + '.cmd'
        os.remove(path_to_delete)
        
    def write_comets_model(self, working_dir : str = None):
        """ writes the COMETS model object to a file 
            
            This writes the current model object in COMETS format to file, either
            in the current working directory or in the directory specified by
            the optional parameter working_dir.  It is mostly used by the comets
            object to run simulations, though a user may with to write these to 
            examine directly or to save. 
            
            Parameters
            ----------
            working_dir : str, optional
                a directory path to put COMETS files.  for example "./data_files/"
                
        """
        path_to_write = ""
        if working_dir is not None:
            path_to_write = working_dir
        path_to_write = path_to_write + self.id + '.cmd'

        # format variables for writing comets model
        bnd = self.reactions.loc[(self.reactions['LB']
                                  != self.default_bounds[0]) |
                                 (self.reactions['UB'] !=
                                  self.default_bounds[1]),
                                 ['ID', 'LB', 'UB']].astype(
                                     str).apply(lambda x: '   '.join(x),
                                                axis=1)
        bnd = '    ' + bnd.astype(str)

        rxn_n = '    ' + self.reactions['REACTION_NAMES'].astype(str)

        met_n = '    ' + self.metabolites.astype(str)

        smat = self.smat.astype(str).apply(lambda x:
                                           '   '.join(x), axis=1)
        smat = '    ' + smat.astype(str)

        exch_r = ' '.join([str(x) for x in
                           self.reactions.loc[self.reactions.EXCH, 'ID']])

        # optional fields (vmax,km, hill)
        if self.vmax_flag:
            Vmax = self.reactions.loc[self.reactions['V_MAX'].notnull(),
                                      ['EXCH_IND', 'V_MAX']]
            Vmax = Vmax.astype(str).apply(lambda x:
                                          '   '.join(x), axis=1)
            Vmax = '    ' + Vmax.astype(str)

        if self.km_flag:
            Km = self.reactions.loc[self.reactions['KM'].notnull(),
                                    ['EXCH_IND', 'KM']]
            Km = Km.astype(str).apply(lambda x:
                                      '   '.join(x), axis=1)
            Km = '    ' + Km.astype(str)

        if self.hill_flag:
            Hill = self.reactions.loc[self.reactions['HILL'].notnull(),
                                      ['EXCH_IND', 'HILL']]
            Hill = Hill.astype(str).apply(lambda x:
                                          '   '.join(x), axis=1)
            Hill = '    ' + Hill.astype(str)

        if os.path.isfile(path_to_write):
            os.remove(path_to_write)

        with open(path_to_write, 'a') as f:

            f.write('SMATRIX  ' + str(len(self.metabolites)) +
                    '  ' + str(len(self.reactions)) + '\n')
            smat.to_csv(f, mode='a', lineterminator = '\n', header=False, index=False)
            f.write(r'//' + '\n')

            f.write('BOUNDS ' +
                    str(self.default_bounds[0]) + ' ' +
                    str(self.default_bounds[1]) + '\n')
            bnd.to_csv(f, mode='a', lineterminator = '\n', header=False, index=False)
            f.write(r'//' + '\n')

            f.write('OBJECTIVE\n' +
                    '    ' + '    '.join([str(rx) for rx in self.objective]) + '\n')
            f.write(r'//' + '\n')

            f.write('METABOLITE_NAMES\n')
            met_n.to_csv(f, mode='a', lineterminator = '\n', header=False, index=False)
            f.write(r'//' + '\n')

            f.write('REACTION_NAMES\n')
            rxn_n.to_csv(f, mode='a', lineterminator = '\n', header=False, index=False)
            f.write(r'//' + '\n')

            f.write('EXCHANGE_REACTIONS\n')
            f.write(' ' + exch_r + '\n')
            f.write(r'//' + '\n')

            if self.vmax_flag:
                f.write('VMAX_VALUES ' +
                        str(self.default_vmax) + '\n')
                Vmax.to_csv(f, mode='a', lineterminator = '\n', header=False, index=False)
                f.write(r'//' + '\n')

            if self.km_flag:
                f.write('KM_VALUES ' +
                        str(self.default_km) + '\n')
                Km.to_csv(f, mode='a', lineterminator = '\n', header=False, index=False)
                f.write(r'//' + '\n')

            if self.hill_flag:
                f.write('HILL_VALUES ' +
                        str(self.default_hill) + '\n')
                Hill.to_csv(f, mode='a', lineterminator = '\n', header=False, index=False)
                f.write(r'//' + '\n')

            if self.light_flag:
                f.write('LIGHT\n')
                for lrxn in self.light:
                    lrxn_ind = str(int(self.reactions.ID[
                        self.reactions['REACTION_NAMES'] == lrxn[0]]))
                    f.write('    {} {} {}\n'.format(lrxn_ind,
                                                    lrxn[1], lrxn[2]))
                f.write(r'//' + '\n')

            if self.signals.size > 0 or self.multitoxins.size > 0:
                f.write('MET_REACTION_SIGNAL\n')
                if self.signals.size > 0:
                    sub_signals = self.signals.drop(['REACTION_NAMES', 'EXCH'],
                                                    axis='columns')
                    col_names = list(self.signals.drop(['REACTION_NAMES',
                                                        'EXCH', 'PARAMETERS'],
                                                       axis='columns').columns)
                    for idx in sub_signals.index:
                        row = sub_signals.drop(['PARAMETERS'], axis='columns').iloc[idx, :]
                        n_parms = len(sub_signals.PARAMETERS[idx])
                        curr_col_names = col_names + [str(i) for i in range(n_parms)]
                        temp_df = pd.DataFrame(columns=curr_col_names)
                        if row.loc['REACTION_NUMBER'] != 'death':
                            row.loc['REACTION_NUMBER'] = str(int(row.loc['REACTION_NUMBER']))
                        temp_df.loc[0, 'REACTION_NUMBER'] = row.loc['REACTION_NUMBER']
                        temp_df.loc[0, 'EXCH_IND'] = row.loc['EXCH_IND']
                        temp_df.loc[0, 'BOUND'] = row.loc['BOUND']
                        temp_df.loc[0, 'FUNCTION'] = row.loc['FUNCTION']
                        for i in range(n_parms):
                            temp_df.loc[0, str(i)] = sub_signals.PARAMETERS[idx][i]
                        temp_df.to_csv(f, mode='a', lineterminator = '\n', sep=' ', header=False, index=False)
                if self.multitoxins.size > 0:
                    for idx in self.multitoxins.index:
                        rxn_num = self.multitoxins.loc[idx,"REACTION_NUMBER"]
                        bound = self.multitoxins.loc[idx,"BOUND"]
                        exch_inds = ','.join([str(ind) for ind in self.multitoxins.loc[idx,"EXCH_INDS"]])
                        vmax = self.multitoxins.loc[idx, 'VMAX']
                        kms = ','.join([str(ind) for ind in self.multitoxins.loc[idx,"KMS"]])
                        hills = ','.join([str(ind) for ind in self.multitoxins.loc[idx,"HILLS"]])
                        curr_line = f"multitoxin {rxn_num} {exch_inds} {bound} {vmax} {kms} {hills}"
                        f.write(curr_line + '\n')
                f.write(r'//' + '\n')

            if self.convection_flag:
                for key, value in self.convection_parameters.items():
                    f.write(key + ' ' + str(value) + '\n')
                    f.write(r'//' + '\n')
            
            if self.multispecies_convmodel_flag:
                for key, value in self.multimodel_parameters.items():
                    f.write(key + ' ' + str(value) + '\n')
                    f.write(r'//' + '\n')

            if self.nonlinear_diffusion_flag:
                for key, value in self.nonlinear_diffusion_parameters.items():
                    f.write(key + ' ' + str(value) + '\n')
                    f.write(r'//' + '\n')
                    
            #H. Shi Feb 2023
            if self.nonlinear_diffusion_ctx_flag:
                for key, value in self.nonlinear_diffusion_ctx_parameters.items():
                    f.write(key + ' ' + str(value) + '\n')
                    f.write(r'//' + '\n')


            if self.noise_variance_flag:
                f.write('noiseVariance' + ' ' +
                        str(self.noise_variance) + '\n')
                f.write(r'//' + '\n')

            if self.neutral_drift_flag:
                f.write("neutralDrift true\n//\n")
                f.write("neutralDriftSigma " + str(self.neutralDriftSigma) + "\n//\n")

            f.write('OBJECTIVE_STYLE\n' + self.obj_style + '\n')
            f.write(r'//' + '\n')

            f.write('OPTIMIZER ' + self.optimizer + '\n')
            f.write(r'//' + '\n')


def _read_file(filename: str) -> str:
    """ helper function to read non-rectangular files.
    """
    f = open(filename, 'r')
    f_lines = f.read()
    f.close()
    return f_lines
