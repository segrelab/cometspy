'''
The Comets module serves as a Python user interface to COMETS.
For more information see https://segrelab.github.io/comets-manual/
'''

import pandas as pd
import os
import numpy as np
import re
import cobra
import io
from comets import read_file


class model:
    def __init__(self, model=None):
        self.initial_pop = [[0, 0, 0.0]]
        self.id = None
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
        self.objective = None
        self.optimizer = 'GUROBI'
        self.obj_style = 'MAXIMIZE_OBJECTIVE_FLUX'

        if model is not None:
            if isinstance(model, cobra.Model):
                self.load_cobra_model(model)
            else:  # assume it is a path
                if model[-3:] == "cmd":
                    self.read_comets_model(model)
                else:
                    self.read_cobra_model(model)

    def get_reaction_names(self):
        return(list(self.reactions['REACTION_NAMES']))

    def add_signal(self, rxn_num, exch_ind, bound,
                   function, parms):

        if str(rxn_num).lower().strip() == 'death':
            rxn_name = 'death'
            rxn_num = 'death'
        else:
            rxn_name = self.reactions.loc[self.reactions.ID == rxn_num+1,
                                          'REACTION_NAMES']
            rxn_num = str(rxn_num)

        exch_name = list(self.get_exchange_metabolites())[exch_ind-1]
        new_row = pd.DataFrame({'REACTION_NUMBER': rxn_num,
                                'EXCH_IND': exch_ind,
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
        """ toggles neutral drift to on (which is in the model file) and
        sets the demographic noise parameter neutralDriftSigma) """
        if not isinstance(neutralDriftSigma, float):
            raise ValueError("neutralDriftSigma must be a float")
        self.neutral_drift_flag = True
        self.neutralDriftSigma = neutralDriftSigma

    def add_nonlinear_diffusion_parameters(self,
                                           convNonlinDiffZero=1.,
                                           convNonlinDiffN=1.,
                                           convNonlinDiffExponent=1.,
                                           convNonlinDiffHillN=10.,
                                           convNonlinDiffHillK=0.9):
        print("Note: for non-linear diffusion parameters to function,\n"
              + "params.all_params['biomassMotionStyle'] = 'ConvNonlin' Diffusion 2D'\n"
              + "must also be set")
        for parm in [convNonlinDiffZero, convNonlinDiffN,
                     convNonlinDiffExponent, convNonlinDiffHillN,
                     convNonlinDiffHillK]:
            if not isinstance(parm, float):
                raise ValueError('all nonlinear diffusion terms must be float')
        self.nonlinear_diffusion_flag = True
        self.nonlinear_diffusion_parameters = {'convNonLinDiffZero': convNonlinDiffZero,
                                               'convNonlinDiffN': convNonlinDiffN,
                                               'convNonlinDiffExponent': convNonlinDiffExponent,
                                               'convNonlinDiffHillN': convNonlinDiffHillN,
                                               'convNonlinDiffHillK': convNonlinDiffHillK}

    def add_light(self, reaction, abs_coefficient, abs_base):
        if (reaction not in self.reactions['REACTION_NAMES'].values):
            raise ValueError('the reaction is not present in the model')
        self.light.append([reaction, abs_coefficient, abs_base])
        self.light_flag = True

    def add_convection_parameters(self, packedDensity=1.,
                                  elasticModulus=1.,
                                  frictionConstant=1.,
                                  convDiffConstant=1.):
        """ running this without named parameters sets default parameters (i.e. 1).
        Named parameters are used to specify how convection works """

        print("Note: for convection parameters to function,\n"
              + "params.all_params['biomassMotionStyle'] = 'Convection 2D'\n"
              + "must also be set")
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

    def add_noise_variance_parameter(self, noiseVariance):
        if not isinstance(noiseVariance, float):
            raise ValueError('noiseVariance must be a float')
        self.noise_variance_flag = True
        self.noise_variance = noiseVariance

    def get_exchange_metabolites(self):
        """ useful for layouts to grab these and get the set of them """
        exchmets = pd.merge(self.reactions.loc[self.reactions['EXCH'], 'ID'],
                            self.smat,
                            left_on='ID', right_on='rxn',
                            how='inner')['metabolite']
        exchmets = self.metabolites.iloc[exchmets-1]
        return(exchmets.METABOLITE_NAMES)

    def change_bounds(self, reaction, lower_bound, upper_bound):
        if reaction not in self.reactions['REACTION_NAMES'].values:
            print('reaction couldnt be found')
            return
        self.reactions.loc[self.reactions['REACTION_NAMES'] == reaction,
                           'LB'] = lower_bound
        self.reactions.loc[self.reactions['REACTION_NAMES'] == reaction,
                           'UB'] = upper_bound

    def get_bounds(self, reaction):
        if reaction not in self.reactions['REACTION_NAMES'].values:
            print('reaction couldnt be found')
            return
        lb = float(self.reactions.loc[self.reactions[
            'REACTION_NAMES'] == reaction, 'LB'])
        ub = float(self.reactions.loc[self.reactions[
            'REACTION_NAMES'] == reaction, 'UB'])
        return((lb, ub))

    def change_vmax(self, reaction, vmax):
        if reaction not in self.reactions['REACTION_NAMES'].values:
            print('reaction couldnt be found')
            return
        self.vmax_flag = True
        self.reactions.loc[self.reactions[
            'REACTION_NAMES'] == reaction, 'V_MAX'] = vmax

    def change_km(self, reaction, km):
        if reaction not in self.reactions['REACTION_NAMES'].values:
            print('reaction couldnt be found')
            return
        self.km_flag = True
        self.reactions.loc[self.reactions[
            'REACTION_NAMES'] == reaction, 'KM'] = km

    def change_hill(self, reaction, hill):
        if reaction not in self.reactions['REACTION_NAMES'].values:
            print('reaction couldnt be found')
            return
        self.hill_flag = True
        self.reactions.loc[self.reactions[
            'REACTION_NAMES'] == reaction, 'HILL'] = hill

    def read_cobra_model(self, path):
        curr_m = cobra.io.read_sbml_model(path)
        self.load_cobra_model(curr_m)

    def load_cobra_model(self, curr_m):
        self.id = curr_m.id
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
            self.smat = pd.concat([self.smat, cdf])

        self.smat = self.smat.sort_values(by=['metabolite', 'rxn'])

        # The rest of stuff
        if hasattr(curr_m, 'default_bounds'):
            self.default_bounds = curr_m.default_bounds

        obj = [str(x).split(':')[0]
               for x in reaction_list
               if x.objective_coefficient != 0][0]
        self.objective = int(self.reactions[self.reactions.
                                            REACTION_NAMES == obj]['ID'])

        if hasattr(curr_m, 'comets_optimizer'):
            self.optimizer = curr_m.comets_optimizer

        if hasattr(curr_m, 'comets_obj_style'):
            self.obj_style = curr_m.comets_obj_style

    def read_comets_model(self, path):
        self.id = os.path.splitext(os.path.basename(path))[0]

        # in this way, its robust to empty lines:
        m_f_lines = [s for s in read_file(path).splitlines() if s]
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

    def write_comets_model(self, working_dir=None):

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
            smat.to_csv(f, mode='a', header=False, index=False)
            f.write(r'//' + '\n')

            f.write('BOUNDS ' +
                    str(self.default_bounds[0]) + ' ' +
                    str(self.default_bounds[1]) + '\n')
            bnd.to_csv(f, mode='a', header=False, index=False)
            f.write(r'//' + '\n')

            f.write('OBJECTIVE\n' +
                    '    ' + str(self.objective) + '\n')
            f.write(r'//' + '\n')

            f.write('METABOLITE_NAMES\n')
            met_n.to_csv(f, mode='a', header=False, index=False)
            f.write(r'//' + '\n')

            f.write('REACTION_NAMES\n')
            rxn_n.to_csv(f, mode='a', header=False, index=False)
            f.write(r'//' + '\n')

            f.write('EXCHANGE_REACTIONS\n')
            f.write(' ' + exch_r + '\n')
            f.write(r'//' + '\n')

            if self.vmax_flag:
                f.write('VMAX_VALUES ' +
                        str(self.default_vmax) + '\n')
                Vmax.to_csv(f, mode='a', header=False, index=False)
                f.write(r'//' + '\n')

            if self.km_flag:
                f.write('KM_VALUES ' +
                        str(self.default_km) + '\n')
                Km.to_csv(f, mode='a', header=False, index=False)
                f.write(r'//' + '\n')

            if self.hill_flag:
                f.write('HILL_VALUES ' +
                        str(self.default_hill) + '\n')
                Hill.to_csv(f, mode='a', header=False, index=False)
                f.write(r'//' + '\n')

            if self.light_flag:
                f.write('LIGHT\n')
                for lrxn in self.light:
                    lrxn_ind = str(int(self.reactions.ID[
                        self.reactions['REACTION_NAMES'] == lrxn[0]]))
                    f.write('    {} {} {}\n'.format(lrxn_ind,
                                                    lrxn[1], lrxn[2]))
                f.write(r'//' + '\n')

            if self.signals.size > 0:
                f.write('MET_REACTION_SIGNAL\n')
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
                    temp_df.loc[0, 'REACTION_NUMBER'] = row.loc['REACTION_NUMBER']
                    temp_df.loc[0, 'EXCH_IND'] = row.loc['EXCH_IND']
                    temp_df.loc[0, 'BOUND'] = row.loc['BOUND']
                    temp_df.loc[0, 'FUNCTION'] = row.loc['FUNCTION']
                    for i in range(n_parms):
                        temp_df.loc[0, str(i)] = sub_signals.PARAMETERS[idx][i]
                    temp_df.to_csv(f, mode='a', sep=' ', header=False, index=False)
                f.write(r'//' + '\n')

            if self.convection_flag:
                for key, value in self.convection_parameters.items():
                    f.write(key + ' ' + str(value) + '\n')
                    f.write(r'//' + '\n')

            if self.nonlinear_diffusion_flag:
                for key, value in self.nonlinear_diffusion_parameters.items():
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
