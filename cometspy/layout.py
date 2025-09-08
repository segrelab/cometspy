import pandas as pd
import os
import numpy as np
import re
import math

from cometspy.model import model


class layout:
    '''
    object containing model and environmental definitions for a COMETS sim.
    
    The layout object is an essential object needed to start a comets
    simulation, along with a params object. For a simulation to run, the layout
    must contain model(s). Typically, the layout will be created after
    creating models, and these models will be given in a list during layout
    creation. Subsequently, methods like 
    set_specific_metabolite() will be used to create environmental media 
    definitions.  Additionally, spatial information can be assigned, such as
    the spatial dimensions in lattice units (layout.grid), spatial barriers,
    etc.
    
    Parameters
    ----------
    
    input_obj : list(cometspy.model) or str, optional
        either a list of cometspy models, or a path to a previously-written
        COMETS layout to load. 
        
    Examples
    --------
    
    >>> import cometspy as c
    >>> import cobra.test
    >>> e_cobra = cobra.test.create_test_model('textbook')
    >>> e = c.model(e_cobra)
    >>> e.open_exchanges() 
    >>> l = c.layout([e])
    >>> # use the media from the textbook in the COMETS media
    >>> for key in e_cobra.medium.keys():
    >>>     l.set_specific_metabolite(key[3:], 10.)
    >>> l.display_current_media()

    Attributes
    ----------
    models : list(cometspy.model)
        a list of cometspy models. do not modify directly. use add_model()
    grid : list of size two
        the number of boxes in the x and y direction. default = [1, 1]
    media : pandas.DataFrame
        info on met initial values and diffusion. see set_specific_metabolite
    refresh : list
        info on met added per time. see set_specific_refresh
    local_media : dict
        info on spatial-explicit media. see set_specific_metabolite_at_location
    local_refresh : dict
        info on met added per time to loc. see set_specific_refresh_at_location
    local_static : dict
        info on constant met value at loc. see set_specific_static_at_location
    default_diff_c : float
        the baseline diffusion constant for all metabolites
    barriers : list(tuple)
        list of impenetrable locations. see add_barriers
    region_map : numpy.array or None
        integer array specifying regions. see set_region_map
    region_parameters : dict
        region-specific default diffusion params. see set_region_parameters
    reactions : list
        list of extracellular reactions. see add_external_reaction
    periodic_media : list
        list of media undergoing periodic cycles. see set_global_periodic_media
    
        


    '''
    def __init__(self, input_obj : list = None):

        # define an empty layout that can be filled later
        self.models = []
        self.grid = [1, 1]
        self.media = pd.DataFrame(columns=['metabolite',
                                           'init_amount',
                                           'diff_c',
                                           'g_static',
                                           'g_static_val',
                                           'g_refresh'])

        # local_media is a dictionary with locations as keys, and as values,
        # another dict with metabolite names as keys and amounts as values
        # this information sets initial, location-specific media amounts.
        self.local_media = {}
        self.refresh = []
        self.local_refresh = {}
        self.local_static = {}
        self.initial_pop_type = "custom"  # JMC not sure purpose of this
        self.initial_pop = []
        self.all_exchanged_mets = []

        self.default_diff_c = 5.0e-6
        self.default_g_static = 0
        self.default_g_static_val = 0
        self.default_g_refresh = 0

        self.barriers = []

        self.reactions = []
        self.periodic_media = []

        self.region_map = None
        self.region_parameters = {}

        self.models_pairs_frictions =[]
        self.models_substrate_frictions =[]

        self.__local_media_flag = False
        self.__diffusion_flag = False
        self.__refresh_flag = False
        self.__static_flag = False
        self.__barrier_flag = False
        self.__region_flag = False
        self.__ext_rxns_flag = False
        self.__periodic_media_flag = False
        self.__models_frictions_flag = False
        self.__models_pairs_frictions_flag = False

        if input_obj is None:
            print('building empty layout model\nmodels will need to be added' +
                  ' with layout.add_model()')
        elif isinstance(input_obj, str):
            if not os.path.isfile(input_obj):
                raise IOError(' when running comets.layout(), input_obj' +
                              ' is a string, and therefore should be a path' +
                              ' to a layout; however, no file could be found' +
                              ' at that path destination')
            self.read_comets_layout(input_obj)
        else:
            if not isinstance(input_obj, list):
                input_obj = [input_obj]  # probably just one cobra model
            self.models = input_obj
            self.update_models()

    def set_region_parameters(self, region : int, 
                              diffusion : list, 
                              friction : float):
        """
        set regions-specific diffusion and friction values
        
        COMETS can have different regions with different substrate 
        diffusivities and frictions.  Here, you set those parameters. 

        This does not affect a simulation unless a region map is also set, 
        using the layout.set_region_map() function.
        
        Parameters
        ----------
        
        region : int
            the region (set by set_region_map) that uses the other parameters
        diffusion : list(float)
            a list containing diffusion constants (cm2/s) for each metabolite
        friction : float
            a single value for the friction observed in this region
            
        Examples
        --------
        
        >>> # this example assumes you have already created a cometspy model
        >>> # called "e" and have imported cometspy as c and numpy as np
        >>> # Here we are making a 2x2 world with two regions in which region
        >>> # 1 metabolites have diffusion constant 1.e-6 cm2/s and in region
        >>> # 2 metabolites have diffusion constant 1.e-7 cm2/s, and in both
        >>> # cases mets have default friction constant 1.
        >>> l = c.layout([e])
        >>> l.grid = [2,2]
        >>> region_map = np.array([[1, 1], [2, 2]], dtype = int)
        >>> l.set_region_map(region_map)
        >>> n_mets = length(l.media) # how many metabolites there are
        >>> l.set_region_parameters(1, [1.e-6] * n_mets, 1.)
        >>> l.set_region_parameters(1, [1.e-7] * n_mets, 1.)
        
        
        """
        if not self.__region_flag:
            print("Warning: You are setting region parameters but a region" +
                  "map has not been set. Use layout.set_region_map() or these" +
                  "parameters will be unused")
        self.region_parameters[region] = [diffusion, friction]

    def set_region_map(self, region_map : np.array):
        """
        set a region_map dictating different regions in a spatial simulation

        COMETS can have different regions with different substrate diffusivities
        and frictions.  Here, you set the map defining the regions. Specifically,
        you provide either:
        1) a numpy array whose shape == layout.grid, or
        2) a list of lists whose first length is grid[0] and second len is grid[1]

        Populating these objects should be integer values, beginning at 1 and
        incrementing only, that define the different grid areas.  These are
        intimately connected to region_parameters, which are set with
        layout.set_region_parameters()
        
        Parameters
        ----------
        region_map : numpy.array(int) or list(list)
            a 2d array of size layout.grid containing integers, or a list of
            lists which can be coerced to an integer array using np.array()
            
        Examples
        --------
        see set_region_parameters for an example
        """
        if isinstance(region_map, list):
            region_map = np.array(region_map)
        if not tuple(self.grid) == region_map.shape:
            raise ValueError("the shape of your region map must be the " +
                             "same as the grid size. specifically, \n" +
                             "tuple(layout.grid) == region_map.shape\n" +
                             "must be True after region_map = np.array(region_map)")
        self.region_map = region_map
        self.__region_flag = True

    def add_external_reaction(self,
                              rxnName : str, 
                              metabolites : list, 
                              stoichiometry : list, **kwargs):
        """
        adds an extracellular reaction to the layout
        
        External reactions are reactions not tied to a growing (living) model.
        They have rates defined in certain kwargs. There is no enyme
        concentration; it is assumed fixed. 

        Parameters
        ----------
        rxnName : str
            name of the external reaction
        metabolites : list(str)
            list of metabolite names (strs)
        stoichiometry : list(float)
            stoichiometry of the metabolites. must be same order as metabolites
        kwargs : ['Kcat' or 'Km' or 'K']
            rate parameters. Either Kcat and Km, or K must be set. 
            
        Examples
        --------
        
        >>> # set an external reaction that converts lactose into glucose + gal
        >>> # assumes import cometspy as c and a model called m has been made
        >>> rxn_name = 'lactase'
        >>> metabolites = ['lcts_e', 'glc__D_e', 'gal_e']
        >>> stoichiometry = [-1., 1., 1.]
        >>> K = 0.2 # mmol / hr
        >>> l = c.layout([m]) # m is a cometspy.model that has exchanges for the mets above
        >>> l.add_external_reaction(rxn_name, metabolites, stoichiometry,
                                    K = K)

        """

        ext_rxn = {'Name': rxnName,
                   'metabolites': metabolites,
                   'stoichiometry': stoichiometry}

        for key, value in kwargs.items():
            if key not in ['Kcat', 'Km', 'K']:
                print('Warning: Parameter ' + key + ' i not recognized and ' +
                      'will be ignored. Please set either Kcat and Km for' +
                      ' enzymatic reactions, or K for non catalyzed ones')
            else:
                ext_rxn[key] = value

        if 'Kcat' in ext_rxn and len([i for i in ext_rxn['stoichiometry']
                                      if i < 0]) > 1:
            print('Warning: Enzymatic reactions are only allowed to have'
                  + 'one reactant')

        self.reactions.append(ext_rxn)
        self.__ext_rxns_flag = True

    def set_global_periodic_media(self,
                                  metabolite : str, function : str,
                                  amplitude : float, period : float, 
                                  phase : float, offset : float):
        """
        set cyclical changes in extracellular nutrient environment

        Parameters
        ----------
        metabolite : str
            name of the metabolite under a cycle
        function : str, one of ['step', 'sin', 'cos', 'half_sin', 'half_cos']
            name of function that defines changing metabolite concentration
        amplitude : float
            amplitude of the function in mmol
        period : float
            duration of the function period in hr
        phase : float
            horizontal phase of the period (hr)
        offset : float
            vertical offset of the function (mmol)


        """

        if (metabolite not in self.media['metabolite'].values):
            raise ValueError('the metabolite is not present in the media')
        if (function not in ['step', 'sin', 'cos', 'half_sin', 'half_cos']):
            raise ValueError(function + ': function unknown')

        self.periodic_media.append([self.media.index[self.media['metabolite']
                                                     == metabolite][0],
                                    function, amplitude, period,
                                    phase, offset])
        self.__periodic_media_flag = True

    def read_comets_layout(self, input_obj : str):
        """
        load a comets layout from a file
        
        If a COMETS layout file has been previously made, this function can
        load that file into the current cometspy layout object.
        
        This is an unusual way of working with cometspy with rare uses.

        Parameters
        ----------
        input_obj : str
            the path to the existing COMETS layout file to be loaded

        """
        # .. load layout file
        f_lines = [s for s in _read_file(input_obj).splitlines() if s]
        filedata_string = os.linesep.join(f_lines)
        end_blocks = []
        for i in range(0, len(f_lines)):
            if '//' in f_lines[i]:
                end_blocks.append(i)

        # '''----------- GRID ------------------------------------------'''
        self.grid = [int(i) for i in f_lines[2].split()[1:]]
        if len(self.grid) < 2:            
            print('\n Warning: Grid must contain only two values, but it \n' +
                  ' currently contains ' + len(self.grid) +
                  '\nCheck your layout file.')

        # '''----------- MODELS ----------------------------------------'''
        '''
        Models can be specified in layout as either comets format models
        or .xml format (sbml cobra compliant)

        '''
        # right now, assume all models in layouts are strings leading to
        # comets model files

        # models need initial pop, so lets grab that first

        # '''----------- INITIAL POPULATION ----------------------------'''
        lin_initpop = re.split('initial_pop',
                               filedata_string)[0].count('\n')
        lin_initpop_end = next(x for x in end_blocks if x > lin_initpop)

        g_initpop = f_lines[lin_initpop].split()[1:]

        # TODO:  I think we should deprecate these, it makes things difficult
        # then, we could just generate these on-the-fly using the py toolbox,
        # and have the initial_pop always appear to be 'custom' type to COMETS
        # DB totally agree

        if (len(g_initpop) > 0 and g_initpop[0] in ['random',
                                                    'random_rect',
                                                    'filled',
                                                    'filled_rect',
                                                    'square']):
            self.initial_pop_type = g_initpop[0]
            self.initial_pop = [float(x) for x in g_initpop[1:]]
        else:
            self.initial_pop_type = 'custom'

            # .. local initial population values
            lin_initpop += 1

            # list of lists of lists. first level per-model, then per-location
            temp_init_pop_for_models = [[] for x in
                                        range(len(f_lines[0].split()[1:]))]

            for i in range(lin_initpop, lin_initpop_end):
                ipop_spec = [float(x) for x in
                             f_lines[i].split()]
                if len(ipop_spec)-2 != len(temp_init_pop_for_models):
                    print('\nWarning: Some initial population lines are corrupt.\n' +
                          'Check your layout file')
                if (ipop_spec[0] >= self.grid[0] or ipop_spec[1] >= self.grid[1]):
                    print('\nWarning: Some initial population values fall outside' +
                          '\nof the defined grid size. Check your layout file.')
                else:
                    for j in range(len(ipop_spec)-2):
                        if ipop_spec[j+2] != 0.0:
                            if len(temp_init_pop_for_models[j]) == 0:
                                temp_init_pop_for_models[j] = [[ipop_spec[0],
                                                                ipop_spec[1],
                                                                ipop_spec[j+2]]]
                            else:
                                temp_init_pop_for_models[j].append([ipop_spec[0],
                                                                    ipop_spec[1],
                                                                    ipop_spec[j+2]])
                                
        models = f_lines[0].split()[1:]
        if len(models) > 0:
            for i, model_path in enumerate(models):
                curr_model = model(model_path)
                # TODO: get the initial pop information for each model, because the models own that info
                curr_model.initial_pop = temp_init_pop_for_models[i]
                self.add_model(curr_model)
                self.update_models()
        else:
            print('Warning: No models in layout')

        # '''----------- MEDIA DESCRIPTION -----------------------------'''
        lin_media = re.split('world_media',
                             filedata_string)[0].count('\n') + 1
        lin_media_end = next(x for x in end_blocks if x > lin_media)

        media_names = []
        media_conc = []
        for i in range(lin_media, lin_media_end):
            metabolite = f_lines[i].split()
            media_names.append(metabolite[0])
            media_conc.append(float(metabolite[1]))

        self.media['metabolite'] = media_names
        self.media['init_amount'] = media_conc

        # '''----------- MEDIA DIFFUSION -------------------------------'''
        self.__diffusion_flag = False
        if 'DIFFUSION' in filedata_string:
            self.__diffusion_flag = True
            lin_diff = re.split('diffusion_constants',
                                filedata_string)[0].count('\n')
            lin_diff_end = next(x for x in end_blocks if x > lin_diff)

            self.default_diff_c = float(re.findall(r'\S+', f_lines[lin_diff].
                                                strip())[1])
            
            for i in range(lin_diff+1, lin_diff_end):
                diff_spec = [float(x) for x in f_lines[i].split()]
                if diff_spec[0] > len(self.media.metabolite)-1:
                    print('\n Warning: Corrupt line ' + str(i) + ' in ' +
                          'diffusion values. \nLine not written. Check your ' +
                          'layout file.')
                else:
                    self.media.loc[int(diff_spec[0]),
                                   'diff_c'] = diff_spec[1]

        self.__local_media_flag = False
        if 'MEDIA' in set(filedata_string.upper().strip().split()):
            self.__local_media_flag = True
            lin_media = [x for x in range(len(f_lines))
                         if f_lines[x].strip().split()[0].upper() ==
                         'MEDIA'][0]+1
            lin_media_end = next(x for x in end_blocks if x > lin_media)

            for i in range(lin_media, lin_media_end):
                media_spec = [float(x) for x in f_lines[i].split()]
                if len(media_spec) != len(self.media.metabolite)+2:
                    print('\nWarning: Some local "media" lines are corrupt\n ' +
                          '(wrong number of entries). Check your layout file.')
                elif (media_spec[0] >= self.grid[0] or
                      media_spec[1] >= self.grid[1]):
                    print('\nWarning: Some local "media" lines are corrupt\n' +
                          '(coordinates outside of the defined grid)\n' +
                          'Check your layout file')                    
                else:
                    loc = (int(media_spec[0]), int(media_spec[1]))
                    self.local_media[loc] = {}
                    media_spec = media_spec[2:]
                    for j in range(len(media_spec)):
                        if media_spec[j] != 0:
                            self.local_media[loc][
                                self.all_exchanged_mets[j]] = media_spec[j]

        # '''----------- MEDIA REFRESH----------------------------------'''
        # .. global refresh values
        self.__refresh_flag = False
        if 'REFRESH' in filedata_string.upper():
            self.__refresh_flag = True
            lin_refr = re.split('REFRESH',
                                filedata_string.upper())[0].count('\n')
            lin_refr_end = next(x for x in end_blocks if x > lin_refr)

            g_refresh = [float(x) for x in f_lines[lin_refr].split()[1:]]

            if len(g_refresh) != len(media_names):
                print('\nWarning: Some local refresh lines are corrupt\n ' +
                      '(wrong number of entries). Check your layout file.')
                
            else:
                self.media['g_refresh'] = g_refresh

            # .. local refresh values
            lin_refr += 1

            for i in range(lin_refr, lin_refr_end):
                refr_spec = [float(x) for x in f_lines[i].split()]
                if len(refr_spec) != len(self.media.metabolite)+2:
                    print('\nWarning: Some local "refresh" lines are corrupt\n ' +
                          '(wrong number of entries). Check your layout file.')
                    
                elif (refr_spec[0] >= self.grid[0] or
                      refr_spec[1] >= self.grid[1]):
                    print('\nWarning: Some local "refresh" lines are corrupt\n' +
                          '(coordinates outside of the defined grid)\n' +
                          'Check your layout file')
                else:
                    loc = (int(refr_spec[0]), int(refr_spec[1]))
                    self.local_refresh[loc] = {}
                    refr_spec = refr_spec[2:]
                    for j in range(len(refr_spec)):
                        if refr_spec[j] != 0:
                            self.local_refresh[loc][
                                self.all_exchanged_mets[j]] = refr_spec[j]

        # region-based information (substrate diffusivity,friction, layout)
        self.__region_flag = False

        if 'SUBSTRATE_LAYOUT' in filedata_string.upper():
            lin_substrate = re.split('SUBSTRATE_LAYOUT',
                                     filedata_string.upper())[0].count('\n')
            lin_substrate_end = next(x for x in end_blocks if x > lin_substrate)
            region_map_data = []
            for i in range(lin_substrate+1, lin_substrate_end):
                region_map_data.append([int(x) for x in f_lines[i].split()])
                region_map_data = np.array(region_map_data, dtype=int)
                if region_map_data.shape != tuple(self.grid):
                    print('\nWarning: Some substrate_layout lines are ' +
                          ' longer or shorter than the grid width, or there are more' +
                          ' lines than the grid length. Check your layout file.' +
                          'Check your layout file.')
                self.__region_flag = True
                self.region_map = region_map_data

        if 'SUBSTRATE_DIFFUSIVITY' in filedata_string.upper():
            lin_substrate = re.split('SUBSTRATE_DIFFUSIVITY',
                                     filedata_string.upper())[0].count('\n')
            lin_substrate_end = next(x for x in end_blocks if x > lin_substrate)
            self.region_parameters = {}
            region = 1
            for i in range(lin_substrate+1, lin_substrate_end):
                self.region_parameters[region] = [None, None]
                self.region_parameters[region][0] = [float(x)
                                                     for x in f_lines[i].split()]
                if len(self.region_parameters[
                        region][0]) != len(self.media.metabolite):
                    print('\nWarning: Some substrate_diffusivity lines are ' +
                          ' longer or shorter than the number of metabolites\n' +
                          'Check your layout file.')
                region += 1

        if 'SUBSTRATE_FRICTION' in filedata_string.upper():
            lin_substrate = re.split('SUBSTRATE_FRICTION',
                                     filedata_string.upper())[0].count('\n')
            lin_substrate_end = next(x for x in end_blocks if x > lin_substrate)
            region = 1
            for i in range(lin_substrate+1, lin_substrate_end):
                self.region_parameters[region][1] = float(f_lines[i].split()[0])
                region += 1

        # '''----------- STATIC MEDIA ----------------------------------'''
        # .. global static values
        self.__static_flag = False
        if 'STATIC' in filedata_string.upper():
            self.__static_flag = True
            lin_static = re.split('STATIC',
                                  filedata_string.upper())[0].count('\n')
            lin_stat_end = next(x for x in end_blocks if x > lin_static)

            g_static = [float(x) for x in f_lines[lin_static].split()[1:]]

            if len(g_static) != 2*len(self.media.metabolite):
                print('\nWarning: Wrong number of global static values .\n' +
                      'Check your layout file')
                
            else:
                self.media.loc[:, 'g_static'] = [int(x)
                                                 for x in g_static[0::2]]
                self.media.loc[:, 'g_static_val'] = [float(x) for x in
                                                     g_static[1::2]]

            # .. local static values
            lin_static += 1
            for i in range(lin_static, lin_stat_end):
                stat_spec = [float(x) for x in f_lines[i].split()]
                if len(stat_spec) != (2*len(self.media.metabolite))+2:
                    print('\nWarning: Wrong number of local static values at some\n' +
                          'lines. Check your layout file.')
                elif (stat_spec[0] >= self.grid[0] or
                      stat_spec[1] >= self.grid[1]):
                    print('\nWarning: Some local "static" lines have coordinates\n' +
                          'that fall outside of the defined grid. Check your layout file')
                else:
                    loc = (int(stat_spec[0]), int(stat_spec[1]))
                    self.local_static[loc] = {}
                    stat_spec = stat_spec[2:]
                    for j in range(int(len(stat_spec)/2)):
                        if stat_spec[j*2] != 0:
                            self.local_static[loc][
                                self.all_exchanged_mets[j]] = stat_spec[j*2+1]
                            
        # '''----------- CHEMOTAXIS COEFFICIENT ----------------------------------'''
        # .. chemotaxis coefficient and nutrient parameter
        #filedata_string is one big string with all the data from the file
        self.__chemotaxis_coeff_flag = False
        if 'MODEL_MEDIA_CHEMOTAXIS_COEFFS' in filedata_string.upper().strip().split():
            #splits at keyword
            self.__chemotaxis_coeff_flag = True
            #there is chemotaxis
            #split at keyword
            lin_ctx = re.split('CHEMOTAXIS',
                                filedata_string.upper())[0].count('\n')
            lin_ctx_end = next(x for x in end_blocks if x > lin_ctx)

            g_ctx = [float(x) for x in f_lines[lin_ctx].split()[1:]]

            if len(g_ctx) != 3:
                print('\nWarning: Some chemotaxis coeff lines are corrupt\n ' +
                      '(wrong number of entries). Check your layout file.')
                
            else:
                #self.chemotaxis = pd.DataFrame(columns =  ['model', 'media', 'coeff', 'n_params'])
        #H. Shi: model, media, coeff, nutrient params
                self.chemotaxis['model'] = [int(x) for x in g_ctx[0::4]] #first element then every fourth element?
                self.chemotaxis['media'] = [int(x) for x in g_ctx[1::4]]
                self.chemotaxis['coeff'] = [float(x) for x in g_ctx[2::4]]
                #self.chemotaxis['n_params'] = [float(x) for x in g_ctx[3:4]]
        '''
                self.chemotaxis['model'] = int(g_ctx[0])
                self.chemotaxis['media'] = int(g_ctx [1])
                self.chemotaxis['coeff'] = g_ctx[2]
                self.chemotaxis['n_params'] =  g_ctx[3]
        '''

        # '''----------- CHEMOTAXIS NUTRIENT PARAMETER ----------------------------------'''
        # .. nutrient parameter
        #filedata_string is one big string with all the data from the file
        self.__chemotaxis_param_flag = False
        if 'MODEL_MEDIA_NUTRIENT_PARAMS' in filedata_string.upper().strip().split():
            #splits at keyword
            self.__chemotaxis_param_flag = True
            #there is chemotaxis
            #split at keyword
            lin_ctx_p = re.split('NUTRIENT',
                                filedata_string.upper())[0].count('\n')
            lin_ctx_p_end = next(x for x in end_blocks if x > lin_ctx_p)

            g_ctx_p = [float(x) for x in f_lines[lin_ctx_p].split()[1:]]

            if len(g_ctx_p) != 3:
                print('\nWarning: Some chemotaxis nutrient params lines are corrupt\n ' +
                      '(wrong number of entries). Check your layout file.')
                
            else:
                #self.chemotaxis = pd.DataFrame(columns =  ['model', 'media', 'coeff', 'n_params'])
        #H. Shi: model, media, coeff, nutrient params
                #self.chemotaxis['model'] = [int(x) for x in g_ctx[0::4]] #first element then every fourth element?
                #self.chemotaxis['media'] = [int(x) for x in g_ctx[1::4]]
                #self.chemotaxis['coeff'] = [float(x) for x in g_ctx[2::4]]
                self.chemotaxis['n_params'] = [float(x) for x in g_ctx_p[2:3]]
        '''
                self.chemotaxis['model'] = int(g_ctx[0])
                self.chemotaxis['media'] = int(g_ctx [1])
                self.chemotaxis['coeff'] = g_ctx[2]
                self.chemotaxis['n_params'] =  g_ctx[3]
        '''
        
    def get_model_ids(self) -> list:
        """
        returns a list of the ids of the models in this layout.

        """
        ids = [x.id for x in self.models]
        return(ids)

    def write_necessary_files(self, working_dir : str, to_append = ""):
        """
        writes the layout and the model files to file
        
        This method is called by comets.run(). Alternatively it may be called
        directly if one wishes to inspect comets layout and model files.
        
        A path must be specified if called directly.
        
        Note: the layout file will be called ".current_layout" and therefore
        be hidden from file browsers by default.

        Parameters
        ----------
        working_dir : str, 
            The directory where the files will be written.
        to_append : str, 
            String to append to written filenames

        """
        self.__check_if_initial_pops_in_range()
        self.write_layout(working_dir, to_append)
        self.write_model_files(working_dir)

    def write_model_files(self, working_dir=""):
        '''writes each model file'''
        for m in self.models:
            m.write_comets_model(working_dir)
            
    def delete_model_files(self, working_dir):
        """
        deletes model files in specified directory
        """
        for m in self.models:
            m.delete_comets_model(working_dir)

    def display_current_media(self):
        """
        a utility function to show current non-zero media amounts

        """
        print(self.media[self.media['init_amount'] != 0.0])

    def add_barriers(self, barriers : list):
        """
        adds impenetrable barriers to the layout at specified locations
        
        Barriers can be used in COMETS to create impenetrable spatial
        locations, for example to simulate rocks or to make a rounded
        simulation. Neither biomass nor metabolites can diffuse through
        barriers.
        
        In COMETS, locations are zero-ordered.

        Parameters
        ----------
        
        barriers : list(tuple)
            A list containing tuples of integers of x,y locations.

        Examples
        --------
        
        >>> # make a diagonal line of barriers
        >>> l = c.layout([model1, model2]) # assumes model1,2 are made
        >>> l.grid = [5,5]
        >>> l.add_barriers([(0,0),(1,1),(2,2),(3,3),(4,4)])

        """
        # first see if they provided only one barrier not in a nested list, and if
        # so, put it into a list
        if len(barriers) == 2:
            if isinstance(barriers[0], int):
                barriers = [barriers]
        # now check each barrier and make sure it has 2 ints that fit within the grid size
        for b in barriers:
            try:
                if len(b) != 2 or b[0] >= self.grid[0] or b[1] >= self.grid[1]:
                    raise ValueError
                self.barriers.append((int(b[0]), int(b[1])))
            except ValueError:
                print('ERROR ADDING BARRIERS in add_barriers\n')
                print("expecting barriers to be a list of tuples of coordinates which"
                      + " fit within the current grid")
                print("  such as  layout.grid = [5,5]")
                print("           barriers = [(0,0),(1,1),(2,2),(4,4)]")
                print("           layout.add_barriers(barriers)")
        if len(self.barriers) > 0:
            self.__barrier_flag = True
            self.barriers = list(set(self.barriers))

    def set_metabolite_diffusion(self, diffusion_constant : float):
        """
        sets the default diffusion constant for all metabolites
        
        Specific diffusion constants can be set subsequent to setting this.

        Parameters
        ----------
        diffusion_constant : float
            the diffusion constant, in cm2 / s, that all metabolites use

        """
        try:
            if not isinstance(diffusion_constant, float):
                raise ValueError
        except ValueError:
            print("ERROR, diffusion_constant must be a float")
            return
        print("Warning: set_metabolite_diffusion overwrites all diffusion constants\nthis should be performed prior to setting specific metabolite diffusion")
        self.default_diff_c = diffusion_constant
        self.media.loc[:, "diff_c"] = diffusion_constant
        self.__diffusion_flag = True
        
    def __check_if_diffusion_flag_should_be_set(self):
        """ Internal function.  If someone directly changes the diff_c in the
        media dataframe, cometspy will be unaware and not set __diffusion_flag = True,
        rendering that change moot.  This function is run when doing write_layout() 
        to flip the flag if necessary """
        if len(np.unique(self.media.diff_c)) > 1:
            self.__diffusion_flag = True
        elif self.media.diff_c[0] != self.default_diff_c:
            self.__diffusion_flag = True
            
    def set_specific_metabolite_diffusion(self, met : str, 
                                          diffusion_constant : float):
        """
        sets the diffusion constant of a specific metabolite

        Parameters
        ----------
        met : str
            name of the metabolite
        diffusion_constant : float
            the diffusion constant, in cm2/s, of the named metabolite

        Examples
        --------
        >>> l = c.layout([m])
        >>> l.set_specific_metabolite_diffusion("ac_e", 5.2e-6)
        >>> l.set_specific_metabolite_diffusion("gal_e", 4.e-6)


        """
        try:
            if not isinstance(diffusion_constant, float):
                raise ValueError
        except ValueError:
            print("ERROR, diffusion_constant must be a float")
            return
        try:
            if len(self.media.loc[self.media.metabolite == met, "diff_c"]) != 1:
                raise ValueError
        except ValueError:
            print("ERROR, met not in layout.media.metabolite list\ncannot change metabolite diffusion constant")
        self.media.loc[self.media.metabolite == met, "diff_c"] = diffusion_constant
        self.__diffusion_flag = True
        
    def set_specific_metabolite(self, met : str, amount : float, 
                                static : bool = False):
        """
        sets the initial (or constant) value for a metabolite across space
        
        This is the most common way to set a metabolite concentration. It 
        sets the initial value of the named metabolite to 'amount' in every
        location. Optionally, if static = True, this concentration is fixed
        during the simulation.

        Parameters
        ----------
        met : str
            name of the metabolite
        amount : float
            initial value for the metabolite in each box, in mmol
        static : bool, optional
            DEFAULT = False. If True, the 'amount' is fixed over time.

        Examples
        --------
        
        >>> l = c.layout([model])
        >>> l.set_specific_metabolite("glc__D_e", 0.005)
        >>> l.set_specific_metabolite("o2_e", 15, static = True)
        >>> # a common thing to do is populate this from the cobra model, as 
        >>> # shown here:
        >>> import cobra.test
        >>> import cometspy as c
        >>> cobra_model = cobra.test.create_test_model("textbook")
        >>> model = c.model(cobra_model)
        >>> model.open_exchanges()
        >>> l = c.layout([model])
        >>> for key, value in cobra_model.medium.items():
        >>>     l.set_specific_metabolite(key[3:], value)


        """
        if met in set(self.media['metabolite']):
            self.media.loc[self.media['metabolite'] == met,
                           'init_amount'] = amount
            if static:
                self.media.loc[self.media['metabolite'] == met,
                               'g_static'] = 1
                self.media.loc[self.media['metabolite'] == met,
                               'g_static_val'] = amount
        else:
            newrow = {'metabolite': met,
                      'g_refresh': self.default_g_refresh,
                      'g_static': 1 if static else self.default_g_static,
                      'g_static_val': amount if static else self.default_g_static_val,
                      'init_amount': amount,
                      'diff_c': self.default_diff_c}
            newrow = pd.DataFrame([newrow], columns=newrow.keys())
            # concat the row only if table isn't empty (avoids FutureWarning)
            if self.media.shape[0] == 0:
                self.media = newrow
            else:
                self.media = pd.concat([self.media,
                                        newrow],
                                        axis=0, sort=False, ignore_index=True)
            print('Warning: The added metabolite (' + met + ') is not ' +
                  'able to be taken up by any of the current models')

    def set_specific_metabolite_at_location(self, met : str, 
                                            location : tuple, 
                                            amount : float):
        """
        set an initial value of a metabolite at a specific spatial location
        
        In contrast to set_specific_metabolite, which sets initial values 
        universally over space, this method sets a specific metabolite 
        initial value at one location only.

        Parameters
        ----------
        met : str
            name of the metabolite
        location : tuple
            a tuple of integers specifying a specific location, e.g. (3,2)
        amount : float
            the initial value for the metabolite at the location, in mmol

        Examples
        --------
        >>> l = c.layout([model])
        >>> l.grid = [10,1]
        >>> l.set_specific_metabolite_at_location("glc__D_e", (9,0), 0.005)

        Notes
        -----
        Remember that COMETS locations are zero-ordered.
        
        """
        if met not in self.all_exchanged_mets:
            raise Exception('met is not in the list of exchangeable mets')
        self.__local_media_flag = True
        if location not in list(self.local_media.keys()):
            self.local_media[location] = {}
        self.local_media[location][met] = amount

    def set_specific_refresh(self, met : str, amount : float):
        """
        set the amount to increment a metabolite continuously through time
        
        metabolites can be "refreshed" or "dripped in" during a COMETS 
        simulation. This sets the amount added each hour, which is divided
        up by the number of time steps. It applies that refresh to all boxes
        in space.
        

        Parameters
        ----------
        met : str
            name of the metabolite refreshed
        amount : float
            the amount added to (or subtracted from) each box, in mmol / hr.


        """
        if met in self.media['metabolite'].values:
            self.media.loc[self.media['metabolite'] == met,
                           'g_refresh'] = amount
            self.__refresh_flag = True
        else:
            raise Exception("the specified metabolite " + met +
                  "is not in the medium; add it first!")

    def set_specific_refresh_at_location(self, met : str, 
                                         location : tuple, 
                                         amount : float):
        """
        sets the amount to adjust a metabolite at a specific location in time
        
        Like set_specific_refresh, but for a specific lacation. This 
        Method is used to specify how much a metabolite is increased (or 
        decreased) each hour in a specific location.

        Parameters
        ----------
        met : str
            the name of the metabolite
        location : tuple
            the (x, y) location to adjust, with x and y of type int
        amount : float
            the amount added to (or subtracted from) the box, in mmol / hr


        """
        if met not in self.all_exchanged_mets:
            raise Exception('met is not in the list of exchangeable mets')
        self.__refresh_flag = True
        if location not in list(self.local_refresh.keys()):
            self.local_refresh[location] = {}
        self.local_refresh[location][met] = amount

    def set_specific_static(self, met : str, amount : float):
        """
        sets a metabolite to a fixed value in every location

        Parameters
        ----------
        met : str
            name of the metabolite
        amount : float
            the amount, in mmol, that is in each box at each time step

        """
        if met in self.media['metabolite'].values:
            self.media.loc[self.media['metabolite'] == met,
                           'g_static'] = 1
            self.media.loc[self.media['metabolite'] == met,
                           'g_static_val'] = amount
            self.__static_flag = True
        else:
            print("the specified metabolite " + met +
                  "is not in the medium; add it first!")

    def set_specific_static_at_location(self, met : str, 
                                        location : tuple, 
                                        amount : float):
        """
        sets a metabolite to a fixed value at a specific location

        Parameters
        ----------
        met : str
            name of the metabolite
        location : tuple
            the (x, y) location of the fixed metabolite amount. x, y are int
        amount : float
            the amount, in mmol, that is in each box at each time step


        """
        if met not in self.all_exchanged_mets:
            raise Exception('met is not in the list of exchangeable mets')
        self.__static_flag = True
        if location not in list(self.local_static.keys()):
            self.local_static[location] = {}
        self.local_static[location][met] = amount

    def set_models_substrate_frictions(self, 
                                        friction : list):
        """
        sets the values of the friction paremeters between the organisms and substrate

        Parameters
        ----------
        model : int
            the order of the model
        location : tuple
            the (x, y) location of the fixed metabolite amount. x, y are int
        amount : float
            the amount, in mmol, that is in each box at each time step


        """
        self.__models_frictions_flag = True
        self.models_substrate_frictions=friction
           
    def set_models_pairs_frictions(self, 
                                        pairs_frictions : list):
        """
        sets the values of the friction paremeters between the organisms and substrate

        Parameters
        ----------
        model : int
            the order of the model
        location : tuple
            the (x, y) location of the fixed metabolite amount. x, y are int
        amount : float
            the amount, in mmol, that is in each box at each time step


        """
        self.__models_pairs_frictions_flag = True
        self.models_pairs_frictions=pairs_frictions

    def add_typical_trace_metabolites(self, amount : float=1000.0, 
                                      static : bool =True):
        """
        adds common BIGG-style metabolites to the environment. This only 
        works if your models use metabolites which are formatted like 'ca2_e'
        and then should still only be used as a starting point. 

        Parameters
        ----------
        amount : float, optional
            the amount, in mmol, of these metabolites. Default is 1000.0.
        static : bool, optional
            The default is True. If false, these metabolites are depletable. 


        """
        trace_metabolites = ['ca2_e',
                             'cl_e',
                             'cobalt2_e',
                             'cu2_e',
                             'fe2_e',
                             'fe3_e',
                             'h_e',
                             'k_e',
                             'h2o_e',
                             'mg2_e',
                             'mn2_e',
                             'mobd_e',
                             'na1_e',
                             'ni2_e',
                             'nh4_e',
                             'o2_e',
                             'pi_e',
                             'so4_e',
                             'zn2_e']

        for met in trace_metabolites:
            if met in set(self.media['metabolite']):
                self.media.loc[self.media['metabolite'] == met,
                               'init_amount'] = amount
                if static:
                    self.media.loc[self.media['metabolite'] == met,
                                   'g_static'] = 1
                    self.media.loc[self.media['metabolite'] == met,
                                   'g_static_val'] = amount
                
            else:
                newrow = {'metabolite': met,
                          'g_refresh': self.default_g_refresh,
                          'g_static': 1 if static else self.default_g_static,
                          'g_static_val': amount if static else self.default_g_static_val,
                          'init_amount': amount,
                          'diff_c': self.default_diff_c}
                newrow = pd.DataFrame([newrow], columns=newrow.keys())
                self.media = pd.concat([self.media,
                                        newrow],
                                       axis=0, sort=False)
        self.media = self.media.reset_index(drop=True)

    def write_layout(self, working_dir : str, to_append = ""):
        """
        writes just the COMETS layout file to the supplied path

        Parameters
        ----------
        working_dir : str
            the path to the directory where .current_layout will be written


        """
        # right now we only check if a user manually set a diff_c.  Ideally,
        # we should check for manual changes to everything. Alternatively,
        # we should print all blocks no matter what. 
        self.__check_if_diffusion_flag_should_be_set()
        outfile = working_dir + ".current_layout" + to_append
        if os.path.isfile(outfile):
            os.remove(outfile)

        lyt = open(outfile, 'a')
        self.__write_models_and_world_grid_chunk(lyt)
        self.__write_media_chunk(lyt)
        self.__write_diffusion_chunk(lyt)
        self.__write_local_media_chunk(lyt)
        self.__write_refresh_chunk(lyt)
        self.__write_static_chunk(lyt)
        self.__write_barrier_chunk(lyt)
        self.__write_regions_chunk(lyt)
        self.__write_periodic_media_chunk(lyt)
        self.__write_models_frictions_chunk(lyt)
        self.__write_models_pairs_frictions_chunk(lyt)
        lyt.write(r'  //' + '\n')

        self.__write_initial_pop_chunk(lyt)
        self.__write_ext_rxns_chunk(lyt)
        lyt.close()

    def __write_models_and_world_grid_chunk(self, lyt):
        """ writes the top 3 lines  to the open lyt file"""

        model_file_line = "{}.cmd".format(".cmd ".join(self.get_model_ids())).split(" ")
        model_file_line = "".join(["./" + _ + " " for _ in model_file_line])
        model_file_line = "model_file " + model_file_line + "\n"
        lyt.write(model_file_line)
        lyt.write('  model_world\n')

        lyt.write('    grid_size ' +
                  ' '.join([str(x) for x in self.grid]) + '\n')

    def __write_media_chunk(self, lyt):
        """ used by write_layout to write the global media information to the
        open lyt file """
        lyt.write('    world_media\n')
        for i in range(0, len(self.media)):
            lyt.write('      ' + self.media.metabolite[i] +
                      ' ' + str(self.media.init_amount[i]) + '\n')
        lyt.write(r'    //' + '\n')

    def __write_local_media_chunk(self, lyt):
        """ used by write_layout to write the location-specific initial
        metabolite data"""
        if self.__local_media_flag:
            lyt.write('    media\n')
            locs = list(self.local_media.keys())
            for loc in locs:
                # this chunk goes in order, not by name, so must get met number
                # for each location, make a list with zeros for each met. Put
                # non-zero numbers where the self.local_media tells us to
                met_amounts_in_order = [0] * len(self.all_exchanged_mets)
                for met in list(self.local_media[loc].keys()):
                    met_amounts_in_order[
                        self.__get_met_number(met)] = self.local_media[loc][met]
                lyt.write('      ')
                lyt.write('{} {} '.format(loc[0], loc[1]))
                lyt.write(' '.join(str(x) for x in met_amounts_in_order))
                lyt.write('\n')
            lyt.write('    //\n')

    def __write_refresh_chunk(self, lyt):
        if self.__refresh_flag:
            lyt.write('    media_refresh ' +
                      ' '.join([str(x) for x in self.media.
                                g_refresh.tolist()]) +
                      '\n')
            locs = list(self.local_refresh.keys())
            if len(locs) > 0:
                for loc in locs:
                    met_amounts_in_order = [0] * len(self.all_exchanged_mets)
                    for met in list(self.local_refresh[loc].keys()):
                        met_amounts_in_order[self.__get_met_number(met)] = self.local_refresh[loc][met]
                    met_amounts_in_order.insert(0, loc[1])
                    met_amounts_in_order.insert(0, loc[0])
                    lyt.write('      ' +
                              ' '.join([str(x) for x in met_amounts_in_order]) +
                              '\n')
            lyt.write(r'    //' + '\n')

    def __write_static_chunk(self, lyt):
        if self.__static_flag:
            g_static_line = [None]*(len(self.media)*2)
            g_static_line[::2] = self.media.g_static
            g_static_line[1::2] = self.media.g_static_val
            lyt.write('    static_media ' +
                      ' '.join([str(x) for x in g_static_line]) + '\n')
            locs = list(self.local_static.keys())
            if len(locs) > 0:
                for loc in locs:
                    # this is 2 * len because there is a pair of values for each met
                    # the first value is a flag--0 if not static, 1 if static
                    # the second value is the amount if it is static
                    met_amounts_in_order = [0] * 2 * len(self.all_exchanged_mets)
                    for met in list(self.local_static[loc].keys()):
                        met_amounts_in_order[self.__get_met_number(met) * 2] = 1
                        met_amounts_in_order[
                            self.__get_met_number(met) * 2 + 1] = self.local_static[
                                loc][met]
                    met_amounts_in_order.insert(0, loc[1])
                    met_amounts_in_order.insert(0, loc[0])
                    lyt.write('      ' +
                              ' '.join([str(x) for x in
                                        met_amounts_in_order]) +
                              '\n')
            lyt.write(r'    //' + '\n')

    def __write_diffusion_chunk(self, lyt):
        """ used by write_layout to write the metab-specific
        diffusion data to the open lyt file """

        if self.__diffusion_flag:
            lyt.write('    diffusion_constants ' +
                      str(self.default_diff_c) +
                      '\n')
            for i in range(0, len(self.media)):
                if not math.isnan(self.media.diff_c[i]):
                    lyt.write('      ' + str(i) + ' ' +
                              str(self.media.diff_c[i]) + '\n')
            lyt.write(r'    //' + '\n')

    def __write_barrier_chunk(self, lyt):
        """ used by write_layout to write the barrier section to the open lyt file """
        if self.__barrier_flag:
            lyt.write('    barrier\n')
            for barrier in self.barriers:
                lyt.write('      {} {}\n'.format(barrier[0], barrier[1]))
            lyt.write('    //\n')

    def __write_ext_rxns_chunk(self, lyt):
        """ used by write_layout to write the external reactions section
        to the open lyt file
        """
        reactants = []
        enzymes = []
        products = []

        if self.__ext_rxns_flag:
            for i, rxn in enumerate(self.reactions):

                current_reactants = [self.media.index[
                    self.media['metabolite'] ==
                    rxn['metabolites'][k]].tolist()[0]+1
                                     for k in range(len(rxn['metabolites']))
                                     if rxn['stoichiometry'][k] < 0]

                current_products = [self.media.index[
                    self.media['metabolite'] ==
                    rxn['metabolites'][k]].tolist()[0]+1
                                     for k in range(len(rxn['metabolites']))
                                     if rxn['stoichiometry'][k] > 0]

                current_react_stoich = [k for k in rxn['stoichiometry']
                                        if k < 0]

                current_prod_stoich = [k for k in rxn['stoichiometry']
                                       if k > 0]

                for ind, k in enumerate(current_reactants):
                    if ind == 0:

                        cl = ('        ' + str(i+1)             # reaction
                              + ' ' + str(k)                     # metabolite
                              + ' ' + str(-current_react_stoich[ind])  # stoich
                              + ' '
                              + str([rxn['K'] if 'K' in rxn else rxn['Km']][0])
                              + '\n')
                        reactants.append(cl)
                    else:
                        cl = ('        ' + str(i+1)
                              + ' ' + str(k)
                              + ' ' + str(-current_react_stoich[ind])
                              + ' ' + '\n')
                        reactants.append(cl)

                for ind, k in enumerate(current_products):
                    cl = ('        ' + str(i+1)
                          + ' ' + str(k)
                          + ' ' + str(current_prod_stoich[ind])
                          + ' ' + '\n')
                    products.append(cl)

                if 'Kcat' in rxn:
                    cl = ('        ' + str(i+1)
                          + ' ' + str(rxn['Kcat'])
                          + '\n')
                    enzymes.append(cl)

            # write the reaction lines
            lyt.write('reactions\n')
            lyt.write('    reactants\n')
            for i in reactants:
                lyt.write(i)

            lyt.write('    enzymes\n')
            for i in enzymes:
                lyt.write(i)

            lyt.write('    products\n')
            for i in products:
                lyt.write(i)
            lyt.write('//\n')

    def __write_periodic_media_chunk(self, lyt):
        """ used by write_layout to write the periodic media
        """
        if self.__periodic_media_flag:
            lyt.write('    periodic_media global\n')
            for media in self.periodic_media:
                lyt.write('        {} {} {} {} {} {}\n'.
                          format(media[0], media[1], media[2],
                                 media[3], media[4], media[5]))
            lyt.write('    //\n')

    def __write_regions_chunk(self, lyt):
        """ used by write_layout to write the regions section to the open lyt file
        specifically this section includes "substrate_diffusivity" "substrate_friction"
        and "substrate_layout".
        """
        if self.__region_flag:
            keys = list(self.region_parameters.keys())
            keys.sort()
            lyt.write('    substrate_diffusivity\n')
            for key in keys:
                diff = [str(x) for x in self.region_parameters[key][0]]
                line = "    " + "    ".join(diff) + "\n"
                lyt.write(line)
            lyt.write("    //\n")
            lyt.write("    substrate_friction\n")
            for key in keys:
                fric = self.region_parameters[key][1]
                line = "    " + str(fric) + "\n"
                lyt.write(line)
            lyt.write("    //\n")
            lyt.write("    substrate_layout\n")
            for i in range(self.region_map.shape[0]):
                for j in range(self.region_map.shape[1]):
                    lyt.write("    ")
                    lyt.write(str(self.region_map[i, j]))
                lyt.write("\n")
            lyt.write("    //\n")

    def __write_initial_pop_chunk(self, lyt):
        """ writes the initial pop to the open
        lyt file and adds the closing //s """
        if (self.initial_pop_type == 'custom'):
            lyt.write('  initial_pop\n')
            for i in self.initial_pop:
                lyt.write('    ' + str(int(i[0])) + ' ' + str(int(i[1])) +
                          ' ' + ' '.join([str(x) for x in i[2:]]) +
                          '\n')
        else:
            # TODO: test this part and fix, probably not functional currently
            lyt.write('  initial_pop ' + self.initial_pop_type +
                      ' '.join([str(x) for x in self.initial_pop]) +
                      '\n')
        lyt.write(r'  //' + '\n')
        lyt.write(r'//' + '\n')

    def __write_models_frictions_chunk(self, lyt):
        """ writes the models frictions to the open
        lyt file and adds the closing //s """
        if self.__models_frictions_flag:
            lyt.write('  modelsFriction\n')
            for i in range(len(self.models_substrate_frictions)):
                lyt.write('    ' + str(i)+' ' +str(self.models_substrate_frictions[i]) + '\n')
            lyt.write(r'  //' + '\n')

    def __write_models_pairs_frictions_chunk(self, lyt):
        """ writes the intermodels pairs frictions to the open
        lyt file and adds the closing //s """
        if self.__models_pairs_frictions_flag:
            lyt.write('  interModelPairsFriction\n')
            for i in range(len(self.models_pairs_frictions)):
                for j in range(len(self.models_pairs_frictions)):
                    lyt.write('    ' + str(i) + ' ' +str(j) + ' ' + str(self.models_pairs_frictions[i][j]) + '\n')
            lyt.write(r'  //' + '\n')

    def update_models(self):
        """ updates layout properties when models are added / removed
        
        This is usually just called internally. However, if you manually 
        change layout.models, then this method should be called to make other
        necessary changes. 
        
        """
        self.__build_initial_pop()
        self.__build_exchanged_mets()
        self.__add_new_mets_to_media()

    def __build_initial_pop(self):
        # This counts how many models there are.  then it goes through
        # each model, and makes a new initial pop line of the right length
        n_models = len(self.models)
        initial_pop = []
        for i, m in enumerate(self.models):
            if not isinstance(m.initial_pop[0], list):
                m.initial_pop = [m.initial_pop]
            for pop in m.initial_pop:
                curr_line = [0] * (n_models + 2)
                curr_line[0] = pop[0]
                curr_line[1] = pop[1]
                curr_line[i+2] = pop[2]
                initial_pop.append(curr_line)
        self.initial_pop = initial_pop
        self.__resolve_initial_pop_line_dups()
        
    def __resolve_initial_pop_line_dups(self):
        # sometimes each model has the same x,y, this resolves that
        init_pop = self.initial_pop
        init_pop_dict = {}
        for row in init_pop:
            if tuple(row[0:2]) in init_pop_dict.keys():
                init_pop_dict[tuple(row[0:2])] += np.array(row[2:])
            else:
                init_pop_dict[tuple(row[0:2])] = np.array(row[2:])
        init_pop_fixed = []
        for key, value in init_pop_dict.items():
            loc = list(key)
            loc.extend(list(value))
            init_pop_fixed.append(loc)
        self.initial_pop = init_pop_fixed

    def __add_new_mets_to_media(self):
        # usually run right after build_exchange mets, to add any new mets
        # to the media data.frame

        for met in self.all_exchanged_mets:
            if met not in self.media['metabolite'].values:
                new_row = pd.DataFrame.from_dict({'metabolite': [met],
                                                  'init_amount': [0.0],
                                                  'diff_c': [self.default_diff_c],
                                                  'g_static': [self.default_g_static],
                                                  'g_static_val': [
                                                      self.default_g_static_val],
                                                  'g_refresh': [self.default_g_refresh]})
                # concat the row only if table isn't empty (avoids FutureWarning)
                if self.media.shape[0] == 0:
                    self.media = new_row
                else:
                    self.media = pd.concat([self.media, new_row],
                                            ignore_index=True, sort=True)

    def __build_exchanged_mets(self):
        # goes through each model, grabs its exchange met names, and bundles
        # them into a single list
        all_exchanged_mets = []
        for m in self.models:
            all_exchanged_mets.extend(m.get_exchange_metabolites())
        all_exchanged_mets = sorted(list(set(list(all_exchanged_mets))))
        self.all_exchanged_mets = all_exchanged_mets


    def add_model(self, model):
        """
        add a cometspy model to this layout
        
        This is the preferred way to add models to existing layout.

        Parameters
        ----------
        
        model : cometspy.model
            a cometspy.model with an initial_pop
            
        Examples
        --------
        
        >>> import cobra.test
        >>> import cometspy as c
        >>> layout = c.layout()
        >>> ecoli = cobra.test.create_test_model("ecoli")
        >>> salmonella = cobra.test.create_test_model("salmonella")
        >>> ecoli = c.model(ecoli)
        >>> salmonella = c.model(salmonella)
        >>> ecoli.open_exchanges()
        >>> salmonella.open_exchanges()
        >>> ecoli.initial_pop = [0, 0, 1.e-7]
        >>> salmonella.initial_pop = [0, 0, 1.e-7]
        >>> layout.add_model(ecoli)
        >>> layout.add_model(salmonella)
        

        """
        self.models.append(model)
        self.update_models()

    def __check_if_initial_pops_in_range(self):
        """
        checks if initial pops for each model are in the grid or raises error
        
        A common error is users setting model.initial_pop but not setting
        layout.grid.  Since they occur in different places, this method is 
        useful for a check before running.

        """
        for m in self.models:
            for founder in m.initial_pop:
                if (founder[0] < 0) or (founder[1] > (self.grid[1]-1)):
                    message = "the initial pop of a model is outside of layout.grid."
                    message += f" Either increase layout.grid or adjust {m.id}'s initial_pop"
                    raise ValueError(message)
                    
    def __get_met_number(self, met):
        """ returns the met number (of the external mets) given a name """
        met_number = [x for x in range(len(self.all_exchanged_mets)) if
                      self.all_exchanged_mets[x] == met][0]
        return(met_number)

def _read_file(filename):
    f = open(filename, 'r')
    f_lines = f.read()
    f.close()
    return f_lines
