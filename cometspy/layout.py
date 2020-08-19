'''
The layout module handles COMETS simulation layouts, including
media and spatial arrangement.
For more information see https://segrelab.github.io/comets-manual/
'''

import pandas as pd
import os
import numpy as np
import re
import math

from cometspy.model import model


def read_file(filename):
    f = open(filename, 'r')
    f_lines = f.read()
    f.close()
    return f_lines


class layout:
    '''
    Generates a COMETS layout either by reading from a file or by building one
    from a list of COBRA models. Or, with no arguments, build an empty layout.

    To read a layout from a file, give the path as a string:

        layout = comets.layout("./path/to/layout/layoutfile.txt")

    To build a layout from a list of models, give the models in a list:
        ijo = cobra.test.load

    '''
    def __init__(self, input_obj=None):

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
        self.global_diff = None
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

        self.__local_media_flag = False
        self.__diffusion_flag = False
        self.__refresh_flag = False
        self.__static_flag = False
        self.__barrier_flag = False
        self.__region_flag = False
        self.__ext_rxns_flag = False
        self.__periodic_media_flag = False

        if input_obj is None:
            print('building empty layout model\nmodels will need to be added' +
                  ' with layout.add_model()')
        elif isinstance(input_obj, str):
            if not os.path.isfile(input_obj):
                raise IOError(' when running comets.layout(), input_obj' +
                              ' is a string, and therefore should be a path' +
                              ' to a layout; however, no file could be found' +
                              ' at that path destionation')
            self.read_comets_layout(input_obj)
        else:
            if not isinstance(input_obj, list):
                input_obj = [input_obj]  # probably just one cobra model
            self.models = input_obj
            self.update_models()

    def set_region_parameters(self, region, diffusion, friction):
        """
        COMETS can have different regions with different substrate diffusivities
        and frictions.  Here, you set those parameters. For example, if a layout
        had three different substrates, and you wanted to define their diffusion
        for region 1, you would use:

            layout.set_region_parameters(1, [1e-6, 1e-6, 1e-6], 1.0)

        This does not affect a simulation unless a region map is also set, using
        the layout.set_region_map() function.
        """
        if not self.__region_flag:
            print("Warning: You are setting region parameters but a region" +
                  "map has not been set. Use layout.set_region_map() or these" +
                  "parameters will be unused")
        self.region_parameters[region] = [diffusion, friction]

    def set_region_map(self, region_map):
        """
        COMETS can have different regions with different substrate diffusivities
        and frictions.  Here, you set the map defining the regions. Specifically,
        you provide either:
            1) a numpy array whose shape == layout.grid, or
            2) a list of lists whose first length is grid[0] and second len is grid[1]

        Populating these objects should be integer values, beginning at 1 and
        incrementing only, that define the different grid areas.  These are
        intimately connected to region_parameters, which are set with
        layout.set_region_parameters()
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
                              rxnName, metabolites, stoichiometry, **kwargs):

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
                                  metabolite, function,
                                  amplitude, period, phase, offset):

        if (metabolite not in self.media['metabolite'].values):
            raise ValueError('the metabolite is not present in the media')
        if (function not in ['step', 'sin', 'cos', 'half_sin', 'half_cos']):
            raise ValueError(function + ': function unknown')

        self.periodic_media.append([self.media.index[self.media['metabolite']
                                                     == metabolite][0],
                                    function, amplitude, period,
                                    phase, offset])
        self.__periodic_media_flag = True

    def read_comets_layout(self, input_obj):

        # .. load layout file
        f_lines = [s for s in read_file(input_obj).splitlines() if s]
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

            self.global_diff = float(re.findall(r'\S+', f_lines[lin_diff].
                                                strip())[1])
            
            for i in range(lin_diff+1, lin_diff_end):
                diff_spec = [float(x) for x in f_lines[i].split()]
                if diff_spec[0] > len(self.media.metabolite)-1:
                    print('\n Warning: Corrupt line ' + str(i) + ' in diffusion' +
                          'values. \nLine not written. Check your layout file.')
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

    def get_model_ids(self):
        ids = [x.id for x in self.models]
        return(ids)

    def write_necessary_files(self, working_dir):
        self.write_layout(working_dir)
        self.write_model_files(working_dir)

    def write_model_files(self, working_dir=""):
        '''writes each model file'''
        for m in self.models:
            m.write_comets_model(working_dir)

    def display_current_media(self):
        print(self.media[self.media['init_amount'] != 0.0])

    def add_barriers(self, barriers):
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

    def set_specific_metabolite(self, met, amount, static=False):
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
            print('Warning: The added metabolite (' + met + ') is not' +
                  'able to be taken up by any of the current models')

    def set_specific_metabolite_at_location(self, met, location, amount):
        """ allows the user to specify a metabolite going to a specific location
        in a specific amount.  useful for generating non-homogenous
        environments. The met should be the met name (e.g. 'o2_e') the
        location should be a tuple (e.g. (0, 5)), and the amount should be
        a float / number"""
        if met not in self.all_exchanged_mets:
            raise Exception('met is not in the list of exchangeable mets')
        self.__local_media_flag = True
        if location not in list(self.local_media.keys()):
            self.local_media[location] = {}
        self.local_media[location][met] = amount

    def set_specific_refresh(self, met, amount):
        if met in self.media['metabolite'].values:
            self.media.loc[self.media['metabolite'] == met,
                           'g_refresh'] = amount
            self.__refresh_flag = True
        else:
            print("the specified metabolite " + met +
                  "is not in the medium; add it first!")

    def set_specific_refresh_at_location(self, met, location, amount):
        if met not in self.all_exchanged_mets:
            raise Exception('met is not in the list of exchangeable mets')
        self.__refresh_flag = True
        if location not in list(self.local_refresh.keys()):
            self.local_refresh[location] = {}
        self.local_refresh[location][met] = amount

    def set_specific_static(self, met, amount):
        if met in self.media['metabolite'].values:
            self.media.loc[self.media['metabolite'] == met,
                           'g_static'] = 1
            self.media.loc[self.media['metabolite'] == met,
                           'g_static_val'] = amount
            self.__static_flag = True
        else:
            print("the specified metabolite " + met +
                  "is not in the medium; add it first!")

    def set_specific_static_at_location(self, met, location, amount):
        if met not in self.all_exchanged_mets:
            raise Exception('met is not in the list of exchangeable mets')
        self.__static_flag = True
        if location not in list(self.local_static.keys()):
            self.local_static[location] = {}
        self.local_static[location][met] = amount

    def add_typical_trace_metabolites(self, amount=1000.0, static=True):
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
                # print('Warning: The added metabolite (' + met + ') is not' +
                #      'able to be taken up by any of the current models')
        self.media = self.media.reset_index(drop=True)

    def write_layout(self, working_dir):
        ''' Write the layout in a file'''
        outfile = working_dir + ".current_layout"
        if os.path.isfile(outfile):
            os.remove(outfile)

        lyt = open(outfile, 'a')
        self.__write_models_and_world_grid_chunk(lyt, working_dir)
        self.__write_media_chunk(lyt)
        self.__write_diffusion_chunk(lyt)
        self.__write_local_media_chunk(lyt)
        self.__write_refresh_chunk(lyt)
        self.__write_static_chunk(lyt)
        self.__write_barrier_chunk(lyt)
        self.__write_regions_chunk(lyt)
        self.__write_periodic_media_chunk(lyt)
        lyt.write(r'  //' + '\n')

        self.__write_initial_pop_chunk(lyt)
        self.__write_ext_rxns_chunk(lyt)
        lyt.close()

    def __write_models_and_world_grid_chunk(self, lyt, working_dir):
        """ writes the top 3 lines  to the open lyt file"""

        model_file_line = "{}.cmd".format(".cmd ".join(self.get_model_ids())).split(" ")
        model_file_line = "".join(["./" + _ + " " for _ in model_file_line])
        #model_file_line = working_dir + working_dir.join(model_file_line)
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
                      str(self.global_diff) +
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

    def update_models(self):
        self.build_initial_pop()
        self.build_exchanged_mets()
        self.add_new_mets_to_media()

    def build_initial_pop(self):
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
        self._resolve_initial_pop_line_dups()
        
    def _resolve_initial_pop_line_dups(self):
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

    def add_new_mets_to_media(self):
        # usually run right after build_exchange mets, to add any new mets
        # to the media data.frame

        for met in self.all_exchanged_mets:
            if met not in self.media['metabolite'].values:
                new_row = pd.DataFrame.from_dict({'metabolite': [met],
                                                  'init_amount': [0],
                                                  'diff_c': [self.default_diff_c],
                                                  'g_static': [self.default_g_static],
                                                  'g_static_val': [
                                                      self.default_g_static_val],
                                                  'g_refresh': [self.default_g_refresh]})
                self.media = pd.concat([self.media, new_row],
                                       ignore_index=True, sort=True)

    def build_exchanged_mets(self):
        # goes through each model, grabs its exchange met names, and bundles
        # them into a single list
        all_exchanged_mets = []
        for m in self.models:
            all_exchanged_mets.extend(m.get_exchange_metabolites())
        all_exchanged_mets = sorted(list(set(list(all_exchanged_mets))))
        self.all_exchanged_mets = all_exchanged_mets

    def update_media(self):
        # TODO: update media with all exchangeable metabolites from all models
        pass

    def add_model(self, model):
        self.models.append(model)
        self.update_models()

    def __get_met_number(self, met):
        """ returns the met number (of the external mets) given a name """
        met_number = [x for x in range(len(self.all_exchanged_mets)) if
                      self.all_exchanged_mets[x] == met][0]
        return(met_number)

