#!/usr/bin/env python

'''
The Comets module serves as a Python user interface to COMETS.
For more information see https://comets-manual.readthedocs.io/en/latest/
'''

import re
import math
import subprocess as sp
import pandas as pd
import os
import cobra
import io
import numpy as np

__author__ = "Djordje Bajic, Jean Vila"
__copyright__ = "Copyright 2019, The COMETS Consortium"
__credits__ = ["Djordje Bajic", "Jean Vila"]
__license__ = "MIT"
__version__ = "0.2.1"
__maintainer__ = "Djordje Bajic"
__email__ = "djordje.bajic@yale.edu"
__status__ = "Beta"


class CorruptLine(Exception):
    pass


class OutOfGrid(Exception):
    pass


class UnallocatedMetabolite(Exception):
    pass


def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


def read_file(filename):
    f = open(filename, 'r')
    f_lines = f.read()
    f.close()
    return f_lines


def readlines_file(filename):
    f = open(filename, 'r')
    f_lines = f.readlines()
    f.close()
    return f_lines


class model:
    def __init__(self, model = None):
        self.initial_pop = [[0, 0, 0.0]]
        self.id = None
        self.reactions = pd.DataFrame(columns=['REACTION_NAMES', 'ID',
                                          'LB', 'UB', 'EXCH', 'EXCH_IND','V_MAX', 'KM', 'HILL'])
        self.smat = pd.DataFrame(columns=['metabolite',
                                     'rxn',
                                     's_coef'])
        self.metabolites = pd.DataFrame(columns=['METABOLITE_NAMES'])
        
        self.vmax_flag = False
        self.km_flag = False
        self.hill_flag = False
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
            else: # assume it is a path
                if model[-3:] == "cmd":
                    self.read_comets_model(model)
                else:
                    self.read_cobra_model(model)

                    
                    
    def get_reaction_names(self):
        return(list(self.reactions['REACTION_NAMES']))
                    
    def get_exchange_metabolites(self):
        """ useful for layouts to grab these and get the set of them"""
        exchmets = pd.merge(self.reactions.loc[self.reactions['EXCH'], 'ID'], self.smat,
                            left_on='ID', right_on='rxn',
                            how='inner')['metabolite']
        exchmets = self.metabolites.iloc[exchmets-1]
        return(exchmets.METABOLITE_NAMES)
        
    def change_bounds(self, reaction, lower_bound, upper_bound):
        if not reaction in self.reactions['REACTION_NAMES'].values:
            print('reaction couldnt be found')
            return
        self.reactions.loc[self.reactions['REACTION_NAMES'] == reaction, 'LB'] = lower_bound
        self.reactions.loc[self.reactions['REACTION_NAMES'] == reaction, 'UB'] = upper_bound

    def get_bounds(self, reaction):
        if not reaction in self.reactions['REACTION_NAMES'].values:
            print('reaction couldnt be found')
            return
        lb = float(self.reactions.loc[self.reactions['REACTION_NAMES'] == reaction, 'LB'])
        ub = float(self.reactions.loc[self.reactions['REACTION_NAMES'] == reaction, 'UB'])
        return((lb,ub))
        
    def change_vmax(self, reaction, vmax):
        if not reaction in self.reactions['REACTION_NAMES'].values:
            print('reaction couldnt be found')
            return
        self.vmax_flag = True
        self.reactions.loc[self.reactions['REACTION_NAMES'] == reaction, 'V_MAX'] = vmax
        
    def change_km(self, reaction, km):
        if not reaction in self.reactions['REACTION_NAMES'].values:
            print('reaction couldnt be found')
            return
        self.km_flag = True
        self.reactions.loc[self.reactions['REACTION_NAMES'] == reaction, 'KM'] = km
        
    def change_hill(self, reaction, hill):
        if not reaction in self.reactions['REACTION_NAMES'].values:
            print('reaction couldnt be found')
            return
        self.hill_flag = True
        self.reactions.loc[self.reactions['REACTION_NAMES'] == reaction, 'HILL'] = hill
        
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
        # assign the dataframes we just built
        self.reactions = reactions
        self.metabolites = metabolites
        
        
    def write_comets_model(self, alternate_path = None):
        
        if alternate_path is not None:
            path_to_write = alternate_path
        else:
            path_to_write = self.id + '.cmd'
        
        # format variables for writing comets model
        bnd = self.reactions.loc[(self.reactions['LB'] != self.default_bounds[0]) |
                            (self.reactions['UB'] != self.default_bounds[1]),
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

            f.write('OBJECTIVE_STYLE\n' + self.obj_style + '\n')
            f.write(r'//' + '\n')

            f.write('OPTIMIZER ' + self.optimizer + '\n')
            f.write(r'//' + '\n')


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
        self.grid = [1,1]
        self.media = pd.DataFrame(columns=['metabolite',
                                           'init_amount',
                                           'diff_c',
                                           'g_static',
                                           'g_static_val',
                                           'g_refresh'])
        # local_media is a dictionary with locations as keys, and as values, 
        # another dictionary with metabolite names as keys and amounts as values
        # this information sets initial, location-specific media amounts.  
        self.local_media = {}
        self.global_diff = None
        self.local_refresh = []
        self.refresh = []
        self.local_refresh2 = {} # wip: do this like local_media? ask djb
        self.local_static = []
        self.local_static2 = {}
        self.initial_pop_type = "custom" # JMC not sure purpose of this 
        self.initial_pop = [] # 
        self.all_exchanged_mets = []
        
        self.default_diff_c = 5.0e-6 
        self.default_g_static = 0
        self.default_g_static_val = 0
        self.default_g_refresh = 0
        
        self.__local_media_flag = False
        self.__diffusion_flag = False
        self.__refresh_flag = False
        self.__static_flag = False
        
        if input_obj is None:
            print("building empty layout model\nmodels will need to be added with layout.add_model()")
        elif isinstance(input_obj, str):
            if not os.path.isfile(input_obj):
                raise IOError(' when running comets.layout(), input_obj is a string, and therefore should be a path to a layout; however, no file could be found at that path destionation')
            self.read_comets_layout(input_obj)
        else:
            if not isinstance(input_obj, list):
                input_obj = [input_obj] # probably just one cobra model 
            self.models = input_obj
            self.update_models()
                


    def read_comets_layout(self, input_obj):
                # .. load layout file
        f_lines = [s for s in read_file(input_obj).splitlines() if s]
        filedata_string = os.linesep.join(f_lines)
        end_blocks = []
        for i in range(0, len(f_lines)):
            if '//' in f_lines[i]:
                end_blocks.append(i)

        # '''----------- MODELS ----------------------------------------'''
        '''
        Models can be specified in layout as either comets format models
        or .xml format (sbml cobra compliant)
        '''            
        models = f_lines[0].split()
        if len(models) > 1:
            self.models = f_lines[0].split()[1:]
            self.update_models()
        else:
            print('Warning: No models in layout')
                
        # '''----------- GRID ------------------------------------------'''
        try:
            self.grid = [int(i) for i in f_lines[2].split()[1:]]
            if len(self.grid) < 2:
                raise CorruptLine
        except CorruptLine:
            print('\n ERROR CorruptLine: Only ' + str(len(self.grid)) +
                  ' dimension(s) specified for world grid')
            
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
            try:
                for i in range(lin_diff+1, lin_diff_end):
                    diff_spec = [float(x) for x in f_lines[i].split()]
                    if diff_spec[0] > len(self.media.metabolite)-1:
                        raise UnallocatedMetabolite
                    else:
                        self.media.loc[int(diff_spec[0]),
                                       'diff_c'] = diff_spec[1]
            except UnallocatedMetabolite:
                print('\n ERROR UnallocatedMetabolite: Some diffusion ' +
                      'values correspond to unallocated metabolites')

        # '''----------- MEDIA REFRESH----------------------------------'''
        # .. global refresh values
        self.__refresh_flag = False
        if 'REFRESH' in filedata_string:
            self.__refresh_flag = True
            lin_refr = re.split('refresh',
                                filedata_string)[0].count('\n')
            lin_refr_end = next(x for x in end_blocks if x > lin_refr)

            g_refresh = [float(x) for x in f_lines[lin_refr].split()[1:]]

            try:
                if len(g_refresh) != len(media_names):
                    raise CorruptLine
                else:
                    self.media['g_refresh'] = g_refresh
            except CorruptLine:
                print('\n ERROR CorruptLine: Number of global refresh ' +
                      'values does not match number of \nmedia ' +
                      'metabolites in provided layout file')

            # .. local refresh values
            lin_refr += 1
            try:
                for i in range(lin_refr, lin_refr_end):
                    refr_spec = [float(x) for x in f_lines[i].split()]
                    if len(refr_spec) != len(self.media.metabolite)+2:
                        raise CorruptLine
                    elif (refr_spec[0] >= self.grid[0] or
                          refr_spec[1] >= self.grid[1]):
                        raise OutOfGrid
                    else:
                        self.local_refresh.append(refr_spec)

            except CorruptLine:
                print('\n ERROR CorruptLine: Some local "refresh" lines ' +
                      'have a wrong number of entries')
            except OutOfGrid:
                print('\n ERROR OutOfGrid: Some local "refresh" lines ' +
                      'have coordinates that fall outside of the ' +
                      '\ndefined ' + 'grid')

        # '''----------- STATIC MEDIA ----------------------------------'''
        # .. global static values
        lin_static = re.split('static',
                              filedata_string)[0].count('\n')
        lin_stat_end = next(x for x in end_blocks if x > lin_static)

        g_static = [float(x) for x in f_lines[lin_static].split()[1:]]
        try:
            if len(g_static) != 2*len(self.media.metabolite):
                raise CorruptLine
            else:
                self.media.loc[:, 'g_static'] = [int(x)
                                                 for x in g_static[0::2]]
                self.media.loc[:, 'g_static_val'] = [float(x) for x in
                                                     g_static[1::2]]
        except CorruptLine:
            print('\nERROR CorruptLine: Wrong number of global ' +
                  'static values')
            
        # .. local static values
        lin_static += 1
        try:
            for i in range(lin_static, lin_stat_end):
                stat_spec = [float(x) for x in f_lines[i].split()]
                if len(stat_spec) != (2*len(self.media.metabolite))+2:
                    raise CorruptLine
                elif (stat_spec[0] >= self.grid[0] or
                      stat_spec[1] >= self.grid[1]):
                    raise OutOfGrid
                else:
                    self.local_static.append(stat_spec)
                    
        except CorruptLine:
            print('\n ERROR CorruptLine: Wrong number of local static ' +
                  'values at some lines')
        except OutOfGrid:
            print('\n ERROR OutOfGrid: Some local "static" lines have ' +
                  ' coordinates that fall outside of the defined grid')

        # '''----------- INITIAL POPULATION ----------------------------'''
        lin_initpop = re.split('initial_pop',
                               filedata_string)[0].count('\n')
        lin_initpop_end = next(x for x in end_blocks if x > lin_initpop)

        g_initpop = f_lines[lin_initpop].split()[1:]
        
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
            try:
                for i in range(lin_initpop, lin_initpop_end):
                    ipop_spec = [float(x) for x in
                                 f_lines[lin_initpop].split()]
                    if len(ipop_spec) != len(self.models)+2:
                        raise CorruptLine
                    if (ipop_spec[0] >= self.grid[0] or
                            ipop_spec[1] >= self.grid[1]):
                        raise OutOfGrid
                    else:
                        self.initial_pop.append(ipop_spec)
            except CorruptLine:
                print('Problem at some initial population lines')
            except OutOfGrid:
                print('Some initial population values' +
                      ' fall outside of the defined grid')
                
          
        
    def get_model_ids(self):
        ids = [x.id for x in self.models]
        return(ids)
        
    def write_necessary_files(self, layoutfile):
        self.write_layout(layoutfile)
        self.write_model_files()
        
    def write_model_files(self):
        '''writes each model file'''
        for model in self.models:
            model.write_comets_model()
            
    def display_current_media(self):
        print(self.media[self.media['init_amount'] != 0.0])
        
        
    def set_specific_metabolite(self, met, amount):
        try:
            self.media.loc[self.media['metabolite'] == met, 'init_amount'] = amount
        except:
            print("the specified metabolite " + met + "is not able to be taken up, not added to media")             
        
    def set_specific_metabolite_at_location(self, met, location, amount):
        """ allows the user to specify a metabolite going to a specific location
        in a specific amount.  useful for generating non-homogenous environments.
        THe met should be the met name (e.g. 'o2_e') the location should be a 
        tuple (e.g. (0, 5)), and the amount should be a float / number"""
        if met not in self.all_exchanged_mets:
            raise Exception('met is not in the list of exchangeable mets')
        self.__local_media_flag = True
        if location not in list(self.local_media.keys()):
            self.local_media[location] = {}
        self.local_media[location][met] = amount
        
    def set_specific_refresh(self, met, amount):
        try:
            self.media.loc[self.media['metabolite'] == met, 'g_refresh'] = amount
            self.__refresh_flag = True
        except:
            print("the specified metabolite " + met + "is not able to be taken up, not added to media")
      
        
    def set_specific_refresh_at_location(self, met, location, amount):
        if met not in self.all_exchanged_mets:
            raise Exception('met is not in the list of exchangeable mets')
        self.__refresh_flag = True
        if location not in list(self.local_refresh2.keys()):
            self.local_refresh2[location] = {}
        self.local_refresh2[location][met] = amount
        
    def set_specific_static(self, met, amount):
        try:
            self.media.loc[self.media['metabolite'] == met, 'g_static'] = 1
            self.media.loc[self.media['metabolite'] == met, 'g_static_val'] = amount
            self.__static_flag = True
        except:
            print("the specified metabolite " + met + "is not able to be taken up, not added to media")

    def set_specific_static_at_location(self, met, location, amount):
        if met not in self.all_exchanged_mets:
            raise Exception('met is not in the list of exchangeable mets')
        self.__static_flag = True
        if location not in list(self.local_static2.keys()):
            self.local_static2[location] = {}
        self.local_static2[location][met] = amount
            
            
    def add_typical_trace_metabolites(self, amount = 1000.0):
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
            self.media.loc[self.media['metabolite'] == met, 'init_amount'] = amount
            
    def write_layout(self, outfile):
        ''' Write the layout in a file'''

        if os.path.isfile(outfile):
            os.remove(outfile)
        
        lyt = open(outfile, 'a')
        self.__write_models_and_world_grid_chunk(lyt)
        self.__write_media_chunk(lyt)
        self.__write_diffusion_chunk(lyt)
        self.__write_local_media_chunk(lyt)
        #self.__write_refresh_chunk(lyt)
        self.__write_refresh2_chunk(lyt)
        #self.__write_static_chunk(lyt)
        self.__write_static2_chunk(lyt)
        lyt.write(r'  //' + '\n')

        self.__write_initial_pop_chunk(lyt)

        lyt.close()
        
    def __write_models_and_world_grid_chunk(self, lyt):
        """ writes the top 3 lines  to the open lyt file"""
        lyt.write('model_file ' +
                  '.cmd '.join(self.get_model_ids()) +
                  '.cmd\n')
        lyt.write('  model_world\n')
        
        lyt.write('    grid_size ' +
                  ' '.join([str(x) for x in self.grid]) + '\n')        
        
    def __write_media_chunk(self, lyt):
        """ used by write_layout to write the global media information to the open lyt file """
        lyt.write('    world_media\n')
        for i in range(0, len(self.media)):
            lyt.write('      ' + self.media.metabolite[i] +
                      ' ' + str(self.media.init_amount[i]) + '\n')
        lyt.write(r'    //' + '\n')
        
    def __write_local_media_chunk(self, lyt):
        """ used by write_layout to write the location-specific initial metabolite data"""
        if self.__local_media_flag:    
            lyt.write('    media\n')
            locs = list(self.local_media.keys())
            for loc in locs:
                # this chunk goes in order, not by name, so must get met number
                # for each location, make a list with zeros for each met.  put 
                # non-zero numbers where the self.local_media tells us to
                met_amounts_in_order = [0] * len(self.all_exchanged_mets)
                for met in list(self.local_media[loc].keys()):
                    met_amounts_in_order[self.__get_met_number(met)] = self.local_media[loc][met]
                lyt.write('      ')
                lyt.write('{} {} '.format(loc[0], loc[1]))
                lyt.write(' '.join(str(x) for x in met_amounts_in_order))
                lyt.write('\n')
            lyt.write('    //\n')
            
                
    def __write_refresh2_chunk(self, lyt):
        if self.__refresh_flag:
            lyt.write('    media_refresh ' +
                      ' '.join([str(x) for x in self.media.
                                g_refresh.tolist()]) +
                      '\n')
            locs = list(self.local_refresh2.keys())
            if len(locs) > 0:
                for loc in locs:
                    met_amounts_in_order = [0] * len(self.all_exchanged_mets)
                    for met in list(self.local_refresh2[loc].keys()):
                        met_amounts_in_order[self.__get_met_number(met)] = self.local_refresh2[loc][met]
                    met_amounts_in_order.insert(0, loc[1])
                    met_amounts_in_order.insert(0, loc[0])
                    lyt.write('      ' +
                          ' '.join([str(x) for x in met_amounts_in_order]) +
                          '\n')
            lyt.write(r'    //' + '\n') 
            
    def __write_static2_chunk(self, lyt):
        if self.__static_flag:
            g_static_line = [None]*(len(self.media)*2)
            g_static_line[::2] = self.media.g_static
            g_static_line[1::2] = self.media.g_static_val
            lyt.write('    static_media ' +
                      ' '.join([str(x) for x in g_static_line]) + '\n')
            locs = list(self.local_static2.keys())
            if len(locs) > 0:
                for loc in locs:
                    # this is 2 * len because there is a pair of values for each met
                    # the first value is a flag--0 if not static, 1 if static
                    # the second value is the amount if it is static
                    met_amounts_in_order = [0] * 2 * len(self.all_exchanged_mets)
                    for met in list(self.local_static2[loc].keys()):
                        met_amounts_in_order[self.__get_met_number(met) * 2] = 1 # the flag
                        met_amounts_in_order[self.__get_met_number(met) * 2 + 1 ] = self.local_static2[loc][met]
                    met_amounts_in_order.insert(0, loc[1])
                    met_amounts_in_order.insert(0, loc[0])
                    lyt.write('      ' +
                          ' '.join([str(x) for x in met_amounts_in_order]) +
                          '\n')
            lyt.write(r'    //' + '\n') 

                    
    def __write_diffusion_chunk(self, lyt):
        """ used by write_layout to write the metab-specific diffusion data to the open lyt file """
        if self.__diffusion_flag:
            lyt.write('    diffusion_constants ' +
                      str(self.global_diff) +
                      '\n')
            for i in range(0, len(self.media)):
                if not math.isnan(self.media.diff_c[i]):
                    lyt.write('      ' + str(i) + ' ' +
                              str(self.media.diff_c[i]) + '\n')
            lyt.write(r'    //' + '\n')        
            
    def __write_refresh_chunk(self, lyt):
        """ used by write_layout to write the global media refresh information to the open lyt file """
        if self.__refresh_flag:
            lyt.write('    media_refresh ' +
                      ' '.join([str(x) for x in self.media.
                                g_refresh.tolist()]) +
                      '\n')
            for i in range(0, len(self.local_refresh)):
                lyt.write('      ' +
                          ' '.join([str(x) for x in self.local_refresh[i]]) +
                          '\n')
            lyt.write(r'    //' + '\n')   
            
    def __write_static_chunk(self, lyt):
        """ used by write_layout to write the global static media information to the open lyt file """        
        g_static_line = [None]*(len(self.media)*2)
        g_static_line[::2] = self.media.g_static
        g_static_line[1::2] = self.media.g_static_val
        lyt.write('    static_media ' +
                  ' '.join([str(x) for x in g_static_line]) + '\n')
        
        for i in range(0, len(self.local_static)):
            lyt.write('      ' +
                      ' '.join([str(x) for x in self.local_static[i]]) +
                      '\n')
        lyt.write(r'    //' + '\n')
        lyt.write(r'  //' + '\n')
        
    def __write_initial_pop_chunk(self, lyt):
        """ writes the initial pop to the open lyt file and adds the closing //s """
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
        for i, model in enumerate(self.models):
            if not isinstance(model.initial_pop[0], list): # in case this wasnt a nested list
                model.initial_pop = [model.initial_pop]
            for pop in model.initial_pop:
                curr_line = [0] * (n_models + 2)
                curr_line[0] = pop[0]
                curr_line[1] = pop[1]
                curr_line[i+2] = pop[2]       
                initial_pop.append(curr_line)
        self.initial_pop = initial_pop
        
    def add_new_mets_to_media(self):
        ## usually run right after build_exchange mets, to add any new mets
        ## to the media data.frame
        
        for met in self.all_exchanged_mets:
            if met not in self.media['metabolite'].values:
                new_row = pd.DataFrame.from_dict({'metabolite': [met],
                                     'init_amount': [0],
                                     'diff_c': [self.default_diff_c],
                                     'g_static': [self.default_g_static],
                                     'g_static_val': [self.default_g_static_val],
                                     'g_refresh': [self.default_g_refresh]})
                self.media = pd.concat([self.media, new_row], 
                                       ignore_index = True)
            
        
    def build_exchanged_mets(self):
        # goes through each model, grabs its exchange met names, and bundles
        # them into a single list
        all_exchanged_mets = []
        for model in self.models:
            all_exchanged_mets.extend(model.get_exchange_metabolites())
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
        met_number = [x for x in range(len(self.all_exchanged_mets)) if self.all_exchanged_mets[x] == met][0]
        return(met_number)
    
        
class params:
    '''
    Class storing COMETS parameters
    '''
    def __init__(self, global_params=None, package_params=None):
        self.all_params = {'BiomassLogName': 'biomass.txt',
                           'BiomassLogRate': 1,
                           'FluxLogName': 'flux_out',
                           'FluxLogRate': 5,
                           'MediaLogName': 'media_out',
                           'MediaLogRate': 5,
                           'TotalBiomassLogName': 'total_biomass_out.txt',
                           'maxCycles': 100,
                           'saveslideshow': False,
                           'totalBiomassLogRate': 1,
                           'useLogNameTimeStamp': False,
                           'writeBiomassLog': False,
                           'writeFluxLog': False,
                           'writeMediaLog': False,
                           'writeTotalBiomassLog': True,
                           'batchDilution': False,
                           'dilFactor': 10,
                           'dilTime': 2,
                           'cellSize': 1e-13,
                           'allowCellOverlap': True,
                           'deathRate': 0,
                           'defaultHill': 1,
                           'defaultKm': 0.01,
                           'defaultVmax': 10,
                           'defaultAlpha': 1,
                           'defaultW': 10,
                           'defaultDiffConst': 1e-5,
                           'exchangestyle': 'Monod Style',
                           'flowDiffRate': 3e-9,
                           'growthDiffRate': 0,
                           'maxSpaceBiomass': 0.1,
                           'minSpaceBiomass': 0.25e-10,
                           'numDiffPerStep': 10,
                           'numRunThreads': 1,
                           'showCycleCount': True,
                           'showCycleTime': False,
                           'spaceWidth': 0.02,
                           'timeStep': 0.1,
                           'toroidalWorld': False,
                           'simulateActivation': False,
                           'activateRate': 0.001,
                           'randomSeed': 0,
                           'colorRelative': True,
                           'slideshowColorRelative': True,
                           'slideshowRate': 1,
                           'slideshowLayer': 0,
                           'slideshowExt': 'png',
                            # 'biomassMotionStyle': 'Diffusion' +
                            # '2D(Crank-Nicolson)', TODO: this not working
                           'numExRxnSubsteps': 5,
                           'costlyGenome': True,
                           'geneFractionalCost': 1e-4,
                           'evolution': False,
                           'mutRate': 1e-5,
                           'addRate': 1e-5}
        self.all_params = dict(sorted(self.all_params.items(),
                                      key=lambda x: x[0]))
        
        self.all_type = {'BiomassLogName': 'global',
                         'BiomassLogRate': 'global',
                         'FluxLogName': 'global',
                         'FluxLogRate': 'global',
                         'MediaLogName': 'global',
                         'MediaLogRate': 'global',
                         'TotalBiomassLogName': 'global',
                         'maxCycles': 'package',
                         'saveslideshow': 'global',
                         'totalBiomassLogRate': 'global',
                         'useLogNameTimeStamp': 'global',
                         'writeBiomassLog': 'global',
                         'writeFluxLog': 'global',
                         'writeMediaLog': 'global',
                         'writeTotalBiomassLog': 'global',
                         'batchDilution': 'global',
                         'dilFactor': 'global',
                         'dilTime': 'global',
                         'cellSize': 'global',
                         'allowCellOverlap': 'package',
                         'deathRate': 'package',
                         'defaultHill': 'package',
                         'defaultKm': 'package',
                         'defaultVmax': 'package',
                         'defaultW': 'package',
                         'defaultAlpha': 'package',
                         'defaultDiffConst': 'package',
                         'exchangestyle': 'package',
                         'flowDiffRate': 'package',
                         'growthDiffRate': 'package',
                         'maxSpaceBiomass': 'package',
                         'minSpaceBiomass': 'package',
                         'numDiffPerStep': 'package',
                         'numRunThreads': 'package',
                         'showCycleCount': 'package',
                         'showCycleTime': 'package',
                         'spaceWidth': 'package',
                         'timeStep': 'package',
                         'toroidalWorld': 'package',
                         'simulateActivation': 'global',
                         'activateRate': 'global',
                         'randomSeed': 'global',
                         'colorRelative': 'global',
                         'slideshowColorRelative': 'global',
                         'slideshowRate': 'global',
                         'slideshowLayer': 'global',
                         'slideshowExt': 'global',
                         'biomassMotionStyle': 'package',
                         'numExRxnSubsteps': 'package',
                         'costlyGenome': 'global',
                         'geneFractionalCost': 'global',
                         'evolution': 'package',
                         'mutRate': 'package',
                         'addRate': 'package'}
        self.all_type = dict(sorted(self.all_type.items(),
                                    key=lambda x: x[0]))

        # .. parse parameters files to python type variables
        if global_params is not None:
            with open(global_params) as f:
                for line in f:
                    if '=' in line:
                        k, v = line.split(' = ')
                        if v.strip() == 'true':
                            self.all_params[k.strip()] = True
                        elif v.strip() == 'false':
                            self.all_params[k.strip()] = False
                        elif v.strip().isdigit():
                            self.all_params[k.strip()] = int(v.strip())
                        elif isfloat(v.strip()):
                            self.all_params[k.strip()] = float(v.strip())
                        else:
                            self.all_params[k.strip()] = v.strip()

        if package_params is not None:
            with open(package_params) as f:
                for line in f:
                    if '=' in line:
                        k, v = line.split(' = ')
                        if v.strip() == 'true':
                            self.all_params[k.strip()] = True
                        elif v.strip() == 'false':
                            self.all_params[k.strip()] = False
                        elif v.strip().isdigit():
                            self.all_params[k.strip()] = int(v.strip())
                        elif isfloat(v.strip()):
                            self.all_params[k.strip()] = float(v.strip())
                        else:
                            self.all_params[k.strip()] = v.strip()

        # Additional processing.
        # If evolution is true, we dont want to write the total biomass log
        if self.all_params['evolution']:
            self.all_params['writeTotalBiomassLog'] = False

    ''' write parameters files; method probably only used by class comets'''
    def write_params(self, out_glb, out_pkg):

        if os.path.isfile(out_glb):
            os.remove(out_glb)

        if os.path.isfile(out_pkg):
            os.remove(out_pkg)

        # convert booleans to java format before writing
        towrite_params = {}
        for k, v in self.all_params.items():
            if v is True:
                towrite_params[k] = 'true'
            elif v is False:
                towrite_params[k] = 'false'
            else:
                towrite_params[k] = str(v)

        with open(out_glb, 'a') as glb, open(out_pkg, 'a') as pkg:
            for k, v in towrite_params.items():
                if self.all_type[k] == 'global':
                    glb.writelines(k + ' = ' + v + '\n')
                else:
                    pkg.writelines(k + ' = ' + v + '\n')


class comets:
    '''
    This class sets up an environment with all necessary for
    a comets simulation to run, runs the simulation, and stores the output
    data from it.
    '''
    def __init__(self, layout, parameters, working_dir=''):
        
        # define instance variables
        self.working_dir = os.getcwd() + '/' + working_dir
        self.GUROBI_HOME = os.environ['GUROBI_HOME']
        self.COMETS_HOME = os.environ['COMETS_HOME']
        
        self.VERSION = 'comets_evo'

        # set default classpaths, which users may change
        self.build_default_classpath_pieces()
        self.build_and_set_classpath()
        self.test_classpath_pieces()
        
        # check to see if user has the libraries where expected

        self.layout = layout
        self.parameters = parameters
        
        # dealing with output files
        self.parameters.all_params['useLogNameTimeStamp'] = False
        self.parameters.all_params['TotalBiomassLogName'] = (
            'total_biomass_log_' + hex(id(self)))
        self.parameters.all_params['BiomassLogName'] = (
            'biomass_log_' + hex(id(self)))
        self.parameters.all_params['FluxLogName'] = (
            'flux_log_' + hex(id(self)))
        self.parameters.all_params['MediaLogName'] = (
            'media_log_' + hex(id(self)))
        
    def build_default_classpath_pieces(self):
        self.classpath_pieces = {}
        self.classpath_pieces['gurobi'] = (self.GUROBI_HOME +
                                           '/lib/gurobi.jar')
        self.classpath_pieces['junit'] = (self.COMETS_HOME +
                                          '/lib/junit/junit-4.12.jar')
        self.classpath_pieces['hamcrest'] = (self.COMETS_HOME +
                                             '/junit/hamcrest-core-1.3.jar')
        self.classpath_pieces['jogl_all'] = (self.COMETS_HOME +
                                             '/lib/jogl/jogamp-all-' +
                                             'platforms/jar/jogl-all.jar')
        self.classpath_pieces['gluegen_rt'] = (self.COMETS_HOME +
                                               '/lib/jogl/jogamp-all-' +
                                               'platforms/jar/gluegen-rt.jar')
        self.classpath_pieces['gluegen'] = (self.COMETS_HOME +
                                            '/lib/jogl/jogamp-all-' +
                                            'platforms/jar/gluegen.jar')
        self.classpath_pieces['gluegen_rt_natives'] = (self.COMETS_HOME +
                                                       '/lib/jogl/jogamp-' +
                                                       'all-platforms/jar/' +
                                                       'gluegen-rt-natives-' +
                                                       'linux-amd64.jar')
        self.classpath_pieces['jogl_all_natives'] = (self.COMETS_HOME +
                                                     '/lib/jogl/' +
                                                     'jogamp-all-platforms/' +
                                                     'jar/jogl-all-natives-' +
                                                     'linux-amd64.jar')
        self.classpath_pieces['jmatio'] = (self.COMETS_HOME +
                                           '/lib/JMatIO/lib/jamtio.jar')
        self.classpath_pieces['jmat'] = (self.COMETS_HOME +
                                         '/lib/JMatIO/JMatIO-041212/' +
                                         'lib/jmat.jar')
        self.classpath_pieces['concurrent'] = (self.COMETS_HOME +
                                               '/lib/colt/lib/concurrent.jar')
        self.classpath_pieces['colt'] = (self.COMETS_HOME +
                                         '/lib/colt/lib/colt.jar')
        self.classpath_pieces['lang3'] = (self.COMETS_HOME +
                                          '/lib/commons-lang3-3.7/' +
                                          'commons-lang3-3.7.jar')
        self.classpath_pieces['bin'] = (self.COMETS_HOME +
                                        '/bin/' + self.VERSION + '.jar')
    
    def build_and_set_classpath(self):
        ''' builds the JAVA_CLASSPATH from the pieces currently in
        self.classpath_pieces '''
        paths = list(self.classpath_pieces.values())
        classpath = ':'.join(paths)
        self.JAVA_CLASSPATH = classpath
    
    def test_classpath_pieces(self):
        ''' checks to see if there is a file at each location in classpath
        pieces. If not, warns the user that comets will not work without the
        libraries. Tells the user to either edit those pieces (if in linux)
        or just set the classpath directly'''
        broken_pieces = self.get_broken_classpath_pieces()
        if len(broken_pieces) == 0:
            pass  # yay! class files are where we hoped
        print('warning:  we cannot find required java class libraries at the' +
              'expected locations')
        print('    specifically, we cannot find the following libraries at ' +
              'these locations:\n')
        print('library common name \t expected path')
        print('___________________ \t _____________')
        for key, value in broken_pieces.items():
            print('{}\t{}'.format(key, value))
        print('\n  You have two options to fix this problem:')
        print('1.  set each class path correctly by doing:')
        print('       comets.set_classpath(libraryname, path)')
        print('       e.g.   comets.set_classpath(\'hamcrest\', \'/home/' +
              'chaco001/comets/junit/hamcrest-core-1.3.jar\')\n')
        print('       note that versions dont always have to exactly match, ' +
              'but you\'re on your own if they don\'t\n')
        print('2.  fully define the classpath yourself by overwriting ' +
              'comets.JAVA_CLASSPATH')
        print('       look at the current comets.JAVA_CLASSPATH to see how ' +
              'this should look.')
                
    def get_broken_classpath_pieces(self):
        ''' checks to see if there is a file at each location in classpath
        pieces. Saves the pieces where there is no file and returns them as a
        dictionary, where the key is the common name of the class library and
        the value is the path '''
        broken_pieces = {}
        for key, value in self.classpath_pieces.items():
            if not os.path.isfile(value):
                broken_pieces[key] = value    
        return(broken_pieces)
        
    def set_classpath(self, libraryname, path):
        ''' tells comets where to find required java libraries
        e.g. comets.set_classpath(\'hamcrest\', \'/home/chaco001/
        comets/junit/hamcrest-core-1.3.jar\')
        Then re-builds the path'''
        self.classpath_pieces[libraryname] = path
        self.build_and_set_classpath()

    def run(self, delete_files = True):
        print('\nRunning COMETS simulation ...')
        # write the files for comets in working_dir
        c_global = self.working_dir + '.current_global'
        c_package = self.working_dir + '.current_package'
        c_script = self.working_dir + '.current_script'

        self.layout.write_necessary_files(self.working_dir + '.current_layout')
        #self.layout.write_layout(self.working_dir + '.current_layout')
        self.parameters.write_params(c_global, c_package)

        if os.path.isfile(c_script):
            os.remove(c_script)
        with open(c_script, 'a') as f:
            f.write('load_comets_parameters ' + c_global + '\n')
            f.writelines('load_package_parameters ' + c_package + '\n')
            f.writelines('load_layout ' + self.working_dir +
                         '.current_layout')
            
        # simulate
        self.cmd = ('java -classpath ' + self.JAVA_CLASSPATH +
                    # ' -Djava.library.path=' + self.D_JAVA_LIB_PATH +
                    ' edu.bu.segrelab.comets.Comets -loader' +
                    ' edu.bu.segrelab.comets.fba.FBACometsLoader' +
                    ' -script ' + c_script)
        
        p = sp.Popen(self.cmd, shell=True, stdout=sp.PIPE, stderr=sp.STDOUT)

        self.run_output, self.run_errors = p.communicate()
        self.run_output = self.run_output.decode()

        if self.run_errors is not None:
            self.run_errors = self.run_errors.decode()
        else:
            self.run_errors = "STDERR empty."
            
        # clean workspace
        if delete_files:
            os.remove(c_global)
            os.remove(c_package)
            os.remove(c_script)
            os.remove('.current_layout')
            os.remove('COMETS_manifest.txt')  # todo: stop writing this in java
        
        # '''----------- READ OUTPUT ---------------------------------------'''

        # Read total biomass output
        if self.parameters.all_params['writeTotalBiomassLog']:
            tbmf = readlines_file(
                self.parameters.all_params['TotalBiomassLogName'])
            self.total_biomass = pd.DataFrame([re.split(r'\t+', x.strip())
                                               for x in tbmf],
                                              columns=['cycle'] +
                                              self.layout.get_model_ids())
            self.total_biomass = self.total_biomass.astype('float')
            if delete_files:
                os.remove(self.parameters.all_params['TotalBiomassLogName'])
            
        # Read flux
        if self.parameters.all_params['writeFluxLog']:
            self.fluxes = []
            tff = readlines_file(
                self.parameters.all_params['FluxLogName'])
            for i in tff:
                self.fluxes.append(re.findall(r'\d+',
                                              re.search(r'\{(.*)\}',
                                                        i).group(1)) +
                                   [float(x)
                                    for x in re.search(r'\[(.*)\]',
                                                       i).group(1).split()])
            if delete_files:
                os.remove(self.parameters.all_params['FluxLogName'])
                
        # Read spatial biomass log
        if self.parameters.all_params['writeBiomassLog']:
            biomass_out_file = 'biomass_log_' + hex(id(self))
            self.biomass = pd.read_csv(biomass_out_file,
                                       header=None, delimiter=r'\s+',
                                       names=['Cycle', 'x', 'y',
                                              'species', 'biomass'])
            if delete_files:
                os.remove(biomass_out_file)
            
        # Read evolution-related logs
        if self.parameters.all_params['evolution']:
            evo_out_file = 'biomass_log_' + hex(id(self))
            self.evolution = pd.read_csv(evo_out_file,
                                         header=None, delimiter=r'\s+',
                                         names=['Cycle', 'x', 'y',
                                                'species', 'biomass'])
            genotypes_out_file = 'GENOTYPES_biomass_log_' + hex(id(self))
            self.genotypes = pd.read_csv(genotypes_out_file,
                                         header=None, delimiter=r'\s+',
                                         names=['Ancestor',
                                                'Mutation',
                                                'Species'])
            
        print('Done!')

# TODO: fix read_comets_layout to always expect text addresses of comets model files
# TODO: make sure layout loading uses the new formats for location-specific media, refresh, etc
# TODO: read media logs (after fixing format in java)
# TODO: read spatial biomass logs
# TODO: remove comets manifest (preferably, dont write it)
# TODO: find quicker reading solution than the pd.read_csv stringIO hack
# TODO: fucntions to generate predefined media, spatial layouts etc
# TODO: write noncustom initial pop in layout
# TODO: add barriers in layout class
# TODO: add units when printing params
# TODO: solve weird rounding errors when reading from comets model
# TODO: include all params in one file (maybe layout?) to avoid file writing
# TODO: update media with all exchangeable metabolites from all models
# TODO: give warning when unknown parameter is set

