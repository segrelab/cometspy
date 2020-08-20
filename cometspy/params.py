
'''
The params module handles COMETS simulation parameters.
For more information see https://segrelab.github.io/comets-manual/
'''

import pandas as pd
import os


def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


class params:
    '''
    Class storing COMETS parameters
    '''
    def show_params(self):
        pdparams = pd.DataFrame({'VALUE': pd.Series(self.all_params),
                                 'UNITS': pd.Series(self.param_units)})
        return pdparams

    def set_param(self, name, value):
        if name in self.all_params:
            self.all_params[name] = value
        else:
            print('Parameter ' + name + ' does not exist')

    def get_param(self, name):
        if name in self.all_params:
            return self.all_params[name]
        else:
            print('Parameter ' + name + ' does not exist')
           
    def __init__(self, global_params=None, package_params=None):
        self.all_params = {'writeSpecificMediaLog': False,
                           'specificMediaLogRate': 1,
                           'specificMedia': 'ac_e',
                           'SpecificMediaLogName': 'specific_media.txt',
                           'BiomassLogName': 'biomass.txt',
                           'BiomassLogRate': 1,
                           'biomassLogFormat': 'COMETS',
                           'FluxLogName': 'flux_out',
                           'FluxLogRate': 5,
                           'fluxLogFormat': 'COMETS',
                           'MediaLogName': 'media_out',
                           'MediaLogRate': 5,
                           'mediaLogFormat': 'COMETS',
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
                           'biomassMotionStyle': 'Diffusion 2D(Crank-Nicolson)',
                           'numExRxnSubsteps': 5,
                           'costlyGenome': False,
                           'geneFractionalCost': 1e-4,
                           'evolution': False,
                           'mutRate': 0,
                           'addRate': 0,
                           'metaboliteDilutionRate': 0.}
        self.all_params = dict(sorted(self.all_params.items(),
                                      key=lambda x: x[0]))

        self.param_units = {'writeSpecificMediaLog': "logical",
                            'specificMediaLogRate': "cycles",
                            'specificMedia': "comma separated string",
                            'SpecificMediaLogName': '',
                            'BiomassLogName': '',
                            'BiomassLogRate': 'cycles',
                            'biomassLogFormat': '',
                            'FluxLogName': '',
                            'FluxLogRate': 'cycles',
                            'fluxLogFormat': '',
                            'MediaLogName': '',
                            'MediaLogRate': 'cycles',
                            'mediaLogFormat': '',
                            'TotalBiomassLogName': '',
                            'maxCycles': 'cycles',
                            'saveslideshow': 'logical',
                            'totalBiomassLogRate': 'cycles',
                            'useLogNameTimeStamp': 'logical',
                            'writeBiomassLog': 'logical',
                            'writeFluxLog': 'logical',
                            'writeMediaLog': 'logical',
                            'writeTotalBiomassLog': 'logical',
                            'batchDilution': 'logical',
                            'dilFactor': '',
                            'dilTime': 'hr.',
                            'cellSize': 'gr. cell dry weight',
                            'allowCellOverlap': 'logical',
                            'deathRate': 'percent/cycle',
                            'defaultHill': 'unitless',
                            'defaultKm': 'M (molar conc.)',
                            'defaultVmax': 'mmol /gr. CDW /hr',
                            'defaultAlpha': 'slope',
                            'defaultW': 'M (molar conc.)',
                            'defaultDiffConst': 'cm2/s',
                            'exchangestyle': 'uptake type',
                            'flowDiffRate': 'cm2/s',
                            'growthDiffRate': 'cm2/s',
                            'maxSpaceBiomass': 'gr. cell dry weight',
                            'minSpaceBiomass': 'gr. cell dry weight',
                            'numDiffPerStep': '',
                            'numRunThreads': '',
                            'showCycleCount': 'logical',
                            'showCycleTime': 'logical',
                            'spaceWidth': 'cm',
                            'timeStep': 'hr.',
                            'toroidalWorld': 'logical',
                            'simulateActivation': 'logical',
                            'activateRate': '',
                            'randomSeed': 'numeric',
                            'colorRelative': 'logical',
                            'slideshowColorRelative': 'logical',
                            'slideshowRate': 'cycles',
                            'slideshowLayer': 'numeric',
                            'slideshowExt': 'extension',
                            'biomassMotionStyle': 'type of motion',
                            'numExRxnSubsteps': 'numeric',
                            'costlyGenome': 'logical',
                            'geneFractionalCost': 'percent growth rate^(1/2) /reaction',
                            'evolution': 'logical',
                            'mutRate': 'deletions /reaction /generation',
                            'addRate': 'additions /generation',
                            'metaboliteDilutionRate': 'percent/cycle'}
        self.param_units = dict(sorted(self.param_units.items(),
                                       key=lambda x: x[0]))

        self.all_type = {'writeSpecificMediaLog': 'global',
                         'specificMediaLogRate': 'global',
                         'specificMedia': 'global',
                         'SpecificMediaLogName': 'global',
                         'BiomassLogName': 'global',
                         'BiomassLogRate': 'global',
                         'biomassLogFormat': 'global',
                         'FluxLogName': 'global',
                         'FluxLogRate': 'global',
                         'fluxLogFormat': 'global',
                         'MediaLogName': 'global',
                         'MediaLogRate': 'global',
                         'mediaLogFormat': 'global',
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
                         'addRate': 'package',
                         'metaboliteDilutionRate': 'package'}
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
            self.all_params['writeBiomassLog'] = True
                    
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

