import pandas as pd
import os


def __isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


class params:
    '''
    object storing simulation-related and some default biological parameters
    
    The params object is an essential object needed to start a comets
    simulation, along with a layout object containing models. The params 
    object stores simulation information such as timeStep and maxCycles, as 
    well as many default biological parameters such as defaultKm and 
    defaultDiffC (the default metabolite diffusion constant). 
    
    All of the parameters are stored in a dictionary called 'all_params' and
    can either be adjusted directly or using set_param(). 
    
    When params are written to a COMETS object for COMETS to process, two
    files are generated ('.current_global' and 'current_package')
    
    Parameters
    ----------
    global_params : str, optional
        optional path to an existing 'global_params' file to load
    package_params : str, optional
        optional path to an existing 'package_params' file to load
        
    Examples
    --------
    
    >>> # make a params object and change two params
    >>> import cometspy as c
    >>> params = c.params()
    >>> params.set_param("timeStep", 0.01) # hours
    >>> params.set_param("spaceWidth", 0.5) # cm
    >>> # view all the params
    >>> params.show_params()
    >>> # access the params directly
    >>> params.all_params # this is a dict
    >>> # use a params object in a comets simulation
    >>> import cobra.test
    >>> model = c.model(cobra.test.create_test_model("textbook"))
    >>> model.initial_pop = [0, 0, 1.e-7]
    >>> model.open_exchanges()
    >>> layout = c.layout([model])
    >>> # normally we'd now alter the media in the layout
    >>> layout.add_typical_trace_metabolites()
    >>> layout.set_specific_metabolite('lcts_e', 0.05) # mmol
    >>> sim = c.comets(layout, params)
    >>> sim.run()
    
    Attributes
    ----------
    
    all_params : dict
        contains every possible param as keys, with their value as the value
    params_units : dict
        dictionary containing the units for each possible parameter

    '''
    def show_params(self):
        """
        utility function to show all possible parameters and their units
        """
        
        pdparams = pd.DataFrame({'VALUE': pd.Series(self.all_params),
                                 'UNITS': pd.Series(self.param_units)})
        return pdparams

    def set_param(self, name : str, value):
        """
        alters a specific parameter value
        
        See show_params() for a list of possible parameters to adjust. 
        
        Parameters
        ----------
        name : str
            the name of a specific parameter to change 
        value : variable
            the type of the value depends on the parameter.

        Examples
        --------
        >>> p = cometspy.params()
        >>> p.set_param("writeBiomassLog", True)
        >>> p.set_param("maxCycles", 100)
        >>> p.set_param("BiomassLogRate",5)
        >>> p.set_param("defaultVmax", 10.24)

        """
        if name in self.all_params:
            self.all_params[name] = value
        elif name.lower() in [n.lower() for n in self.all_params.keys()]:
            match = [n for n in self.all_params.keys() if n.lower() == name.lower()][0]
            self.all_params[match] = value
            print('Warning: inputted parameter ' + name + ' is properly spelled as ' + match)
        else:
            print('Warning: inputted parameter ' + name + ' does not exist')

    def get_param(self, name : str) -> object:
        """
        returns the value associated with the given parameter name

        Parameters
        ----------
        name : str
            the name of the parameter whose value is desired

        Returns
        -------
        object
            the type of object depends on the parameter requested

        """
        if name in self.all_params:
            return self.all_params[name]
        else:
            print('Parameter ' + name + ' does not exist')
           
    def __init__(self, global_params : str =None, 
                 package_params : str =None):
        self.all_params = {'writeSpecificMediaLog': False,
                           'specificMediaLogRate': 1,
                           'specificMedia': 'ac_e',
                           'SpecificMediaLogName': 'specific_media',
                           'BiomassLogName': 'biomass',
                           'BiomassLogRate': 1,
                           'biomassLogFormat': 'COMETS',
                           'FluxLogName': 'flux_out',
                           'FluxLogRate': 5,
                           'fluxLogFormat': 'COMETS',
                           'MediaLogName': 'media_out',
                           'MediaLogRate': 5,
                           'mediaLogFormat': 'COMETS',
                           'velocityMultiConvLogName': 'velocity_out',
                           'velocityMultiConvLogRate': 5,
                           'velocityMultiConvLogFormat': 'COMETS',
                           'TotalBiomassLogName': 'total_biomass_out',
                           'maxCycles': 100,
                           'saveslideshow': False,
                           'totalBiomassLogRate': 1,
                           'useLogNameTimeStamp': False,
                           'writeBiomassLog': False,
                           'writeFluxLog': False,
                           'writeMediaLog': False,
                           'writeVelocityMultiConvLog': False,
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
                            'velocityMultiConvLogName': '',
                            'velocityMultiConvLogRate': 'cycles',
                            'velocityMultiConvLogFormat': '',
                            'TotalBiomassLogName': '',
                            'maxCycles': 'cycles',
                            'saveslideshow': 'logical',
                            'totalBiomassLogRate': 'cycles',
                            'useLogNameTimeStamp': 'logical',
                            'writeBiomassLog': 'logical',
                            'writeFluxLog': 'logical',
                            'writeMediaLog': 'logical',
                            'writeVelocityMultiConvLog': 'logical',
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
                         'velocityMultiConvLogName': 'global',
                         'velocityMultiConvLogRate': 'global',
                         'velocityMultiConvLogFormat': 'global',
                         'TotalBiomassLogName': 'global',
                         'maxCycles': 'package',
                         'saveslideshow': 'global',
                         'totalBiomassLogRate': 'global',
                         'useLogNameTimeStamp': 'global',
                         'writeBiomassLog': 'global',
                         'writeFluxLog': 'global',
                         'writeMediaLog': 'global',
                         'writeVelocityMultiConvLog': 'global',
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
                        elif __isfloat(v.strip()):
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
                        elif __isfloat(v.strip()):
                            self.all_params[k.strip()] = float(v.strip())
                        else:
                            self.all_params[k.strip()] = v.strip()
                
        # Additional processing.
        # If evolution is true, we dont want to write the total biomass log
        if self.all_params['evolution']:
            self.all_params['writeTotalBiomassLog'] = False
            self.all_params['writeBiomassLog'] = True
                    
    def write_params(self, out_glb : str, out_pkg : str):
        """
        writes params data to the specified files
        
        Usually this is only used internally, though a user could use it to 
        pre-generate params files.

        Parameters
        ----------
        out_glb : str
            the path and name of the 'global' parameters file to be written
        out_pkg : str
            the path and name of the 'package' parameters file to be written


        """
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

