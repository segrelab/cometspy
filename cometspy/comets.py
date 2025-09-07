
'''
The comets module runs COMETS simulations and stores output.

Generally, a comets object is created just before calling run(). Afterwards,
any saved data (e.g. total_biomass) can be accessed from the object.



'''

import re
import subprocess as sp
import pandas as pd
import os
import glob
import numpy as np
import platform

__author__ = "Djordje Bajic, Jean Vila, Jeremy Chacon, Ilija Dukovski"
__copyright__ = "Copyright 2024, The COMETS Consortium"
__credits__ = ["Djordje Bajic", "Jean Vila", "Jeremy Chacon", "Ilija Dukovski"]
__license__ = "MIT"
__version__ = "0.6.1"
__comets_compatibility__ = "2.12.4" # version of comets this was tested with (except signaling)
__comets_compatibility_signaling__ = "2.12.4" # version signaling was tested with

__maintainer__ = "Djordje Bajic"
__email__ = "djordje.bajic@yale.edu"
__status__ = "Beta"


def _readlines_file(filename):
    f = open(filename, 'r')
    f_lines = f.readlines()
    f.close()
    return f_lines


class comets:
    """
    the main simulation object to run COMETS

    The comets class stores a previously-made layout and params, and
    interacts with COMETS to run the simulation. It also stores any
    generated output, which could include total_biomass, biomass,
    media, or fluxes, as set in the params object. std_out from the COMETS
    simulation is saved in the attribute run_output and can be useful to
    examine to debug errors.

    When creating a comets object, the optional relative_dir path is useful
    when one is to run multiple simulations simultaneously, otherwise
    temporary files may overwrite each other.

    Parameters
    ----------

    layout : layout
        a cometspy.layout containing models and media information
    parameters : params
        a cometspy.params containing specified simulation parameters
    relative_dir : str, optional
        a directory to place temporary simulation files.

    Attributes
    ----------

    layout : cometspy.layout
        the layout containing cometspy.model[s] and media information
    params : cometspy.params
        the params containing simulation and default biological parameters
    working_dir : str
        the directory at which to save temporary sim files
    GUROBI_HOME : str
        the directory where GUROBI exists on the system
    COMETS_HOME : str
        the directory where COMETS exists on the system
    VERSION : str
        the version of comets, read from the files at COMETS_HOME
    classpath_pieces : dict
        classpath separated into library name (key) and location (value)
    JAVA_CLASSPATH : str
        a generated (overwritable) string containing the java classpath
    run_output : str
        generated object containing text from COMETS sim's std_out
    run_errors : str
        generated object containing text from COMETS sim's std_err
    total_biomass : pandas.DataFrame
        generated object containing total biomass from simulation
    biomass : pandas.DataFrame
        generated object containing spatially-explicit biomass from sim
    velocity : pandas.DataFrame
        generated object containing spatially-explicit velocity from sim
    specific_media : pandas.DataFrame
        generated object containing information on selected media from sim
    media : pandas.Dataframe
        generated object containing spatially-explicit media from sim
    fluxes_by_species : dict{model_id : pandas.DataFrame}
        generated object containing each species' spatial-explicit fluxes
    genotypes : pandas.DataFrame
        generated object containing genotypes if an evolution sim was run

    Examples
    --------

    >>> # assume layout and params have already been created
    >>> # and cometspy was imported as c
    >>> sim = c.comets(layout, params)
    >>> sim.run() # this could take from seconds to hours
    >>> print(sim.run_output) # to see the std_out
    >>> sim.total_biomass.plot(x = "cycle")
    >>> # assume params.all_params["writeBiomassLog"] == True, and that
    >>> # params.all_params["BiomassLogRate"] = 1000, and that
    >>> # params.all_params["maxCycles"] = 5000
    >>> im = sim.get_biomass_image(layout.models[0].id, 4000)
    >>> from matplotlib import pyplot as plt
    >>> plt.imshow(im / np.max(im))


    """

    def __init__(self, layout,
                 parameters, relative_dir : str =''):

        # define instance variables
        self.working_dir = os.getcwd() + '/' + relative_dir
        try:
            self.GUROBI_HOME = os.environ['GUROBI_COMETS_HOME']
            os.environ['GUROBI_HOME'] = self.GUROBI_HOME
        except:
            try:
                self.GUROBI_HOME = os.environ['GUROBI_HOME']
            except:
                try:
                    self.GUROBI_HOME = os.environ['COMETS_GUROBI_HOME']
                except:
                    self.GUROBI_HOME = ''
                    print("could not find environmental variable GUROBI_COMETS_HOME or GUROBI_HOME or COMETS_GUROBI_HOME")
                    print("COMETS will not work with GUROBI until this is solved. ")
                    print("Here is a solution:")
                    print("    1. import os and set os.environ['GUROBI_HOME'] then try to make a comets object again")
                    print("       e.g.   import os")
                    print("              os.environ['GUROBI_HOME'] = 'C:\\\\gurobi902\\\\win64'")
        self.COMETS_HOME = os.environ['COMETS_HOME']
        self.VERSION = os.path.splitext(os.listdir(os.environ['COMETS_HOME'] +
                                                   '/bin')[0])[0]

        # set default classpaths, which users may change
        self.__build_default_classpath_pieces()
        self.__build_and_set_classpath()
        self.__test_classpath_pieces()

        # check to see if user has the libraries where expected

        self.layout = layout
        self.parameters = parameters

        # dealing with output files
        self.parameters.set_param("useLogNameTimeStamp", False)
        self.parameters.set_param("TotalBiomassLogName",
                                  "totalbiomasslog" + '_' + hex(id(self)))
        self.parameters.set_param("BiomassLogName",
                                  "biomasslog" + '_' + hex(id(self)))
        self.parameters.set_param("FluxLogName",
                                  "fluxlog" + '_' + hex(id(self)))
        self.parameters.set_param("MediaLogName",
                                  "medialog" + '_' + hex(id(self)))
        self.parameters.set_param("velocityMultiConvLogName",
                                  "velocitymulticonvlog" + '_' + hex(id(self)))


    def __build_default_classpath_pieces(self):
        """
        sets up what it thinks the classpath should be
        """
        self.classpath_pieces = {}
        self.classpath_pieces['gurobi'] = (self.GUROBI_HOME +
                                           '/lib/gurobi.jar')
	
        self.classpath_pieces['or_tools_java'] = (self.COMETS_HOME +
                                           '/lib/or-tools/9.4.1874/' + 'ortools-java-9.4.1874.jar')

        self.classpath_pieces['or_tools_linux'] = (self.COMETS_HOME +
                                           '/lib/or-tools/9.4.1874/' + 'ortools-linux-x86-64-9.4.1874.jar')

        self.classpath_pieces['junit'] = glob.glob(self.COMETS_HOME +
                                                   '/lib/junit' + '/**/*junit*',
                                                   recursive=True)[0]
        self.classpath_pieces['hamcrest'] = glob.glob(self.COMETS_HOME +
                                                      '/lib' + '/**/*hamcrest*',
                                                      recursive=True)[0]

        self.classpath_pieces['jogl_all'] = glob.glob(self.COMETS_HOME +
                                                      '/lib' + '/**/jogl-all.jar',
                                                      recursive=True)[0]
        self.classpath_pieces['gluegen_rt'] = glob.glob(self.COMETS_HOME +
                                                      '/lib' + '/**/gluegen-rt.jar',
                                                      recursive=True)[0]
        self.classpath_pieces['gluegen'] = glob.glob(self.COMETS_HOME +
                                                      '/lib' + '/**/gluegen.jar',
                                                      recursive=True)[0]
        self.classpath_pieces['gluegen_rt_natives'] = glob.glob(
            self.COMETS_HOME + '/lib' + '/**/gluegen-rt-natives-linux-amd64.jar',
            recursive=True)[0]
        self.classpath_pieces['jogl_all_natives'] = glob.glob(
            self.COMETS_HOME + '/lib' + '/**/jogl-all-natives-linux-amd64.jar',
            recursive=True)[0]

        self.classpath_pieces['jmatio'] = glob.glob(self.COMETS_HOME +
                                                      '/lib' + '/**/jmatio.jar',
                                                      recursive=True)[0]
        self.classpath_pieces['jmat'] = glob.glob(self.COMETS_HOME +
                                                      '/lib' + '/**/jmatio.jar',
                                                      recursive=True)[0]

        self.classpath_pieces['concurrent'] = glob.glob(self.COMETS_HOME +
                                                      '/lib' + '/**/concurrent.jar',
                                                      recursive=True)[0]
        self.classpath_pieces['colt'] = glob.glob(self.COMETS_HOME +
                                                      '/lib' + '/**/colt.jar',
                                                      recursive=True)[0]

        __lang3 = glob.glob(self.COMETS_HOME +
                            '/lib' + '/**/commons-lang3*jar',
                            recursive=True)
        self.classpath_pieces['lang3'] = [i for i in __lang3
                                          if 'test' not in i
                                          and 'sources' not in i][0]

        __math3 = glob.glob(self.COMETS_HOME +
                            '/lib' + '/**/commons-math3*jar',
                            recursive=True)
        self.classpath_pieces['math3'] = [i for i in __math3
                                          if 'test' not in i
                                          and 'sources' not in i
                                          and 'tools' not in i
                                          and 'javadoc' not in i][0]

        self.classpath_pieces['jdistlib'] = glob.glob(self.COMETS_HOME +
                                                      '/lib' + '/**/*jdistlib*',
                                                      recursive=True)[0]
        self.classpath_pieces['bin'] = (self.COMETS_HOME +
                                        '/bin/' + self.VERSION + '.jar')


    def __build_and_set_classpath(self):
        ''' builds the JAVA_CLASSPATH from the pieces currently in
        self.classpath_pieces '''
        paths = list(self.classpath_pieces.values())
        if platform.system() == 'Windows':
            classpath = ';'.join(paths)
            classpath = '\"' + classpath + '\"'
            self.JAVA_LIB = '\"' + self.GUROBI_HOME + '/lib;' + self.GUROBI_HOME+ '/bin;'+ self.COMETS_HOME+ '/lib/jogl/jogamp-all-platforms/lib'+ '\"'

        else:
            classpath = ':'.join(paths)
        self.JAVA_CLASSPATH = classpath

    def __test_classpath_pieces(self):
        ''' checks to see if there is a file at each location in classpath
        pieces. If not, warns the user that comets will not work without the
        libraries. Tells the user to either edit those pieces (if in linux)
        or just set the classpath directly'''
        if platform.system() == 'Windows':
            return # Windows uses the script, so classpath doesn't matter as long as env variable set
        broken_pieces = self.__get_broken_classpath_pieces()
        if len(broken_pieces) == 0:
            pass  # yay! class files are where we hoped
        else:
            print('Warning: java class libraries cannot be found')
            print('These are the expected locations for dependencies:')

            print('Dependency \t\t\t expected path')
            print('__________ \t\t\t _____________')
            for key, value in broken_pieces.items():
                print('{}\t\t\t{}'.format(key, value))
            print('\n  You have two options to fix this problem:')
            print('1.  set each class path correctly by doing:')
            print('    comets.set_classpath(libraryname, path)')
            print('    e.g.   comets.set_classpath(\'hamcrest\', ' +
                  '\'/home/chaco001/comets/junit/hamcrest-core-1.3.' +
                  'jar\')\n')
            print('    note that versions dont always have to ' +
                  'exactly match, but you\'re on your own if they ' +
                  'don\'t\n')
            print('2.  fully define the classpath yourself by ' +
                  'overwriting comets.JAVA_CLASSPATH')
            print('       look at the current comets.JAVA_CLASSPATH ' +
                  'to see how this should look.')

    def __get_broken_classpath_pieces(self):
        ''' checks to see if there is a file at each location in classpath
        pieces. Saves the pieces where there is no file and returns them as a
        dictionary, where the key is the common name of the class library and
        the value is the path '''  #

        broken_pieces = {}         #
        for key, value in self.classpath_pieces.items():
            if not os.path.isfile(value):  #
                broken_pieces[key] = value
        return(broken_pieces)

    def set_classpath(self, libraryname : str, path : str):
        """
        sets the full path to a required specified java library

        This can be used to set non-default classpaths. Note that currently,
        it does not work for windows, because windows runs COMETS slightly
        differently than Unix.

        Parameters
        ----------

        libraryname : str
            the library for which the new path is being supplied.
        path : str
            the full path, including file name, to the library

        Examples
        --------

        >>> sim = c.comets(layout, params)
        >>> sim.set_classpath("jmatio", "/opt/jmatio/jmatio.jar")

        """
        self.classpath_pieces[libraryname] = path
        self.__build_and_set_classpath()

    def run(self, delete_files : bool = True, progress : bool = True):
        """
        run a COMETS simulation

        This runs the COMETS simulation specified by the layout and params
        objects supplied to this comets object. It creates the files needed
        for COMETS into the current directory or the optional one specified
        when this object was created.

        Once complete (or if an error has occurred), this tries to read
        simulation logs including data, as well as the std_out in the
        run_output attribute.

        If the optional delete_files is set to False, then temporary files and
        data log files are not deleted. They are deleted by default.

        Parameters
        ----------

        delete_files : bool, optional
            Whether to delete simulation and log files. The default is True.

        progress : bool, optional
            Whether to display a progress bar via tqdm. The default is True.

        Examples
        --------

        >>> sim = c.comets(layout, params)
        >>> sim.run(delete_files = True)
        >>> print(sim.run_output)
        >>> print(sim.total_biomass)

        """
        print('\nRunning COMETS simulation ...')
        #print('\nDebug Here ...')

        # If evolution is true, write the biomass but not the total biomass log
        if self.parameters.all_params['evolution']:
            self.parameters.all_params['writeTotalBiomassLog'] = False
            self.parameters.all_params['writeBiomassLog'] = True

        to_append = '_' + hex(id(self))
        # write the files for comets in working_dir
        c_global = self.working_dir + '.current_global' + to_append
        c_package = self.working_dir + '.current_package' + to_append
        c_script = self.working_dir + '.current_script' + to_append

        self.layout.write_necessary_files(self.working_dir, to_append)

        # self.layout.write_layout(self.working_dir + '.current_layout')
        self.parameters.write_params(c_global, c_package)

        if os.path.isfile(c_script):
            os.remove(c_script)
        with open(c_script, 'a') as f:
            f.write('load_comets_parameters ' + '.current_global' + to_append + '\n')
            f.writelines('load_package_parameters ' + '.current_package' + to_append + '\n')
            f.writelines('load_layout ' + '.current_layout' + to_append)

        if platform.system() == 'Windows':
            self.cmd = ('\"' + self.COMETS_HOME +
                        '\\comets_scr' + '\" \"' +
                        c_script +
                        '\"')
        else:
            # simulate
            self.cmd = ('java -classpath ' + self.JAVA_CLASSPATH +
                        # ' -Djava.library.path=' + self.D_JAVA_LIB_PATH +
                        ' edu.bu.segrelab.comets.Comets -loader' +
                        ' edu.bu.segrelab.comets.fba.FBACometsLoader' +
                        ' -script "' + c_script + '"')

        process = sp.Popen(self.cmd,
                     cwd = self.working_dir,
                     shell=True, stdout=sp.PIPE, stderr=sp.STDOUT)

        # progress bar
        if progress:
            from tqdm.auto import tqdm # auto-detects whether to use terminal progress bar or notebook-style one
            
            prog = tqdm(range(self.parameters.all_params["maxCycles"]),
                        desc = "Progress",
                        unit = "cycle")

            self.run_output = ""
            
            with process.stdout:
                for line in iter(process.stdout.readline, b''):
                    line = line.decode("ascii")
                    self.run_output += line
                    if(line.strip().startswith("Cycle ") & (line.strip() != "Cycle 1")):
                        prog.update(1)
            process.wait()

            self.run_errors = process.stderr
        
        else:
            self.run_output, self.run_errors = process.communicate()
            self.run_output = self.run_output.decode('ascii','ignore')

        if self.run_errors is not None:
            self.run_errors = self.run_errors.decode('ascii','ignore')
        else:
            self.run_errors = "STDERR empty."

        # Raise RuntimeError if simulation had nonzero exit
        self.__analyze_run_output()

        # '''----------- READ OUTPUT ---------------------------------------'''
        # Read total biomass output
        if self.parameters.all_params['writeTotalBiomassLog']:
            tbmf = _readlines_file(
                self.working_dir + self.parameters.all_params['TotalBiomassLogName'])
            tbmf = [x.replace(",",".") for x in tbmf] # for systems that use commas as decimal place
            self.total_biomass = pd.DataFrame([re.split(r'\t+', x.strip())
                                               for x in tbmf],
                                              columns=['cycle'] +
                                              self.layout.get_model_ids())
            self.total_biomass = self.total_biomass.astype('float')
            self.total_biomass.cycle = self.total_biomass.cycle.astype('int')
            if delete_files:
                os.remove(self.working_dir + self.parameters.all_params['TotalBiomassLogName'])

        # Read flux
        if self.parameters.all_params['writeFluxLog']:

            max_rows = 4 + max([len(m.reactions) for m in self.layout.models])

            self.fluxes = pd.read_csv(self.working_dir + self.parameters.all_params['FluxLogName'],
                                      sep="\\s+",
                                      header=None, names=range(max_rows))
            # deal with commas-as-decimals
            if any([isinstance(self.fluxes.iloc[0,i], str) for i in range(self.fluxes.shape[1])]):
                self.fluxes = pd.read_csv(self.working_dir + self.parameters.all_params['FluxLogName'],
                                      decimal = ",", sep="\\s+",
                                      header=None, names=range(max_rows))
            if delete_files:
                os.remove(self.working_dir + self.parameters.all_params['FluxLogName'])
            self.__build_readable_flux_object()

        # Read media logs
        if self.parameters.all_params['writeMediaLog']:
            self.media = pd.read_csv(self.working_dir + self.parameters.all_params[
                'MediaLogName'], sep="\\s+", names=('metabolite',
                                                               'cycle', 'x',
                                                               'y',
                                                               'conc_mmol'))
            # deal with commas-as-decimals
            if isinstance(self.media.loc[0, "conc_mmol"], str):
                self.media = pd.read_csv(self.working_dir + self.parameters.all_params[
                'MediaLogName'],
                    decimal = ",", sep="\\s+", names=('metabolite',
                                                               'cycle', 'x',
                                                               'y',
                                                               'conc_mmol'))
            if delete_files:
                os.remove(self.working_dir + self.parameters.all_params['MediaLogName'])

        # Read spatial biomass log
        if self.parameters.all_params['writeBiomassLog']:
            self.biomass = pd.read_csv(self.working_dir + self.parameters.all_params[
                'BiomassLogName'], header=None, delimiter=r'\s+', names=['cycle', 'x', 'y',
                                                                         'species', 'biomass'])
            # deal with commas-as-decimals
            if isinstance(self.biomass.loc[0,"biomass"], str):
                self.biomass = pd.read_csv(self.working_dir + self.parameters.all_params[
                'BiomassLogName'], header=None,
                    decimal = ",",
                    delimiter=r'\s+', names=['cycle', 'x', 'y','species', 'biomass'])

            # cut off extension added by toolbox
            self.biomass['species'] = [sp[:-4] if '.cmd' in sp else sp for sp in self.biomass.species]

            if delete_files:
                os.remove(self.working_dir + self.parameters.all_params['BiomassLogName'])

        # Read spatial velocity log
        if self.parameters.all_params['writeVelocityMultiConvLog']:
            self.velocity = pd.read_csv(self.working_dir + self.parameters.all_params[
                'velocityMultiConvLogName'], header=None, delimiter=r'\s+', names=['cycle', 'species', 'x', 'y',
                                                                            'velocityX','velocityY'])
            # deal with commas-as-decimals
            #if isinstance(self.velocity.loc[0,"velocity"], str):
            #    self.velocity = pd.read_csv(self.working_dir + self.parameters.all_params[
            #    'velocityMultiConvLogName'], header=None,
            #        decimal = ",",
            #        delimiter=r'\s+', names=['cycle', 'x', 'y','species', 'velocityX','velocityY'])

            # cut off extension added by toolbox
            self.velocity['species'] = [sp[:-4] if '.cmd' in sp else sp for sp in self.velocity.species]

            if delete_files:
                os.remove(self.working_dir + self.parameters.all_params['velocityMultiConvLogName'])

        # Read evolution-related logs
        if 'evolution' in list(self.parameters.all_params.keys()):
            if self.parameters.all_params['evolution']:
                genotypes_out_file = 'GENOTYPES_' + self.parameters.all_params[
                    'BiomassLogName']
                self.genotypes = pd.read_csv(genotypes_out_file,
                                             header=None, delimiter=r'\s+',
                                             names=['Ancestor',
                                                    'Mutation',
                                                    'Species'])
                if delete_files:
                    os.remove(self.working_dir + genotypes_out_file)

        # Read specific media output
        if self.parameters.all_params['writeSpecificMediaLog']:
            spec_med_file = self.working_dir + self.parameters.all_params['SpecificMediaLogName']
            self.specific_media = pd.read_csv(spec_med_file, delimiter=r'\s+')
            # deal with commas-as-decimals
            if any([isinstance(self.specific_media.iloc[0,i], str) for i in range(3, self.specific_media.shape[1])]):
                self.specific_media = pd.read_csv(spec_med_file, decimal = ",",delimiter=r'\s+')

            if delete_files:
                os.remove(self.working_dir + self.parameters.all_params['SpecificMediaLogName'])

        # clean workspace
        if delete_files:
            self.layout.delete_model_files(self.working_dir)
            os.remove(c_global)
            os.remove(c_package)
            os.remove(c_script)
            os.remove(self.working_dir + '.current_layout' + to_append)
            os.remove(self.working_dir + 'COMETS_manifest.txt')  # todo: stop writing this in java
        print('Done!')

    def __build_readable_flux_object(self):
        """ comets.fluxes is an odd beast, where the column position has a
        different meaning depending on what model the row is about. Therefore,
        this function creates separate dataframes, stored in a dictionary with
        model_id as a key, that are much more human-readable."""

        self.fluxes_by_species = {}
        for i in range(len(self.layout.models)):
            model_num = i + 1

            model_id = self.layout.models[model_num - 1].id
            model_rxn_names = list(self.layout.models[
                model_num - 1].reactions.REACTION_NAMES)
            model_rxn_len = len(model_rxn_names)

            sub_df = self.fluxes.loc[self.fluxes[3] == model_num]

            # this tosses extraneous columns and the model num column
            sub_df = sub_df.drop(sub_df.columns[model_rxn_len+4: len(sub_df.columns)],
                                 axis=1)
            sub_df = sub_df.drop(sub_df.columns[3], axis=1)
            sub_df.columns = ["cycle", "x", "y"] + model_rxn_names
            self.fluxes_by_species[model_id] = sub_df

    def __analyze_run_output(self):
        if "End of simulation" in self.run_output:
            return
        else:
            print("Error: COMETS simulation did not complete\n")
            print("     examine comets.run_output for the full java trace\n")
            print("     if we detect a common reason, it will be stated in the RuntimeError at the bottom")

            if "Could not find or load main class edu.bu.segrelab.comets.Comets" in self.run_output:
                message = "Could not find or load main class edu.bu.segrelab.comets.Comets\n"
                message += "check if comets.version and comets.classpath_pieces['bin'] \n"
                message += "point to an actual comets.jar file\n"
                message += "this problem may be associated with a malformed\n"
                message += "os.environ['COMETS_HOME'] environmental variable\n"
                message += "that can be overwritten by, for example, \n"
                message += ">>> import os\n"
                message += ">>> os.environ['COMETS_HOME'] = '/home/comets/'"
                raise RuntimeError(f"COMETS simulation did not complete:\n {message}")

            loc = self.run_output.find("NoClassDefFoundError")
            if loc != -1:
                error_string = self.run_output[(loc+22):(loc+100)]
                missing_class = error_string.split("\n")[0]
                if missing_class[0:6] == "gurobi":
                    message = "JAVA could not find gurobi.\n"
                    message += "try the following: \n"
                    message += ">>> import os\n"
                    message += ">>> os.environ['GUROBI_COMETS_HOME']\n"
                    message += "if there is nothing there try setting that variable\n"
                    message += "to the location of gurobi.jar, for example:\n"
                    message += ">>> os.environ['GUROBI_COMETS_HOME'] = '/opt/gurobi900/linux64'"
                else:
                    message = f"JAVA could not find a needed class: {missing_class}\n"
                    message += "make sure it is in your java classpath\n"
                    message += "this can be changed with comets.set_classpath()\n"
                    message += "if in Unix. In Windows, it suggests that something changed\n"
                    message += "with the dependencies installed alongside COMETS"
                raise RuntimeError(f"COMETS simulation did not complete:\n {message}")

            message = "undetected reason. examine comets.run_output for JAVA trace"
            raise RuntimeError(f"COMETS simulation did not complete:\n {message}")


    def get_metabolite_image(self, met : str, cycle : int) -> np.array:
        """
        returns an image of metabolite concentrations at a given cycle

        This will only work if media was saved at the current cycle. This
        requires the following parameters to have been set:

            params.set_param("writeMediaLog",True)
            params.set_param("MediaLogRate", n) # n is in cycles

        Notes
        -----

        There may be a bug where cycles are + 1 what they should be. We will
        fix this soon but for now be aware you may need to +1 your desired
        cycle. The bug does not affect anything within the simulation.

        Parameters
        ----------

        met : str
            the name of the metabolite
        cycle : int
            the cycle to get the metabolite data

        Returns
        -------

            A 2d numpy array which can be visualized like an image.

        Examples
        --------

        >>> sim.run()
        >>> # assume sim.params.all_params["MediaLogRate"] = 100 and that
        >>> # sim.params.all_params["maxCycles"] = 1000
        >>> im = sim.get_metabolite_image("ac_e", 500)
        >>> from matplotlib import pyplot as plt # may need to be installed
        >>> plt.imshow(im)

        """
        if not self.parameters.all_params['writeMediaLog']:
            raise ValueError("media log was not recorded during simulation")
        if met not in list(self.layout.media.metabolite):
            raise NameError("met " + met + " is not in layout.media.metabolite")
        if cycle not in list(np.unique(self.media['cycle'])):
            raise ValueError('media was not saved at the desired cycle. try another.')
        im = np.zeros((self.layout.grid[0], self.layout.grid[1]))
        aux = self.media.loc[np.logical_and(self.media['cycle'] == cycle,
                                            self.media['metabolite'] == met)]
        for index, row in aux.iterrows():
            im[int(row['x']-1), int(row['y']-1)] = row['conc_mmol']
        return(im)

    def get_biomass_image(self, model_id : str, cycle : int) -> np.array:
        """
        returns an image of biomass concentrations at a given cycle

        This will only work if biomass was saved at the current cycle. This
        requires the following parameters to have been set:
        >>> params.set_param("writeBiomassLog",True)
        >>> params.set_param("BiomassLogRate", n) # n is in cycles


        Parameters
        ----------

        model_id : str
            the id of the model to get biomass data on
        cycle : int
            the cycle to get the biomass data

        Returns
        -------

            A 2d numpy array which can be visualized like an image.

        Examples
        --------

        >>> sim.run()
        >>> # assume sim.params.all_params["BiomassLogRate"] = 100 and that
        >>> # sim.params.all_params["maxCycles"] = 1000
        >>> im = sim.get_biomass_image("iJO1366", 500) # e.g. the ecoli model
        >>> from matplotlib import pyplot as plt # may need to be installed
        >>> plt.imshow(im)

        """
        if not self.parameters.all_params['writeBiomassLog']:
            raise ValueError("biomass log was not recorded during simulation")
        if model_id not in list(np.unique(self.biomass['species'])):
            raise NameError("model " + model_id + " is not one of the model ids")
        if cycle not in list(np.unique(self.biomass['cycle'])):
            raise ValueError('biomass was not saved at the desired cycle. try another.')
        im = np.zeros((self.layout.grid[0], self.layout.grid[1]))
        aux = self.biomass.loc[np.logical_and(self.biomass['cycle'] == cycle,
                                    self.biomass['species'] == model_id), :]
        for index, row in aux.iterrows():
            im[int(row['x']-1), int(row['y']-1)] = row['biomass']
        return(im)

    def get_flux_image(self, model_id : str,
                       reaction_id : str, cycle : int) -> np.array:
        """
        returns a 2d numpy array showing fluxes at a given cycle

        This will only work if flux was saved at the current cycle. This
        requires the following parameters to have been set:

            params.set_param("writeFluxLog",True)
            params.set_param("FluxLogRate", n) # n is in cycles

        Parameters
        ----------

        model_id : str
            the id of the model about which to get fluxes
        reaction_id : str
            the id of the reaction about which to get fluxes
        cycle : int
            the cycle at which to get fluxes

        Returns
        -------
        a 2d numpy array which can be visualized like an image

        Examples
        --------

        >>> sim.run()
        >>> # assume sim.params.all_params["FluxLogRate"] = 100 and that
        >>> # sim.params.all_params["maxCycles"] = 1000
        >>> # assume a model used was iJO1366
        >>> im = sim.get_flux_image("iJO1366", "EX_ac_e", 500)
        >>> from matplotlib import pyplot as plt # may need to be installed
        >>> plt.imshow(im)

        """
        if not self.parameters.all_params['writeFluxLog']:
            raise ValueError("flux log was not recorded during simulation")
        if model_id not in [m.id for m in self.layout.models]:
            raise NameError("model " + model_id + " is not one of the model ids")
        im = np.zeros((self.layout.grid[0], self.layout.grid[1]))
        temp_fluxes = self.fluxes_by_species[model_id]
        if cycle not in list(np.unique(temp_fluxes['cycle'])):
            raise ValueError('flux was not saved at the desired cycle. try another.')
        if reaction_id not in list(temp_fluxes.columns):
            raise NameError("reaction_id " + reaction_id +
                            " is not a reaction in the desired model")
        aux = temp_fluxes.loc[temp_fluxes['cycle'] == cycle, :]
        for index, row in aux.iterrows():
            im[int(row['x']-1), int(row['y']-1)] = row[reaction_id]
        return(im)

    def get_metabolite_time_series(self, upper_threshold : float = 1000.) -> pd.DataFrame:
        """
        returns a pandas DataFrame containing extracellular metabolite time series

        Parameters
        ----------

        upper_threshold : float (optional)
            metabolites ever above this are not returned
        """
        total_media = self.media.groupby(by = ["metabolite", "cycle"]).agg(func = sum).reset_index().drop(columns = ["x", "y"])
        total_media = total_media.pivot(columns = "metabolite", values = "conc_mmol", index = ["cycle"]).reset_index().fillna(0.)
        exceeded_threshold = [x for x in total_media.min().index[total_media.min() > upper_threshold] if x != "cycle"]
        total_media = total_media.drop(columns = exceeded_threshold)
        return(total_media)

    def get_species_exchange_fluxes(self, model_id : str, threshold : float = 0.) -> pd.DataFrame:
        """
        returns a pandas DataFrame containing one model's exchange flux time series

        Parameters
        ----------

        model_id : str
            id of the model
        threshold : float (optional)
            abs(flux) must exceed this to be returned
        """
        fluxes = self.fluxes_by_species[model_id].copy()
        fluxes = fluxes.groupby(by = "cycle").agg(func = sum).reset_index().drop(columns = ["x", "y"])
        not_exch = [x for x in fluxes.columns if "EX_" not in x and x != "cycle"]
        fluxes = fluxes.drop(columns =not_exch)
        didnt_exceed_threshold = [x for x in np.abs(fluxes).max().index[np.abs(fluxes).max() < threshold] if x != "cycle"]
        fluxes = fluxes.drop(columns = didnt_exceed_threshold)
        return(fluxes)

# TODO: check for manual changes within layout that may not have triggered flags. See layout.write_layout for details
# TODO: fix read_comets_layout to always expect text addresses of comets model files
# TODO: remove comets manifest (preferably, dont write it)
# TODO: find quicker reading solution than the pd.read_csv stringIO hack
# TODO: solve weird rounding errors when reading from comets model
# TODO: give warning when unknown parameter is set
# TODO: write parameters in single file
# TODO: model biomass should be added in the layout "add_model" method, and not as a model class field
# TODO: adding models seems to remove media that has been previously set up
# TODO: write documentation properly for all functions and classes
