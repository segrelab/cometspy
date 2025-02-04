# Installation 


The COMETS Python toolbox (cometspy) is available from the package manager PyPI cometspy using the pip command. To install, run: 

```
pip install cometspy 
```

Warning! The current version will work only with Pandas version 1.5.0 or higer. 
Pandas changed the variable line_terminator to lineterminator in the 'csv' modules.
More at: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.to_csv.html
"Changed in version 1.5.0: Previously was line_terminator, changed for consistency with read_csv and the standard library ‘csv’ module."

If you get an error message about line_terminator, update your pandas package. 

Of course, to use cometspy you first need to install COMETS on your computer. See section below. 

## Installation of COMETS 
There are two ways to install COMETS: using the COMETS installer or by unpacking the .tar.gz file. The easiest is to use the installer, especially recommended for individual use on a laptop or desktop. The installer can be downloaded from: https://comets.bu.edu The users are required to register, after which they can obtain the installer appropriate for their system. The installer guides the user through a standard GUI installation procedure that includes accepting the license agreement, choosing the directory where COMETS will be installed (the default directory is recommended), the option to create a desktop shortcut etc. The installer is available for the Windows, MacOS and Linux systems. 

In addition to the GUI installer, we provide the comets_x.x.x.tar.gz file for custom installation, typically on a Linux system. The file should be unpacked in the directory were COMETS will be installed with:

```
$tar -xzvf comets_x.x.x.tar.gz   ./
```

This will create the comets installation directory.

In Unix systems (Linux or MacOS) the user also needs to specify the COMETS_HOME environment variable, which has to point to the comets installation folder. In Linux systems, this is done by adding the following line to the .bashrc file located in the home folder: 
```
export COMETS_HOME = "/your/comets/installation/folder" 
```

<!-- ### Installation of the COMETS MATLAB Toolbox -->
<!-- The COMETS toolbox for MATLAB can be downloaded from https://github.com/segrelab/comets-toolbox. You may download the toolbox as an archive from the GitHub repository, or execute the following command from the command line (the folder ./comets-toolbox will be created in the working directory): -->

<!-- ``` -->
<!-- git clone https://github.com/segrelab/comets-toolbox.git comets-toolbox -->
<!-- ``` -->

<!-- A prerequisite to install the toolbox this way is to have installed git which can be found here: https://git-scm.com/. Once this folder has been created, run the following commands in MATLAB to add the toolbox and its subfolders to the MATLAB path: -->

<!-- ```matlab -->
<!-- >> addpath(genpath("comets-toolbox"),"-end"); -->
<!-- >> savepath(); -->
<!-- ``` -->
<!-- where comets-toolbox is the full path to the directory where the toolbox was installed. On a Windows system, for example, this path may be:  -->

<!-- ``` -->
<!-- C:\Users\username\comets-toolbox -->
<!-- ``` -->

<!-- where username is replaced with the specific one for the user. -->
   
<!-- In addition, this toolbox requires the installation of the COBRA toolbox, available at https://opencobra.github.io/ 40. Many functions of the COMETS toolbox will not work before loading the COBRA toolbox using the initCobraToolbox() command. The detailed instructions for installing the COBRA toolbox can be found here: https://opencobra.github.io/cobratoolbox/stable/installation.html. A prerequisite to install the COBRA toolbox is to have installed git which can be found here: https://git-scm.com/. Once git is installed, the toolbox can be installed by running:  -->

<!-- ``` -->
<!-- git clone --depth=1 https://github.com/opencobra/cobratoolbox.git cobratoolbox -->
<!-- ``` -->

<!-- Once this folder has been created, run the following commands in MATLAB to add the toolbox and its subfolders to the MATLAB path: -->

<!-- ```matlab -->
<!-- >> addpath(genpath("cobratoolbox"),"-end"); -->
<!-- >> savepath(); -->
<!-- ``` -->

<!-- where cobratoolbox is the full path to the directory where the toolbox was installed. On a Windows system, for example, this path may be:  -->

<!-- ``` -->
<!-- C:\Users\username\cobratoolbox -->
<!-- ``` -->

<!-- where username is replaced with the one specific for the user. -->
