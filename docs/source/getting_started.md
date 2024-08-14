# Getting started

The Comets python module is intended to offer a programatic, easy and intuitive interface to `COMETS`. While it internally uses the same `COMETS` Java engine as always, it replaces the legacy parameter files, simulation layouts or models, bash scripts and output files, by python objects in a single environment, which users deal with to perform simulations and analyze the results. 

Any comets simulation starts from a **layout** and a set of **parameters**. The layout specifies the environment (media metabolites, refresh values, periodic dilutions) and the species present in it, that is, the **models**. The parameters specify many simulation characteristics, such as number of iterations, timestep, type of metabolite exchange or whether to record different output logs. Therefore, the two main types of objects are `layout` and `params`. These two are passed to the `comets` class, which perform simulations and contain their output. 

In this section, we will walk through the basic functionalities of `COMETS` using the Python Toolbox, and more specific examples of usage will be provided in the next sections. 

## Preparing a model for `COMETS`
The class `model` is used to store the genome-scale metabolic models used in `COMETS` simulations. Most frequently, we will first load a model using `COBRAPy`. Then, we can pass it to the `COMETS` `model` class, which allows us to change `COMETS`-specific model parameters, such as initial population sizes.  


```python
import cobra
import cobra.test
import cometspy as c

# Load a textbook example model using the COBRAPy toolbox 
test_model = cobra.test.create_test_model('textbook')

# Use the above model to create a COMETS model
test_model = c.model(test_model)

# Change comets specific parameters, e.g. the initial biomass of the model
# Notre 
test_model.initial_pop = [0, 0, 1e-7] 
```

    Using license file /home/djordje/gurobi.lic
    Academic license - for non-commercial use only


## Setting `COMETS` simulation parameters
`COMETS` simulation parameters are stored in the `params` class, which contains just a `dict` object with the parameter names and values. If we initialize the class without arguments, it will contain the default parameter values (see [here]()). Once loaded, the parameter values can be visualized and modified as desired. 


```python
# Create a parameters object with default values 
my_params = c.params()

# Change the value of a parameter, for example number of simulation cycles
my_params.set_param('maxCycles', 100)

# Set some writeTotalBiomassLog parameter to True, in order to save the output
my_params.set_param('writeTotalBiomassLog', True)

# See avaliable parameters and their values
my_params.show_params()
```


<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>VALUE</th>
      <th>UNITS</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>BiomassLogName</th>
      <td>biomass.txt</td>
      <td></td>
    </tr>
    <tr>
      <th>BiomassLogRate</th>
      <td>1</td>
      <td>cycles</td>
    </tr>
    <tr>
      <th>FluxLogName</th>
      <td>flux_out</td>
      <td></td>
    </tr>
    <tr>
      <th>FluxLogRate</th>
      <td>5</td>
      <td>cycles</td>
    </tr>
    <tr>
      <th>MediaLogName</th>
      <td>media_out</td>
      <td></td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>writeBiomassLog</th>
      <td>False</td>
      <td>logical</td>
    </tr>
    <tr>
      <th>writeFluxLog</th>
      <td>False</td>
      <td>logical</td>
    </tr>
    <tr>
      <th>writeMediaLog</th>
      <td>False</td>
      <td>logical</td>
    </tr>
    <tr>
      <th>writeSpecificMediaLog</th>
      <td>False</td>
      <td>logical</td>
    </tr>
    <tr>
      <th>writeTotalBiomassLog</th>
      <td>True</td>
      <td>logical</td>
    </tr>
  </tbody>
</table>
<p>62 rows × 2 columns</p>
</div>



### Preparing a `COMETS` simulation layout
The layout class describes the characteristics of the environment, i.e. the "world", including which species (models) are in it. It can be instantiated in empty or using `COMETS` models: 

 * If instantiated without arguments (as `my_layout = c.layout()`), an empty layout is created with all necessary fields that have to be populated.
 * If a layout is instantiated passing a `model` (or several models), it will generate a layout with all metabolites those models can exchange with the environment at zero concentration, plus metals and ions at unlimited concentration (default -1000).  

To examine the different parts of a Comets `layout`, let\'s first create one from the above loaded textbook model: 


```python
my_layout = c.layout(test_model)
```

The layout stores information about the species (`my_layout.models`) and spatial structure (`my_layout.grid`) in the environment. In this case, the model is only the textbook one, and the grid is the default one, which is $1 \times 1$ i.e. only one cell.

The layout stores also information about the **media** as a pandas dataframe. In this case, no amount of any media component is present. 


```python
my_layout.media
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>diff_c</th>
      <th>g_refresh</th>
      <th>g_static</th>
      <th>g_static_val</th>
      <th>init_amount</th>
      <th>metabolite</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0.000005</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>ac_e</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0.000005</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>acald_e</td>
    </tr>
    <tr>
      <th>2</th>
      <td>0.000005</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>akg_e</td>
    </tr>
    <tr>
      <th>3</th>
      <td>0.000005</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>co2_e</td>
    </tr>
    <tr>
      <th>4</th>
      <td>0.000005</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>etoh_e</td>
    </tr>
    <tr>
      <th>5</th>
      <td>0.000005</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>for_e</td>
    </tr>
    <tr>
      <th>6</th>
      <td>0.000005</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>fru_e</td>
    </tr>
    <tr>
      <th>7</th>
      <td>0.000005</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>fum_e</td>
    </tr>
    <tr>
      <th>8</th>
      <td>0.000005</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>glc__D_e</td>
    </tr>
    <tr>
      <th>9</th>
      <td>0.000005</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>gln__L_e</td>
    </tr>
    <tr>
      <th>10</th>
      <td>0.000005</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>glu__L_e</td>
    </tr>
    <tr>
      <th>11</th>
      <td>0.000005</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>h2o_e</td>
    </tr>
    <tr>
      <th>12</th>
      <td>0.000005</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>h_e</td>
    </tr>
    <tr>
      <th>13</th>
      <td>0.000005</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>lac__D_e</td>
    </tr>
    <tr>
      <th>14</th>
      <td>0.000005</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>mal__L_e</td>
    </tr>
    <tr>
      <th>15</th>
      <td>0.000005</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>nh4_e</td>
    </tr>
    <tr>
      <th>16</th>
      <td>0.000005</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>o2_e</td>
    </tr>
    <tr>
      <th>17</th>
      <td>0.000005</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>pi_e</td>
    </tr>
    <tr>
      <th>18</th>
      <td>0.000005</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>pyr_e</td>
    </tr>
    <tr>
      <th>19</th>
      <td>0.000005</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>succ_e</td>
    </tr>
  </tbody>
</table>
</div>



| metabolite | init_amount | diff_c | g_static | g_static_val | g_refresh_val |
|------------|-------------|--------|----------|--------------|---------------|
| ca2_e      | 1000        | NaN    | 1        | 1000         | 0             |
| cbl_e      | 1000        | NaN    | 1        | 1000         | 0             |
| cl_e       | 1000        | NaN    | 1        | 1000         | 0             |
| . . .      | . . .       | . . .  | . . .    | . . .        | . . .         |


When initiated from models, the media compounds that can be in the environment are all those for which there is an exchange reaction in at least one of the models. The media, shown in the table above, is a `pandas` dataframe where several pieces of information are stored:

* `init_amount` is the initial amount to be added to each cell of the simulation grid (in mmol). 
* `diff_c` indicates whether the molecule has a diffusion constant different than the default one (stored in `ec_layout.global_diff`)
* `g_static` indicates whether the component should remain at a static value, i.e. without change due to consumption and other effects of the simulation. This is useful for example, for setting some nutrients as unlimited. 
* `g_static_val` indicates at which value shouold the nutrient remain static, if the previous coulmn value is 1. 
* `g_refresh_val` indicates the amount of the metabolite that should be added after each simulation cycle to each cell of the grid. 


In addition, we can set local `static` and `refresh` values, specific to a cell of the simulation grid. 

When a media component is `static`, this means that its concentration is returned in each cycle to the set static value. This is used when we want a media component to remain virtually unlimited during a simulation. 

When a media component has a `refresh` value, this means it will be replenished by adding the set amount at every simulation cycle. 

Local refresh values are stored in a list, `my_layout.local_refresh`, where each element of the list is itself a list with the form `[ x  y  m1_r  m2_r  m3_r ... ]`. The first two elements `x` and `y` represent the coordinates, and are followed by the refresh values for all metabolites, in the same order as in `media`.

Local static values are stored in a similar way. Each element of the `my_layout.local_static` list is itself a list with the form `[ x  y  m1_s  m1_s_v m2_s  m2_s_v ... ]`. The difference here is that for each metabolite, there are two values, one defining whether the molecule is to be static at that coordinate (`m1_s`, `m2_s`, ... ) and another with the value at which it should be kept (`m1_s_v`, `m2_s_v`, ... ). 

Note that both `local_refresh` and `local_static` can be empty (the default), or contain only entries for the coordinates where there is at least one nonzero refresh or static value, respectively. 

Finally, the layout also contains information about the starting biomass of each model. This information is stored in the `initial_pop` list. Each component of `initial_pop` is itself a list with the format `[x y biomass_1 biomass_2 ...]`specifying the amount of biomass of each model in each coordinate. 

### Running a `COMETS` simulation 

The `comets` class uses a layout object and a parameters object to run simulations and store the output. Running a comets simulation is pretty straightforward. We firstly define the `comets` object by passing it a `layout` and a `params` objects as arguments. Then, we `run()` the simulation: 


```python
my_simulation = c.comets(my_layout, my_params)
my_simulation.run()
```
    
    Running COMETS simulation ...
    Done!


### Checking simulation output and possible errors
In the background, this command invokes the `COMETS` Java engine in a console, giving a standard output (stdout) and standard error (stderr) logs. Both can be acessed through the fields `run_outputs` and `run_errors`, respectively. 


```python
print(my_simulation.run_output)
```

    -script
    running script file: /home/djordje/Dropbox/projects/comets_paper/Methods_Paper_Materials/COMETS_Examples/COMETS_getting_started/Python/.current_script
    Current Java version: 11.0.7
    Loading layout file '/home/djordje/Dropbox/projects/comets_paper/Methods_Paper_Materials/COMETS_Examples/COMETS_getting_started/Python/.current_layout'...
    Found 1 model files!
    Loading '/home/djordje/Dropbox/projects/comets_paper/Methods_Paper_Materials/COMETS_Examples/COMETS_getting_started/Python/e_coli_core.cmd' ...
    Loading '/home/djordje/Dropbox/projects/comets_paper/Methods_Paper_Materials/COMETS_Examples/COMETS_getting_started/Python/e_coli_core.cmd' ...
    Academic license - for non-commercial use only
    Academic license - for non-commercial use only
    Done!
     Testing default parameters...
    Done!
    Optimizer status code = 5 (looks ok!)
    objective solution = 0.8739215069684305
    Constructing world...
    Done!
    medialist	ac_e	acald_e	akg_e	co2_e	etoh_e	for_e	fru_e	fum_e	glc__D_e	gln__L_e	glu__L_e	h2o_e	h_e	lac__D_e	mal__L_e	nh4_e	o2_e	pi_e	pyr_e	succ_e
    Cycle 1
    Total biomass:
    Model e_coli_core.cmd: 1.0E-7
    Cycle 2
    Total biomass:
    Model e_coli_core.cmd: 1.0E-7
    ...
    Total time = 0.312s
    



```python
print(my_simulation.run_errors)
```

    STDERR empty.


### Accessing the results of the simulation
The results of the successful simulation are stored in several fields in the `comets` object, depending on whether the parameters `writeTotalBiomasslog`,  `writeBiomassLog`, `writeFluxLog` and `writeMediaLog` were set to `True`. 

* The field `total_biomass`  stores the total biomass (summed up over all coordinates) for each timepoint and species.
* The field `biomass` stores detailed biomass values for each timepoint, coordinate and species. 
* The field `media` stores the composition of the media at each timepoint.
* The field `fluxes` stores the metabolic fluxes for each species, coordinate and timepoint. 

Additionally, specific comets modes will have additional output fields; for instance, if we run an evolution simulation, the field `genotypes` will store information about each species such as its ancestor and which mutation it suffered. 

All of the output files ae `pandas` dataframes which can be further analyzed or plotted using standard Python tools.  


```python
my_simulation.total_biomass
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>cycle</th>
      <th>e_coli_core</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0.0</td>
      <td>1.000000e-07</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1.0</td>
      <td>1.000000e-07</td>
    </tr>
    <tr>
      <th>2</th>
      <td>2.0</td>
      <td>1.000000e-07</td>
    </tr>
    <tr>
      <th>3</th>
      <td>3.0</td>
      <td>1.000000e-07</td>
    </tr>
    <tr>
      <th>4</th>
      <td>4.0</td>
      <td>1.000000e-07</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>96</th>
      <td>96.0</td>
      <td>1.000000e-07</td>
    </tr>
    <tr>
      <th>97</th>
      <td>97.0</td>
      <td>1.000000e-07</td>
    </tr>
    <tr>
      <th>98</th>
      <td>98.0</td>
      <td>1.000000e-07</td>
    </tr>
    <tr>
      <th>99</th>
      <td>99.0</td>
      <td>1.000000e-07</td>
    </tr>
    <tr>
      <th>100</th>
      <td>100.0</td>
      <td>1.000000e-07</td>
    </tr>
  </tbody>
</table>
<p>101 rows × 2 columns</p>
</div>


