# Growth in a test tube
This simple example illustrates the basic workflow of COMETS, including how to build the simulation layout, to specify parameters, load the model and plot the data once the simulation is finished.

The core of the COMETS methodology is the Dynamic Flux Balance Analysis algorithm (Madahevan et al 2002). One of the first successful simulations of the time dynamics of bacterial metabolism was the classical study of *Escherichia coli* batch culture by Varma and Palsson (1994). Here, we reproduce one of the results in that study, the anaerobic fermentation in minimal media with glucose as the only carbon source.

We will simulate a test tube by setting a well-mixed virtual container with $1cm3$ of media, which we will inoculate with $5\times10^{−6}$ grams of E. coli biomass. We will set the initial composition of the substrate to 11mM of glucose and unlimited amounts of ammonia and phosphate. For the nutrient uptake, we will use standard Michaelis-Menten kinetics, using the experimentally measured Monod parameter for anaerobic uptake of glucose by E. coli.

For this simple example, we use the rudimentary "core" model of E. coli (Orth et al. 2009), which can be downloaded from http://bigg.ucsd.edu/models/e_coli_core or loaded, as we do here, using a function built in CobraPy. This model represents an E. coli metabolism simplified to its core functions (glycolysis, tricarboxylic cycle, pentose phosphate shunt, etc).


## Loading the COMETS Python Toolbox

We first import the COMETS Python toolbox, which will also load all the dependencies, such as CobraPy or pandas.


```python
import cometspy as c
import cobra.test
import matplotlib.pyplot as plt
```

## Creating a test tube

We are now ready to create a "layout" for our simulation. By instantiating the class layout without arguments, we will create the default simulation layout, represents an empty, well mixed space (called "cell") with volume 1cm^3. We can then modify this layout according to our needs, in this case adding the media composition in the desired concentration.


```python
# Create empty 1x1 layout
test_tube = c.layout()

# Add 11mM glucose and remove o2
test_tube.set_specific_metabolite('glc__D_e', 0.011)
test_tube.set_specific_metabolite('o2_e', 0)

# Add the rest of nutrients unlimited (ammonia, phosphate, water and protons)
test_tube.set_specific_metabolite('nh4_e',1000);
test_tube.set_specific_metabolite('pi_e',1000);
test_tube.set_specific_metabolite('h2o_e',1000);
test_tube.set_specific_metabolite('h_e',1000);
```

    building empty layout model
    models will need to be added with layout.add_model()
    Warning: The added metabolite (glc__D_e) is notable to be taken up by any of the current models
    Warning: The added metabolite (o2_e) is notable to be taken up by any of the current models
    Warning: The added metabolite (nh4_e) is notable to be taken up by any of the current models
    Warning: The added metabolite (pi_e) is notable to be taken up by any of the current models
    Warning: The added metabolite (h2o_e) is notable to be taken up by any of the current models
    Warning: The added metabolite (h_e) is notable to be taken up by any of the current models


## Loading a model

Next, we have to load the model and add it to the layout (or "inoculate it in our test tube").

We will instantiate the comets model class using a loaded CobraPy model as input.

Note that we remove the bounds on glucose import, which will be set dynamically by COMETS during the simulation according to the dynamically changing external glucose concentration. We will set the initial biomass of our model at $10^{-6}$ grams.


```python
# create the model using CobraPy functionality
e_coli_cobra = cobra.test.create_test_model('textbook')

# use the loaded model to build a comets model
e_coli = c.model(e_coli_cobra)

# remove the bounds from glucose import (will be set dynamically by COMETS)
e_coli.change_bounds('EX_glc__D_e', -1000, 1000)

# set the model's initial biomass
e_coli.initial_pop = [0, 0, 5e-6]

# add it to the test_tube
test_tube.add_model(e_coli)
```

    Using license file /home/djordje/gurobi.lic
    Academic license - for non-commercial use only


## Setting the simulation parameters

We next instantiate the params class, which generates a set of parameters for the COMETS simulation with the [TODO LINK TO DEF VALS] default values for all of them. All of the parameters are contained in the all_params field which is a Python dict object, making it easy to change the value of the desired parameters.


```python
# Set the parameters that are different from the default
sim_params = c.params()
```


```python
sim_params.set_param('defaultVmax', 18.5)
sim_params.set_param('defaultKm', 0.000015)
sim_params.set_param('maxCycles', 1000)
sim_params.set_param('timeStep', 0.01)
sim_params.set_param('spaceWidth', 1)
sim_params.set_param('maxSpaceBiomass', 10)
sim_params.set_param('minSpaceBiomass', 1e-11)
sim_params.set_param('writeMediaLog', True)
```

## Running the simulation

With all set up, we can now instantiate the comets class by passing the layout (containing the model) and the params objects we just created.


```python
experiment = c.comets(test_tube, sim_params)
```


Finally, we can run the simulation as:


```python
experiment.run()
```

    
    Running COMETS simulation ...
    Done!


## Analyzing the results

The results of our simulation are stored in several pandas data frames contained in the comets object that we just simulated. The growth of the simulated model can be seen by plotting the total_biomass field.


```python
ax = experiment.total_biomass.plot(x = 'cycle')
ax.set_ylabel("Biomass (gr.)")
```



![](img/test_tube_1.png)


Similarly, we can plot composition of the media. In this case, we will limit the plot to those components that are not added to the layout in unlimited amounts (“static” compounds, e.g. ammonia, phosphate, water, etc in this simulation). In this case, we do this by limiting the plot to compounds with concentration lower than 900mM.


```python
media = experiment.media.copy()
media = media[media.conc_mmol<900]

fig, ax = plt.subplots()
media.groupby('metabolite').plot(x='cycle', ax =ax, y='conc_mmol')
ax.legend(('acetate','ethanol', 'formate', 'glucose'))
ax.set_ylabel("Concentration (mmol)")
```




 ![](img/test_tube_2.png)

