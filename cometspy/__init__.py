"""
*a python package for making and running COMETS simulations.*

cometspy is a straight-forward python interface to the COMETS program. 
Use for dynamic flux-balance analysis, with multiple models, heterogenous 
spatio-temporal environments, and evolution. 

**To use this package, you must also have the actual COMETS program installed.**
 
* For more information on COMETS, see https://runcomets.org

* For COMETS development, see https://github.com/segrelab/comets

* For development of this package, see http://github.com/segrelab/cometspy

**cometspy workflow:**

models -> layout, layout + params -> comets, comets.run()

**cometspy hello.world**

>>> import cobra.test
>>> import cometspy as c
>>> # make a model from a cobra model, open exchange reactions, and give a pop
>>> tb = cobra.test.create_test_model("textbook")
>>> m = c.model(tb)
>>> m.initial_pop = [0, 0, 1.e-4]
>>> m.open_exchanges()
>>> # make a layout with the model, and give it three nutrients
>>> l = c.layout([m])
>>> l.set_specific_metabolite("glc__D_e", 0.01)
>>> l.set_specific_metabolite("nh4_e", 1000, static = True)
>>> l.set_specific_metabolite("pi_e", 1000, static = True)
>>> # make params and change one default
>>> p = c.params()
>>> p.set_param("maxCycles", 100)
>>> # make a simulation and run it!
>>> sim = c.comets(l, p)
>>> sim.run()
>>> print(sim.total_biomass)
>>> # saved data objects are pandas.DataFrames, and therefore can be plotted
>>> # sim.total_biomass.plot(x = "cycle")

"""
from cometspy.comets import comets
from cometspy.model import model
from cometspy.layout import layout
from cometspy.params import params
