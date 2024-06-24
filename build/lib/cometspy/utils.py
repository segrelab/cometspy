#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

The utils module contains helper functions generating spatial patterns. 
"""

from .layout import layout
from .params import params
import random

def pick_random_locations(n : int, 
                          xrange : tuple, yrange : tuple, 
                          forbidden_locs : set = set()) -> list:
    """
    returns a list of n x,y tuples corresponding to locations in the range given

    Parameters
    ----------

    n : int
        number of locations desired
    xrange : tuple
        the x-range (min, max) of the x range possible
    yrange : tuple
        the y-range (min, max) of the y range possible
    forbidden_locs : set, optional
        A list of tuples that cannot be chosen. 

    Returns
    -------

    list
        a list of (x,y) values. 
        
    Examples
    --------

    >>> from cometspy.utils import pick_random_locations
    >>> locs = pick_random_locations(3, (0, 10), (0,10))
    >>> locs

    """
    pickable_locs = [(x,y) for x in range(xrange[0], xrange[1]) for y in range(yrange[0], yrange[1])]
    pickable_locs = list(set(pickable_locs).difference(set(forbidden_locs)))
    if len(pickable_locs) < n:
        print("there are fewer available locations than n, returning all available locs")
        locs = pickable_locs
    else:
        locs = random.sample(pickable_locs, n)
       
    return(locs)
    
def grow_rocks(n : int, 
               xrange : tuple, 
               yrange : tuple, 
               mean_size: int) -> list:
    """
    grows simple simulated rocks by adding random adjacent points from seeds
    
    n number of seed points are generated first with pick_random_locations. 
    Then, mean_size * n - n additional points are added to these seed locations. 
    For each new point, one random location out of the set of all possible 
    unoccupied locations next to occupied locations are chosen. Only 
    lattice points directly to the NSEW are considered. This process is
    repeated until all new points are assigned. This usually results in
    'rocks' of different shapes and sizes. 
    
    This function can be very slow (> 1 min) when n * mean_size > 2000

    Parameters
    ----------
    n : int
        number of seed points to generate rocks from
    xrange : tuple
        x range possible, e.g. (0, 5)
    yrange : tuple
        y range possible, e.g. (0, 10)
    mean_size : int
        average size in lattice units of a generated rock

    Returns
    -------
    list
        list of all points generated including seed points

    """
    grow_points = n * mean_size - n
    locs = pick_random_locations(n, xrange, yrange)
    if (n*mean_size) > ((xrange[1]-xrange[0]) * (yrange[1]-yrange[0])):
        print("more rocks asked for than space possible, try again")
        return
    for i in range(grow_points):
        adjacents = _find_unoccupied_adjacent(locs, xrange, yrange)
        locs.append(random.sample(adjacents, 1)[0])
    return(locs)

def _find_unoccupied_adjacent(locs, xrange, yrange):
    """
    returns the unoccupied adjacent locations to the given locations
    
    this is not the set, so unoccupied locations adjacent to > 1 occupied
    location will have more spots
    """ 
    locs = [(loc[0], loc[1]) for loc in locs]
    adjacent = []
    for loc in locs:
        if ((loc[0]-1, loc[1]) not in locs) and ((loc[0]-1) >= xrange[0]):
            adjacent.append((loc[0]-1, loc[1]))
        if ((loc[0]+1, loc[1]) not in locs) and ((loc[0]+1) <= xrange[1]):
            adjacent.append((loc[0]+1, loc[1]))
        if ((loc[0], loc[1]-1) not in locs) and ((loc[1]-1) >= yrange[0]):
            adjacent.append((loc[0], loc[1]-1))
        if ((loc[0], loc[1]+1) not in locs) and ((loc[1]+1) <= yrange[1]):
            adjacent.append((loc[0], loc[1]+1))
    #adjacent = list(set(adjacent))
    return(adjacent)

def chemostat(models : list, 
              reservoir_media : dict, 
              dilution_rate : float) -> tuple:
    """
    helper function to let a user skip some steps when generating a chemostat
    
    This sets relevant simulation parameters (e.g. deathRate, 
    metaboliteDilutionRate) and layout values (e.g. refresh media) based upon
    a "reservoir" definition and a dilution rate.
    
    It generates a layout that has the reservoir_media as the initial values,
    as well as set it to drip in / out based upon the dilution rate. 
    
    The returned layout and params can be further modified before supplying
    to a comets object if desired. 

    Parameters
    ----------
    
    models : list(cometspy.model)
        list of cometspy.model(s) with initial_pop set to use in the sim
    reservoir_media : dict
        media definition with metabolite names as keys as mmol amounts as values
    dilution_rate : float
        the dilution rate of the chemostat, in 1/hr

    Returns
    -------
    
    tuple (layout, params)
        a cometspy.layout object and a cometspy.params object
        
    Examples
    --------
    
    >>> import cobra.test
    >>> import cometspy as c
    >>> from cometspy.utils import chemostat
    >>> # make a model from a cobra model, open exchange reactions, and give a pop
    >>> tb = cobra.test.create_test_model("textbook")
    >>> m = c.model(tb)
    >>> m.initial_pop = [0, 0, 1.e-4]
    >>> m.open_exchanges()
    >>> reservoir = {'glc__D_e' : 0.01, 'nh4_e' : 1000., 'pi_e' : 1000.}
    >>> layout, params = chemostat([m], reservoir, 0.1)
    >>> params.set_param("maxCycles", 100)
    >>> sim = c.comets(layout, params)
    >>> sim.run()
    >>> print(sim.total_biomass)

    """
    mylayout = layout(models)
    for key, value in reservoir_media.items():
        mylayout.set_specific_metabolite(key, value)
        mylayout.set_specific_refresh(key, value * dilution_rate)

    parameters = params()
    parameters.all_params['metaboliteDilutionRate'] = dilution_rate
    parameters.all_params['deathRate'] = dilution_rate
    return((mylayout, parameters))
