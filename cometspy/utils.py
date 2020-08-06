#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

The utils module contains helper functions for running COMETS simulations
For more information see https://segrelab.github.io/comets-manual/
"""

from .layout import layout
from .params import params
import random

def pick_random_locations(n, xrange, yrange, forbidden_locs = set()):
    pickable_locs = [(x,y) for x in range(xrange[0], xrange[1]) for y in range(yrange[0], yrange[1])]
    pickable_locs = list(set(pickable_locs).difference(set(forbidden_locs)))
    if len(pickable_locs) < n:
        print("there are fewer available locations than n, returning all available locs")
        locs = pickable_locs
    else:
        locs = random.sample(pickable_locs, n)
       
    return(locs)
    
def grow_rocks(n, xrange, yrange, mean_size):
    grow_points = n * mean_size - n
    locs = pick_random_locations(n, xrange, yrange)
    if (n*mean_size) > ((xrange[1]-xrange[0]) * (yrange[1]-yrange[0])):
        print("more rocks asked for than space possible, try again")
        return
    while grow_points > 0:
        new_locs = []
        for loc in locs:
            for dir in range(4):
                if dir == 0:
                    if loc[0]+1 >= xrange[0]:
                        continue
                    new_loc = (loc[0]+1, loc[1])
                elif dir == 1:
                    if loc[0]-1 <= xrange[0]:
                        continue
                    new_loc = (loc[0]-1, loc[1])
                elif dir == 2:
                    if loc[1]-1 <= xrange[1]:
                        continue
                    new_loc = (loc[0], loc[1]-1)
                elif dir == 3:
                    if loc[1]+1 >= xrange[1]:
                        continue
                    new_loc = (loc[0], loc[1]+1)
                if new_loc not in locs:
                    new_locs.append(new_loc)
        locs.append(random.sample(new_locs,1)[0])
        grow_points -= 1
    return(locs)

def chemostat(models, reservoir_media, dilution_rate):
    mylayout = layout(models)
    for key, value in reservoir_media.items():
        mylayout.set_specific_metabolite(key, value)
        mylayout.set_specific_refresh(key, value * dilution_rate)

    parameters = params()
    parameters.all_params['metaboliteDilutionRate'] = dilution_rate
    parameters.all_params['deathRate'] = dilution_rate
    return((mylayout, parameters))
