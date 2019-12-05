#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 10:37:31 2019

@author: chaco001
"""

import comets as c
# import os
# import pandas as pd
# import cobra as cb
# import cobra.test
# import copy

b_params = c.params()
b_params.all_params['maxCycles'] = 400

p = c.model('test_models/iJO1366.xml')
p.initial_pop = [0, 0, 1e-8]
layout = c.layout()
layout.add_model(p)
layout.add_typical_trace_metabolites()
layout.set_specific_metabolite('glc__D_e', 10.)
layout.display_current_media()
batch_test = c.comets(layout, b_params)
#batch_test.set_classpath("lang3", "/opt/commons-lang3-3.8.1/commons-lang3-3.8.1.jar")
#batch_test.set_classpath("concurrent", "/opt/colt/lib/concurrent.jar")
#batch_test.set_classpath("colt", "/opt/colt/lib/colt.jar")
batch_test.run()
batch_test.run_output.decode()
myplot = batch_test.total_biomass.plot(x = 'cycle')


layout = c.layout()
layout.add_model(o)
p.id = 'modified_model'
p.change_vmax('EX_glc__D_e', 8)
layout.add_model(p)
layout.add_typical_trace_metabolites()
layout.set_specific_metabolite('glc__D_e', 10.)
layout.display_current_media()
batch_test = c.comets(layout, b_params)
batch_test.set_classpath("lang3", "/opt/commons-lang3-3.8.1/commons-lang3-3.8.1.jar")
batch_test.set_classpath("concurrent", "/opt/colt/lib/concurrent.jar")
batch_test.set_classpath("colt", "/opt/colt/lib/colt.jar")
batch_test.run()
batch_test.run_output
batch_test.total_biomass.plot(x = 'cycle')

