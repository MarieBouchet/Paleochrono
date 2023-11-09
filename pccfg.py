#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  7 11:52:46 2019
Global variables used across all modules.
@author: parrenif
"""
import sys
import os
import yaml

##Default Parameters
list_sites = []
opt_method = 'trf'  #leastsq, leastsq-parallel, none
is_parallel = False
is_analytical_jacobian = True
jacobian = 'semi_adjoint'
nb_nodes = 6         #Number of nodes for the leastsq-parallel mode
datadir = './'
color_obs = 'r'       #color for the observations
color_opt = 'k'       #color for the posterior scenario
color_mod = 'b'       #color for the prior scenario
color_ci = '0.8'      #color for the confidence intervals
color_sigma = 'm'     #color for the uncertainty
color_di = 'g'        #color for the dated intervals
color_resolution = 'y' #Color for the resolution
show_initial = False  #always put to False for now
color_init = 'c'      #always put to 'c' for now
scale_ageci = 10.     #scaling of the confidence interval in the ice and air age figures
show_figures = False  #whether to show or not the figures at the end of the run
show_airlayerthick = False #whether to show the air layer thickness figure (buggy on anaconda)
tol = 1e-6      #Tolerance for the termination.
tr_solver = 'lsmr'
age_unit = 'yr'
age_unit_ref = 'B1950'
#nb_runs = 0


def read_parameters():
    global list_sites
    global list_drillings
    global datadir

    ###Reading parameters directory
    datadir = sys.argv[1]
    if datadir[-1] != '/':
        datadir = datadir+'/'
    print('Parameters directory is: ', datadir)
    #os.chdir(datadir)
    
    filename = datadir+'parameters.yml'
    if os.path.isfile(filename):
        data = yaml.load(open(filename).read(), Loader=yaml.FullLoader)
        if data != None:
            globals().update(data)
    else:
        print('Python format for the global parameters file is deprecated.'
              'Use YAML format instead.')
        exec(open(datadir+'/parameters.py').read())
        globals().update(locals())
        
    try:
        list_sites = list_drillings
    except NameError:
        pass