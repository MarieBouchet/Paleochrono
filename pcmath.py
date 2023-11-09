"""
Created on Thu Feb  7 10:00:12 2019
Some mathematical functions for paleochrono.
@author: parrenif
"""
import numpy as np
import math as m

def interp_lin_aver(x_out, x_in, y_in):
    """Return a linear interpolation of a (x_in,y_in) series at x_out abscissas with averaging."""
    y_out = np.nan*np.zeros(np.size(x_out)-1)
    if x_out[0] < min(x_in):
        x_mod = np.concatenate((np.array([x_out[0]]), x_in))
        y_mod = np.concatenate((np.array([y_in[0]]), y_in))
    else:
        x_mod = x_in+0
        y_mod = y_in+0
    if x_out[-1] > max(x_in):
        x_mod = np.concatenate((x_mod, np.array([x_out[-1]])))
        y_mod = np.concatenate((y_mod, np.array([y_in[-1]])))
    for i in range(np.size(x_out)-1):
        x_loc = x_mod[np.where(np.logical_and(x_mod > x_out[i], x_mod < x_out[i+1]))]
        x_loc = np.concatenate((np.array([x_out[i]]), x_loc, np.array([x_out[i+1]])))
        y_loc = np.interp(x_loc, x_mod, y_mod)
        y_out[i] = np.sum((y_loc[1:]+y_loc[:-1])/2*(x_loc[1:]-x_loc[:-1]))/(x_out[i+1]-x_out[i])
    return y_out

def interp_stair_aver(x_out, x_in, y_in):
    """Return a staircase interpolation of a (x_in,y_in) series at x_out abscissas with averaging.
    """
    x_mod = x_in+0
    y_mod = y_in+0
    if x_out[0] < x_in[0]:
        x_mod = np.concatenate((np.array([x_out[0]]), x_mod))
        y_mod = np.concatenate((np.array([y_in[0]]), y_mod))
    if x_out[-1] > x_in[-1]:
        x_mod = np.concatenate((x_mod, np.array([x_out[-1]])))
        y_mod = np.concatenate((y_mod, np.array([y_in[-1]])))
    y_int = np.cumsum(np.concatenate((np.array([0]), y_mod[:-1]*(x_mod[1:]-x_mod[:-1]))))
#Maybe this is suboptimal since we compute twice g(xp[i]):
    y_out = (np.interp(x_out[1:], x_mod, y_int)-np.interp(x_out[:-1], x_mod, y_int))/\
            (x_out[1:]-x_out[:-1])
    return y_out


def gaussian(x_in):
    """Return the value of the gaussian function (no multiplicative constant)
    at a given x_in abscissa."""
    return np.exp(-x_in**2/2)

def grid(para):
    start = para['start']
    end = para['end']
    try:
        nb_steps = para['nb_steps']
    except KeyError:
        resolution = para['resolution']
        nb_steps = m.floor((end-start)/resolution)
        end = start + resolution * nb_steps
    if para['type'] == 'regular':
        eps = (end-start)/nb_steps/2
        grid = np.arange(start, end+eps, (end-start)/nb_steps)
    elif para['type'] == 'linear':
        ratio = para['ratio']
        if ratio == None:
            ratio = 2/(nb_steps+1)
        eps = (1.-ratio)/nb_steps
        grid = np.arange(ratio, 2.-ratio+eps, (2.-2*ratio)/(nb_steps-1))
        grid = grid * (end-start)/nb_steps
        grid = np.cumsum(np.concatenate((np.array([start]), grid)))
    else:
        print('Type of grid not recognized.')
    try:
        inverted = para['inverted']
    except KeyError:
        inverted = False
    if inverted:
        grid = grid[::-1]
        grid = grid[:-1]-grid[1:]
        grid = np.cumsum(np.concatenate((np.array([start]), grid)))
    return grid

def truncation(grid, inf, sup):
    if inf == None:
        inf = grid[0]
    if sup == None:
        sup = grid[-1]
    grid = grid[np.logical_and(grid>=inf, grid<=sup)]
    return grid

def stretch(grid, start, end):
    grid = start + (grid-grid[0])/(grid[-1]-grid[0])*(end-start)
    return grid