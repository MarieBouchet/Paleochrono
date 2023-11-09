#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  7 16:55:54 2019
Module for the SitePair class.
@author: parrenif
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as mpl
from matplotlib.backends.backend_pdf import PdfPages
from scipy.linalg import lu_factor, lu_solve
from scipy.linalg import cholesky
import pccfg

class SitePair(object):
    """Class for a pair of sites."""

    def __init__(self, site1, site2):
        self.site1 = site1
        self.site2 = site2
        self.label = self.site1.label+'-'+self.site2.label

        self.age_age_label = self.site1.age_label+self.site2.age_label
        if len(self.age_age_label)>0:
            self.age_age_label = self.age_age_label + '_'
        self.age_age2_label = self.site1.age_label+self.site2.age2_label + '_'
        self.age2_age_label = self.site1.age2_label+self.site2.age_label + '_'
        self.age2_age2_label = self.site1.age2_label+self.site2.age2_label + '_'
        

#TODO: allow to have either dlabel1+'-'dlabel2 or dlbel2+'-'dlabel1 as directory
        filename = pccfg.datadir+self.site1.label+'-'+self.site2.label+\
            '/'+self.age_age_label+'synchro_horizons.txt'
        if not os.path.isfile(filename):
            filename = pccfg.datadir+self.site1.label+'-'+self.site2.label+'/ice_depth.txt'
        if os.path.isfile(filename) and open(filename).read():
            readarray = np.loadtxt(filename)
            if np.size(readarray) == np.shape(readarray)[0]:
                readarray.resize(1, np.size(readarray))
            self.iceicehorizons_depth1 = readarray[:, 0]
            self.iceicehorizons_depth2 = readarray[:, 1]
            self.iceicehorizons_sigma = readarray[:, 2]
        else:
            self.iceicehorizons_depth1 = np.array([])
            self.iceicehorizons_depth2 = np.array([])
            self.iceicehorizons_sigma = np.array([])
        self.iceicehorizons_correlation = np.diag(np.ones(np.size(self.iceicehorizons_depth1)))

        if self.site1.archive == 'icecore' and self.site2.archive == 'icecore':
            filename = pccfg.datadir+self.site1.label+'-'+self.site2.label+\
                       '/'+self.age2_age2_label+'synchro_horizons.txt'
            if not os.path.isfile(filename):
                filename = pccfg.datadir+self.site1.label+'-'+self.site2.label+'/air_depth.txt'
            if os.path.isfile(filename) and open(filename).read():
                readarray = np.loadtxt(filename)
                if np.size(readarray) == np.shape(readarray)[0]:
                    readarray.resize(1, np.size(readarray))
                self.airairhorizons_depth1 = readarray[:, 0]
                self.airairhorizons_depth2 = readarray[:, 1]
                self.airairhorizons_sigma = readarray[:, 2]
            else:
                self.airairhorizons_depth1 = np.array([])
                self.airairhorizons_depth2 = np.array([])
                self.airairhorizons_sigma = np.array([])
            self.airairhorizons_correlation = np.diag(np.ones(np.size(self.airairhorizons_depth1)))

        if self.site2.archive == 'icecore':
            filename = pccfg.datadir+self.site1.label+'-'+\
                            self.site2.label+'/'+self.age_age2_label+'synchro_horizons.txt'
            if not os.path.isfile(filename):
                filename = pccfg.datadir+self.site1.label+'-'+self.site2.label+'/iceair_depth.txt'
            if os.path.isfile(filename) and open(filename).read():
                readarray = np.loadtxt(filename)
                if np.size(readarray) == np.shape(readarray)[0]:
                    readarray.resize(1, np.size(readarray))
                self.iceairhorizons_depth1 = readarray[:, 0]
                self.iceairhorizons_depth2 = readarray[:, 1]
                self.iceairhorizons_sigma = readarray[:, 2]
            else:
                self.iceairhorizons_depth1 = np.array([])
                self.iceairhorizons_depth2 = np.array([])
                self.iceairhorizons_sigma = np.array([])
            self.iceairhorizons_correlation = np.diag(np.ones(np.size(self.iceairhorizons_depth1)))

        if self.site1.archive == 'icecore':
            filename = pccfg.datadir+self.site1.label+'-'+\
                        self.site2.label+'/'+self.age2_age_label+'synchro_horizons.txt'
            if not os.path.isfile(filename):
                filename = pccfg.datadir+self.site1.label+'-'+self.site2.label+'/airice_depth.txt'
            if os.path.isfile(filename) and open(filename).read():
                readarray = np.loadtxt(filename)
                if np.size(readarray) == np.shape(readarray)[0]:
                    readarray.resize(1, np.size(readarray))
                self.airicehorizons_depth1 = readarray[:, 0]
                self.airicehorizons_depth2 = readarray[:, 1]
                self.airicehorizons_sigma = readarray[:, 2]
            else:
                self.airicehorizons_depth1 = np.array([])
                self.airicehorizons_depth2 = np.array([])
                self.airicehorizons_sigma = np.array([])
            self.airicehorizons_correlation = np.diag(np.ones(np.size(self.airicehorizons_depth1)))


        filename1 = pccfg.datadir+'/parameters_covariance_observations_all_site_pairs.py'
        filename2 = pccfg.datadir+'/parameters-CovarianceObservations-AllDrillings.py'
        if os.path.isfile(filename1):
            exec(open(filename1).read())
        elif os.path.isfile(filename2):
            exec(open(filename2).read())
        filename3 = pccfg.datadir+self.label+'/parameters_covariance_observations.py'
        filename4 = pccfg.datadir+self.label+'/parameters-CovarianceObservations.py'
        if os.path.isfile(filename3):
            exec(open(filename3).read())
        elif os.path.isfile(filename4):
            exec(open(filename4).read())
            
        if ((os.path.isfile(filename1) or os.path.isfile(filename2) or os.path.isfile(filename3)\
            or os.path.isfile(filename4)) and (pccfg.jacobian=='analytical' or \
            pccfg.jacobian=='semi_adjoint' or pccfg.jacobian=='adjoint')):
            print('Covariance on observations not implemented for analytical Jacobian. Exiting.')
            sys.exit()
            
        if np.any(self.iceicehorizons_correlation != \
                  np.diag(np.ones(np.size(self.iceicehorizons_depth1)))):
            self.iceicehorizons_correlation_bool = True
            self.iceicehorizons_chol = cholesky(self.iceicehorizons_correlation)
            self.iceicehorizons_lu_piv = lu_factor(self.iceicehorizons_chol)
        else:
            self.iceicehorizons_correlation_bool = False
            
        if self.site1.archive == 'icecore' and self.site2.archive == 'icecore':
            if np.any(self.airairhorizons_correlation != \
                  np.diag(np.ones(np.size(self.airairhorizons_depth1)))):
                self.airairhorizons_correlation_bool = True
                self.airairhorizons_chol = cholesky(self.airairhorizons_correlation)
                self.airairhorizons_lu_piv = lu_factor(self.airairhorizons_chol)
            else:
                self.airairhorizons_correlation_bool = False
                
        if self.site2.archive == 'icecore':
            if np.any(self.iceairhorizons_correlation != \
                  np.diag(np.ones(np.size(self.iceairhorizons_depth1)))):
                self.iceairhorizons_correlation_bool = True
                self.iceairhorizons_chol = cholesky(self.iceairhorizons_correlation)
                self.iceairhorizons_lu_piv = lu_factor(self.iceairhorizons_chol)
            else:
                self.iceairhorizons_correlation_bool = False
                
        if self.site1.archive == 'icecore':
            if np.any(self.airicehorizons_correlation != \
                  np.diag(np.ones(np.size(self.airicehorizons_depth1)))):
                self.airicehorizons_correlation_bool = True
                self.airicehorizons_chol = cholesky(self.airicehorizons_correlation)
                self.airicehorizons_lu_piv = lu_factor(self.airicehorizons_chol)
            else:
                self.airicehorizons_correlation_bool = False

    def residuals(self):
        """Calculate the residual terms of a pair of sites."""

        if np.size(self.iceicehorizons_depth1) > 0:
            resi_iceice = (self.site1.fct_age(self.iceicehorizons_depth1)-\
                           self.site2.fct_age(self.iceicehorizons_depth2))/self.iceicehorizons_sigma
            if self.iceicehorizons_correlation_bool:
                resi_iceice = lu_solve(self.iceicehorizons_lu_piv, resi_iceice)
            resi = [resi_iceice]
        else:
            resi = [np.array([])]

        if self.site1.archive == 'icecore' and self.site2.archive == 'icecore' and \
            np.size(self.airairhorizons_depth1) > 0:
            resi_airair = (self.site1.fct_airage(self.airairhorizons_depth1)-\
                          self.site2.fct_airage(self.airairhorizons_depth2))/\
                          self.airairhorizons_sigma
            if self.airairhorizons_correlation_bool:
                resi_airair = lu_solve(self.airairhorizons_lu_piv, resi_airair)
            resi.append(resi_airair)

        if self.site2.archive == 'icecore' and np.size(self.iceairhorizons_depth1) > 0:
            resi_iceair = (self.site1.fct_age(self.iceairhorizons_depth1)-\
                          self.site2.fct_airage(self.iceairhorizons_depth2))/\
                          self.iceairhorizons_sigma
            if self.iceairhorizons_correlation_bool:
                resi_iceair = lu_solve(self.iceairhorizons_lu_piv, resi_iceair)
            resi.append(resi_iceair)

        if self.site1.archive == 'icecore' and np.size(self.airicehorizons_depth1) > 0:
            resi_airice = (self.site1.fct_airage(self.airicehorizons_depth1)-\
                           self.site2.fct_age(self.airicehorizons_depth2))/self.airicehorizons_sigma
            if self.airicehorizons_correlation_bool:
                resi_airice = lu_solve(self.airicehorizons_lu_piv, resi_airice)
            resi.append(resi_airice)

        return np.concatenate(resi)


    def residuals_jacobian1(self):
#        if np.size(self.iceicehorizons_depth1) > 0:
#            print(np.shape(self.site1.fct_age_jac(self.iceicehorizons_depth1)),
#                  np.shape(self.site2.fct_age(self.iceicehorizons_depth2)),
#                  np.shape(self.iceicehorizons_sigma)
        resi_iceice = self.site1.fct_age_jac(self.iceicehorizons_depth1)/self.iceicehorizons_sigma
        resi = [resi_iceice]
        if self.site1.archive == 'icecore' and self.site2.archive == 'icecore':
            resi_airair = self.site1.fct_airage_jac(self.airairhorizons_depth1)/\
                            self.airairhorizons_sigma
            resi.append(resi_airair)
        if self.site2.archive == 'icecore':
            resi_iceair = self.site1.fct_age_jac(self.iceairhorizons_depth1)/\
                            self.iceairhorizons_sigma
            resi.append(resi_iceair)
        if self.site1.archive == 'icecore':
            resi_airice = self.site1.fct_airage_jac(self.airicehorizons_depth1)/\
                            self.airicehorizons_sigma
            resi.append(resi_airice)
#       else:
#            resi = [np.array([])]
        return np.concatenate(resi, axis=1)

    def residuals_jacobian2(self):
#        if np.size(self.iceicehorizons_depth1) > 0:
        resi_iceice = -self.site2.fct_age_jac(self.iceicehorizons_depth2)/self.iceicehorizons_sigma
        resi = [resi_iceice]
        if self.site1.archive == 'icecore' and self.site2.archive == 'icecore':
            resi_airair = -self.site2.fct_airage_jac(self.airairhorizons_depth2)/\
                            self.airairhorizons_sigma
            resi.append(resi_airair)
        if self.site2.archive == 'icecore':
            resi_iceair = -self.site2.fct_airage_jac(self.iceairhorizons_depth2)/\
                            self.iceairhorizons_sigma
            resi.append(resi_iceair)
        if self.site1.archive == 'icecore':
            resi_airice = -self.site2.fct_age_jac(self.airicehorizons_depth2)/\
                            self.airicehorizons_sigma
            resi.append(resi_airice)
#        else:
#            resi = [np.array([])]
        return np.concatenate(resi, axis=1)

    def residuals_delta(self):
        resi_iceice = (self.site1.fct_age_delta(self.iceicehorizons_depth1)-\
            self.site2.fct_age_delta(self.iceicehorizons_depth2))/self.iceicehorizons_sigma
        return resi_iceice

    def figures(self):
        
        """Build the figures related to a pair of sites."""
        if np.size(self.iceicehorizons_depth1)>0:
            fig, ax = mpl.subplots()
            mpl.xlabel(self.site1.label+' '+self.site1.age_labelsp+'age ('+pccfg.age_unit+' '+pccfg.age_unit_ref+')')
            mpl.ylabel(self.site2.label+' '+self.site2.age_labelsp+'age ('+pccfg.age_unit+' '+pccfg.age_unit_ref+')')
            if np.size(self.iceicehorizons_depth1) > 0:
                if pccfg.show_initial:
                    mpl.plot(self.site1.fct_age_init(self.iceicehorizons_depth1),
                                 self.site2.fct_age_init(self.iceicehorizons_depth2),
                                 color=pccfg.color_init, linestyle='', marker='o', markersize=2,
                                 label="Initial")
                mpl.plot(self.site1.fct_age_model(self.iceicehorizons_depth1),
                             self.site2.fct_age_model(self.iceicehorizons_depth2),
                             color=pccfg.color_mod, linestyle='', marker='o', markersize=2,
                             label="Prior")
                mpl.errorbar(self.site1.fct_age(self.iceicehorizons_depth1),
                             self.site2.fct_age(self.iceicehorizons_depth2), color=pccfg.color_opt,
                             ecolor=pccfg.color_ci,
                             xerr=np.zeros(np.size(self.iceicehorizons_depth1)),
                             linestyle='', marker='o', markersize=2,
                             label="Posterior $\pm\sigma$")
                xstart = self.site1.fct_age(self.iceicehorizons_depth1)-self.iceicehorizons_sigma/2
                ystart = self.site2.fct_age(self.iceicehorizons_depth2)+self.iceicehorizons_sigma/2
                for i in range(np.size(self.iceicehorizons_depth1)):
                    mpl.arrow(xstart[i], ystart[i], self.iceicehorizons_sigma[i],
                              -self.iceicehorizons_sigma[i], color=pccfg.color_ci,
                              width=0.0, head_length=0.0, head_width=0.0)
            x_low, x_up, y_low, y_up = mpl.axis()
#            x_low = self.site1.age_top
#            y_low = self.site2.age_top
#            mpl.axis((x_low, x_up, y_low, y_up))
            rangefig = np.array([min(x_low, y_low), max(x_up, y_up)])
            mpl.plot(rangefig, rangefig, color=pccfg.color_obs, label='Perfect agreement', zorder=0)
            mpl.legend(loc="best")
            ax.set_aspect('equal')
            printed_page = PdfPages(pccfg.datadir+self.label+'/'+self.age_age_label+'synchro.pdf')
            printed_page.savefig(fig)
            printed_page.close()
            if not pccfg.show_figures:
                mpl.close()

        if self.site1.archive == 'icecore' and self.site2.archive == 'icecore':
            if np.size(self.airairhorizons_depth1)>0:
                fig, ax = mpl.subplots()
                mpl.xlabel(self.site1.label+' '+self.site1.age2_labelsp+'age ('+pccfg.age_unit+' '+
                           pccfg.age_unit_ref+')')
                mpl.ylabel(self.site2.label+' '+self.site2.age2_labelsp+'age ('+pccfg.age_unit+' '+
                           pccfg.age_unit_ref+')')
                if np.size(self.airairhorizons_depth1) > 0:
                    if pccfg.show_initial:
                        mpl.plot(self.site1.fct_airage_init(self.airairhorizons_depth1),
                                     self.site2.fct_airage_init(self.airairhorizons_depth2),
                                     color=pccfg.color_init,
                                     linestyle='',
                                     marker='o', markersize=2, label="Initial")
                    mpl.plot(self.site1.fct_airage_model(self.airairhorizons_depth1),
                                 self.site2.fct_airage_model(self.airairhorizons_depth2),
                                 color=pccfg.color_mod,
                                 linestyle='', marker='o', markersize=2,
                                 label="Prior")
                    mpl.errorbar(self.site1.fct_airage(self.airairhorizons_depth1),
                                 self.site2.fct_airage(self.airairhorizons_depth2),
                                 color=pccfg.color_opt, ecolor=pccfg.color_ci,
                                 xerr=np.zeros_like(self.airairhorizons_sigma),
                                 linestyle='', marker='o', markersize=2,
                                 label="Posterior $\pm\sigma$")
                    xstart = self.site1.fct_airage(self.airairhorizons_depth1)-\
                                 self.airairhorizons_sigma/2
                    ystart = self.site2.fct_airage(self.airairhorizons_depth2)+\
                                 self.airairhorizons_sigma/2
                    for i in range(np.size(self.airairhorizons_depth1)):
                        mpl.arrow(xstart[i], ystart[i], self.airairhorizons_sigma[i],
                                  -self.airairhorizons_sigma[i], color=pccfg.color_ci,
                                  width=0.0, head_length=0.0, head_width=0.0)
                x_low, x_up, y_low, y_up = mpl.axis()
#                x_low = self.site1.age_top
#                y_low = self.site2.age_top
#                mpl.axis((x_low, x_up, y_low, y_up))
                rangefig = np.array([min(x_low, y_low), max(x_up, y_up)])
                mpl.plot(rangefig, rangefig, color=pccfg.color_obs, label='Perfect agreement',
                         zorder=0)
                mpl.legend(loc="best")
                ax.set_aspect('equal')
                printed_page = PdfPages(pccfg.datadir+self.label+'/'+self.age2_age2_label+
                                        'synchro.pdf')
                printed_page.savefig(fig)
                printed_page.close()
                if not pccfg.show_figures:
                    mpl.close()

        if self.site2.archive == 'icecore':
            if np.size(self.iceairhorizons_depth1)>0:
                fig, ax = mpl.subplots()
                mpl.xlabel(self.site1.label+' '+self.site1.age_labelsp+'age ('+pccfg.age_unit+' '+
                           pccfg.age_unit_ref+')')
                mpl.ylabel(self.site2.label+' '+self.site2.age2_labelsp+'age ('+pccfg.age_unit+' '+
                           pccfg.age_unit_ref+')')
                if np.size(self.iceairhorizons_depth1) > 0:
                    if pccfg.show_initial:
                        mpl.plot(self.site1.fct_age_init(self.iceairhorizons_depth1),
                                     self.site2.fct_airage_init(self.iceairhorizons_depth2),
                                     color=pccfg.color_init,
                                     linestyle='',
                                     marker='o', markersize=2, label="Initial")
                    mpl.plot(self.site1.fct_age_model(self.iceairhorizons_depth1),
                                 self.site2.fct_airage_model(self.iceairhorizons_depth2),
                                 color=pccfg.color_mod,
                                 linestyle='', marker='o', markersize=2,
                                 label="Prior")
                    mpl.errorbar(self.site1.fct_age(self.iceairhorizons_depth1),
                                 self.site2.fct_airage(self.iceairhorizons_depth2),
                                 color=pccfg.color_opt, ecolor=pccfg.color_ci,
                                 xerr=np.zeros_like(self.iceairhorizons_sigma),
                                 linestyle='', marker='o', markersize=2,
                                 label="Posterior $\pm\sigma$")
                    xstart = self.site1.fct_age(self.iceairhorizons_depth1)-\
                                 self.iceairhorizons_sigma/2
                    ystart = self.site2.fct_airage(self.iceairhorizons_depth2)+\
                                 self.iceairhorizons_sigma/2
                    for i in range(np.size(self.iceairhorizons_depth1)):
                        mpl.arrow(xstart[i], ystart[i], self.iceairhorizons_sigma[i],
                                  -self.iceairhorizons_sigma[i], color=pccfg.color_ci,
                                  width=0.0, head_length=0.0, head_width=0.0)                    
                x_low, x_up, y_low, y_up = mpl.axis()
#                x_low = self.site1.age_top
#                y_low = self.site2.age_top
#                mpl.axis((x_low, x_up, y_low, y_up))
                rangefig = np.array([min(x_low, y_low), max(x_up, y_up)])
                mpl.plot(rangefig, rangefig, color=pccfg.color_obs, label='Perfect agreement',
                         zorder=0)
                mpl.legend(loc="best")
                ax.set_aspect('equal')
                printed_page = PdfPages(pccfg.datadir+self.label+'/'+self.age_age2_label+
                                        'synchro.pdf')
                printed_page.savefig(fig)
                printed_page.close()
                if not pccfg.show_figures:
                    mpl.close()

        if self.site1.archive == 'icecore':
            if np.size(self.airicehorizons_depth1)>0:
                fig, ax = mpl.subplots()
                mpl.xlabel(self.site1.label+' '+self.site1.age2_labelsp+'age ('+pccfg.age_unit+' '+
                           pccfg.age_unit_ref+')')
                mpl.ylabel(self.site2.label+' '+self.site2.age_labelsp+'age ('+pccfg.age_unit+' '+
                           pccfg.age_unit_ref+')')
                if np.size(self.airicehorizons_depth1) > 0:
                    if pccfg.show_initial:
                        mpl.plot(self.site1.fct_airage_init(self.airicehorizons_depth1),
                                     self.site2.fct_age_init(self.airicehorizons_depth2),
                                     color=pccfg.color_init,
                                     linestyle='', marker='o', markersize=2, label="Initial")
                    mpl.plot(self.site1.fct_airage_model(self.airicehorizons_depth1),
                                 self.site2.fct_age_model(self.airicehorizons_depth2),
                                 color=pccfg.color_mod,
                                 linestyle='', marker='o', markersize=2,
                                 label="Prior")
                    mpl.errorbar(self.site1.fct_airage(self.airicehorizons_depth1),
                                 self.site2.fct_age(self.airicehorizons_depth2),
                                 color=pccfg.color_opt, ecolor=pccfg.color_ci,
                                 xerr=np.zeros_like(self.airicehorizons_sigma),
                                 linestyle='', marker='o', markersize=2,
                                 label="Posterior $\pm\sigma$")
                    xstart = self.site1.fct_airage(self.airicehorizons_depth1)-\
                                 self.airicehorizons_sigma/2
                    ystart = self.site2.fct_age(self.airicehorizons_depth2)+\
                                 self.airicehorizons_sigma/2
                    for i in range(np.size(self.airicehorizons_depth1)):
                        mpl.arrow(xstart[i], ystart[i], self.airicehorizons_sigma[i],
                                  -self.airicehorizons_sigma[i], color=pccfg.color_ci,
                                  width=0.0, head_length=0.0, head_width=0.0)
                x_low, x_up, y_low, y_up = mpl.axis()
#                x_low = self.site1.age_top
#                y_low = self.site2.age_top
#                mpl.axis((x_low, x_up, y_low, y_up))
                rangefig = np.array([min(x_low, y_low), max(x_up, y_up)])
                mpl.plot(rangefig, rangefig, color=pccfg.color_obs, label='Perfect agreement')
                mpl.legend(loc="best")
                ax.set_aspect('equal')
                printed_page = PdfPages(pccfg.datadir+self.label+'/'+self.age2_age_label+
                                        'synchro.pdf')
                printed_page.savefig(fig)
                printed_page.close()
                if not pccfg.show_figures:
                    mpl.close()
