#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import numpy as np
from matplotlib import pyplot as plt

#options for plots
plt.rcParams["font.size"]= 16
plt.rcParams['lines.linewidth'] = 3

savefig_options=dict(bbox_inches='tight')

marker_style = dict( marker='o',markersize=8, markerfacecolor='white', linestyle='',markeredgewidth=1.2)

path_data = '../data'
path_plot = '../plots'

# data k0
################################
k0=0
EOCL1cont_k0 = np.loadtxt(path_data+'/kmax='+str(k0)+'/EOCL1cont_k0.txt')
EOCL1dis_k0 = np.loadtxt(path_data+'/kmax='+str(k0)+'/EOCL1dis_k0.txt')
################################


#data k1
################################
k1=1
EOCL1cont_k1 = np.loadtxt(path_data+'/kmax='+str(k1)+'/EOCL1cont_k1.txt')
EOCL1dis_k1 = np.loadtxt(path_data+'/kmax='+str(k1)+'/EOCL1dis_k1.txt')
################################

#data k2
################################
k2=2
EOCL1cont_k2 = np.loadtxt(path_data+'/kmax='+str(k2)+'/EOCL1cont_k2.txt')
EOCL1dis_k2 = np.loadtxt(path_data+'/kmax='+str(k2)+'/EOCL1dis_k2.txt')
###############################

#data k3
################################
k3=3
EOCL1cont_k3 = np.loadtxt(path_data+'/kmax='+str(k3)+'/EOCL1cont_k3.txt')
EOCL1dis_k3 = np.loadtxt(path_data+'/kmax='+str(k3)+'/EOCL1dis_k3.txt')
###############################

#figure for order of convergence continuous L1 norm
plt.figure(1)
plt.loglog(EOCL1cont_k0[:,0],EOCL1cont_k0[:,1],'o',c='black',label=r'$k=0$')  
plt.loglog(EOCL1cont_k1[:,0],EOCL1cont_k1[:,1],'o',c='C3',label=r'$k=1$')
plt.loglog(EOCL1cont_k2[:,0],EOCL1cont_k2[:,1],'o',c='C1',label=r'$k=2$')
plt.loglog(EOCL1cont_k3[:,0],EOCL1cont_k3[:,1],'o',c='C2',label=r'$k=3$') 

plt.loglog(EOCL1cont_k0[1:,0],1*EOCL1cont_k0[1:,0]**(-1),':',c='black')
plt.text(6,0.2,r'$\propto N_{\mathrm{bins}/\mathrm{dec}}^{-1}$',color='black')

plt.loglog(EOCL1cont_k1[1:,0],0.7*EOCL1cont_k1[1:,0]**(-2),':',c='C3') 
plt.text(6,0.02,r'$\propto N_{\mathrm{bins}/\mathrm{dec}}^{-2}$',color='C3')

plt.loglog(EOCL1cont_k2[1:,0],0.5*EOCL1cont_k2[1:,0]**(-3),':',c='C1')
plt.text(6,0.0025,r'$\propto N_{\mathrm{bins}/\mathrm{dec}}^{-3}$',color='C1')

plt.loglog(EOCL1cont_k3[1:,0],0.1*EOCL1cont_k3[1:,0]**(-4),':',c='C2') 
plt.text(6,1*10**(-4),r'$\propto N_{\mathrm{bins}/\mathrm{dec}}^{-4}$',color='C2')

plt.xlabel(r'$N_{\mathrm{bins}/\mathrm{decade}}$')
plt.ylabel(r'$e_{\mathrm{c},N}$')
plt.xlim(xmax=17)
plt.legend(loc='lower left')
plt.tight_layout()
# plt.savefig(path_plot+'/EOCL1cont_kconst.png',dpi=192,**savefig_options)


#figure for order of convergence discrete L1 norm
plt.figure(2)
plt.loglog(EOCL1dis_k0[:,0],EOCL1dis_k0[:,1],'o',c='black',label=r'$k=0$')  
plt.loglog(EOCL1dis_k1[:,0],EOCL1dis_k1[:,1],'o',c='C3',label=r'$k=1$')
plt.loglog(EOCL1dis_k2[:,0],EOCL1dis_k2[:,1],'o',c='C1',label=r'$k=2$')
plt.loglog(EOCL1dis_k3[:,0],EOCL1dis_k3[:,1],'o',c='C2',label=r'$k=3$') 

plt.loglog(EOCL1dis_k0[1:,0],1.2*EOCL1dis_k0[1:,0]**(-2),':',c='black')
plt.text(6,0.04,r'$\propto N_{\mathrm{bins}/\mathrm{dec}}^{-2}$',color='black')

plt.loglog(EOCL1dis_k1[1:,0],0.4*EOCL1dis_k1[1:,0]**(-2),':',c='C3') 
plt.text(6,0.0015,r'$\propto N_{\mathrm{bins}/\mathrm{dec}}^{-2}$',color='C3')

plt.loglog(EOCL1dis_k2[1:,0],0.5*EOCL1dis_k2[1:,0]**(-4),':',c='C1')
plt.text(6,4*10**(-4),r'$\propto N_{\mathrm{bins}/\mathrm{dec}}^{-4}$',color='C1')

plt.loglog(EOCL1dis_k3[1:,0],0.1*EOCL1dis_k3[1:,0]**(-4),':',c='C2') 
plt.text(3,2*10**(-5),r'$\propto N_{\mathrm{bins}/\mathrm{dec}}^{-4}$',color='C2')

plt.xlabel(r'$N_{\mathrm{bins}/\mathrm{decade}}$')
plt.ylabel(r'$e_{\mathrm{d},N}$')
plt.xlim(xmax=17)
plt.legend(loc='lower left')
plt.tight_layout()
# plt.savefig(path_plot+'/EOCL1dis_kconst.png',dpi=192,**savefig_options)

plt.show()