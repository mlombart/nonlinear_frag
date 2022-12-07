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

nbins=20

path_data = '../data'
path_plot = '../plots'


# data k0
################################
k0=0
abserrM1_k0 = np.loadtxt(path_data+'/kmax='+str(k0)+'/abserrM1.txt')
################################


#data k1
################################
k1=1
abserrM1_k1 = np.loadtxt(path_data+'/kmax='+str(k1)+'/abserrM1.txt')
################################

#data k2
################################
k2=2
abserrM1_k2 = np.loadtxt(path_data+'/kmax='+str(k2)+'/abserrM1.txt')
###############################

#data k3
################################
k3=3
abserrM1_k3 = np.loadtxt(path_data+'/kmax='+str(k3)+'/abserrM1.txt')
###############################

xmin = abserrM1_k0[0,0]
xmax = abserrM1_k0[-1,0]

yminloglog = 10**(-8)
ymaxloglog = 10**(-3)

#figure total mass versus time
plt.figure(1)
plt.loglog(abserrM1_k0[:,0],abserrM1_k0[:,1],linestyle='dotted',color='black',label=r'$k=0$')
plt.loglog(abserrM1_k1[:,0],abserrM1_k1[:,1],linestyle='dashdot',color='C3',label=r'$k=1$')
plt.loglog(abserrM1_k2[:,0],abserrM1_k2[:,1],linestyle='dashed',color='C1',label=r'$k=2$')
plt.loglog(abserrM1_k3[:,0],abserrM1_k3[:,1],linestyle='solid',color='C2',label=r'$k=3$')
plt.legend(loc='upper center',ncol=2)
plt.xlabel(r'time $\tau$')
plt.ylabel(r'numerical error $e_{M_1,N}$');
plt.xlim(xmin,xmax)
plt.ylim(yminloglog,ymaxloglog)
plt.tight_layout()

# plt.savefig(path_plot+'/abserrM1_kmul.png',**savefig_options)

plt.show()
