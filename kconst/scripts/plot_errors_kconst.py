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
errL1cont_k0 = np.loadtxt(path_data+'/kmax='+str(k0)+'/errL1cont.txt')
errL1dis_k0 = np.loadtxt(path_data+'/kmax='+str(k0)+'/errL1dis.txt')
################################


#data k1
################################
k1=1
errL1cont_k1 = np.loadtxt(path_data+'/kmax='+str(k1)+'/errL1cont.txt')
errL1dis_k1 = np.loadtxt(path_data+'/kmax='+str(k1)+'/errL1dis.txt')
################################

#data k2
################################
k2=2
errL1cont_k2 = np.loadtxt(path_data+'/kmax='+str(k2)+'/errL1cont.txt')
errL1dis_k2 = np.loadtxt(path_data+'/kmax='+str(k2)+'/errL1dis.txt')
###############################

#data k3
################################
k3=3
errL1cont_k3 = np.loadtxt(path_data+'/kmax='+str(k3)+'/errL1cont.txt')
errL1dis_k3 = np.loadtxt(path_data+'/kmax='+str(k3)+'/errL1dis.txt')
###############################

xmin = np.min(errL1dis_k0[:,0])
xmax = np.max(errL1dis_k0[:,0])

yminloglog = 10**(-5)
ymaxloglog = 10**(1)

#figure continuous erreur L1 norm 
plt.figure(1)
plt.loglog(errL1cont_k0[:,0],errL1cont_k0[:,1],linestyle='dotted',color='black',label=r'$k=0$')
plt.loglog(errL1cont_k1[:,0],errL1cont_k1[:,1],linestyle='dashdot',color='C3',label=r'$k=1$')
plt.loglog(errL1cont_k2[:,0],errL1cont_k2[:,1],linestyle='dashed',color='C1',label=r'$k=2$')
plt.loglog(errL1cont_k3[:,0],errL1cont_k3[:,1],linestyle='solid',color='C2',label=r'$k=3$')
plt.legend(loc='lower center',ncol=2)
plt.xlabel(r'time $\tau$')
plt.ylabel(r'numerical error $e_{c,N}$');
plt.xlim(xmin,xmax)
plt.ylim(yminloglog,ymaxloglog)
plt.tight_layout()
# plt.savefig(path_plot+'/errL1_cont_kconst.png',dpi=192,**savefig_options)

#figure discrete erreur L1 norm 
plt.figure(2)
plt.loglog(errL1dis_k0[:,0],errL1dis_k0[:,1],linestyle='dotted',color='black',label=r'$k=0$')
plt.loglog(errL1dis_k1[:,0],errL1dis_k1[:,1],linestyle='dashdot',color='C3',label=r'$k=1$')
plt.loglog(errL1dis_k2[:,0],errL1dis_k2[:,1],linestyle='dashed',color='C1',label=r'$k=2$')
plt.loglog(errL1dis_k3[:,0],errL1dis_k3[:,1],linestyle='solid',color='C2',label=r'$k=3$')
plt.legend(loc='lower center',ncol=2)
plt.xlabel(r'time $\tau$')
plt.ylabel(r'numerical error $e_{d,N}$');
plt.xlim(xmin,xmax)
plt.ylim(yminloglog,ymaxloglog)
plt.tight_layout()
# plt.savefig(path_plot+'/errL1_dis_kconst.png',dpi=192,**savefig_options)

plt.show()