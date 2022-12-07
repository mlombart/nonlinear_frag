#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import scipy.special as special
import scipy.integrate as integrate
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset

#options for plots
plt.rcParams["font.size"]= 16
plt.rcParams['lines.linewidth'] = 3
plt.rcParams["legend.columnspacing"] = 0.5

marker_style = dict( marker='o',markersize=12, markerfacecolor='white', linestyle='',markeredgewidth=2)


savefig_options=dict(bbox_inches='tight')



#solution kconst g(x,0)=x exp(-x)
def solkconstDL(x,t):
   f = np.exp(-x)*np.exp(-t) + np.sqrt(2.*t)*np.exp(-t)*integrate.quad(lambda y: (special.iv(1,2.*np.sqrt(2.*t*np.log(y/x))) * np.exp(-y))/(y*np.sqrt(np.log(y/x))),x,np.inf)[0]
   return x*f


nbins=20

path_data = '../data'
path_plot = '../plots'

#data k0
################################

k0=0
massgridk0 = np.loadtxt(path_data+'/kmax='+str(k0)+'/massgrid.txt')
xmeanlogk0 = np.loadtxt(path_data+'/kmax='+str(k0)+'/xmeanlog.txt')


gt0_xmeanlog_k0 = np.genfromtxt(path_data+'/kmax='+str(k0)+'/gt0_xmeanlog.txt')
gtend_xmeanlog_k0 = np.genfromtxt(path_data+'/kmax='+str(k0)+'/gtend_xmeanlog.txt')

timek0 = np.loadtxt(path_data+'/kmax='+str(k0)+'/time.txt')


#data k1
################################
k1=1
massgridk1 = np.loadtxt(path_data+'/kmax='+str(k1)+'/massgrid.txt')
xmeanlogk1 = np.loadtxt(path_data+'/kmax='+str(k1)+'/xmeanlog.txt')


gt0_xmeanlog_k1 = np.genfromtxt(path_data+'/kmax='+str(k1)+'/gt0_xmeanlog.txt')
gtend_xmeanlog_k1 = np.genfromtxt(path_data+'/kmax='+str(k1)+'/gtend_xmeanlog.txt')

timek1 = np.loadtxt(path_data+'/kmax='+str(k1)+'/time.txt')

#data k2
################################
k2=2
massgridk2 = np.loadtxt(path_data+'/kmax='+str(k2)+'/massgrid.txt')
xmeanlogk2 = np.loadtxt(path_data+'/kmax='+str(k2)+'/xmeanlog.txt')


gt0_xmeanlog_k2 = np.genfromtxt(path_data+'/kmax='+str(k2)+'/gt0_xmeanlog.txt')
gtend_xmeanlog_k2 = np.genfromtxt(path_data+'/kmax='+str(k2)+'/gtend_xmeanlog.txt')

timek2 = np.loadtxt(path_data+'/kmax='+str(k2)+'/time.txt')


#data k3
################################
k3=3
massgridk3 = np.loadtxt(path_data+'/kmax='+str(k3)+'/massgrid.txt')
xmeanlogk3 = np.loadtxt(path_data+'/kmax='+str(k3)+'/xmeanlog.txt')


gt0_xmeanlog_k3 = np.genfromtxt(path_data+'/kmax='+str(k3)+'/gt0_xmeanlog.txt')
gtend_xmeanlog_k3 = np.genfromtxt(path_data+'/kmax='+str(k3)+'/gtend_xmeanlog.txt')

timek3 = np.loadtxt(path_data+'/kmax='+str(k3)+'/time.txt')



xmin = massgridk1[0]
xmax = massgridk1[-1]

yminloglog = 10**(-17)
ymaxloglog = 10**(3)

x=np.logspace(np.log10(xmin),np.log10(xmax),num=1000,dtype=np.float64)

#load data for exact solution at tau=tend
solkconstDLtend = [solkconstDL(x[i],timek0[-1]) for i in range(len(x))]
# np.savetxt(path_data+"/data_solkconstDL_tend.txt",solkconstDLtend)
# solkconstDLtend = np.loadtxt(path_data+"/data_solkconstDL_tend.txt")

#figure in log-log scale
fig,ax = plt.subplots()
ax.set_ylim(yminloglog,10**15)
ax.set_xlim(xmin,xmax)
ax.set_xscale('log')
ax.set_yscale('log')
ax.loglog(x,solkconstDLtend,'--',c='C0',label='Analytic')
ax.plot(xmeanlogk0,gtend_xmeanlog_k0,markeredgecolor='black',label=r'$k=0$',**marker_style)
ax.plot(xmeanlogk1,gtend_xmeanlog_k1,markeredgecolor='C3',label=r'$k=1$',**marker_style)
ax.plot(xmeanlogk2,gtend_xmeanlog_k2,markeredgecolor='C1',label=r'$k=2$',**marker_style)
ax.plot(xmeanlogk3,gtend_xmeanlog_k3,markeredgecolor='C2',label=r'$k=3$',**marker_style)


#zoomin around the peak in lin-log scale
axins2 = zoomed_inset_axes(ax, 3, loc=1)
axins2.semilogx(x,solkconstDLtend,'--',c='C0',label='Analytic')
axins2.semilogx(xmeanlogk0,gtend_xmeanlog_k0,markeredgecolor='black',label=r'$k=0$',**marker_style)
axins2.semilogx(xmeanlogk1,gtend_xmeanlog_k1,markeredgecolor='C3',label=r'$k=1$',**marker_style)
axins2.semilogx(xmeanlogk2,gtend_xmeanlog_k2,markeredgecolor='C1',label=r'$k=2$',**marker_style)
axins2.semilogx(xmeanlogk3,gtend_xmeanlog_k3,markeredgecolor='C2',label=r'$k=3$',**marker_style)
index_maxvalue = np.where(gtend_xmeanlog_k1==np.max(gtend_xmeanlog_k1))[0][0]
xlim_r = xmeanlogk0[index_maxvalue]*50
xlim_l = xmeanlogk0[index_maxvalue]/10.
axins2.set_xlim(xlim_l, xlim_r)
ylim_up = gtend_xmeanlog_k1[index_maxvalue]*1.2
ylim_down = gtend_xmeanlog_k1[index_maxvalue]/1000.
axins2.set_ylim(ylim_down, ylim_up)
axins2.yaxis.tick_left()
axins2.tick_params(labelsize=10)
plt.setp(axins2.get_yticklabels(), visible=True)


ax.set_xlabel(r'mass $x$')
ax.set_ylabel(r"mass density $g(x,\tau')$")
ax.set_title(r"$\tau'=%.1f$" %(timek0[-1]))
ax.legend(loc='lower left',ncol=2)
mark_inset(ax, axins2, loc1=2, loc2=4, fc="none", ec="0.5")


# plt.savefig(path_plot+'/kconst_loglog.png',dpi=192,**savefig_options)



plt.show()
