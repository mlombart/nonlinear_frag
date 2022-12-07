#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sympy as sym
import numpy as np
from scipy.special import legendre,iv
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches

#options for plots
plt.rcParams["font.size"]= 16
plt.rcParams['lines.linewidth'] = 3

savefig_options=dict(bbox_inches='tight')


#functions to reconstruct the polynomials approximation
#legendre polynomials
def LegendreP(i,x):
   coeffs=legendre(i)
   res = 0
   for j in range(0,i+1):
      res = res+ coeffs[j]*x**j

   return res

#gtilde
def gtilde(massgrid,massbins,gij,theta,k,j,x):
   xij = 2/(massgrid[j+1]-massgrid[j])*(x-massbins[j])
   res1 = 0 
   res2 = 0
   if k==0:
      res2 = gij[j]
   else:
      for i in range(k+1):
         res1 = res1+ gij[j,i]*LegendreP(i,xij)   

      res2 = theta[j]*(res1-gij[j,0])+gij[j,0]
   return res2

def gk0(gij,j,x):
   xij = x
   return gij[j]*LegendreP(0,xij)


def I(massgrid,j):
   res= np.logspace(np.log10(massgrid[j]),np.log10(massgrid[j+1]),num=100)
   return res


#solution kmul g(x,0)=x exp(-x)
def solkmulDL(x,tau):
   res = x*(1+tau)**2*np.exp(-x*(1+tau))
   return res


marker_style = dict( marker='o',markersize=8, markerfacecolor='white', linestyle='',markeredgewidth=1.2)


path_data = '../data'
path_plot = '../plots'

nbins=20

#data k0
################################
k0=0
massgridk0 = np.loadtxt(path_data+'/kmax='+str(k0)+'/massgrid.txt')
massbinsk0 = np.loadtxt(path_data+'/kmax='+str(k0)+'/massbins.txt')

gij_t0_k0 = np.genfromtxt(path_data+'/kmax='+str(k0)+'/gij_t0.txt')
gij_tend_k0 = np.genfromtxt(path_data+'/kmax='+str(k0)+'/gij_tend.txt')

timek0 = np.loadtxt(path_data+'/kmax='+str(k0)+'/time.txt')


################################


#data k1
################################
k1=1
massgridk1 = np.loadtxt(path_data+'/kmax='+str(k1)+'/massgrid.txt')
massbinsk1 = np.loadtxt(path_data+'/kmax='+str(k1)+'/massbins.txt')

gij_t0_k1 = np.genfromtxt(path_data+'/kmax='+str(k1)+'/gij_t0.txt')
gij_tend_k1 = np.genfromtxt(path_data+'/kmax='+str(k1)+'/gij_tend.txt')
theta_k1 = np.loadtxt(path_data+'/kmax='+str(k1)+'/theta.txt')

timek1 = np.loadtxt(path_data+'/kmax='+str(k1)+'/time.txt')

gij_t0_k1 = np.reshape(gij_t0_k1,(nbins,k1+1))
gij_tend_k1 = np.reshape(gij_tend_k1,(nbins,k1+1))


################################

#data k2
################################
k2=2
massgridk2 = np.loadtxt(path_data+'/kmax='+str(k2)+'/massgrid.txt')
massbinsk2 = np.loadtxt(path_data+'/kmax='+str(k2)+'/massbins.txt')

gij_t0_k2 = np.genfromtxt(path_data+'/kmax='+str(k2)+'/gij_t0.txt')
gij_tend_k2 = np.genfromtxt(path_data+'/kmax='+str(k2)+'/gij_tend.txt')

theta_k2 = np.loadtxt(path_data+'/kmax='+str(k2)+'/theta.txt')

timek2 = np.loadtxt(path_data+'/kmax='+str(k2)+'/time.txt')

gij_t0_k2 = np.reshape(gij_t0_k2,(nbins,k2+1))
gij_tend_k2 = np.reshape(gij_tend_k2,(nbins,k2+1))

###############################

#data k3
################################
k3=3
massgridk3 = np.loadtxt(path_data+'/kmax='+str(k3)+'/massgrid.txt')
massbinsk3 = np.loadtxt(path_data+'/kmax='+str(k3)+'/massbins.txt')

gij_t0_k3 = np.genfromtxt(path_data+'/kmax='+str(k3)+'/gij_t0.txt')
gij_tend_k3 = np.genfromtxt(path_data+'/kmax='+str(k3)+'/gij_tend.txt')

theta_k3 = np.loadtxt(path_data+'/kmax='+str(k3)+'/theta.txt')

timek3 = np.loadtxt(path_data+'/kmax='+str(k3)+'/time.txt')

gij_t0_k3 = np.reshape(gij_t0_k3,(nbins,k3+1))
gij_tend_k3 = np.reshape(gij_tend_k3,(nbins,k3+1))

###############################

xmin = massgridk0[0]
xmax = massgridk0[-1]


ymint0linlog=-0.01
ymaxt0linlog=0.4

ymintendlinlog=-10**(-5)
ymaxtendlinlog=10** (-5)

yminloglog = 10**(-16)
ymaxloglog = 1

x=np.logspace(np.log10(xmin),np.log10(xmax),num=1000)


#grid plot in lin-log scale
fig, axes = plt.subplots(4,2,figsize=(10,12),sharex='col', gridspec_kw={'hspace': 0, 'wspace': 0.05})

#add grey lines to highlight bins
for (m,n), subplot in np.ndenumerate(axes):
   axes[m,n].ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
   axes[m,n].set_xlim(xmin,xmax)
   axes[m,n].autoscale(enable=True, axis="y", tight=False)
   axes[m,n].axvline(massgridk0[0],c='grey',alpha=0.3)
   axes[m,n].axhline(0,c='grey',alpha=0.3,linestyle='--')
   axes[m,1].yaxis.tick_right()
   for j in range(nbins):
      axes[m,n].axvline(massgridk0[j+1],ymin=-0.01,c='grey',alpha=0.3)


#add dashed line for peak of the curve at tau=0 + arrow to highlight the movement of the curve to small masses (fragmentation)
arrow0 = mpatches.FancyArrowPatch((1.1, 150), (0.005, 150),mutation_scale=20,color="C4")
arrow1 = mpatches.FancyArrowPatch((1.1, 150), (0.005, 150),mutation_scale=20,color="C4")
arrow2 = mpatches.FancyArrowPatch((1.1, 150), (0.005, 150),mutation_scale=20,color="C4")
arrow3 = mpatches.FancyArrowPatch((1.1, 150), (0.005, 150),mutation_scale=20,color="C4")

axes[0,1].axvline(1.,c='C4',linestyle='dashdot')
axes[0,1].text(1.5, 140,r'$\tau=0$',rotation='vertical')
axes[0,1].add_patch(arrow0)
axes[1,1].axvline(1.,c='C4',linestyle='dashdot')
axes[1,1].text(1.5, 140,r'$\tau=0$',rotation='vertical')
axes[1,1].add_patch(arrow1)
axes[2,1].axvline(1.,c='C4',linestyle='dashdot')
axes[2,1].text(1.5, 140,r'$\tau=0$',rotation='vertical')
axes[2,1].add_patch(arrow2)
axes[3,1].axvline(1.,c='C4',linestyle='dashdot')
axes[3,1].text(1.5, 140,r'$\tau=0$',rotation='vertical')
axes[3,1].add_patch(arrow3)


#add analytic solution in each plot
axes[0,0].semilogx(x,solkmulDL(x,timek0[0]),'--',c='C0',label='Analytic')
axes[0,1].semilogx(x,solkmulDL(x,timek0[-1]),'--',c='C0',label='Analytic')
axes[1,0].semilogx(x,solkmulDL(x,timek1[0]),'--',c='C0',label='Analytic')
axes[1,1].semilogx(x,solkmulDL(x,timek1[-1]),'--',c='C0',label='Analytic')
axes[2,0].semilogx(x,solkmulDL(x,timek2[0]),'--',c='C0',label='Analytic')
axes[2,1].semilogx(x,solkmulDL(x,timek2[-1]),'--',c='C0',label='Analytic')
axes[3,0].semilogx(x,solkmulDL(x,timek3[0]),'--',c='C0',label='Analytic')
axes[3,1].semilogx(x,solkmulDL(x,timek3[-1]),'--',c='C0',label='Analytic')

#add numerical solution
for j in range(nbins):
   axes[0,0].plot(I(massgridk0,j),gk0(gij_t0_k0,j,I(massgridk0,j)),c='black',label=(r'$p_j(x,\tau)$' if j==0 else '_'))
   axes[0,1].plot(I(massgridk0,j),gk0(gij_tend_k0,j,I(massgridk0,j)),c='black',label=(r'$p_j(x,\tau)$' if j==0 else '_'))
   axes[1,0].plot(I(massgridk1,j),gtilde(massgridk1,massbinsk1,gij_t0_k1,theta_k1[0,:],k1,j,I(massgridk0,j)),c='C3',label=(r'$p_j(x,\tau)$' if j==0 else '_'))
   axes[1,1].plot(I(massgridk1,j),gtilde(massgridk1,massbinsk1,gij_tend_k1,theta_k1[-1,:],k1,j,I(massgridk1,j)),c='C3',label=(r'$p_j(x,\tau)$' if j==0 else '_'))
   axes[2,0].plot(I(massgridk2,j),gtilde(massgridk2,massbinsk2,gij_t0_k2,theta_k2[0,:],k2,j,I(massgridk2,j)),c='C1',label=(r'$p_j(x,\tau)$' if j==0 else '_'))
   axes[2,1].plot(I(massgridk2,j),gtilde(massgridk2,massbinsk2,gij_tend_k2,theta_k2[-1,:],k2,j,I(massgridk2,j)),c='C1',label=(r'$p_j(x,\tau)$' if j==0 else '_'))
   axes[3,0].plot(I(massgridk3,j),gtilde(massgridk3,massbinsk3,gij_t0_k3,theta_k3[0,:],k3,j,I(massgridk3,j)),c='C2',label=(r'$p_j(x,\tau)$' if j==0 else '_'))
   axes[3,1].plot(I(massgridk3,j),gtilde(massgridk3,massbinsk3,gij_tend_k3,theta_k3[-1,:],k3,j,I(massgridk3,j)),c='C2',label=(r'$p_j(x,\tau)$' if j==0 else '_'))
   

axes[0,0].plot([], [], ' ', label=r'$k=0$')
axes[1,0].plot([], [], ' ', label=r'$k=1$')
axes[2,0].plot([], [], ' ', label=r'$k=2$')
axes[3,0].plot([], [], ' ', label=r'$k=3$')

axes[0,0].legend()
axes[1,0].legend()
axes[2,0].legend()
axes[3,0].legend()

axes[0,0].set_title(r'$\tau=%d$' %(timek0[0]))
axes[0,1].set_title(r'$\tau=%d$' %(timek0[-1]))

axes[1,0].yaxis.get_offset_text().set_visible(False)
axes[1,1].yaxis.get_offset_text().set_visible(False)
axes[2,0].yaxis.get_offset_text().set_visible(False)
axes[2,1].yaxis.get_offset_text().set_visible(False)
axes[3,0].yaxis.get_offset_text().set_visible(False)
axes[3,1].yaxis.get_offset_text().set_visible(False)

for j in range(4):
   axes[j,0].set_ylabel(r'mass density $g$')
axes[3,0].set_xlabel(r'mass $x$')
axes[3,1].set_xlabel(r'mass $x$')

# plt.savefig(path_plot+'/kmul_linlog.png',**savefig_options)

plt.show()

