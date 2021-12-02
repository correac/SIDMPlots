###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2020 Camila Correa
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
##############################################################################

import matplotlib

matplotlib.use("Agg")
from pylab import *
import scipy.stats as stat
from scipy.interpolate import interp1d
import h5py
import numpy as np
import glob
import os.path
from scipy.optimize import curve_fit
import matplotlib.ticker as ticker
from matplotlib.ticker import FuncFormatter

def tick_function(X):
    # X : km/s
    G = 4.30091e-3 #pc Msun^-1 (km/s)^2
    rho = 2.7753e11 * 0.6777**2 #Msun/Mpc^3
    rho /= (1e6)**3  #Msun/pc^3
    M = X**2 / G
    M *= (3 / (4. * np.pi* 200 * rho))**(1/3)
    M = M**(3/2)
    M = np.log10(M)
    return ["%.1f" % z for z in M]

def notation(x, pos):
    power = int(np.log10(x))
    mantissa = x/(10**power)
    return r'$%i$' % (int(x))

def momentum_transfer(x,w0):
    f = 2. * w0**4 / x**4
    f *= (2*np.log(1+0.5*x**2/w0**2)-np.log(1+x**2/w0**2))
    return f

def sigma(x,mx,mphi,alpha):
    w0 = 300. * (mphi/10.) * (10./mx)
    sigma0 = 274.85 * (alpha/0.01)**2 * (mx/10.) * (10./mphi)**4
    sigmam = sigma0
    y = momentum_transfer(x,w0)
    y *= sigmam #cm^2/gr
    y *= 2
    return y

def velocity_dependence(x,w0):
    f = 1. / (1. + x**2 / w0**2)
    return f

def sigma2(x,mx,mphi,alpha):
    w0 = 300. * (mphi/10.) * (10./mx)
    sigma0 = 274.85 * (alpha/0.01)**2 * (mx/10.) * (10./mphi)**4
    y = sigma0 * velocity_dependence(x,w0)
    return y


# Plot parameters
params = {
    "font.size": 12,
    "font.family":"Times",
    "text.usetex": True,
    "figure.figsize": (5, 3.5),
    "figure.subplot.left": 0.15,
    "figure.subplot.right": 0.95,
    "figure.subplot.bottom": 0.13,
    "figure.subplot.top": 0.85,
    "figure.subplot.wspace": 0.25,
    "figure.subplot.hspace": 0.25,
    "lines.markersize": 6,
    "lines.linewidth": 1.0,
    "figure.max_open_warning": 0,
}
rcParams.update(params)
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
plt.grid(True)

xrange = np.arange(0,3,0.05)
xrange = 10**xrange

xdata = np.array([10,11,12,15,20,80,90])
ydata = np.array([100,98,94,86,70,3,1])
popt, pcov = curve_fit(sigma, xdata, ydata)
print("SigmaVel100 model")
print("Momentum-transfer")
print('mx',popt[0])
print('mphi',popt[1])
print('alpha',popt[2])
print(sigma(10,*popt),sigma(90,*popt))
print('---')

#ax1.plot(xdata,ydata,'o',color='tab:grey')
label = r'$m_{x}{=}$%.3f'%popt[0]+'GeV, $m_{\phi}{=}$%.3f'%popt[1]+r'MeV, $\alpha{=}$%.3e'%popt[2]
ax1.plot(xrange,sigma(xrange,*popt),'-',color='tab:purple',label=label)

# popt, pcov = curve_fit(sigma2, xdata, ydata)
# print("Velocity-dep")
# print('mx',popt[0])
# print('mphi',popt[1])
# print('alpha',popt[2])
# print('---')
# ax1.plot(xrange,sigma2(xrange,*popt),'--',color='tab:purple')#,label=r'$m_{x}{=}3.391$GeV, $m_{\phi}{=}0.377$MeV, $\alpha{=}1.53e-5$')


xdata = np.array([10,11,12,15,20,80,90])
ydata = np.array([60,60,58,55,45,3,1])
popt, pcov = curve_fit(sigma, xdata, ydata)
print("SigmaVel60 model")
print('mx',popt[0])
print('mphi',popt[1])
print('alpha',popt[2])
print(sigma(10,*popt),sigma(90,*popt))
print('---')

#ax1.plot(xdata,ydata,'o',color='tab:grey')
label = r'$m_{x}{=}$%.3f'%popt[0]+'GeV, $m_{\phi}{=}$%.3f'%popt[1]+r'MeV, $\alpha{=}$%.3e'%popt[2]
ax1.plot(xrange,sigma(xrange,*popt),'-',color='tab:red',label=label)

# popt, pcov = curve_fit(sigma2, xdata, ydata)
# print('mx',popt[0])
# print('mphi',popt[1])
# print('alpha',popt[2])
# print('---')
# ax1.plot(xrange,sigma2(xrange,*popt),'--',color='tab:red')#,label=r'$m_{x}{=}3.604$GeV, $m_{\phi}{=}0.450$MeV, $\alpha{=}1.63e{-}5$')


xdata = np.array([10,11,12,15,20,80,90])
ydata = np.array([45,45,43,40,35,3,1])
popt, pcov = curve_fit(sigma, xdata, ydata)
print("SigmaVel45 model")
print('mx',popt[0])
print('mphi',popt[1])
print('alpha',popt[2])
print(sigma(10,*popt),sigma(90,*popt))
print('---')

#ax1.plot(xdata,ydata,'o',color='tab:grey')
label = r'$m_{x}{=}$%.3f'%popt[0]+'GeV, $m_{\phi}{=}$%.3f'%popt[1]+r'MeV, $\alpha{=}$%.3e'%popt[2]
ax1.plot(xrange,sigma(xrange,*popt),'-',color='tab:orange',label=label)

# popt, pcov = curve_fit(sigma2, xdata, ydata)
# print('mx',popt[0])
# print('mphi',popt[1])
# print('alpha',popt[2])
# print('---')
# ax1.plot(xrange,sigma2(xrange,*popt),'--',color='tab:orange')#,label=r'$m_{x}{=}3.089$GeV, $m_{\phi}{=}0.396$MeV, $\alpha{=}1.17e{-}5$')

xdata = np.array([10,11,12,15,20,80,90])
ydata = np.array([35,35,33,30,25,4,2])
popt, pcov = curve_fit(sigma, xdata, ydata)
print("SigmaVel35 model")
print('mx',popt[0])
print('mphi',popt[1])
print('alpha',popt[2])
print(sigma(10,*popt),sigma(90,*popt))
print('---')

#ax1.plot(xdata,ydata,'o',color='tab:grey')
label = r'$m_{x}{=}$%.3f'%popt[0]+'GeV, $m_{\phi}{=}$%.3f'%popt[1]+r'MeV, $\alpha{=}$%.3e'%popt[2]
ax1.plot(xrange,sigma(xrange,*popt),'-',color='tab:blue',label=label)

# popt, pcov = curve_fit(sigma2, xdata, ydata)
# print('mx',popt[0])
# print('mphi',popt[1])
# print('alpha',popt[2])
# print('---')
# ax1.plot(xrange,sigma2(xrange,*popt),'--',color='tab:blue')#,label=r'$m_{x}{=}3.634$GeV, $m_{\phi}{=}0.316$MeV, $\alpha{=}6.38e{-}6$')

xdata = np.array([10,11,12,15,20,80,90])
ydata = np.array([21,20,19,18,15,3,1])
popt, pcov = curve_fit(sigma, xdata, ydata)
print("SigmaVel20 model")
print('mx',popt[0])
print('mphi',popt[1])
print('alpha',popt[2])
print(sigma(10,*popt),sigma(90,*popt))
print('---')

#ax1.plot(xdata,ydata,'o',color='tab:grey')
label = r'$m_{x}{=}$%.3f'%popt[0]+'GeV, $m_{\phi}{=}$%.3f'%popt[1]+r'MeV, $\alpha{=}$%.3e'%popt[2]
ax1.plot(xrange,sigma(xrange,*popt),'-',color='tab:green',label=label)

# popt, pcov = curve_fit(sigma2, xdata, ydata)
# print('mx',popt[0])
# print('mphi',popt[1])
# print('alpha',popt[2])
# print('---')
# ax1.plot(xrange,sigma2(xrange,*popt),'--',color='tab:green')#,label=r'$m_{x}{=}4.195$GeV, $m_{\phi}{=}0.385$MeV, $\alpha{=}6.68e-6$')

# ax1.plot(xrange,sigma(xrange,10,0.5,6.74e-6),'-',lw=1,color='tab:purple',label=r'$m_{x}{=}10$GeV, $m_{\phi}{=}0.5$MeV, $\alpha{=}6.74e{-}6$')
# ax1.plot(xrange,sigma(xrange,5,0.5,6.74e-6),'-',lw=1,color='tab:green',label=r'$m_{x}{=}5$GeV, $m_{\phi}{=}0.5$MeV, $\alpha{=}6.74e{-}6$')
# ax1.plot(xrange,sigma(xrange,3,0.3,6.74e-6),'-',lw=1,color='tab:red',label=r'$m_{x}{=}3$GeV, $m_{\phi}{=}0.3$MeV, $\alpha{=}6.74e{-}6$')
# ax1.plot(xrange,sigma(xrange,1,0.3,6.74e-6),'-',lw=1,color='tab:blue',label=r'$m_{x}{=}1$GeV, $m_{\phi}{=}0.3$MeV, $\alpha{=}6.74e{-}6$')

ax1.plot(xrange,np.ones(len(xrange)),':',lw=1,color='tab:grey')
ax1.plot(xrange,10*np.ones(len(xrange)),':',lw=1,color='tab:grey')

ax1.legend(loc="lower left",labelspacing=0.2,handlelength=1,handletextpad=0.4,frameon=False,fontsize=10.5)
ax1.set_xlabel('$v$ [km/s]')
ax1.set_ylabel('$\sigma_{T}/m_{x}$ [cm$^{2}$g$^{-1}$]')
ax1.set_xlim([10,200])
ax1.set_ylim([1e-1,140])
ax1.set_xscale('log')
ax1.set_yscale('log')

new_tick_locations = np.array([10,50,100,150,200])
print(tick_function(new_tick_locations))
ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(new_tick_locations)
ax2.set_xticklabels(tick_function(new_tick_locations))
ax2.set_xlabel(r"M$_{200}$ [log$_{10}$ M$_{\odot}$]")
ax1.tick_params(direction='in', axis='both', which='both', pad=4.5)
ax2.tick_params(direction='in', axis='both', which='both', pad=4.5)
plt.savefig('find_cross_section_2_models.png', dpi=200)


