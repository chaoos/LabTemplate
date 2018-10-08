# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy import interpolate
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline
from scipy.optimize import fsolve
from decimal import Decimal
import pandas as pd
from uncertainties import ufloat, unumpy
from uncertainties.umath import *
from subprocess import call
import os
#from helpers import *
import helpers as hp

font = {'family' : 'normal',
        'size'   :  16}
matplotlib.rc('font', **font)

# Custom code starts here

# some initial contants with their errors
T_bs = np.array([30.0, 50.0, 80.0, 110.0]) + 273.15 # °K
P_1 = ufloat(200.0, 0.25*200.0) # Ohm
R_1 = ufloat(1200.0, 120.0) # Ohm
R_2 = ufloat(389*10**-3, 0.1*389*10**-3) # Ohm
I_Shunt	= ufloat(20.0, 0.1) # A
V_Shunt = ufloat(50.0*10**-3, 0.1*50.0*10**-3) # V
L_1 = 0.05 # m
L_2 = 0.05 # m
d_1 = 2.03*10**-3 # m
d_2 = 7.04*10**-3 # m
rho_1 = lambda T: np.polyval(np.polyfit(np.array([293.0, 353.0]), np.array([1.68, 2.06])*10**-8, 1), T) # linear polyfit
rho_2 = lambda T: 44.0*10**-8 # Ohm meters
k_1 = lambda T: np.polyval(np.polyfit([250.0, 400.0], [401.0, 391.0], 1), T) # W/mK
k_2 = lambda T: np.polyval(np.polyfit([275.0, 400.0], [21.9, 26.6], 1), T) # W/mK
F_1 = np.pi*(0.5*d_1)**2# m^2
F_2 = np.pi*(0.5*d_2)**2# m^2
deltaI_T = 8.5*10**-6

# linear polyfit of the reference table of the type k thermocouple
C = np.arange(-10, 21, 1) # °C
V = np.array([-0.392, -0.353, -0.314, -0.275, -0.236, -0.197, -0.157, -0.118, -0.079, -0.039, 0.000, 0.039, 0.079, 0.119, 0.158, 0.198, 0.238, 0.277, 0.317, 0.357, 0.397, 0.437, 0.477, 0.517, 0.557, 0.597, 0.637, 0.677, 0.718, 0.758, 0.798]) * 10**-3 # V
coeffs = hp.upolyfit(V, C, 1)
p = lambda x: np.polyval(coeffs, x)
f = lambda x: np.exp(2*x) + np.sin(2*np.pi*x^2) + R_2

# create plots for all T's separately
for T_b in T_bs:

	Temp = str(int(T_b-273.15)) # T in °C

	# fetch the measured data from the csv file
	V_S = hp.fetch('data/example_data_' + Temp + '.csv', 'V_Shunt [V]', 1.0*10**-3)
	I_T = hp.fetch('data/example_data_' + Temp + '.csv', 'I_T [A]', deltaI_T)
	V_p = hp.fetch('data/example_data_' + Temp + '.csv', 'V_p [V]', 1.0*10**-3)

	# do some calculation
	n = int(V_p.size/2) # n=4, half of the number of measurements
	V_T = R_2*I_T
	I = V_S*I_Shunt/V_Shunt
	V_p_mean = 0.5*(np.abs(np.flipud(V_p[0:n])) + np.abs(V_p[n:]))
	I_mean = 0.5*(np.abs(np.flipud(I[0:n])) + np.abs(I[n:]))
	deltaT = p(V_T)
	deltaTminus = np.flipud(deltaT[0:n])
	deltaTplus = deltaT[n:]
	R_tot = L_1*rho_1(T_b)/F_1 + L_2*rho_2(T_b)/F_2

	Pi12V = 0.5*V_p_mean*(deltaTplus - deltaTminus)/(deltaTplus + deltaTminus)
	Pi12I = ((k_1(T_b)*F_1)/L_1 + (k_2(T_b)*F_2)/L_2)*(deltaTplus - deltaTminus)/(2*I_mean)

	### PLOTS ###

	# I - Pi12V
	plt.figure(0)
	plt.errorbar(hp.nominal(I_mean), hp.nominal(Pi12V)*10**3, xerr=hp.stddev(I_mean), yerr=hp.stddev(Pi12V)*10**3, label=r'' + Temp + '°C')
	plt.ylabel(r'$\Pi_{12}^{(V)}$ [mV]')

	# I - Pi12I
	plt.figure(1)
	plt.errorbar(hp.nominal(I_mean), hp.nominal(Pi12I)*10**3, xerr=hp.stddev(I_mean), yerr=hp.stddev(Pi12I)*10**3, label=r'' + Temp + '°C')
	plt.ylabel(r'$\Pi_{12}^{(I)}$ [mV]')

	# I - dT
	plt.figure(2)
	plt.errorbar(hp.nominal(I), hp.nominal(deltaT), xerr=hp.stddev(I), yerr=hp.stddev(deltaT), label=r'' + Temp + '°C')
	plt.ylabel(r'$\Delta T^{\pm}$ [°K]')

	# I - (dT+ - dT-)
	plt.figure(3)
	plt.errorbar(hp.nominal(I_mean), hp.nominal(deltaTplus - deltaTminus), xerr=hp.stddev(I_mean), yerr=hp.stddev(deltaTplus - deltaTminus), label=r'' + Temp + '°C')
	plt.ylabel(r'$\Delta T^{+} - \Delta T^{-}$ [°K]')

	# I - (dT+ + dT-)
	plt.figure(4)
	plt.errorbar(hp.nominal(I_mean), hp.nominal(deltaTplus + deltaTminus), xerr=hp.stddev(I_mean), yerr=hp.stddev(deltaTplus + deltaTminus), label=r'' + Temp + '°C')
	plt.ylabel(r'$\Delta T^{+} + \Delta T^{-}$ [°K]')

	# I - V_p
	plt.figure(5)
	plt.errorbar(hp.nominal(I), hp.nominal(V_p), xerr=hp.stddev(I), yerr=hp.stddev(V_p), label=r'' + Temp + '°C')
	plt.ylabel(r'$V_p$ [V]')

	# attributes of all plots
	for i in [0, 1, 2, 3, 4, 5]:
		plt.figure(i)
		plt.legend(loc="best")
		plt.grid(True)
		plt.xlabel(r'$I$ [A]')

	V_S = np.concatenate((np.array([r"$V_{Shunt} [mV]$"]), V_S*10**3))
	I_T = np.concatenate((np.array([r"$I_T [micro A]$"]), I_T*10**6))
	V_p = np.concatenate((np.array([r"$V_p [mV]$"]), V_p*10**3))
	I   = np.concatenate((np.array([r"$I [A]$"]), I))
	arr = np.array([I, V_S, I_T, V_p]).T
	hp.replace("table"+Temp, arr)

# save the plots im format .fmt
fmt = "png"
plt.figure(0)
plt.savefig("plots/Pi12V_vs_I." + fmt)

plt.figure(1)
plt.savefig("plots/Pi12I_vs_I." + fmt)

plt.figure(2)
plt.savefig("plots/dT_vs_I." + fmt)

plt.figure(3)
plt.savefig("plots/dTp-dTm_vs_I." + fmt)

plt.figure(4)
plt.savefig("plots/dTp+dTm_vs_I." + fmt)

plt.figure(5)
plt.savefig("plots/Vp_vs_I." + fmt)

# replace the markers in main.tex with the formatted variables/fitlines/tables etc...
hp.replace("EQ1", "f(x) = " + hp.fitline([2, 5], 1))
hp.replace("variable1", ufloat(12345678923.3, 345678923))
hp.replace("variable2", ufloat(69.911212, 6.721212))
hp.replace("variable3", ufloat(0.00992123123, 0.00095123123))
hp.replace("variable4", 0.00992123123)
hp.replace("variable5", ufloat(1, 10**-5))
hp.replace("variable6", 0.0)
hp.replace("variable7", 0)
hp.replace("variable8", 23)
hp.replace("variable9", ufloat(10, 0))
hp.replace("variable10", 10**23)
hp.replace("polyfit1", "f(x) = " + hp.fitline(coeffs, 1))
hp.replace("polyfit2", hp.fitline(coeffs, 2))
hp.replace("polyfit3", hp.fitline(coeffs, 4))
hp.replace("polyfit4", hp.fitline(coeffs))
hp.replace("Name", "Roman Gruber")
hp.replace("Experiment", "Peltier Effect")
hp.replace("deltaI_T", deltaI_T)

hp.compile()

#plt.show()