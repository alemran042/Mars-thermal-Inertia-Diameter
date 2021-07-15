#####################################################################################
## Program for estimating particle diameters (grain-size) of surface materials of Mars
## from THEMIS thermal nighttime thermal inertia data at different pressure condition
## Prepared by: Al Emran, Space and Planetary Sciences, U of A, [alemran@uark.edu]
#####################################################################################

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
from matplotlib.ticker import ScalarFormatter, NullFormatter
plt.rc('font', family='serif')
plt.rc('mathtext',fontset='cm')
import pdb
from scipy.special import factorial

###################
def TI2Diameter(mi, mx, Pmin, Pmax, meu = 1250, SH = 800, B= 8.1*10**4, C = 0.0015): 
	
	''''
	#Input parameters:
	mi = minimum thermal inertia (TIU)
	mx = maximum thermal inertia (TIU)
	P(min and max) = pressure range at martian surface (in torr)
	meu = typical particle density (ρ)
	SH = specific heat (usally 800 Jkg-1K-1)
	C and B = constants
	
	#Output parameters:
	DI = an array of diamter
	TI = an array of thermal inertia at different pressure
	'''

	TI = np.arange(mi, mx)
	DI = []
	P = np.arange(Pmin, Pmax, .5)
	for ii in P:
		DI.append(((TI**2/(meu*SH))/(C*ii**.6))**(1/(-.11*math.log10(ii/B))))

	return [DI, TI]

data = TI2Diameter(200, 426, 3, 7)

print(len(data[0])) # This shows the number of arrary for the pressure condition

#### Plot results
fig, ax = plt.subplots()
plt.subplots_adjust(left = 0.17, bottom = 0.15, right = 0.97, top = 0.97)

fnt = 12
l1, = ax.plot(data[0][0], data[1], 'k', label = '$Pressure - 3.0  torr $')
l2, = ax.plot(data[0][1], data[1], 'b', label = '$Pressure - 3.5  torr $')
l3, = ax.plot(data[0][2], data[1], 'r', label = '$Pressure - 4.0  torr $')
l4, = ax.plot(data[0][3], data[1], 'g', label = '$Pressure - 4.5  torr $')
l5, = ax.plot(data[0][4], data[1], 'c', label = '$Pressure - 5.0  torr $')
l6, = ax.plot(data[0][5], data[1], 'm', label = '$Pressure - 5.5  torr $')
l7, = ax.plot(data[0][6], data[1], 'lawngreen', label = '$Pressure - 6.0  torr $')
l8, = ax.plot(data[0][7], data[1], 'y', label = '$Pressure - 6.5  torr $')

ax.grid(True)
plt.tick_params(which='both',axis='both',direction='in',top='False',right='True',\
				labelsize=14)
plt.tick_params(which='major',length=10)
plt.minorticks_on()
plt.tick_params(which='minor',length=5)

ax.legend(handles = [l1, l2, l3, l4, l5, l6, l7, l8], loc = 'best', frameon = True,\
		 fontsize = 'medium')
plt.xlabel( 'Diameter ($\mu m$)',fontsize=fnt)
plt.ylabel( 'Thermal inertia ($TI$)',fontsize=fnt)

#plt.savefig('DI vs TI.eps', format='eps')
plt.savefig('DI vs TI.png')
plt.show()

pdb.set_trace()  ### STOP
