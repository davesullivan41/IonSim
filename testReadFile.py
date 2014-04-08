from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
import os.path

satellite = np.loadtxt('satellite.txt')

plt.figure(1)
plt.plot(satellite[:,0],satellite[:,1])


planetsRaw = np.loadtxt('solarSystem.txt')




totalLength = len(planetsRaw)
numPlanets = 6
pLength = totalLength/numPlanets

rComp = np.array([np.array([planetsRaw[0],planetsRaw[pLength],planetsRaw[2*pLength],
	planetsRaw[3*pLength],planetsRaw[4*pLength],planetsRaw[5*pLength]])])

# rComp = np.array([np.array()])
for i in range(1,int(pLength)):
	rComp = np.append(rComp,np.array([np.array([planetsRaw[0+i],planetsRaw[pLength+i],planetsRaw[2*pLength+i],
		planetsRaw[3*pLength+i],planetsRaw[4*pLength+i],planetsRaw[5*pLength+i]])]),axis=0)
print(rComp)
# tempLength = len(planetsRAW)/5
# planets = np.array([np.array([planetsRAW[0:tempLength]])])
# for i in range(1,5):
# 	planets = np.append(planets,np.array([planetsRAW[i*tempLength:(i+1)*tempLength]]),axis=0)
# plt.figure(2)
# plt.plot(planets[0,:,0],planets[0,:,1])
# plt.show()