from __future__ import print_function, division
import numpy as np
import Euler_Methods as em
 

# orbit - Program to compute the orbit of a comet.
#clear all;  help orbit;  % Clear memory and print header

# Set initial position and velocity of the comet.
# r0 = float(raw_input("Enter initial radial distance (AU): "))
# v0 = float(raw_input("Enter initial tangential velocity (AU/yr): "))
# matrices are ordered in radial distance from the sun (i.e. the sun first ,jupiter last)
# TODO: Read these from file that contains values in kg/m/s units, and translate
# r = np.array([[1.0022,0.],[0.,0.],[1.,0.],[1.0024,0.],[1.52,0.],[5.2,0.]])
# v = np.array([[0.,6.493185],[0.,0.],[0.,6.283185],[0.,6.0],[0.,5.09634],[0.,2.75536]])


# decrease units by ~9 orders of magnitude
au = 1.496e11
earthV = 2.978e4
solarMass = 1.989e30
velConversion = 2.108e-4
names = ['Sun','Venus','Earth','Moon','Mars','Jupiter']
# r = np.array([[au+5e7,0.],[0.,0.],[au,0.],[au+3.633e8,0.],[2.289e11,0.],[7.785e11,0.]])
# v = np.array([[0.,earthV + 6e3],[0.,0.],[0.,earthV],[0.,earthV + 1.022e3],[0.,2.408e4],[0.,1.307e4]])
r = np.array([[0.,0.],[1.0894e11,0.],[au,0.],[au+3.633e8,0.],[2.289e11,0.],[7.785e11,0.]])
v = np.array([[0.,0.],[0.,3.5e4],[0.,earthV],[0.,earthV + 1.022e3],[0.,2.408e4],[0.,1.307e4]])
# convert from m,s to au,year
r = r/au
v = v * velConversion
# 6.493185
# manually create the sun
state = np.array([[ r[0,0], r[0,1], v[0,0], v[0,1] ]])
pCount = 1
totalPlanets = len(r)
while pCount < totalPlanets:
  state = np.append(state,np.array([[r[pCount,0], r[pCount,1], v[pCount,0],v[pCount,1]]]),0)
  pCount += 1

#Set physical parameters (mass, G*M)
GM = 4*np.pi**2      # Grav. const. * Mass of Sun (au^3/yr^2)
# GM = 6.674e-11         # Grav. const. 
# mass in fraction of sun mass
# mass = np.array([3.7e-14,1.,3.e-6,3.7e-8,3.3e-7,9.5e-4])
mass = np.array([1.989e30,4.8676e24,5.972e24,7.348e22,6.418e23,1.899e27])
# conver to solar masses
mass = mass/solarMass
#mass = np.array([(1./3.)*1.e6,1.,316.+(2./3.)])
adaptErr = 10.# Error parameter used by adaptive Runge-Kutta
time = 0.0

#%* Loop over desired number of steps using specified
#%  numerical method.
tau =  float(raw_input("Enter time step (yr): "))
NumericalMethod = 4
maxSteps = int(3./tau)
orbits = 0
maxOrbits = 2
## TODO: Make loop stop at predetermined time, not after running out of moves
print('Loading...')
for istep in range(0,maxSteps):

  # Progress bar
  if(istep%(maxSteps/10) == 0):
    print('|',end="")
  elif(istep%(maxSteps/50) == 0):
    print('-',end="")
  #%* Record position and energy for plotting.
  # Initially set the arrays for the first step
  if istep == 0:
    # record both r components
    rComp = np.array(np.array([r]))
    tplot = time
  else:
    rComp = np.append(rComp,np.array([r]),axis=0)
    tplot = np.append(tplot,time)

  # orbit count break
  # if(istep > 2 and thplot[-2,4] < 0 and thplot[-1,4] >= 0):
    # break  

  [state,mass] = em.EulerCromer(state,mass,time,tau,em.gravmatrk)
  r = np.copy(state[:,0:2])
  v = np.copy(state[:,2:4])
  # print(v)
  # time = time + tau
print('|')
# modelName = 'model_'+str(tau)
# print(rComp[:,0,0])
# np.savetxt(modelName,(rComp[:,0,0],rComp[:,1]))
# np.savetxt('time_'+str(tau),tplot)