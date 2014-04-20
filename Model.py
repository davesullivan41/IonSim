from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
import Euler_Methods as em
import os.path

def orbitPlot(center,moons,r,figure,names):
  plt.plot(0.,0.,'o')
  for index in moons:
    plt.plot(r[:,index,0]-r[:,center,0],r[:,index,1]-r[:,center,1])

  plt.title(names[center]+' and satellites')

def orbitPlotSat(center,outer):
  plt.plot(0.,0.,'o')
  plt.plot(outer[:,0]-center[:,0],outer[:,1]-center[:,1])
 

# orbit - Program to compute the orbit of a comet.
#clear all;  help orbit;  % Clear memory and print header


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
# manually create the sun
state = np.array([[ r[0,0], r[0,1], v[0,0], v[0,1] ]])
pCount = 1
totalPlanets = len(r)
while pCount < totalPlanets:
  state = np.append(state,np.array([[r[pCount,0], r[pCount,1], v[pCount,0],v[pCount,1]]]),0)
  pCount += 1

#Set physical parameters (mass, G*M)
GM = 4*np.pi**2      # Grav. const. * Mass of Sun (au^3/yr^2)
# mass in fraction of sun mass
mass = np.array([1.989e30,4.8676e24,5.972e24,7.348e22,6.418e23,1.899e27])
# conver to solar masses
mass = mass/solarMass
adaptErr = 1.# Error parameter used by adaptive Runge-Kutta
time = 0.0

####################################################################
## Build the solar system model that the spacecraft can be placed ##
## in. Never recalculate data that already exists                 ##
####################################################################

#%* Loop over desired number of steps using specified
#%  numerical method.
tau =  float(raw_input("Enter time step (yr): "))
duration = float(raw_input("Enter duration (yr): "))
NumericalMethod = 4
maxSteps = int(duration/tau)
# Note -- only use longest data collection period
filename = 'dat/solarSystem_'+str(maxSteps)+'_'+str(tau)+'.txt'
if(os.path.isfile(filename)):
  print('Loading file from data store')
  planetsRaw = np.loadtxt(filename)
  totalLength = len(planetsRaw)
  planetsRaw = planetsRaw.tolist()
  numPlanets = 6
  pLength = int(totalLength/numPlanets)

  # rComp = np.array([np.array([planetsRaw[0],planetsRaw[pLength],planetsRaw[2*pLength],
  #   planetsRaw[3*pLength],planetsRaw[4*pLength],planetsRaw[5*pLength]])])
  rComp = [[planetsRaw[0],planetsRaw[pLength],planetsRaw[2*pLength],
    planetsRaw[3*pLength],planetsRaw[4*pLength],planetsRaw[5*pLength]]]

  # rComp = np.array([np.array()])
  for i in range(1,int(pLength)):
    rComp.append([planetsRaw[0+i],planetsRaw[pLength+i],planetsRaw[2*pLength+i],
      planetsRaw[3*pLength+i],planetsRaw[4*pLength+i],planetsRaw[5*pLength+i]])
  rComp = np.array(rComp)
  # print(rComp[2])
  tplot = np.arange(pLength)*tau+time
else:
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
    #%* Record position for plotting.
    # Initially set the arrays for the first step
    if istep == 0:
      # record both r components
      rComp = np.array(np.array([r]))
      tplot = time
    else:
      rComp = np.append(rComp,np.array([r]),axis=0)
      tplot = np.append(tplot,time)

    ##### This is the line that calculates the new positions of each planet ####
    state = em.EulerCromer(state,mass,time,tau,em.gravmatrk)
    r = np.copy(state[:,0:2])
  print('|')

  with file(filename,'w') as outfile:
    for planet in range(0,totalPlanets):
      outfile.write('# New planet\n')
      np.savetxt(outfile,rComp[:,planet])
print('Planet model built')



####################################################################
## Now that the solar system has been built, place the spacecraft ##
## into the system and have it accelerate away from starting point##
####################################################################
## Probably need to decrease time step to ~ minutes or seconds
r = [(au+4.2164e7)/au,0.]
v = [0.,(earthV-3.075e3)*velConversion]
satMass = 3.

satMass = satMass/solarMass
massPlot = [satMass]
state = np.array([r[0],r[1],v[0],v[1]])
step = 0

filename = 'dat/satellite_'+str(maxSteps)+'_'+str(tau)+'.txt'
massFilename = 'dat/mass_'+str(maxSteps)+'_'+str(tau)+'.txt'
# Check if data for this run has already been calculated
if(os.path.isfile(filename) and os.path.isfile(massFilename)):
  # Load data from appropriate files
  rSComp = np.loadtxt(filename)
  massPlot = np.loadtxt(massFilename)
  print('Satellite path loaded from data store')
else:
  totalSteps = len(tplot)
  print('Calculating satellite path')
  # rSComp = np.array(np.array([r]))
  rSComp = [r]
  vSComp = [v]
  # points to interp over
  planets = np.copy(rComp[step:step+3])
  [state,satMass] = em.EulerCromerSat(state,planets,satMass,mass,tplot[1],tau,em.gravSat)
  r = np.copy(state[0:2])
  # v = np.copy(state[2:4])
  step += 1
  print('|',end='')
  # Time step when the satellite runs out of fuel
  outOfFuel =0
  for time in tplot[2:-1]:
    # rSComp = np.append(rSComp,np.array([r]),axis=0)
    rSComp.append(r.tolist())
    massPlot.append(satMass)
    # vSComp.append(np.copy(state[2:4]).tolist())
    # points to interp over
    planets = np.copy(rComp[step-1:step+2])
    #### This is the line that calculates the new position and mass for the satellite #####
    [state,satMass] = em.EulerCromerSat(state,planets,satMass,mass,tplot[time],tau,em.gravSat)
    r = np.copy(state[0:2])

    if (satMass == massPlot[-1] and outOfFuel == 0):
      outOfFuel = 1
      print('Out of fuel on step '+str(step))

    # v = np.copy(state[2:4])
    # Progress bar
    step += 1
    if(step%(int(totalSteps/10))==0):
      print('|',end='')
    elif(step%(int(totalSteps/50))==0):
      print('-',end='')
  print('|')
  rSComp = np.array(rSComp)
  massPlot = np.array(massPlot)
  np.savetxt(filename,rSComp)
  np.savetxt(massFilename,massPlot)
  print('Satellite path calculated')




# v = np.linalg.norm(state[2:4])/ velConversion

# Calculate Energy
try:
  rNorm = np.sqrt(rSComp[:,0]**2 + rSComp[:,1]**2)
  vSquared = [np.linalg.norm(v)**2]
  for i in range(1,len(rSComp)):
    vSquared.append(((rSComp[i,0]-rSComp[i-1,0])/(tau*velConversion))**2+((rSComp[i,1]-rSComp[i-1,1])/(tau*velConversion))**2)

  kEnergy = 0.5*massPlot*vSquared
  # G*M = 1.327e20 m^3 /s^2
  pEnergy = (1.327e20*massPlot)/(rNorm*au)
  # calculate the time step when the satellite reached escape velocity
  escapeStep = 0
  for i in range(0,len(kEnergy)):
    if kEnergy[i] > pEnergy[i]:
      escapeStep = i
      break
  print('Final velocity of satellite:'+str(round(np.sqrt(vSquared[-1]),2))+'m/s')
  print('Final mass of satellite:'+str(round(massPlot[-1]*solarMass,3))+'kg')
except NameError:
  print('Attempted to use old data file for satellite -- mass values not recorded')


# Calculate when teh satellite crossed jupiter orbit -- for jupiter flyby
jupiterStep = 0
for i in range(0,len(rNorm)):
  if rNorm[i] > rComp[0,-1,0]:
    jupiterStep = i
    break

# print('Step when satellite crossed jupiter orbit :'+str(jupiterStep))
# print('S Coords:'+str(rSComp[jupiterStep,0])+', '+str(rSComp[jupiterStep,1]))
# print('J coords:'+str(rComp[jupiterStep,-1,0])+', '+str(rComp[jupiterStep,-1,1]))


# Plot the whole system
fig = plt.figure(1); plt.clf()  #Clear figure 1 window and bring forward
orbitPlot(0,[1,2,3,4,5],rComp,fig,names)
plt.plot(rSComp[:,0],rSComp[:,1])
plt.plot(rSComp[escapeStep,0],rSComp[escapeStep,1],'o')
# plt.plot(rSComp[jupiterStep,0],rSComp[jupiterStep,1],'o')
plt.xlabel('Distance (AU)')
plt.grid('on')

# The following two plot statements are only used if you want to just observe the satellite
# Satellite trajectory from Earth frame
# fig = plt.figure(2)
# # orbitPlot(2,[3],rComp,fig,names)
# plt.plot(rSComp[:,0] - rComp[1:-1,2,0],rSComp[:,1]-rComp[1:-1,2,1])

# Satellite trajectory from Solar frame
# fig = plt.figure(3)
# plt.plot(rComp[:,2,0],rComp[:,2,1],rSComp[:,0],rSComp[:,1])

# Show the satellites distance from earth
fig = plt.figure(4)
satOrbitR = np.sqrt((rSComp[:,0]-rComp[1:-1,2,0])**2 + (rSComp[:,1]-rComp[1:-1,2,1])**2)
plt.plot(satOrbitR)
plt.title('Satellite distance from Earth')

# Show the energy so we know exactly when the satellite escaped the sun's gravity well
fig = plt.figure(5)
plt.plot(kEnergy)
plt.plot(pEnergy)
plt.legend(['Kinetic Energy','Potential Energy'])

plt.show()