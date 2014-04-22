from __future__ import print_function, division
import numpy as np

GM = 4*np.pi**2      # Grav. const. * Mass of Sun (au^3/yr^2)

def intrpf(xi,x,y):
  #'''
  #% Function to interpolate between data points
  #% using Lagrange polynomial (quadratic)
  #% Inputs
  #%   x    Vector of x coordinates of data points (3 values)
  #%   y    Vector of y coordinates of data points (3 values)
  #%   xi   The x value where interpolation is computed
  #% Output
  #%   yi   The interpolation polynomial evaluated at xi
  #'''

  # %* Calculate yi = p(xi) using Lagrange polynomial
  yi = (xi-x[1])*(xi-x[2])/((x[0]-x[1])*(x[0]-x[2]))*y[0]\
  + (xi-x[0])*(xi-x[2])/((x[1]-x[0])*(x[1]-x[2]))*y[1] \
  + (xi-x[0])*(xi-x[1])/((x[2]-x[0])*(x[2]-x[1]))*y[2];
  return (yi);

# Euler-Cromer approximation method using derivs to find the new velocity and acceleration
def EulerCromer(x,mass,t,tau,derivs):
  dx = derivs(x,mass,t,tau)
  x[:,2:4] = x[:,2:4] + dx[:,2:4]*tau
  x[:,0:2] = x[:,0:2] + x[:,2:4]*tau
  return x

# Returns the new velocity and acceleration due to the gravitational force of each object in s on each other object in s
def gravMat(s,mass,t,tau):
  count = 0
  # Loop through each planet
  for body1 in s:
    # r1 = body1[0:2]
    # v1 = body1[2:4]
    accel = np.array([0.0,0.0])
    massCount = 0
    if count == 0:
      for body2 in s:
        # r2 = body2[0:2]
        # v2 = body2[2:4]
        # Distance to body2
        diff = np.linalg.norm(body1[0:2]-body2[0:2])
        if(diff != 0):
          # force from planets
          accel += - GM*mass[massCount]*(body1[0:2]-body2[0:2])/diff**3
        massCount += 1
      derivs = np.array([[body1[2], body1[3], accel[0], accel[1]]])
    else:
      for body2 in s:
        # r2 = body2[0:2]
        # v2 = body2[2:4]
        # Distance to body2
        diff = np.linalg.norm(body1[0:2]-body2[0:2])
        if(diff != 0):
          # force form planets
          accel +=  - GM*mass[massCount]*(body1[0:2]-body2[0:2])/diff**3
        massCount += 1
      derivs = np.append(derivs,np.array([[body1[2], body1[3], accel[0],accel[1]]]),0)
    count += 1
  return derivs

# Euler-Cromer approximation method that interpolates new planet coordinates between the coordinates given in planets
# Used to model the spaceship
def EulerCromerSat(x,planets,mass,planetMass,engineCount,t,tau,derivs):
  # Number of steps satellite takes between each planet step
  scale = 100.
  newMass = mass
  # Decrease tau
  smallTau = tau/scale
  # Create scale between t-tau and t
  counter = np.arange(scale) * smallTau + (t-tau)
  # time vector to be used by interp
  interpTime = [t-tau,t,t+tau]
  for time in counter:
    # Interpolate position of each planet at time
    planetLoc = intrpf(time,interpTime,planets)
    # Calculate acceleration of satellite
    [dx,newMass] = derivs(x,planetLoc,newMass,planetMass,engineCount,time,smallTau)
    # EulerCromer
    x[2:4] += dx[2:4]*smallTau
    x[0:2] += x[2:4]*smallTau
  return x,newMass

# Returns the new velocity and acceleration on teh spaceship due to the gravitational force from all the planets
def gravSat(s,planets,mass,planetMass,engineCount,t,tau):
  body = 0
  accel = np.array([0.,0.])	
  # Loop over each planet
  for planet in planets:
    # distance to planet
    diff = np.linalg.norm(s[0:2]-planet)
    if(diff != 0.):
      # acceleration from planet
      accel = accel - GM*planetMass[body]*(s[0:2]-planet)/diff**3
    body += 1
  # acceleration from the satellite's engines
  [satAccel,newMass] = shipEngine(s[0:2],s[2:4],mass,engineCount,tau)
  accel += satAccel
  derivs = np.array([s[2],s[3],accel[0],accel[1]])
  return derivs,newMass

# This is the function that implements simulating an ion engine
def shipEngine(position,velocity,mass,engineCount,tau):
  # massLoss = 3.333e-6 # kg/s
  # massLoss = 1.68e-32 # solar masses per second
  # exhaustVelocity = 30000 # m/s
  # exhaustVelocity = 6.324 # au/year
  # massPerEngine = 100g = 5.03e-32 solar masses
  # engineCount = 12.
  if mass > (2.51e-31+engineCount*5.03e-32):
    newMass = mass-engineCount*1.68e-32*tau
    dv = (velocity/np.linalg.norm(velocity))*6.324*np.log(mass/newMass)/tau
    return [dv,newMass]
  else:
    return [[0.,0.],mass]
  
