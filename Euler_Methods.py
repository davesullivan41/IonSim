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

def EulerCromer(x,mass,t,tau,derivs):
  dx = derivs(x,mass,t,tau)
  x[:,2:4] = x[:,2:4] + dx[:,2:4]*tau
  x[:,0:2] = x[:,0:2] + x[:,2:4]*tau
  return x

def EulerCromerSat(x,planets,mass,planetMass,t,tau,derivs):
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
    [dx,newMass] = derivs(x,planetLoc,newMass,planetMass,time,smallTau)
    # EulerCromer
    x[2:4] += dx[2:4]*smallTau
    x[0:2] += x[2:4]*smallTau
  return x,newMass

def gravSat(s,planets,mass,planetMass,t,tau):
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
  [satAccel,newMass] = shipEngine(s[0:2],s[2:4],mass,tau)
  accel += satAccel
  derivs = np.array([s[2],s[3],accel[0],accel[1]])
  return derivs,newMass

def gravmatrk(s,mass,t,tau):
  #%  Returns right-hand side of Kepler ODE; used by Runge-Kutta routines
  #%  Inputs
  #%    s      State matrix [[r(1) r(2) v(1) v(2)],[r(1) r(2) v(1) v(2)]]
  #     mass   mass vector [m(1),m(2)]
  #%    t      Time (not used)
  #%  Output
  #%    deriv  Derivatives [[dr(1)/dt dr(2)/dt dv(1)/dt dv(2)/dt],[dr(1)/dt dr(2)/dt dv(1)/dt dv(2)/dt]]
  count = 0
  # Loop through each planet
  for body1 in s:
    # r1 = body1[0:2]
    # v1 = body1[2:4]
    accel = np.array([0.0,0.0])
    massCount = 0
    # Possibly keep sun motionless
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

def shipEngine(position,velocity,mass,tau):
  # massLoss = 3.333e-9 # kg/s
  # massLoss = 3.06e-32 # solar masses per second
  # exhaustVelocity = 30000 # m/s
  # exhaustVelocity = 6.324 # au/year
  if mass > 2.51e-31:
    newMass = mass-6.*3.06e-32*tau
    dv = (velocity/np.linalg.norm(velocity))*6.324*np.log(mass/newMass)/tau
    return [dv,newMass]
  else:
    return [[0.,0.],mass]
  
