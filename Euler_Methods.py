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
  # accel  = dx[:,2:4]
  # dv = dx[:,0:2]
  x[:,2:4] = x[:,2:4] + dx[:,2:4]*tau
  x[:,0:2] = x[:,0:2] + x[:,2:4]*tau
  # x = x+dx*tau
  return x

def EulerCromerSat(x,planets,mass,planetMass,t,tau,derivs):
  scale = 100.
  smallTau = tau/scale
  counter = np.arange(scale) * smallTau + (t-tau)
  interpTime = [t-tau,t,t+tau]
  # print(planetMass)
  for time in counter:
    planetLoc = intrpf(time,interpTime,planets)
    # print(planetLoc)
    [dx,mass] = derivs(x,planetLoc,mass,planetMass,time,smallTau)
    # print(dx[2:4])
    x[2:4] += dx[2:4]*smallTau
    x[0:2] += x[2:4]*smallTau
  # accel  = dx[:,2:4]
  # dv = dx[:,0:2]
  # x = x+dx*tau
  return x,mass

def gravSat(s,planets,mass,planetMass,t,tau):
  body = 0
  accel = np.array([0.,0.])	
  # print(planets)
  for planet in planets:
    # print(planet)
    diff = np.linalg.norm(s[0:2]-planet)
    # print(planet,planetMass[body],diff)
    if(diff != 0.):
      accel = accel - GM*planetMass[body]*(s[0:2]-planet)/diff**3
      # print(planetMass[body],GM*planetMass[body]*(s[0:2]-planet)/diff**3)
    body += 1
    # [satAccel,mass] = shipEngine(s[0:2],s[2:4],mass,tau)
    # accel += satAccel
  [satAccel,mass] = shipEngine(s[0:2],s[2:4],mass,tau)
  accel += satAccel
  derivs = np.array([s[2],s[3],accel[0],accel[1]])
  return derivs,mass

def gravmatrk(s,mass,t,tau):
  #%  Returns right-hand side of Kepler ODE; used by Runge-Kutta routines
  #%  Inputs
  #%    s      State matrix [[r(1) r(2) v(1) v(2)],[r(1) r(2) v(1) v(2)]]
  #     mass   mass vector [m(1),m(2)]
  #%    t      Time (not used)
  #%  Output
  #%    deriv  Derivatives [[dr(1)/dt dr(2)/dt dv(1)/dt dv(2)/dt],[dr(1)/dt dr(2)/dt dv(1)/dt dv(2)/dt]]
  count = 0
  shipDeltaV = 0
  for body1 in s:
    r1 = body1[0:2]
    v1 = body1[2:4]
    accel = np.array([0.0,0.0])
    #mass1 = mass[count]
    massCount = 0
    # count == 0 means that this is the ship --- I could probably do eomthing more sophisticated
    # by modifying the state vector...
    if count == 0:
      for body2 in s:
        r2 = body2[0:2]
        v2 = body2[2:4]
        #mass2 = mass[massCount]

        diff = np.linalg.norm(r1-r2)
        if(diff != 0):
          # force from planets
          accel += - GM*mass[massCount]*(r1-r2)/diff**3
          #print(count,massCount-1,accel)
        massCount += 1
      # [shipDeltaV,shipMass] = shipEngine(r1,v1,mass[count],tau,1)
      # mass[0] = shipMass
      derivs = np.array([[v1[0], v1[1], accel[0], accel[1]]])
      # derivs = np.array([[accel[0]*tau, accel[1]*tau, accel[0]+shipDeltaV[0]/tau, accel[1]+shipDeltaV[1]/tau]])
      # derivs = np.array([[accel[0]*tau, accel[1]*tau, accel[0], accel[1]]])
    else:
      for body2 in s:
        r2 = body2[0:2]
        v2 = body2[2:4]
        diff = np.linalg.norm(r1-r2)
        if(diff != 0):
          accel +=  - GM*mass[massCount]*(r1-r2)/diff**3
          #print(count,massCount-1,accel)
        massCount += 1
      derivs = np.append(derivs,np.array([[v1[0], v1[1], accel[0],accel[1]]]),0)
      # derivs = np.append(derivs,np.array([[accel[0]*tau,accel[1]*tau, accel[0],accel[1]]]),0)
    count += 1
  #print(derivs)
  # because the sun is a jerk and needs to learn its place.
  # derivs[1] = np.array([[0.,0.,0.,0.]])
  return derivs

def shipEngine(position,velocity,mass,tau):
  massLossConstant = 0.99 # number chosen at random
  exhaustVelocity = .0001 # also chosen at random
  newMass = mass#*massLossConstant**tau
  dv = 0.01*velocity#np.array([0.,0.])#exhaustVelocity * np.linalg.norm(mass/newMass) * (velocity/np.linalg.norm(velocity))
  return [dv,newMass]
