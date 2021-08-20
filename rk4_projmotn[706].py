# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 07:20:43 2021

@author: massam
"""

""" Section 1:  Start by importing relevant libraries
#---------------------------------------------------------------------------"""
import numpy as np
import matplotlib.pyplot as plt


""" Section 2:  Define functions for the script
#---------------------------------------------------------------------------"""

# This function steps through each of the four RK4 steps according to the
# equations of motion (given by f), where q is the set of coords (x,y,vx,vy)
def RK4step(f,q,dt):
    
    k1 = dt*f(q)
    k2 = dt*f(q+0.5*k1)
    k3 = dt*f(q+0.5*k2)
    k4 = dt*f(q+k3)
    
    return q + (k1+2*k2+2*k3+k4)/6

# this function embodies the 2D equations of motion - specifically here it is
# projectile motion with viscous drag, F_v = -cv
def motionEqns(v):
    x,y,vx,vy = v
    return np.array([vx,vy, -c/m * vx, -c/m * vy + g])

""" Section 3:  Main - define variables
#---------------------------------------------------------------------------"""

# input parameters for projectile motion
c = 0.5
m = 2
g = -9.8

# set up time steps
dt = 0.1
T = np.arange(0,10,dt)

# Q0 is a 4-vector set coordinates for initial conditions: x0, y0, vx0, vy0
Q0 = np.array([0, 0, 173, 100])
# Q is the array of coordinates for each time step
Q = np.zeros([T.size,Q0.size])
Q[0,:] = Q0

# for loop to calculate coordinates 
for i in range(T.size-1):
    Q[i+1,:] = RK4step(motionEqns, Q[i], dt)

x,y,vx,vy = Q.T     # the '.T' here transposes Q

fig = plt.figure(figsize=(9,9))
ax1 = fig.add_subplot(2, 3, 1)
ax2 = fig.add_subplot(2, 3, 2)
ax3 = fig.add_subplot(2, 3, 3)
ax4 = fig.add_subplot(2, 3, 4)
ax5 = fig.add_subplot(2, 3, 5)
ax1.plot(T,x)
ax2.plot(T,y)
ax3.plot(T,vx)
ax4.plot(T,vy)
ax5.plot(x,y,'r')

plt.show()