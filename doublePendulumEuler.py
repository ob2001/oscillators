# -*- coding: utf-8 -*-
"""
Created on Tue May 25 13:36:41 2021

@author: Obutler
"""

"Section 1: Import libraries"
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

"Section 2a: Define classes"
# Define a new Mass class. This allows us to keep all variables associated
# with one object in an easily accessible place. This will also increase readability
# of the code.
# Any quantities specified in the class body are the defaults for objects of the class.
# The values of any particular object of a class can be modified independently of the
# other objects of that class.
class Mass:
    l = 0 # Length of pendulum in metres
    m = 0 # Mass of bob at end of pendulum in kilograms
    
    # For angular measurements, positive values are counter-clockwise
    t = 0 # Angle of pendulum with the vertical in radians
    w = 0 # Angular velocity
    a = 0 # Angular acceleration
    
    # Define an initialization function to simplify the creation of new "Mass" objects
    # __init__() is a special function that is called whenever a new object is instantiated
    def __init__(self, L = 0, M = 0, T = 0, W = 0, A = 0):
        self.l = L
        self.m = M
        self.t = T
        self.w = W
        self.a = A

"Section 2b: Define functions"
# Streak function is used in animation to update plotted points and lines for each frame
# This function is called each frame and updates the position of objects being drawn to
# the screen.
def streak(t, r1, r2):
    # Update animation
    bob1.set_data(r1[0][t], r1[1])
    bob2.set_data(r2[0][t], r2[1][t])
    string1.set_data([0, r1[0][t]], [0, r1[1]])
    string2.set_data([r1[0][t], r2[0][t]], [r1[1], r2[1][t]])
    
    return string1, string2, bob1, bob2 # Returning bobs after strings allows bobs to
                                        # be drawn on top of strings

"Section 3: Main body of code"
t_sim = 1000 # Simulation time
dt = 0.01 # Size of time step in seconds
T = np.arange(0, t_sim, dt) # Time array
tsteps = T.size # Number of time steps

g = 9.81 # Acceleration due to gravity

# Order of initialization variables: Length, Mass, Angle, Velocity, Acceleration
# Initialize kinematic quantities with arrays for use in Eler's method
m1 = Mass(1, 1, np.zeros(tsteps), np.zeros(tsteps), np.zeros(tsteps)) # Initialize first pendulum
m2 = Mass(1, 1, np.zeros(tsteps), np.zeros(tsteps), np.zeros(tsteps)) # Initialize second pendulum

# You may change the first entries in any of the kinematic arrays to change the initial
# conditions of the double pendulum. We reference values within these different arrays with
# <name of mass object>.<kinematic array>[time step]
m1.t[0] = np.pi/2
m2.t[0] = np.pi

# Perform Euler's method on the kinematic arrays with a for loop to simulate the motion
# of the double pendulum.
for i in range(tsteps - 1):
    num = -g*(2*m1.m + m2.m)*np.sin(m1.t[i]) - m2.m*g*np.sin(m1.t[i] - 2*m2.t[i]) - 2*np.sin(m1.t[i] - m2.t[i])*m2.m*((m2.w[i]**2)*m2.l + (m1.w[i]**2)*m1.l*np.cos(m1.t[i] - m2.t[i]))
    denom = m1.l*(2*m1.m + m2.m - m2.m*np.cos(2*m1.t[i] - 2*m2.t[i]))
    m1.a[i + 1] = dt*num/denom
    m1.w[i + 1] = m1.w[i] + m1.a[i + 1] * dt
    m1.t[i + 1] = m1.t[i] + m1.w[i + 1] * dt
    
    num = 2*np.sin(m1.t[i] - m2.t[i])*(((m1.w[i]**2)*m1.l*(m1.m + m2.m)) + g*(m1.m + m2.m)*np.cos(m1.t[i]) + (m2.w[i]**2)*m2.l*m2.m*np.cos(m1.t[i] - m2.t[i]))
    denom = m2.l*(2*m1.m + m2.m - m2.m*np.cos(2*m1.t[i] - 2*m2.t[i]))
    m2.a[i + 1] = dt*num/denom
    m2.w[i + 1] = m2.w[i] + m2.a[i + 1] * dt
    m2.t[i + 1] = m2.t[i] + m2.w[i + 1] * dt

# Trig/geometry to determine the polar position of the second bob given the
# position of the first bob and the second bob's relative position to it.
L = np.sqrt(m1.l**2 + m2.l**2 + 2*m1.l*m2.l*np.cos(m2.t - m1.t))
phi = m1.t - np.arcsin(m2.l*np.sin(m1.t - m2.t)/L)

"Section 4: Plotting and Animation"
fig = plt.figure(figsize = (12, 9))

# Make sure we plot on a polar coordinate system because our kinematic quantities
# are for angles
ax1 = fig.add_subplot(111, projection = 'polar')
ax1.set_ylim(0, m1.l + m2.l)

# Set the 0 degree mark to face downwards (south), rather than to the right
ax1.set_theta_zero_location('S')

# Initialize quantities used in animation
string1, = ax1.plot([0, m1.t[0]], [0, m1.l], 'k') # string 1
string2, = ax1.plot([m1.t[0], phi[0]], [m1.l, L[0]], 'k') # string 2
bob1, = ax1.plot(m1.t[0], m1.l, 'bo', markersize = 12 + m1.m) # bob 1
bob2, = ax1.plot(phi[0], L[0], 'bo', markersize = 12 + m2.m) # bob 2

anispeed = 0 # Time between frames (ms)
framefactor = 10 # Number of frames to skip when animating.
                # Useful when animation is very long/slow

# Animate double pendulum
ani = animation.FuncAnimation(fig, streak, frames = int(np.floor(tsteps/framefactor)), repeat = False,
                              fargs = ([m1.t[::framefactor], m1.l], [phi[::framefactor], L[::framefactor]]),
                              interval = anispeed, blit = True)

plt.show()