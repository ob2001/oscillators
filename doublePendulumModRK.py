# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 13:22:53 2021

@author: Obutler
"""

"Section 1: Import libraries"
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

"Section 2: Define functions"
# Streak function is used in animation to update plotted points and lines for each frame
# This function is called each frame and updates the position of objects being drawn to
# the screen.
def evolve(t, r1, r2):
    # Update animation
    bob1.set_data(r1[0][t], r1[1])
    bob2.set_data(r2[0][t], r2[1][t])
    string1.set_data([0, r1[0][t]], [0, r1[1]])
    string2.set_data([r1[0][t], r2[0][t]], [r1[1], r2[1][t]])
    
    return string1, string2, bob1, bob2 # Returning bobs after strings allows bobs to
                                        # be drawn on top of strings

# Define acceleration function to calculate the acceleration of the bobs given
# their angular displacement and velocity.
def acc(v):
    t1, t2, w1, w2 = v
    num1 = -g*(2*m1 + m2)*np.sin(t1) - m2*g*np.sin(t1 - 2*t2) - 2*np.sin(t1 - t2)*m2*(w2**2*m2 + w1**2*m1*np.cos(t1 - t2))
    denom1 = m1*(2*m1 + m2 - m2*np.cos(2*t1 - 2*t2))

    num2 = 2*np.sin(t1 - t2)*(w1**2*m1*(m1 + m2) + g*(m1 + m2)*np.cos(t1) + w2**2*m2*m2*np.cos(t1 - t2))
    denom2 = m2*(2*m1 + m2 - m2*np.cos(2*t1 - 2*t2))
    
    return np.array([w1, w2, num1/denom1, num2/denom2])

# Define a function to run the RK4 algorithm to update the positions of the
# bobs at each time step.
def RK4(func, s, dt):
    k1 = dt*func(s) # k1 and all subsequent k values are actually vectors
                    # because func() both takes and returns a vector
    k2 = dt*func(s+0.5*k1)
    k3 = dt*func(s+0.5*k2)
    k4 = dt*func(s+k3)
    
    return s + (k1 + 2*k2 + 2*k3 + k4)/6

"Section 3: Main body of code"
t_sim = 100 # Simulation time
dt = 0.01 # Size of time step in seconds
T = np.arange(0, t_sim, dt) # Time array
tsteps = T.size # Number of time steps

g = 9.81 # Acceleration due to gravity

m1, m2 = 1, 1 # Masses of bobs
l1, l2 = 1, 1 # Lengths of pendula
X = np.zeros([T.size, 4]) # Kinematic array

# Set first row of kinematic array to be initial conditions
X[0, :] = [np.pi, np.pi - 0.01, 0, 0]

# Perform RK4 method on the kinematic arrays with a for loop to simulate the motion
# of the double pendulum.
for i in range(tsteps - 1):
    X[i + 1, :] = RK4(acc, X[i], dt)

# Calculate kinetic, potential, and total energy of double pendulum.
# To be plotted vs time to see whether energy is conserved.
K1 = 0.5*m1*l1**2*X[:, 2]**2
K2 = 0.5*m2*(l1**2*X[:, 2]**2 + l2**2*X[:, 3]**2 + 2*l1*l2*X[:, 2]*X[:, 3]*np.cos(X[:, 0] - X[:, 1]))
K = K1 + K2
V1 = -m1*g*l1*np.cos(X[:, 0])
V2 = -m2*g*(l1*np.cos(X[:, 0]) + l2*np.cos(X[:, 1]))
V = V1 + V2

ETot = K + V

# Trig/geometry to determine the polar position of the second bob given the
# position of the first bob and the second bob's relative position to it.
L = np.sqrt(l1**2 + l2**2 + 2*l1*l2*np.cos(X[:, 1] - X[:, 0]))
phi = X[:, 0] - np.arcsin(l2*np.sin(X[:, 0] - X[:, 1])/L)

"Section 4: Plotting and Animation"
if(True):
    fig = plt.figure(figsize = (12, 9))
    
    # Make sure we plot on a polar coordinate system because our kinematic quantities
    # are for angles.
    # Also set the 0 degree mark to face downwards (south), rather than to the right.
    ax1 = fig.add_subplot(111, projection = 'polar')
    ax1.set_ylim(0, l1 + l2)
    ax1.set_theta_zero_location('S')
    
    # Initialize quantities used in animation
    string1, = ax1.plot([0, X[0][0]], [0, l1], 'k') # string 1
    string2, = ax1.plot([X[0][0], phi[0]], [l1, L[0]], 'k') # string 2
    bob1, = ax1.plot(X[0][0], l1, 'bo', markersize = 10 + 5*m1) # bob 1
    bob2, = ax1.plot(phi[0], L[0], 'bo', markersize = 10 + 5*m2) # bob 2
    
    anispeed = 1 # Time between frames (ms)
    framefactor = 1 # Number of frames to skip when animating.
                    # Useful when animation is very long/slow
    
    # Animate double pendulum
    ani = animation.FuncAnimation(fig, evolve, frames = int(np.floor(tsteps/framefactor)), repeat = False,
                                  fargs = ([X[::framefactor ,0], l1], [phi[::framefactor], L[::framefactor]]),
                                  interval = anispeed, blit = True)

# Plot kinetic, potential, and total energy versus time.
if(False):
    fig2 = plt.figure(figsize = (12, 9))
    ax1 = fig2.add_subplot(111)
    ax1.plot(T, K, label = 'Kinetic energy')
    ax1.plot(T, V, label = 'Potential energy')
    ax1.plot(T, ETot, label = 'Total energy')
    ax1.legend()

plt.show()