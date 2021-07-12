# -*- coding: utf-8 -*-
"""
Created on Wed May 12 15:46:06 2021

@author: Obutler
"""

"1. Import libraries"
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

"2. Define functions"
# Define a periodic function pushing the oscillator
def push(t, A = 1, w = 1, phi = 0):
    result = A*np.cos(t*w - phi)
    return result

# Used in the animation, updates the position of the oscillator and the length of the spring
def streak(t, r1, r2):
    front1.set_data(r1[0][t], r1[1])
    front2.set_data([0, r2[0][t]], [0, 0])

    return front2, front1,

# WIP
# Algorithm to find the amplitude of oscillations after equilibration
def amplitude(position):
    return

"3. Main body of code"
t_sim = 10 # Simulation time in seconds
dt = 0.01 # Time step for Euler's method
t = np.arange(0, t_sim, dt) # Array of times
tsteps = t.size # Number of time steps

m = 1 # Mass on spring
c = 1 # Damping constant
k = 10 # Spring constant
w0 = np.sqrt(k/m) # Natural frequency
w2 = np.sqrt(w0**2 - 2*c**2) # True critical frequency

# Print out the natural frequency of the oscillator
print(f"w0: {w0}")
print(f"w2: {w2}")

# Setup arrays of kinematic quantities
A, V, X = np.zeros(tsteps), np.zeros(tsteps), np.zeros(tsteps)
A[0], V[0], X[0] = 0, 0, 10

# Perform Euler's method on the kinematic arrays for the specified time range
for i in range(1, tsteps):
    A[i] = (-2*c*V[i - 1] - k*X[i - 1] + push(t[i], A = 10, w = w2))/m
    V[i] = V[i - 1] + dt*A[i]
    X[i] = X[i - 1] + dt*V[i]

"4. Plotting and animation"
# Create figure for graphing kinematic quantities
fig = plt.figure(figsize = (19, 9))

# Set up axes
ax1 = fig.add_subplot(131, xlabel = 'Time', ylabel = 'Acceleration')
ax2 = fig.add_subplot(132, xlabel = 'Time', ylabel = 'Velocity')
ax3 = fig.add_subplot(133, xlabel = 'Time', ylabel = 'Position')

# Plot kinematic quantities over time on their own axes
ax1.plot(t, A)
ax2.plot(t, V)
ax3.plot(t, X)

animate = True
if animate:
    # Create second figure to display animation
    fig2 = plt.figure(figsize = (12, 9))
    ax4 = fig2.add_subplot(111, xlim = (-max(np.abs(X)), max(np.abs(X))))

    # Used in animation
    front1, = ax4.plot(X[0], 0, 'bo', markersize = 12) # Circle showing position of mass
    front2, = ax4.plot([0, X[0]], [0, 0]) # Line from equilibrium point to mass

    aniSpeed = 1 # ms
    frameFactor = int((tsteps*aniSpeed)/(1000)) # Forces animation to be 10 seconds

    ani = animation.FuncAnimation(fig2, streak, frames = tsteps, repeat = False,
                              fargs = ([X[::frameFactor], 0], [X[::frameFactor], 0]), interval = aniSpeed, blit = True)


plt.show()
