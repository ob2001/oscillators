# -*- coding: utf-8 -*-
"""
Created on Fri May 14 13:35:01 2021

@author: Obutler
"""

"1. Import libraries"
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

"2. Define functions"
# Define a periodic function pushing the oscillator
def push(t, A = 1, w = 1, phi = 0):
    if(w != 0):
        return A*np.cos((w*t) - phi)
    else:
        return 0

# Used in the animation, updates the position of the oscillator and the length of the spring
def streak(t, r1, r2):
    front1.set_data(r1[0][t], r1[1])
    front2.set_data([0, r2[0][t]], [0, r2[1]])
    
    return front2, front1,

"3. Main body of code"
t_sim = 25 # Simulation time in seconds
dt = 0.01 # Time step for Euler's method
t = np.arange(0, t_sim, dt) # Array of times
tsteps = t.size # Number of time steps

m = 1 # Mass on pendulum in kilograms
g = 9.81 # Acceleration due to gravity
L = 1 # Length of pendulum in metres

crit = 2*np.sqrt(g*m/L) # Critical damping value

b = 0 # Damping constant

w0 = np.sqrt(g/L) # Natuaral frequency of undamped oscillator
z = b/crit # Damping ratio
if(np.abs(z) < 1):
    w2 = w0*np.sqrt(1-z**2) # Angular frequency of underdamped oscillator
else:
    w2 = 0

# Print out the natural frequencies and damping ratio of the oscillator
print(f"w0: {w0}")
print(f"z: {z}")
print(f"w2: {w2}")

# Setup arrays of kinematic quantities
A, V, Theta = np.zeros(tsteps), np.zeros(tsteps), np.zeros(tsteps)
A[0], V[0], Theta[0] = 0, 0, np.pi - 0.01

# Perform Euler's method on the kinematic arrays
for i in range(1, tsteps):
    A[i] = -b*V[i - 1]/m - g*np.sin(Theta[i - 1])/L + 0*push(t[i], A = 0, w = w2)
    V[i] = V[i - 1] + dt*A[i]
    Theta[i] = Theta[i - 1] + dt*V[i]

"4. Plotting and animation"
# Create figure for graphing kinematic quantities
fig = plt.figure(figsize = (19, 9))

# Set up axes
ax1 = fig.add_subplot(131, xlabel = 'Time', ylabel = 'Angular Acceleration')
ax2 = fig.add_subplot(132, xlabel = 'Time', ylabel = 'Angular Velocity')
ax3 = fig.add_subplot(133, xlabel = 'Time', ylabel = 'Angle')

# Plot kinematic quantities over time on their own axes
ax1.plot(t, A)
ax2.plot(t, V)
ax3.plot(t, Theta)

animate = True
if animate:
    # Create second figure to display animation
    fig2 = plt.figure(figsize = (12, 9))
    ax4 = fig2.add_subplot(111, polar = True)
    ax4.set_ylim(0, L + 0.1*L)
    ax4.set_theta_zero_location("S")

    # Used in animation
    front1, = ax4.plot(Theta[0], L, 'bo', markersize = 12 + m) # Circle showing position of mass
    front2, = ax4.plot([0, 0], [Theta[0], L], 'k') # Line from origin to mass
    aniSpeed = 1 # ms
    
    ani = animation.FuncAnimation(fig2, streak, frames = tsteps, repeat = False,
                              fargs = ([Theta, L], [Theta, L]), interval = aniSpeed, blit = True)
    
plt.show()