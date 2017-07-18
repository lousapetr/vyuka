#!/usr/bin/env python

from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt


# theta'(t) = omega(t)
# omega'(t) = -b*omega(t) - c*sin(theta(t))

def example_pendulum():
    def pendulum(y, t, b, c):
        theta, omega = y
        dydt = [omega, -b*omega - c*np.sin(theta)]
        return dydt
    
    y0 = [np.pi - 0.00001, 0.0]
    b = 0.25
    c = 5.0
    
    sol = odeint(pendulum, y0, t, args=(b,c))
    
    plt.plot(t, sol[:,0], 'b', label='theta(t)')
    plt.plot(t, sol[:,1], 'g', label='omega(t)')
    plt.legend(loc='best')
    plt.xlabel('t')
    plt.show()
    

def lotka_volterra(y, t, a, b, c, d):
    prey, predators = y
    d_prey = a * prey - b * prey * predators
    d_predators = d * prey * predators - c * predators
    dydt = [ d_prey, d_predators ]
    return dydt

t = np.linspace(0,30,1000)

y0 = [0.9, 0.3]
a, b, c, d = 0.66, 1.33, 1.0, 1.0
sol = odeint(lotka_volterra, y0, t, args=(a,b,c,d))

# plt.plot(t, sol[:,0], 'b', label='prey')
# plt.plot(t, sol[:,1], 'g', label='predators')
plt.plot(sol[:,0], sol[:,1], label='phase-space')
plt.legend(loc='best')
plt.xlabel('t')
plt.show()


# example_pendulum()
