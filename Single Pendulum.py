# -*- coding: utf-8 -*-
import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import argrelextrema

   
def realSolution(theta0, D, t):
    
    alpha = np.sqrt(1 - (D**2)/4)
    trueValue = np.exp(-D*t/2) * (theta0*np.cos(alpha*t) + ((D*theta0)/(2*alpha))*np.sin(alpha*t))
    trueDerivative = (-D/2) * trueValue + np.exp(-D*t/2) * ((-theta0*alpha)*np.sin(alpha*t) + ((D*theta0)/2)*np.cos(alpha*t))
    trueEnergy = (trueDerivative**2)/2 + (trueValue**2)/2
    
    return trueValue, trueEnergy
    
def stabilityCheck(v, theta, t, trueEnergies):
    
    energies = np.zeros(t.size)
    for i in range(t.size):
        currentEnergy = (v[i]**2)/2 + (theta[i]**2)/2
        energies[i] = currentEnergy
        
    energyChanges = [0]
    for i in range(1, len(energies)):
        energyChange = energies[i] - energies[i-1]
        energyChanges.append(energyChange)
       
    peaks = argrelextrema(energies, np.greater)
    
    first5Peaks = [energies[0]]
    
    if len(peaks[0]) < 4:
        pass
    else:
        for i in range(4):
            first5Peaks.append(energies[peaks[0][i]])
        
    highestInitialEnergy = max(first5Peaks)
    
    maximumEnergy = max(energies)
    
    percentage = ((float(maximumEnergy) - float(highestInitialEnergy))/float(highestInitialEnergy)) *100.
    
    if percentage > 0.5:
        print 'This solution is unstable'
        print 'Maximum Energy:', maximumEnergy
        print 'Highest Energy in first 4 Energy Oscillations:', highestInitialEnergy
        print 'Percentage Difference:', percentage
        
    else:
        print 'This solution is stable'
        print 'Maximum Energy:', maximumEnergy
        print 'Highest Energy in first 4 Energy Oscillations:', highestInitialEnergy
        print 'Percentage Difference:', percentage
        
    plt.figure()
    plt.plot(t, energies)
    #plt.plot(t, trueEnergies, '--')
    plt.title('Energy over time')
    plt.xlabel('Time')
    plt.ylabel('Energy')
    plt.grid() 
    
    plt.figure()
    #plt.axes(xlim = 700, ylim=1 )
    plt.plot(t, energyChanges, 'y')
    plt.plot(t, np.zeros(t.size), 'r')
    plt.title('Energy Change per step interval')
    plt.xlabel('Time')
    plt.ylabel('Energy Change')
    plt.grid()

def euler(h, D, tMax):
    
    t = np.arange(0, tMax, h)
    v = np.zeros(t.size)
    theta = np.zeros(t.size)
    theta[0] = 0.1
    
    for i in range(1, t.size):
        v[i] = v[i-1] + (-theta[i-1] - D*v[i-1]) * h
        theta[i] = theta[i-1] + (v[i-1]) * h
        
    trueValues, trueEnergies = realSolution(theta[0], D, t)
    
    stabilityCheck(v, theta, t, trueEnergies)
    
    plt.figure()
    plt.plot(t, theta)
    plt.plot(t, trueValues, '--')
    plt.title('Explicit Euler with h=' +str(h) + ' and D=' +str(D))
    plt.xlabel('Time')
    plt.ylabel('Angle')
    plt.grid()
    
def leapfrog(h, D, tMax):
    
    t = np.arange(0, tMax, h)
    v = np.zeros(t.size)
    theta = np.zeros(t.size)
    theta[0] = 0.1
    
    v[1] = v[0] + (-theta[0] - D*v[0]) * h
    theta[1] = theta[0] + (v[0]) * h
    
    for i in range(2, t.size):
        v[i] = v[i-2] + 2 * (-theta[i-1] - D*v[i-1]) * h
        theta[i] = theta[i-2] + 2 * (v[i-1]) * h
        
    trueValues, trueEnergies = realSolution(theta[0], D, t)
    
    stabilityCheck(v, theta, t, trueEnergies)
        
    plt.figure()
    plt.plot(t, theta)
    plt.plot(t, trueValues, '--')
    plt.title('Leapfrog with h=' +str(h) + ' and D=' +str(D))
    plt.xlabel('Time')
    plt.ylabel('Angle')
    plt.grid()
    
   
def RK4(h, D, tMax):
    
    t = np.arange(0, tMax, h)
    v = np.zeros(t.size)
    theta = np.zeros(t.size)
    theta[0] = 0.1
    
    for i in range(1, t.size):
        kv1 = -theta[i-1] - D * v[i-1]
        ktheta1 = v[i-1]
        v1 = v[i-1] + kv1 * h/2.
        theta1 = theta[i-1] + ktheta1 * h/2.
        kv2 = -theta1 - D * v1
        ktheta2 = v1
        v2 = v[i-1] + kv2 * h/2.
        theta2 = theta[i-1] + ktheta2 * h/2.
        kv3 = -theta2 - D * v2
        ktheta3 = v2
        v3 = v[i-1] + kv3 * h
        theta3 = theta[i-1] + ktheta3 * h
        kv4 = -theta3 - D * v3
        ktheta4 = v3
        
        v[i] = v[i-1] + h/6. * (kv1 + 2*kv2 + 2*kv3 + kv4)
        theta[i] = theta[i-1] + h/6. * (ktheta1 + 2*ktheta2 + 2*ktheta3 + ktheta4)
        
    trueValues, trueEnergies = realSolution(theta[0], D, t)
    
    stabilityCheck(v, theta, t, trueEnergies)
        
    plt.figure()
    plt.plot(t, theta)
    #plt.plot(t, trueValues, '--')
    plt.title('RK4 with h=' +str(h) + ' and D=' +str(D))
    plt.xlabel('Time')
    plt.ylabel('Angle')
    plt.grid()

def implicit(h, D, tMax):
    
    t = np.arange(0, tMax, h)
    v = np.zeros(t.size)
    theta = np.zeros(t.size)
    theta[0] = 0.1
    
    denominator = 1 / (h**2 + D*h +1)
    for i in range(1, t.size):
        v[i] = denominator * (-h*theta[i-1] + v[i-1])
        #theta[i] = denominator * ((1 + h*D)*theta[i-1] + h*v[i-1])
        theta[i] = theta[i-1] + v[i]*h
        
    trueValues, trueEnergies = realSolution(theta[0], D, t)
    
    stabilityCheck(v, theta, t, trueEnergies)
        
    plt.figure()
    plt.plot(t, theta)
    plt.plot(t, trueValues, '--')
    plt.title('Implicit Euler with h=' +str(h) + ' and D=' +str(D))
    plt.xlabel('Time')
    plt.ylabel('Angle')
    plt.grid()
    #
def largeEuler(h, D, tMax):
    
    t = np.arange(0, tMax, h)
    v = np.zeros(t.size)
    theta = np.zeros(t.size)
    theta[0] = 0.75 * np.pi
    
    for i in range(1, t.size):
        v[i] = v[i-1] + (-np.sin(theta[i-1]) - D*v[i-1]) * h
        theta[i] = theta[i-1] + (v[i-1]) * h
        
    #trueValues, trueEnergies = realSolution(theta[0], D, t)
    #
    #stabilityCheck(v, theta, t, trueEnergies)
    
    plt.figure()
    plt.plot(t, theta)
    #plt.plot(t, trueValues, '--')
    plt.title('Explicit Euler (Large Angle) with h=' +str(h) + ' and D=' +str(D))
    plt.xlabel('Time')
    plt.ylabel('Angle')
    plt.grid()
    
    energies = []
    for i in range(t.size):
        currentEnergy = (v[i]**2)/2 + (1 - np.cos(theta[i]))
        energies.append(currentEnergy)
        
    energyChanges = [0]
    for i in range(1, len(energies)):
        energyChange = energies[i] - energies[i-1]
        energyChanges.append(energyChange)
        
    plt.figure()
    plt.plot(t, energies)
    #plt.plot(t, trueEnergies, '--')
    plt.title('Energy over time')
    plt.xlabel('Time')
    plt.ylabel('Energy')
    plt.grid() 
    
    plt.figure()
    #plt.axes(xlim = 700, ylim=1 )
    plt.plot(t, energyChanges, 'y')
    plt.plot(t, np.zeros(t.size), 'r')
    plt.title('Energy Change per step interval')
    plt.xlabel('Time')
    plt.ylabel('Energy Change')
    plt.grid()
    
    
#def RK4matrix(h, D, tMax):
#    
#    t = np.arange(0, tMax, h)
#    v = np.zeros(t.size)
#    theta = np.full(t.size, 0.01)
#    
#    L = np.matrix([[0, 1], [-1, -D]])
#    
#    for i in range(1, t.size):
#        vector1 = np.matrix([[theta[i-1]], [v[i-1]]])
#        k1 = h * np.dot(L, vector1)
#        k2 = h * np.dot(L, vector1 + 0.5*h*k1)
#        k3 = h * np.dot(L, vector1 + 0.5*h*k2)
#        k4 = h * np.dot(L, vector1 + h*k3)
#        
#        vector2 = vector1 + 1/6 * (k1 + 2*k2 + 2*k3 + k4)
#        theta[i] = vector2[0]
#        v[i] = vector2[1]
#        
#    plt.figure()
#    plt.plot(t, theta)
#    #plt.plot(t, trueValues, '--')
#    plt.title('RK4 MATRIX')
#    plt.xlabel('Time')
#    plt.ylabel('Angle')
#    plt.grid()