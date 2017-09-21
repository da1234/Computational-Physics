import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import argrelextrema

   
def realSolution(theta0, D, t):
    
    alpha = np.sqrt(1 - (D**2)/4)
    trueValue = np.exp(-D*t/2) * (theta0*np.cos(alpha*t) + ((D*theta0)/(2*alpha))*np.sin(alpha*t))
    trueDerivative = (-D/2) * trueValue + np.exp(-D*t/2) * ((-theta0*alpha)*np.sin(alpha*t) + ((D*theta0)/2)*np.cos(alpha*t))
    trueEnergy = (trueDerivative**2)/2 + (trueValue**2)/2
    
    return trueValue, trueEnergy
    
def stabilityCheck(v, theta, t, trueEnergies, fig):
    
    KEs = np.zeros(t.size)
    PEs = np.zeros(t.size)
    energies = np.zeros(t.size)
    
    if theta[0] < 1:
        for i in range(t.size):
            KE = (v[i]**2)/2.
            PE = (theta[i]**2)/2.
            
            energies[i] = KE + PE
            KEs[i] = KE
            PEs[i] = PE
    else:
        for i in range(t.size):
            KE = (v[i]**2)/2.
            PE = (1 - np.cos(theta[i]))
            
            energies[i] = KE + PE
            KEs[i] = KE
            PEs[i] = PE
       
    peaks = argrelextrema(energies, np.greater)
    
    first20Peaks = [energies[0]]
    
    if len(peaks[0]) < 20:
        for i in range(1, 21):
            first20Peaks.append(energies[i])
    else:
        for i in range(20):
            first20Peaks.append(energies[peaks[0][i]])
        
    highestInitialEnergy = max(first20Peaks)
    
    maximumEnergy = max(energies)
    
    percentage = ((float(maximumEnergy) - float(highestInitialEnergy))/float(highestInitialEnergy)) *100.
    
    if percentage > 1.0:
        print 'This solution is unstable'
        print 'Maximum Energy:', maximumEnergy
        print 'Highest Energy in first 20 Energy Oscillations:', highestInitialEnergy
        print 'Percentage Difference:', percentage
        
    else:
        print 'This solution is stable'
        print 'Maximum Energy:', maximumEnergy
        print 'Highest Energy in first 20 Energy Oscillations:', highestInitialEnergy
        print 'Percentage Difference:', percentage
        
    ax2 = fig.add_subplot(212)
    ax2.plot(t, energies, 'g', label='Total Energy')
    ax2.plot(t, KEs, 'b', label='Kinetic Energy')
    ax2.plot(t, PEs, 'r', label='Potential Energy')
    plt.ylim(0, max(energies)+0.0005)
    #ax2.plot(t, trueEnergies, 'y--', label='Analytic Energy')
    plt.title('Energy over time')
    plt.xlabel('Time')
    plt.ylabel('Energy')
    plt.legend()
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
    
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax1.plot(t, theta, label='Theta (Numerial Solution)')
    ax1.plot(t, trueValues, '--', label='Theta (Analytic Solution)')
    plt.title('Explicit Euler for Simple Pendulum with h=' +str(h) + ' and D=' +str(D))
    plt.xlabel('Time')
    plt.ylabel('Angle')
    plt.legend()
    plt.grid()
    
    stabilityCheck(v, theta, t, trueEnergies, fig)
    
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
        
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax1.plot(t, theta, label='Theta (Numerial Solution)')
    ax1.plot(t, trueValues, '--', label='Theta (Analytic Solution)')
    plt.title('Leapfrog for Simple Pendulum with h=' +str(h) + ' and D=' +str(D))
    plt.xlabel('Time')
    plt.ylabel('Angle')
    plt.legend()
    plt.grid()
    
    stabilityCheck(v, theta, t, trueEnergies, fig)
    
   
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
        
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax1.plot(t, theta, label='Theta (Numerial Solution)')
    ax1.plot(t, trueValues, '--', label='Theta (Analytic Solution)')
    plt.title('RK4 for Simple Pendulum with h=' +str(h) + ' and D=' +str(D))
    plt.xlabel('Time')
    plt.ylabel('Angle')
    plt.legend()
    plt.grid()
    
    stabilityCheck(v, theta, t, trueEnergies, fig)

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
        
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax1.plot(t, theta, label='Theta (Numerial Solution)')
    ax1.plot(t, trueValues, '--', label='Theta (Analytic Solution)')
    plt.title('Implicit Euler for Simple Pendulum with h=' +str(h) + ' and D=' +str(D))
    plt.xlabel('Time')
    plt.ylabel('Angle')
    plt.legend()
    plt.grid()
    
    stabilityCheck(v, theta, t, trueEnergies, fig)

def largeEuler(h, D, tMax):
    
    t = np.arange(0, tMax, h)
    v = np.zeros(t.size)
    theta = np.zeros(t.size)
    theta[0] = 0.75 * np.pi
    
    for i in range(1, t.size):
        v[i] = v[i-1] + (-np.sin(theta[i-1]) - D*v[i-1]) * h
        theta[i] = theta[i-1] + (v[i-1]) * h
    
    trueValues, trueEnergies = realSolution(theta[0], D, t)
    
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax1.plot(t, theta, label='Theta (Numerial Solution)')
    ax1.plot(t, trueValues, '--', label='Theta (Analytic Solution)')
    plt.title('Explicit Euler (Large Angle) for Simple Pendulum with h=' +str(h) + ' and D=' +str(D))
    plt.xlabel('Time')
    plt.ylabel('Angle')
    plt.legend()
    plt.grid()
    plt.show()
    
    stabilityCheck(v, theta, t, trueEnergies, fig)
    
largeEuler(1.0,0.2,150)