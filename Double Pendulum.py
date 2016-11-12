import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import argrelextrema


def RKDouble(h, R, G, tMax):
    
    t = np.arange(0, tMax, h)
    theta = np.zeros(t.size)
    phi = np.zeros(t.size)
    w = np.zeros(t.size)
    v = np.zeros(t.size)
    energies = np.zeros(t.size)
    KEs = np.zeros(t.size)
    PEs = np.zeros(t.size)
    theta[0] = 0.1
    
    for i in range(1, t.size):
        
        kv1 = (R+1) * theta[i-1] - (R+1) * phi[i-1] + G * (1-(R**-1)) * w[i-1] - (G*R**-1) * v[i-1]
        kw1 = -(R+1) * theta[i-1] + R * phi[i-1] - G * w[i-1]
        kphi1 = v[i-1]
        ktheta1 = w[i-1]
        
        v1 = v[i-1] + kv1 * h/2.
        w1 = w[i-1] + kw1 * h/2.
        phi1 = phi[i-1] + kphi1 * h/2.
        theta1 = theta[i-1] + ktheta1 * h/2.
        
        kv2 = (R+1) * theta1 - (R+1) * phi1 + G * (1-(R**-1)) * w1 - (G*R**-1) * v1
        kw2 = -(R+1) * theta1 + R * phi1 - G * w1
        kphi2 = v1
        ktheta2 = w1
        
        v2 = v[i-1] + kv2 * h/2.
        w2 = w[i-1] + kw2 * h/2.
        phi2 = phi[i-1] + kphi2 * h/2.
        theta2 = theta[i-1] + ktheta2 * h/2.
        
        kv3 = (R+1) * theta2 - (R+1) * phi2 + G * (1-(R**-1)) * w2 - (G*R**-1) * v2
        kw3 = -(R+1) * theta2 + R * phi2 - G * w2
        kphi3 = v2
        ktheta3 = w2
        
        v3 = v[i-1] + kv3 * h
        w3 = w[i-1] + kw3 * h
        phi3 = phi[i-1] + kphi3 * h
        theta3 = theta[i-1] + ktheta3 * h
        
        kv4 = (R+1) * theta3 - (R+1) * phi3 + G * (1-(R**-1)) * w3 - (G*R**-1) * v3
        kw4 = -(R+1) * theta3 + R * phi3 - G * w3
        kphi4 = v3
        ktheta4 = w3
        
        v[i] = v[i-1] + h/6. * (kv1 + 2*kv2 + 2*kv3 + kv4)
        w[i] = w[i-1] + h/6. * (kw1 + 2*kw2 + 2*kw3 + kw4)
        phi[i] = phi[i-1] + h/6. * (kphi1 + 2*kphi2 + 2*kphi3 + kphi4)
        theta[i] = theta[i-1] + h/6. * (ktheta1 + 2*ktheta2 + 2*ktheta3 + ktheta4)
    
    for i in range(t.size):
        
        KE = 1/2. * w[i]**2 + 1/2. * R * (w[i]**2 + v[i]**2 + 2*w[i]*v[i])
        PE = (theta[i]**2)/2. + R*(theta[i]**2)/2. + R*(phi[i]**2)/2.
        
        energies[i] = KE + PE
        KEs[i] = KE
        PEs[i] = PE
        
    peaks = argrelextrema(energies, np.greater)
    #print peaks[0]
    
    #for i in range(20):
    #    print i, energies[i]
    
    first20Peaks = [energies[0]]
    
    if len(peaks[0]) < 20:
        for i in range(1, 20):
            first20Peaks.append(energies[i])
    else:
        for i in range(20):
            first20Peaks.append(energies[peaks[0][i]])
            
    #print first20Peaks
        
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
                   
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax1.plot(t, theta, 'b', label='Theta')
    ax1.plot(t, phi, 'r', label='Phi')
    plt.title('RK4 for Double Pendulum with R = ' +str(R) + ' and G = ' +str(G))
    plt.xlabel('Time')
    plt.ylabel('Angle')
    plt.legend()
    plt.grid()
    plt.show()

    ax2 = fig.add_subplot(212)
    ax2.plot(t, energies, 'g', label='Total Energy')
    ax2.plot(t, KEs, 'b', label='Kinetic Energy')
    ax2.plot(t, PEs, 'r', label='Potential Energy')
    plt.ylim(0, max(energies)+0.0005)
    plt.title('Energy of Double Pendulum with R = ' +str(R) + ' and G = ' +str(G))
    plt.xlabel('Time')
    plt.ylabel('Energy')
    plt.legend()
    plt.grid()
    plt.show()