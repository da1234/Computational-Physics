# -*- coding: utf-8 -*-
import numpy as np
from matplotlib import pyplot as plt

   
def euler(h, D, tMin, tMax):
    
    t = np.arange(tMin, tMax, h)
    v = np.full(t.size, 0.01)
    theta = np.full(t.size, 0.01)
    
    for i in range(1, t.size):
        v[i] = v[i-1] + (-theta[i-1] - D*v[i-1]) * h
        theta[i] = theta[i-1] + (v[i-1]) * h
        
    plt.plot(t, theta)
    plt.xlabel('Time')
    plt.ylabel('Angle')
    
def leapfrog(h, D, tMin, tMax):
    
    t = np.arange(tMin, tMax, h)
    v = np.full(t.size, 0.01)
    theta = np.full(t.size, 0.01)
    
    
    v[1] = v[0] + (-theta[0] - D*v[0]) * h
    theta[1] = theta[0] + (v[0]) * h
    
    for i in range(2, t.size):
        v[i] = v[i-2] + 2 * (-theta[i-1] - D*v[i-1]) * h
        theta[i] = theta[i-2] + 2 * (v[i-1]) * h
        
    plt.plot(t, theta)
    plt.xlabel('Time')
    plt.ylabel('Angle')
   
def RK4(h, D, tMin, tMax):
    
    t = np.arange(tMin, tMax, h)
    v = np.full(t.size, 0.01)
    theta = np.full(t.size, 0.01)
    
    for i in range(1, t.size):
        kv1 = -theta[i-1] - D * v[i-1]
        ktheta1 = v[i-1]
        v1 = v[i-1] + kv1 * h/2
        theta1 = theta[i-1] + ktheta1 * h/2
        kv2 = -theta1 - D * v1
        ktheta2 = v1
        v2 = v[i-1] + kv2 * h/2
        theta2 = theta[i-1] + ktheta2 * h/2
        kv3 = -theta2 - D * v2
        ktheta3 = v2
        v3 = v[i-1] + kv3 * h
        theta3 = theta[i-1] + ktheta3 * h
        kv4 = -theta3 - D * v3
        ktheta4 = v3
        
        v[i] = v[i-1] + h/6 * (kv1 + 2*kv2 + 2*kv3 + kv4)
        theta[i] = theta[i-1] + h/6 * (ktheta1 + 2*ktheta2 + 2*ktheta3 + ktheta4)
        
    plt.plot(t, theta)
    plt.xlabel('Time')
    plt.ylabel('Angle')

def implicit(h, D, tMin, tMax):
    
    t = np.arange(tMin, tMax, h)
    v = np.full(t.size, 0.01)
    theta = np.full(t.size, 0.01)
    
    for i in range(1, t.size):
        denominator = 1 / (h**2 + D*h +1)
        v[i] = denominator * (-h*theta[i-1] + v[i-1])
        theta[i] = denominator * ((1 + h*D)*theta[i-1] + h*v[i-1])
        
    plt.plot(t, theta)
    plt.xlabel('Time')
    plt.ylabel('Angle')