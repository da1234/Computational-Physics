#using Euler forward 


from matplotlib import pyplot as plt 
import numpy as np 

def Euler():
    D = 0.0
    h =0.05
    t= np.arange(0,50,h)
    u= np.zeros(t.size)
    v= np.ones(t.size)
    
    
    
    for i in range(1, t.size):
        
        v[i] = v[i-1] +u[i-1]*h
        u[i]= u[i-1] - h*(v[i-1] + D*(u[i-1]))
    
    
    plt.plot(t,v,'.-')
    plt.show()


#implemenmting leap frog 

def leap():
    
    D = 0.0
    h =0.05
    t= np.arange(0,50,h)
    u= np.zeros(t.size)
    v= np.ones(t.size)
    
    
    
    for i in range(1, t.size):
        
        v[i] = v[i-2] +u[i-1]*2*h
        u[i]= u[i-2] - 2*h*(v[i-1] + D*(u[i-1]))
    
    
    plt.plot(t,v,'.-')
    plt.show()
    
##implementing RK4 

def RK4():
    
    D = 0.0
    h =0.05
    t= np.arange(0,50,h)
    u= np.zeros(t.size)
    v= np.ones(t.size)
    
    
    for i in range(1, t.size):
        
        fv1 = u[i-1]
        fu1 = -(v[i-1] + D*u[i-1])
        v1 = v[i-1] + fv1*h/2
        u1 = u[i-1] - fu1*h/2
        fv2 = u1 
        fu2 =-(v1 + D*u1)
        v2 = v[i-1] + fv2*h/2
        u2 = u[i-1] + fu2*h/2
        fv3 = u2
        fu3 = -(v2 + D*u2)
        v3 = v[i-1] + fv3*h
        u3 = u[i-1] + fu3*h
        fv4 = u3 
        fu4 = -(v3 + D*u3)
        
        v[i] = v[i-1] + (fv1 + 2*(fv2 +fv3)+ fv4) /6. * h 
        u[i] = u[i-1] + (fu1 + 2*(fu2 +fu3)+ fu4) /6. * h 
        
        
        
        
        
        
    
    plt.plot(t,v,'.-')
    plt.show()    

##  implementing Euler backwards 
def Euler_backwrds():
    D = 0.2
    h =0.05
    t= np.arange(0,100,h)
    u= np.zeros(t.size)
    v= np.ones(t.size)
    B = 1/ (-h**2 -(1+ h*D))
    
    
    
    for i in range(1, t.size):
        
        u[i] = B*(h*v[i-1] - u[i-1])
        v[i]= B*(-v[i-1]*(1+h*D)- u[i-1]*h)
    
    
    plt.plot(t,v,'.-')
    plt.show()   
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
     
    
Euler_backwrds()