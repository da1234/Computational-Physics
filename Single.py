
from matplotlib import pyplot as plt 
import numpy as np 



class Single:
    
    
    def __init__(self,method,Egraph,damping,step,tmax):
        
        if method == "Euler":
            self.Euler(damping,step,tmax,Egraph)
            
        elif method == "leap":
            self.leap(damping,step,tmax,Egraph)
            
        elif method =="Implicit":
            self.Euler_backwrds(damping,step,tmax,Egraph)
            
        elif method =="RK4":
            self.RK4(damping,step,tmax,Egraph)
            
        elif method =="Euler_large":
            self.RK4(damping,step,tmax,Egraph)
        
        
        
    

    def Energy_plot(self,u,v,t,fig,Egraph):
        
        energies = []
        
        for i in range(t.size):
            
            energy = 0.5*(u[i]**2 + v[i]**2)
            energies.append(energy)
            
        energyChanges = [0]
        for i in range(1, len(energies)):
            energyChange = energies[i] - energies[i-1]
            energyChanges.append(energyChange)
            
        
        
        if Egraph == "energy_changes":
          
            ax2 = fig.add_subplot(313)
            ax2.plot(t, energyChanges, 'g', label='Energy Changes')
            plt.ylim(0, max(energyChanges))
            plt.title('Energy Change over time')
            plt.xlabel('Time')
            plt.ylabel('Energy Change')
            plt.legend()
            plt.grid() 
        
        
        elif Egraph == "energy_values":
         
            ax2 = fig.add_subplot(313)
            ax2.plot(t, energies, 'g', label='Energy')
            plt.ylim(0, max(energies))
            plt.title('Energy over time')
            plt.xlabel('Time')
            plt.ylabel('Energy')
            plt.legend()
            plt.grid() 
    
        
    def Energy_plot2(self,u,v,t,fig,Egraph):
        
        energies = []
        
        for i in range(t.size):
            
            
            KE = 0.5*u[i]**2
            PE= (1 - np.cos(v[i]))
            energy = KE +PE
            
            energies.append(energy)
            
        energyChanges = [0]
        for i in range(1, len(energies)):
            energyChange = energies[i] - energies[i-1]
            energyChanges.append(energyChange)
            
            
        
        if Egraph == "energy_changes":
          
            ax2 = fig.add_subplot(313)
            ax2.plot(t, energyChanges, 'g', label='Energy Changes')
            plt.ylim(0, max(energyChanges))
            plt.title('Energy Change over time')
            plt.xlabel('Time')
            plt.ylabel('Energy Change')
            plt.legend()
            plt.grid() 
        
        
        elif Egraph == "energy_values":
         
            ax2 = fig.add_subplot(313)
            ax2.plot(t, energies, 'g', label='Energy')
            plt.ylim(0, max(energies))
            plt.title('Energy over time')
            plt.xlabel('Time')
            plt.ylabel('Energy')
            plt.legend()
            plt.grid() 
        
    
    
    def Euler(self,damping,step,time,Egraph):
        
        D = damping	
        h =step
        t= np.arange(0,time,h)
        u = np.full(t.size, 0.01)
        v = np.full(t.size, 0.01)
        
        
        
        for i in range(1, t.size):
            
            v[i] = v[i-1] +u[i-1]*h
            u[i]= u[i-1] - h*(v[i-1] + D*(u[i-1]))
            
        
        fig = plt.figure()
        self.Energy_plot(u,v,t,fig,Egraph)
        ax1 = fig.add_subplot(311)
        ax1.plot(t, v, label='Theta (Numerial Solution)')
        plt.ylim(0, max(v)+0.01)
        plt.title('Explicit Euler (small Angle) for Simple Pendulum with h=' +str(h) + ' and D=' +str(D))
        plt.xlabel('Time')
        plt.ylabel('Angle')
        plt.legend()
        plt.grid()
        plt.show()
    
    
    #implemenmting leap frog 
    
    def leap(self,damping,step,time,Egraph):
        
        D = damping
        h =step
        t= np.arange(0,time,h)
        u = np.full(t.size, 0.01)
        v = np.full(t.size, 0.01)
        
        u[1] = u[0] + h*(-v[0] - D*u[0])
        v[1] = v[0] + (u[0])*h
        
        
        for i in range(1, t.size):
            
            v[i] = v[i-2] +u[i-1]*2*h
            u[i]= u[i-2] - 2*h*(v[i-1] + D*(u[i-1]))
        
        
        fig = plt.figure()
        self.Energy_plot2(u,v,t,fig,Egraph)
        ax1 = fig.add_subplot(311)
        ax1.plot(t, v, label='Theta (Numerial Solution)')
        plt.title('Leapfrog for Simple Pendulum with h=' +str(h) + ' and D=' +str(D))
        plt.xlabel('Time')
        plt.ylabel('Angle')
        plt.legend()
        plt.grid()
        plt.show()
    
        
    ##implementing RK4     
    def RK4(self,damping,step,time,Egraph):
        
        D = damping
        h =step
        t= np.arange(0,time,h)
        u = np.full(t.size, 0.01)
        v = np.full(t.size, 0.01)
        
        for i in range(1, t.size):
            ku1 = -v[i-1] - D * u[i-1]
            kv1 = u[i-1]
            u1 = u[i-1] + ku1 * h/2
            v1 = v[i-1] + kv1 * h/2
            ku2 = -v1 - D * u1
            kv2 = u1
            u2 = u[i-1] + ku2 * h/2
            v2 = v[i-1] + kv2 * h/2
            ku3 = -v2 - D * u2
            kv3 = u2
            u3 = u[i-1] + ku3 * h
            v3 = v[i-1] + kv3 * h
            ku4 = -v3 - D * u3
            kv4 = u3
            
            v[i] = v[i-1] + h/6 * (kv1 + 2*kv2 + 2*kv3 + kv4)
            u[i] = u[i-1] + h/6 * (ku1 + 2*ku2 + 2*ku3 + ku4)
            
        fig = plt.figure()
        self.Energy_plot(u,v,t,fig,Egraph)
        ax1 = fig.add_subplot(311)
        ax1.plot(t, v, label='Theta (Numerial Solution)')
        plt.title('Rk4 for Simple Pendulum with h=' +str(h) + ' and D=' +str(D))
        plt.xlabel('Time')
        plt.ylabel('Angle')
        plt.legend()
        plt.grid()
        plt.show()
        
    
    
            
        
    
    ##  implementing Euler backwards 
    def Euler_backwrds(self,damping,step,time,Egraph):
        
        D = damping
        h =step
        t= np.arange(0,time,h)
        u = np.full(t.size, 0.01)
        v = np.full(t.size, 0.01)
        B = 1/ (-h**2 -(1+ h*D))
        
        
        
        for i in range(1, t.size):
            
            u[i] = B*(h*v[i-1] - u[i-1])
            v[i]= B*(-v[i-1]*(1+h*D)- u[i-1]*h)
        
    
        fig = plt.figure()
        self.Energy_plot(u,v,t,fig,Egraph)
        ax1 = fig.add_subplot(311)
        ax1.plot(t, v, label='Theta (Numerial Solution)')
        plt.title('Implicit Euler for Simple Pendulum with h=' +str(h) + ' and D=' +str(D))
        plt.xlabel('Time')
        plt.ylabel('Angle')
        plt.legend()
        plt.grid()
        plt.show()
    
        
    def Euler_large(self,damping,step,time,Egraph):
    
        D = damping
        h =step
        t= np.arange(0,time,h)
        u= np.zeros(t.size)
        v= np.full(t.size, 0.75*np.pi)
        
        
        
        for i in range(1, t.size):
            
            u[i]= u[i-1] - h*(np.sin(v[i-1]) + D*(u[i-1]))
            v[i] = v[i-1] +u[i-1]*h
        
        
        fig = plt.figure()
        self.Energy_plot2(u,v,t,fig,Egraph)
        ax1 = fig.add_subplot(311)
        ax1.plot(t, v, label='Theta (Numerial Solution)')
        plt.title('Explicit Euler (Large Angle) for Simple Pendulum with h=' +str(h) + ' and D=' +str(D))
        plt.xlabel('Time')
        plt.ylabel('Angle')
        plt.legend()
        plt.grid()
        plt.show()
    

