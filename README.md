# Comphys-


This repository contains a plethora of tools for solving partial differential euqations using computational methods.

The Euler forward, Runge-Kutta methods and other Finite Difference Methods are covered, all within the context of solving the single and double pendulum. 



Enjoy!!!

Follow the instructions bellow, you'll need python 2


SinglePendulum 

import Single.py 

make an instance of Single 

takes 5 arguments:

		method : Euler, leap, Implicit, RK4,or Euler_large (as a string)
		
		Egraph: energy_values or energy_changes 

		damping
		
		step size

		max time 


DoublePendulum 

import Double.py 

make an instance of Double 

takes 5 arguments:

		R: ratio of bottom mass to topmost pendulum

		damping
		
		step size 
		
		max time  

		Egraph: energy_values or energy_changes 
