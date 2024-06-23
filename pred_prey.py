"""
mite model with local hunting on [-1,1]^2
u_t - nu * div(grad(u)) + div(u*grad(w)) = (alpha * w - beta) u
w_t - mu * div(grad(w)) = (gamma - delta * u) w 

at t = 0 
predators are distributed uniformly u = 20 on whole domain
prey are distributed in a Gaussian with peak 20 at (0,0)
"""

from fenics import *
from dolfin import *
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
import sys

set_log_level(50)

tol = 1E-14

argv = sys.argv
folder_name = str(argv[0]) + '_kappa=' + str(argv[1]) + '_nu=' + str(argv[2]) + '_mu=' + str(argv[3])

##create mesh and define function space
nx = ny = 50
mesh = RectangleMesh(Point(-1,-1),Point(1,1), nx, ny)

#discretise time
T = 1.0 # final time
num_steps = 100 # number of time steps
#T = 12.0 # final time
#num_steps = 6144 # number of time steps
dt = T / num_steps # time step size

#setting the coefficients
#mu = 0.000125 #Diffusion coefficient of the prey
mu = float(argv[3])
gamma = 0.8 #Birth rate of the prey
delta = 2 #Predation rate of the prey
#nu = 0.0002 #Diffusion coefficient of the predator
nu = float(argv[2])
alpha = 2.0 #Birth rate of predator due to feeding
beta = 0.8 #Death rate of predator
#kappa = 30.0 #Impact of prey taxis
kappa = float(argv[1])

# Define function space for system of concentrations
P1 = FiniteElement('P', triangle, 1)
P2 = FiniteElement('P', triangle, 2)
element = MixedElement([P1, P2])
V = FunctionSpace(mesh, element)

#Define initial value
uw_n = Expression(('4*exp(-30*pow(x[0]+0.6,2) - 30*pow(x[1]-0.6,2))','2*exp(-9*pow(x[0]+0.4,2) - 9*pow(x[1]+0.5,2)) + 2*exp(-9*pow(x[0]-0.5,2) - 9*pow(x[1]-0.4,2))'), degree=2)
uw_n = interpolate(uw_n, V)

#Define variational problem
uw = Function(V)
v_1, v_2 = TestFunctions(V)

# Split system functions to access components
u, w = split(uw)
u_n, w_n = split(uw_n)

#Define variational problem
F = (u-u_n-dt*alpha*w*u+dt*beta*u)*v_1*dx + dt*nu*dot(grad(u),grad(v_1))*dx\
	- kappa * dt*u*dot(grad(w_n),grad(v_1))*dx\
	+ (w-w_n-dt*gamma*w+dt*delta*u*w)*v_2*dx + dt*mu*dot(grad(w),grad(v_2))*dx 

#create VTK plot file
vtkfile_u = File('./'+ folder_name +'/predator.pvd')
vtkfile_w = File('./'+ folder_name +'/prey.pvd')

_u, _w = uw_n.split()
vtkfile_u << (_u,0)
vtkfile_w << (_w,0)
File('./'+ folder_name +'/pred_prey'+str(0)+'.xml') << uw_n

#Initiate list of density values for u and w
u_tot = []
w_tot = []
u_tot.append(assemble(_u*dx(mesh)))
w_tot.append(assemble(_w*dx(mesh)))

#Time-stepping
t=0
for i in tqdm(range(num_steps)):

	#Update current time
	t += dt
    
	#solve variational problem for time step
	solve(F == 0, uw)
    
	#Update previous solution
	uw_n.assign(uw)
    
	
	_u, _w = uw_n.split(True)

	#Project nodal values on the positive half space
	_u_arr = _u.vector().get_local()
	_w_arr = _w.vector().get_local()
	dof_P_1 = V.sub(0).dim()
	dof_P_2 = V.sub(1).dim()
	for j in range(dof_P_1):
		if _u_arr[j] <= 0 + tol:
			_u_arr[j] = 0 + tol
	for l in range(dof_P_2):
		if _w_arr[l] <= 0 + tol:
			_w_arr[l] = 0 + tol


	_u.vector()[:] = _u_arr
	_w.vector()[:] = _w_arr

	assign(uw_n.sub(0), _u)
	assign(uw_n.sub(1), _w)

	#Save solution to file (VTK)
	_u, _w = uw_n.split()
	if i < 1000 or (i/100).is_integer():
		vtkfile_u << (_u,t)
		vtkfile_w << (_w,t)
		File('./'+ folder_name +'/pred_prey'+str(t)+'.xml') << uw_n


    
	#Evaluate integral of u and w over whole domanin
	u_tot.append(assemble(_u*dx(mesh))) 
	w_tot.append(assemble(_w*dx(mesh)))

#Plot the densities
fig = plt.figure()
time = [i*dt for i in range(len(u_tot))]
plt.plot(time,u_tot,label = 'predator')
plt.plot(time,w_tot, label='prey')
plt.ylabel('Mite densities')
plt.xlabel('Time')
plt.title('Mite densities over time')
plt.legend()
plt.savefig('./'+ folder_name +'/mite_densities_hunt.png')
plt.close(fig) 
