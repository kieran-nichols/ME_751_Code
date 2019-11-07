import numpy as np
import matplotlib.pyplot as plt

phi_d_steady = 0 # at steady state and is a function of space where phi = phi(x)
Delta = 0
u = 1 # u(x,t)
phi_bar = 0
D = 1 # D(x,t)
Delta_phi_bar = 0
s = 1 # s(x,t)
x0 = 0

# Laplace: u=0 and s!=0
# Diffusion u=0 and s=0

n = 200
m = 1000
itr = n + 1
itr_m = m + 1

x = np.linspace(0., 1., num=itr)
time = np.linspace(0., 5., num=itr_m)
dt = time[1]
phi = np.zeros([itr])
#A = np.zeros([itr,itr,itr])
#b = np.zeros([itr,itr])
A = np.zeros([itr,itr])
b = np.zeros([itr])
dx = x[1]
#phi = np.linspace(0., 1., num=itr) # want to solve for phi using linear equation A*phi = b
# what is A and b?

#s = phi_d + np.matmul(Delta,(np.matmul(u,phi_bar) - Delta_phi_bar))
A[0,0] = 1; b[0] = 0 #A[0,1] = - 2*D/dx**2; A[0,2] = - u/(2*dx) - D/dx**2 ;
A[n,n] = - D/dx**2; A[n,n-1] = 2*D/dx**2; A[n,n-2] = - D/dx**2#A[n,n-2] = 1;
b[n] = s 
Phi = np.zeros([itr,itr_m])
#Phi[0,] = 0
#f = 0

# Using a 1D Analysis in the x direction
# want to solve for vectore phi
for j in np.arange(1,m+1,1):
#    f = (Phi[j] - Phi[j-1])/dt
    for i in np.arange(1,n,1):
        b[i] = s + Phi[i,j-1]/dt
        A[i,i-1] = - u/(2*dx) - D/dx**2        
        A[i,i] = 2*D/dx**2 + 1/dt  
        A[i,i+1] =  u/(2*dx) - D/dx**2 
    phi = np.linalg.solve(A,b)
    Phi[:,j] = phi 
    
#A[0,0] = 1; b[0] = 0 #A[0,1] = - 2*D/dx**2; A[0,2] = - u/(2*dx) - D/dx**2 ;
#A[n,n] = - D/dx**2; A[n,n-1] = 2*D/dx**2; A[n,n-2] = - D/dx**2#A[n,n-2] = 1;
#b[n] = s 
#phi = np.linalg.solve(A,b)

#for j in np.arange(1,n,1):
#for i in np.arange(1,n,1):
#    b[i-1,i] = s + 1/dt
#    b[i,i] = s - 1/dt
#    A[:,i,i-1] = - u/(2*dx) - D/dx**2        
#    A[:,i,i] = 2*D/dx**2  
#    A[:,i,i+1] =  u/(2*dx) - D/dx**2 

#    
#A[:,0,0] = 1; b[:,0] = 0 #A[0,1] = - 2*D/dx**2; A[0,2] = - u/(2*dx) - D/dx**2 ;
##A[n,n] = 3; A[n,n-1] = -4; A[n,n-2] = 1;
#A[:,n,n] = - D/dx**2; A[:,n,n-1] = 2*D/dx**2; A[:,n,n-2] = - D/dx**2#A[n,n-2] = 1;
#b[:,n] = s #2/n**2

plt.figure(1,figsize=(6, 4),dpi=50)
plt.plot(x,Phi[:,itr_m-2], label='phi')
#plt.plot(x,Phi, label='phi')
plt.title("Phi")
plt.xlabel('x_n')
plt.ylabel('Phi(x_n)')
#plt.legend()
#plt.show()


#plt.figure(1,figsize=(6, 4),dpi=50)
#plt.plot(center,phi1, label='phi')
#plt.ylim([-2, 2])
#plt.title("Phi")
#plt.xlabel('midpoint of segment')
#plt.ylabel('Phi(x_n)')
#plt.legend()
#plt.show()
