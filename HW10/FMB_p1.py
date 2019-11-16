import numpy as np
import matplotlib.pyplot as plt

u = 1 # u(x,t)
D = 1 # D(x,t)
s = 1 # s(x,t)

# Laplace: u=0 and s!=0
# Diffusion u=0 and s=0

n = 20
itr = n + 1

x = np.linspace(0., 1., num=itr)
phi = np.zeros([itr])
psi1 = np.zeros([itr])
psi2 = np.zeros([itr])
dpsi1 = np.zeros([itr])
dpsi2 = np.zeros([itr])
A = np.zeros([itr,itr])
b = np.zeros([itr])
c = np.zeros([itr])
d = np.zeros([itr])
r = np.zeros([itr])
#M = np.zeros([itr,itr])
C = np.zeros([itr,itr])
#D = np.zeros([itr,itr])
r_full_const = np.zeros([itr])
r_full_var = np.zeros([itr])
dx = x[1]
#phi = np.linspace(0., 1., num=itr) # want to solve for phi using linear equation A*phi = b
# what is A and b?

#s = phi_d + np.matmul(Delta,(np.matmul(u,phi_bar) - Delta_phi_bar))

# Using a 1D Analysis in the x direction
# want to solve for vectore phi

for i in np.arange(1,n,1):
#    print(i)
    
    psi1 = 1 #(x[i+1]-x)/(x[i+1]-x[i])
    psi2 = 1 #(-x[i]+x)/(x[i+1]-x[i])
    dpsi1 = -1/(x[i+1]-x[i])
    dpsi2 = 1/(x[i+1]-x[i]) 
    
    r_full_const[i] = s*psi1 #- psi1*(u*phi[i] - D*dphi1[i])
    r_full_var[i] = s*psi1
    
#    M[1,1] = psi1[i]*psi1[i]
#    M[1,2] = psi1[i]*psi2[i]
#    M[2,1] = psi2[i]*psi1[i]
#    M[2,2] = psi2[i]*psi2[i]
    
    C11 = u*psi1*dpsi1
    C12 = D*psi1*dpsi2
#    C21 = u*psi2*dpsi1
#    C22 = D*psi2*dpsi2
#    
    D11 = dpsi1*D*dpsi1
    D12 = dpsi2*D*dpsi2
#    D[i+1,i] = dpsi2[i+1]*D*dpsi1[i]
#    D[i+1,i+1] = dpsi2[i+1]*D*dpsi2[i+1]

    A[i,i-1] = C11/2 - D11 #+ u*dpsi2       
    A[i,i] = 2*D*D11
    A[i,i+1] =  C12/2 - D12
    
#    A[i,i-1] = u*dpsi1[i]/2 - D*dpsi1[i]**2 #+ u*dpsi2       
#    A[i,i] = 2*D*dpsi1[i]**2
#    A[i,i+1] =  u*dpsi2[i]/2 - D*dpsi2[i]**2
    
#    A[i,i-1] = - u/(2*dx) - D/dx**2         
#    A[i,i] = 2*D/dx**2
#    A[i,i+1] =  u/(2*dx) - D/dx**2 

    
A[0,0] = 1; r_full_const[0] = 0 #A[0,1] = - 2*D/dx**2; A[0,2] = - u/(2*dx) - D/dx**2 ;
#A[n,n] = 3; A[n,n-1] = -4; A[n,n-2] = 1;
A[n,n] = -D/dx**2; A[n,n-1] = 2*D/dx**2; A[n,n-2] = -D/dx**2#A[n,n-2] = 1;
r_full_const[n] = s #2/n**2

phi = np.linalg.solve(A,r_full_const)

plt.figure(1,figsize=(6, 4),dpi=50)
plt.plot(x,phi, label='phi')
plt.title("Phi")
plt.xlabel('x_n')
plt.ylabel('Phi(x_n)')
#plt.legend()
#plt.show()
#
#V = np.ones(itr)/n
#width = V[1]
#edge = np.linspace(0., 1., num=itr) # edges of segments
#center = np.linspace(0, 1., num=itr) + 1/2/n # edges of segments
#center = center[:-1]
#phi1 = np.zeros([n])
#A1 = np.zeros([n,n])
##c2f = np.zeros([n,n])
##c2c = np.zeros([n,n])
#b1 = np.zeros([n])
##phi = np.linspace(0., 1., num=itr) # want to solve for phi using linear equation A*phi = b
#
## Using a 1D Analysis in the x direction
## want to solve for vectore phi
#
#for i in np.arange(1,n-1,1):
##    print(i)
#    b1[i] = V[i] # magnitude of Vi multiplied by si(t); volume in this case is time invariant
#    c2f_forward = edge[i+1] - center[i]
#    c2f_backward = center[i] - edge[i]
#    c2c_forward =  center[i+1] - center[i]
#    c2c_backward = center[i] - center[i-1]
#    A1[i,i-1] = (-u*(c2f_backward/c2c_backward) - D/c2c_backward)      
#    A1[i,i] = (2*D/c2c_backward) 
#    A1[i,i+1] = (u*(c2f_forward/c2c_forward) - D/c2c_forward)
#    
#A1[0,0] = 1; b1[0] = 0 
#A1[n-1,n-1] = - D/c2c_forward; A1[n-1,n-2] = 2*D/c2c_backward; A1[n-1,n-3] = -D/c2c_backward;#A[n,n-2] = 1;
#b1[n-1] = V[n-1] 
#
#phi1 = np.linalg.solve(A1,b1)
#
#plt.figure(1,figsize=(6, 4),dpi=50)
#plt.plot(center,phi1, label='phi')
##plt.ylim([-2, 2])
#plt.title("Phi")
#plt.xlabel('midpoint of segment')
#plt.ylabel('Phi(x_n)')
##plt.legend()
##plt.show()
