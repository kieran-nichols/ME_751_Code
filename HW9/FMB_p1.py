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

n = 500
itr = n + 1

x = np.linspace(0., 1., num=itr)
phi = np.zeros([itr])
A = np.zeros([itr,itr])
b = np.zeros([itr])
dx = x[1]
#phi = np.linspace(0., 1., num=itr) # want to solve for phi using linear equation A*phi = b
# what is A and b?

#s = phi_d + np.matmul(Delta,(np.matmul(u,phi_bar) - Delta_phi_bar))

# Using a 1D Analysis in the x direction
# want to solve for vectore phi

for i in np.arange(1,n,1):

    b[i] = s
    A[i,i-1] = - u/(2*dx) - D/dx**2         
    A[i,i] = 2*D/dx**2
    A[i,i+1] =  u/(2*dx) - D/dx**2 
    
A[0,0] = 1; b[0] = 0 
A[n,n] = - D/dx**2; A[n,n-1] = 2*D/dx**2; A[n,n-2] = - D/dx**2#A[n,n-2] = 1;
b[n] = s 

phi = np.linalg.solve(A,b)

plt.figure(1,figsize=(6, 4),dpi=50)
plt.plot(x,phi, label='phi')


V = np.ones(itr)/n
width = V[1]
edge = np.linspace(0., 1., num=itr) # edges of segments
center = np.linspace(0, 1., num=itr) + 1/2/n # edges of segments
center = center[:-1]
phi1 = np.zeros([n])
A1 = np.zeros([n,n])
#c2f = np.zeros([n,n])
#c2c = np.zeros([n,n])
b1 = np.zeros([n])
#phi = np.linspace(0., 1., num=itr) # want to solve for phi using linear equation A*phi = b

# Using a 1D Analysis in the x direction
# want to solve for vectore phi

for i in np.arange(1,n-1,1):
#    print(i)
    b1[i] = V[i] # magnitude of Vi multiplied by si(t); volume in this case is time invariant
    c2f_forward = edge[i+1] - center[i]
    c2f_backward = center[i] - edge[i]
    c2c_forward =  center[i+1] - center[i]
    c2c_backward = center[i] - center[i-1]
    A1[i,i-1] = (-u*(c2f_backward/c2c_backward) - D/c2c_backward)      
    A1[i,i] = (2*D/c2c_backward) 
    A1[i,i+1] = (u*(c2f_forward/c2c_forward) - D/c2c_forward)
    
A1[0,0] = 1; b1[0] = 0 
A1[n-1,n-1] = - D/c2c_forward; A1[n-1,n-2] = 2*D/c2c_backward; A1[n-1,n-3] = -D/c2c_backward;#A[n,n-2] = 1;
b1[n-1] = V[n-1] 

phi1 = np.linalg.solve(A1,b1)

plt.figure(1,figsize=(6, 4),dpi=50)
plt.plot(center,phi1, label='phi')
#plt.ylim([-2, 2])
plt.title("Phi")
plt.xlabel('midpoint of segment')
plt.ylabel('Phi(x_n)')
plt.legend(['Finite Diff','Finite Volume'])
#plt.show()
