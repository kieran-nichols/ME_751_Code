import numpy as np
import matplotlib.pyplot as plt

phi_d_steady = 0 # at steady state and is a function of space where phi = phi(x)
u = 1 # u(x,t)
D = 1 # D(x,t)
s = 1 # s(x,t)
x0 = 0

# Laplace: u=0 and s!=0
# Diffusion u=0 and s=0

n = 20
itr = n + 1
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
    b1[i] = V[i]/width # magnitude of Vi multiplied by si(t); volume in this case is time invariant
#    c2f_forward = edge[i+1] - center[i]
#    c2f_backward = center[i] - edge[i]
    c2c_forward =  center[i+1] - center[i]
    c2c_backward = center[i] - center[i-1]
#    A[i,i-1] = -u/2/c2c_backward - D/c2f_backward**2      
#    A[i,i] = 2*D/c2c_backward**2 
#    A[i,i+1] = u/2/c2c_forward - D/c2f_forward**2 
    A1[i,i-1] = (-u/2/c2c_backward - D/c2c_backward**2)      
    A1[i,i] = (2*D/c2c_backward**2) 
    A1[i,i+1] = (u/2/c2c_forward - D/c2c_forward**2)
    
#b[n] = 2/itr**2
A1[0,0] = 1; b1[0] = 0 #A[0,1] = - 2*D/dx**2; A[0,2] = - u/(2*dx) - D/dx**2 ;
#A[n,n] = 3; A[n,n-1] = -4; A[n,n-2] = 1;
A1[n-1,n-1] = - D/c2c_forward**2; A1[n-1,n-2] = 2*D/c2c_backward**2; A1[n-1,n-3] = -D/c2c_backward**2;#A[n,n-2] = 1;
#A[n-1,n-1] = A[n-1,n-2]
b1[n-1] = V[n-1]/width #2/n**2 # need to fix this

phi1 = np.linalg.solve(A1,b1)

plt.figure(1,figsize=(6, 4),dpi=50)
plt.plot(center,phi1, label='phi')
#plt.ylim([-2, 2])
plt.title("Phi")
plt.xlabel('midpoint of segment')
plt.ylabel('Phi(x_n)')
#plt.legend()
#plt.show()
