import numpy as np
import matplotlib.pyplot as plt

u = 1 # u(x,t)
D = 1 # D(x,t)
s = 1 # s(x,t)

# Laplace: u=0 and s!=0
# Diffusion u=0 and s=0

n = 20
itr = n + 5

x = np.linspace(0., 1., num=itr)
phi = np.zeros([itr])
W = np.zeros([4])
psi2 = np.zeros([itr])
A = np.zeros([itr,itr])
b = np.zeros([itr])
dx = x[1]

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
#plt.legend()
#plt.show()