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
W = np.zeros([3])
DW = np.zeros([itr])
A = np.zeros([itr,itr])
b = np.zeros([itr])
h = x[1]

V = np.ones(itr)/n

# not sure what to put for rij
rij = h
eij = rij
q = rij/h
W[0] = 1/(6*h)*(2-q)**2 - 4*(1-q)**2
W[1] =  1/(6*h)*(2-q)**2
W[2] =  1/(6*h)*(2-q)**2 

DW[0] = -3*eij/(6*h**2)*(2-q)**2 - 4*(1-q)**2
DW[1] = -3*eij/(6*h**2)*(2-q)**2
DW[2] = -3*eij/(6*h**2)*(2-q)**2 

# Using a 1D Analysis in the x direction
# want to solve for vectore phi

for i in np.arange(2,itr-2,1):
#    print(i)
    b[i] = V[i] # ?
    A[i,i-1] = V[i]*W[0] + DW[0]  
    A[i,i] = V[i]*W[1] + DW[1]  
    A[i,i+1] = V[i]*W[2] + DW[2]  

# The phi[0] and phi[1] should be equal and opposite to phi[3] and phi[4] respectively
A[2,0] = 1; b[0] = 0
A[1,0] = 1; b[0] = 0
A[0,0] = 1; b[0] = 0
# The phi[itr-1] and phi[itr-2] may be equal to 0
A[n-1,n-1] = 1; 
A[n-1,n-2] = 1; 
A[n-1,n-3] = 1; # these ones are not correct; not sure what number to put
b[n-1] = V[n-1] 
b[n-2] = V[n-3] 
b[n-3] = V[n-3] 

# couldn't solve the problem yet
#phi1 = np.linalg.solve(A,b)

plt.figure(1,figsize=(6, 4),dpi=50)
plt.plot(x,phi, label='phi')
#plt.ylim([-2, 2])
plt.title("Phi")
plt.xlabel('x location')
plt.ylabel('Phi(x_n)')
#plt.legend()
#plt.show()