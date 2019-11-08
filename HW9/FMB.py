import numpy as np
import matplotlib.pyplot as plt

def finite_diff(A,b,x,phi,u,D,s):
    for i in np.arange(1,n,1):
        dx = x[1]
        b[i] = s
        A[i,i-1] = - u/(2*dx) - D/dx**2         
        A[i,i] = 2*D/dx**2
        A[i,i+1] =  u/(2*dx) - D/dx**2 
        
    A[0,0] = 1; b[0] = 0 
    A[n,n] = - D/dx**2; A[n,n-1] = 2*D/dx**2; A[n,n-2] = - D/dx**2
    b[n] = s 
    
    phi = np.linalg.solve(A,b)
    return phi

def finite_volume(A1,b1,phi1,V,center,edge,u,D,s):
    for i in np.arange(1,n-1,1):
        b1[i] = V[i] 
        c2f_forward = edge[i+1] - center[i]
        c2f_backward = center[i] - edge[i]
        c2c_forward =  center[i+1] - center[i]
        c2c_backward = center[i] - center[i-1]
        A1[i,i-1] = (-u*(c2f_backward/c2c_backward) - D/c2c_backward)      
        A1[i,i] = (2*D/c2c_backward) 
        A1[i,i+1] = (u*(c2f_forward/c2c_forward) - D/c2c_forward)
        
    A1[0,0] = 1; b1[0] = 0 
    A1[n-1,n-1] = - D/c2c_forward; A1[n-1,n-2] = 2*D/c2c_backward; A1[n-1,n-3] = -D/c2c_backward;
    b1[n-1] = V[n-1] 

    phi = np.linalg.solve(A1,b1)
    return phi

u = 1 # u(x,t)
D = 1 # D(x,t)
s = 1 # s(x,t)

# Laplace: u=0 and s!=0
# Diffusion u=0 and s=0

n = 1000
itr = n + 1

x = np.linspace(0., 1., num=itr)
phi = np.zeros([itr])
A = np.zeros([itr,itr])
b = np.zeros([itr])

phi = finite_diff(A,b,x,phi,u,D,s)

V = np.ones(itr)/n
edge = np.linspace(0., 1., num=itr) # edges of segments
center = np.linspace(0, 1., num=itr) + 1/2/n # edges of segments
center = center[:-1]
phi1 = np.zeros([n])
A1 = np.zeros([n,n])
b1 = np.zeros([n])

phi1 = finite_volume(A1,b1,phi1,V,center,edge,u,D,s)

plt.figure(1,figsize=(6, 4),dpi=50)
plt.plot(x,phi, label='phi')
plt.plot(center,phi1, label='phi')
plt.title("Phi")
plt.xlabel('x location (midpoint for volume method)')
plt.ylabel('Phi(x_n)')
plt.legend(['Finite Diff','Finite Volume'])
#plt.show()

## 
# integral of d^2phi/dx^2 = s  == phi = s*x^2/2
# integral of dphi/dx = s  == phi = s*x

#phi_analytical1 = np.exp(-s*np.sqrt(x)/D)
#phi_analytical1 = np.exp(s*(x)**2/D)
phi_analytical1 = -s*(x)**2/D
phi_analytical2 = u*x
#phi_a = phi_analytical2 - phi_analytical1

plt.figure(2,figsize=(6, 4),dpi=50)
plt.plot(x,phi_analytical1, label='phi')
plt.plot(x,phi_analytical2, label='phi')
#plt.plot(x,phi, label='phi')
#plt.plot(x,phi_a, label='phi')
#plt.ylim([0,2])
plt.title("Phi")
plt.xlabel('x location (midpoint for volume method)')
plt.ylabel('Phi(x_n)')
plt.legend(['Finite Diff','Finite Volume'])

