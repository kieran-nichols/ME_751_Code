import numpy as np
import matplotlib.pyplot as plt

u = 1 # u(x,t)
D = 1 # D(x,t)
s = 1 # s(x,t)

n = 100
m = 500
itr = n + 1
itr_m = m + 1

x = np.linspace(0., 1., num=itr)
time = np.linspace(0., 5., num=itr_m)
dt = time[1]
phi = np.zeros([itr])
A = np.zeros([itr,itr])
b = np.zeros([itr])
dx = x[1]

A[0,0] = 1; b[0] = 0 #A[0,1] = - 2*D/dx**2; A[0,2] = - u/(2*dx) - D/dx**2 ;
A[n,n] = - D/dx**2; A[n,n-1] = 2*D/dx**2; A[n,n-2] = - D/dx**2#A[n,n-2] = 1;
b[n] = s 
Phi = np.zeros([itr,itr_m])

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

plt.figure(1,figsize=(6, 4),dpi=50)
plt.plot(x,Phi[:,itr_m-2], label='phi')
#plt.plot(x,Phi, label='phi')
plt.title("Phi")
plt.xlabel('x_n')
plt.ylabel('Phi(x_n)')

V = np.ones(itr)/n
width = V[1]
edge = np.linspace(0., 1., num=itr) # edges of segments
center = np.linspace(0, 1., num=itr) + 1/2/n # edges of segments
center = center[:-1]
phi1 = np.zeros([n])
A1 = np.zeros([n,n])
b1 = np.zeros([n])
Phi1 = np.zeros([itr-1,itr_m-1])

c2f_forward = edge[2] - center[1]
c2f_backward = center[1] - edge[1]
c2c_forward =  center[2] - center[1]
c2c_backward = center[1] - center[0]

A1[0,0] = 1; b1[0] = 0 
A1[n-1,n-1] = - D/c2c_forward; A1[n-1,n-2] = 2*D/c2c_backward; A1[n-1,n-3] = -D/c2c_backward;
b1[n-1] = V[n-1] 

# Using a 1D Analysis in the x direction
# want to solve for vectore phi
for j in np.arange(1,m,1):
    for i in np.arange(1,n-1,1):
        b1[i] = V[i] + Phi[i,j-1]/dt # how do I change these equations
        c2f_forward = edge[i+1] - center[i]
        c2f_backward = center[i] - edge[i]
        c2c_forward =  center[i+1] - center[i]
        c2c_backward = center[i] - center[i-1]
        A1[i,i-1] = (-u*(c2f_backward/c2c_backward) - D/c2c_backward)      
        A1[i,i] = (2*D/c2c_backward) + 1/dt 
        A1[i,i+1] = (u*(c2f_forward/c2c_forward) - D/c2c_forward)
    phi1 = np.linalg.solve(A1,b1)
    Phi1[:,j] = phi1 
    
plt.figure(2,figsize=(6, 4),dpi=50)
plt.plot(center,Phi1[:,itr_m-2], label='phi')
#plt.ylim([-2, 2])
plt.title("Phi")
plt.xlabel('midpoint of segment')
plt.ylabel('Phi(x_n)')
#plt.legend(['Finite Diff','Finite Volume'])
#plt.show()

# what do I do which j=1?
for j in np.arange(2,m,1):
    for i in np.arange(1,n-1,1):
        b1[i] = V[i] + 4*Phi[i,j-1]/dt/2 - Phi[i,j-2]/dt/2 # how do I change these equations
        c2f_forward = edge[i+1] - center[i]
        c2f_backward = center[i] - edge[i]
        c2c_forward =  center[i+1] - center[i]
        c2c_backward = center[i] - center[i-1]
        A1[i,i-1] = (-u*(c2f_backward/c2c_backward) - D/c2c_backward)      
        A1[i,i] = (2*D/c2c_backward) + 3/dt/2 
        A1[i,i+1] = (u*(c2f_forward/c2c_forward) - D/c2c_forward)
    phi1 = np.linalg.solve(A1,b1)
    Phi1[:,j] = phi1 
    
plt.figure(3,figsize=(6, 4),dpi=50)
plt.plot(center,Phi1[:,itr_m-2], label='phi')
#plt.ylim([-2, 2])
plt.title("Phi")
plt.xlabel('midpoint of segment')
plt.ylabel('Phi(x_n)')
#plt.legend(['Finite Diff','Finite Volume'])
#plt.show()