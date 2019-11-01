import numpy as np
import os, shutil, time
import matplotlib.pyplot as plt

from GCons.DP1 import *
from GCons.DP2 import *
from GCons.CD import *
from GCons.D import *
from simEngine3D import *

def save2file(arr, fout):
    with open(fout, 'w') as f:
        for i in range(len(arr)):
            f.write(str(arr[i,0]) + "\t" + str(arr[i,1]) + \
            "\t" + str(arr[i,2]) + "\n")

def plot_from_file(fin, rva):
    x, y, z = [], [], []
    with open(fin, 'r') as f:
        for line in f:
            x.append(float(line.split()[0]))
            y.append(float(line.split()[1]))
            z.append(float(line.split()[2]))
    time = np.arange(len(x))*1e-3
    plt.plot(time, x, label="x", linewidth=1.5)
    plt.plot(time, y, label="y", linewidth=1.5)
    plt.plot(time, z, label="z", linewidth=1.5)
    plt.xlabel('Time (Sec)',fontsize=14)
    plt.ylabel(rva,fontsize=14)
    plt.legend(prop={'size':12})
    plt.grid(True)
    plt.show()

t0 = time.time()
# define parameters
t = 0.
L = 2.
width = 0.05
rho   = 7800.
m     = 2.* L * width**2 * rho
g     = np.array([0, 0, -m * 9.81])

# body j is ground
j_ground = True

# set up G-RF
x = np.array([1,0,0])
y = np.array([0,1,0])
z = np.array([0,0,1])

# set up L-RF
x_p = np.array([1,0,0])
y_p = np.array([0,1,0])
z_p = np.array([0,0,1])

# settings for CD constraints
# position of Q in G-RF, and it's fixed
s_j_Q_bar = np.array([0,0,0])
r_j = np.array([0,0,0])
p_j = np.array([1,0,0,0])
r_j_dot = np.array([0,0,0])
p_j_dot = np.array([0,0,0,0])

# position of Q in L-RF, consider theta = pi/4 as initial position
# distance between O' and Q are fixed
s_i_P_bar_1 = np.array([-L,0,0])
s_i_P_bar_2 = np.array([-L/2,0,0])
r_i_1 = np.array([0, np.sin(np.pi/4.) * L, -np.cos(np.pi/4.) * L])
r_i_1_2 = np.array([0, np.sin(np.pi/4.) * L/2, -np.cos(np.pi/4.) * L/2])
r_i_dot_1 = np.array([0,0,0])
r_i_dot_2 = np.array([0,0,0])
A = np.array([[0, 0, 1],
              [ np.sin(np.pi/4.), np.cos(np.pi/4.), 0],
              [-np.cos(np.pi/4.), np.sin(np.pi/4.), 0]])
p_i = get_p(A)
p_i_dot = 0.001*np.ones([4])

# call functions
phi_1   = getPhi(x, y, z, x_p, z_p,
               r_i, p_i, s_i_P_bar,
               r_j, p_j, s_j_Q_bar, t, j_ground)
phi_2   = getPhi(x, y, z, x_p, z_p,
               r_i_1, p_i_1, s_i_P_bar_1,
               r_j, p_j, s_j_Q_bar, t, j_ground)
phi_q = getPhi_q(x, y, z, x_p, z_p,
                 r_i, p_i, s_i_P_bar,
                 r_j, p_j, s_j_Q_bar, t, j_ground)
# guess of gamma
#gamma    = getGamma(x, y, z, x_p, z_p,
#                    r_i, r_i_dot, p_i, p_i_dot, s_i_P_bar,
#                    r_j, r_j_dot, p_j, p_j_dot, s_j_Q_bar, t, j_ground)
#print(gamma)


# inverse dynamics settings
n_constraints = 6
t_start, t_end = 0, 10
dt = 1e-3
num_steps = int((t_end-t_start) / dt)
step = 0
max_iter = 10
tol = 1e-15
beta = 2/3
#    z = []
G = []
delta = 0.
Psi = []
q_prev_1 = np.zeros([7]); q_prev_2 = np.zeros([7]); q_new = np.zeros([7])
qd_prev_1 = np.zeros([7]); qd_prev_2 = np.zeros([7]); qd_new = np.zeros([7])
q_dd = 0.001*np.ones([7])
Zn = 0.001*np.ones([15]) # 3rddot 4pddot 1lambdap 7lambda

# create arrays to save time evolution results
pos_1 = np.zeros([num_steps,3])
ang_vel_1 = np.zeros([num_steps,3])
nu_norm_1 = np.zeros([num_steps,3])
pos_2 = np.zeros([num_steps,3])
ang_vel_2 = np.zeros([num_steps,3])
nu_norm_2 = np.zeros([num_steps,3])
torque_2 = np.zeros((num_steps, 3))

q_prev_1 = np.concatenate((r_i, p_i))

#while step < 1:
while step < num_steps:
        # get real time associated
    t = step * dt
#tn = 0
    iter = 0
    
    while delta > tol and iter < max_iter:
        t = tn
        q_dd = Zn[:7]
        lamb = Zn[8:]
        lamb_p = Zn[7]
        
        if step == 0:
            # Newton order 1
            q_new_1 = q_prev_1_1 + beta**2*dt**2*q_dd_1
            qd_new_1 = qd_prev_1_1 + dt*beta*q_dd_1 
            q_new_2 = q_prev_1_2 + beta**2*dt**2*q_dd_2
            qd_new_2 = qd_prev_1_2 + dt*beta*q_dd_2 
        else:
            # Newton order 2
            q_new = 4/3*q_prev_1[:3] - 1/3*q_prev_2[:3] - beta**2*dt**2*Zn
            qd_new = 4/3*qd_prev_1[:3] - 1/3*qd_prev_2[:3]- beta*dt*Zn
        
        phi   = getPhi(x, y, z, x_p, z_p,
                       q_new[:3], q_new[3:], s_i_P_bar,
                       r_j, p_j, s_j_Q_bar, t, j_ground)
        phi_q = getPhi_q(x, y, z, x_p, z_p,
                         q_new[:3], q_new[3:], s_i_P_bar,
                         r_j, p_j, s_j_Q_bar, t, j_ground)
        #    
        #            # step 1: getting accelerations
        #            # acceleration analysis
        F = g
        T_hat = np.zeros([4])
        P = np.transpose(q_new[3:])
        gamma_p = np.matmul(P,q_dd[3:])
        gamma    = getGamma(x, y, z, x_p, z_p,
                            q_new[:3], qd_new[:3], q_new[3:], qd_new[3:], s_i_P_bar,
                            r_j, r_j_dot, p_j, p_j_dot, s_j_Q_bar, t, j_ground)
        #    print(gamma_p.shape,gamma.shape)
        RHS = np.hstack((F, T_hat, np.array([gamma_p]), gamma)).reshape(([14,1]))
        
        
        # LHS set up
        x2 = np.zeros([3,4])
        x3 = np.zeros([3])
        x6 = np.zeros([3,1])
        x4 = np.zeros([n_constraints])
        x5 = np.zeros([n_constraints,n_constraints])
        x6 = np.zeros([1])
        
        # Setting up the 4 LHS terms
        M = m*np.eye(3)
        J_p = getJp(q_new[3:], m, width, 2.* L)
        # for some reason when I use the x1,..,x6 variables the arrays don't concat
        LHS1 = np.hstack((M, x2, np.zeros([3,1]), np.transpose(phi_q[:,:3])))
        LHS2 = np.hstack((np.transpose(x2), J_p, np.transpose(P.reshape((1,4))),  np.transpose(phi_q[:,3:])))
        LHS3 = np.hstack((np.transpose(x3), P, x6, np.transpose(x4)))
        LHS4 = np.hstack((phi_q[:,:3], phi_q[:,3:], np.zeros([6,1]), x5))
        #    print(LHS1.shape,LHS2.shape,LHS3.reshape([1,14]).shape,LHS4.shape)
        LHS = np.concatenate([LHS1,LHS2,LHS3.reshape([1,14]),LHS4],axis=0)  
        #print(RHS)
        
        # Solve for accel and lamda
        # Tested with an offset due to original singular matrix result; gave tiny perturbation
        # to originally 0 number guesses; still giving singular matrix for offset of 0
        offset = -1
        Zn = np.linalg.solve(LHS[:offset,:offset], RHS[:offset])
        print(Zn)
        #Zn = np.linalg.solve(LHS, RHS)
        ##            
        G = np.array(([-g + np.matmul(M, Zn[:3]) + np.matmul(np.transpose(phi_q[:,:3]),Zn[8:])],
                     [np.matmul(-getJp(q_new[3:], m, width, 2.* L), Zn[3:7]) + np.matmul(np.transpose(phi_q[:,3:]),Zn[8:]) + np.transpose(P)*Zn[8:]],
                     [1/(beta**2*dt**2)*phi_q[:,3:]],
                     [1/(beta**2*dt**2)*phi]))
        
        if iter == 0:
            #
            x2 = np.zeros([1,3])
            x3 = np.zeros([1])
            x4 = np.zeros([7,1])
            x5 = np.zeros([7,7])
            
            # Assume Psi terms approx by the 0 order terms?
            phirlamb_rr = 0.
            Fr= 0.
            F_rdot = 0.
            phiLamb_rp = 0.
            F_p = 0.
            F_pdot = 0.
            phiLamb_pr = 0.
            TauHat_r = 0.
            TauHat_rdot = 0.
            JpPddot_p = 0.
            PLambP_p = 0.
            phiLamb_pp = 0.
            TauHat_p = 0.
            TauHat_pdot = 0.
            h = dt
            # Setting up the 4 Psi terms
            Psi_11 = m*np.eye(3) + h**2*beta**2*phirlamb_rr - h**2*beta**2*Fr -h*beta*F_rdot
            Psi_12 = h**2*beta**2*phiLamb_rp - h**2*beta**2*F_p - h*beta*F_pdot
            Psi_21 = np.zeros([4,1]) #h**2*beta**2*phiLamb_pr - h**2*beta**2*TauHat_r - h*beta*TauHat_rdot
            Psi_22 = getJp(q_new[3:], m, width, 2.* L) + h**2*beta**2*JpPddot_p + h**2*beta**2*PLambP_p + h**2*beta**2*phiLamb_pp - h**2*beta**2*TauHat_p - h*beta*TauHat_pdot
            #                print(phi_q[:,:3].shape,phi_q[:,3:].shape,x4.shape,x5.shape)
            Psi = np.array(([Psi_11, Psi_12, np.transpose(x2), np.transpose(phi_q[:,:3])],
                           [Psi_21, Psi_22, np.transpose(P),  np.transpose(phi_q[:,3:])],
                           [x2, P, x3, np.transpose(x4)],
                           [phi_q[:,:3], phi_q[:,3:], x4, x5]))
            #            print(Psi)
        
        # Calculate the Newton correction
        delta = np.linalg.solve(Psi,-G)
        #delta = 0.0001
        
        Zn = Zn + delta
        
        # Improve quality of solution
        q_prev_2 = q_prev_1
        q_prev_1   = q_new
        iter += 1
    pos_1[step,:] = q_new_1[:3]
    ang_vel_1[step,:] = qd_new_1[:3]
    nu_norm_1[step,1] = np.linalg.norm(getNu_1(t))
    pos_2[step,:] = q_new_2[:3]
    ang_vel_2[step,:] = qd_new_2[:3]
    step += 1     
t1 = time.time()    
print('finished')  
######################     

# save to files
if os.path.exists("data"):
    shutil.rmtree("data")
os.makedirs("data")

save2file(pos_1, "data/Pos_1.txt")
save2file(ang_vel_1, "data/Vel_1.txt")
save2file(nu_norm_1, "data/Norm_1.txt")
save2file(pos_2, "data/Pos_2.txt")
save2file(ang_vel_2, "data/Vel_2.txt")

#main()
print("Simulation finished in ",t1-t0, "seconds.")

# plotting
plot_from_file("data/Pos_1.txt", "Translation_1 (m)")
plot_from_file("data/Vel_1.txt", "Ang_Vel_1 (rad/s)")
plot_from_file("data/Norm_1.txt", "Norm_1 ()")
plot_from_file("data/Pos_2.txt", "Translation_2 (m)")
plot_from_file("data/Vel_2.txt", "Ang_Vel_2 (rad/s)")

