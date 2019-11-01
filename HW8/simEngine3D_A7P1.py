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

def main():
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
    s_i_P_bar = np.array([-L,0,0])
    r_i = np.array([0, np.sin(np.pi/4.) * L, -np.cos(np.pi/4.) * L])
    r_i_dot = np.array([0,0,0])
    A = np.array([[0, 0, 1],
                  [ np.sin(np.pi/4.), np.cos(np.pi/4.), 0],
                  [-np.cos(np.pi/4.), np.sin(np.pi/4.), 0]])
    p_i = get_p(A)
    p_i_dot = np.array([0,0,0,0])

    # call functions
    phi   = getPhi(x, y, z, x_p, z_p,
                   r_i, p_i, s_i_P_bar,
                   r_j, p_j, s_j_Q_bar, t, j_ground)
    phi_q = getPhi_q(x, y, z, x_p, z_p,
                     r_i, p_i, s_i_P_bar,
                     r_j, p_j, s_j_Q_bar, t, j_ground)
#    print(phi)
    

    # inverse dynamics settings
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
    q_prev_1 = []; q_prev_2 = []; q_new = np.zeros([7])
    Zn = np.zeros([15]) # 3rddot 4pddot 1lambdap 7lambda
    # create arrays to save time evolution results
    # torque
    torque = np.zeros((num_steps, 3))

    q_prev_1 = np.concatenate((r_i, p_i))

    while step < 1:
#    while step < num_steps:
        # get real time associated
        t = step * dt


######################    
        # Seeding for initial conditions
        if step < 1:
            while diff > tol and iter < 1:#max_iter:
                q_new = q_prev_1 - dt*Zn[:7]
                q_i_dot = - dt*beta*Zn[:7] 
                
                phi   = getPhi(x, y, z, x_p, z_p,
                               q_new[:3], q_new[3:], s_i_P_bar,
                               r_j, p_j, s_j_Q_bar, t, j_ground)
                phi_q = getPhi_q(x, y, z, x_p, z_p,
                                 q_new[:3], q_new[3:], s_i_P_bar,
                                 r_j, p_j, s_j_Q_bar, t, j_ground)
    
                # velocity analysis
                nu      = getNu(t)
    #            print(phi_q,'\n',nu)
                q_i_dot = np.linalg.solve(phi_q, nu)
    #    
    #            # step 1: getting accelerations
    #            # acceleration analysis
                gamma    = getGamma(x, y, z, x_p, z_p,
                                    q_new[:3], q_i_dot[:3], q_new[3:], q_i_dot[3:], s_i_P_bar,
                                    r_j, r_j_dot, p_j, p_j_dot, s_j_Q_bar, t, j_ground)
                q_i_ddot = np.linalg.solve(phi_q, gamma)
    ##            
    ##           # Get r_ddot, p_ddot, lambda_p, lambda        
    ##
                rhs = np.concatenate((g - np.matmul(m * np.eye(3), q_i_ddot[:3]), # r
                                      np.matmul(-getJp(q_new[3:], m, width, 2.* L), q_i_ddot[3:]))) # p
                lamb = np.linalg.solve(np.transpose(phi_q), rhs) # lambda
                lamb_p = 
                Pi = 0.5 * np.matmul(phi_q[:,3:], np.transpose(getE(q_new[3:]))) # norm
                # z matrix = [rddot; pddot; lamb_p; lamb]
                Zn = np.concatenate((q_i_ddot, lamb))
                print(lamb)
                q_prev_1 = q_new
                iter += 1 
            
######################     
        diff = float("inf")
        iter = 0

        for iter in range(1):
#        while diff > tol and iter < max_iter:
            
            # Computation of pos and vel using 4th order BDF and most recent accel
            q_new[:3] = 4/3*q_prev_1[:3] - 1/3*q_prev_2[:3] - beta**2*dt**2*Zn[:3]
            q_new[3:7] = - beta**2*dt**2*Zn[3:7] 
#            print(iter,q_new)
                
#            q_new = q_old - np.matmul(np.linalg.inv(phi_q), phi)
#            diff  = np.linalg.norm(q_new - q_old)
#            q_old = q_new
            
            # Find constraint and jacobian
            phi   = getPhi(x, y, z, x_p, z_p,
                           q_new[:3], q_new[3:], s_i_P_bar,
                           r_j, p_j, s_j_Q_bar, t, j_ground)
            phi_q = getPhi_q(x, y, z, x_p, z_p,
                             q_new[:3], q_new[3:], s_i_P_bar,
                             r_j, p_j, s_j_Q_bar, t, j_ground)

            # velocity analysis
#            nu      = getNu(t)
#            print(phi_q,'\n',nu)
#            q_i_dot = np.linalg.solve(phi_q, nu)
#    
#            # step 1: getting accelerations
#            # acceleration analysis
            gamma    = getGamma(x, y, z, x_p, z_p,
                                q_new[:3], q_i_dot[:3], q_new[3:], q_i_dot[3:], s_i_P_bar,
                                r_j, r_j_dot, p_j, p_j_dot, s_j_Q_bar, t, j_ground)
            q_i_ddot = np.linalg.solve(phi_q, gamma)
##            
##           # Get r_ddot, p_ddot, lambda_p, lambda        
##
            rhs = np.concatenate((g - np.matmul(m * np.eye(3), q_i_ddot[:3]), # r
                                  np.matmul(-getJp(q_new[3:], m, width, 2.* L), q_i_ddot[3:]))) # p
            lamb = np.linalg.solve(np.transpose(phi_q), rhs) # lambda
            Pi = 0.5 * np.matmul(phi_q[:,3:], np.transpose(getE(q_new[3:]))) # norm
            P = np.array([q_new[3:]])
            
            # z matrix = [rddot; pddot; lamb_p; lamb]
            Zn = np.concatenate((q_i_ddot, lamb))
            print(Zn)

            G = np.array(([-g + np.matmul(m * np.eye(3), q_i_ddot[:3]) + np.matmul(np.transpose(phi_q[:,:3]),lamb)],
                         [np.matmul(-getJp(q_new[3:], m, width, 2.* L), q_i_ddot[3:]) + np.matmul(np.transpose(phi_q[:,3:]),lamb) + np.transpose(P)*lamb[3:]],
                         [1/(beta**2*dt**2)*phi_q[:,3:]],
                         [1/(beta**2*dt**2)*phi]))
#            print(G)
            h = dt
            if iter == 0:
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
#            delta = np.linalg.solve(Psi,-G)
            delta = 0.0001
            
            Zn = Zn + delta
            # Improve quality of solution
            q_prev_2 = q_prev_1
            q_prev_1   = q_new
            
            iter += 1          
            # end of while loop
            print(delta)
        # step 3: populate torque for x, y, z in G-RF
#        torque[step,0] = -Pi[5,0] * lamb[5]
#        torque[step,1] = -Pi[5,1] * lamb[3]
#        torque[step,2] = -Pi[5,0] * lamb[4]

        step += 1

    # save to files
    save2file(torque, "data/Q_torque.txt")

if __name__ == "__main__":
    if os.path.exists("data"):
        shutil.rmtree("data")
    os.makedirs("data")
    t0 = time.time()
    main()
    t1 = time.time()
    print("Simulation finished in ",t1-t0, "seconds.")

    # plotting
    plot_from_file("data/Q_torque.txt", "Reaction Torque (Nm)")
