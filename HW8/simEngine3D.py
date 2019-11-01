import numpy as np

from GCons.DP1 import *
from GCons.DP2 import *
from GCons.CD import *
from GCons.D import *

def f(t):
    return 0 #np.sin(np.pi/4. * np.cos(2.* t))

def df(t):
    return 0#-np.pi/2. * np.sin(2.* t) * np.cos(np.pi/4.* np.cos(2.* t))

def ddf(t):
    return 0#-np.pi/4. * (4.* np.cos(2.* t) * np.cos(np.pi/4.* np.cos(2* t)) + \
            #np.pi * np.sin(2.* t) * np.sin(2.* t) * np.sin(np.pi/4.* np.cos(2.* t)))

# calculate p given A
def get_p(A):
    e0 = np.sqrt(0.25 * (np.trace(A) + 1.))
    e1 = np.sqrt(0.25 * (-np.trace(A) + 1. + 2.* A[0][0]))
    e2 = np.sqrt(0.25 * (-np.trace(A) + 1. + 2.* A[1][1]))
    e3 = np.sqrt(0.25 * (-np.trace(A) + 1. + 2.* A[2][2]))
    return np.array([e0, e1, e2, e3])

def getE(p):
    e0, e1, e2, e3 = p[0], p[1], p[2], p[3]
    return np.array([[-e1, e0, -e3, e2],
                     [-e2, e3, e0, -e1],
                     [-e3, -e2, e1, e0]])

def getG(p):
    e0, e1, e2, e3 = p[0], p[1], p[2], p[3]
    return np.array([[-e1, e0, e3, -e2],
                     [-e2, -e3, e0, e1],
                     [-e3, e2, -e1, e0]])

def getJp(p, m, w, L):
    # here L is the length of pendulum
    # in L-RF
    J_bar = 1./12. * m * np.array([2.* w**2, L**2 + w**2, L**2 + w**2]) * np.eye(3)
    G = getG(p)

    return 4. * np.matmul(np.transpose(G), np.matmul(J_bar, G))


def getPhi(x, y, z, x_p, z_p,
           r_i, p_i, s_i_P_bar,
           r_j, p_j, s_j_Q_bar, t, j_ground):
    # CD constraints at point Q
    phi_CD_x = CD_phi(x, r_i, p_i, s_i_P_bar, r_j, p_j, s_j_Q_bar, 0.)
    phi_CD_y = CD_phi(y, r_i, p_i, s_i_P_bar, r_j, p_j, s_j_Q_bar, 0.)
    phi_CD_z = CD_phi(z, r_i, p_i, s_i_P_bar, r_j, p_j, s_j_Q_bar, 0.)
#    print(phi_CD_z)
    # DP1 constraints
    phi_DP1_y = DP1_phi(p_i, z_p, p_j, y, 0.)
    phi_DP1_z = DP1_phi(p_i, z_p, p_j, z, 0.)

    # driving constraint
    phi_DP1_drive = DP1_phi(p_i, x_p, p_j, y, f(t))

    # euler parameter normalization constraint
    phi_p = 0.5 * (np.dot(p_i, p_i) - 1.)

    return np.array([phi_CD_x, phi_CD_y, phi_CD_z,
                     phi_DP1_y, phi_DP1_z,
#                     phi_DP1_drive,
                     phi_p])

def getNu(t):
    nu_CD_x  = CD_nu(0.)
    nu_CD_y  = CD_nu(0.)
    nu_CD_z  = CD_nu(0.)
    nu_DP1_y = DP1_nu(0.)
    nu_DP1_z = DP1_nu(0.)
    nu_DP1_drive = DP1_nu(df(t))
    nu_p = 0.
    return np.array([nu_CD_x, nu_CD_y, nu_CD_z,
                     nu_DP1_y, nu_DP1_z,
#                     nu_DP1_drive,
                     nu_p])

def getPhi_q(x, y, z, x_p, z_p,
             r_i, p_i, s_i_P_bar,
             r_j, p_j, s_j_Q_bar, t, j_ground):
    # get phi_r for 7 constraints
    phi_r_CD_x  = CD_phi_r(x, j_ground)
    phi_r_CD_y  = CD_phi_r(y, j_ground)
    phi_r_CD_z  = CD_phi_r(z, j_ground)
    phi_r_DP1_y = DP1_phi_r(j_ground)
    phi_r_DP1_z = DP1_phi_r(j_ground)
    phi_r_DP1_drive = DP1_phi_r(j_ground)
    phi_r_p = np.zeros(3)

    # get phi_p for 7 constraints
    phi_p_CD_x = CD_phi_p(x, p_i, s_i_P_bar, p_j, s_j_Q_bar, j_ground)
    phi_p_CD_y = CD_phi_p(y, p_i, s_i_P_bar, p_j, s_j_Q_bar, j_ground)
    phi_p_CD_z = CD_phi_p(z, p_i, s_i_P_bar, p_j, s_j_Q_bar, j_ground)
    phi_p_DP1_y = DP1_phi_p(p_i, z_p, p_j, y, j_ground)
    phi_p_DP1_z = DP1_phi_p(p_i, z_p, p_j, z, j_ground)
    phi_p_DP1_drive = DP1_phi_p(p_i, x_p, p_j, y, j_ground)
    phi_p_p = p_i

    return np.array([np.concatenate((phi_r_CD_x, phi_p_CD_x)),
                     np.concatenate((phi_r_CD_y, phi_p_CD_y)),
                     np.concatenate((phi_r_CD_z, phi_p_CD_z)),
                     np.concatenate((phi_r_DP1_y, phi_p_DP1_y)),
                     np.concatenate((phi_r_DP1_z, phi_p_DP1_z)),
#                     np.concatenate((phi_r_DP1_drive, phi_p_DP1_drive)),
                     np.concatenate((phi_r_p, phi_p_p))])

def getGamma(x, y, z, x_p, z_p,
             r_i, r_i_dot, p_i, p_i_dot, s_i_P_bar,
             r_j, r_j_dot, p_j, p_j_dot, s_j_Q_bar, t, j_ground):
    gamma_CD_x = CD_gamma(x, p_i_dot, s_i_P_bar, p_j_dot, s_j_Q_bar, 0.)
    gamma_CD_y = CD_gamma(y, p_i_dot, s_i_P_bar, p_j_dot, s_j_Q_bar, 0.)
    gamma_CD_z = CD_gamma(z, p_i_dot, s_i_P_bar, p_j_dot, s_j_Q_bar, 0.)

    gamma_DP1_y = DP1_gamma(p_i, p_i_dot, z_p, p_j, p_j_dot, y, 0.)
    gamma_DP1_z = DP1_gamma(p_i, p_i_dot, z_p, p_j, p_j_dot, z, 0.)
    gamma_DP1_drive = DP1_gamma(p_i, p_i_dot, x_p, p_j, p_j_dot, y, ddf(t))
    gamma_DP1_p = -np.dot(p_i_dot, p_i_dot)

    return np.array([gamma_CD_x, gamma_CD_y, gamma_CD_z,
                     gamma_DP1_y, gamma_DP1_z,
#                     gamma_DP1_drive,
                     gamma_DP1_p])
