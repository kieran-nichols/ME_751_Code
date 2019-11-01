import numpy as np

# function tilde, getA, and getB are repeated in each GCon

# calculate a_tilde
def tilde(a):
    ax, ay, az = a[0], a[1], a[2]
    t = np.array([[0, -az, ay],
                  [az, 0, -ax],
                  [-ay, ax, 0]])
    return t

# get A matrix from euler parameters
def getA(p):
    e0, e1, e2, e3 = p[0], p[1], p[2], p[3]
    A = 2 * np.array([[e0**2 + e1**2 - 0.5, e1*e2 - e0*e3, e1*e3 + e0*e2],
                      [e1*e2 + e0*e3, e0**2 + e2**2 - 0.5, e2*e3 - e0*e1],
                      [e1*e3 - e0*e2, e2*e3 + e0*e1, e0**2 + e3**2 - 0.5]])
    return A

# get B matrix from p and a_bar
def getB(p, a_bar):
    B     = np.zeros((3, 4))
    e0, e = p[0], p[1:]
    # first column of B matrix
    c = e0 * np.identity(3) + tilde(e)
    # 3*3 matrix inside B
    m = np.matmul(c, tilde(a_bar))

    B[:, 0]  = np.dot(c, a_bar)
    B[:, 1:] = np.outer(e, a_bar) - m
    B *= 2.
    return B

# DP2 constraint
def DP2_phi(a_i_bar, r_i, p_i, s_i_P_bar,
                     r_j, p_j, s_j_Q_bar, f):
    A_i   = getA(p_i)
    A_j   = getA(p_j)
    s_i_P = np.matmul(A_i, s_i_P_bar)
    s_j_Q = np.matmul(A_j, s_j_Q_bar)
    d_ij  = r_j + s_j_Q - r_i - s_i_P

    return np.dot(np.matmul(A_i, a_i_bar), d_ij) - f

def DP2_phi_r(a_i_bar, p_i, j_ground=False):
    res = np.zeros(6)
    A_i = getA(p_i)
    a_i = np.matmul(A_i, a_i_bar)

    res[:3] = -np.transpose(a_i)
    res[3:] = np.transpose(a_i)
    if j_ground:
        return res[:3]
    else:
        return res

def DP2_phi_p(a_i_bar, r_i, p_i, s_i_P_bar,
              a_j_bar, r_j, p_j, s_j_Q_bar, j_ground=False):
    res = np.zeros(8)
    A_i = getA(p_i)
    A_j = getA(p_j)

    a_i   = np.matmul(A_i, a_i_bar)
    s_i_P = np.matmul(A_i, s_i_P_bar)
    s_j_Q = np.matmul(A_j, s_j_Q_bar)
    d_ij  = r_j + s_j_Q - r_i - s_i_P

    # note that in slide 33 of lecture10 there is a typo
    # for the calculation of phi_p_i
    res[:4] = np.dot(d_ij, getB(p_i, a_i_bar)) - \
              np.dot(a_i,  getB(p_i, s_i_P_bar))
    res[4:] = np.dot(a_i,  getB(p_j, s_j_Q_bar))
    if j_ground:
        return res[:4]
    else:
        return res

def DP2_nu(df):
    return df

def DP2_gamma(a_i_bar, r_i, p_i, s_i_P_bar, r_i_dot, p_i_dot,
                       r_j, p_j, s_j_Q_bar, r_j_dot, p_j_dot, ddf):
    A_i      = getA(p_i)
    A_j      = getA(p_j)
    a_i      = np.matmul(A_i, a_i_bar)
    a_i_dot  = np.dot(getB(p_i, a_i_bar), p_i_dot)
    s_i_P    = np.matmul(A_i, s_i_P_bar)
    s_j_Q    = np.matmul(A_j, s_j_Q_bar)
    d_ij     = r_j + s_j_Q - r_i - s_i_P
    d_ij_dot = r_j_dot + np.matmul(getB(p_j,s_j_Q_bar),p_j_dot) - \
               r_i_dot - np.matmul(getB(p_i,s_i_P_bar),p_i_dot)
    gamma = -np.dot(a_i,  np.matmul(getB(p_j_dot, s_j_Q_bar), p_j_dot)) + \
             np.dot(a_i,  np.matmul(getB(p_i_dot, s_i_P_bar), p_i_dot)) - \
             np.dot(d_ij, np.matmul(getB(p_i_dot, a_i_bar), p_i_dot))  - \
             2. * np.dot(a_i_dot, d_ij_dot) + ddf
    return gamma
