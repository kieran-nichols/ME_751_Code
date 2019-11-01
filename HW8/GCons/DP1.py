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

# DP1 constraint
def DP1_phi(p_i, a_i_bar, p_j, a_j_bar, f):
    A_i = getA(p_i)
    A_j = getA(p_j)
    a_i = np.matmul(A_i, a_i_bar)
    a_j = np.matmul(A_j, a_j_bar)
    return np.dot(a_i, a_j) - f

def DP1_phi_r(j_ground=False):
    if j_ground:
        return np.zeros(3)
    else:
        return np.zeros(6)

def DP1_phi_p(p_i, a_i_bar, p_j, a_j_bar, j_ground=False):
    res = np.zeros(8)
    A_i = getA(p_i)
    A_j = getA(p_j)
    a_i = np.dot(A_i, a_i_bar)
    a_j = np.dot(A_j, a_j_bar)

    res[:4] = np.dot(a_j, getB(p_i, a_i_bar))
    res[4:] = np.dot(a_i, getB(p_j, a_j_bar))
    if j_ground:
        return res[:4]
    else:
        return res

def DP1_nu(df):
    return df

def DP1_gamma(p_i, p_i_dot, a_i_bar, p_j, p_j_dot, a_j_bar, ddf):
    a_i_dot = np.dot(getB(p_i, a_i_bar), p_i_dot)
    a_j_dot = np.dot(getB(p_j, a_j_bar), p_j_dot)
    A_i     = getA(p_i)
    A_j     = getA(p_j)
    a_i     = np.matmul(A_i, a_i_bar)
    a_j     = np.matmul(A_j, a_j_bar)
    gamma   = -np.dot(a_i, np.matmul(getB(p_j_dot, a_j_bar), p_j_dot)) - \
               np.dot(a_j, np.matmul(getB(p_i_dot, a_i_bar), p_i_dot)) - \
               2. * np.dot(a_i_dot, a_j_dot) + ddf
    return gamma
