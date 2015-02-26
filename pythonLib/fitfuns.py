# Anthony Ho, ahho@stanford.edu, 2/3/2015
# Last update 2/4/2015
"""Library containing some commonly used mathematical functions and their derivatives"""


import numpy as np
import scipy.special


# ------ Linear equation ------ #


# Compute A*x + B
# where params[0] = A
#       params[1] = B
def linear(params, x):
    return params[0]*x + params[1]


# Compute the residual of A*x + B
# where params[0] = A
#       params[1] = B
def linearResidual(params, x, y):
    return params[0]*x + params[1] - y


# ------ Single exponential with tau ------ #


# Compute A*exp(-t/tau) + C
# where params[0] = A
#       params[1] = tau
#       params[2] = C
def singleExp(params, x):
    return params[0]*np.exp(-x/params[1]) + params[2]


# Compute the Jacobian of singleExp
# Returns a 2D-ndarray of shape (len(x), len(params))
def singleExpPrime(params, x):
    partial_p0s = np.exp(-x/params[1])
    partial_p1s = params[0]*x*np.exp(-x/params[1])/(params[1]**2)
    partial_p2s = np.ones(len(x))
    return np.column_stack((partial_p0s, partial_p1s, partial_p2s))


# Compute the residual of A*exp(-t/tau) + C - y
# where params[0] = A
#       params[1] = tau
#       params[2] = C
def singleExpResidual(params, x, y):
    return params[0]*np.exp(-x/params[1]) + params[2] - y


# ------ Single exponential with k ------ #


# Compute A*exp(-kt) + C
# where params[0] = A
#       params[1] = k
#       params[2] = C
def singleExpK(params, x):
    return params[0]*np.exp(-params[1]*x) + params[2]


# Compute the Jacobian of singleExpK
# Returns a 2D-ndarray of shape (len(x), len(params))
def singleExpKPrime(params, x):
    partial_p0s = np.exp(-params[1]*x)
    partial_p1s = -params[0]*x*np.exp(-params[1]*x)
    partial_p2s = np.ones(len(x))
    return np.column_stack((partial_p0s, partial_p1s, partial_p2s))


# Compute the residual of A*exp(-kt) + C - y
# where params[0] = A
#       params[1] = k
#       params[2] = C
def singleExpKResidual(params, x, y):
    return params[0]*np.exp(-params[1]*x) + params[2] - y


# ------ Double exponential with taus ------ #


# Compute (A1*exp(-t/tau1) + A2)*exp(-t/tau2) + C
# where params[0] = A1
#       params[1] = tau1
#       params[2] = A2
#       params[3] = tau2
#       params[4] = C
def doubleExp(params, x):
    return params[0]*np.exp(-x*(1/params[1]+1/params[3])) + params[2]*np.exp(-x/params[3]) + params[4]


# Compute the Jacobian of doubleExp
# Returns a 2D-ndarray of shape (len(x), len(params))
def doubleExpPrime(params, x):
    partial_p0s = np.exp(-x*(1/params[1]+1/params[3]))
    partial_p1s = params[0]*x*np.exp(-x*(1/params[1]+1/params[3]))/(params[1]**2)
    partial_p2s = np.exp(-x/params[3])
    partial_p3s = params[2]*x*np.exp(-x/params[3])/(params[3]**2)
    partial_p4s = np.ones(len(x))
    return np.column_stack((partial_p0s, partial_p1s, partial_p2s, partial_p3s, partial_p4s))


# Compute the residual of (A1*exp(-t/tau1) + A2)*exp(-t/tau2) + C - y
# where params[0] = A1
#       params[1] = tau1
#       params[2] = A2
#       params[3] = tau2
#       params[4] = C
def doubleExpResidual(params, x, y):
    return params[0]*np.exp(-x*(1/params[1]+1/params[3])) + params[2]*np.exp(-x/params[3]) + params[4] - y


# ------ Double exponential with k's ------ #


# Compute (A1*exp(-k1t) + A2)*exp(-k2t) + C
# where params[0] = A1
#       params[1] = k1
#       params[2] = A2
#       params[3] = k2
#       params[4] = C
def doubleExpK(params, x):
    return params[0]*np.exp(-(params[1]+params[3])*x) + params[2]*np.exp(-params[3]*x) + params[4]


# Compute the Jacobian of doubleExp
# Returns a 2D-ndarray of shape (len(x), len(params))
def doubleExpKPrime(params, x):
    partial_p0s = np.exp(-(params[1]+params[3])*x)
    partial_p1s = -params[0]*x*np.exp(-(params[1]+params[3])*x)
    partial_p2s = np.exp(-params[3]*x)
    partial_p3s = -params[2]*x*np.exp(-params[3]*x)
    partial_p4s = np.ones(len(x))
    return np.column_stack((partial_p0s, partial_p1s, partial_p2s, partial_p3s, partial_p4s))


# Compute the residual of (A1*exp(-k1t) + A2)*exp(-k2t) + C - y
# where params[0] = A1
#       params[1] = k1
#       params[2] = A2
#       params[3] = k2
#       params[4] = C
def doubleExpKResidual(params, x, y):
    return params[0]*np.exp(-(params[1]+params[3])*x) + params[2]*np.exp(-params[3]*x) + params[4] - y


# ------ Double exponential with taus, constant tau2 ------ #


# Compute (A1*exp(-t/tau1) + A2)*exp(-t/tau2) + C
# where params[0] = A1
#       params[1] = tau1
#       params[2] = A2
#       params[3] = C
#       constants[0] = tau2
def doubleExpConstT2(params, x, constants):
    return params[0]*np.exp(-x*(1/params[1]+1/constants[0])) + params[2]*np.exp(-x/constants[0]) + params[3]


# Compute the Jacobian of doubleExpConstT2
# Returns a 2D-ndarray of shape (len(x), len(params))
def doubleExpConstT2Prime(params, x, constants):
    partial_p0s = np.exp(-x*(1/params[1]+1/constants[0]))
    partial_p1s = params[0]*x*np.exp(-x*(1/params[1]+1/constants[0]))/(params[1]**2)
    partial_p2s = np.exp(-x/constants[0])
    partial_p3s = np.ones(len(x))
    return np.column_stack((partial_p0s, partial_p1s, partial_p2s, partial_p3s))


# Compute the residual of (A1*exp(-t/tau1) + A2)*exp(-t/tau2) + C - y
# where params[0] = A1
#       params[1] = tau1
#       params[2] = A2
#       params[3] = C
#       constants[0] = tau2
def doubleExpConstT2Residual(params, x, y, constants):
    return params[0]*np.exp(-x*(1/params[1]+1/constants[0])) + params[2]*np.exp(-x/constants[0]) + params[3] - y


# ------ Double exponential with taus, constant tau2 and C ------ #


# Compute (A1*exp(-t/tau1) + A2)*exp(-t/tau2) + C
# where params[0] = A1
#       params[1] = tau1
#       params[2] = A2
#       constants[0] = tau2
#       constants[1] = C
def doubleExpConstT2C(params, x, constants):
    return params[0]*np.exp(-x*(1/params[1]+1/constants[0])) + params[2]*np.exp(-x/constants[0]) + constants[1]


# Compute the Jacobian of doubleExpConstT2C
# Returns a 2D-ndarray of shape (len(x), len(params))
def doubleExpConstT2CPrime(params, x, constants):
    partial_p0s = np.exp(-x*(1/params[1]+1/constants[0]))
    partial_p1s = params[0]*x*np.exp(-x*(1/params[1]+1/constants[0]))/(params[1]**2)
    partial_p2s = np.exp(-x/constants[0])
    return np.column_stack((partial_p0s, partial_p1s, partial_p2s))


# Compute the residual of (A1*exp(-t/tau1) + A2)*exp(-t/tau2) + C - y
# where params[0] = A1
#       params[1] = tau1
#       params[2] = A2
#       constants[0] = tau2
#       constants[1] = C
def doubleExpConstT2CResidual(params, x, y, constants):
    return params[0]*np.exp(-x*(1/params[1]+1/constants[0])) + params[2]*np.exp(-x/constants[0]) + constants[1] - y


# ------ Distributions ------ #


# Compute the pmf of Poisson distribution, pre-normalized
def poisson(params, x):
    return params[0]*(params[1]**x)*np.exp(-params[1])/scipy.special.gamma(x+1)
