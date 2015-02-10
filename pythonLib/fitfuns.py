# Anthony Ho, ahho@stanford.edu, 2/3/2015
# Last update 2/4/2015
""" Library containing some commonly used mathematical functions and their derivatives """


import numpy as np



## ------ ------ ##

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



## ------ ------ ##

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



## ------ ------ ##

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



## ------ ------ ##

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



## ------ ------ ##

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
