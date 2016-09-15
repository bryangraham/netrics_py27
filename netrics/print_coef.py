# Ensure "normal" division
from __future__ import division

# Load library dependencies
import numpy as np

# Define print_coef() function
#-----------------------------------------------------------------------------#

def print_coef(beta, vcov, var_names=None):
    
    """
    AUTHOR: Bryan S. Graham, bgraham@econ.berkeley.edu, June 2016
    
    This function prints out a list of variable names, coefficient estimates and standard errors
    in a unified format. Assume beta and var_names are conformably ordered
    
    INPUTS
    -------
    beta        : K vector of estimated coefficients
    vcov        : K x K estimated variance-covariance matrix (2d numpy array)
    var_names   : list of length K with variable names
    
    
    OUTPUTS
    -------
    This functionr returns "None". All "output" is to screen.
        
    """

    # if var_names is None then label variables X0, X1...etc.
    
    if not var_names:
        var_names = []
        for k in range(0,len(beta)):
            var_names.append("X_" + str(k))
            
    print ""
    print "Independent variable       Coef.    ( Std. Err.) "
    print "-------------------------------------------------------------------------------------------"
        
    c = 0
    for names in var_names:
        print names.ljust(25) + "%10.6f" % beta[c] + \
                         " (" + "%10.6f" % np.sqrt(vcov[c,c]) + ")"
        c += 1
                
    print ""
    print "-------------------------------------------------------------------------------------------"
       
    
    return None