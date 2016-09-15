# Ensure "normal" division
from __future__ import division

# Load library dependencies
import numpy as np
import scipy as sp
import scipy.optimize
import itertools as it

from numba import jit
import numexpr as ne

from logit import logit

from print_coef import print_coef
from helpers import generate_dyad_to_tetrads_dict, generate_tetrad_indices, \
                    organize_data_tetrad_logit, tetrad_logit_score_proj

# Define tetrad_logit() function
#-----------------------------------------------------------------------------#

def tetrad_logit(D, W, dtcon=None, silent=False, W_names=None):
    
    """
    AUTHOR: Bryan S. Graham, bgraham@econ.berkeley.edu, June 2016
    
    This function computes the Tetrad Logit estimator introduced in Graham (2014, NBER No. 20341) -- "An Econometric
    Model of Link Formation with Degree Heterogeneity". The implementation is as described in the paper. Notation
    attempts to follow that used in the paper.
    
    INPUTS
    ------
    D                 : N x N undirected adjacency matrix
    W                 : List with elements consisting of N x N 2d numpy arrays of dyad-specific 
                        covariates such that W[k][i,j] gives the k-th covariate for dyad ij
    dtcon             : Dyad and tetrad concordance (dtcon) List with elements [tetrad_to_dyads_indices, 
                        dyad_to_tetrads_dict]. If dtcon=None, then construct it using generate_tetrad_indices() 
                        function. See header to generate_tetrad_indices() for more information.
    silent            : If True then suppress all optimization and estimation output, show output otherwise.  
    W_names           : List of K strings giving names of the columns of W_tilde. If silent=False then use
                        these in presentation of estimation output.
                   
    OUTPUTS
    -------
    beta_TL          :  K x 1 vector of tetrad logit point estimates of beta
    vcov_beta_TL     :  K x K asymptotic-variance matrix for beta_TL 
                        NOTE: vcov_beta_TL is already "divided by n" (just take square root of diagonal for std. errs.)
    tetrad_frac_TL   :  Fraction of tetrads that contribute to Tetrad Logit criterion function                    
    success          :  corresponds to success component of OptimizeResult associated with Scipy's minimize function;
                        success = True if the tetrad logit optimizer exited successfully
            
    
    CALLS:           : ...logit()...
                       ...organize_data_tetrad_logit()...
    ------
    
    """
    
    # ------------------------------------------------------- #
    # - STEP 1: Prepare data for estimation                 - #
    # ------------------------------------------------------- #
    
    # compute size of the network and dimension of regressor matrix
    K        = len(W)                    # Number of dyad-specific covariates
    N        = np.shape(D)[0]            # Number of agents in the network
    n        = N*(N-1) // 2              # Number of dyads in network  
    Nchoose4 = N*(N-1)*(N-2)*(N-3) // 24 # Number of tetrads in network
    
    # organize data for input into logit maximizer
    [S, W_tilde, tetrad_frac_TL, proj_tetrads_dict] = organize_data_tetrad_logit(D, W, dtcon)
    
    # NOTES:  S  is a 6 (N choose 4) x 1 vector of -1, 0, 1 according to the configuration of the tetrad. The
    #         6 (N choose 4) elements correspond to the six dyads in each of the N choose 4 tetrads.
    #         W_tilde is a 6 (N choose 4) x K matrix of regressors corresponding to S. This is as described in the 
    #         paper.
    
    # Drop all rows where S = 0 (i.e., no identifying content & not part of the criterion function)
    g       = S.nonzero()                                      # row & column indices of for all elements of S 
                                                               # that are non-zero 
                                                               # (i.e., equal to -1 or 1)
    Y_trim  = (0*(S[g[0],:]==-1) + 1*(S[g[0],:]==1))           # Make outcome binary
    W_trim  = W_tilde[g[0],:]                   
            
    # ------------------------------------------------------- #
    # - STEP 2: Compute Tetrad Logit Point Estimate         - #
    # ------------------------------------------------------- #
    
    [beta_TL, hess_TL, success_TL] = logit(Y_trim, W_trim, nocons=True, \
                                               silent=silent)                            # TL estimates of beta
    beta_TL                        = np.reshape(beta_TL,(-1,1))                          # Put beta_TL into 
                                                                                         # 2-dimensional form
        
    # ------------------------------------------------------- #
    # - STEP 3: Compute Variance-Covariance Matrix          - #
    # ------------------------------------------------------- #
   
    # ------------------------------------- #
    # - Compute covariance of the "score" - #
    # ------------------------------------- #
   
    # place full "score" vector, including non-contributing tetrads, into
    # a 6 (N choose 4) x 1 two dimensional scipy sparse matrix
    score = sp.sparse.coo_matrix((np.ravel(Y_trim - (1 + np.exp(-np.dot(W_trim,beta_TL)))**-1),\
                                 (g[0],g[1])), shape = (6*Nchoose4, 1), dtype='float64').tocsr()
    
    # Compute n x K matrix containing all six components of each dyad's contribution to the
    # projection of the "score". These correspond to the six non-redudant permutation of ij, kl
    # which enter the criterion function as described in Graham (2014, NBER). Note the use of
    # list comprehensions in what follows
    
    # -------------------------------------------------------------- #
    # - Compute n x K proj_score matrix using a list comprehension - #
    # -------------------------------------------------------------- #
    
    #---------------------------#
    #- Serial implementation   -#
    #---------------------------#
    # pass in relevant score components to compute projection for each dyad
    proj_score = np.array([tetrad_logit_score_proj([score[tetrads,:].toarray(), W_tilde[tetrads,:]]) \
                           for dyad, tetrads in proj_tetrads_dict.iteritems()])/(n - 2*(N-1) + 1)
    
    # compute the covariance matrix of the score projection
    OMEGA_hat = np.cov(proj_score, rowvar=False)*((n-1)/(n-K))
           
    # ------------------------------------- #
    # - Compute inverse hessian matrix    - #
    # ------------------------------------- #
    
    iGAMMA_hat = np.linalg.inv(-hess_TL/Nchoose4)
    
    # Sandwich variance-covariance matrix estimate for beta_TL
    vcov_beta_TL = 36*np.dot(np.dot(iGAMMA_hat,OMEGA_hat),iGAMMA_hat)/n
    
    # ------------------------------------------------------- #
    # - STEP 4: Report estimation results                   - #
    # ------------------------------------------------------- #
    
    if not silent:
        print ""
        print "-------------------------------------------------------------------------------------------"
        print "- TETRAD LOGIT ESTIMATION RESULTS                                                         -"
        print "-------------------------------------------------------------------------------------------"
        print ""
        print "Number of agents,           N : " + "{:>15,.0f}".format(N)
        print "Number of dyads,            n : " + "{:>15,.0f}".format(n)
        print "Number of tetrads             : " + "{:>15,.0f}".format(Nchoose4) 
        print "Number identifying tetrads    : " + "{:>15,.0f}".format(tetrad_frac_TL*Nchoose4)
        print "Fraction identifying tetrads  : " + "{:>15.6f}".format(tetrad_frac_TL)
        print ""
        print "-------------------------------------------------------------------------------------------"
        print_coef(beta_TL, vcov_beta_TL, W_names)
    
    
    return [beta_TL, vcov_beta_TL, tetrad_frac_TL, success_TL]