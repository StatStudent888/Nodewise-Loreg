import sys
import os
import datetime
import random

from simulate_data import *

sys.path.append('../algorithms')
from neighborhood_continuation import *
from MIO import *


# Runs an example to illustrate using our code.

## Generate synthetic data. Simulation parameters

# Type of covariance matrix of covariates [see simulate_data.py file in the same folder for details] (Another option for type_Sigma = 2)
type_Sigma = 2

N, P = 50, 50
k0   = 7
rho  = 0.2
SNR  = 2

X, l2_X, real_beta, eps_train, eps_val, y_train, y_val = simulate_data(type_Sigma, N, P, k0, rho, SNR)


## MIO: Run the MIO solver to obtain solutions to the L0-LQ Problem
MIO_L0_LQ('l1',   X, y_train, 7, 0.1, time_limit=30)
##MIO_L0_LQ('l2^2', X, y_train, 7, 0.1, time_limit=30)

