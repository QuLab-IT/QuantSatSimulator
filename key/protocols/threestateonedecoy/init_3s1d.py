import numpy as np

__all__ = ['get_x_bounds','arrange_output','x0_rand','check_constraints',
           'bool_constraints']

###############################################################################

def get_x_bounds(opt_dict,num_min):
    """
    Returns the intial values and upper & lower bounds for the optimised 
    parameters.

    Parameters
    ----------
    opt_dict : dict
        Dictionary of parameters related to optimisation.
    mu3 : float
        Intensity of pulse 3 (vacuum).
    num_min : float
        An arbitrarily small number.

    Returns
    -------
    x : float, array
        Optimised parameters initial values.
    xb : float, array
        Optimised parameters upper & lower bounds.

    """
    xb = np.zeros((4,2))
    xb[0,:] = opt_dict['Px'][3:None]
    xb[1,:] = opt_dict['P1'][3:None]
    xb[2,:] = opt_dict['mu1'][3:None]
    xb[3,:] = opt_dict['mu2'][3:None]

    x = np.zeros((4,))
    if opt_dict['Px'][1]:
        # Initial value for Px is specified
        x[0] = opt_dict['Px'][2]
    else:
        # Initial value for Px is randomly allocated
        x[0] = np.random.uniform(xb[0,0]+num_min , xb[0,1]-num_min)
    if opt_dict['P1'][1]:
        # Initial values for both P1 are specified
        x[1] = opt_dict['P1'][2]
    else:
        # P1 is randomly allocated
        x[1] = np.random.uniform(xb[1,0]+num_min , xb[1,1]-num_min)
    if opt_dict['mu1'][1]:
        # Initial value for mu1 is specified
        x[2] = opt_dict['mu1'][2]
    else:
        # Initial value for mu1 is randomly allocated
        x[2] = np.random.uniform(xb[2,0]+num_min , xb[2,1]-num_min)
    if opt_dict['mu2'][1]:
        # Initial value for mu2 is specified
        x[3] = opt_dict['mu2'][2]
    else:
        # Initial value for mu2 is randomly allocated
        x[3] = np.random.uniform(xb[3,0]+num_min , min(xb[3,1],x[2])-num_min)
    return x, xb

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def arrange_output(SysLoss,ls,dt,mu3,QBERI,Pec,Pbx,x,SKLdata,sys_params):
    """
    Put calculation output parameters into an ordered list.

    Parameters
    ----------
    sysLoss : float
        The nominal loss or system loss metric (dB). For symmetric transmission
        windows this is the system loss at zenith.
    ls : float
        Excess system loss (dB).
    dt : int
        Transmission window half-width (s).
    mu3 : float
        Fixed protocol parameter. The intensity of the third pulse.
    QBERI : float
        Intrinsic quantum bit error rate.
    Pec : float
        Extraneous count rate/probability.
    theta_max : float
        Maximum elevation of satellite overpass (rad).
    x : float, array-like
        The initial set of parameters to use.
    SKLdata : float, array-like
        Calculation output parameters.
    sys_params : dict
        Additional system parameters to be included in fulldata.

    Returns
    -------
    list
        Ordered list of data to write out.

    """
    # dt (s),ls (dB),QBERI,Pec,maxElev (deg),SKL (b),QBER,phiX,nX,nZ,mX,lambdaEC,
    # sX0,sX1,vZ1,sZ1,mean photon no.,PxA,PxB,P1,P2,mu1,mu2,mu3,eps_c,eps_s,Pap,NoPass,
    # fs (Hz),minElev (deg),shiftElev (deg),SysLoss (dB)
    return [dt,ls,QBERI,Pec,sys_params[0],*SKLdata,x[0],Pbx,x[1],1-x[1],
            x[2],x[3],mu3,*sys_params[1:],SysLoss]

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def x0_rand(mu3,xb,num_min):
    """
    Randomly initialise the 5 protocol parameters using the specified bounds.
    Parameters and bounds should be specified in the order {Px,pk1,pk2,mu1,mu2}.

    Parameters
    ----------
    mu3 : float
        Intensity of pulse 3 (vacuum).
    xb : float, array-like
        Upper and lower bounds for the protocol parameters. (5,2)
    num_min : float
        An arbitrarily small number.

    Returns
    -------
    x0 : float, array
        Randomly initialised protocol parameters.

    """
    Px_i  = np.random.uniform(xb[0,0]+num_min , xb[0,1]-num_min)
    pk1_i = np.random.uniform(xb[1,0]+num_min , xb[1,1]-num_min)
    mu1_i = np.random.uniform(xb[2,0]+num_min , xb[2,1]-num_min)
    mu2_i = np.random.uniform(xb[3,0]+num_min , min(xb[3,1],mu1_i)-num_min)
    return np.array([Px_i,pk1_i,mu1_i,mu2_i])

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def check_constraints(Px,pk1,pk2,mu1,mu2,mu3):
    """
    Check that the parameters are within the bounds and constraints of the
    asymmetric BB84 protocol with weak coherent pulses with 2 'decoy' states.
    Stops the script if any bounds or constraints are violated.

    Parameters
    ----------
    Px : float
        Asymmetric polarisation probability.
    pk1 : float
        Probability Alice sends pulse intensity 1.
    pk2 : float
        Probability Alice sends pulse intensity 2.
    mu1 : float
        Intensity of pulse 1.
    mu2 : float
        Intensity of pulse 2.
    mu3 : float
        Intensity of pulse 3.

    Returns
    -------
    None.

    """
    # Constraint 1: Check polarisation basis probabilities are valid.
    if (Px >= 1.0 or Px <= 0.0):
        print("Error! Constraint 1 < Px < 0: ", Px)
        exit(1)
    # Constraint 2: Check probability of pulse with intensity 1 is in bounds.
    if (pk1 >= 1.0 or pk1 <= 0.0):
        print("Error! Constraint 1 < pk1 < 0: ", pk1)
        exit(1)
    # Constraint 3: Check probability of pulse with intensity 2 is in bounds.
    if (pk2 >= 1.0 or pk2 <= 0.0):
        print("Error! Constraint 1 < pk2 < 0: ", pk2)
        exit(1)
    # Constraint 4: Check sum of probabilities for intensity 1 & 2 are less
    # than unity.
    if ((pk1 + pk2) >= 1.0):
        print("Error! Constraint (pk1 + pk2) < 1: ", pk1 + pk2)
        exit(1)
    # Constraint 5: Check value of intensity 1 is in bounds.
    if (mu1 >= 1.0 or mu1 <= 0.0):
        print("Error! Constraint 1 < mu1 < 0: ", mu1)
        exit(1)
    # Constraint 6: Check value of intensity 2 is in bounds.
    if (mu2 >= 1.0 or mu2 <= 0.0):
        print("Error! Constraint 1 < mu2 < 0: ", mu2)
        exit(1)
    # Constraint 7: Check values of all intensities are in bounds.
    if ((mu1 - mu3) <= mu2):
        print("Error! Constraint (mu1-mu3) > mu2: ", (mu1-mu3), mu2)
        exit(1)
    # Constraint 8: Check values of intensities 2 & 3 are in bounds.
    if (mu2 <= mu3):
        print("Error! Constraint mu2 > mu3: ", mu2, mu3)
        exit(1)
    return None

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def bool_constraints(Px,pk1,mu1,mu2,mu3):
    """
    Check that the parameters are within the bounds and constraints of the
    asymmetric BB84 protocol with weak coherent pulses with 2 'decoy' states.
    Returns a boolean array corresponding to each of the constraints.

    Parameters
    ----------
    Px : float
        Asymmetric polarisation probability.
    pk1 : float
        Probability Alice sends pulse intensity 1.
    pk2 : float
        Probability Alice sends pulse intensity 2.
    mu1 : float
        Intensity of pulse 1.
    mu2 : float
        Intensity of pulse 2.
    mu3 : float
        Intensity of pulse 3.

    Returns
    -------
    C : boolean, array-like.
        Do the parameters satisfy the constraints? True or False

    """
    C = np.array([1,1,1,1,1,1,1,1], dtype=bool) # Initialise array as True
    
    # Constraint 1: Check polarisation basis probabilities are valid.
    if (Px >= 1.0 or Px <= 0.0):
        C[0] = False
    # Constraint 2: Check probability of pulse with intensity 1 is in bounds.
    if (pk1 >= 1.0 or pk1 <= 0.0):
        C[1] = False
    # Constraint 3: Check value of intensity 1 is in bounds.
    if (mu1 >= 1.0 or mu1 <= 0.0):
        C[4] = False
    # Constraint 6: Check value of intensity 2 is in bounds.
    if (mu2 >= 1.0 or mu2 <= 0.0):
        C[5] = False
    return C