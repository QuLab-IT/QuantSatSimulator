import numpy as np

__all__ = ['DRate_j','error_j','nxz','mxz','mk_j','nXZpm','nXZpm_HB','nXZpm_inf',
           'tau','s0','s1','vxz1','mean_photon_a','Qber']

###############################################################################

def DRate_j(Pap,Pec,exp_loss_jt):
    """
    Calculates the expected detection rate including afterpulse contributions
    for each intensity and time slot.
    Defined as R_k in Sec. IV of [1].

    Parameters
    ----------
    eta : float
        Excess loss parameter.
    Pap : float
        Probability of an afterpulse event.
    Pec : float
        Extraneous count probability.
    exp_loss_jt : float, array
        Loss, per intensity per time slot, decay function.

    Returns
    -------
    float, array
        Expected detection rate.

    """
    return (1 + Pap)*(1 - (1 - 2*Pec)*exp_loss_jt)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def error_j(Dj,Pap,Pec,QBERI,exp_loss_jt):
    """
    Calculates the conditional probability for a pulse of intensity mu_j
    to cause an error, after sifting, in the time slot t.
    Defined as e_k in Sec. IV of [1].

    Parameters
    ----------
    Dj : float, array
        Expected detection rate.
    Pap : float
        Probability of an afterpulse event.
    Pec : float
        Extraneous count probability.
    QBERI : float
        Intrinsic Quantum Bit Error Rate.
    exp_loss_jt : float, array
        Loss, per intensity per time slot, decay function.

    Returns
    -------
    float, array
        Error rate per intensity per time slot.

    """
    return Pec + (0.5*Pap*Dj) + QBERI*(1 - exp_loss_jt)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def nxz(PAxz,PBxz,Rate,dt,P_times_Dj):
    """
    Calculates the number of events in the X or Z sifted basis per pulse
    intensity per time slot

    Parameters
    ----------
    PAxz : float
        Probability of Alice preparing a state in the X/Z basis.
    PBxz : float
        Probability of Bob measuring a state in the X/Z basis.
    Npulse : integer/float
        Number of pulses sent by Alice.
    P_times_Dj : float, array
        The element-wise multiplication of the intensity probability array P
        with the expected detection rate per time slot array Dj.

    Returns
    -------
    float, array
        The number of events in the sifted X/Z basis.

    """
    return PAxz*PBxz*P_times_Dj*Rate*dt

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def mxz(n,P_dot_Dj,P_dot_ej):
    """
    Calculates the number of errors in the sifted X or Z basis for each time 
    slot.
    ----------
    n : float
        The number of events in the sifted X/Z basis per time
        slot.
    P_dot_Dj : float
        The dot product of the intensity probability array P with the expected 
        detection rate per time slot array Dj.
    P_dot_ej : float
        The dot product of the intensity probability array P with the 
        conditional probability for a pulse with a given intensity and time 
        slot to create an error array ej.

    Returns
    -------
    float
        The number of errors in the sifted X/Z basis per time slot.

    """
    return P_dot_ej*n / P_dot_Dj

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def mk_j(mj,P_dot_Dj,P_times_Dj):
    """
    Calculates the number of errors in the sifted X or Z basis for a pulse
    with a given intensity in a particular time slot.
    ----------
    mj : float, array
        The number of errors in the sifted X/Z basis per time slot.
    P_dot_Dj : float
        The dot product of the intensity probability array P with the expected 
        detection rate per time slot array Dj.
    P_times_Dj : float, array
        The element-wise multiplication of the intensity probability array P
        with the expected detection rate per time slot array Dj.

    Returns
    -------
    float, array
        Number of errors in sifted X/Z basis per intensity per time slot.

    """


    return np.divide( np.multiply(mj,P_times_Dj) , P_dot_Dj )

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def nXZpm(mu,P,nxz_mu,eps):
    """
    Calculates the upper and lower bounds on the number of events in the 
    sifted X or Z basis, for each pulse intensity using, the Chernoff bounds.
    Defined after Eq. (2) in [1].
        nXplus[j] and nXmin[j], or nZplus[j] and nZmin[j];  j = {1:3}

    Parameters
    ----------
    mu : float, array
        Pulse intensities.
    P : float, array
        Probability of Alice preparing a pulse intensity.
    nxz_mu : float, array
        Number of events per intensity for the sifted X/Z basis.
    eps_s : float
        The secrecy error; the key is eps_s secret.

    Returns
    -------
    nXZmin : float, array
        Lower bound on the expected number of events per intensity in the
        sifted X/Z basis.
    nXZplus : float, array
        Upper bound on the expected number of events per intensity in the
        sifted X/Z basis.

    """
    log_21es = np.math.log(1 / eps)
    term_m   = 0.5*log_21es + np.sqrt(2*nxz_mu*log_21es + 0.25*log_21es**2)
    term_p   = log_21es + np.sqrt(2*nxz_mu*log_21es + log_21es**2)

    nXZmin   = np.divide(np.multiply(np.exp(mu), nxz_mu - term_m), P)
    nXZplus  = np.divide(np.multiply(np.exp(mu), nxz_mu + term_p), P)
    return nXZmin, nXZplus

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def nXZpm_HB(mu,P,nxz_mu,nXZ,eps):
    """
    Calculates the upper and lower bounds on the number of events in the 
    sifted X or Z basis, for each pulse intensity using, the Hoeffding bound.
    Defined after Eq. (2) in [1].
        nXplus[j] and nXmin[j], or nZplus[j] and nZmin[j];  j = {1:3}

    Parameters
    ----------
    mu : float, array
        Pulse intensities.
    P : float, array
        Probability of Alice preparing a pulse intensity.
    nxz_mu : float, array
        Number of events per intensity for the sifted X/Z basis.
    nXZ : float, array
        Number of events per intensity per time slot for the sifted X/Z basis.
    eps_s : float
        The secrecy error; the key is eps_s secret.

    Returns
    -------
    nXZmin : float, array
        Lower bound on the expected number of events per intensity in the
        sifted X/Z basis.
    nXZplus : float, array
        Upper bound on the expected number of events per intensity in the
        sifted X/Z basis.

    """
    term2   = np.sqrt(0.5*nXZ * np.math.log(1 / eps))
    nXZmin  = np.divide(np.multiply(np.exp(mu), nxz_mu - term2), P)
    nXZplus = np.divide(np.multiply(np.exp(mu), nxz_mu + term2), P)
    return nXZmin, nXZplus

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def nXZpm_inf(mu,P,nxz_mu):
    """
    Calculates the number of events in the  sifted X or Z basis, for each pulse 
    intensity in the asymptotic limit.
    Defined after Eq. (2) in [1].
        nXi[j], or nZi[j];  j = {1:3}

    Parameters
    ----------
    mu : float, array
        Pulse intensities.
    P : float, array
        Probability of Alice preparing a pulse intensity.
    nxz_mu : float, array
        Number of events per intensity for the sifted X/Z basis.

    Returns
    -------
    nXZi : float, array
        The expected number of events per intensity in the
        sifted X/Z basis.

    """
    nXZi = np.divide(np.multiply(np.exp(mu), nxz_mu), P)
    return nXZi, nXZi

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def tau(n,mu,P):
    """
    Calculates the total probability that Alice prepares an n-photon state.
    Defined after Eq. (2) in [1].

    Parameters
    ----------
    n : integer
        Number of photons in the state.
    mu : float, array
        Intensities of the weak coherent pulses.
    P : float, array
        Probabilities that Alice prepares a particular intensity.

    Returns
    -------
    tau : float
        Total probability of an n-photon state.

    """
    if (not isinstance(n, int)): # Python 3.
        #if (not isinstance(n, (int, long))): # Python 2.
        print("Error! n must be an integer: ", n)
        exit(1)
    tau = 0
    for jj in range(len(mu)):
        tau += np.math.exp(-mu[jj]) * mu[jj]**n * P[jj]
    tau = tau / np.math.factorial(n)
    return tau

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def s0(mu,P,nMin):
    """
    Calculates the approximate number of vacuum events in the sifted X or Z 
    basis.
    See Eq. (2) in [1].

    Parameters
    ----------
    mu : float, array
        Pulse intensities.
    P : float, array
        Probabilities that Alice prepares a particular intensity.
    nMin : float, array
        Lower bound on the expected number of events per intensity in the
        sifted X/Z basis.

    Returns
    -------
    float
        The number of vacuum events in the sifted X or Z basis.

    """
    return tau(0,mu,P) * (mu[0]*nMin[1] - mu[1]*nMin[0]) / (mu[0] - mu[1])

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def s1(mu,P,nMin,nPlus,s0=None):
    """
    Calculates the number of single photon events in the sifted X or Z basis.
    See Eq. (3) in [1].

    Parameters
    ----------
    mu : float, array
        Pulse intensities.
    P : float, array
        Probabilities that Alice prepares a particular intensity.
    nMin : float, array
        Lower bound on the expected number of events per intensity in the
        sifted X/Z basis.
    nPlus : float, array
        Upper bound on the expected number of events per intensity in the
        sifted X/Z basis.
    s0 : float, optional
        The number of vacuum events in the sifted X or Z basis. 
        The default is None.

    Returns
    -------
    float
        The number of single photon events in the sifted X or Z basis.

    """
    if (s0):
        # Use the value of s0 provided
        return tau(1,mu,P)*mu[0]* ( nMin[1] - (mu[1]**2 / mu[0]**2)*nPlus[0] - \
                                   (mu[0]**2 - mu[1]**2) / mu[0]**2 * s0/tau(0,mu,P) ) \
                                   / (mu[1] *(mu[0]-mu[1]) )
    else:
        # Calculate s0
        return tau(1,mu,P)*mu[0]* ( nMin[1] - (mu[1]**2 / mu[0]**2)*nPlus[0] - \
                                   (mu[0]**2 - mu[1]**2) / mu[0]**2 * s0(mu,P,nMin)/tau(0,mu,P) ) \
                                   / (mu[1] *(mu[0]-mu[1]) )

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def vxz1(mu,P,mXZmin,mXZplus):
    """
    Calculates the upper bound to the number of bit errors associated with 
    single photon events in the sifted X or Z basis.
    See Eq. (4) in [1].

    Parameters
    ----------
    mu : float, array
        Pulse intensities.
    P : float, array
        Probability of Alice preparing a pulse intensity.
    mXZmin : float, array
        Lower bound on the expected number of errors per intensity in the
        sifted X/Z basis.
    mXZplus : float, array
        Upper bound on the expected number of errors per intensity in the
        sifted X/Z basis.

    Returns
    -------
    float
        Upper bound to the number of bit errors associated with single photon 
        events in the sifted X/Z basis.

    """
    return tau(1,mu,P)*(mXZplus[0] - mXZmin[1]) / (mu[0] - mu[1])
        
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def mean_photon_a(P,mu):
    """
    Calculate the mean photon number for a signal sent by Alice.
    This function uses arrays.

    Parameters
    ----------
    P : float, array
        Probability Alice sends a signal intensity.
    mu : float, array
        Intensity of the pulses.

    Returns
    -------
    float
        Mean siganl photon number.

    """
    return np.dot(P, mu)
        
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def Qber(PA,PB,nZ,mX,nXZ,nZX):
    paz = 1-PA
    pbz = 1-PB
    cb = mX/(PA*PB) - nXZ/(paz*PB) - nZX/(PA*pbz) + 2*nZ/(paz*PB)
    return min(0.5*(paz*pbz/nZ)*( mX/(PA*PB) + max(0, cb) ),1)