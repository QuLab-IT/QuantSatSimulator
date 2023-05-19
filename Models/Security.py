import numpy as np
from scipy.stats import binom
from sys import float_info

# Numerical value to use when denominator values are potentially zero
# Extract the smallest float that the current system can:
ZERO = float_info.epsilon # round (relative error due to rounding)
# Extract the largest float that the current system can represent
INF = float_info.max

MU3 = 0.0

class protocol():

    def __init__(self, states , decoys , parameters):
        self.states = states
        self.decoys = decoys
        if not paramcheck(decoys, parameters):
            exit(1)
        self.params = parameters
        self.boundM = 'chernoff'
        self.ECfunc = 'logm'
        self.eta = 0
        self.dt = 0

    def nstate(self,n):
        self.states = n

    def ndecoy(self,n):
        self.decoys = n

    def param(self,parameters):
        if paramcheck(parameters):
            self.params = parameters

    def bound_method(self,method):
        if method in ['Chernoff','chernoff','Hoeffding','hoeffding','Asymptotic','asymptotic']:
            self.boundM = method
        else:
            return

    def error_func(self,func):
        if func in ['logM','logm','Block','block','mX','mx','None','none']:
            self.ECfunc = func
        else:
            return

    def System_state(self,eta,dt):
        self.eta = eta
        self.dt = dt

    def key_length(self,var):
        """
        Returns the secure key length for an asymmetric BB84 protocol with weak
        coherent pulses for the specified States-Decoy number.

        Parameters
        ----------
        var : float, array/tuple
            var[0] = Asymmetric basis choice probability - Pax
            var[1] = Asymmetric basis choice probability - Pbx
            var[2] = Weak coherent pulse 1 probability - pk_1
            var[3] = Weak coherent pulse 2 probability - pk_2
            var[4] = Weak coherent pulse 1 intensity - mu_1
            var[5] = Weak coherent pulse 2 intensity - mu_2
        self.params : float, array/tuple
            'pax' = Alice's asymmetric basis choice probability
            'pbx' = Bob's asymmetric basis choice probability
            'pk1' = Weak coherent pulse 1 probability
            'pk2' = Weak coherent pulse 2 probability
            'mu1' = Weak coherent pulse 1 intensity
            'mu2' = Weak coherent pulse 2 intensity
            'mu3' = Weak coherent pulse 3 intensity
            'QBERI' = Intrinsic Quantum Bit Error Rate
            'pap' = Probability of an afterpulse event
            'pec' = Extraneous count probability
            'srate' = Alice's source rate (number of pulses per second)
            'eps_s' = The secrecy error; the key is eps_s secret
            'eps_c' = Correctness error in the secure key
        self.eta : float
            Total loss of the optical channel
        self.dt : float
            Time slot (in seconds)

        Returns
        -------
        l : float
            Secure key length (in bits).
        QBERx : float
            The Quantum Bit Error rate for the X basis.
        phi_x : float
            The phase error rate for the X basis.
        nX : float
            The total number of measured events in the X basis.
        nZ : float
            The total number of measured events in the Z basis.
        mXtot : float
            The total number of measurement errors in the X basis.
        lambdaEC : float
            The estimated number of bits used for error correction.
        sx0 : float
            The number of vacuum events in the X basis.
        sx1 : float
            The number of single-photon events in the X basis.
        vz1 : float
            The number of bit errors associated with single-photon events in the Z 
            basis.
        sz1 : float
            The number of single-photon events in the Z basis.
        mpn : float
            The mean photon number of the signal.
        
        """

        ###########################################################################
        ## Store input as named arrays/parameters for ease of use
        ###########################################################################
        # Asymmetric basis choice probability
        PAx = var[0]
        PBx = var[1]
        if self.decoys == 1:
            # Probability of Alice sending a pulse of intensity mu_k
            P = np.array([var[2],max(1-var[2],ZERO)])
            # Weak coherent pulse intensities
            mu = np.array([var[4],var[5]])
        elif self.decoys == 2:
            # Probability of Alice sending a pulse of intensity mu_k
            P = np.array([var[2],var[3],max(1-var[2]-var[3],ZERO)])
            # Weak coherent pulse intensities
            mu = np.array([var[4],var[5],MU3])
        
        ###########################################################################
        ## Estimate count and error rates for signals
        ###########################################################################
        # Exponential loss decay function
        loss_j = np.exp(-self.eta*mu)
        # Expected detection rate (including afterpulse contributions)
        # Returns the (3) array Dj
        Dj = Dk_j(self.params.get('pap'),self.params.get('pec'),loss_j)
        # Probability of having a bit error per intensity for a given time slot
        # Returns the (3) array ej
        ej = error_j(Dj,self.params.get('pap'),self.params.get('pec'),self.params.get('QBERI'),loss_j)

        ###########################################################################
        # Store some useful array products
        ###########################################################################
        # Take the dot product of the (3,) array P and the (3,) array ej.
        # Returns the (1,) array P_dot_ej
        P_dot_ej   = np.dot(P, ej)
        # Take the dot product of the (3,) array P and the (3,t) array Dj.
        # Returns the (1,) array P_dot_Dj
        P_dot_Dj   = np.dot(P, Dj)
        # Do an element-wise multiplication of the (3,) array P and the (3,) 
        # array Dj.
        # Returns the (3,) array P_times_Dj
        P_times_Dj = np.multiply(Dj, P)

        ###########################################################################
        ## Estimate count statistics
        ###########################################################################
        # Number of events in the sifted X basis for each time slot and intensity
        # Returns the (3,) array nx_j
        nx_mu = n_j(PAx, PBx, self.params.get('srate'), self.dt, P_times_Dj)
        # Number of events in the sifted Z basis for each time slot and intensity
        # Returns the (3,) array nz_j
        nz_mu = n_j(1 - PAx, 1 - PBx, self.params.get('srate'), self.dt, P_times_Dj)
        # Total number of events in the sifted X basis
        nX = np.sum(nx_mu)
        # Total number of events in the sifted Z basis
        nZ = np.sum(nz_mu)
        
        # Total number of errors in the sifted X and Z basis
        # Returns the (1,) array mX
        mX = m_j(nX, P_dot_Dj, P_dot_ej)
        mZ = m_j(nZ, P_dot_Dj, P_dot_ej)
        # Number of errors in the sifted X basis for each time slot and intensity
        # Returns the (3,) array mx
        mx_mu = mk_j(mX, P_dot_Dj, P_times_Dj)
        mz_mu = mk_j(mZ, P_dot_Dj, P_times_Dj)

        ###########################################################################
        ## Estimate bounds on count estimates due to statistical fluctuations
        ###########################################################################
        # Security of the protocol
        eps = self.params.get('eps_s') + self.params.get('eps_c')
        # Upper and lower bounds used to estimate the ideal number of X and Z basis
        # events accounting for statistical fluctuations  
        if self.boundM in ['Chernoff','chernoff']:
            # Use Chernoff bounds
            # Returns the (3,) arrays nXmin, nXplus
            nXmin, nXplus = CHbound(mu,P,nx_mu,eps)
            # Returns the (3,) arrays nZmin, nZplus
            nZmin, nZplus = CHbound(mu,P,nz_mu,eps)
            # Returns the (3,) arrays mZmin and mZplus
            mZmin, mZplus = CHbound(mu,P,mz_mu,eps)
            if self.decoys == 1:
                mXmin, mXplus = CHbound(mu,P,mx_mu,eps)
                mXmax = np.sum(mXplus)
                mZmax = np.sum(mZplus)
        elif self.boundM in ['Hoeffding','hoeffding']:
            # Use Hoeffding bound
            # Returns the (3,) arrays nXmin, nXplus
            nXmin, nXplus = HObound(mu,P,nx_mu,nX,eps)
            # Returns the (3,) arrays nZmin, nZplus
            nZmin, nZplus = HObound(mu,P,nz_mu,nZ,eps)
            # Returns the (3,) arrays mZmin and mZplus
            mZmin, mZplus = HObound(mu,P,mz_mu,mZ,eps)
            if self.decoys == 1:
                mXmin, mXplus = HObound(mu,P,mx_mu,mX,eps)
                mXmax = np.sum(mXplus)
                mZmax = np.sum(mZplus)
        elif self.boundM in ['Asymptotic','asymptotic']:
            # Use asymptotic bounds - no bounds
            # Returns the (3,) arrays nXmin, nXplus
            nXmin, nXplus = ASbound(mu,P,nx_mu)
            # Returns the (3,) arrays nZmin, nZplus
            nZmin, nZplus = ASbound(mu,P,nz_mu)
            # Returns the (3,) arrays mZmin and mZplus
            mZmin, mZplus = ASbound(mu,P,mz_mu)
            if self.decoys == 1:
                mXmin, mXplus = ASbound(mu,P,mx_mu)
                mXmax = np.sum(mXplus)
                mZmax = np.sum(mZplus)
        else:
            # Use Chernoff bounds
            # Returns the (3,) arrays nXmin, nXplus
            nXmin, nXplus = CHbound(mu,P,nx_mu,eps)
            # Returns the (3,) arrays nZmin, nZplus
            nZmin, nZplus = CHbound(mu,P,nz_mu,eps)
            # Returns the (3,) arrays mZmin and mZplus
            mZmin, mZplus = CHbound(mu,P,mz_mu,eps)
            if self.decoys == 1:
                mXmin, mXplus = CHbound(mu,P,mx_mu,eps)
                mXmax = np.sum(mXplus)
                mZmax = np.sum(mZplus)
            
        ###########################################################################
        ## Calculate the number of n-photon events
        ###########################################################################
        # Estimated number of n-photon for X and Z basis
        if self.decoys == 1:
            # For 1 decoy
            # Number of vacuum events in the sifted X basis
            sx0 = max(s0_1D(mu,P,nXmin,nXplus), ZERO)
            # Number of vacuum events in the sifted Z basis
            sz0 = max(s0_1D(mu,P,nZmin,nZplus), ZERO)

            if self.boundM in ['Chernoff','chernoff']:
                # Number of single photon events in the sifted X basis
                sx1 = max(s1_1D(mu,P,nXmin,nXplus,CHs0(mu,P,mXmax,nX,eps)), ZERO)
                # Number of single photon events in the sifted Z basis
                sz1 = max(s1_1D(mu,P,nZmin,nZplus,CHs0(mu,P,mZmax,nZ,eps)), ZERO)
            elif self.boundM in ['Hoeffding','hoeffding']:
                # Number of single photon events in the sifted X basis
                sx1 = max(s1_1D(mu,P,nXmin,nXplus,HOs0(mu,P,mXmax,nX,eps)), ZERO)
                # Number of single photon events in the sifted Z basis
                sz1 = max(s1_1D(mu,P,nZmin,nZplus,HOs0(mu,P,mZmax,nZ,eps)), ZERO)
            elif self.boundM in ['Asymptotic','asymptotic']:
                # Number of single photon events in the sifted X basis
                sx1 = max(s1_1D(mu,P,nXmin,nXplus,ASs0(mu,P,mXmax)), ZERO)
                # Number of single photon events in the sifted Z basis
                sz1 = max(s1_1D(mu,P,nZmin,nZplus,ASs0(mu,P,mZmax)), ZERO)
            else:
                # Number of single photon events in the sifted X basis
                sx1 = max(s1_1D(mu,P,nXmin,nXplus,CHs0(mu,P,mXmax,nx_mu,eps)), ZERO)
                # Number of single photon events in the sifted Z basis
                sz1 = max(s1_1D(mu,P,nZmin,nZplus,CHs0(mu,P,mZmax,nz_mu,eps)), ZERO)

        elif self.decoys == 2:
            # For 2 decoy
            # Number of vacuum events in the sifted X basis
            sx0 = max(s0_2D(mu,P,nXmin,nXplus), ZERO)
            # Number of vacuum events in the sifted Z basis
            sz0 = max(s0_2D(mu,P,nZmin,nZplus), ZERO)

            # Number of single photon events in the sifted X basis
            sx1 = max(s1_2D(mu,P,nXmin,nXplus,sx0), ZERO)
            # Number of single photon events in the sifted Z basis
            sz1 = max(s1_2D(mu,P,nZmin,nZplus,sz0), ZERO)

        ###########################################################################
        ## Calculate metrics such as the error rate, QBER and mean photon number
        ###########################################################################
        # Number of bit errors associated with single photon events in the sifted
        # Z basis 
        if self.decoys == 1:
            vz1   = min(max(v1_1D(mu,P,mZmin,mZplus), ZERO), mZ)
        elif self.decoys == 2:
            vz1   = min(max(v1_2D(mu,P,mZmin,mZplus), ZERO), mZ)
        # Ratio between the number of bit errors and the number of events
        # associated with single photons in the sifted Z basis.
        ratio = min(vz1 / sz1, 1 - ZERO)
        # The quantum bit error rate in the sifted X basis
        if self.states == 3:
            QBERx = QBERsimple()
        elif self.states == 4:
            if nX == 0:
                QBERx = 0
            else:
                QBERx = mX / nX
        # Calculate the mean photon number
        mpn   = mean_photon_a(P,mu)

        ###########################################################################
        ## Estimate the number of bits sacrificed for error correction
        ###########################################################################
        if self.ECfunc in ['logM','logm']:
            lambdaEC = logM(nX, QBERx, self.params.get('eps_c'))
        elif self.ECfunc in ['Block','block']:
            lambdaEC = 1.16 * nX * h(QBERx)
        elif self.ECfunc in ['mXtot','mxtot']:
            lambdaEC = 1.16 * mX
        elif self.ECfunc in ['None','none']:
            lambdaEC = 0
        else:
            lambdaEC = 0

        ###########################################################################
        ## Calculate the approximate length of the secret key per pass
        ###########################################################################
        if self.boundM in ['Asymptotic','asymptotic']:
            # Asymptotic phase error rate for single photon events in the sifted X
            # basis
            phi_x = min(ratio, 0.5)
            # Secret key length in the asymptotic regime
            l = max((sx0 + sx1 * (1 - h(phi_x)) - lambdaEC) , 0.0)
        else:
            # Phase error rate for single photon events in the sifted X basis
            phi_x = min(ratio + gamma(self.params.get('eps_s'),ratio,sz1,sx1), 0.5)
            if self.decoys == 1:
                a = 6
                b = 19
            elif self.decoys == 2:
                a = 6
                b = 21
            # Secret key length in the finite regime
            l = max(sx0 + sx1 * (1 - h(phi_x)) - lambdaEC -
                a*np.math.log(b / self.params.get('eps_s'), 2) - np.math.log(2.0 / self.params.get('eps_c'), 2) , 0.0)

        return l, QBERx, phi_x, nX, nZ, mX, lambdaEC, sx0, sx1, vz1, sz1, mpn, Dj, ej

    def key_length_inv(self,var):
        """
        Returns the inverse of the secure key length for an asymmetric BB84 
        protocol with weak coherent pulses and 2 'decoy' states. The intensity of 
        the weak coherent pulse 3, mu_3, is assumed to be a pre-defined global 
        parameter.

        Parameters
        ----------
        var : float, array/tuple
            var[0] = Asymmetric basis choice probability - Pax
            var[1] = Asymmetric basis choice probability - Pbx
            var[2] = Weak coherent pulse 1 probability - pk_1
            var[3] = Weak coherent pulse 2 probability - pk_2
            var[4] = Weak coherent pulse 1 intensity - mu_1
            var[5] = Weak coherent pulse 2 intensity - mu_2

        Returns
        -------
        1/l : float
            Inverse of the secure key length (in bits).

        """
        
        # Calculate the secure key length (normalised by NoPass)
        l, _, _, _, _, _, _, _, _, _, _, _, _, _ = self.key_length(var)

        # Safety check that all parameters are positive
        if (np.any(var[var < 0])):
            return INF  # ~1/0
        # Safety check that the parameters satisfy the constraints
        C = bool_constraints(self.decoys,var[0],var[1],var[2],var[3],var[4],var[5],MU3)
        if (not np.all(C)):
            return INF  # ~1/0

        if (l > 0):
            return (1.0/l)  # Inverse key length
        elif (l == 0):
            return INF  # ~1/0
        else:
            #return num_max  # ~1/0
            print("Error! Unexpected key length:", l)
            #return l        # Negative key length, NaN? --Useful for troubleshooting
            return INF  # ~1/0 --Stops calculations grinding to a halt

###########################################################################
## System Parameter Dependence Functions
###########################################################################

def Dk_j(Pap,Pec,dLoss):
    """
    Calculates the expected detection rate including afterpulse contributions
    for each intensity and time slot.
    ----------
    Pap : float
        Probability of an afterpulse event.
    Pec : float
        Extraneous count probability.
    dLoss : float, array
        Loss, per intensity per time slot, decay function.

    Returns
    -------
    float, array
        Expected detection rate.

    """
    return (1 + Pap)*(1 - (1 - 2*Pec)*dLoss)

def error_j(Dj,Pap,Pec,QBERI,dLoss):
    """
    Calculates the conditional probability for a pulse of intensity mu_j
    to cause an error, after sifting, in the time slot t.
    ----------
    Dj : float, array
        Expected detection rate.
    Pap : float
        Probability of an afterpulse event.
    Pec : float
        Extraneous count probability.
    QBERI : float
        Intrinsic Quantum Bit Error Rate.
    dLoss : float, array
        Loss, per intensity per time slot, decay function.

    Returns
    -------
    float, array
        Error rate per intensity per time slot.

    """
    return Pec + (0.5*Pap*Dj) + QBERI*(1 - dLoss)

def n_j(PA,PB,Fs,dt,P_times_Dj):
    """
    Calculates the number of events in the X or Z sifted basis per pulse
    intensity per time slot.
    ----------
    PA : float
        Probability of Alice preparing a state in the X/Z basis.
    PB : float
        Probability of Bob measuring a state in the X/Z basis.
    Fs : integer/float
        Number of pulses sent by Alice.
    P_times_Dj : float, array
        The element-wise multiplication of the intensity probability array P
        with the expected detection rate per time slot array Dj.

    Returns
    -------
    float, array
        The number of events in the sifted X/Z basis.

    """
    return PA*PB*Fs*dt*P_times_Dj

def m_j(n,P_dot_Dj,P_dot_ej):
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

###########################################################################
## Statistical Fluctuations Bound Functions
###########################################################################

def CHbound(mu,P,x_mu,eps):
    """
    Calculates the upper and lower bounds  for each pulse intensity using, the Chernoff bounds.
    ----------
    mu : float, array
        Pulse intensities.
    P : float, array
        Probability of Alice preparing a pulse intensity.
    x_mu : float, array
        Number of events per intensity for the sifted X/Z basis.
    eps_s : float
        The secrecy error; the key is eps_s secret.

    Returns
    -------
    Xmin : float, array
        Lower bound on the expected number of events per intensity in the
        sifted X/Z basis.
    Xplus : float, array
        Upper bound on the expected number of events per intensity in the
        sifted X/Z basis.

    """
    log_es = np.math.log(1 / eps)
    term_m  = 0.5*log_es + np.sqrt(2*x_mu*log_es + 0.25*log_es**2)
    term_p  = log_es + np.sqrt(2*x_mu*log_es + log_es**2)
    aux_P = np.array([max(ZERO,x) for x in P])
    Xmin  = np.divide(np.multiply(np.exp(mu), x_mu - term_m), aux_P)
    Xplus = np.divide(np.multiply(np.exp(mu), x_mu + term_p), aux_P)
    return Xmin, Xplus

def HObound(mu,P,x_mu,X,eps):
    """
    Calculates the upper and lower bounds for each pulse intensity using, the Hoeffding bound.
    ----------
    mu : float, array
        Pulse intensities.
    P : float, array
        Probability of Alice preparing a pulse intensity.
    x_mu : float, array
        Number of events per intensity for the sifted X/Z basis.
    X : float, array
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
    term   = np.sqrt(0.5*X * np.math.log(1 / eps))
    Xmin  = np.divide(np.multiply(np.exp(mu), x_mu - term), P)
    Xplus = np.divide(np.multiply(np.exp(mu), x_mu + term), P)
    return Xmin, Xplus

def ASbound(mu,P,x_mu):
    """
    Calculates the number of events in the  sifted X or Z basis, for each pulse 
    intensity in the asymptotic limit.
    ----------
    mu : float, array
        Pulse intensities.
    P : float, array
        Probability of Alice preparing a pulse intensity.
    x_mu : float, array
        Number of events per intensity for the sifted X/Z basis.

    Returns
    -------
    X : float, array
        The expected number of events per intensity in the
        sifted X/Z basis.

    """
    X = np.divide(np.multiply(np.exp(mu), x_mu), P)
    return X, X

###########################################################################
## Number of n-photon events Functions
###########################################################################

def tau(n,mu,P):
    """
    Calculates the total probability that Alice prepares an n-photon state.
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
    tau = 0
    for k in range(len(mu)):
        tau += np.math.exp(-mu[k]) * mu[k]**n * P[k]
    return tau / np.math.factorial(n)

def s0_2D(mu,P,nMin,nPlus):
    """
    Calculates the approximate number of vacuum events in the sifted X or Z 
    basis.
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
    return tau(0,mu,P) * (mu[1]*nMin[2] - mu[2]*nPlus[1]) / (mu[1] - mu[2])

def s1_2D(mu,P,nMin,nPlus,s0=None):
    """
    Calculates the number of single photon events in the sifted X or Z basis.
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
    if s0 is not None:
        # Use the value of s0 provided
        return tau(1,mu,P)*mu[0]* ( nMin[1] - nPlus[2] - \
                                   (mu[1]**2 - mu[2]**2) / mu[0]**2 * \
                                   (nPlus[0] - s0 / tau(0,mu,P)) ) / \
                            (mu[0] * (mu[1] - mu[2]) - mu[1]**2 + mu[2]**2)
    else:
        # Calculate s0
        return tau(1,mu,P)*mu[0]* ( nMin[1] - nPlus[2] - \
                                   (mu[1]**2 - mu[2]**2) / mu[0]**2 * \
                                   (nPlus[0] - s0_2D(mu,P,nMin,nPlus) / tau(0,mu,P)) ) / \
                            (mu[0] * (mu[1] - mu[2]) - mu[1]**2 + mu[2]**2)

def s0_1D(mu,P,nMin,nPlus):
    """
    Calculates the approximate number of vacuum events in the sifted X or Z 
    basis.
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
    return tau(0,mu,P) * (mu[0]*nMin[1] - mu[1]*nPlus[0]) / (mu[0] - mu[1])

def CHs0(mu,P,mMax,n_mu,eps):
    """
    Calculates the approximate number of vacuum events in the sifted X or Z 
    basis, using the Chernoff bound.
    ----------
    mu : float, array
        Pulse intensities.
    P : float, array
        Probabilities that Alice prepares a particular intensity.
    nMax : float, array
        Upper bound on the expected number of events per intensity in the
        sifted X/Z basis.
    n_mu : float, array
        Number of events per intensity for the sifted X/Z basis.
    eps_s : float
        The secrecy error; the key is eps_s secret.

    Returns
    -------
    float
        The number of vacuum events in the sifted X or Z basis.

    """
    log_es = np.math.log(1 / eps)
    term  = log_es + np.sqrt(2*n_mu*log_es + log_es**2)
    return 2*(tau(0,mu,P)*mMax + term)

def HOs0(mu,P,mMax,X,eps):
    """
    Calculates the approximate number of vacuum events in the sifted X or Z 
    basis, using the Hoeffding bound.
    ----------
    mu : float, array
        Pulse intensities.
    P : float, array
        Probabilities that Alice prepares a particular intensity.
    nMin : float, array
        Lower bound on the expected number of events per intensity in the
        sifted X/Z basis.
    X : float, array
        Number of events for the sifted X/Z basis.
    eps_s : float
        The secrecy error; the key is eps_s secret.

    Returns
    -------
    float
        The number of vacuum events in the sifted X or Z basis.

    """
    term   = np.sqrt(0.5*X * np.math.log(1 / eps))
    return 2*(tau(0,mu,P)*mMax + term)

def ASs0(mu,P,mMax):
    """
    Calculates the approximate number of vacuum events in the sifted X or Z 
    basis, using the Asymptotic case.
    ----------
    mu : float, array
        Pulse intensities.
    P : float, array
        Probabilities that Alice prepares a particular intensity.
    nMax : float, array
        expected number of events per intensity in the
        sifted X/Z basis.

    Returns
    -------
    float
        The number of vacuum events in the sifted X or Z basis.

    """
    return 2*(tau(0,mu,P)*mMax)

def s1_1D(mu,P,nMin,nPlus,s0):
    """
    Calculates the number of single photon events in the sifted X or Z basis.
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
    s0 : float
        The number of vacuum events in the sifted X or Z basis. 
        The default is None.

    Returns
    -------
    float
        The number of single photon events in the sifted X or Z basis.

    """
    
    return tau(1,mu,P)*mu[0]* ( nMin[1] - (mu[1]**2 / mu[0]**2)*nPlus[0] - \
                                (mu[0]**2 - mu[1]**2) / mu[0]**2 * \
                                s0 / tau(0,mu,P) ) / (mu[1] *(mu[0]-mu[1]) )

###########################################################################
## Quantum Bit Error Rate Functions
###########################################################################

def v1_2D(mu,P,mMin,mPlus):
    """
    Calculates the upper bound to the number of bit errors associated with 
    single photon events in the sifted X or Z basis.
    ----------
    mu : float, array
        Pulse intensities.
    P : float, array
        Probability of Alice preparing a pulse intensity.
    mMin : float, array
        Lower bound on the expected number of errors per intensity in the
        sifted X/Z basis.
    mPlus : float, array
        Upper bound on the expected number of errors per intensity in the
        sifted X/Z basis.

    Returns
    -------
    float
        Upper bound to the number of bit errors associated with single photon 
        events in the sifted X/Z basis.

    """
    return tau(1,mu,P)*(mPlus[1] - mMin[2]) / (mu[1] - mu[2])

def v1_1D(mu,P,mMin,mPlus):
    """
    Calculates the upper bound to the number of bit errors associated with 
    single photon events in the sifted X or Z basis.
    ----------
    mu : float, array
        Pulse intensities.
    P : float, array
        Probability of Alice preparing a pulse intensity.
    mMin : float, array
        Lower bound on the expected number of errors per intensity in the
        sifted X/Z basis.
    mPlus : float, array
        Upper bound on the expected number of errors per intensity in the
        sifted X/Z basis.

    Returns
    -------
    float
        Upper bound to the number of bit errors associated with single photon 
        events in the sifted X/Z basis.

    """
    return tau(1,mu,P)*(mPlus[0] - mMin[1]) / (mu[0] - mu[1])

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

def QBERsimple():
    return 0

def logM(nX, QBERx, eps_c):
    """
    This function os used to approximate the amount of secure key that is
    required to perform the error correction in bits.
    See Eq. (6) in [3].

    Parameters
    ----------
    nX : float
        The total number of measured events in the sifted X basis.
    QBERx : float
        The quantum bit error rate in the sifted X basis.
    eps_c : float
        The correctness error.

    Returns
    -------
    lM : float
        Number of bits for error correction.

    """
    if QBERx == 1:
        lM = nX * h(QBERx) + (nX * (1.0 - QBERx) - FInv(int(nX), QBERx, eps_c) - 1) \
             * (-INF) - 0.5*np.math.log(nX) - np.math.log(1.0 / eps_c)
    elif QBERx == 0:
        lM = 0.0
    else:
        lM = nX * h(QBERx) + (nX * (1.0 - QBERx) - FInv(int(nX), QBERx, eps_c) - 1) \
             * np.math.log((1.0 - QBERx) / QBERx) - 0.5*np.math.log(nX) - \
               np.math.log(1.0 / eps_c)
    return lM

def FInv(nX, QBERx, eps_c):
    """
    Calculates the quantile function, or inverse cumulative distribution
    function for a binomial distribution, nCk p^k (1 - p)^(n-k), where
    n = floor(nX), p = 1-QBERx, k \propto eps_c

    Parameters
    ----------
    nX : float/integer
        Number of events in sifted X basis.
    QBERx : float
        Quantum bit error rate in sifted X basis.
    eps_c : float
        Correctness error in the secure key.

    Returns
    -------
    float
        Quartile function of binomial distribution.

    """
    return binom.ppf(eps_c * (1. + 1.0 / np.sqrt(nX)), int(nX), 1.0 - QBERx)

def gamma(a,b,c,d):
    """
    Gamma Function
    ----------
    a : float
        Argument 1.
    b : float
        Argument 2.
    c : float
        Argument 3.
    d : float
        Argument 4.

    Returns
    -------
    g : float
        Output value.

    """
    g1 = max((c + d) * (1 - b) * b / (c*d * np.math.log(2)), 0.0)
    g2 = max((c + d) * 21**2 / (c*d * (1 - b) * b*a**2), 1.0)
    g  = np.math.sqrt(g1 * np.math.log(g2, 2))
    return g

def h(x):
    """
    Binary entropy function.
    ----------
    x : float
        Function argument.

    Returns
    -------
    h : float
        Binary entropy.

    """
    if x == 0 or x ==1:
        h = 0
    else:
        h = -x*np.math.log(x, 2) - (1 - x)*np.math.log(1 - x, 2)
    return h

###########################################################################
## Other Functions
###########################################################################

def paramcheck(decoys,param):
    """
    Parameters
    ----------
    param : float, array/tuple
        'pax' = Alice's asymmetric basis choice probability
        'pbx' = Bob's asymmetric basis choice probability
        'pk1' = Weak coherent pulse 1 probability
        'pk2' = Weak coherent pulse 2 probability
        'pk3' = Weak coherent pulse 3 probability
        'mu1' = Weak coherent pulse 1 intensity
        'mu2' = Weak coherent pulse 2 intensity
        'mu3' = Weak coherent pulse 2 intensity
        'QBERI' = Intrinsic Quantum Bit Error Rate
        'pap' = Probability of an afterpulse event
        'pec' = Extraneous count probability
        'srate' = Alice's source rate (number of pulses per second)
    """
    flag = True

    if not isinstance(param,dict):
        print('Error! Parameters are not a dictionary')
        return flag

    if decoys == 1:
        C = ['pax','pbx','pk1','pk2','mu1','mu2','QBERI','pap','pec']
    else:
        C = ['pax','pbx','pk1','pk2','pk3','mu1','mu2','mu3','QBERI','pap','pec']
    
    for p in C:
        if p in param:
            var = param.get(p)
            if isinstance(var,float):
                if (var > 1.0 or var < 0.0):
                    print('Error! Constraint 0 < {} < 1: {}'.format(p,var))
                    flag = False
            else:
                print('Error! Type of parameter {} is not a float'.format(p))
                flag = False

    if 'srate' in param:
        if param.get('srate') < 0:
            print('Error! Constraint srate > 0: {}'.format(param.get('srate')))
            flag = False

    if decoys == 1:
        if param.get('pk1') > 1:
            print("Error! Constraint pk_1 < 1: ", param.get('pk1'))
            flag = False
        if (param.get('mu1') <= param.get('mu2')):
            print("Error! Constraint mu1 > mu2: ", param.get('mu1'), param.get('mu2'))
            flag = False
    else:
        if param.get('pk1') + param.get('pk2') > 1:
            print("Error! Constraint pk_1 + pk_2 < 1: ", param.get('pk1') + param.get('pk2'))
            flag = False
        if ((param.get('mu1') - MU3) <= param.get('mu2')):
            print("Error! Constraint (mu1-mu3) > mu2: ", (param.get('mu1') - MU3), param.get('mu2'))
            flag = False
        if (param.get('mu2') <= MU3):
            print("Error! Constraint mu2 > mu3: ", param.get('mu2'), MU3)
            flag = False

    return flag

def bool_constraints(decoys,Pax,Pbx,pk1,pk2,mu1,mu2,mu3):

    C = np.array([1,1,1,1,1,1,1,1,1], dtype=bool) # Initialise array as True

    if decoys > 1:
        # Constraint 1: Check polarisation basis probabilities are valid.
        if (Pax >= 1.0 or Pax <= 0.0):
            C[0] = False
        # Constraint 2: Check polarisation basis probabilities are valid.
        if (Pbx >= 1.0 or Pbx <= 0.0):
            C[1] = False
        # Constraint 2: Check probability of pulse with intensity 1 is in bounds.
        if (pk1 >= 1.0 or pk1 <= 0.0):
            C[2] = False
        # Constraint 3: Check probability of pulse with intensity 2 is in bounds.
        if (pk2 >= 1.0 or pk2 <= 0.0):
            C[3] = False
        # Constraint 4: Check sum of probabilities for intensity 1 & 2 are less
        # than unity.
        if ((pk1 + pk2) >= 1.0):
            C[4] = False
        # Constraint 5: Check value of intensity 1 is in bounds.
        if (mu1 >= 1.0 or mu1 <= 0.0):
            C[5] = False
        # Constraint 6: Check value of intensity 2 is in bounds.
        if (mu2 >= 1.0 or mu2 <= 0.0):
            C[6] = False
        # Constraint 7: Check values of all intensities are in bounds.
        if ((mu1 - mu3) <= mu2):
            C[7] = False
        # Constraint 8: Check values of intensities 2 & 3 are in bounds.
        if (mu2 <= mu3):
            C[8] = False
    else:
        # Constraint 1: Check polarisation basis probabilities are valid.
        if (Pax >= 1.0 or Pax <= 0.0):
            C[0] = False
        # Constraint 2: Check polarisation basis probabilities are valid.
        if (Pbx >= 1.0 or Pbx <= 0.0):
            C[1] = False
        # Constraint 2: Check probability of pulse with intensity 1 is in bounds.
        if (pk1 >= 1.0 or pk1 <= 0.0):
            C[2] = False
        # Constraint 5: Check value of intensity 1 is in bounds.
        if (mu1 >= 1.0 or mu1 <= 0.0):
            C[5] = False
        # Constraint 6: Check value of intensity 2 is in bounds.
        if (mu2 >= 1.0 or mu2 <= 0.0):
            C[6] = False
        # Constraint 7: Check values of all intensities are in bounds.
        if (mu1 <= mu2):
            C[7] = False
    return C