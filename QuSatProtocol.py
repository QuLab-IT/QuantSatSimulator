import numpy as np
import Security as sc
import OptSet as opts
import os
import warnings
from tqdm import tqdm
from scipy.optimize import minimize
from sys import float_info
from time import perf_counter

#Directories:
CD = os.getcwd()
DDATA = os.path.join(CD,"Data")
OUTFILE = 'Pedro'

# Physical Constants:
RT = 6378000            # Earth Radius in [m]

# Numerical value to use when denominator values are potentially zero
# Extract the smallest float that the current system can:
ZERO = float_info.epsilon # round (relative error due to rounding)
# Extract the largest float that the current system can represent
INF = float_info.max

MU3 = 0.0

# Micius Data:
M_H = 600000
M_DT = 0.5
M_DR = 0.3
M_W = 0.89*M_DT
M_WL = 785*10**(-9)

M_DETCLOSS = 2.3
M_OPTLOSS = 3
M_INT = M_DETCLOSS + M_OPTLOSS

F_or_T = [False, True] # List used to switch between False or True values

#warnings.filterwarnings("ignore") # Surpress printing warnings

# Printing Options
BAR_SIZE = 60

class Parameters():
    def __init__(self):
        self.PParam = { 'pax':   None,
                        'pbx':   None,
                        'pk1':   None,
                        'pk2':   None,
                        'mu1':   None,
                        'mu2':   None,
                        'mu3':   None,
                        'QBERI': None,
                        'pap':   None,
                        'pec':   None,
                        'srate': None,
                        'eps_s': None,
                        'eps_c': None }

        self.SParam = { 'states': None,
                        'decoys': None,
                        'h':      None,
                        'dt':     None,
                        'dr':     None,
                        'w':      None,
                        'wl':     None }

        self.AdvParam = { 'boundM': None,
                          'ECfunc': None,
                          'intr':   None,
                          'optM':   None,
                          'optP':   None,
                          'R_pax':  None,
                          'R_pbx':  None,
                          'R_pk1':  None,
                          'R_pk2':  None,
                          'R_mu1':  None,
                          'R_mu2':  None,
                          'R_mu3':  None }

    def print(self,sparam = True, pparam = True, advparam = True):
        if sparam:
            print('System Parameters:')
            for item in list(self.SParam.items()):
                print(' '*2 + '> {}: {}'.format(item[0],item[1]))
        if pparam:
            print('Protocol Parameters:')
            for item in list(self.PParam.items()):
                print(' '*2 + '> {}: {}'.format(item[0],item[1]))
        if advparam:
            print('Advanced Parameters:')
            for item in list(self.AdvParam.items()):
                print(' '*2 + '> {}: {}'.format(item[0],item[1]))

    def check (self):
        x = [0,0,0]
        if None in self.PParam.values():
            x[0] = -1
        elif None in self.SParam.values():
            x[1] = -1
        elif None in self.AdvParam.values():
            x[2] = -1
        
        return x
        
    def get(self):
        if -1 in self.check():
            return None, None, None
        else:
            return self.SParam, self.PParam, self.AdvParam

    def PAx (self, value, range = None):
        if not ratio_check(value):
            return -1
        self.PParam['pax'] = value

        if range is not None:
            if not (ratio_check(range[0]) and ratio_check(range[1])):
                return -1
            elif range[0] > range[1]:
                return -1
            elif range[0] > value or range[1] < value:
                return -1
            self.AdvParam['R_pax'] = range
        else:
            self.AdvParam['R_pax'] = [value,value]
        
        return 0

    def PBx (self, value, range = None):
        if not ratio_check(value):
            return -1
        self.PParam['pbx'] = value

        if range is not None:
            if not (ratio_check(range[0]) and ratio_check(range[1])):
                return -1
            elif range[0] > range[1]:
                return -1
            elif range[0] > value or range[1] < value:
                return -1
            self.AdvParam['R_pbx'] = range
        else:
            self.AdvParam['R_pbx'] = [value,value]
            
        return 0

    def PK1 (self, value, range = None):
        if not ratio_check(value):
            return -1
        self.PParam['pk1'] = value

        if range is not None:
            if not (ratio_check(range[0]) and ratio_check(range[1])):
                return -1
            elif range[0] > range[1]:
                return -1
            elif range[0] > value or range[1] < value:
                return -1
            self.AdvParam['R_pk1'] = range
        else:
            self.AdvParam['R_pk1'] = [value,value]
        
        return 0

    def PK2 (self, value, range = None):
        if self.PParam['pk1'] is None:
            return -1
        elif not ratio_check(value):
            return -1
        elif self.PParam['pk1'] < value:
            return -1
        self.PParam['pk2'] = value

        if range is not None:
            if not (ratio_check(range[0]) and ratio_check(range[1])):
                return -1
            elif range[0] > range[1]:
                return -1
            elif range[0] > value or range[1] < value:
                return -1
            self.AdvParam['R_pk2'] = range
        else:
            self.AdvParam['R_pk2'] = [value,value]
        
        return 0

    def Mu1 (self, value, range = None):
        if not ratio_check(value):
            return -1
        self.PParam['mu1'] = value

        if range is not None:
            if not (ratio_check(range[0]) and ratio_check(range[1])):
                return -1
            elif range[0] > range[1]:
                return -1
            elif range[0] > value or range[1] < value:
                return -1
            self.AdvParam['R_mu1'] = range
        else:
            self.AdvParam['R_mu1'] = [value,value]
        
        return 0

    def Mu2 (self, value, range = None):
        if self.SParam['states'] is None:
            return -1
        if self.PParam['mu1'] is None:
            return -1
        elif not ratio_check(value):
            return -1
        elif self.PParam['mu1'] < value:
            return -1
        elif self.SParam['decoys'] < 1:
            return -1
        self.PParam['mu2'] = value

        if range is not None:
            if not (ratio_check(range[0]) and ratio_check(range[1])):
                return -1
            elif range[0] > range[1]:
                return -1
            elif range[0] > value or range[1] < value:
                return -1
            self.AdvParam['R_mu2'] = range
        else:
            self.AdvParam['R_mu2'] = [value,value]
        
        return 0

    def Mu3 (self, value, range = None):
        if self.SParam['states'] is None:
            return -1
        if self.PParam['mu2'] is None:
            return -1
        elif not ratio_check(value):
            return -1
        elif self.PParam['mu2'] < value:
            return -1
        elif self.PParam['mu1']-value <= self.PParam['mu2']:
            return -1
        self.PParam['mu3'] = value

        if range is not None:
            if not (ratio_check(range[0]) and ratio_check(range[1])):
                return -1
            elif range[0] > range[1]:
                return -1
            elif range[0] > value or range[1] < value:
                return -1
            self.AdvParam['R_mu3'] = range
        else:
            self.AdvParam['R_mu3'] = [value,value]
        
        return 0

    def QBERI (self, value):
        if not ratio_check(value):
            return -1
        self.PParam['QBERI'] = value
        return 0
    
    def Pap (self, value):
        if not ratio_check(value):
            return -1
        self.PParam['pap'] = value
        return 0

    def Pec (self, value):
        if not ratio_check(value):
            return -1
        self.PParam['pec'] = value
        return 0

    def SRate (self, value):
        if value < 0:
            return -1
        self.PParam['srate'] = value
        return 0

    def Eps_c (self, value):
        if not ratio_check(value):
            return -1
        self.PParam['eps_c'] = value
        return 0

    def Eps_s (self, value):
        if not ratio_check(value):
            return -1
        self.PParam['eps_s'] = value
        return 0

    def States (self,value):
        if value < 3 or value > 4:
            return -1
        self.SParam['states'] = value
        return 0

    def Decoys (self,value):
        if value < 1 or value > 2:
            return -1
        self.SParam['decoys'] = value
        return 0

    def Alt (self,value):
        if value < 0:
            return -1
        self.SParam['h'] = value
        return 0

    def DTransmitter (self,value):
        if value < 0:
            return -1
        self.SParam['dt'] = value
        return 0

    def DReceiver (self,value):
        if value < 0:
            return -1
        self.SParam['dr'] = value
        return 0

    def Waist (self,value):
        if value < 0:
            return -1
        self.SParam['w'] = value
        return 0

    def Wavelenght (self,value):
        if value < 0:
            return -1
        self.SParam['wl'] = value
        return 0

    def BoundM (self,value):
        if value in ['Chernoff','chernoff','Hoeffding','hoeffding','Asymptotic','asymptotic']:
            self.AdvParam['boundM'] = value
            return 0
        return -1

    def ECfunc (self,value):
        if value in ['logM','logm','Block','block','mX','mx','None','none']:
            self.AdvParam['ECfunc'] = value
            return 0
        return -1

    def Intr (self,value):
        if value < 0:
            return -1
        self.AdvParam['intr'] = value
        return 0

    def OptM (self,value):
        if value in ['COBYLA','SLSQP']:
            self.AdvParam['optM'] = value
            return 0
        return -1

    def OptP (self,value):
        if self.AdvParam['optM'] is None:
            return -1
        default = opts.opt_method(self.AdvParam['optM'])
        flag = 0
        for key in list(default.param_list().keys()):
            if key not in list(value.keys()):
                flag = -1
        if flag == -1:
            return -1
        self.AdvParam['optP'] = value
        return 0

def loadData(path,filename,usecols = None):
    return np.loadtxt(os.path.join(path,filename), delimiter=',', skiprows=1, usecols=usecols)

def writeDataCSV(data,outpath,outfile,out_head=None):
    """
    Write out data to a CSV file

    Parameters
    ----------
    data : float, array-like
        Data array containing parameters, SKL, and protocol metrics.
    outpath : string
        Path for output file.
    outfile : string
        Name for output file.
    out_head : string, optional
        Header for data file
    message : string, optional
        Data description for print command, default = 'data'.
    Returns
    -------
    None.

    """
    if (out_head is not None):
        nhead = len(out_head.split(',')) # Split header at every comma
        if (data.shape[1] != nhead):
            print('Warning: No. of fields does not match number of headings in', 
                  'output file:',outfile+'.csv')
            print('No. fields =',data.shape[1],', No. headings =',nhead)
    filename = os.path.join(outpath, outfile + '.csv')
    np.savetxt(filename,data,delimiter=',',header=out_head) 
    return None

def Opt_Param_find(eff,delta,protocol,AdvParam):

    decoys = protocol.decoys

    #******************************************************************************
    # Initialise calculation parameters
    #******************************************************************************
    # Store optimisation parameters
    method = AdvParam['optM']
    optparam = opts.opt_method(method)
    for item in list(AdvParam['optP'].items()):
        optparam.edit(item[0],item[1])
    O = opts.bounds(optparam)
    order = ['pax', 'pbx', 'pk1', 'pk2', 'mu1', 'mu2']
    for n in order:
        O.add_var(AdvParam['R_'+n],n)
    bounds,cons,options = O.opt_param(decoys)
    Nmax = O.method.param['Nmax']
    NoptMin = O.method.param['NoptMin']
    tStopBetter = O.method.param['tStopBetter']
    tStopZero = O.method.param['tStopZero']
    # Store initial/fixed protocol parameters as an array
    x0 = x0_rand(O.lowerb(),O.upperb(),ZERO,decoys)

    # Use the optimised parameters from the previous shift angle calculation
    protocol.System_state(eff,delta)
    res = minimize(protocol.key_length_inv,x0,args=(),method=method,
                    jac=None,hess=None,hessp=None,bounds=bounds, 
                    constraints=cons,tol=None,callback=None, 
                    options=options)
    Ntot = res.nfev # Initilaise total No. of function evaluations
    Nopt = 1        # Number of optimisation calls
    # Re-run optimization until Nmax function evaluations
    # have been used. Take a copy of initial results to compare.
    x0c  = x0
    SKLc = int(1.0 / res.fun)
    resc = res
    Nzero = 0 # Number of times we get SKL == 0
    while Ntot < Nmax or Nopt < NoptMin:
        Nopt += 1
        # Randomly assign new initial parameters
        x0 = x0_rand(O.lowerb(),O.upperb(),ZERO,decoys)
        # Calculate optimised SKL
        res = minimize(protocol.key_length_inv,x0,args=(),method=method,
                        jac=None,hess=None,hessp=None,bounds=bounds, 
                        constraints=cons,tol=None,callback=None, 
                        options=options)
        if int(1.0 / res.fun) > 0:
            if int(1.0 / res.fun) > SKLc:
                if Nopt >= NoptMin and tStopBetter:
                    break # A better value was found!
                else:
                    # Store new set of best parameters
                    x0c  = x0
                    resc = res
                    SKLc = int(1.0 / res.fun)
            else:
                # Reset to the best parameters
                x0  = x0c
                res = resc
        else:
            # SKL = 0. Reset to the 'best' parameters,
            # (may still give SKL = 0).
            Nzero += 1
            if Nopt > NoptMin:
                if Nzero / (Nopt - 1) == 1:
                    # We get SKL = 0 every time.
                    if tStopZero:
                        break
            x0  = x0c
            res = resc
        Ntot += res.nfev

    return res.x

def Fixed_Sim(SParam, PParam, AdvParam, Print=None, GUI=None, outfile=OUTFILE, opt=None):
    if None in [SParam, PParam, AdvParam]:
        if Print is not None:
            print("Error: Some parameters are not defined:")
            print('System Parameters:')
            for p in SParam.items():
                print(' '*2 + '> ' + p[0] + ': ' + str(p[1]),end="")
                if p[1] is None:
                    print(' (!)')
                else:
                    print('\n')
        return

    pass_data = loadData(os.path.join(DDATA,'Micius Data'), "Micius.txt", usecols=(0,1,2))
    time = pass_data[:,0]
    elevation = pass_data[:,1]

    states = SParam['states']
    decoys = SParam['decoys']
    P = sc.protocol(states,decoys,PParam)
    P.bound_method(AdvParam['boundM'])
    P.error_func(AdvParam['ECfunc'])

    h         = SParam['h']
    dtrans    = SParam['dt']
    dreceiv   = SParam['dr']
    w         = SParam['w']
    wl        = SParam['wl']

    # Losses
    dif  = diffraction(elevation, h, dtrans, dreceiv, w, wl)
    atm  = atmospheric(elevation, pass_data[:,2])
    intr = np.full(len(elevation), AdvParam['intr'])
    loss_dB = dif + atm + intr - 16
    eff = 10**( -loss_dB*0.1)

    #******************************************************************************
    # Initialise loop ranges
    #******************************************************************************
    Total_c = len(eff)
    #******************************************************************************
    # Initialise calculation parameters
    #******************************************************************************
    # Store initial/fixed protocol parameters as an array
    order = ['pax', 'pbx', 'pk1', 'pk2', 'mu1', 'mu2', 'mu3']
    if opt is not None:
        i = np.argmin(eff)
        x0 = Opt_Param_find(eff[i],np.abs(time[i]-time[i-1]),P,AdvParam)
    else:
        x0 = np.array([PParam[i] for i in order])
    #******************************************************************************
    # Initialise output data storage and headers
    #******************************************************************************
    # Header for CSV file: Data columns
    header = "SysLoss [dB],DiffLoss,AtmLoss,IntrLoss,time,dt,elev,SKL,QBERx,phiX,nX,nZ,mX,lambdaEC,sX0,sX1,vZ1,sZ1," + \
             "mean photon no.,D1,D2,D3,e1,e2,e3,QBERI,Pec,Pap,Rrate,eps_c,eps_s," + \
             "PAx,PBx,P1,P2,P3,mu1,mu2,mu3"
    # Initialise a data storage array: shape(No. of data runs, No. of metrics)   
    fulldata = np.empty((Total_c,len(header.split(","))))

    if GUI is not None:
        return
    elif Print is not None:
        print('_'*60,'\n')
        print('-'*60)
        print('QuSatProtocol v1.0')
        print('-'*60)
        print('Efficient BB84 security protocol')
        print('For {} States and {} Decoys'.format(P.states,P.decoys))
        print('Using Fixed protocol parameters: Px, pk, mu')
        print('Using {} bounds for statistical fluctuations'.format(P.boundM))
        if P.ECfunc in ['logM','logm']:
            print('Error correction model: logM(nX, QBERx, eps_c)')
        elif P.ECfunc in ['Block','block']:
            print('Error correction model: 1.16 * nX * h(QBERx)')
        elif P.ECfunc in ['mXtot','mxtot']:
            print('Error correction model: 1.16 * mXtot')
        elif P.ECfunc in ['None','none']:
            print('Error correction term not included')
        print('-'*60)
        print('_'*60,'\n')

    tc00 = perf_counter()       # Start clock timer
    tc = 0
    count = 0                   # Initialise calculation counter
    for i in tqdm(range(Total_c),ncols=BAR_SIZE,disable=not Print):

        if GUI is not None:
            return
        
        if i >= len(eff)-1 :
            delta = np.abs(time[i]-time[i-1])
        else:
            delta = np.abs(time[i+1]-time[i])
        
        # Use the optimised parameters from the previous shift angle calculation
        P.System_state(eff[count],delta)
        # Get final parameters from standard key length function
        SKL,QBERx,phi_x,nX,nZ,mX,lambdaEC,sX0,sX1,vz1,sZ1,mpn,Dj,ej = P.key_length(x0)
        # Store calculation parameters
        fulldata[count] = [loss_dB[count],dif[count],atm[count],intr[count],
                           time[count],delta,elevation[count],SKL,QBERx,phi_x,nX,nZ,mX,lambdaEC,sX0,sX1,
                           vz1,sZ1,mpn,Dj[0],Dj[1],dim_check(Dj,2),ej[0],ej[1],dim_check(ej,2),
                           PParam['QBERI'],PParam['pec'],PParam['pap'],
                           PParam['srate'],PParam['eps_c'],PParam['eps_s'],x0[0],
                           x0[1],x0[2],x0[3],1-x0[2]-x0[3],x0[4],
                           x0[5],MU3]
        #******************************************************************************
        # Calculation timings
        #******************************************************************************
        count += 1           # Increment calculation counter
    
    #******************************************************************************
    #******************************************************************************
    # Sort and output data
    #******************************************************************************
    #******************************************************************************
    # Write out full data in CSV format
    writeDataCSV(fulldata,os.path.join(DDATA,'RawData'),outfile,header)
    #******************************************************************************
    # Print the calculation timings
    #******************************************************************************
    tc11 = perf_counter() # Stop clock timer
    tc   = tc11-tc00      # Calculation duration from clock
    if GUI is not None:
        return
    elif Print is not None:
        print('')
        print('Final clock timer (s): ' + format_time(tc))
        print('All done!')
    return

def Pass_Sim(SParam, PParam, AdvParam, Print=False, GUI=None, outfile=OUTFILE):

    if None in [SParam, PParam, AdvParam]:
        if Print is not None:
            print("Error: Some parameters are not defined:")
            print('System Parameters:')
            for p in SParam.items():
                print(' '*2 + '> ' + p[0] + ': ' + str(p[1]),end="")
                if p[1] is None:
                    print(' (!)')
                else:
                    print('\n')
        return

    pass_data = loadData(os.path.join(DDATA,'Micius Data'), "Micius.txt", usecols=(0,1,2))
    time = pass_data[:,0]
    elevation = pass_data[:,1]

    states = SParam['states']
    decoys = SParam['decoys']
    P = sc.protocol(states,decoys,PParam)
    P.bound_method(AdvParam['boundM'])
    P.error_func(AdvParam['ECfunc'])

    h         = SParam['h']
    dtrans    = SParam['dt']
    dreceiv   = SParam['dr']
    w         = SParam['w']
    wl        = SParam['wl']

    # Losses
    dif  = diffraction(elevation, h, dtrans, dreceiv, w, wl)
    atm  = atmospheric(elevation, pass_data[:,2])
    intr = np.full(len(elevation), AdvParam['intr'])
    loss_dB = dif + atm + intr - 16
    eff = 10**( -loss_dB*0.1)

    #******************************************************************************
    # Initialise loop ranges
    #******************************************************************************
    Total_c = len(eff)

    #******************************************************************************
    # Initialise calculation parameters
    #******************************************************************************
    # Store initial/fixed protocol parameters as an array
    order = ['pax', 'pbx', 'pk1', 'pk2', 'mu1', 'mu2']
    x0 = np.array([PParam[i] for i in order])
    x0i = np.empty((len(order),Total_c))
    # Store optimisation parameters
    method = AdvParam['optM']
    optparam = opts.opt_method(method)
    for item in list(AdvParam['optP'].items()):
        optparam.edit(item[0],item[1])
    O = opts.bounds(optparam)
    for n in order:
        O.add_var(AdvParam['R_'+n],n)
    bounds,cons,options = O.opt_param(decoys)
    Nmax = O.method.param['Nmax']
    NoptMin = O.method.param['NoptMin']
    tStopBetter = O.method.param['tStopBetter']
    tStopZero = O.method.param['tStopZero']

    #******************************************************************************
    # Initialise output data storage and headers
    #******************************************************************************
    # Header for CSV file: Data columns
    header = "SysLoss [dB],DiffLoss,AtmLoss,IntrLoss,time,dt,elev,SKL,QBERx,phiX,nX,nZ,mX,lambdaEC,sX0,sX1,vZ1,sZ1," + \
             "mean photon no.,D1,D2,D3,e1,e2,e3,QBERI,Pec,Pap,Rrate,eps_c,eps_s," + \
             "PAx,PBx,P1,P2,P3,mu1,mu2,mu3"
    # Initialise a data storage array: shape(No. of data runs, No. of metrics)   
    fulldata = np.empty((Total_c,len(header.split(","))))
    if method == 'COBYLA':
        opt_head = "Nopt,Ntot,x0i,x1i,x2i,x3i,x4i,x5i,x0,x1,x2,x3,x4,x5,1/fun," + \
                   "status,success,nfev,maxcv"
    elif method == 'SLSQP':
        opt_head = "Nopt,Ntot,x0i,x1i,x2i,x3i,x4i,x5i,x0,x1,x2,x3,x4,x5,1/fun," + \
                   "status,success,nfev,njev,Nit"
    else:
        opt_head = ""
    # Initialise a data storage array for optimiser metrics
    optdata = np.empty((Total_c,len(opt_head.split(","))))

    if GUI is not None:
        return
    elif Print is True:
        print('_'*BAR_SIZE,'\n')
        print('-'*BAR_SIZE)
        print('QuSatProtocol v1.0')
        print('-'*BAR_SIZE)
        print('Efficient BB84 security protocol')
        print('For {} States and {} Decoys'.format(P.states,P.decoys))
        print('Using optimised protocol parameters: Px, pk, mu')
        print('Using {} bounds for statistical fluctuations'.format(P.boundM))
        if P.ECfunc in ['logM','logm']:
            print('Error correction model: logM(nX, QBERx, eps_c)')
        elif P.ECfunc in ['Block','block']:
            print('Error correction model: 1.16 * nX * h(QBERx)')
        elif P.ECfunc in ['mXtot','mxtot']:
            print('Error correction model: 1.16 * mXtot')
        elif P.ECfunc in ['None','none']:
            print('Error correction term not included')
        print('-'*BAR_SIZE)
        print('_'*BAR_SIZE,'\n')

    tc00 = perf_counter()       # Start clock timer
    tc = 0
    count = 0                   # Initialise calculation counter
    Err = []
    for i in tqdm(range(Total_c),ncols=BAR_SIZE,disable=not Print):

        if GUI is not None:
            return
        
        if i >= len(eff)-1 :
            delta = np.abs(time[i]-time[i-1])
        else:
            delta = np.abs(time[i+1]-time[i])

        if count > 0:
            x0 = x0i[:,count-1]
        
        # Use the optimised parameters from the previous shift angle calculation
        P.System_state(eff[i],delta)
        res = minimize(P.key_length_inv,x0,args=(),method=method,
                       jac=None,hess=None,hessp=None,bounds=bounds, 
                       constraints=cons,tol=None,callback=None, 
                       options=options)
        Ntot = res.nfev # Initilaise total No. of function evaluations
        Nopt = 1        # Number of optimisation calls
        # Re-run optimization until Nmax function evaluations
        # have been used. Take a copy of initial results to compare.
        x0c  = x0
        SKLc = int(1.0 / res.fun)
        resc = res
        Nzero = 0 # Number of times we get SKL == 0
        while Ntot < Nmax or Nopt < NoptMin:
            Nopt += 1
            # Randomly assign new initial parameters
            x0 = x0_rand(O.lowerb(),O.upperb(),ZERO,decoys)
            # Calculate optimised SKL
            res = minimize(P.key_length_inv,x0,args=(),method=method,
                            jac=None,hess=None,hessp=None,bounds=bounds, 
                            constraints=cons,tol=None,callback=None, 
                            options=options)
            if int(1.0 / res.fun) > 0:
                if int(1.0 / res.fun) > SKLc:
                    if Nopt >= NoptMin and tStopBetter:
                        break # A better value was found!
                    else:
                        # Store new set of best parameters
                        x0c  = x0
                        resc = res
                        SKLc = int(1.0 / res.fun)
                else:
                    # Reset to the best parameters
                    x0  = x0c
                    res = resc
            else:
                # SKL = 0. Reset to the 'best' parameters,
                # (may still give SKL = 0).
                Nzero += 1
                if Nopt > NoptMin:
                    if Nzero / (Nopt - 1) == 1:
                        # We get SKL = 0 every time.
                        if tStopZero:
                            break
                x0  = x0c
                res = resc
            Ntot += res.nfev

        if not res.success:
            aux = "Optimiser status = {}: {}".format(res.status, res.message)
            if aux in [entry[0] for entry in Err]:
                for i in range(len(Err)):
                    if Err[i][0] == aux:
                        Err[i][1] += 1
            else:
                Err += [["Optimiser status = {}: {}".format(res.status, res.message),1]]
        pp = {}
        for i in range(len(order)):
            pp[order[i]] = res.x[i]
        sc.paramcheck(decoys,pp)
        # Get final parameters from standard key length function
        SKL,QBERx,phi_x,nX,nZ,mX,lambdaEC,sX0,sX1,vz1,sZ1,mpn,Dj,ej = P.key_length(res.x)
        # Store calculation parameters
        fulldata[count] = [loss_dB[count],dif[count],atm[count],intr[count],
                           time[count],delta,elevation[count],SKL,QBERx,phi_x,nX,nZ,mX,lambdaEC,sX0,sX1,
                           vz1,sZ1,mpn,Dj[0],Dj[1],dim_check(Dj,2),ej[0],ej[1],dim_check(ej,2),
                           PParam['QBERI'],PParam['pec'],PParam['pap'],
                           PParam['srate'],PParam['eps_c'],PParam['eps_s'],res.x[0],
                           res.x[1],res.x[2],res.x[3],1-res.x[2]-res.x[3],res.x[4],
                           res.x[5],MU3]
        # Store optimiser metrics
        optdata[count] = opts.getOptData(Nopt,Ntot,x0,res,method)
        # Store protocol parameters to initialise calculations
        if np.isnan(int(1.0 / res.fun)) or np.isinf(int(1.0 / res.fun)):
            x0i[:,count] = x0_rand(O.lowerb(),O.upperb(),ZERO,decoys)
        else:
            if int(1.0 / res.fun) > 0:
                x0i[:,count] = res.x
            else:
                x0i[:,count] = x0_rand(O.lowerb(),O.upperb(),ZERO,decoys)
        #******************************************************************************
        # Calculation timings
        #******************************************************************************
        count += 1           # Increment calculation counter
    
    #******************************************************************************
    #******************************************************************************
    # Sort and output data
    #******************************************************************************
    #******************************************************************************
    # Write out full data in CSV format
    writeDataCSV(fulldata,os.path.join(DDATA,'RawData'),outfile,header)
    # Write optimiser metrics
    #print(res.keys()) # Prints available outputs for the object res
    writeDataCSV(optdata,os.path.join(DDATA,'RawData'),outfile+'_metrics',opt_head)
    #******************************************************************************
    # Print the calculation timings
    #******************************************************************************
    tc11 = perf_counter() # Stop clock timer
    tc   = tc11-tc00      # Calculation duration from clock
    if GUI is not None:
        return
    elif Print is True:
        print('')
        print('Error Messages:')
        for entry in Err:
            print('({}) > '.format(entry[1]) + entry[0])
        print('')
        print('Final clock timer (s): ' + format_time(tc))
        print('All done!')
    return

def diffraction(elevation, altitude, dtransmitter, dreceiver, beamwaist, wavelenght):
    d = 0.5 * ( 2*RT*np.cos(elevation + np.pi/2) + np.sqrt( 4*altitude**2 + 8*altitude*RT + (2*RT*np.cos(elevation + np.pi/2))**2 ) )
    pc = 2.44*wavelenght/dtransmitter
    div = wavelenght/(np.pi*beamwaist)
    m = 1-np.abs(pc-div)/pc

    return -20*np.log10( dreceiver*dtransmitter / ( dtransmitter**2 + 2.44*m*d*wavelenght ) )

def atmospheric(elev,eff):
    dif = diffraction(elev,M_H,M_DT,M_DR,M_W,M_WL)
        
    return -10*np.log10(eff) - dif - M_INT

def x0_rand(lb,ub,num_min,decoys):
    """
    Randomly initialise the 7 protocol parameters using the specified bounds.
    Parameters and bounds should be specified in the order {Px,pk1,pk2,mu1,mu2}.

    Parameters
    ----------
    xb : float, array-like
        Upper and lower bounds for the protocol parameters. (5,2)
    num_min : float
        An arbitrarily small number.

    Returns
    -------
    x0 : float, array
        Randomly initialised protocol parameters.

    """
    # Pax_i  = np.random.rand() * (ub[0] - lb[0] - 2*num_min) + lb[0] + num_min
    # Pbx_i  = np.random.rand() * (ub[1] - lb[1] - 2*num_min) + lb[1] + num_min
    Pax_i  = np.random.uniform(lb[0],ub[0])
    Pbx_i  = np.random.uniform(lb[1],ub[1])
    if decoys == 1:
        # pk1_i = np.random.rand() * (ub[2] - lb[2] - 2*num_min) + lb[2] + num_min
        pk1_i = np.random.uniform(lb[2],ub[2])
        pk2_i = 1-pk1_i
    elif decoys == 2:
        pk1_i, pk2_i = 1.0, 1.0
        while (pk1_i+pk2_i >= 1.0):
            pk1_i = np.random.rand() * (ub[2] - lb[2] - 2*num_min) + lb[2] + num_min
            pk2_i = np.random.rand() * (min(ub[3],1-pk1_i) - lb[3] - 2*num_min) + lb[3] + num_min
    #mu3_i = np.random.rand() * (ub[6] - lb[6] - 2*num_min) + lb[6] + num_min
    mu3_i = MU3
    mu1_i = np.random.rand() * (ub[4] - max(lb[4],2*mu3_i) - 2*num_min) + max(lb[4],2*mu3_i) + num_min
    mu2_i = np.random.rand() * (min(ub[5],mu1_i) - max(lb[5],mu3_i) - 2*num_min) + max(lb[5],mu3_i) + num_min
    return np.array([Pax_i,Pbx_i,pk1_i,pk2_i,mu1_i,mu2_i])

def format_time(s):
    ms = int(s*1000)
    min, ms = divmod(ms,60000)
    sec, ms = divmod(ms,1000)
    return "{:02d}:{:02d}:{:03d}".format(min,sec,ms)

def ratio_check(value):
    if isinstance(value,float):
        if (value > 1.0 or value < 0.0):
            return False
    else:
        return False
    return True

def dim_check(arr,index):
    if index > len(arr)-1:
        return -1
    else:
        return arr[index]
