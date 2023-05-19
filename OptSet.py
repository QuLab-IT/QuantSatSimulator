import numpy as np
from sys import float_info
from scipy.optimize import Bounds

# Numerical value to use when denominator values are potentially zero
# Extract the smallest float that the current system can:
ZERO = float_info.epsilon # round (relative error due to rounding)
# Extract the largest float that the current system can represent
INF = float_info.max

MU3 = 0.0

F_or_T = [False, True]       # List used to switch between False or True values
OPT    = ['COBYLA','SLSQP']  # List of optimisation methods

class opt_method():

    def __init__(self,method):
        
        if method in OPT:
            self.name = method
            self.param = { 'NoptMin': 10,            # Minimum No. of optimisations to strive for
                           'NoptMax': 10000,         # Maximum No. of optimisations (not used)
                           'tStopZero': F_or_T[1],   # Stop optimizing if the first NoptMin return SKL = 0?
                           'tStopBetter': F_or_T[1], # Stop after NoptMin optimizations if SKL improved?
                           }

            if method == 'COBYLA':
                # Optimiser options: minimize, method='COBYLA'
                self.param['Nmax']   = 1000       # Max No. of iterations
                self.param['ctol']   = 1.0e-12    # Constraint absolute tolerance.
                self.param['rhobeg'] = 0.002      # Reasonable initial changes to the variables.
            else:
                # Optimiser options: minimize, method='SLSQP'
                self.param['Nmax'] = 1000         # Max No. of iterations
                self.param['ftol'] = 1.0e-12      # Precision goal for the value of f in the stopping criterion.
                self.param['eps']  = 1.0e-7       # Step size used for numerical approximation of the Jacobian.
        else:
            print('Error: {} does not exist in the Optimisation Methods list:'.format(method))
            for m in OPT:
                print(' '*3 + '> ' + m)
            return None    

    def param_list(self, print_list=None):
        if print_list is not None:
            print("Optimisation parameter list:")
            for key,value in self.param.items():
                print(key + ': ' + str(value))
        return self.param

    def edit(self,key,value):
        if key in self.param:
            self.param[key] = value
        else:
            print('Warning! {} does not exist in the parameter list:'.format(key))
            self.param_list(print_list = F_or_T[1])
        return

class bounds():

    def __init__(self,method):
        self.method = method
        self.name = {}
        self.fixed = {}
        self.order = ['pax', 'pbx', 'pk1', 'pk2', 'mu1', 'mu2']

    def add_var(self,bounds,name):
        self.name[name] = bounds
        if bounds[0] == bounds[1]:
            self.fixed[name] = 1
        else:
            self.fixed[name] = 0

    def pop_var(self,name):
        if name in self.name.keys():
            self.name.popitem(name)
        else:
            return

    def lowerb(self):
        lb = []
        for n in self.order:
            if self.name[n][0] == self.name[n][1]: #Fixed
                if self.name[n][0] == 0:
                    lb += [ZERO]
                elif self.name[n][0] == 1:
                    lb += [1 - ZERO]
                else:
                    lb += [self.name[n][0] - digit_order(self.name[n][0])/100]
            else:
                #lb += [self.name[n][0] + ZERO]
                lb += [self.name[n][0]]

        return lb

    def upperb(self):
        ub = []
        for n in self.order:
            if self.name[n][0] == self.name[n][1]: #Fixed
                if self.name[n][0] == 0:
                    ub += [2*ZERO]
                elif self.name[n][0] == 1:
                    ub += [1]
                else:
                    ub += [self.name[n][0] + digit_order(self.name[n][0])/100]
            else:
                #ub += [self.name[n][1] - ZERO]
                ub += [self.name[n][1]]

        return ub

    def fix_rep(self):
        f = []
        for n in self.order:
            f += [self.fixed[n]]
        return f

    def opt_param(self,decoys):
        lb = self.lowerb()
        ub = self.upperb()
        # Store the parameters bounds as an object, as required by minimize() for
        cons_type = 'ineq' # C_j[x] >= 0
        if decoys == 1:
            lin_cons = {'type' : cons_type,
                            'fun'  : cons_1D,
                            'jac'  : cons_1D_jac}
        elif decoys == 2:
            lin_cons = {'type' : cons_type,
                            'fun'  : cons_2D,
                            'jac'  : cons_2D_jac}
        # method='COBYLA'
        if self.method.name == 'COBYLA':
            b = [(lb[i],ub[i]) for i in range(len(lb))]
            # Set specific optimiser options
            options = {'rhobeg': self.method.param['rhobeg'],
                       'maxiter': self.method.param['Nmax'], 
                       'disp': False, 
                       'catol': self.method.param['ctol']}
        elif self.method.name == 'SLSQP':
            b = Bounds(lb, ub)
            # Set specific optimiser options
            options = {'maxiter': self.method.param['Nmax'],
                       'ftol': self.method.param['ftol'],
                       'iprint': 1,
                       'disp': False,
                       'eps': self.method.param['eps'],
                       'finite_diff_rel_step': None}
        return b, lin_cons, options

            
def cons_1D(x):
    # Linear constraint inequality function:
    #   (1)         1 - pk1 >= 0,
    #   (2) mu1 - mu2 - eps >= 0,
    #   (3)       mu2 - eps >= 0,
    # where eps is an arbitrarily small number
    return np.array([1 - x[2],
                     x[4] - x[5] - ZERO,
                     x[5] - ZERO])

def cons_1D_jac(x):
    # Jacobian of the linear constraint inequality function
    return np.array([[0,0,-1,0,0,0],
                     [0,0,0,0,1,-1],
                     [0,0,0,0,0,1]])

def cons_2D(x):
    # Linear constraint inequality function:
    #   (1)         1 - pk1 - pk2 >= 0,
    #   (2) mu1 - mu2 - mu3 - eps >= 0,
    #   (3)       mu2 - mu3 - eps >= 0,
    # where eps is an arbitrarily small number
    return np.array([1 - x[2] - x[3],
                     x[4] - x[5] - MU3 - ZERO,
                     x[5] - MU3 - ZERO])

def cons_2D_jac(x):
    # Jacobian of the linear constraint inequality function
    return np.array([[0,0,-1,-1,0,0],
                     [0,0,0,0,1,-1],
                     [0,0,0,0,0,1]])

def getOptData(Nopt,Ntot,x0,res,method):
    """
    Returns a list of output metrics from the scipy.optimize results object res.

    Parameters
    ----------
    Nopt : integer
        Number of optimisations performed.
    Ntot : integer
        Total number of function evaluations.
    x0 : float, array-like
        Initial protocol parameters
    res : object, dictionary
        Optimisation results.
    method : string
        Optimization method.

    Returns
    -------
    optData : list
        List of optimisation metrics and data.

    """
    if method == 'COBYLA':
        # 'Nopt' 'Ntot' 'x0', 'x', 'fun', 'status', 'success', 'nfev', 'maxcv'
        return [Nopt,Ntot,*x0,*res.x,res.fun,res.status,res.success,res.nfev,
                res.maxcv]
    elif method == 'SLSQP':
        # 'Nopt' 'Ntot' 'x0', 'x', 'fun', 'status', 'success', 'nfev', 'njev',
        # 'nit'
        return [Nopt,Ntot,*x0,*res.x,res.fun,res.status,res.success,res.nfev,
                res.njev,res.nit]
    else:
        return []

def digit_order(num):
    str_num = '{:.15g}'.format(num)
    return 10**-(len(str_num)-str_num.index('.')-1)