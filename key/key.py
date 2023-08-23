from sys import float_info

import numpy as np
from scipy.optimize import minimize
from tqdm import tqdm

from channel.time_dependent_loss import get_losses
#from key.protocols.fourstatetwodecoy.init_4s2d import arrange_output
#from key.protocols.fourstatetwodecoy.init_4s2d import (x0_rand, check_constraints)
#from key.protocols.fourstatetwodecoy.key_4s2d import (set_params, key_length, key_length_inv, key_length_sim)
#from key.protocols.threestateonedecoy.init_3s1d import (x0_rand, check_constraints)
#rom key.protocols.threestateonedecoy.key_3s1d import (set_params, key_length, key_length_inv, key_length_sim)
from output.outputs import (getOptData, writeDataCSV, write_data, writeMultiData,
                            get_timings, format_time)

__all__ = ['optimise_loop_sim','optimise_loop_inv','check_opt','store_data','optimise_key',
           'stack_arrays','cut_arrays','SKL_opt_loop_loss_and_time',
           'SKL_loop_loss_and_time','SKL_sys_loop','sys_param_list',
           'args_fixed_list','SKL_main_loop']

###############################################################################

x0_random           = None
set_parameters      = None
key_length_func     = None
key_length_inverse  = None
key_length_simetric = None
arrange_out         = None

###############################################################################

num_min = float_info.epsilon # A very small number

BAR_SIZE = 80

list_protocols = ['aBB84-WCP','3S1D']

###############################################################################


def optimise_loop_sim(tInit,x,mu3,xb,args,bounds,cons,options,opt_params):
    """
    Execute a loop of optimisations of the main protocol parameters limited by
    either the total number of optimisations or the aggregate number of function
    evaluations. The initial parameters for optimisation are taken from
    previous successful (non-zero) optimisations, otherwise they are generated
    randomly.

    Parameters
    ----------
    tInit : bool
        Has an initial set of parameters been specified.
    x : float, array-like
        The initial set of parameters to use.
    x0i : float, array-like
        The final optimised parameters from previous calculations.
    ci : int, array-like
        Calculation loop counter array.
    mu3 : float
        Fixed protocol parameter. The intensity of the third pulse.
    xb : float, array-like
        The upper and lower bounds for the optimised parameters.
    args : dict
        Dictionary of extra optimiser arguments.
    bounds : obj
        Scipy bounds object containing optimised parameter bounds.
    cons : obj
        Scipy object containing the optimisation constraints.
    options : dict
        Dictionary of optimiser parameters.
    opt_params : dict
        Dictionary of parameters related to optimisation.

    Returns
    -------
    res : dict
        Dictionary of optimisation results.
    x0 : float, array-like
        Initial parameters used for final optimisation.
    Nopt : int
        Number of optimisations performed.
    Ntot : int
        Number of function evaluations.

    """
    # Re-set initial parameters (if required)
    Ninit = 0
    if tInit:
        # Store initial/fixed protocol parameters as an array
        x0 = x
    else:
        x0 = x0_random(mu3,xb,num_min)
        
    # Initial optimisation of SKL
    res = minimize(key_length_simetric,x0,args=(args,),method=opt_params['method'],
                   jac=None,hess=None,hessp=None,bounds=bounds, 
                   constraints=cons,tol=None,callback=None, 
                   options=options)
    Ntot = res.nfev # Initilaise total No. of function evaluations
    Nopt = 1        # Number of optimisation calls
    # Re-run optimization until Nmax function evaluations
    # have been used. Take a copy of initial results to compare.
    x0_   = x0
    SKL_  = int(-res.fun)
    res_  = res
    Nzero = 0 # Number of times we get SKL == 0
    while Nopt < opt_params['NoptMin'] or Ntot < opt_params['Nmax']:
        # Initialise the optimised parameters
        x0 = x0_random(mu3,xb,num_min)
        # Calculate optimised SKL
        res = minimize(key_length_simetric,x0,args=(args,),method=opt_params['method'],
                       jac=None,hess=None,hessp=None,bounds=bounds, 
                       constraints=cons,tol=None,callback=None, 
                       options=options)
        Nopt += 1 # Increment optimisation counter
        if int(-res.fun) > 0:
            if int(-res.fun) > SKL_:
                if Nopt >= opt_params['NoptMin'] and opt_params['tStopBetter']:
                    break # A better value was found!
                else:
                    # Store new set of best parameters
                    x0_  = x0
                    res_ = res
                    SKL_ = int(-res.fun)
            else:
                # Reset to the best parameters
                x0  = x0_
                res = res_
        else:
            # SKL = 0. Reset to the 'best' parameters,
            # (may still give SKL = 0).
            Nzero += 1
            if Nopt > opt_params['NoptMin']:
                if Nzero / (Nopt - 1) == 1:
                    # We get SKL = 0 every time.
                    if opt_params['tStopZero']:
                        break
                    x0  = x0_
                    res = res_
        Ntot += res.nfev

    return res, x0, Nopt, Ntot

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def optimise_loop_inv(tInit,x,mu3,xb,args,bounds,cons,options,opt_params):
    """
    Execute a loop of optimisations of the main protocol parameters limited by
    either the total number of optimisations or the aggregate number of function
    evaluations. The initial parameters for optimisation are taken from
    previous successful (non-zero) optimisations, otherwise they are generated
    randomly.

    Parameters
    ----------
    tInit : bool
        Has an initial set of parameters been specified.
    x : float, array-like
        The initial set of parameters to use.
    x0i : float, array-like
        The final optimised parameters from previous calculations.
    ci : int, array-like
        Calculation loop counter array.
    mu3 : float
        Fixed protocol parameter. The intensity of the third pulse.
    xb : float, array-like
        The upper and lower bounds for the optimised parameters.
    args : dict
        Dictionary of extra optimiser arguments.
    bounds : obj
        Scipy bounds object containing optimised parameter bounds.
    cons : obj
        Scipy object containing the optimisation constraints.
    options : dict
        Dictionary of optimiser parameters.
    opt_params : dict
        Dictionary of parameters related to optimisation.

    Returns
    -------
    res : dict
        Dictionary of optimisation results.
    x0 : float, array-like
        Initial parameters used for final optimisation.
    Nopt : int
        Number of optimisations performed.
    Ntot : int
        Number of function evaluations.

    """
    # Re-set initial parameters (if required)
    Ninit = 0
    if tInit:
        # Store initial/fixed protocol parameters as an array
        x0 = x
    else:
        x0 = x0_random(mu3,xb,num_min)
        
    # Initial optimisation of SKL
    res = minimize(key_length_inverse,x0,args=(args,),method=opt_params['method'],
                   jac=None,hess=None,hessp=None,bounds=bounds, 
                   constraints=cons,tol=None,callback=None, 
                   options=options)
    Ntot = res.nfev # Initilaise total No. of function evaluations
    Nopt = 1        # Number of optimisation calls
    # Re-run optimization until Nmax function evaluations
    # have been used. Take a copy of initial results to compare.
    x0_   = x0
    SKL_  = int(1.0 / res.fun)
    res_  = res
    Nzero = 0 # Number of times we get SKL == 0
    while Nopt < opt_params['NoptMin'] or Ntot < opt_params['Nmax']:
        # Initialise the optimised parameters
        x0 = x0_random(mu3,xb,num_min)
        # Calculate optimised SKL
        res = minimize(key_length_inverse,x0,args=(args,),method=opt_params['method'],
                       jac=None,hess=None,hessp=None,bounds=bounds, 
                       constraints=cons,tol=None,callback=None, 
                       options=options)
        Nopt += 1 # Increment optimisation counter
        if int(1.0 / res.fun) > 0:
            if int(1.0 / res.fun) > SKL_:
                if Nopt >= opt_params['NoptMin'] and opt_params['tStopBetter']:
                    break # A better value was found!
                else:
                    # Store new set of best parameters
                    x0_  = x0
                    res_ = res
                    SKL_ = int(1.0 / res.fun)
            else:
                # Reset to the best parameters
                x0  = x0_
                res = res_
        else:
            # SKL = 0. Reset to the 'best' parameters,
            # (may still give SKL = 0).
            Nzero += 1
            if Nopt > opt_params['NoptMin']:
                if Nzero / (Nopt - 1) == 1:
                    # We get SKL = 0 every time.
                    if opt_params['tStopZero']:
                        break
                    x0  = x0_
                    res = res_
        Ntot += res.nfev

    return res, x0, Nopt, Ntot

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def check_opt(res,mu3):
    """
    Check optimiser output is within bounds and constraints of protocol.

    Parameters
    ----------
    res : dict
        Dictionary of optimisation results.
    mu3 : float
        Intensity of pulse three (vacuum).
    method : str
        Optimisation method.
    tPrint : bool
        Print output to std out?

    Returns
    -------
    None.

    """
    Err = None
    if not res.success:
        Err = [res.status, res.message]
    # Check if optimised parameters satisfy the constraints
    #check_constraints(res.x[0],res.x[1],res.x[2],res.x[3],res.x[4],mu3)
    return Err

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def store_data_sim(res,Nopt,Ntot,x0,method,mu3,xb,args):
    """
    Store results, calculation parameters and optimisation metrics as arrays.

    Parameters
    ----------
    res : dict
        Dictionary of optimisation results.
    Nopt : int
        Number of optimisations performed.
    Ntot : int
        Number of function evaluations.
    x0 : float, array-like
        Initial parameters used for final optimisation.
    method : str
        Optimisation method.
    mu3 : float
        Intensity of pulse three (vacuum).
    xb : float, array-like
        The upper and lower bounds for the optimised parameters.
    args : dict
        Dictionary of extra optimiser arguments.

    Returns
    -------
    fulldata : float, array-like
        Calculation output parameters.
    optdata : list
        Optimisation metrics and parameters.
    x0i : float, array-like
        The final optimised protocol parameters including this calculation.

    """
    # Get final parameters from standard key length function
    SKL,QBERx,phi_x,nX,nZ,mX,lambdaEC,sX0,sX1,vz1,sZ1,mpn = key_length_func(res.x,args)
    # Store calculation parameters
    fulldata = [SKL,QBERx,phi_x,nX,nZ,mX,lambdaEC,sX0,sX1,vz1,sZ1,mpn]
    # Store optimiser metrics
    optdata = getOptData(Nopt,Ntot,x0,res,method)
    # Store protocol parameters to initialise calculations
    if np.isnan(int(-res.fun)) or np.isinf(int(-res.fun)):
        # SKL = 0 or was invalid. Store a random set of optimised params for 
        # intialising future calculations
        x0i = x0_random(mu3,xb,num_min)
    else:
        if int(-res.fun) > 0:
            # SKL > 0. Store optimised parameters to intialise future 
            # calculations
            x0i = res.x
        else:
            # Unexpected SKL result. Store a random set of optimised params for 
            # intialising future calculations
            x0i = x0_random(mu3,xb,num_min)
    return fulldata, optdata, x0i

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def store_data_inv(res,Nopt,Ntot,x0,method,mu3,xb,args):
    """
    Store results, calculation parameters and optimisation metrics as arrays.

    Parameters
    ----------
    res : dict
        Dictionary of optimisation results.
    Nopt : int
        Number of optimisations performed.
    Ntot : int
        Number of function evaluations.
    x0 : float, array-like
        Initial parameters used for final optimisation.
    method : str
        Optimisation method.
    mu3 : float
        Intensity of pulse three (vacuum).
    xb : float, array-like
        The upper and lower bounds for the optimised parameters.
    args : dict
        Dictionary of extra optimiser arguments.

    Returns
    -------
    fulldata : float, array-like
        Calculation output parameters.
    optdata : list
        Optimisation metrics and parameters.
    x0i : float, array-like
        The final optimised protocol parameters including this calculation.

    """
    # Get final parameters from standard key length function
    SKL,QBERx,phi_x,nX,nZ,mX,lambdaEC,sX0,sX1,vz1,sZ1,mpn = key_length_func(res.x,args)
    # Store calculation parameters
    fulldata = [SKL,QBERx,phi_x,nX,nZ,mX,lambdaEC,sX0,sX1,vz1,sZ1,mpn]
    # Store optimiser metrics
    optdata = getOptData(Nopt,Ntot,x0,res,method)
    # Store protocol parameters to initialise calculations
    if np.isnan(int(1.0 / res.fun)) or np.isinf(int(1.0 / res.fun)):
        # SKL = 0 or was invalid. Store a random set of optimised params for 
        # intialising future calculations
        x0i = x0_random(mu3,xb,num_min)
    else:
        if int(1.0 / res.fun) > 0:
            # SKL > 0. Store optimised parameters to intialise future 
            # calculations
            x0i = res.x
        else:
            # Unexpected SKL result. Store a random set of optimised params for 
            # intialising future calculations
            x0i = x0_random(mu3,xb,num_min)
    return fulldata, optdata, x0i

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def optimise_key(tInit,x,mu3,xb,args,bounds,cons,options,opt_params):
    """
    Find the optimised parameters, by optimising repeteadly, then check the
    parameters and then store the output data in arrays.

    Parameters
    ----------
    tInit : bool
        Have initial values for the protocol parameters been specified?
    x : float, array-like
        The initial set of parameters to use.
    x0i : float, array-like
        The final optimised protocol parameters from previous calculations.
    ci : int, array-like
        Calculation loop counter array.
    mu3 : float
        Fixed protocol parameter. The intensity of the third pulse
    xb : float, array-like
        The upper and lower bounds for the optimised parameters.
    args : dict
        Dictionary of extra optimiser arguments.
    bounds : obj
        Scipy bounds object containing optimised parameter bounds.
    cons : obj
        Scipy object containing the optimisation constraints.
    options : dict
        Dictionary of optimiser parameters.
    opt_params : dict
        Dictionary of parameters related to optimisation.
    tPrint : bool
        Print output to std out?

    Returns
    -------
    res : dict
        Dictionary of optimisation results.
    fulldata : float, array-like
        Calculation output parameters.
    optdata : list
        Optimisation metrics and parameters.
    x0i : float, array-like
        The final optimised protocol parameters including this calculation.

    """
    # Set Objective funtion
    if opt_params['obj_fun'] == 'Symetric':
        optimise_loop =  optimise_loop_sim
        store_data    =  store_data_sim
    elif opt_params['obj_fun'] == 'Inverse':
        optimise_loop =  optimise_loop_inv
        store_data    =  store_data_inv
    else:
        optimise_loop =  optimise_loop_sim
        store_data    =  store_data_sim
    # Run optimisation loop
    res, x0, Nopt, Ntot = optimise_loop(tInit,x,mu3,xb,args,bounds,cons,
                                        options,opt_params)
    # Check optimised output
    Err = check_opt(res,mu3)
    # Retrive data, metrics, and initial values
    fulldata, optdata, x0i = store_data(res,Nopt,Ntot,x0,opt_params['method'],
                                        mu3,xb,args)
    return res, fulldata, optdata, Err

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def stack_arrays(arr):
    rev = arr[::-1]
    stack = np.vstack((arr,rev[1:-1]))
    return stack

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def cut_arrays(Array1,Array2):
    zipp = zip(Array1,Array2)
    Arr = list(sorted(zipp, key=lambda x: x[0]))
    inds = []
    for i in range(1,len(Arr)):
        if Arr[i][0] == Arr[i-1][0]:
            inds += [i]
    Filter = [value for index, value in enumerate(Arr) if index not in inds]
    Ar1, Ar2 = zip(*Filter)

    return np.array(Ar1) , np.array(Ar2)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def SKL_opt_loop_loss_and_time(theta_max,Pec,QBERI,dt,tInit,x,mu3,xb,args_fixed,bounds,
                               cons,options,opt_params,fulldata,
                               optdata,sys_params,Syseff,Sysloss):
    """
    Perfom secret key optimisations for iterated values of the transmission
    window half-width (dt) and the excess system loss (ls).

    Parameters
    ----------
    count : int
        Absolute calculation counter.
    ci : int, array-like
        Calculation loop counter.
    ni : int, array-like
        The number of iterations in each loop.
    theta_max : float
        Maximum elevation of satellite overpass (rad).
    Pec : float
        Extraneous count rate/probability.
    QBERI : float
        Intrinsic quantum bit error rate.
    ls_range : int, array-like (3)
        Start, stop, and No. of step values for the excess loss loop.
    dt_range : int, array-like (3)
        Start, stop, and step values for the transmission window (half-width) 
        loop (s).
    tInit : bool
        Have initial values for the protocol parameters been specified?
    x : float, array-like
        The initial set of parameters to use.
    x0i : float, array-like
        The initial parameters to use.
    ci : int, array-like
        Calculation loop counter array.
    mu3 : float
        Fixed protocol parameter. The intensity of the third pulse.
    xb : float, array-like
        The upper and lower bounds for the optimised parameters.
    args_fixed : dict
        Dictionary of arguments required by key_length functions.
    bounds : obj
        Scipy bounds object containing optimised parameter bounds.
    cons : obj
        Scipy object containing the optimisation constraints.
    options : dict
        Dictionary of optimiser parameters.
    opt_params : dict
        Dictionary of parameters related to optimisation.
    tPrint : bool
        Print output to std out?
    fulldata : float, array-like
        Calculation output parameters.
    optdata : list
        Optimisation metrics and parameters.
    sys_params : dict
        Additional system parameters to be included in fulldata.
    sysLoss : float
        The nominal loss or system loss metric (dB). For symmetric transmission
        windows this is the system loss at zenith.

    Returns
    -------
    fulldata : float, array-like
        Calculation output parameters.
    optdata : list
        Optimisation metrics and parameters.
    x0i : float, array-like
        The final optimised protocol parameters including this calculation.
    count : int
        Updated absolute calculation counter.

    """
    sr_eff, sr_delta = cut_arrays(Syseff,dt)

    Err = []
    xi = x
    # Perform the optimization of the secure key length
    for i in tqdm(range(len(sr_eff)),ncols=BAR_SIZE,position=0):
        # Store key params in a dict
        args = set_parameters(Pec,QBERI,sr_eff[i],sr_delta[i],*args_fixed)
        # Optimise key length
        res, SKLdata, optimise, aux = optimise_key(tInit,xi,mu3,xb,args,bounds,cons,options,
                            opt_params)
        # Get 0
        if SKLdata[0] == 0 or aux is not None:
            lim = opt_params['loop_range']
            for j in tqdm(range(lim),ncols=BAR_SIZE,position=1,leave=False):
                # Optimise key length
                res, SKLdata, optimise, aux = optimise_key(tInit,xi,mu3,xb,args,bounds,cons,options,
                                                            opt_params)
                if SKLdata[0] > 0 and aux is not None:
                    j = lim
        # Check Error Mensages
        if aux is not None:
            aux = "Optimiser status = {}: {}".format(aux[0], aux[1])
            if aux in [entry[0] for entry in Err]:
                for w in range(len(Err)):
                    if Err[w][0] == aux:
                        Err[w][1] += 1
            else:
                Err += [[aux,1]] 
        # Start with previous initial parameters
        xi = res.x
        # Store output data
        optdata  = np.vstack((optdata,optimise))
        fulldata = np.vstack((fulldata,arrange_out(-10*(np.math.log10(sr_eff[i])),0,dt[i],
                                                          mu3,QBERI,Pec,
                                                          args['Pbx'],res.x,
                                                          SKLdata,
                                                          sys_params)))
    fulldata = stack_arrays(fulldata)
    return fulldata, optdata, Err

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def SKL_loop_loss_and_time(theta_max,Pec,QBERI,x0,dt,mu3,args_fixed,
                            fulldata,sys_params,Syseff,Sysloss):
    """
    Perfom secret key (non-optimised) calculations for iterated values of the 
    transmission window half-width (dt) and the excess system loss (ls).

    Parameters
    ----------
    count : int
        Absolute calculation counter.
    ci : int, array-like
        Calculation loop counter.
    ni : int, array-like
        The number of iterations in each loop.
    theta_max : float
        Maximum elevation of satellite overpass (rad).
    Pec : float
        Extraneous count rate/probability.
    QBERI : float
        Intrinsic quantum bit error rate.
    ls_range : int, array-like (3)
        Start, stop, and No. of step values for the excess loss loop.
    dt_range : int, array-like (3)
        Start, stop, and step values for the transmission window (half-width) 
        loop (s).
    x0 : float, array-like
        The set of protocol parameters to use.
    mu3 : float
        Fixed protocol parameter. The intensity of the third pulse.
    args_fixed : dict
        Dictionary of arguments required by key_length functions.
    tPrint : bool
        Print output to std out?
    fulldata : float, array-like
        Calculation output parameters.
    sys_params : dict
        Additional system parameters to be included in fulldata.
    sysLoss : float
        The nominal loss or system loss metric (dB). For symmetric transmission
        windows this is the system loss at zenith.

    Returns
    -------
    fulldata : float, array-like
        Calculation output parameters.
    count : int
        Updated absolute calculation counter.

    """
    
    # Perform the optimization of the secure key length
    for i in range(len(Syseff)):
            # Store key params in a dict
            args = set_parameters(Pec,QBERI,Syseff[i],dt[i],*args_fixed)
            # SKL calculation
            SKL,QBERx,phi_x,nX,nZ,mX,lambdaEC,sX0,sX1,vz1,sZ1,mpn = key_length_func(x0,args)
            # Make a list of output data
            SKLdata = [SKL,QBERx,phi_x,nX,nZ,mX,lambdaEC,sX0,sX1,vz1,sZ1,mpn]
            # Store output data
            fulldata = np.vstack((fulldata,arrange_out(Sysloss[i],0,dt[i],
                                                          mu3,QBERI,Pec,
                                                          theta_max,x0,
                                                          SKLdata,
                                                          sys_params)))
            
    return fulldata

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def SKL_sys_loop(x,xb,theta_max,main_params,opt_params,
                 args_fixed,bounds,cons,options,header,opt_head,fulldata,optdata,
                 multidata,sys_params,Syseff,sysLoss,dt):
    """
    Calculate the SKL over the main iterated parameter loops. 

    Parameters
    ----------
    count : int
        Absolute calculation counter.
    ci : int, array-like
        Calculation loop counter.
    ni : int, array-like
        The number of iterations in each loop.
    x : float, array-like
        The initial set of protocol parameters to use.
    x0i : float, array-like
        The final protocol parameters from previous calculations.
    xb : float, array-like
        The upper and lower bounds for the optimised parameters.
    theta_max : float
        Maximum elevation of satellite overpass (rad).
    dt_range : int, array-like (3)
        Start, stop, and step values for the transmission window (half-width) 
    main_params : dict
        Dictionary of main calculation parameters.
    opt_params : dict
        Dictionary of parameters related to optimisation.
    args_fixed : dict
        Dictionary of arguments required by key_length functions.
    bounds : obj
        Scipy bounds object containing optimised parameter bounds.
    cons : obj
        Scipy object containing the optimisation constraints.
    options : dict
        Dictionary of optimiser parameters.
    header : str
        Comma separated string of output parameters corresponding to columns
        in fulldata. Used as header when writing output to file.
    opt_head : str
        Comma separated string of optimiser metrics and parameters corresponding
        to columns of optdata. Used as a header when writing output to file.
    fulldata : float, array-like
        Calculation output parameters.
    optdata : list
        Optimisation metrics and parameters.
    multidata : float, array-like
        Optimal dt output parameters for all calculations.
    sys_params : dict
        Additional system parameters to be included in fulldata.
    sysLoss : float
        The nominal loss or system loss metric (dB). For symmetric transmission
        windows this is the system loss at zenith.

    Returns
    -------
    multidata : float, array-like
        Updated optimal dt output parameters for all calculations.
    count : int
        Updated absolute calculation counter.

    """
    # Extract parameters for the calculation
    tOptimise = main_params['opt']['tOptimise'] # Perform an optimised calculation?
    tInit     = main_params['opt']['tInit']     # Use specified intial protocol paraneters for optimised calculation?
    mu3       = main_params['fixed']['mu3']     # The intensity of pulse 3 (BB84 decoy states).
    Pec       = main_params['iter']['Pec'][0]
    QBERI     = main_params['iter']['QBERI'][0]

    Err = None

    outfile = main_params['out']['out_base'] + \
            '_th_m_{:5.2f}_Pec_{}_QBERI_{}'.format(np.degrees(theta_max),
                                                    Pec,QBERI)


    if tOptimise:
        # Find the optimised SKL and protocol parameters
        fulldata, optdata, Err = \
            SKL_opt_loop_loss_and_time(theta_max,Pec,QBERI,dt,tInit,x,
                                    mu3,xb,args_fixed,bounds,cons,
                                    options,opt_params,fulldata,optdata,
                                    sys_params,Syseff,sysLoss)
    else:
        # Calculate the SKL for a given set of protocol parameters
        fulldata = SKL_loop_loss_and_time(theta_max,Pec,QBERI,x,dt,mu3,args_fixed,
                                     fulldata,sys_params,Syseff,sysLoss)
    
    #*******************************************************************
    #*******************************************************************
    # Sort and output data
    #*******************************************************************
    #*******************************************************************
    if (main_params['out']['tWriteFiles']):
        multidata = write_data(main_params['out'],tOptimise,
                                header,opt_head,outfile,fulldata,
                                multidata,optdata)
    return multidata, Err
  
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def sys_param_list(fixed_params,max_elev):
    """
    Generate a list of fixed system parameters.

    Parameters
    ----------
    fixed_params : dict
        Dictionary of fixed arguments required by key_length functions.
    max_elev : float
        Maximum elevation of the staellite overpass (rad).

    Returns
    -------
    sys_params : dict
        Dictionary of fixed system parameters.

    """
    sys_params = []
    sys_params.append(max_elev)                # The maximum elevation of the satellite overpass
    sys_params.append(fixed_params['eps_c'])   # The correctness parameter
    sys_params.append(fixed_params['eps_s'])   # The secrecy parameter
    sys_params.append(fixed_params['Pap'])     # Probability of an afterpulse event
    sys_params.append(fixed_params['NoPass'])  # Number of satellite overpasses
    sys_params.append(fixed_params['Rrate'])   # Source repetition rate
    sys_params.append(fixed_params['minElev']) # The minimum satellite elevation used for transmission
    sys_params.append(fixed_params['shift0'])  # Number of time steps to shift the t=0 point from zenith
    return sys_params

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def args_fixed_list(fixed_params,calc_params,FSeff,time0pos):
    """
    Extract fixed parameters from a dictionary and return as a list.

    Parameters
    ----------
    fixed_params : dict
        Dictionary of fixed system parameters.
    calc_params : dict
        Dictionary of calculation parameters.
    FSeff : float, array-like
        The free space transmission losses (as efficiency) as a function of time.
    time0pos : int
        Index of the t=0 point of the transmission window in the FSeff array.

    Returns
    -------
    list
        Fixed system parameters.

    """
    mu3         = fixed_params['mu3']     # Intensity of pulse 3 (vacuum)
    Pap         = fixed_params['Pap']     # Probability of an afterpulse event
    rate        = fixed_params['Rrate']   # Repetition rate of the source in Hz
    pbx         = fixed_params['Pbx']     # Probability Bob measures an X basis signal
    eps_c       = fixed_params['eps_c']   # The correctness parameter
    eps_s       = fixed_params['eps_s']   # The secrecy parameter
    NoPass      = fixed_params['NoPass']  # Number of satellite overpasses
    boundFunc   = calc_params['bounds']   # Type of tail bounds
    errcorrFunc = calc_params['funcEC']   # Error correction function
    fEC         = calc_params['fEC']      # Error correction factor   
    num_zero    = calc_params['num_zero'] # Small value used to approximate zero
    return [mu3,rate,Pap,boundFunc,pbx,eps_c,eps_s,num_zero,
            errcorrFunc,fEC,NoPass]

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def SKL_main_loop(protocol,main_params,adv_params,x,xb,f_atm,bounds,cons,
                  options,header,opt_head,fulldata,optdata,multidata):
    """
    Calculate the SKL using either optimised or specified protocol parameters
    for iterated values of:
        theta_max : Maximum elevation angle (rad)
        Pec : Probability of extraneous counts
        QBERI : Intrinsic quantum bit error rate
        ls : Excess system loss (dB)
        dt: Transmission window half-width (s)

    Parameters
    ----------
    main_params : dict
        Dictionary of general calculation parameters.
    adv_params : dict
        Dictionary of advanced calculation parameters.
    x : float, array-like
        The initial set of protocol parameters to use.
    x0i : float, array-like
        The final protocol parameters from previous calculations.
    xb : float, array-like
        The upper and lower bounds for the optimised parameters.
    ci : int, array-like
        Calculation loop counter.
    ni : int, array-like
        The number of iterations in each loop.
    f_atm : function
        Atmospheric transmissivity vs elevation angle (rad) function.
    bounds : obj
        Scipy bounds object containing optimised parameter bounds.
    cons : obj
        Scipy object containing the optimisation constraints.
    options : dict
        Dictionary of optimiser parameters.
    header : str
        Comma separated string of output parameters corresponding to columns
        in fulldata. Used as header when writing output to file.
    opt_head : str
        Comma separated string of optimiser metrics and parameters corresponding
        to columns of optdata. Used as a header when writing output to file.
    fulldata : float, array-like
        Calculation output parameters.
    optdata : list
        Optimisation metrics and parameters.
    multidata : float, array-like
        Optimal dt output parameters for all calculations.

    Returns
    -------
    None.

    """
    # Retrieve loss data array
    loss_data = get_losses(main_params['iter']['theta_max'],
                            main_params['loss'],f_atm,
                            main_params['out']['tPrint'],
                            main_params['out']['out_path'])

    time = loss_data[:,0]
    dt = [np.abs(time[i]-time[i-1]) if i >= len(time)-1 else np.abs(time[i+1]-time[i])
             for i in range(len(time))]

    # Free space efficiency and loss in db
    # Returns the (t,) array FSeff, where t is the total number of time-slots
    Syseff  = loss_data[:,2]
    Sysloss = -10*np.log10(Syseff)

    theta_max = main_params['iter']['theta_max']

    # Find the time slot at the centre of the pass where t = 0.
    time0pos   = np.where(loss_data[:,0] == 0)[0][0]
    time0elev  = loss_data[time0pos,1] # Elevation angle at t = 0 (rads).
    
    # Maximum elevation angle (degs) of satellite pass
    max_elev = np.degrees(time0elev)

    min_elev   = main_params['fixed']['minElev'] # Reset value
    # Find first dt value corresponding to an elevation greater than min_elev
    minElevpos = np.where(loss_data[:,1] >= np.radians(min_elev))[0][0] # Index of first value
    dt_elev    = loss_data[minElevpos,0] # Max value of dt less than, or equal to, the

    # Store list of fixed params for output data file
    sys_params = sys_param_list(main_params['fixed'],max_elev)
    # Store list of optimiser arguments
    args_fixed = args_fixed_list(main_params['fixed'],adv_params['calc'],
                                    Syseff,time0pos)

    global x0_random
    global set_parameters
    global key_length_func 
    global key_length_inverse  
    global key_length_simetric
    global arrange_out

    if protocol.lower() == list_protocols[0].lower():
        from key.protocols.fourstatetwodecoy.init_4s2d import (x0_rand,arrange_output)
        from key.protocols.fourstatetwodecoy.key_4s2d import (set_params, key_length, key_length_inv, key_length_sim)
    elif protocol.lower() == list_protocols[1].lower():
        from key.protocols.threestateonedecoy.init_3s1d import (x0_rand,arrange_output)
        from key.protocols.threestateonedecoy.key_3s1d import (set_params, key_length, key_length_inv, key_length_sim)

    x0_random           = x0_rand
    set_parameters      = set_params
    key_length_func     = key_length
    key_length_inverse  = key_length_inv
    key_length_simetric = key_length_sim
    arrange_out         = arrange_output

    # Run main SKL loops
    multidata, Err = SKL_sys_loop(x,xb,theta_max[0],main_params,
                        adv_params['opt'],args_fixed,bounds,cons,options,
                        header,opt_head,fulldata,optdata,multidata,
                        sys_params,Syseff,Sysloss,dt)

    return Err