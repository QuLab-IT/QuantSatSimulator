from os.path import join

from input.read_input import get_input
from input.write_input import (print_header, print_input)
from channel.atmosphere.atmos_data import get_f_atm
from optimize.optimiser import (set_constraints, opt_arrays, out_heads)
from key.key import SKL_main_loop
from output.outputs import (out_arrays, get_timings, format_time)

###############################################################################

# Set the inpufile path and filename(s);
# Path to input files, empty string means location of execution
inputpath    = 'inputfiles'

# File name for main input parameters
filename     = 'input.txt'

# File name for advanced/optional input parameters (optional)
filename_adv = 'input-adv.txt'

###############################################################################

# Print the SatQuMA header
print_header()

# Start clock and CPU timer
tc00, tp00 = get_timings()

# Error Variable: Error Mensages
Err = None

# Retrieve input parameters
protocol, main_params, adv_params = get_input(join(inputpath,filename),
                                              join(inputpath,filename_adv))

# Print input parameters
if main_params['out']['tPrint']:
    print('-'*80)
    print('Security protocol:',protocol)
    print('-'*80)
    print_input(main_params,'Main')
    print('-'*80)
    print_input(adv_params,'Adv')
    print('-'*80)
    
###############################################################################

# Get atmospheric function (if required), else a dummy function is returned
f_atm = get_f_atm(main_params['loss'])

# Get bounds/constraints for the optimiser
bounds, cons, options, x, xb = set_constraints(main_params['opt'],
                                               main_params['fixed'],
                                               adv_params['opt'],protocol)

# Generate output file headers based on input parameters
header, opt_head = out_heads(protocol,main_params['opt']['tOptimise'],
                             adv_params['opt']['method'])

# Initialise output data arrays - arrays are sized based on the length of the 
# headers generated above
fulldata, optdata, multidata = out_arrays(header,opt_head,
                                          main_params['opt']['tOptimise'],
                                          main_params['out']['tdtOptData'])

# Get key!!!
Err = SKL_main_loop(protocol,main_params,adv_params,x,xb,f_atm,bounds,cons,
                    options,header,opt_head,fulldata,optdata,multidata)

###############################################################################

# Stop clock and CPU timers
tc11, tp11 = get_timings()

tc   = tc11-tc00      # Calculation duration from clock
tp   = tp11-tp00      # Calculation duration from CPU

if len(Err) > 0:
    print('Error Messages:')
    for entry in Err:
        print('({}) > '.format(entry[1]) + entry[0])

print('\nFinal clock timer:',format_time(tc))
print('Final CPU timer:  ',format_time(tp),'\n')
print('All done!')