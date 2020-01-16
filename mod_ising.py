#################################################################################
#
# This file contains the functions most likely to be modified by all students.
#
# See the lab handout for more details.
#
# First method:
#  set_input:          Define allowed input parameters and defaults
#
# There are five generators, each named 'gen_<name>' that provide 
#   five separate inputs for the main loop that calculated values at each 
#   final temperature. They are:
#  
#  gen_T: Tempature at each step
#  gen_B: Magnetic field at each step
#  gen_collect_EM:  When to sample the average E and M values of each site
#  gen_collect_SC:  When to sample the spin correlation
#  gen_prog_update: This is purely cosmetic and constrolls update outputs to screen
#                   If in doubt, it's ok to just always yield False
#
# The remaining functions are:
#  run_ising_lattice:  Run the simulation for each step
#  make_T_array:       Make the default array of temperatures on which
#                      to run the simulation. Default is evenly spaced from
#                      T_min to T_max with T_spacing
#
#################################################################################

import sys
import time
import numpy as np
import matplotlib.pyplot as plt
from   sys import exit, argv
from   library_interface import IsingLattice
from   fn_ising import *
import curses

def set_input(cmd_line_args, parse_cmd_line_args):
    """Parse command-line parameters into an input dictionary.
    Provide default values.
    Add new parameters as you wish.

    note: This fills a Python dictionary called 'inp' which is passed
          around throughout the Python part of the program.

          Feel free to add your own input parameters.
    
    input:   sys.argv
             use syntax of keyword:value on command line
    return:  dict[key] = value
    note:    Any value which can be turned into a number will be a float 
             if it has a '.', otherwise it will be an int.

    """

    inp = dict()
    #----------------------------------------
    # Lattice size: i.e. lattice = N^2 points
    #----------------------------------------
    inp['N']          = 20

    #-------------------------------------------------------------------------------------------------
    # For temperatures, either use:
    # (1) T_min, T_max, and T_spacing which will cause it to run through all temperatures in the range
    #  or 
    # (2) T_input_file and T_input_line, which will run the single temperature on line(T_input_line) 
    #     of file(T_input_file)
    #
    # The program will default to option(2) unless T_input_file != '', in which case it will try
    # to use option 2.
    #
    # The advantage to using option (2) is that slurm can be used to run jobs over all the temperarues
    # simultaneously.
    #-------------------------------------------------------------------------------------------------

    # Option (1) of above
    inp['T_min']      = 2.0    # minimum temperature -> Not used if reading temperatures from 
    inp['T_max']      = 4.3    # maximum temperature
    inp['T_spacing']  = 0.1    # step size from min to max temperature

    # Option (2) of above
    inp['T_input_file']  = ''
    inp['T_input_line']  = 0


    #----------------------------------------------
    # Options used for annealing and burnin
    #----------------------------------------------
    inp['T0_anneal']         = 4.0        # start temperature (arbitrary; feel free to change)
    inp['steps_anneal']      = int(10000)   # number of lattice steps in simulation
    inp['steps_burnin']      = int(10000)   # optional parameter, used as naive default

    #----------------------------------------------
    # Options for statistics collected
    #----------------------------------------------
    inp['EM_samples']        = int(10000)
    inp['EM_sample_spacing'] = 10
    inp['SC_samples']        = int(500)
    inp['SC_algorithm']      = 1 # 0 for legacy version, 1 for <ab>-<a><b>

    # Group the EM data into this many groups, take the std of each group, and then
    # take the std of these standard deviations.
    inp['EM_n_groups_stdstd'] = 100 # number of sub groups from which to calc std(std(EM))

    # strongly recommended to use option 1 here. The legacy version is undocumented.
    inp['SC_algorithm']      = 1 # 0 for legacy version, 1 for <ab>-<a><b>

    # Should probably use B=0, and then modify B in the gen_B function for magnetic annealing
    inp['B']                 = 0.0    # magnetic field strength
    inp['r_flip']            = 0.02   # Very cool result from PJ and Elizabeth
                                      # showing autocorrelation effects starting about 0.05 in a 100x100 matrix
                                      # recommended to ultimately use smaller value than 0.05

    #----------------------------------------------
    # output options
    #----------------------------------------------
    inp['data_dir']          = 'data' # output directory for all files
    inp['output_dir']        = ''     # directory for results of this run (defaults to name like
                                      # s2_N100 for 'Set 2'  'Lattice 100x100'
                                      # will use 'file_prefix' and 'run_dir_date' if they are set)
                                      # Important:
                                      # If set to a specific value, it will output new pickles to the same place
    # -- options used in picking output_dir name. These are ignored if output_dir is provided --
    inp['output_date']       = False
    inp['output_prefix']     = ''

    inp['print_inp']         = False  # optional flag to print input to command line
    inp['plots']             = False  # whether or not plots are generated

    inp['use_cpp']           = True   # the c++ implementation is *much* faster
    inp['multiprocess']      = True   # will be turned off automatically if using T_input_file
    inp['curses']            = True   # '                                                    '
    inp['draw_lattice']      = False  # '                                                    '

    #-------------------------------
    # options for what to pickle
    #-------------------------------
    inp['pickle_final_lattice_state'] = True 
    inp['pickle_full_E']          = False # warning: this may take a lot of memory
    inp['pickle_full_M']          = False # warning: this may take a lot of memory
    inp['pickle_full_SC']         = False # warning: this may take a lot of memory
    inp['unpickle_to_csv']        = True

    inp['seed_offset']                = 0 # Important to use when submiting many jobs at the same time
                                          # - if 0 will default to T_input_line if T_input_file is used 
                                          # - will set the seed exactly if < 0 (for an exact repeat run)

    #function call from support_methods
    parse_cmd_line_args(inp, cmd_line_args)
    set_pickle_dir(inp)

    if not inp['T_input_file'] == '':
        if not os.path.isfile(inp['T_input_file']):
            exit(f'Fatal error: cannot time T_input_file "{inp["T_input_file"]}"')
        try:
            lines = open(inp['T_input_file']).readlines()
            inp['T_min'] = float(lines[inp['T_input_line']])
            inp['T_max'] = inp['T_min'] - 1.0
        except:
            exit(f'Fatal error in trying to read input line {inp["T_input_line"]}')
        
        if inp["seed_offset"] == 0:
            inp["seed_offset"] = inp['T_input_line']

        inp['multiprocess']    = False
        inp['curses']          = False
        inp['draw_lattice']    = False
        inp['plots']           = False
        inp['unpickle_to_csv'] = False

    if inp['use_cpp'] == False and inp['SC_algorithm'] == 1:
        print('Warning: Not using C++ implementation of the Ising Lattice.')
        print('         Therefore defaulting the SC_algorithm to option 0 (vs 1).\n')
        inp['SC_algorithm'] = 0
        
    if inp['print_inp']:
        print('Printed list of input keys:')
        for key in sorted(inp.keys()):
            print('%-20s'%key,' ',inp[key])

    return inp


def gen_T(inp, T_final):
    '''Yield values of temperature.
    Default implementation: 
    (1) Start at T0_anneal and linearly drop to T_final in steps anneal
    (2) Yield T_final for all remaining steps.'''

    for T in np.linspace(inp['T0_anneal'],T_final,inp['steps_anneal']):
        yield T

    for T in range(inp['steps_burnin'] +
                   inp['EM_samples']*inp['EM_sample_spacing']):
        yield T_final

def gen_B(inp, T_final=None):
    '''Yield magnetic field values.
    Default implementation: Always yield inp['B'].'''

    # default annealing: 0.1 decreasing linearly with tempearture
    # to the final value of inp['B']

    B_init  = 0.1;
    B_final = inp['B'];

    for B in np.linspace(B_init,B_final,inp['steps_anneal']):
        yield B

    for B in range(inp['steps_burnin'] +
                   inp['EM_samples']*inp['EM_sample_spacing']):
        yield B_final

    # while True:
        # yield inp['B']

def gen_collect_EM(inp, T_final=None):
    '''Yield when to sample E and M values.
    Default implementation: Don't sample during steps_anneal
      and steps_burning. Then sample EM_samples separated by
      EM_sample_spacing-1 between each sample.'''
    for _ in range(inp['steps_anneal']+inp['steps_burnin']):
        yield False
    
    for _ in range(inp['EM_samples']):
        for __ in range(inp['EM_sample_spacing']-1):
            yield False
        yield True

def gen_collect_SC(inp, T_final=None):
    '''Yield when to collect Spin Correlation values.
    Default implementation: Collect EM_samples evenly spaced 
             throughout where E and M values are sampled.'''
    for _ in range(inp['steps_anneal']+inp['steps_burnin']):
        yield False
    
    space_corr = int( (inp['EM_samples']*inp['EM_sample_spacing'])
                      /inp['SC_samples'])
    if space_corr == 0:
        print('Fatal error.',
         'Collecting more SC_samples than EM_samples*EM_sample_spacing.')
        exit(1)
    if space_corr == 1:
        space_corr = 2

    for _ in range(inp['SC_samples']):
        for __ in range(space_corr-1):
            yield False
        yield True

    while True:
        yield False

def gen_prog_update(inp):
    '''Yield when to update the screen. 
    Do not use with multiprocess.'''

    #if using multiprocess, always return False
    if inp['multiprocess']:
        while True:
            yield False

    n_total = (inp['steps_anneal']+inp['steps_burnin']
              +inp['EM_samples']*inp['EM_sample_spacing'])
    n_spacing = int( n_total/100);
    # print (f'\n\n {n_spacing} \n')
    if n_spacing == 0:
        n_spacing = 1
    n_complete = 1

    while True:
        for _ in range(n_spacing-1):
            n_complete += 1
            yield False
            
        n_complete += 1
        yield (n_complete, n_total)

def run_ising_lattice(inp, T_final, seed_offset=0, updates=True):
    '''Run a 2-D Ising model on a lattice. Use the above
    five generators for loop values.
    Return mean and std values for Energy, Magnetization, and Spin Correlations
    '''
    seed_offset += inp['seed_offset']
    if updates:
        inp['print'].start_Tprogress(T_final)

    lattice = IsingLattice(inp, seed_offset)
    if inp['SC_algorithm'] == 0:
        get_SC = lattice.get_SC_v0;
    else:
        get_SC = lattice.get_SC_v1;


    try:
        E_collect = []
        M_collect = []
        SC_collect = []
        
        for T, B, sample_EM, sample_SC, prog_update in zip(
            gen_T(inp, T_final),
            gen_B(inp),
            gen_collect_EM(inp, T_final),
            gen_collect_SC(inp, T_final),
            gen_prog_update(inp)
        ):
            lattice.step(T,B)

            if sample_EM:
                E_collect.append(lattice.get_E())
                M_collect.append(lattice.get_M())
            if sample_SC:
                SC_collect.append( get_SC() )
            if updates and prog_update:
                inp['print'].update_Tprogress(prog_update)

        # Done with calculating all E, M, and SpinCorr data
        # The remainder of this function deals with displaying and saving
        # the results.
        if updates:
            inp['print'].update_Tprogress(prog_update,finished=True)
        
        if inp['draw_lattice'] and not inp['multiprocess']:
            inp['print'].draw_lattice(lattice, T=T_final)

        SC_data = np.array(SC_collect)
        E_data  = np.array(E_collect)
        M_data  = np.array(M_collect)

        # add the E_stdstd values
        np.random.shuffle(E_data) 
        np.random.shuffle(M_data) 

        n_groups = inp['EM_n_groups_stdstd']
        group_size = len(E_data)//n_groups

        # pickle the results
        pickle_data = {
            'T'   : T_final,
            'E_mean'   : E_data.mean(),
            'E_std'    : E_data.std(),
            'E_stdstd' : np.std(np.array( [ np.std(E_data[i*group_size:(i+1)*group_size]) for i in range(n_groups) ] ) ),
            'E_len'    : len(E_data),
            'M_mean'   : M_data.mean(),
            'M_std'    : M_data.std(),
            'M_stdstd' : np.std(np.array( [ np.std(M_data[i*group_size:(i+1)*group_size]) for i in range(n_groups) ] ) ),
            'M_len'    : M_data.std(),
            'SC_mean'  : SC_data.mean(0),
            'SC_std'   : SC_data.std(0),
            'SC_len'   : np.shape(SC_data)[0]
        }     

        # make a copy of inp data -- except the printer class
        inp_pickle = inp.copy()
        del inp_pickle['print']
        pickle_data['inp'] = inp_pickle

        if inp['pickle_final_lattice_state']:
            pickle_data['lattice_spins'] = lattice.get_numpy_spin_matrix()

        lattice.free_memory()

        if inp['pickle_full_E']:
            pickle_data['E_data']  = E_data
        if inp['pickle_full_M']:
            pickle_data['M_data']  = M_data
        if inp['pickle_full_SC']:
            pickle_data['SC_data'] = SC_data

        with open(get_pickle_name(inp,T_final),'wb') as f_out:
            pickle.dump(pickle_data, f_out)

    except KeyboardInterrupt:
        try:
            lattice.free_memory()
        except:
            pass
        print("\n\nProgram terminated by keyboard. Good Bye!")
        sys.exit()

def make_T_array(inp):
    '''Function that returns 'array' (may be singly valued) of temperatures
    for which to run the simluation.'''

    if inp['T_max'] <= inp['T_min']:
        return [inp['T_min'],]
    else:
        # for when T_max-T_min is evenly divisible by T_spacing. 
        return np.repeat(np.arange(inp['T_min'], inp['T_max']+1E-5, inp['T_spacing']), 2)
