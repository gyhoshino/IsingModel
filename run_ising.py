import multiprocessing as mp
from time import sleep
from mod_ising import *
from time import process_time
import time

def multi_lattice(inp, T,seed_offset):
    print(f'Starting Temp {T:.3}')
    time_start = time.time()
    # rval = run_ising_lattice(inp,T,seed_offset, updates=False)
    run_ising_lattice(inp,T,seed_offset, updates=False)
    time_pass_str = time.strftime('%H:%M:%S', time.gmtime(time.time()-time_start));
    print(f'Finished Temp {T:.3}. Total time: {time_pass_str}')
    return True
    # return rval

def run_multi_core(inp):
    print("\n2D Ising Model Simulation; multi-core\n")
    T_array = make_T_array(inp)
    pool = mp.Pool(mp.cpu_count())

    jobs = []
    seed_offset = 1
    it_0 = time.time();
    for T in T_array:
        jobs.append(pool.apply_async(multi_lattice,(inp,T,seed_offset)))
        seed_offset += 1
    pool.close()
    pool.join()

def run_single_core(inp):
    # sequentially run through the desired temperatures and collect the output for each temperature
    # EM_data = []
    # SC_data = []
    inp['print'].start_singlecorr()
    for T in make_T_array(inp):
        run_ising_lattice(inp,T)

    # a1 lion
        # EM_vals, SC_vals = run_ising_lattice(inp, T)
        # EM_data.append( EM_vals )
        # SC_data.append( SC_vals )
    # print_results(inp, EM_data, SC_data)

    # if inp['plots']:
        # plot_graphs(EM_data)

class run_wrapper:
    '''Wrap execution either in, or out, of a python cursor terminal ("curses")'''
    def __init__(self,inp):
        self.inp = inp
    def run(self, stdscr=None):
        self.inp['print'] = printer(self.inp, stdscr)
        if self.inp['multiprocess']:
            run_multi_core(self.inp)
        if not self.inp['multiprocess']:
            run_single_core(self.inp)

        # if desired, un-pickle and make a csv file
        # note: For now: 
        #       If there is data from two separate runs in the pickle directory for a
        #       given temperatures, then either:
        #         (a) the full data is present, in which case the two numpy arrays are appended
        #         (b) only the most recent file is used.
        #       If is possible to improve this in the future, by saving number of samples, and sum data^2
        #       but this is not yet implemented.
        if inp['unpickle_to_csv']:
            unpickle_to_csv(inp['output_dir'])

            # a0 lion


            # print_results(inp, EM_data, SC_data)

        self.inp['print'].print_done()

if __name__ == "__main__":
    """Main program: run Ising Lattice here"""
    inp = set_input(argv, parse_cmd_line_args)
    # print('a0')
    print(inp)
    # print('a1')

    # using ncurses with multiprocess will cause the program to fail
    if inp['multiprocess']:
        inp['curses'] = False

    wrap = run_wrapper(inp)
    
    if inp['curses']:
        #this is the call that provides stdscr to run_wrapper.run
        curses.wrapper(wrap.run) 
    else:
        wrap.run()
