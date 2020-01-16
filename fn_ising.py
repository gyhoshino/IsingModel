import os
import time
import csv
import numpy as np
import logging
import curses
import pickle
from glob import glob
'''This is a file with methods (and classes) called by mod_ising.py.
You are welcome to modify these. However, it is not anticipated that you
should need to do so.
'''

# colors for the ncurses screen
'''Up ()  down []'''
UP_FG = 16
UP_BG = 205
DN_FG = 16
DN_BG = 63
i_UP  = 10
i_DN  = 11

def set_pickle_dir(inp):
    '''
    Make the pickle directory under data_dir/output_dir/pickles
    data_dir is set by user
    output_dir may be set by user. If not, make it.
    No return value.
    '''
    data_dir = inp['data_dir']
    name = data_dir+'/'

    if not len(inp['output_dir']) == 0:
        name += inp['output_dir']
    else:
        if inp['output_prefix']:
            name += inp['output_prefix']+'_'
        name += f"s%i_N{inp['N']}"
        if inp['output_date']:
            name += str(time.strftime("_%Y%m%d-%H%M%S"))

        set_number = 0
        while os.path.isdir(name % set_number):
            set_number += 1

        name = name % set_number

    inp['output_dir'] = name

    if not os.path.isdir(inp['output_dir']):
        os.makedirs(inp['output_dir'])

    inp['pickles'] = inp['output_dir']+'/'+'pickles'

    if not os.path.isdir(inp['pickles']):
        os.makedirs(inp['pickles'])

        with open(inp['pickles']+'/readme','w') as fout:
            fout.write('''
This directory contains Python3 pickles of run output for various temperatures.
The name of each pickle is the set and temperature:
    exmaple:
        s0_2.3.pickle for the first time T=2.3 was run
        s1_2.3.pickle for the second time T=2.3 was run
        etc...
The pickle contians either summaries of the EM and Spin Correlation data,
or the full EM and spin corr data sets.

These objects can be unpickled by Python and compiled into CSV files.
Use the Python script unpickle_to_csv.py to do this.
''')

def get_pickle_name(inp,T):
    '''Return pickle name like: s0_2.345.pickle, for T=2.345, set 0. 
    Increment s0 up to sN to find a unique name.'''
    i = 0
    name = inp['pickles']+f'/s{i}_{T:.9}.pickle'
    while os.path.isdir(name):
        i += 1
        name = inp['pickles']+f'/s{i}_{T:.9}.pickle'
    return name

# def get_filenames(inp): #make data folder if doesn't exist, then specify filename
#     '''Generate the output file names for the EM (energy and megnetism) and SC (spin correlation) files'''
#     try:
#         dir_out = inp['dir_out']
#         prefix  = inp['file_prefix']
#         if inp['date_output']:
#             dir_out += str(time.strftime("_%Y%m%d-%H%M%S"))

#         if not os.path.isdir(dir_out):
#             os.makedirs(dir_out)

#         # file name = [file_prefix]##_EM_v#.csv if only one temperature (example: runA_4.20_EM_v0.csv)
#         #             [file_prefix]##T##_EM_v#.csv if there are two temperatures (example: runA_4.2T5.3_EM_v0.csv) 
#         # the other file name is identical, but with "SC" (for spin correlation)) instead of EM
#         if inp['T_max'] <= inp['T_min']:
#             t_name = '%.4f'%inp['T_min']
#         else:
#             t_name = '%.4fT%.2f'%(inp['T_min'],inp['T_max'])

#         v = 0
#         while True:
#             EM_file = f'{dir_out}/{prefix}{t_name}_EM_v{v}.csv'
#             SC_file = f'{dir_out}/{prefix}{t_name}_SC_v{v}.csv'
#             if not (os.path.isfile(EM_file) or os.path.isfile(SC_file)):
#                 break
#             v += 1
#         return EM_file, SC_file

#     except:
#         print ('fatal: Failed to make output file names')
#         sys.exit()

# def get_pickle_jar_name(inp):
#     if 'pickle_jar' in inp:
#         return inp['pickle_jar']
#     else:
#         return f"jar_{inp['N']}"

# def pickle_results(inp, T, E, M):
#     data = {}
#     if not os.path.isdir(jar_dir):
#         os.mkdir(jar_dir)
#     if os.path.isfile(pickle_name):
#         data = pickle.load(open(pickle_name,'rb'))
#         data['nEM']  += len(E)
#         data['sumE'] += np.sum(E)
#         data['sumM'] += np.sum(M)
        
#     else:
#         data = {'n_EM':0,
#                 'sumE':3.4}
     

def parse_cmd_line_args(inp, cmd_line_args):
    '''Parse arguments from the command line into input dictionary inp'''
    for x in cmd_line_args[1:]:
        if ':' in x:
            try:
                key, val = x.split(':')
                if key not in inp:
                    print(f'fatal error: invalid command-line input "{key}"')
                    print('  If you want to add a new input, add it in\n'
                          '  the "set_input" function in mod_ising.py')
                    raise KeyError
                try:
                    if '.' in val:
                        inp[key] = float(val)
                        print('%-20s'%('inp["%s"]'%key),'set to float  ',inp[key])
                    elif val.lower() == 'false' or val.lower() == 'f':
                        inp[key] = False
                    elif val.lower() == 'true' or val.lower() == 't':
                        inp[key] = True
                    else:
                        inp[key] = int(val)
                        print('%-20s'%('inp["%s"]'%key),'set to int    ',inp[key])
                except:
                    inp[key] = val
                    print('%-20s'%('inp["%s"]'%key),'set to string ',inp[key])
            except KeyError:
                exit(2)
            except:
                print('warning: input "%s" not added to arguments'%x)
        else:
            print(f'fatal error: invalid command-line format "{x}"')
            print( '             proper format is key:value key:value etc...')
            exit(2)

class checkered_nums:
    '''Class to write increasing digits (always two columns limit) to column and row.
       Print in checkered colors: black-on-white then white-on-black'''
    def __init__(self, stdscr_in):
        self.stdscr = stdscr_in
        self.i = 1
    def print(self, row, col):
        if col+2 > curses.COLS or row >= curses.LINES:
            return
        num = self.i
        if num > 100:
            num = num % 100
        num = '%2i'%num
        # self.stdscr.addstr(col,row,'%i'%num, curses.color_pair(5))
        self.stdscr.addstr(row,col,num,curses.color_pair(5+self.i%2))
        self.i += 1

class printer:
    '''class to print to the screen'''
    def __init__(self, inp, stdscr=None):
        ''' Use stdscr (curses terminal screen) if present.  Otherwise use stdout ''' 
        self.start_time = time.time()
        self.inp = inp
        self.stdscr = stdscr

        if self.stdscr:
            self.col = 0
            self.row = 1
            curses.start_color()
            curses.use_default_colors()
            curses.init_pair(1, curses.COLOR_RED,   254) #curses.COLOR_WHITE)
            curses.init_pair(2, curses.COLOR_WHITE, curses.COLOR_MAGENTA)
            curses.init_pair(3, curses.COLOR_GREEN, 16) #curses.COLOR_BLACK)
            curses.init_pair(4, curses.COLOR_WHITE, curses.COLOR_CYAN)
            curses.init_pair(5, curses.COLOR_WHITE, 245) #16) #curses.COLOR_BLACK)
            curses.init_pair(6, curses.COLOR_BLACK, curses.COLOR_WHITE)
            curses.init_pair(7, curses.COLOR_RED, 16)
            curses.init_pair(8, curses.COLOR_WHITE, 16)
            curses.init_pair(i_UP, UP_FG, UP_BG)
            curses.init_pair(i_DN, DN_FG, DN_BG)

    def start_singlecorr(self):
        '''Write start of program message to screen'''
        if self.stdscr:
            self.stdscr.addstr(0,0, f'{self.inp["N"]}', curses.color_pair(7))
            self.stdscr.addstr('x')
            self.stdscr.addstr(f'{self.inp["N"]}', curses.color_pair(7))
            self.stdscr.addstr(' 2D Ising Lattice: Single Core')
                    # 2D Ising Lattice: Single Core',
                   # curses.color_pair(7) )
            self.stdscr.refresh()
        else:
            print(f'{self.inp["N"]}x{self.inp["N"]} 2D Ising Lattice: Single Core')

    def start_Tprogress(self, Tfinal):
        # get number of finished steps
        '''Write initial progress line to screen for a given temperature'''
        self.Tfinal = Tfinal;
        self.Tstart_time = time.time()
        self.fmt_prog = ('T:{Tfinal:>5.4} '
                'steps: {prog[0]:>7}/{prog[1]:>7}, {ratio:3}% '
                   'time: {time_pass} '
                   'est.time-to-go: {est_time}'
                   )
        self.fmt_finish = ('T:{Tfinal:>.2} '
                   'steps: {prog[1]:>7}/{prog[1]:>7}, 100% '
                   'time: {time_pass} '
                   'est.time-to-go: done!'
                   )
        if self.stdscr:
            pass
        else:
            print(f'T:{self.Tfinal:>.2}',end='\r')

    def update_Tprogress(self, prog,finished=False):
        time_pass = time.time() - self.Tstart_time 
        time_pass_str = time.strftime('%H:%M:%S', time.gmtime(time_pass));

        # if finished is True, prog might = False (if the iteration
        # didn't end on an update)
        if prog:
            self.last_prog = prog
        else:
            prog = (self.last_prog[1],self.last_prog[1])

        try:
            ratio = float(prog[0])/prog[1]
        except:
            print(f'\n MY ERROR prog: {prog}')
        est_time  = (1-ratio)*time_pass/ratio

        est_time_str = time.strftime('%H:%M:%S', time.gmtime(est_time));
        steps = list(prog)
        if finished:
            est_time_str = 'done!    '
            steps[0] = steps[1]
        msg = self.fmt_prog.format(prog=steps,ratio=(int(100*ratio)),
                                   time_pass=time_pass_str,
                                   est_time=est_time_str,
                                   Tfinal=self.Tfinal)
        if self.stdscr:
            self.stdscr.addstr(self.row,0, msg)
            self.stdscr.refresh()
            if finished:
                self.row_plus1()
        else:
            if finished:
                print(msg)
            else:
                print(msg,end='\r')
    def finish_Tprogress(self):
        pass

    def row_plus1(self):
        '''progress the column count. Loop if necessary'''
        self.row += 1
        if self.row > curses.LINES-3:
            self.stdscr.addstr(0,0,'... output looping from bottom of screen ...')
            self.stdscr.clrtoeol()
            self.row = 1

    #             print(string)
    def print_done(self):
        time_pass_str = time.strftime('%H:%M:%S', time.gmtime(time.time()-self.start_time));
        if self.stdscr:
            self.stdscr.clrtobot()
            self.stdscr.refresh()
            self.stdscr.addstr(self.row,0,f'Total run time: {time_pass_str} seconds.')
            self.row_plus1()
            self.stdscr.addstr(self.row, 0, '  Finished ising model program. Press any key to continue.  ',
                    curses.color_pair(1))
            self.stdscr.clrtoeol()
            # curses.flash()
            self.stdscr.getkey()
        else:
            print('\n---- Program Finished ----')
            print(f'Total run time: {time_pass_str} seconds.')

    def draw_lattice(self, lattice, step=0, nstep=0, phase='',T='', B=''):
        '''
        Print the lattice.
        If not in curses windows, use lattice.print_spins().
        If in curses window, print as:
       0         1         2         3         4         5         6         7
       01234567890123456789012345678901234567890123456789012345678901234567890123456789
        200x200 Step   12000/100000   [Anneal/Burning/Analysis]
        T 1.34    B 9.13   |E| 12.34   |M| -0.123  
                             1 2 3 4 5 6 7
        0 1 0 1 0 1 0  or  1()[]()[]()[]()
        1 1 0 0 0 0 0      2()()()()()()()
        1 0 1 0 0 1 1      3()[][][][]()()
        1 0 1 0 1 0 1      4()[]()[]()[]()
        0 0 0 1 0 1 0      5[][][]()[]()[]
        1 0 1 0 0 1 0      6()[]()[][]()[]
        0 1 0 1 0 0 1      7[]()[]()[][]()
        '''
        if not self.stdscr:
            print('Not using n-curses screen. Will not print 2D Map')
            return

        # check for screen large enough
        N = self.inp['N']
        # if curses.COLS < 2*(N+1):
        #     self.stdscr.addstr(1, 1, 'Screen not wide enough to show 2D-map',
        #             curses.color_pair(1))
        #     return

        self.stdscr.clear()
        self.stdscr.addstr(0,0,curses.COLS*' ',curses.color_pair(3))
        self.stdscr.addstr(1,0,curses.COLS*' ',curses.color_pair(3))
        self.stdscr.addstr(0,1,'%ix%i'%(N,N), curses.color_pair(3))
        if step and nstep:
            self.stdscr.addstr(0,9,'Step %i/%i'%(step,nstep), curses.color_pair(3))
        if phase:
            self.stdscr.addstr(0,31,phase, curses.color_pair(3))
        if T:
            self.stdscr.addstr(1,1,'T %.2f'%T, curses.color_pair(3))
        if B:
            self.stdscr.addstr(1,11,'B %.2f'%B, curses.color_pair(3))

        self.stdscr.addstr(1,20,'|E| %.2f'%lattice.get_E(),curses.color_pair(3))
        self.stdscr.addstr(1,32,'|M| %.2f'%lattice.get_M(),curses.color_pair(3))
        values=lattice.get_numpy_spin_matrix()
        self.stdscr.addstr(3,0,'%i %i'%(np.min(curses.LINES-2),N))

        # plot column indices
        row_labels = checkered_nums(self.stdscr)
        col_labels = checkered_nums(self.stdscr)
        flip = 0
        for row in range(np.min([curses.LINES-3,N])):
            flip += 1
            # if flip%2:
            #     curses.init_pair(2,curses.COLOR_RED, curses.COLOR_WHITE)
            # else:
            #     curses.init_pair(2,curses.COLOR_BLUE, curses.COLOR_BLACK)
            row_labels.print(row+3,0)
            for col in range(np.min([int(curses.COLS/2)-2, N])):
                if row == 1:
                    col_labels.print(2, col*2+2)
                if values[row,col]==1:
                    self.stdscr.addstr(row+3,col*2+2,'()',curses.color_pair(i_UP))
                else:
                    self.stdscr.addstr(row+3,col*2+2,'[]',curses.color_pair(i_DN))

        self.stdscr.addstr(curses.LINES-1, 0, 'press a key to continue',
                curses.color_pair(1))
        # self.stdscr.addstr(5,0,'%s'%pairs)
        self.stdscr.refresh()
        self.stdscr.getkey()  

def plot_graphs(data): #T,E_mean,E_std,M_mean,M_std): #plot graphs at end
    dat = np.array(data)
    plt.figure(1)
    # plt.ylim(0,1)
    plt.errorbar(dat[:,0], dat[:,1], yerr=dat[:,2], fmt='o')
    plt.xlabel('Temperature')
    plt.ylabel('Average Site Energy')
    plt.figure(2)
    plt.errorbar(dat[:,0], np.absolute(dat[:,3]), yerr=dat[:,4], uplims=True, lolims=True,fmt='o')
    plt.xlabel('Temperature')
    plt.ylabel('Aveage Site Magnetization')
    plt.show()

def unpickle_to_csv(output_dir):
    print(f'\nWriting cvs files in {output_dir}/\n')
    with open(f'{output_dir}/inp_readme.csv','w') as readme:
     with open(f'{output_dir}/EM.csv','w') as EM:
      with open(f'{output_dir}/SpinCorr.csv','w') as SC:
        readme.write("Making CSV files EM.csv and SpinCorr.csv from the data in pickles/*.pickle\n")
        # readme.write("Note that for now, if there are more than one file for each temperature, then the most recent file only is read.\n")
        readme.write("Options from 'inp' in each pickle are given here.\n") 
        readme_writer = csv.writer(readme, delimiter=',', lineterminator='\n')
        EM_writer = csv.writer(EM, delimiter=',', lineterminator='\n')
        SC_writer = csv.writer(SC, delimiter=',', lineterminator='\n')

        T_vals = dict()
        for pfile in glob(f'{output_dir}/pickles/*.pickle'):
            # make a list of temperatures and sort them
            val = float(pfile.split('_')[-1][:-7])
            T_vals[val] = pfile
        keys = list(T_vals.keys())
        keys.sort()

    # EMata = []
    # SCata = []

        EM.write('E and M values from pickles/*.pickle\n')
        EM.write('For inp values for each row see inp_readme.csv\n' '\n')
        EM_writer.writerow('Temp n_samples E_mean E_std M_mean M_std E_stdstd M_stdstd'.split())

        SC.write('Spin Correlations values from pickles/*.pickle\n')
        SC.write('For inp values for each row see inp_readme.csv\n' '\n')

        first_entry = True
        
        for T in keys:
            # print ("t: " , T, ' ' , output_dir)
            with open(T_vals[T],'rb') as pickle_file:
                data = pickle.load(pickle_file)
                inp = data['inp']

                first_keys = 'N T_min T_max T_spacing r_flip EM_samples SC_samples'.split()
                first_vals = [inp[K] for K in first_keys]
                rest_keys = [K for K in inp.keys() if K not in first_keys]
                rest_keys.sort()
                rest_vals = [inp[K] for K in rest_keys]

                readme_writer.writerow(('File:', T_vals[T], 'Temperature:', T))
                readme_writer.writerow(first_keys+rest_keys)
                readme_writer.writerow(first_vals+rest_vals)
                readme_writer.writerow(tuple(' '))

                if first_entry:
                    first_entry = False
                    SC_writer.writerow(['Temp','n_samples']
                                   +['mean_d%i'%i for i in range(1,len(data['SC_mean'])+1)]
                                   +['std_d%i'%i  for i in range(1,len(data['SC_mean'])+1)]
                    )

                # write the data 
                EM_writer.writerow( tuple((T, data['E_len'], data['E_mean'], data['E_std'], data['M_mean'], data['M_std'], data['E_stdstd'], data['M_stdstd'])) )
                SC_writer.writerow( np.hstack((T, data['SC_len'], data['SC_mean'], data['SC_std'])) )

                if 'E_data' in data:
                    data['E_data'].tofile(f'{output_dir}/{T}_E_vals.csv','\n','%s')
                if 'M_data' in data:
                    data['M_data'].tofile(f'{output_dir}/{T}_M_vals.csv','\n','%s')

                # ok, this needs to be done, but not just now...
                # if 'SC_data' in data:
                    # with open (f'{output_dir}/{T}_SC_vals.csv','w') as f_SC:
                        # writer = csv.writer
                        # for SC in data['SC_data']:
                    # data['SC_data'].tofile(f'{output_dir}/{T}_SC_vals.csv',',','%s')
            
                
                #print the inp to the readme file
                # EM_data.append( tuple(T, data['E_len'], data['E_mean'], data['E_std'], data['M_mean'], data['M_std']) )
                # SC_data.append( tuple(T, data['SC_len'], data['SC_mean'], data['SC_std']) )

        # assume that the std 
                # print("Temp: " , T, " E_mean: ",data['E_mean'])
    

# def print_results(inp, EM_data, SC_data):
#     data_filename, corr_filename = get_filenames(inp)
#     with open(data_filename,'w') as f_out:
#         writer = csv.writer(f_out, delimiter=',', lineterminator='\n')

#         # first echo all the 'inp' dictionary, with a few 'important' settings first
#         # followed by the remaining settings
#         first_keys = 'N T_min T_max T_spacing r_flip EM_samples SC_samples'.split()
#         first_vals = [inp[K] for K in first_keys]
#         rest_keys = [K for K in inp.keys() if K not in first_keys]
#         rest_keys.remove('print') # print is a method passed by convenience. Don't print it.
#         rest_keys.sort()
#         rest_vals = [inp[K] for K in rest_keys]
#         writer.writerow(first_keys+rest_keys)
#         writer.writerow(first_vals+rest_vals)

#         # write an empty row
#         writer.writerow([' '])
#         writer.writerow('Temp n_samples E_mean E_std M_mean M_std'.split())
#         for entry in EM_data:
#             writer.writerow(entry)

#     with open(corr_filename,'w') as f_out:
#         writer = csv.writer(f_out, delimiter=',', lineterminator='\n')
#         first_keys = 'N T_min T_max T_spacing r_flip EM_samples SC_samples'.split()
#         first_vals = [inp[K] for K in first_keys]
#         rest_keys = [K for K in inp.keys() if K not in first_keys]
#         rest_keys.remove('print')
#         rest_keys.sort()
#         rest_vals = [inp[K] for K in rest_keys]
#         writer.writerow(first_keys+rest_keys)
#         writer.writerow(first_vals+rest_vals)

#         writer.writerow([' '])
#         writer.writerow(['Temp','n_samples']
#                        +['mean_d%i'%i for i in range(1,len(SC_data[0][2])+1)]
#                        +['std_d%i'%i  for i in range(1,len(SC_data[0][2])+1)]
#         )
#         for entry in SC_data:
#             writer.writerow(np.hstack(entry))

