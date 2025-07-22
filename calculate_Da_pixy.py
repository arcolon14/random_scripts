#!/usr/bin/env python3
import sys, os, argparse
from datetime import datetime

# Some constants
PROG = sys.argv[0].split('/')[-1]
DEF_PREFIX = 'pixy'

def parse_args(prog=PROG):
    '''Set and verify command line options.'''
    p = argparse.ArgumentParser()
    p.add_argument('-p', '--pi', required=True, 
                   help='(str) Path to the pi output from Pixy.')
    p.add_argument('-d', '--dxy', required=True, 
                   help='(str) Path to the dxy output from Pixy.')
    p.add_argument('-o', '--out-dir', required=False, default='.',
                help='(str) Path to output directory [default=./].')
    p.add_argument('-u', '--out-prefix', required=False, default=DEF_PREFIX,
                   help=f'(str) Optional prefix for output file(s). [default={DEF_PREFIX}]')
    # Check inputs
    args = p.parse_args()
    assert os.path.exists(args.pi)
    assert os.path.exists(args.dxy)
    args.out_dir = args.out_dir.rstrip('/')
    assert os.path.exists(args.out_dir)
    # Set some constants
    return args

def date() -> str:
    '''Print the current date in YYYY-MM-DD format.'''
    return datetime.now().strftime("%Y-%m-%d")

def time() -> str:
    '''Print the current time in HH:MM:SS format.'''
    return datetime.now().strftime("%H:%M:%S")

class PiWindow:
    '''
    Class to store the Pi values from Pixy.
    '''
    def __init__(self, chromosome:str, window_pop_1:int, window_pop_2:int,
                 avg_pi:float, population:str):
        self.chrom = chromosome
        self.start = window_pop_1
        self.end   = window_pop_2
        self.pi    = avg_pi
        self.pop   = population
    def __str__(self):
        return f'{self.pop}\t{self.chrom}\t{self.start}\t{self.end}\t{self.pi}'

def parse_pi(pi_f:str)->dict:
    '''
    Parse a the Pi output file from Pixy.
    Args:
        pi_f: (str) Path to the Pixy Pi file.
    Returns:
        pi_vals: (dict) Dictionary of Pi values per population, per window.
            pi_vals = { pop1 : { window1 : PiWindow, window2 : PiWindow }, 
                        pop2 : { window1 : PiWindow, window2 : PiWindow },
                        pop3 : { window1 : PiWindow, window2 : PiWindow } }
    '''
    print(f'\nLoading Pi values from file:\n    {pi_f}', flush=True)
    pi_vals = dict()
    pop_tally = dict()
    records = 0
    with open(pi_f) as fh:
        for i, line in enumerate(fh):
            line = line.strip('\n')
            if line.startswith('#') or len(line)==0:
                continue
            # Skip the header line
            if line.startswith('pop\tchromosome\twindow_pos'):
                continue
            # Process the rest of the records
            fields = line.split('\t')
            records += 1
            pop = fields[0]
            chromosome = fields[1]
            window_pos_1 = int(fields[2])
            window_pos_2 = int(fields[3])
            # Create a window id based on chromosome and positions.
            # This will be used to match windows across different files.
            window_id = f'{chromosome}:{window_pos_1}-{window_pos_2}'
            # Pi values can be NAs, so handle accordingly
            avg_pi = fields[4]
            try:
                avg_pi = float(fields[4])
            except ValueError:
                avg_pi = None
            # Store all the needed values as a PiWindow
            pi_window = PiWindow(chromosome, window_pos_1, 
                                 window_pos_2, avg_pi, pop)
            # Add this to the output dictionary
            pi_vals.setdefault(pop, dict())
            pi_vals[pop][window_id] = pi_window
            # Add to the tally dictionary, but only for 
            # non NA values.
            pop_tally.setdefault(pop, 0)
            if avg_pi is not None:
                pop_tally[pop] += 1
    # Report to log
    print(f'\n    Read {records:,} records from input Pi file.')
    print(f'    Loaded Pi values for {len(pop_tally):,} populations.')
    print(f'    Retained the following windows with non-NA values:')
    for pop in pop_tally:
        tally = pop_tally[pop]
        print(f'        {pop} : {tally:,}')
    return pi_vals

class DxyWindow:
    '''
    Class to store the Dxy values from Pixy.
    '''
    def __init__(self, chromosome:str, window_pop_1:int, window_pop_2:int,
                 avg_dxy:float, population1:str, population2:str):
        self.chrom = chromosome
        self.start = window_pop_1
        self.end   = window_pop_2
        self.dxy   = avg_dxy
        self.pop1  = population1
        self.pop2  = population2
    def __str__(self):
        return f'{self.pop1}\t{self.pop2}\t{self.chrom}\t{self.start}\t{self.end}\t{self.dxy}'

def parse_dxy(dxy_f:str)->dict:
    '''
    Parse a the Dxy output file from Pixy.
    Args:
        dxy_f: (str) Path to the Pixy Dxy file.
    Returns:
        dxy_vals: (dict) Dictionary of Dxy values per pairwise population comparison, per window.
            dxy_vals = { pop1-pop2 : { window1 : DxyWindow, window2 : DxyWindow },
                         pop1-pop3 : { window1 : DxyWindow, window2 : DxyWindow }, 
                         pop2-pop3 : { window1 : DxyWindow, window2 : DxyWindow } } 
    '''
    print(f'\nLoading dxy values from file:\n    {dxy_f}', flush=True)
    dxy_vals = dict()
    comp_tally = dict()
    records = 0
    with open(dxy_f) as fh:
        for i, line in enumerate(fh):
            line = line.strip('\n')
            if line.startswith('#') or len(line)==0:
                continue
            # Skip the header line
            if line.startswith('pop1\tpop2\tchromosome\twindow_pos'):
                continue
            # Process the rest of the records
            fields = line.split('\t')
            records += 1
            pop1 = fields[0]
            pop2 = fields[1]
            # Create a comparison ID values on the pop1 and pop2 values
            comparison = f'{pop1}-{pop2}'
            chromosome = fields[2]
            window_pos_1 = int(fields[3])
            window_pos_2 = int(fields[4])
            # Create a window id based on chromosome and positions.
            # This will be used to match windows across different files.
            window_id = f'{chromosome}:{window_pos_1}-{window_pos_2}'
            # Dxy values can be NAs, so handle accordingly
            avg_dxy = fields[5]
            try:
                avg_dxy = float(fields[5])
            except ValueError:
                avg_dxy = None
            # Create a DxyWindow object to store this window
            dxy_window = DxyWindow(chromosome, window_pos_1,
                                   window_pos_2, avg_dxy, pop1, pop2)
            # Add this to the output dictionary
            dxy_vals.setdefault(comparison, dict())
            dxy_vals[comparison][window_id] = dxy_window
            # Add to the tally dictionary, but only for 
            # non NA values.
            comp_tally.setdefault(comparison, 0)
            if avg_dxy is not None:
                comp_tally[comparison] += 1

    # Report to log
    print(f'\n    Read {records:,} records from input Dxy file.')
    print(f'    Loaded Dxy values for {len(comp_tally):,} pairwise population comparisons.')
    print(f'    Retained the following windows with non-NA values:')
    for pop in comp_tally:
        tally = comp_tally[pop]
        print(f'        {pop} : {tally:,}')
    return dxy_vals

def get_pi(pi_vals:dict, pop_id:str,
           window_id:str)->PiWindow:
    '''
    Extract a given Pi window from the pi_vals
    dictionary and report necessary errors.
    Args:
        pi_vals: (dict) Dictionary of Pi values per population, per window.
        pop_id: (str) ID of target population.
        window_id: (str) Coordinate IDs of target window.
    Returns:
        pop_pi: (PiWindow) Pi values for a given window.
    '''
    # Process population
    pop_pi_windows = pi_vals.get(pop_id, None)
    if pop_pi_windows is None:
        sys.exit(f'Error: population {pop_id} not found in Pi windows.')
    pop_pi = pop_pi_windows[window_id]
    if pop_pi is None:
        sys.exit(f'Error: Window {window_id} not found for population {pop_id} in Pi windows.')
    assert isinstance(pop_pi, PiWindow)
    return pop_pi

def calculate_da(d_xy:float, 
                 pi_x:float,
                 pi_y:float)->float:
    '''
    Calculate Da (Nei's unbiased genetic distance) for a given 
    window based on Dxy and Pi valuies.
        Da = Dxy_12-((Pi_1+Pi_2)/2)
    Args:
        d_xy: (float) Absolute divergence between populations X and Y
        pi_x: (float) Nucleotide diversity for population X
        pi_y: (float) Nucleotide diversity for population Y
    Returns:
        da: (float) Genetic distance between populations X and Y
    '''
    if d_xy is None:
        return None
    if pi_x is None:
        return None
    if pi_y is None:
        return None
    da = d_xy - ((pi_x + pi_y)/2)
    return da

def write_output(pi_vals:dict, dxy_vals:dict,
                 out_dir:str='.',
                 out_prefix:str=DEF_PREFIX)->None:
    f'''
    Calculate Da from the Pi and Dxy values and save to 
    the selected output file.
    Args:
        pi_vals: (dict) Dictionary of Pi values per population, per window.
        dxy_vals: (dict) Dictionary of Dxy values per pairwise population comparison, per window.
        out_dir: (str) Path to output directory [default=.].
        out_prefix: (str) Prefix for output files [default={DEF_PREFIX}].
    Returns:
        None
    '''
    kept = 0
    outf = f'{out_dir}/{out_prefix}_Da.txt'
    with open(outf, 'w') as fh:
        # Write the header of the output file
        header = 'pop1\tpop2\tchromosome\twindow_pos_1\twindow_pos_2\tavg_da\n'
        fh.write(header)
        # Loop over the Dxy windows first, since the have the 
        # pairwise populations that will be used to get the Pi
        # values later.
        for comparison in dxy_vals:
            dxy_windows = dxy_vals[comparison].keys()
            for window_id in dxy_windows:
                dxy_window = dxy_vals[comparison][window_id]
                assert isinstance(dxy_window, DxyWindow)
                # Sqelect the corresponding Pi windows based on the
                # population and window ID
                pop1 = dxy_window.pop1
                pop2 = dxy_window.pop2
                window_id = f'{dxy_window.chrom}:{dxy_window.start}-{dxy_window.end}'
                # Process population1
                pop1_pi = get_pi(pi_vals, pop1, window_id)
                # Process population2
                pop2_pi = get_pi(pi_vals, pop2, window_id)
                # Calculate Da
                da = calculate_da(dxy_window.dxy, pop1_pi.pi, pop2_pi.pi)
                da_str = None
                if da is None:
                    da_str = 'NA'
                else:
                    da_str = f'{da:0.8f}'
                    kept += 1
                # Write output
                chrom = dxy_window.chrom
                pos1 = dxy_window.start
                pos2 = dxy_window.end
                row = f'{pop1}\t{pop2}\t{chrom}\t{pos1}\t{pos2}\t{da_str}\n'
                fh.write(row)
    # Print to log.
    print(f'\nCalculated Da for {kept:,} windows with non-NA values.')

def main():
    print(f'{PROG} started on {date()} {time()}.')
    # Parse args
    args = parse_args()
    # Load the Pi values
    pi_vals = parse_pi(args.pi)
    # Load the Dxy values
    dxy_vals = parse_dxy(args.dxy)
    # Generate output
    write_output(pi_vals, dxy_vals, args.out_dir, args.out_prefix)

    # Done!
    print(f'\n{PROG} finished on {date()} {time()}.')


# Run Code
if __name__ == '__main__':
    main()
