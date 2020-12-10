import os, shutil, glob, re
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from subprocess import Popen

'''
In the paper they suggest to use a polynomial expression to counter the
effect of TDS of the form alpha = a*x^2 + b*x^3 as a function of
sin(theta)/lambda. The program will set-up a 12 scale factor refinement
for a given set of .mas, .inp and .hkl files and refine three cycles.
After every cycle the polynomial is fitted to the refined scale factors,
optimal parameters for a, b are determined, the polynomial with those
parameters is used to write a new ‘corrected’ .hkl file. This
nested-intervals approach is necessary as the model needs to relax to
the new more correct data but three cycles seemed sufficient, the
problem is however that the overall correction factor is not directly
obtainable. Therefore, the overall sum of the iteratively applied a, b
parameters are then used to write the final .hkl file and a final
refinement is performed. The iterative correction differs marginally
from the combined correction, however with the combined approach we can
report the actual parameters used to correct the data.

How to:
All you need to do is put your latest .inp, .hkl and .mas files into a
folder together with the TDS.py program and run the program using python3.
The Masterfile should be set-up to at least refine the non-Hydrogen Uij
and monopoles of all atoms, multipoles and Kappa would help in determining
better parameters but take a lot more time to refine, please perform
different TDS refinements and check for consistency!
The parameters a and b used to correct an hkl file ARE WRITTEN TO THE FILE
and the line: ! TDS CORRECTION FACTOR: a=x.xxx, b=x.xxx
is added and should remain there for later reference!
XD will correctly ignore the line.

One big note:
It seems that there is an issue with the very high resolution data as
for that data the Fo/Fc vs. resolution distribution does not follow the
expected TDS slope, e.g. the curve falls off towards unity again
instead of keeping to increase continuously. I had to cut the data (in
the .mas file:
SKIP   obs  0.0E+00 0.1E+11 *sigobs  3. 0.1E+07 *sinthl  0.000 1.200)
and I am not sure if this is needed for your remaining data as well.
Please check the DRK plot for the slope and set the sintl cut-off to
where (if) the slope starts to decrease.

Output (after the program completed!):
xd00.* -> Initial 1-scale refinement
xd01.* -> Initial 12-scale refinement
xd02.* -> 12-scale cycles
xd0x.* -> cycles
xd.* -> final corrected .hkl and 1-scale .mas and .inp
'''

path_xdlsm = r'c:\xd2016\bin\xdlsm.exe'
scalefac_num = 12
max_cycles = 2
ext_to_save = ['hkl','fco','res','inp','fou','mas']

# polynomial to fit
def func(x,a,b):
    return 1+(a*x**2+b*x**3)

def save_copies(ext_to_save, num):
    for ext in ext_to_save:
        fname = f"xd.{ext}"
        if os.path.exists(fname):
            shutil.copy(fname, f"xd{num}.{ext}")

# do it
def main():
    # get the next available number
    num = f"{len(glob.glob('xd??.hkl')):>02}"
    
    # print an info header
    print(f"\n {'>':>^66}")
    print(f" >>> {f'Initital 1 Scalefactor Refinement':^58} >>>")
    # run xdlsm
    p = Popen(path_xdlsm)
    p.wait()
    
    # save copies
    save_copies(ext_to_save, num)
    
    # read fco data
    df = pd.DataFrame(np.loadtxt(f"xd.fco", skiprows=26, usecols=(0,1,2,3,4,6,7)),
                      columns=['h','k','l','F2c','F2o','stl','xdr'], dtype=np.float)
    
    # reject unused
    df.query("xdr == 0", inplace=True)
    df.drop(['xdr'], axis=1, inplace=True)
    
    # read hkl header and data
    with open(f"xd.hkl") as of:
        header = of.readline()
        df_hkl = pd.DataFrame(np.loadtxt(of, usecols=(0,1,2,3,4,5), comments='!'),
                              columns=['h','k','l','bn','hkl_F2o','hkl_F2s'],
                              dtype=np.float)
    
    # merge hkl and fco data (get sintl)
    df = pd.merge(df, df_hkl, on=['h','k','l'])
    
    # apply types
    df[['h','k','l','bn']] = df[['h','k','l','bn']].astype(np.int8)
    
    # bin it
    binned, bins = pd.cut(df['stl'], scalefac_num, retbins=True)
    # calculate the centers of the bins
    bins = (bins[1:]+bins[:-1])/2
    # batch number = bin number
    df['bn'] = df.groupby(binned)['bn'].ngroup() + 1
    
    # write new hkl
    with open(f"xd.hkl", 'w') as wf:
        wf.write(header)
        df.drop(['F2c','F2o','stl'], axis=1, inplace=True)
        df.to_string(wf, index=False, header=False,
                     columns=['h','k','l','bn','hkl_F2o','hkl_F2s'],
                     formatters=['{:4.0f}'.format,'{:3.0f}'.format,
                                 '{:3.0f}'.format,'{:1.0f}'.format,
                                 '{:12.3f}'.format,'{:12.3f}'.format])
    
    # edit xd input
    with open(f"xd.inp") as of:
        rf = of.read()
    # update the USAGE line
    rf = re.sub('(USAGE(?:\s+\d+){5})\s+(?:\d+)((?:\s+\d+){8})',
                f"\g<1>{scalefac_num:4}\g<2>", rf)
    # add the new scalefactor statring values
    rf = re.sub('^(\s+\d+\.\d+E(?:\-|\+)\d[1-9])$',
                ('\g<1>'*6+'\n')*(scalefac_num//6)+'\g<1>'*(scalefac_num%6),
                rf, flags=re.MULTILINE)
    # write new .inp
    with open(f"xd.inp", 'w') as wf:
        wf.write(rf)
    
    # edit xd master
    with open(f"xd.mas") as of:
        rf = of.read()
    # add the new number of scalefactors
    rf = re.sub('(SCALE\s+)(\d)+', '\g<1>'+'1'*scalefac_num, rf)
    # write new .mas
    with open(f"xd.mas", 'w') as wf:
        wf.write(rf)
    
    # print an info header
    print(f"\n {'>':>^66}")
    print(f" >>> {f'Initital {scalefac_num:>2} Scalefactors Refinement':^58} >>>")
    # run xdlsm
    p = Popen(path_xdlsm)
    p.wait()
    
    # get the next available number
    num = f"{len(glob.glob('xd??.hkl')):>02}"
    
    # save copies
    save_copies(ext_to_save, num)
    
    #############################
    ## iterate it (if needed!) ##
    #############################
    corr_fact = np.array([0.0,0.0])
    for _ in range(max_cycles):
        
        # read xd result
        with open(f"xd.res") as of:
            rf = of.read()
        
        # parse new scalefactors from .res
        sfacs = np.array(re.findall('(\d+\.\d+E(?:\-|\+)\d[1-9])', rf,
                         flags=re.MULTILINE)[-scalefac_num:]).astype(np.float)
        # we need to square the scalefactors (K)
        # Fo**2 = Fc**2 * K**2
        sfacs = sfacs**2
        # for some reason this scales it future problem!
        sfacs /= sfacs.max()
        
        # fit it
        popt, pcov = curve_fit(func, bins, sfacs, p0=(0.0,0.0))
        
        # update overall correction factor
        corr_fact += popt
        
        # read fco data
        df = pd.DataFrame(np.loadtxt(f"xd.fco", skiprows=26, usecols=(0,1,2,3,4,6,7)),
                          columns=['h','k','l','F2c','F2o','stl','xdr'], dtype=np.float)
        # reject unused
        df.query("xdr == 0", inplace=True)
        df.drop(['xdr'], axis=1, inplace=True)
        
        # read hkl header and data
        with open(f"xd.hkl") as of:
            header = of.readline()
            df_hkl = pd.DataFrame(np.loadtxt(of, usecols=(0,1,2,3,4,5), comments='!'),
                                  columns=['h','k','l','bn','hkl_F2o','hkl_F2s'],
                                  dtype=np.float)
        
        # merge hkl and fco data (get sin(t)/l)
        df = pd.merge(df, df_hkl, on=['h','k','l'])
        
        # apply types
        df[['h','k','l','bn']] = df[['h','k','l','bn']].astype(np.int8)
        
        # apply correction
        df['hkl_F2o'] = df['hkl_F2o'] / df['stl'].apply(func, args=(*popt,))
        
        # write new hkl
        with open(f"xd.hkl", 'w') as wf:
            wf.write(header)
            wf.write(f'!TDS CORRECTION FACTOR: a={popt[0]:.3f}, b={popt[1]:.3f}\n')
            df.drop(['F2c','F2o','stl'], axis=1, inplace=True)
            df.to_string(wf, index=False, header=False,
                         columns=['h','k','l','bn','hkl_F2o','hkl_F2s'],
                         formatters=['{:4.0f}'.format,'{:3.0f}'.format,
                                     '{:3.0f}'.format,'{:1.0f}'.format,
                                     '{:12f}'.format,'{:12f}'.format])
        
        # print an info header
        print(f"\n {'>':>^66}")
        print(f" >>> {f'Cycle {num} {popt}':^58} >>>")
        
        # run xdlsm
        p = Popen(path_xdlsm)
        p.wait()

        # get the next available number
        num = f"{len(glob.glob('xd??.hkl')):>02}"
        
        # save copies
        save_copies(ext_to_save, num)
    
    #################################
    ## do 1 scalefactor refinement ##
    ##       with the overall      ##
    ##      correction factor      ##
    #################################
    # read fco data
    df = pd.DataFrame(np.loadtxt(f"xd00.fco", skiprows=26, usecols=(0,1,2,3,4,6,7)),
                      columns=['h','k','l','F2c','F2o','stl','xdr'], dtype=np.float)
    # reject unused
    df.query("xdr == 0", inplace=True)
    df.drop(['xdr'], axis=1, inplace=True)
    
    # read hkl header and data
    with open(f"xd00.hkl") as of:
        header = of.readline()
        df_hkl = pd.DataFrame(np.loadtxt(of, usecols=(0,1,2,3,4,5), comments='!'),
                              columns=['h','k','l','bn','hkl_F2o','hkl_F2s'],
                              dtype=np.float)
    
    # merge hkl and fco data (get sin(t)/l)
    df = pd.merge(df, df_hkl, on=['h','k','l'])
    
    # apply types
    df[['h','k','l','bn']] = df[['h','k','l','bn']].astype(np.int8)
    
    # apply overall correction
    df['hkl_F2o'] = df['hkl_F2o'] / df['stl'].apply(func, args=(*corr_fact,))
    
    # write new hkl
    with open(f"xd.hkl", 'w') as wf:
        wf.write(header)
        wf.write(f'!TDS CORRECTION FACTOR: a={corr_fact[0]:.3f}, b={corr_fact[1]:.3f}\n')
        df.drop(['F2c','F2o','stl'], axis=1, inplace=True)
        df.to_string(wf, index=False, header=False,
                    columns=['h','k','l','bn','hkl_F2o','hkl_F2s'],
                    formatters=['{:4.0f}'.format,'{:3.0f}'.format,
                                '{:3.0f}'.format,'{:1.0f}'.format,
                                '{:12f}'.format,'{:12f}'.format])
    
    shutil.copy(f"xd00.mas", f"xd.mas")
    shutil.copy(f"xd00.inp", f"xd.inp")
    
    # print an info header
    print(f"\n {'>':>^66}")
    print(f" >>> {f'Final 1 Scalefactor Refinement':^58} >>>")
    print(f" >>> {f'{corr_fact}':^58} >>>")
    
    # run xdlsm
    p = Popen(path_xdlsm)
    p.wait()
        
if __name__ == '__main__':
    main()