from .common import *
from multipos_rustpy import integrate_rp
import time

def runintegration_rp(cbfdir, ponidir, tthmin:float, tthmax:float,tthbins:int, chimin:float, chimax:float, 
                      chibins:int, pfactor:float=0.85, maskfile=None, savcakes=False,outsubdir='cakes',
                        cakemaskfile=None, ponipattern = '*.poni'):
    '''
    a function for collecting poni and cbf files. Interpolating in Python (SciPy), then integrating and averaging with 
    a Rust extension. Requires installation of the multipos_rustpy library
    '''
    t0 = time.time()
    plist = getponilist(ponidir,ponipattern)
    mf = getIPlist(cbfdir,plist)
    ponisubdir = 'ponis_rp'
    mf.saveAllPonis(ponisubdir)
    
    newponidir = f'{cbfdir}/{ponisubdir}'
    print(f'ponis interpolated and saved to {newponidir}. Elapsed time {time.time()-t0:.2f} s')
    print('running integration with Rust extension')
    integrate_rp(cbfdir, newponidir, tthmin, tthmax,tthbins, chimin,chimax,chibins,pfactor, maskfile, savcakes,
                 outsubdir, cakemaskfile)
    print(f'finished. Total time elapsed: {time.time()-t0:.2f} s')