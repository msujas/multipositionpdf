from .common import *
from multipos_rustpy import integrate_rp
import time

def runintegration_rp(cbfdir, ponidir, tthmin:float, tthmax:float,tthbins:int, chimin:float, chimax:float, 
                      chibins:int, pfactor:float=0.85, maskfile=None, savecakes=False,outsubdir='cakes',
                      cakemaskfile=None, ponipattern = '*.poni', outponisubdir = 'ponis_rp'):
    '''
    a function for collecting poni and cbf files. Interpolating in Python (SciPy), then integrating and averaging with 
    a Rust extension. Requires installation of the multipos_rustpy library.
    cbfdir - directory of cbf images
    ponidir - directory for original ponis - i.e. small number to be interpolated
    pfactor - polarization factor (default 0.85)
    cakemaskfile - file for the cake mask - must be same shape as the integrated cake files (will ignore if not)
    savecakes - bool, save a cake for each image (default False, average cake always saved). Could be useful in some cases, but takes time
                and disk space
    outsubdir - subdirectory for integrated cake files (default - cakes)
    outponisubdir - subdirectory to save all the interpolated ponis (one for each image, default - ponis_rp)
    '''
    t0 = time.time()
    plist = getponilist(ponidir,ponipattern)
    mf = getIPlist(cbfdir,plist)
    mf.saveAllPonis(outponisubdir)
    
    newponidir = f'{cbfdir}/{outponisubdir}'
    print(f'ponis interpolated and saved to {newponidir}. Elapsed time {time.time()-t0:.2f} s')
    print('running integration with Rust extension')
    integrate_rp(cbfdir, newponidir, tthmin, tthmax,tthbins, chimin,chimax,chibins,pfactor, maskfile, savecakes,
                 outsubdir, cakemaskfile)
    print(f'finished. Total time elapsed: {time.time()-t0:.2f} s')