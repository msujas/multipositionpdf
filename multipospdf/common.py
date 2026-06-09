from .functions import MultiFile, ImagePoni, PoniYZ, PoniList
import os
from glob import glob

def getyz(file): #assuming files are in format: x_dty125.40_dtz135.00_...
    '''
    function for getting y and z values from file name in common file structure.
    Assumes file names are in format ..._dty124.52_dtz260.10_...
    '''
    basefile = os.path.splitext(os.path.basename(file))[0]
    filesplit = basefile.split('_')

    ypart = [f for f in filesplit if 'dty' in f][0]
    zpart = [f for f in filesplit if 'dtz' in f][0]
    ypos = float(ypart.replace('dty',''))
    zpos = float(zpart.replace('dtz',''))
    return ypos,zpos


def getponilist(ponidir, ponipattern = '*.poni'):
    files = glob(f'{ponidir}/{ponipattern}')
    ponilist= []
    for file in files:
        try:
            y,z = getyz(file)
        except IndexError:
            print(f"couldn't extract y and z positions for {file}. Skipping")
            continue
        ponilist.append(PoniYZ(file, y,z))
    return PoniList(ponilist)

def getIPlist(cbfdir,ponilist:PoniList, maskfile=None, pfactor=0.85, wavelength = None):
    files = glob(f'{cbfdir}/*.cbf')
    iplist = []
    for file in files:
        try:
            y,z = getyz(file)
        except:
            print(f"couldn't extract y and z positions for {file}. Skipping")
            continue
        iplist.append(ImagePoni(file, y,z,ponilist, maskfile, pfactor, wavelength=wavelength))
    return MultiFile(iplist)
