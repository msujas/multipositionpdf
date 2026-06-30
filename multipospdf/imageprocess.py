from cryio.cbfimage import CbfImage
import numpy as np
from pyFAI.geometry import Geometry
from .common import getyz
from .functions import PoniYZ
import math, os, fabio
from fabio.edfimage import EdfImage
from glob import glob


def gainCorrection(avim,gainArray):
    avimGain = avim/gainArray
    avimGain = np.where(gainArray <0, -1, avimGain)
    return avimGain

def estponi(poniyz:PoniYZ, posy,posz)->Geometry:
    '''
    function for estimating ponis for an image
    '''
    poni = poniyz.poni

    geo = Geometry()
    geo.load(poni)
    ydiff = posy - poniyz.ypos
    newponi2 = geo.poni2 + ydiff/1000
    rot1 = geo.rot1
    rot2 = geo.rot2
    dist = geo.dist
    zdiff = posz - poniyz.zpos
    newponi1 = geo.poni1 + zdiff*math.cos(rot1)/1000
    newdist = dist + zdiff*math.sin(rot1)/1000

    newgeo = Geometry(dist=newdist, poni1=newponi1,poni2=newponi2, rot1=rot1,rot2=rot2, 
                                  detector=geo.detector, wavelength=geo.wavelength)
    return newgeo

def correctImage(imagefile, poniyz:PoniYZ, pfactor=0.85, outsubdir='polsaCorrected') -> np.ndarray:  
    y,z = getyz(imagefile)
    array = CbfImage(imagefile).array
    geo:Geometry = estponi(poniyz, y,z)
    correctedarray = array/(geo.solidAngleArray() * geo.polarization(factor=pfactor))
    if outsubdir:
        outdir = f'{os.path.dirname(imagefile)}/{outsubdir}'
        os.makedirs(outdir,exist_ok=True)
        outfile = f'{outdir}/{os.path.basename(imagefile)}'.replace('.cbf','.edf')
        im = EdfImage(data=correctedarray)
        im.save(outfile)
    return correctedarray


def correctdir(dirname, poniyz:PoniYZ, pfactor=0.85, outsubdir='polsaCorrected' ):
    files = glob(f'{dirname}/*.cbf')
    for file in files:
        correctImage(file, poniyz, pfactor=pfactor, outsubdir=outsubdir)


def gaincorrectdir(dirname,gainfile, ext = 'cbf', outsubdir = 'GC'):
    gainarray = fabio.open(gainfile).data
    files = glob(f'{dirname}/*.{ext}')
    newdir = f'{dirname}/{outsubdir}'
    os.makedirs(newdir,exist_ok=True)
    for file in files:
        array = CbfImage(file).array
        arraygc = gainCorrection(array,gainarray)
        basename = os.path.splitext(os.path.basename(file))[0]
        newfile = f'{newdir}/{basename}.edf'
        im = EdfImage(data=arraygc)
        im.save(newfile)

def correctdir_polsagain(dirname, gainfile, poniyz:PoniYZ, pfactor=0.85, outsubdir='polsaGC'):
    files = glob(f'{dirname}/*.cbf')
    newdir = f'{dirname}/{outsubdir}'
    os.makedirs(newdir,exist_ok=True)
    gainarray = fabio.open(gainfile).data
    for file in files:
        arraycorrected = correctImage(file,poniyz, pfactor=pfactor, outsubdir=None)
        arraycorrected = gainCorrection(arraycorrected, gainarray)
        newfile = f'{newdir}/{os.path.basename(file)}'.replace('.cbf','.edf')
        im = EdfImage(arraycorrected)
        im.save(newfile)

