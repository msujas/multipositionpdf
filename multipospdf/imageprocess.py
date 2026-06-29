from cryio.cbfimage import CbfImage
import numpy as np
from pyFAI.geometry import Geometry
from .common import getyz
from .functions import PoniYZ
import math, os
from fabio.edfimage import EdfImage

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

def correctImage(imagefile, poniyz, pfactor=0.85, outsubdir='polsaCorrected') -> np.ndarray:  
    y,z = getyz(imagefile)
    array = CbfImage(imagefile).array
    geo:Geometry = estponi(poniyz, y,z)
    correctedarray = array/(geo.solidAngleArray() * geo.polarization(factor=pfactor))
    if outsubdir:
        outdir = f'{os.path.dirname(imagefile)}/{outsubdir}'
        os.makedirs(outdir,exist_ok=True)
        outfile = f'{outdir}/{os.path.basename(imagefile)}'
        im = EdfImage(data=correctedarray)
        im.save(outfile)
    return correctedarray


