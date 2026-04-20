import numpy as np
import pyFAI
from fabio.edfimage import EdfImage
from cryio.cbfimage import CbfImage
import os
from glob import glob
import matplotlib.pyplot as plt
from multipospdf import *



ponidir = r'D:\beamlineData\April2026\multipositions\Si'
direc = r'D:\beamlineData\April2026\multipositions\C60'
maskfile = r'D:\beamlineData\April2026\multipositions/baseMask_April2026.edf' # mask file to apply to all cbf images. None to not use
cakemask = r'D:\beamlineData\April2026\multipositions/cakemask.edf' #mask file to apply to final merged cake. Must generate initial cake first. None to not use
tthmin = 0.75
tthmax = 58
datadir = f'{direc}/multiPDF0'

def readfiles(direc):
    files = glob(f'{direc}/*.cbf')
    for i,file in enumerate(files):
        im = CbfImage(file)
        array = im.array
        header = im.header
        if i==0:
            dataset = np.empty(shape=(*array.shape,len(files)))
        dataset[:,:,i] = array

    return dataset


def readponis(direc, ypositions = None) -> PoniList:
    ponis = glob(f'{direc}/*MD.poni')
    ponis.sort()
    ponilist = PoniList()
    for p in ponis:
        basename= os.path.basename(p)
        index = int(basename.split('_')[0].replace('pos',''))
        if ypositions:
            ypos = ypositions[index-1]
        else:
            ypos = index
        ponilist.append(PoniData(pyFAI.load(p),ypos, zpos=0))
    return ponilist

def plotponi(direc, sifiles:dict):
    ypositions = list(sifiles.keys())
    ponidct = readponis(direc,ypositions)
    y = []
    p1s = []
    p2s = []
    r1s = []
    r2s = []
    dists = []
    for p in ponidct:
        print(p)
        poni1 = ponidct[p].poni1
        poni2 = ponidct[p].poni2
        r1 = ponidct[p].rot1
        r2 = ponidct[p].rot2
        dist = ponidct[p].dist
        y.append(ponidct[p].ypos)
        p1s.append(poni1)
        p2s.append(poni2)
        r1s.append(r1)
        r2s.append(r2)
        dists.append(dist)
    plt.figure()
    plt.plot(y,p1s, 'o-', label = 'poni1')
    plt.plot(y,p2s,'o-', label = 'poni2')
    plt.plot(y,r1s,'o-',label= 'rot1')
    plt.plot(y,r2s,'o-',label = 'rot2')
    plt.plot(y,dists,'o-',label = 'distance')
    plt.legend()
    plt.show()

def readData(direc):
    files = glob(f'{direc}/*.cbf')
    files.sort()
    ypositions = [float(os.path.basename(f).split('_')[1].replace('dty','')) for f in files]
    filedct = {y:f for y,f in zip(ypositions,files)}
    return filedct
        

datadct = readData(datadir)
ypositions = list(datadct.keys())
files = list(datadct.values())
ponilist = readponis(ponidir, ypositions)
#plotponi(direc, sidct)

flist = [FilePoni(f,y,ponilist) for f,y in zip(files,ypositions)]
multi = MultiFile(flist)

def plotFilePonis(fileponis:list[FilePoni]):
    ys = [f.ypos for f in fileponis]
    p1s = [f.poni1 for f in fileponis]
    p2s = [f.poni2 for f in fileponis]
    dists = [f.dist for f in fileponis]
    r1s = [f.rot1 for f in fileponis]
    r2s = [f.rot2 for f in fileponis]
    plt.figure()
    plt.plot(ys,p1s,'o-', label = 'poni1')
    plt.plot(ys,p2s,'o-',label = 'poni2')
    plt.plot(ys,dists,'o-',label = 'distances')
    plt.plot(ys,r1s,'o-', label = 'rot1')
    plt.plot(ys,r2s,'o-',label = 'rot2')
    plt.legend()
    plt.show()

def generateMaps(fplist:list[FilePoni], destdir:str):
    os.makedirs(destdir,exist_ok=True)
    os.makedirs(f'{destdir}/tth',exist_ok=True)
    os.makedirs(f'{destdir}/pol',exist_ok=True)
    os.makedirs(f'{destdir}/sa',exist_ok=True)
    os.makedirs(f'{destdir}/chi',exist_ok=True)
    for f in fplist:
        basename = os.path.splitext(os.path.basename(f.fname))[0]
        outfile = f'{destdir}/{basename}'
        tth = f.geometry.twoThetaArray()
        pol = f.geometry.polarization(factor=0.99)
        sa = f.geometry.solidAngleArray()
        chi = f.geometry.chiArray()
        im = EdfImage(tth)
        im.save(f'{destdir}/tth/{basename}.edf')
        im = EdfImage(pol)
        im.save(f'{destdir}/pol/{basename}.edf')
        im = EdfImage(sa)
        im.save(f'{destdir}/sa/{basename}.edf')
        im = EdfImage(chi)
        im.save(f'{destdir}/chi/{basename}.edf')
        
    
def getbinarrays(flist:list[FilePoni], tthbins, chibins, tthmax):
    for f in flist:
        f.arrayBins(tthbins, chibins, tthmax)


#ponilist.plot1d()

print('interpolating ponis, integrating 1d')
x,y = multi.average1d(tthmin,tthmax,5000,basemask=maskfile)
np.savetxt(f'{direc}/av1d.xy',np.array([x,y]).transpose(),fmt='%.6f')
#multi.plotAll1d()
print('calculating cakes')
av2d = multi.regrid2d(tthmin,tthmax,5000, outdir=f'{datadir}/cake', cakemask=cakemask)

#multi.saveEDF_noheader(direc) #use this to generate a cake file with no header modification, so it can be read with silx/pyFAI for making mask
#plt.imshow(av2d,aspect='auto')
#plt.show()