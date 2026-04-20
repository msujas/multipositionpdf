import numpy as np
import pyFAI
from fabio.edfimage import EdfImage
from cryio.cbfimage import CbfImage
import os, fabio
from glob import glob
from pyFAI.integrator.azimuthal import AzimuthalIntegrator
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d, interpn, RegularGridInterpolator, LinearNDInterpolator
from pyFAI.geometry import Geometry

def bubbleHeader(file2d,array2d, tth, eta, y, e):
    xye = np.array([tth,y,e]).transpose().flatten()
    xyestring = ' '.join([str(i) for i in xye])
    header = {
    'Bubble_cake_version' : 3,
    'Bubble_cake' : f'{tth[0]} {tth[-1]} {eta[0]} {eta[-1]}',
    'Bubble_normalized': 1 ,
    'Bubble_pattern': xyestring
    }
    f = EdfImage(data = array2d[::-1,:], header = header)
    f.write(file2d)

class PoniData():
    def __init__(self, poni:AzimuthalIntegrator, ypos:float, zpos:float = None):
        if isinstance(poni,AzimuthalIntegrator):
            self.poni:AzimuthalIntegrator = poni
        elif isinstance(poni,str):
            self.poni = pyFAI.load(poni)
        else:
            raise ValueError('poni must be AzimuthalIntegrator or str type')
        self.ypos:float = ypos
        self.zpos:float = zpos
        self.poni1 = self.poni.poni1
        self.poni2 = self.poni.poni2
        self.rot1 = self.poni.rot1
        self.rot2 = self.poni.rot2
        self.rot3 = self.poni.rot3
        self.dist = self.poni.dist
        self.wavelength = self.poni.wavelength
        self.detector = self.poni.detector

class PoniList():
    def __init__(self, alist:list[PoniData] = None):
        self.__list__ : list[PoniData] = []
        if alist:
            self.__list__ = alist
            self.getValues()
    def __getitem__(self, key) -> PoniData :
        return self.__list__[key]
    def __contains__(self, item):
        return item in self.__list__
    def __iter__(self):
        return iter(self.__list__)
    def append(self,item):
        self.__list__.append(item)
        self.getValues()
    def getValues(self):
        self.ypositions = [p.ypos for p in self.__list__]
        self.zpositions = [p.zpos for p in self.__list__]
        self.poni1s = [p.poni1 for p in self.__list__]
        self.poni2s = [p.poni2 for p in self.__list__]
        self.rot1s = [p.rot1 for p in self.__list__]
        self.rot2s = [p.rot2 for p in self.__list__]
        self.rot3s = [p.rot3 for p in self.__list__]
        self.distances = [p.dist for p in self.__list__]
        self.wavelengths = [p.wavelength for p in self.__list__]
        self.detectors = [p.detector for p in self.__list__]
        self.interpolationFunctions()
        allsame = (np.array(self.zpositions) == self.zpositions[0]).all()
        self.poniinterpolation = '1d'
        if not None in self.zpositions and not allsame:
            self.poniinterpolation = '2d'
            self.interpolationFunctions2d()
    def __setitem__(self, key, value):
        self.__list__[key] = value
        self.getValues()
    def interpolationFunctions(self):
        self.poni1int = interp1d(self.ypositions, self.poni1s)
        self.poni2int = interp1d(self.ypositions,self.poni2s)
        self.distint = interp1d(self.ypositions,self.distances)
        self.rot1int =interp1d(self.ypositions,self.rot1s)
        self.rot2int = interp1d(self.ypositions, self.rot2s)       
    def interpolationFunctions2d(self):
        self.poni1int2d = LinearNDInterpolator(list(zip(self.ypositions,self.zpositions)), self.poni1s)
        self.poni2int2d = LinearNDInterpolator(list(zip(self.ypositions,self.zpositions)),self.poni2s)
        self.distint2d = LinearNDInterpolator(list(zip(self.ypositions,self.zpositions)),self.distances)
        self.rot1int2d =LinearNDInterpolator(list(zip(self.ypositions,self.zpositions)),self.rot1s)
        self.rot2int2d = LinearNDInterpolator(list(zip(self.ypositions,self.zpositions)), self.rot2s)
        
class FilePoni():
    def __init__(self, fname:str, ypos:float, ponilist:PoniList, zpos:float = None):
        self.fname = fname
        self.ypos = ypos
        self.zpos = zpos
        self.ponilist = ponilist
        self.aiexample:AzimuthalIntegrator = ponilist[0].poni
        self.detector = self.aiexample.detector
        self.wavelength = self.aiexample.wavelength
        self.array = CbfImage(fname).array
        
        self.integrated = False
        self.poni1s = self.ponilist.poni1s 
        self.poni2s = self.ponilist.poni2s 
        self.dists = self.ponilist.distances 
        self.rot1s = self.ponilist.rot1s 
        self.rot2s = self.ponilist.rot2s 
        self.ypositions = self.ponilist.ypositions 
        self.zpositions = self.ponilist.zpositions 

    def interpolatePoni(self, basemask:str = None):
        self.poni1 = self.ponilist.poni1int(self.ypos)
        self.poni2 = self.ponilist.poni2int(self.ypos)
        self.dist = self.ponilist.distint(self.ypos)
        self.rot1 = self.ponilist.rot1int(self.ypos)
        self.rot2 = self.ponilist.rot2int(self.ypos)
        self.rot3=0
        self.integrate(basemask)
    def interpolatePoni2D(self,basemask:str=None):

        self.poni1 = self.ponilist.poni1int2d(self.ypos,self.zpos)
        self.poni2 = self.ponilist.poni2int2d(self.ypos,self.zpos)
        self.dist = self.ponilist.distint2d(self.ypos,self.zpos)
        self.rot1 = self.ponilist.rot1int2d(self.ypos,self.zpos)
        self.rot2 = self.ponilist.rot2int2d(self.ypos,self.zpos)
        self.rot3 = 0
        self.integrate(basemask)

    def integrate(self,basemask:str = None):
        self.ai = AzimuthalIntegrator(dist=self.dist, poni1=self.poni1, poni2=self.poni2, rot1=self.rot1,
                                      rot2=self.rot2, rot3=self.rot3, detector=self.detector,wavelength=self.wavelength)
        self.geometry = Geometry(dist=self.dist, poni1=self.poni1, poni2=self.poni2, rot1=self.rot1,
                                      rot2=self.rot2, rot3=self.rot3, detector=self.detector,wavelength=self.wavelength)
        
        self.ttharray:np.ndarray = self.geometry.twoThetaArray()
        self.chiarray:np.ndarray = self.geometry.chiArray()
        self.pol:np.ndarray = self.geometry.polarization(factor=0.99)
        self.sa :np.ndarray = self.geometry.solidAngleArray()
        self.sinchi2:np.ndarray = np.sin(self.chiarray)**2
        self.arrayCorrected = self.array/(self.pol*self.sa)
        mask = np.where(self.array < 0,1,0)
        if basemask:
            mask = fabio.open(basemask).data
        dirname = os.path.dirname(self.fname)
        outdir = f'{dirname}/cake'
        outdir1d = f'{dirname}/xye'
        basename = os.path.basename(self.fname).replace('.cbf','')
        outfile = f'{outdir}/{basename}.edf'
        outfile1d = f'{outdir1d}/{basename}.xye'
        os.makedirs(outdir,exist_ok=True)
        os.makedirs(outdir1d,exist_ok=True)
        
        self.x,self.y,self.e = self.ai.integrate1d(self.array,5000,correctSolidAngle=True, polarization_factor=0.99,method='bbox',
                                            unit='2th_deg', mask=mask, error_model='poisson')
        self.result2d = self.ai.integrate2d(self.array,5000,360,correctSolidAngle=True, polarization_factor=0.99,method='bbox',
                                            unit='2th_deg', mask=mask, error_model='poisson')
        bubbleHeader(outfile, *self.result2d[:3],self.y,self.e)
        np.savetxt(outfile1d, np.array([self.x,self.y,self.e]).transpose(), fmt = '%.6f')
        self.array2d = self.result2d[0]
        self.tth = self.result2d[1]
        self.chi = self.result2d[2]
        self.integrated= True
        

    def arrayBins(self,tthbins:int, chibins:int, tthmax:float):
        tthmax = tthmax*np.pi/180
        self.tthbinarray = (self.ttharray*tthbins/tthmax).astype(int)
        self.chibinarray = (self.sinchi2*chibins).astype(int)
        self.combinedbinarray = self.tthbinarray+1j*self.chibinarray
        return self.combinedbinarray

class MultiFile():
    def __init__(self,alist:list[FilePoni]):
        self.list = alist
    def integrateAll(self):
        for file in self.list:
            file.interpolatePoni()
    def average1d(self,x0,xend,npoints, basemask:str=None):
        self.x = np.linspace(x0,xend,npoints)
        avarray = np.empty(shape=(len(self.x), len(self.list)))
        
        poniinterpolation = self.list[0].ponilist.poniinterpolation

        if not self.list[0].integrated:
            print(f'interpolating ponis in {poniinterpolation} and integrating files')
            for file in self.list:
                match poniinterpolation:
                    case '1d': file.interpolatePoni(basemask=basemask)
                    case '2d': file.interpolatePoni2D(basemask=basemask)
                    case _: raise ValueError('poniinterpolation must be "1d" or "2d"')
        for i,file in enumerate(self.list):
            gridfunc = interp1d(file.x,file.y)
            avarray[:,i] = gridfunc(self.x)
        self.yav = avarray.mean(axis=1)
        return self.x,self.yav
    def plotAll1d(self):
        plt.figure()
        for file in self.list:
            label = os.path.basename(file.fname)
            plt.plot(file.x,file.y,label = label)
        plt.plot(self.x,self.yav,label = 'average', linewidth = 3)
        plt.legend()
        plt.xlabel('2theta (°)')
        plt.ylabel('intensity')
        plt.show()
    def regrid2d(self,tth0, tthend, tthpoints, chi0=-178, chiend=178, chipoints=354, nstdevs = 3, medianFilter = 4, 
                 outdir = None, cakemask:str = None):
        tthgrid = np.linspace(tth0,tthend, tthpoints)
        chigrid = np.linspace(chi0,chiend,chipoints)
        mesh = np.empty(shape=(len(chigrid),len(tthgrid),2))
        for i in range(len(chigrid)):
            for j in range(len(tthgrid)):
                mesh[i,j] = [chigrid[i],tthgrid[j]]
        ndlist = np.empty(shape = (*mesh.shape[:2],len(self.list)))
        for i,item in enumerate(self.list):
            interparray = np.where(item.array2d<=0, np.nan, item.array2d)
            rgf = RegularGridInterpolator((item.chi,item.tth), interparray)
            newdata = rgf(mesh)
            ndlist[:,:,i] = newdata
        stdev = np.nanstd(ndlist,axis=2)
        median = np.nanmedian(ndlist,axis = 2)
        masks = np.zeros(shape = ndlist.shape)
        if cakemask:
            cakemask = fabio.open(cakemask).data
        else:
            cakemask = 0
            
        for i in range(ndlist.shape[2]):
            masks[:,:,i] = np.where((ndlist[:,:,i]< median-nstdevs*stdev)|(ndlist[:,:,i]> median+nstdevs*stdev) | (ndlist[:,:,i]<=0) |
                                    (ndlist[:,:,i] > medianFilter*median), 1, cakemask )
        maskeddata = np.where(masks == 1, np.nan, ndlist)
        y2 = 0
        for i in range(maskeddata.shape[2]):
            y2 += np.nanmean(maskeddata[:,:,i],axis = 0)
        self.y2 = y2/maskeddata.shape[2]
        self.avmasked = np.nanmean(maskeddata,axis= 2)
        self.x = tthgrid
        self.y = np.nanmean(self.avmasked,axis = 0)
        e = self.y**0.5
        if outdir:
            os.makedirs(outdir,exist_ok=True)
            bubbleHeader(f'{outdir}/av2d.edf', np.where(np.isnan(self.avmasked),0,self.avmasked), tthgrid, chigrid,self.y,e)
            np.savetxt(f'{outdir}/av2d.xy',np.array([self.x,self.y2]).transpose(),fmt='%.6f')
        return self.avmasked
    def saveEDF_noheader(self,dirname):
        '''
        so files can be read with silx/pyfai to make masks
        '''
        im = EdfImage(self.avmasked)
        im.save(f'{dirname}/av2d_noheader.edf')
