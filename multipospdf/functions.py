import numpy as np
import pyFAI
from fabio.edfimage import EdfImage
from cryio.cbfimage import CbfImage, CbfHeader
import os, fabio
from glob import glob
from pyFAI.integrator.azimuthal import AzimuthalIntegrator
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d, LinearNDInterpolator
from pyFAI.geometry import Geometry
from typing import Literal
try:
    from fluoCorrectionPilatus import fluoSub_integrated_base, getpolcakebase, FluosubCake
except ImportError:
    print("couldn't import fluo correction library. Fluo correction not available. Try installing it if needed")



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

def gainCorrection(avim,gainArray):
    avimGain = avim/gainArray
    avimGain = np.where(gainArray <0, -1, avimGain)
    return avimGain

class PoniYZ():
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
    def __init__(self, alist:list[PoniYZ] = None):
        self.__list__ : list[PoniYZ] = []
        if not all(isinstance(a, PoniYZ) for a in alist):
            raise ValueError('all input list values must be PoniYZ type')
        if alist:
            self.__list__ = alist
            self.getValues()
    def __getitem__(self, key) -> PoniYZ :
        return self.__list__[key]
    def __contains__(self, item):
        return item in self.__list__
    def __iter__(self):
        return iter(self.__list__)
    def append(self,item):
        if not isinstance(item,PoniYZ):
            raise ValueError('can only append PoniYZ values')
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
        self.interpolationDimension = 1
        if not None in self.zpositions and not allsame:
            self.interpolationDimension = 2
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
    def plot1d(self):
        plt.figure(dpi = 150)
        plt.plot(self.ypositions, self.poni1s,'o-',label = 'poni1')
        plt.plot(self.ypositions,self.poni2s,'o-',label = 'poni2')
        plt.plot(self.ypositions,self.rot1s,'o-',label = 'rot1')
        plt.plot(self.ypositions,self.rot2s,'o-',label = 'rot2')
        plt.plot(self.ypositions,self.rot3s,'o-',label = 'rot3')
        plt.plot(self.ypositions,self.distances,'o-',label = 'distance')
        plt.legend()
        plt.xlabel('y-position')
        plt.show()
    def plot2d(self):
        ygrid = np.linspace(min(self.ypositions),max(self.ypositions), 20)
        zgrid = np.linspace(min(self.zpositions),max(self.zpositions),20)
        y,z = np.meshgrid(ygrid,zgrid)
        poni1s = self.poni1int2d(y,z)
        poni2s = self.poni2int2d(y,z)
        rot1s = self.rot1int2d(y,z)
        rot2s = self.rot2int2d(y,z)
        distances = self.distint2d(y,z)
        markersize=  5
        plotdct = {'poni1':poni1s,
                   'poni2':poni2s,
                   'rot1':rot1s,
                   'rot2':rot2s,
                   'distances':distances}
        fig,ax = plt.subplots(2,3,dpi=150,figsize=(6.4*1.5,4.8))
        for i,item in enumerate(plotdct):
            xplot = i%3
            yplot = i//3
            im = ax[yplot,xplot].pcolormesh(y,z,plotdct[item],shading='auto')
            ax[yplot,xplot].scatter(self.ypositions,self.zpositions,c='red',s=markersize)
            fig.colorbar(im, ax=ax[yplot,xplot])
            ax[yplot,xplot].set_title(item)
        lasti = i
        for i in range(lasti+1,6):
            yplot = (i)//3
            xplot = (i)%3
            fig.delaxes(ax[yplot,xplot])
            plt.tight_layout()
        plt.show()

        
class ImagePoni():
    def __init__(self, fname:str, ypos:float,  zpos:float, ponilist:PoniList, maskfile:str = None, pfactor=0.85):
        self.fname = fname
        self.dirname = os.path.dirname(self.fname)
        self.basename = os.path.splitext(os.path.basename(self.fname))[0]
        header = CbfHeader(self.fname).header
        self.flux = header['Flux']
        exptime = header['Exposure_time']
        self.fluxrate = self.flux/exptime
        self.ypos = ypos
        self.zpos = zpos
        self.ponilist = ponilist
        self.aiexample:AzimuthalIntegrator = ponilist[0].poni
        self.detector = self.aiexample.detector
        self.wavelength = self.aiexample.wavelength
        self.array = CbfImage(fname).array
        self.pfactor = pfactor

        self.integrated = False
        self.mapscalculated = False
        self.poni1s = self.ponilist.poni1s 
        self.poni2s = self.ponilist.poni2s 
        self.dists = self.ponilist.distances 
        self.rot1s = self.ponilist.rot1s 
        self.rot2s = self.ponilist.rot2s 
        self.ypositions = self.ponilist.ypositions 
        self.zpositions = self.ponilist.zpositions 
        self.maskfile = maskfile
        self.solidanglescale = 10**6
        match self.ponilist.interpolationDimension:
            case 1: self.interpolatePoni()
            case 2: self.interpolatePoni2D()
        self.geometry = Geometry(dist=self.dist, poni1=self.poni1, poni2=self.poni2, rot1=self.rot1,
                                      rot2=self.rot2, rot3=self.rot3, detector=self.detector,wavelength=self.wavelength)
        self.ai = AzimuthalIntegrator(dist=self.dist, poni1=self.poni1, poni2=self.poni2, rot1=self.rot1,
                                      rot2=self.rot2, rot3=self.rot3, detector=self.detector,wavelength=self.wavelength)
    def interpolatePoni(self):
        self.poni1 = self.ponilist.poni1int(self.ypos)
        self.poni2 = self.ponilist.poni2int(self.ypos)
        self.dist = self.ponilist.distint(self.ypos)
        self.rot1 = self.ponilist.rot1int(self.ypos)
        self.rot2 = self.ponilist.rot2int(self.ypos)
        self.rot3=0

    def interpolatePoni2D(self):
        self.poni1 = self.ponilist.poni1int2d(self.ypos,self.zpos)
        self.poni2 = self.ponilist.poni2int2d(self.ypos,self.zpos)
        self.dist = self.ponilist.distint2d(self.ypos,self.zpos)
        self.rot1 = self.ponilist.rot1int2d(self.ypos,self.zpos)
        self.rot2 = self.ponilist.rot2int2d(self.ypos,self.zpos)
        self.rot3 = 0

    def integrate(self, tthmin, tthmax,tthbins=5000,  chimin = -178, chimax=178, chibins = 354, gainfile = None, xyedir = None, 
                  cakedir = None,  scale = 10**5, unit:Literal["2th_deg", "q_A^-1"] = '2th_deg'):
        mask = np.where(self.array < 0,1,0)
        if self.maskfile:
            mask = fabio.open(self.maskfile).data

        correctedarray = self.array.copy()
        if gainfile:
            gainmap = fabio.open(gainfile).data
            correctedarray = gainCorrection(correctedarray,gainmap)

        correctSolidAngle = False

        self.absSolidAnglemap = self.geometry.solidAngleArray(absolute=True) #have to correct with absolute solid angle, as it can change image to image
        correctedarray = (correctedarray/(self.absSolidAnglemap*self.solidanglescale))
        correctedarray = correctedarray*scale/self.flux

        self.x,self.y,self.e = self.ai.integrate1d(correctedarray,tthbins,correctSolidAngle=correctSolidAngle, polarization_factor=self.pfactor,
                                                   method='bbox',unit=unit, mask=mask, error_model='poisson', radial_range=(tthmin,tthmax))
        self.result2d = self.ai.integrate2d(correctedarray,tthbins,correctSolidAngle=correctSolidAngle, polarization_factor=self.pfactor,
                                            method='bbox',unit=unit, mask=mask, error_model='poisson',radial_range=(tthmin,tthmax), 
                                            azimuth_range=(chimin, chimax), npt_azim=chibins)

        self.array2d = self.result2d[0]
        self.tth = self.result2d[1]
        self.chi = self.result2d[2]
        self.integrated= True
        if cakedir:
            self.saveCake(cakedir=cakedir)
        if xyedir:
            self.save1d(xyedir=xyedir)
    def fluosub(self, fluoK):
        polcake = getpolcakebase(self.tth, self.chi, self.pfactor)
        self.cake_fluosub = fluoSub_integrated_base(self.array2d, polcake, fluoK)
    def saveCake(self,cakedir = 'cake'):
        outdir = f'{self.dirname}/{cakedir}'
        os.makedirs(outdir,exist_ok=True)
        outfile = f'{outdir}/{self.basename}.edf'
        bubbleHeader(outfile,self.array2d,self.tth,self.chi,self.y,self.e)
    def save1d(self,xyedir = 'xye'):
        outdir = f'{self.dirname}/{xyedir}'
        outfile = f'{self.dirname}/{self.basename}.xye'
        os.makedirs(outdir,exist_ok=True)
        np.savetxt(outfile, np.array([self.x,self.y,self.e]).transpose(), fmt = '%.6f')

    def getMaps(self):

        self.tthmap:np.ndarray = self.geometry.twoThetaArray()*180/np.pi
        self.chimap:np.ndarray = self.geometry.chiArray()
        self.polmap:np.ndarray = self.geometry.polarization(factor=self.pfactor)
        self.samap :np.ndarray = self.geometry.solidAngleArray()
        dety,detx,_ = self.geometry.detector.calc_cartesian_positions()
        ponidist = ((dety-self.poni1)**2 + (detx - self.poni2)**2)**0.5
        self.sampledistmap = (ponidist**2 + self.dist**2)**0.5
        self.sinchi2:np.ndarray = np.sin(self.chimap)**2
        
        self.absSolidAnglemap = self.geometry.solidAngleArray(absolute=True)
        self.mapscalculated = True
        self.arrayCorrected = self.array/(self.polmap*self.absSolidAnglemap*self.solidanglescale)
        

    def saveMaps(self,dirname):
        if not self.mapscalculated:
            self.getMaps(self.pfactor)
        basename = os.path.splitext(os.path.basename(self.fname))[0]
        basename += '.edf'
        os.makedirs(f'{dirname}/tth',exist_ok=True)
        os.makedirs(f'{dirname}/sa',exist_ok=True)
        os.makedirs(f'{dirname}/pol',exist_ok=True)
        os.makedirs(f'{dirname}/chi',exist_ok=True)
        os.makedirs(f'{dirname}/sampleDist',exist_ok=True)
        im = EdfImage(self.tthmap)
        im.save(f'{dirname}/tth/{basename}')
        im = EdfImage(self.chimap)
        im.save(f'{dirname}/chi/{basename}')
        im=EdfImage(self.polmap)
        im.save(f'{dirname}/pol/{basename}')
        im=EdfImage(self.samap)
        im.save(f'{dirname}/sa/{basename}')
        im = EdfImage(self.sampledistmap)
        im.save(f'{dirname}/sampleDist/{basename}')


    def arrayBins(self,tthbins:int, chibins:int, tthmax:float):
        tthmax = tthmax*np.pi/180
        self.tthbinarray = (self.tthmap*tthbins/tthmax).astype(int)
        self.chibinarray = (self.sinchi2*chibins).astype(int)
        self.combinedbinarray = self.tthbinarray+1j*self.chibinarray
        return self.combinedbinarray
    def saveponi(self, subdir):
        dirname = f'{self.dirname}/{subdir}'
        filename = f'{dirname}/{self.basename}.poni'
        os.makedirs(dirname, exist_ok=True)
        self.ai.save(filename)
        self.ponifile = filename

class MultiFile():
    def __init__(self,alist:list[ImagePoni]):
        self.list:list[ImagePoni] = []
        for i in alist:
            if not isinstance(i, ImagePoni):
                raise ValueError('input list must contain only ImagePoni values')
            if i.fluxrate < 100:
                print(f'image {i.fname} did not have enough flux, ignoring')
            elif not i.poni1:
                print(f'no valid poni for {i.fname}')
            else:
                self.list.append(i)
                
    def saveMaps(self,dirname):
        for f in self.list:
            f.saveMaps(dirname)
        
    def average1d(self,tthmin,tthmax,tthbins, chimin=-178, chimax=178, chibins=354,   outsubdir= 'xye', 
                  fname='', unit:Literal["2th_deg","q_A^-1"] = '2th_deg', **kwargs):
        '''
        run the interpolations (if not done already) and regrid and average all 1d patterns
        x0 - 2theta0
        xend - 2theta end
        npoints - number of points in regrid
        outsubdir - subdirectory for averaged 1d file
        kwargs - kwargs for interpolation
        '''
        print(f'integrating images')
        for file in self.list:
            print(file.fname)
            file.integrate(tthmin,tthmax,tthbins=tthbins, chimin=chimin, chimax=chimax, chibins=chibins, unit = unit,**kwargs)
        
        self.pfactor = file.pfactor
        self.x = self.list[0].x
        avarray = np.empty(shape=(len(self.x), len(self.list)))
        for i,file in enumerate(self.list):
            avarray[:,i] = np.where(file.y <=0, np.nan, file.y)
        self.yav = np.nanmean(avarray,axis=1)
        if outsubdir:
            outdir = f'{os.path.dirname(self.list[0].fname)}/{outsubdir}'
            os.makedirs(outdir,exist_ok=True)
            outfile = f'{outdir}/{fname}av1d.xy'
            np.savetxt(outfile,np.array([self.x,self.yav]).transpose(),fmt='%.6f')
        return self.x,self.yav
    def average2d(self,outsubdir='cake',fname='', nstdevs=3, medianfilter=4, cakemask=None, fluoK = 0):
        print('averaging cakes')
        if cakemask:
            cakemask = fabio.open(cakemask).data
        else:
            cakemask=0
        if fluoK:
            outsubdir += 'fluosub'
        for i,item in enumerate(self.list):
            if i == 0:
                allarrays = np.empty(shape=(*item.array2d.shape,len(self.list)))
            if fluoK:
                item.fluosub(fluoK )
                allarrays[:,:,i] = np.where(item.cake_fluosub <= 0 ,np.nan, item.cake_fluosub)
            else:
                allarrays[:,:,i] = np.where(item.array2d <= 0 ,np.nan, item.array2d)
        masks = self._getmasks(allarrays, cakemask, nstdevs=nstdevs,medianfilter=medianfilter)
        allarrays = np.where(masks >0, np.nan, allarrays)
        self.avarray = np.nanmean(allarrays,axis=2)
        self.ycake = np.nanmean(self.avarray,axis=0)
        self.avarray = np.where(np.isnan(self.avarray),0, self.avarray)
        self.tth = self.list[0].tth
        self.chi = self.list[0].chi
        self.ycake2 = np.empty(shape=(len(self.tth), len(self.list)))
        for i in range(allarrays.shape[2]):
            self.ycake2[:,i] = np.nanmean(allarrays[:,:,i],axis= 0)
        self.ycake2 = np.nanmean(self.ycake2,axis=1)
        if outsubdir:
            basedir = os.path.dirname(self.list[0].fname)
            outdir = f'{basedir}/{outsubdir}'
            os.makedirs(outdir,exist_ok=True)
            bubbleHeader(f'{outdir}/{fname}av2d.edf', self.avarray , self.tth, self.chi, self.ycake2,self.ycake**0.5)
            np.savetxt(f'{outdir}/{fname}av2d.xy',np.array([self.tth,self.ycake2]).transpose(),fmt = '%.6f')
            np.savetxt(f'{outdir}/{fname}av2d_2.xy',np.array([self.tth,self.ycake]).transpose(),fmt = '%.6f')
    def fluosubav(self,k0):
        fc = FluosubCake(pfactor=self.pfactor)
        self.polcake = getpolcakebase(self.tth, self.chi, self.pfactor)
        self.avcakefluosub = fc.optimise_fluoIntegrated(self.avarray, self.polcake, k0)
        self.fluoK = fc.kopt
    def average2d_optimise_rerun(self, k0, outsubdir, fname,saveindividual=False, **kwargs):
        '''
        slow - runs the integration twice, once to get the average, then optimise fluo correction, 
        then rerun with fluo subtraction
        '''
        self.average2d(outsubdir=outsubdir, fname=fname,**kwargs)
        self.fluosubav(k0)
        basedir = os.path.dirname(self.list[0].fname)
        outdir = f'{basedir}/{outsubdir}fluosub'
        os.makedirs(outdir,exist_ok=True)
        outfile = f'{outdir}/{fname}_fluosub1.edf'
        y= np.nanmean(np.where(self.avcakefluosub<=0, np.nan, self.avcakefluosub), axis=0)
        bubbleHeader(outfile, self.avcakefluosub, self.tth, self.chi, y, y**0.5)
        np.savetxt(f'{outdir}/{fname}av1d_1.xy', np.array([self.tth, y]).transpose(), fmt = '%.6f')
        self.average2d(fluoK = self.fluoK, outsubdir=outsubdir, fname=fname, **kwargs)
        if saveindividual:
            for item in self.list:
                y = np.nanmean(np.where(item.cake_fluosub<=0, np.nan, item.cake_fluosub), axis=0)
                bubbleHeader(f'{outdir}/{item.basename}.edf', item.cake_fluosub, self.tth,self.chi, y, y**0.5)
    def _getmasks(self,data:np.ndarray,cakemask:np.ndarray|int = 0, nstdevs=3,medianfilter = 4):
        '''
        for applying cosmic masking and mask for final cake
        '''
        masks = np.zeros(shape = data.shape)
        median = np.nanmedian(data,axis=2)
        stdev = np.nanstd(data,axis= 2)
        for i in range(data.shape[2]):
            masks[:,:,i] = np.where((data[:,:,i]< median-nstdevs*stdev)|(data[:,:,i]> median+nstdevs*stdev) | (data[:,:,i]<=0) |
                                    (data[:,:,i] > medianfilter*median)|(data[:,:,i] < median/medianfilter), 1, cakemask)
        return masks

    def plotAll1d(self):
        plt.figure()
        for file in self.list:
            label = os.path.basename(file.fname)
            plt.plot(file.tth ,file.y,label = label)
        plt.plot(self.x,self.yav,label = 'average', linewidth = 3)
        plt.legend()
        plt.xlabel('2theta (°)')
        plt.ylabel('intensity')
        plt.show()

    def saveEDF_noheader(self,dirname):
        '''
        so files can be read with silx/pyfai to make masks
        '''
        im = EdfImage(self.avarray)
        im.save(f'{dirname}/av2d_noheader.edf')
    def saveAllPonis(self,subdir = 'poni'):
        for item in self.list:
            item.saveponi(subdir)

    def __getitem__(self, key):
        return self.list[key]
    def __contains__(self, item):
        return item in self.list
    def __iter__(self):
        return iter(self.list)

