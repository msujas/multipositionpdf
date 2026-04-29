import numpy as np
import pyFAI
from fabio.edfimage import EdfImage
from cryio.cbfimage import CbfImage, CbfHeader
import os, fabio
from glob import glob
from pyFAI.integrator.azimuthal import AzimuthalIntegrator
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d, RegularGridInterpolator, LinearNDInterpolator
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

def gainCorrection(avim,gainArray):
    avimGain = avim/gainArray
    avimGain = np.where(gainArray <0, -1, avimGain)
    return avimGain

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

        
class FilePoni():
    def __init__(self, fname:str, ypos:float, ponilist:PoniList, zpos:float = None, maskfile:str = None):
        self.fname = fname
        self.dirname = os.path.dirname(self.fname)
        self.basename = os.path.splitext(os.path.basename(self.fname))[0]
        self.flux = CbfHeader(self.fname).header['Flux']
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
        self.maskfile = maskfile
        match self.ponilist.interpolationDimension:
            case 1: self.interpolatePoni()
            case 2: self.interpolatePoni2D()

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

    def integrate(self, tthmin, tthmax,tthbins=5000,  chimin = -178, chimax=178, chibins = 354, gainfile = None, xyedir = 'xye', 
                  cakedir = 'cake', polarization_factor = 0.85, scale = 10**5):
        self.ai = AzimuthalIntegrator(dist=self.dist, poni1=self.poni1, poni2=self.poni2, rot1=self.rot1,
                                      rot2=self.rot2, rot3=self.rot3, detector=self.detector,wavelength=self.wavelength)

        mask = np.where(self.array < 0,1,0)
        if self.maskfile:
            mask = fabio.open(self.maskfile).data


        if gainfile:
            gainmap = fabio.open(gainfile).data
            self.array = gainCorrection(self.array,gainmap)
        
        self.x,self.y,self.e = self.ai.integrate1d(self.array,tthbins,correctSolidAngle=True, polarization_factor=polarization_factor,method='bbox',
                                            unit='2th_deg', mask=mask, error_model='poisson', radial_range=(tthmin,tthmax))
        self.result2d = self.ai.integrate2d(self.array,tthbins,correctSolidAngle=True, polarization_factor=polarization_factor,method='bbox',
                                            unit='2th_deg', mask=mask, error_model='poisson',radial_range=(tthmin,tthmax), 
                                            azimuth_range=(chimin, chimax), npt_azim=chibins)

        self.array2d = (self.result2d[0])*scale/self.flux
        self.y = self.y*scale/self.flux
        self.tth = self.result2d[1]
        self.chi = self.result2d[2]
        self.integrated= True
        if cakedir:
            self.saveCake(cakedir=cakedir)
        if xyedir:
            self.save1d(xyedir=xyedir)

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
        
    def saveMaps(self,dirname, polarization_factor):
        self.geometry = Geometry(dist=self.dist, poni1=self.poni1, poni2=self.poni2, rot1=self.rot1,
                                      rot2=self.rot2, rot3=self.rot3, detector=self.detector,wavelength=self.wavelength)
        
        self.ttharray:np.ndarray = self.geometry.twoThetaArray()
        self.chiarray:np.ndarray = self.geometry.chiArray()
        self.pol:np.ndarray = self.geometry.polarization(factor=polarization_factor)
        self.sa :np.ndarray = self.geometry.solidAngleArray()
        self.sinchi2:np.ndarray = np.sin(self.chiarray)**2
        self.arrayCorrected = self.array/(self.pol*self.sa)
        basename = os.path.splitext(os.path.basename(self.fname))[0]
        basename += '.edf'
        os.makedirs(f'{dirname}/tth',exist_ok=True)
        os.makedirs(f'{dirname}/sa',exist_ok=True)
        os.makedirs(f'{dirname}/pol',exist_ok=True)
        os.makedirs(f'{dirname}/chi',exist_ok=True)
        im = EdfImage(self.ttharray)
        im.save(f'{dirname}/tth/{basename}')
        im = EdfImage(self.chiarray)
        im.save(f'{dirname}/chi/{basename}')
        im=EdfImage(self.pol)
        im.save(f'{dirname}/pol/{basename}')
        im=EdfImage(self.sa)
        im.save(f'{dirname}/sa/{basename}')

    def arrayBins(self,tthbins:int, chibins:int, tthmax:float):
        tthmax = tthmax*np.pi/180
        self.tthbinarray = (self.ttharray*tthbins/tthmax).astype(int)
        self.chibinarray = (self.sinchi2*chibins).astype(int)
        self.combinedbinarray = self.tthbinarray+1j*self.chibinarray
        return self.combinedbinarray

class MultiFile():
    def __init__(self,alist:list[FilePoni]):
        self.list = alist
    def saveMaps(self,dirname):
        for f in self.list:
            f.saveMaps(dirname)
    def integrateAll(self,tthmin,tthmax):
        poniinterp = self.list[0].ponilist.interpolationDimension
        print(f'integrating with ponis interpolated in {poniinterp}d and integrating')
        for file in self.list:
            file.integrate(tthmin,tthmax)
    def average1d(self,x0,xend,npoints,  outsubdir= 'xye', fname='', **kwargs):
        '''
        run the interpolations (if not done already) and regrid and average all 1d patterns
        x0 - 2theta0
        xend - 2theta end
        npoints - number of points in regrid
        outsubdir - subdirectory for averaged 1d file
        kwargs - kwargs for interpolation
        '''
        poniinterpolation = self.list[0].ponilist.interpolationDimension

        if not self.list[0].integrated:
            print(f'interpolating ponis in {poniinterpolation}d and integrating files')

            for file in self.list:
                file.integrate(x0,xend,tthbins=npoints,**kwargs)
        
        self.x = self.list[0].x
        avarray = np.empty(shape=(len(self.x), len(self.list)))
        for i,file in enumerate(self.list):
            #gridfunc = interp1d(file.x,file.y, fill_value=np.nan, bounds_error=False)
            avarray[:,i] = np.where(file.y <=0, np.nan, file.y)
        self.yav = np.nanmean(avarray,axis=1)
        if outsubdir:
            outdir = f'{os.path.dirname(self.list[0].fname)}/{outsubdir}'
            os.makedirs(outdir,exist_ok=True)
            outfile = f'{outdir}/{fname}av1d.xy'
            np.savetxt(outfile,np.array([self.x,self.yav]).transpose(),fmt='%.6f')
        return self.x,self.yav
    def average2d(self,outsubdir='cake',fname='', nstdevs=3, medianfilter=4, cakemask=None):
        print('averaging cakes')
        if cakemask:
            cakemask = fabio.open(cakemask).data
        else:
            cakemask=0
        for i,item in enumerate(self.list):
            if i == 0:
                allarrays = np.empty(shape=(*item.array2d.shape,len(self.list)))
            allarrays[:,:,i] = np.where(item.array2d <= 0 ,np.nan, item.array2d)
        #allarrays = np.where((allarrays <=0)| (cakemask == 1), np.nan, allarrays)
        masks = self.getmasks(allarrays, cakemask, nstdevs=nstdevs,medianfilter=medianfilter)
        allarrays = np.where(masks == 1, np.nan, allarrays)
        self.avarray = np.nanmean(allarrays,axis=2)
        self.avarray = np.where(np.isnan(self.avarray),0, self.avarray)
        self.tth = self.list[0].tth
        self.chi = self.list[0].chi
        self.ycake = np.nanmean(self.avarray,axis=0)
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
    def getmasks(self,data:np.ndarray,cakemask = None, nstdevs=3,medianfilter = 4):
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
            plt.plot(file.x,file.y,label = label)
        plt.plot(self.x,self.yav,label = 'average', linewidth = 3)
        plt.legend()
        plt.xlabel('2theta (°)')
        plt.ylabel('intensity')
        plt.show()

    def regrid2d(self,tth0, tthend, tthpoints, chi0=-178, chiend=178, chipoints=354, nstdevs = 3, medianFilter = 4, 
                 outsubdir = 'cake', cakemask:str = None):
        print('regridding and merging cakes')
        tthgrid = np.linspace(tth0,tthend, tthpoints)
        chigrid = np.linspace(chi0,chiend,chipoints)
        mesh = np.empty(shape=(len(chigrid),len(tthgrid),2))
        for i in range(len(chigrid)):
            for j in range(len(tthgrid)):
                mesh[i,j] = [chigrid[i],tthgrid[j]]
        ndlist = np.empty(shape = (*mesh.shape[:2],len(self.list)))
        for i,item in enumerate(self.list):
            interparray = np.where(item.array2d<=0, np.nan, item.array2d)
            rgf = RegularGridInterpolator((item.chi,item.tth), interparray, fill_value=np.nan, bounds_error=False)
            newdata = rgf(mesh)
            ndlist[:,:,i] = newdata

        if cakemask:
            cakemask = fabio.open(cakemask).data
        else:
            cakemask = 0
        masks = self.getmasks(ndlist, cakemask, nstdevs=nstdevs,medianfilter=medianFilter)

        maskeddata = np.where(masks == 1, np.nan, ndlist)
        y2 = np.empty(shape=(tthpoints,maskeddata.shape[2]))
        for i in range(maskeddata.shape[2]):
            y2[:,i] = np.nanmean(maskeddata[:,:,i],axis = 0)
        self.y2 = np.nanmean(y2,axis=1)
        self.avarray = np.nanmean(maskeddata,axis= 2)
        self.x = tthgrid
        self.y = np.nanmean(self.avarray,axis = 0)
        e = self.y**0.5
        baseoutdir = os.path.dirname(self.list[0].fname)
        outsubdir = f'{baseoutdir}/{outsubdir}'
        if outsubdir:
            os.makedirs(outsubdir,exist_ok=True)
            bubbleHeader(f'{outsubdir}/av2d.edf', np.where(np.isnan(self.avarray),0,self.avarray), tthgrid, chigrid,self.y,e)
            np.savetxt(f'{outsubdir}/av2d.xy',np.array([self.x,self.y2]).transpose(),fmt='%.6f')
            np.savetxt(f'{outsubdir}/cake2d.xy',np.array([self.x,self.y]).transpose(),fmt='%.6f')
        return self.avarray
    def saveEDF_noheader(self,dirname):
        '''
        so files can be read with silx/pyfai to make masks
        '''
        im = EdfImage(self.avarray)
        im.save(f'{dirname}/av2d_noheader.edf')

    def __getitem__(self, key):
        return self.list[key]
    def __contains__(self, item):
        return item in self.list

