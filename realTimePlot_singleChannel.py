#-*- encoding: utf-8 -*-

"""
realTimePlot_singleChannel.py 

Show an animation "chunk by chunk" in real-time of both
              1)   an EEG signal 
              2)   the spectrogram of the EEG

             For a single channel of a EEG that is extracted from an EDF file (not with LSL streaming yet)

WARNING: THE SIGNAL OBJEcT of EEG_EDFsignal class is doomed to be thoroughly modify in a newer script that must support multi-channel features. 

"""
## to override an abstract method of an abstract class : It uses a decorator. 
from abc import ABCMeta, abstractmethod    
import time 
import sys
from collections import deque               

import numpy as np
from numpy import pi, sin, cos, zeros, ones, log2, ceil

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.image as Image
from matplotlib.lines import Line2D   

from scipy.signal import spectrogram #, stft, periodogram, welch, cwt, morlet

## # personal class for signals management: used to import EEG signal in EDF format 
import signalClass as sC  


 
### ================================================================================
    
### deque handling : 
def concatenateDeque(D, ny):
    """
    D : deque of arrays nx x ny 
    Default : concatenation along nx 
    
    nD = len(D)
    np.array(D).shape = (nD, nx, ny)
    A.shape =  (nD*nx, ny)
    """
    #if ny == 1 :
    ## convert in a arrays with shape = (nx,)   instead of (1,nx)
    concat =np.concatenate(np.array(D))       
    return concat    
    
    #elif ny > 1 :
    #    concat =np.concatenate(np.array(D))  
    #    return concat



class baseRealTimePlot(metaclass=ABCMeta):          
    ### ================================================================================
    ###  A abstract class. Cannot be instanced directly. The sendDataToUpdate() method 
    ### is abstract (empty) and is defined only in derived classes.      
    ### ================================================================================
    
    """
    Draw the animation of a signal  chunk by chunk, for chunk of arbitrary size
    And for signal chunk OR image chunk 
    
    data = (tchunk, ychunk)  
    ychunk size :  (nt, ny)
    
    When we plot a signal, ny = 1, ychunk is an 1D array  of float or a list of float. 
    Whe we plot/show an image :  ny >1 , ychunk is a 2D array 
    
    In any cases: tchunk is an 1D array  of float or a list of float
    

    Attributes: 
        maxt
        dt
        nt
        ny
        ymin
        ymax
        tdata
        ydata
        fig
        ax
        line

    Methods: 
        setDim(self, dim)
        updateDeque(self, tdata, ydata, data)
        initializeDeque(self,tchunk, ychunk )
        getRightChunkLeftPoint(self, tdata, nt)
        
        showAnimation()     : used                   
        update( data)       An abstract method
        sendDataToUpdate()   An abstract method 
        
        (abstract methods in python : Equivalent to C++ virtual functions )
        About abstract methods :   https://stackoverflow.com/questions/5856963/abstract-methods-in-python
        
    """
    def __init__(self, dim = (100,1), maxt=100, dt=1, ymin = -1, ymax = 1):  
        '''
        Everything that is common to cases ny==1 and ny>1
        
        What is not common is placed in the__init__() of derived classes : 
            for ny==1 :  
            
            for ny>1  :
        '''
        ## when yChunk is 1-dim:  line object: 
        ## case ny==1:  create a Line2D object to be added to ax object . 
        ##              tdata, ydata    
        ## case ny > 1  create an AxesImage object with pcolorfast , and we add it to the ax object. 
        ##              tdata,  fdata, Sdata  
        
        
        ## set self.ny and self.nt
        self.setDim(dim)
        self.maxt = maxt 
        self.dt = dt
        self.ymin = ymin
        self.ymax = ymax    
        
        ## we require the chunks to enter a entire number of time in the total width of the plot. 
        assert maxt/dt%self.nt ==0 
           
        ## figure 
        self.fig, self.ax = plt.subplots()   
        self.ax.set_ylim(self.ymin , self.ymax)
        self.ax.set_xlim(0, self.maxt)        
        
        ## titles and axis labels have to be set in derived classes         

    def setDim(self, dim):
        
        ## data dimensions 
        if type(dim) == int: 
            ## data is a  1D signal  
            self.nt = dim
            self.ny =1       
        elif len(dim)==1:    
            self.nt = dim[0] 
            self.ny =1                   
        elif len(dim)==2:    
            ## data is an image
            self.nt = dim[0]            
            self.ny = dim[1] 
        else:
            print('format error  !! ')

    def initializeDeque(self,tchunk, ychunk ):
        '''   
        In each derived class, ydata may have a different structure 
        
        tchunk and ychunk must be time-reversed order. 
        '''
        tdata = deque() 
        ydata = deque()          
        tdata.appendleft(tchunk)
        ydata.appendleft(ychunk)
        return tdata, ydata
    
    
    def updateDeque(self, tdata, ydata, data):
        '''
        tdata and ydata are not member of the basic class : they are received as parameter and then returned. 
        '''
        
        if len(data) == 2 :
            tChunk, yChunk = data            
        else: print('format error !!  ')

        t0 = tdata[-1][-1] 
        t1 = tdata[0][0]

        if t1 - t0 >=self.maxt:
            assert len(ydata)==len(tdata ) 
            tdata.pop()            
            ydata.pop()
        
        tdata.appendleft(tChunk)   
        ydata.appendleft(yChunk)              
    
        ##  AT this point: everything in deque is reversed, but deque itself is not reversed again (i.e. deque is NOT in increasing order )
        ## REM:  self.tdata[-1]  corresponds to the first data  pushed in the deque (a chunk or a single point).
        if self.nt==1 :
            self.ax.set_xlim(tdata[-1], tdata[-1] + self.maxt)         
        elif self.nt > 1 :        
            self.ax.set_xlim(tdata[-1][-1], tdata[-1][-1] + self.maxt)   #  self.tdata[-1] : the most left chunk :
        self.ax.set_ylim(self.ymin, self.ymax)        
        
        return tdata, ydata

        
    def getRightChunkLeftPoint(self, tdata, nt=None):
        '''  return the most left point  in the last chunk (i.e. the most right ), in case nt==1  or nt > 1
        
        By default, nt is self.nt : it is only valid in the signal case. 
        In the image case, nt must be the time length of the image chunk, not of the signal chunk.
        '''
        if nt is None: nt =self.nt
        
        if nt>1:
            ## Each element of tdata deque should have time reversed.
            ## tdata[0] : the last  ('most right') chunk that have been pushed in the deque. It is a list or a 1D array, that contains nt points.
            
            # Check that time is reversed in the chunk
            assert tdata[0][-1]  <= tdata[0][0]    #  
                
            return tdata[0][0]
            
        if nt==1:
            ## tdata[0] : the last point in the deque, a scalar. 
            ## index of the NEXT point to pick to push it in the deque. 
            return tdata[0]  
        

    @abstractmethod 
    def update(self, data):      
        ''' To be override 
        Data is get from sendDataToUpdate throught the FuncAnimation function
        '''
        pass

        
    @abstractmethod 
    def sendDataToUpdate(self):
        ''' To be override '''
        pass
    
    def showAnimation(self):
        
        ## a generator is passed from "sendDataToUpdate" to the update function. 
        ani = animation.FuncAnimation(self.fig, self.update, self.sendDataToUpdate, interval=self.dt )     
        plt.show()  
    
 
### ================================================================================

##===========================================================================
class realTimePlot():     
    """ draw the animation of a signal point by point (and not chunk by chunk)    """
    def __init__(self, maxt=100, dt=0.1, ymin = -1, ymax = 1):  
        nt=1
        super().__init__(nt, maxt, dt, ymin,ymax)
        
  
### ================================================================================
class realTimePlot_fct( realTimePlot):    
        
    def sendDataToUpdate(self):
        '''return the next point (t, y(t)) of myFunction()
        Essentially for testing.
        '''
        
        t=self.tdata[0]
        while True:
            t += self.dt
            y = myFunction(t)  
            yield t, y   # return the generator for a 2-tuple 

def myFunction(x):
    return sin(x)  ## example 

### ======================================================

class realTimePlotSignal( baseRealTimePlot): 
    """
     For signal (  For example:  EEG from EDF file )
     Show in real-time the signal, either point by point (nt=1) or chunk by chunk (nt>1)
     
     Can support signals with multiple channels, but only once is displayed at a time. 
     
     The object Signal is provided via parameters :  
     Its required attributes are tSignal, ySignal, 
    """
    
    def __init__(self,Signal, nt=50, maxt=5000):
        
        ## baseRealTimePlot initializer
        super().__init__(nt, maxt)   
        
        assert hasattr(Signal, 'tSignal') and hasattr(Signal, 'ySignal')
        
        self.Signal = Signal
        tSignal =Signal.tSignal  
        #ySignal =Signal.ySignal
        
        ## Assuming that the time is the same between each points.         
        ## initialization of dt here: 
        self.dt = tSignal[2]-tSignal[1]

        ##  chunk with increasing time  (timenot reversed )
        tchunk = np.arange(self.nt)    #  # because matplotlib does not support generators 
        ychunk = [0.0]*self.nt                         # a list [0.0, 0.0, ..... 0.0], with nt times 0.0
        
        
        
        ## In parent class : 
        self.tdata, self.ydata = self.initializeDeque( tchunk[::-1], ychunk[::-1] ) # chunks are placed in deque in reversed order
        
        ## 
        self.line = Line2D(self.tdata, self.ydata)
        self.ax.add_line(self.line)
        
    def update(self, data):
        ## override an abstract method in the basic abstract class 
        
        ##  data is get from sendDataToUpdate throught the FuncAnimation function  , data = (t,y)
        ## The deque tdata contains maxt/dt points, that means maxt/(dt*nt) chunks of nt points       
    
        assert self.ny == 1   # case for 1D signal 

        ## maximum is taken over ALL the elements of the lists or array that are piped in the deque ydata             
        ymaxNew = np.max(self.ydata)   
        yminNew = np.min(self.ydata)            
    
        if ymaxNew  > self.ymax:  self.ymax = ymaxNew   
        if yminNew  < self.ymin:  self.ymin = yminNew 
        
        self.tdata, self.ydata = self.updateDeque(self.tdata, self.ydata, data)

        self.line.set_data(self.tdata, self.ydata)
        return self.line   
    
    
        
    def sendDataToUpdate(self):
        '''return the next signal value  when the entire signal is already available.
        
        ySignal must be a unique time serie (ny=1), (i.e. from a single channel), and not an image or multiple channels. 
        
        '''
        assert self.ny==1
        tSignal =self.Signal.tSignal  
        ySignal =self.Signal.ySignal
        
        ## return a signal data point (case nt=1) or a signal chunk (case nt > 1)
        while True:            
            ## The deque tdata contains maxt/dt points, that means maxt/(dt*nt) chunks of nt points
            ## index :  index in tSignal that represents the beginning of the next element in the deque. 
            ## index = t in tSignal
            ## An element in deque is either a signal data point (case nt=1) or a signal chunk (case nt > 1)         
            
            ## Reverse the chunk along time axis : to be in the same order than the chunks in the deque
            index = int(  self.getRightChunkLeftPoint(self.tdata)/self.dt ) 
            
            tChunk = tSignal[index + self.nt : index:-1 ] 
            yChunk = ySignal[index + self.nt : index:-1 ]
            yield tChunk, yChunk   # "yield" returns a generator for a 2-tuple 




### ==============================================================================================

class realTimePlotSignal_EDF( realTimePlotSignal):       
    """
    It works well for signals that have only a single channel. 
    
    For version that deals with multi-channel signals: see RTmultiChannels.py file. 
    """

    def __init__(self, nt=50, maxt=5000, iChan =0,fileEDF = 'Subject00_1.edf', datapath = 'C:/Users/Moi2/Documents/DONNEES/EEG_physionet/EEG_arithmeticTasks/'):
       
        ## EEG_EDFsignal instance : An attribute of the class realTimePlot_EDF
        Signal = sC.EEG_EDFsignal(iChan,fileEDF ,datapath )    
        super().__init__(Signal, nt, maxt )

    
### ==============================================================================================
### ==============================================================================================
class RTspectgram(baseRealTimePlot):            
    """
    for any signal, not only EEG EDF file. 
    
    From a chunk of single-channel signal :( i.e.  ySignalChunk  )
    Compute spectrogram ----->  imageChunk to draw : shape = (nt, ny)
    Animate the display of ImageChunk 
    
    
    self.ntIm  is only introduced with this class. It is not included in baseRealTimePlot
    (ntIm is to be use instead of nt in several places. )
    ntIm :  length in time of spectrogram image 
    
    """
  
    def __init__(self, Signal, nt=600, ny = 512, maxt=6000, Nwin =2):
        '''
        nt :  signal chunk length , Not the length of image chunk
        ntSpec:  length of image chunk (i.e. chunk of spectrogram image)
        
        ny :  height of the spectrogram image.  
        
        Signal: an object that has :  tSignal, ySignal, selectedChannel attributes
        Signal.tSignal and Signal.ySignal are 1-D arrays 
        '''
    
        ## Check the attributes of the Signal object. 
        assert hasattr(Signal, 'tSignal') and hasattr(Signal, 'ySignal')
        self.Signal = Signal
        
        ## adjust NPERSEG and ny = NFFT so that NFFT is a power of 2, and NFFT >= NPERSEG.
        self.NPERSEG, self.NFFT = self.adjustNFFT(nt, ny, Nwin)   
        
        ##  for display of frequency ( i.e. y) axis : 
        ymin = 0.5 # Hz
        ymax = 30  # Hz 
        dt = 1

        ## baseRealTimePlot initializer : initialize: 
        dim = [nt,ny]         
        super().__init__(dim, maxt, dt, ymin, ymax ) 
        
        ## initialization
        self.ntIm  =2  

        ## signal chunk at index =0. Time is in ms in tSignalChunk and in t. 
        tSignalChunk, ySignalChunk= Signal.returnSignalChunk(0, self.nt, 'non-reversed')
        t, f, Sxx = self.computeSpectgram_scipy(tSignalChunk, ySignalChunk, PropOverlap=0.5)
        
        ## initialization of ny, fmax, fmin  
        #assert self.ntIm == len(t)
        self.ny= len(f)     
        assert f[0] == min(f)
        assert f[-1] == max(f)        
        self.fmax = max(f)
        self.fmin = min(f)
        self.fdata = f
        
        ## pcolorfast is faster than imshow or pcolormesh, but it does not have shading or interpolation 
        ## returns an image.AxesImage object 
        self.image = self.ax.pcolorfast(t,f , Sxx )    
        self.ax.set_xlabel('t        [ms]')
        self.ax.set_ylabel(' frequency     [Hz]')
        self.ax.set_title('Channel: %s'% Signal.selectedChannel ) 

        ## the extent of the data we want to display (ymin and ymax : apriori not the same as self.fmin and self.fmax)
        self.ax.set_ylim(self.ymin, self.ymax)          
        self.ax.set_xlim( 0,  self.maxt )   #        
    
        ##  extend of actual data.  fmin = min(f) , and fmax = max(f)  
        extent = ( 0, self.ntIm*dt, self.fmin, self.fmax )  
        self.image.set_extent(extent) # extent is data axes (left, right, bottom, top) for making image plots    
                
        ##  Initialization after computation of first chunk spectrogram , ny et ntIm are previously defined, but NOT ymax and ymin
        ## We reverse the chunks only after taking the spectrogram: 
        t2 = t[::-1]
        f2 = f[::-1]
        S2 = Sxx[:,::-1]      
        
        ## initialization
        ## After reversing each chunks separately,  deque themselves have then to be reversed so that tdata, ydata,  are in increasing order: 
        self.tdata, self.Sdata = self.initializeDeque(t2, S2.T)   # S2 was not already transposed;  but t2 and S2 are reversed in time. 
                


    def adjustNFFT(self, nt, ny, Nwin):
        ''' set  NFFT NPERSEG so that: 
        1)  self.ny = len(f) = nfft/2 +1, because the FFT is one-side, we have NFFT = 2*(ny -1)
        2)   NPERSEG <=NFFT    
        3)   NFFT is a power of 2 
    
    
         nt :  time length of signal chunk (not of spectrogram chunk)
         Nwin : number of FFT windows desired by signal time chunk. : (However: not exact when windows overlap )
        '''
    
        NFFT = 2*(ny -1)
        NPERSEG= int(nt/Nwin )  

        NFFT = max(NFFT, NPERSEG)
        exponent = int(ceil(log2(NFFT)) )
        NFFT = int(2**exponent)  # with zero padding         # number of samples in the FFT,i.e. in frequency.  --->  NFFT/2 +1 
        
        return NPERSEG, NFFT
        

    def computeSpectgram_scipy(self, tSignalChunk, ySignalChunk, PropOverlap = 0.5):
        """
        Using the scipy.signal library :
            from scipy.signal import spectrogram 
        Alternative possibility: 
            from matplotlib.Axes import specgram
            
            
        NFFT:  Length of the FFT used, (eventually zero padded )
        """
        
        ## Sampling frequency of the t time series:  in order to have frequency in Hertz while time is in ms: 
        fs = 1000/self.dt     # Defaults to 1.0.   [Hz]  --->  1kHz
             

        ## The method adjustNFFT() has set NPERSEG and NFFT in the following way: 
        ## 1)  self.ny = len(f) = nfft/2 +1, because the FFT is one-side, we have NFFT = 2*(ny -1)
        ## 2)   NPERSEG <=NFFT  
        ## 3)   NFFT is a power of 2 
        NFFT = self.NFFT   
        NPERSEG = self.NPERSEG
        
        window0=('tukey', PropOverlap)   # i.e. Tukey window with 1/4 of a window’s length overlap at each end, i.e. overall : 1/2 overlapped   
        SCALING = 'spectrum' #scaling : { 'density', 'spectrum' }
               
        ##  REM  tSignalChunk is not directly involved itself in spectrogram function. 
        f, t, Sxx  = spectrogram(ySignalChunk, fs, window=window0, nperseg=NPERSEG,nfft = NFFT, detrend='constant', return_onesided=True, scaling=SCALING, axis=-1, mode='psd')
                
        ## because time chunk in scipy spectrogram is  in sec , since frequency is in Hz. We want it in ms. 
        t = t*1000
        return t, f, Sxx
        
   
    def computeSpectgram_mpl(self, tSignalChunk, ySignalChunk, Nwin, PropOverlap=0):
        ## NOT USED YET... but we could do it to compare.     
        """Axes.specgram(x, NFFT=None, Fs=None, Fc=None, detrend=None, window=None, noverlap=None, cmap=None, xextent=None, pad_to=None, sides=None, scale_by_freq=None, mode=None, scale=None, vmin=None, vmax=None, *, data=None, **kwargs)
        https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.specgram.html
        """        
        pass
    
    def itIsTheEnd(self):
        ''' 
        Return a boolean value that says if the end of the signal is reach. 
        
        Criterion 1 : is the two extremities of the last signal chunk is the same point ?
        '''
        criterion1 = self.tdata[0][-1]  == self.tdata[0][0] 
        
        ## criterion2 = Signal.isEndReach() 
        
        return criterion1
    
        
    def sendDataToUpdate(self):   
        '''      
        
        tSignalChunk :  time of original signal chunk 
        tImChunk     :  time of spectrogram image chunk
        '''
        ntIm = self.ntIm
        while True:  
            
            ## A criterion fo End detection is called
            if self.itIsTheEnd(): 
                ## pause for 2 min before closing
                print('----------  End of signal is reach: it will be closed in 2 minutes  -------------')
                
                time.sleep(120)
                sys.exit()
                #break 
            
            
            index = int(  self.getRightChunkLeftPoint(self.tdata, ntIm)/self.dt ) 
            #print(index) #   advance by step of ntIm
 
            tSignalChunk, ySignalChunk= self.Signal.returnSignalChunk(index, self.nt, 'non-reversed')
            tImChunk, fImChunk, SChunk = self.computeSpectgram_scipy( tSignalChunk, ySignalChunk)
            
            ##  For the last image, in case the total length of ySignal is not a multiple of ntIm
            if len(tImChunk) < ntIm:             
                new_tImChunk = zeros([ntIm])          # shape = (ntIm,)
                new_SChunk = zeros([self.ny, ntIm])   # shape = (ny, ntIm)
               
                new_tImChunk[: len(tImChunk)]  = tImChunk
                new_SChunk[:, :len(tImChunk)]  = SChunk

                tImChunk = new_tImChunk
                SChunk = new_SChunk
                
            ## "yield" returns a generator for a 3-tuple to make them an image object 
            yield tImChunk, fImChunk, SChunk 
    

    def update(self, data):         
        '''  
        data is get from sendDataToUpdate throught the FuncAnimation function  , data = tChunk, fChunk, SChunk
        The deque tdata contains maxt/dt points, that means maxt/(dt*nt) chunks of nt points    
        We assume   Chunks to be in reversed-time order 
        '''
        
        ## format check 
        if(len(data) == 3) : 
            tChunk, fChunk, SChunk = data  

            ## must be the same 
            #print(self.fdata)            
            #print(fChunk[::-1])
        else: 
            print('format error !! ')  # should raise an error more properly *****
         
        ## dimension check  
        if tChunk.shape[0] !=self.ntIm:   
            tChunk = tChunk.T
        if SChunk.shape[0] !=self.ntIm:   
            SChunk = SChunk.T      
            assert SChunk.shape[1]==self.ny
            
        ## check the time -reversed order. 
        if tChunk[-1]> tChunk[0]: 
            tChunk = tChunk[::-1]
            SChunk = SChunk[::-1, :]  

        ## In the case of the spectrogram : each tChunk time begins to 0, because the spectrogram does not take tChunk into account but only fs, i.e. ds.    
        ## self.tdata[0][0] is the max value of the previous tChunk. The current tChunk is added after it. 
    
        ## There is a bug if I try:  tChunk += self.tdata[0][0],   tdata itself was modified by incrementing tChunk !!
        ## No bug with this simple syntax... 
        tChunk = tChunk + self.tdata[0][0]         
        
        ## axes adjustments   
        self.tdata, self.Sdata = self.updateDeque(self.tdata, self.Sdata, [tChunk,SChunk] )     #*  SHOULD BE MOVED *****

        ## concatenation         
        tConcat = concatenateDeque(self.tdata, 1)
        SConcat = concatenateDeque(self.Sdata, self.ny)
        
        ## reversed wrt time : 
        tConcat = tConcat[::-1]
        SConcat = SConcat[::-1, :]

        ntIm = self.ntIm 
        dt = self.dt
        assert len(tConcat)%ntIm == 0
              
        self.ax.set_xlim( tConcat[0], tConcat[0]+ self.maxt)         
        self.ax.set_ylim(self.ymin, self.ymax)  
        self.image.set_data(  SConcat.T )    
        extent = ( tConcat[0], tConcat[-1], self.fmin, self.fmax )
        self.image.set_extent(extent)   
        
        return self.image   

    
    def test_update(self):
        ''' 
        Only for testing 
        Test a single step of the update loop
        '''
        index = int(  self.getRightChunkLeftPoint(self.tdata, self.ntIm)/self.dt ) 
        tSignalChunk, ySignalChunk= Signal.returnSignalChunk(index, self.nt, 'non-reversed')
        tImChunk, fImChunk, SChunk = self.computeSpectgram_scipy( tSignalChunk, ySignalChunk)
        
        data = (tImChunk, fImChunk, SChunk )
        self.update(data)
        plt.show()        

## ====================================================================================
##   Test functions 
## ====================================================================================
def test_displayFigure(tSignal, ySignal,t, f, Sxx):
    '''
    Only for testing 
    It simply displays signal and spectrogram for a static first chunk of signal data. 
     Time axis of signal and spectre are not properly aligned: very ugly but not a problem for a simple test 
    '''
    ## y limit in the plot
    ymin =0.5 # Hz
    ymax =50  # Hz 
    #fs = 1000/dt # Sampling frequency of the t time series. Defaults to 1.0.   [Hz]  --->  1kHz 


    plt.close()
    plt.clf 

    plt.figure(1)    
    plt.plot(tSignal, ySignal)
    
    plt.figure(2)

    ## the one we used in update() for now
    ## return an AxisImage object
    ## no interpolation or shading (which would be prettier)
    ax1=plt.subplot(311)
    extent = ( t[0], t[-1], f[0], f[-1] )
    #print(extent)
    image = ax1.pcolorfast(Sxx )  
    ax1.set_ylim(ymin,ymax )     # and not f[0], f[-1]
    image.set_extent(extent)       
    
    ## returns a collections.QuadMesh object:  
    ## works, , but difficult to adapt since it DOES NOT WORK with AxisImage
    ax2=plt.subplot(312)
    plt.pcolormesh(t, f, Sxx, shading='gouraud' )
    plt.ylim(ymin,ymax ) 
 
    
    ## return an AxisImage object; and is slower than pcolorfast()
    ## But it allows interpolation
    # affichage le long de l'axes des y ne MARCHE PAS !!                  ******************  
    ax3=plt.subplot(313)
    extent = ( t[0], t[-1], f[0], f[-1] )
    print(extent)
    image = ax3.imshow(Sxx, extent = extent, interpolation='gaussian')  
    
    ax3.set_xlim(t[0], t[-1])
    ax3.set_ylim(ymin,ymax )    # and not f[0], f[-1]
    image.set_extent(extent)   
    
    plt.show()


def testMainSpectDisplay():
    '''
    It doesnt use spectro.update method 
    To test display with imshow, or pcolormesh, or pcolorfast
    '''
    maxt = 6000 ; ntSignal = 600;  ny = 513
    Signal = returnOurSignal()     
    spectro = RTspectgram(Signal, ntSignal, ny, maxt)
    
    tSignalChunk, ySignalChunk= Signal.returnSignalChunk(0, ntSignal, 'non-reversed')
    tImChunk, fImChunk, SChunk = spectro.computeSpectgram_scipy( tSignalChunk, ySignalChunk)
    
    ## to show:  the signal , and the spectrogram with imshow, or pcolormesh, or pcolorfast
    test_displayFigure(tSignalChunk, ySignalChunk,tImChunk, fImChunk, SChunk)  

def testSinus():
    ''' Test plotting of a simple fuction (a sinus) '''
    dt = 0.1
    maxt = 4*pi
    scopeSinus = realTimePlot_fct(maxt, dt)
    scopeSinus.showAnimation()


    
def testMain_update():
    ''' To test update(), without using animation '''
    maxt = 6000 ; ntSignal = 600;  ny = 513
    Signal = returnOurSignal()
    spectro = RTspectgram(Signal, ntSignal, ny, maxt)
    spectro.test_update()


    
### =========================================================================================
###  Signal  , see also the signalClass file 
### =========================================================================================
def returnOurSignal():
    '''return the channel #0 of an EEG extracted from an EDF file that has been picked on physionet.org 
       Return a Signal object, with tSignal and ySignal attributes.  
    
    '''
    
    fileEDF = 'Subject00_1.edf'
    datapath = 'C:/Users/Moi2/Documents/DONNEES/EEG_physionet/EEG_arithmeticTasks/'
    iChan =0
    Signal = sC.EEG_EDFsignal(iChan,fileEDF, datapath )
    return Signal

### =========================================================================================
###  To run a monitor : with the use of multiprocessing or multithreading 
### =========================================================================================
def animProcessPool(obj):
    '''  
    To lauch showAnimation() for an object obj by using a process pool. 
    It does not work if this function id placed inside the def animate_Processes()
    '''
    obj.showAnimation()

def animate_Processes(objList, mode):
    '''
    Use separate processes to run objects contained in objList. 
    
    objList : list of multiple object that are assumed of compatible type.
    '''
    from multiprocessing import freeze_support, Pool, Process
        
    freeze_support()
    if mode == 'Pool':
        ## using a pool of maximum 4 processes to treat object in objList
        Npool = min( 4, len(objList) )
        
        ## run processes that run animation for each object. 
        with Pool(Npool) as p:
            #processList = p.map(animProcessPool, objList)  
            processList = p.map(lambda o: o.showAnimation(), objList)  
            
        
    if mode == 'Process':
        
        ## list of processes that run the animation of each object 
        processList = [Process(target= o.showAnimation ) for o in objList ]

        for p in processList:   p.start()
        for p in processList:   p.join()
  
  
def animate_Threads(objList, mode):
    '''
    Use separate threads to run objects contained in objList. 
    
    objList : list of multiple object that are assumed of compatible type.
    '''
    import threading as th
    
    pass 


### =========================================================================================
###  Main 
### =========================================================================================

def mainPlotEDFSignal():
    '''  Display an EDF signal in real-time, 
         Use double queue data structures '''
    
    ##   maxt = 6000  ;  ntSignal = 600  ##
    
    ## Signal is returned internally in the __init__ of the class 
    ## and default parameters are the ones of this test 
    signalPlot = realTimePlotSignal_EDF()
    signalPlot.showAnimation()
    
def mainSpectro():
    ''' 
    Use 
    - Double queues as data structures 
    - Spectrogram from scipy
    - pcolorfast from matplotlib
    '''
    maxt = 6000 ; ntSignal = 600;  ny = 513
    Signal = returnOurSignal()
    spectro = RTspectgram(Signal, ntSignal, ny, maxt)
    spectro.showAnimation()  

def main():
    maxt = 6000 ; ntSignal = 600;  ny = 513
    nt=ntSignal # 60
    Signal = returnOurSignal()

    ## two tasks to run either on separate threads or on separate processes.
    ## Here we use multithreading because multiprocessing will be used to deal with the 4 different channels: (one process per channel to treat)
    
    signalPlot=realTimePlotSignal(Signal, nt, maxt)
    spectro = RTspectgram(Signal, ntSignal, ny, maxt)
    
    ## With multiprocessing :  it works as well with Pool than with Process classes. 
    
    mode ='Pool'  # 'Process'
    animate_Processes([signalPlot,spectro], mode)  
    
    
    #animate_Threads([signalPlot,spectro], mode) 
    
if __name__ == '__main__' :
    #mainPlotEDFSignal()
    
    #mainSpectro()
    
    main()
    

   
