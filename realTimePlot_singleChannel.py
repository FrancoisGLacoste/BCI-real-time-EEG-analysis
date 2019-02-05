# -*- encoding: utf-8 -*-

"""
test_realTimeEEG_EDFfile.py 

Show an animation chunk by chunk (and hopefully in real-time) of 
              1)   an EEG signal that is extracted from an EDF file. 
              2)   the spectrogram of the EEG, chunk by chunk in (hopefully) real-time   ********************

THIS PROGRAM WORKS ONLY FOR EEG FILE of EDF FORMAT.  NOT WITH LSL yet. 
THE REAL_TIME ANIMATION WORKS ONLY FOR SIGNAL WITH SINGLE CHANNEL   --------------->  See the newer python script

THE SIGNAL OBJEcT of EEG_EDFsignal class is doomed to be thoroughly modify in a newer script that must support multi-channel features. 



"""
from abc import ABCMeta, abstractmethod    # to override an abstract method of an abstract class : It uses a decorator. 
#from time import sleep
from collections import deque               

import numpy as np
from numpy import pi, sin, cos, zeros, ones, log2, ceil

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.image as Image
from matplotlib.lines import Line2D


from scipy.signal import spectrogram  #, stft, periodogram, welch, cwt, morlet



##==========================================================================
### EEG  signal from EDF file
import pyedflib

class EEG_EDFsignal():       
    """
    A class to manage EEG signals found in EDF format 
    
    format edf :  european data format 
    About pyEDFlib :   https://pypi.org/project/pyEDFlib/  (install with Anaconda:  conda install -c conda-forge pyedflib)
    https://pyedflib.readthedocs.io/en/latest/
    REM: MNE software can also open EDF files    
    """

    def __init__(self, iChan=0, fileEDF = 'Subject00_1.edf',datapath = 'C:/Users/Moi2/Documents/DONNEES/EEG_physionet/EEG_arithmeticTasks/'):
        '''
        Default :  a specific edf file and its path are selected. 
        
        attributes: 
        f
        signals
        tSignal
        ySignal
        iChan
        selectedChannel
        '''
    
        ## open EDF file and keep the file identifier as an attribute of the class
        self.f = pyedflib.EdfReader(datapath+fileEDF ) #  pyedflib.edfreader.EdfReader object 
        self.loadSignalEDF(iChan)
        
    def setEDFsignals(self):
        '''
        Use self.f 
        And set : 
            self.tSignal 
            self.signals  : including ALL channels. 
            self.signal_labels
        '''           
        nChannels = self.f.signals_in_file   # number of signal channels in the file. (ex 21 , 16, 8, 4  for EEG)
    
        ## We assume that all channel signals have the same number of samples.             
        Nsamples = self.f.getNSamples()[0]   # length (in samples) of the 0-th channel signal, which is assumed to be the same for all signals.        
        signals = zeros((nChannels, Nsamples))      
        
        ## The entire signal is already available:  
        for iChannel in range(nChannels):      # i-th channel considered : 0, 1,....  (nChannels-1) 
            signals[iChannel, :] = self.f.readSignal(iChannel)            
        assert nChannels==signals.shape[0]  
        assert Nsamples==signals.shape[1]  

        self.signals =signals 
        timeMax =Nsamples  # 91000  ONLY if sample_width= 1 ms or 1 s, i.e. fsampling =1 kHz  or 1 Hz (depending on units)
        self.tSignal  = np.arange(timeMax)    
        
        self.signalLabels = self.f.getSignalLabels()
           
    
    def loadSignalEDF(self, iChan):  
        '''  
        Import a complete EEG signal from EDF file. 
        
        Set 
            Selected signal :        self.tSignal, self.ySignal ,self.selectedChannel , self.iChan
            All channel signals:     self.signals,  self.signalLabels
        '''
        ## set self.tSignal ; self.signals; self.signal_labels        
        self.setEDFsignals()            

        ## selection of the channels that are considered 
        self.iChan = iChan 
        self.ySignal =  self.signals[iChan, :]
        self.selectedChannel = self.signalLabels[iChan] 
        print(" The channel %s is loaded"% self.selectedChannel)


    def returnSignalChunk(self, index, nt, mode):
        '''  CORRECT           
        mode can be: 'reversed' ,  'non-reversed'. 
        
        USED in spectrogram, and SHOULD ALSO be used in the the signal plot 
        '''''
        
        
        if mode =='non-reversed':
            tChunk = self.tSignal[index : index +nt ] 
            yChunk = self.ySignal[index : index +nt]
            
        elif mode =='reversed':    
            if index > 0:
                tChunk = self.tSignal[index + nt -1 : index-1:-1 ] 
                yChunk = self.ySignal[index + nt -1 : index-1:-1 ]
            elif index==0: 
                tChunk = self.tSignal[index + nt -1 ::-1 ] 
                yChunk = self.ySignal[index + nt -1 ::-1 ]
                
        else:  
            return 'GROS CACA'                       #   RAISE ERROR **********
        return tChunk, yChunk  

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
    if ny == 1 :
        ## must convert in a arrays with shape = (nx,)   instead of (1,nx)
        #print(np.array(D).shape)                  # 
        #print( np.array(D)         )              #  
        concat =np.concatenate(np.array(D))       #  la concat est OK puor ny==1
        #print(concat)
        return concat    
    elif ny > 1 :
        
        print(ny)        
        print( D[0])
        print( D[1])
        
        #   S2 = Sxx[:,::-1] 
        print(type(D[0]))
        print(len(D))                                   # 1                             2
        print(D[0].shape)                            #(2,513)                  (le dernier "chunk" ) (2, 513)
        print(D[1].shape)                            #(2,513)  
        print(np.array(D)) 
        print(np.array(D).shape)                     # (2, 2, 513)  OK
        
        concat =np.concatenate(np.array(D))  
        print(concat.shape)                        #  (4,513),   i.e. : (nD*ntIm, ny)
        return concat
        #return np.concatenate(np.array(D))  

"""
Pour tdata
np.array(D)[0].shape (2, 2)   //  concatenateDeque.shape = (2,)

pour Sdata : 
np.array(D)[0].shape  = (2,)   !!!!!?????
"""
### ================================================================================
###  A abstract class. Cannot be instanced directly. The sendDataToUpdate() method is abstract (empty) and is defined only in derived classes.      
class baseRealTimePlot(metaclass=ABCMeta):          #         
    """
    An abstract class.  
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

        getRightChunkLeftPoint()  
        update( data)         AN abstract method
        showAnimation()                  
        sendDataToUpdate()  :    abstract method required to send data to update method via animation.FuncAnimation()
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
        '''   Called from the derived classes.
        In each derived class, ydata may have a different structure 
        
        tchunk and ychunk should have time reversed. 
        '''
        tdata = deque() 
        ydata = deque()          
        tdata.appendleft(tchunk)
        ydata.appendleft(ychunk)
        return tdata, ydata
    
    
    def updateDeque(self, tdata, ydata, data):
        '''tdata and ydata are not member of the basic class : So they are returned 
        
          ***** SERAIT PEUT ETRE PLUS EFFICACE DE DECLARER tdata et ydata dans la classe de base      ****^^????
                                        et d'utiliser self.tdata au lieu de tdata etc... 
        '''
        
        if len(data) == 2 :
            tChunk, yChunk = data            
        else: print('format error !!  ')

        t0 = tdata[-1][-1] 
        t1 = tdata[0][0]

        #if self.dt * len(tdata ) * nt >= self.maxt:   # NE FONCTIONNE PAS avec dt = dtIm et nt = ntIm .  ****        
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
            
            assert tdata[0][-1]  <= tdata[0][0]    #  
            
            assert tdata[0][-1]  < tdata[0][0]    # Check that time is reversed in the chunk  
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
class realTimePlot():       # WORKS   
    """ draw the animation of a signal point by point (and not chunk by chunk)    """
    def __init__(self, maxt=100, dt=0.1, ymin = -1, ymax = 1):  
        nt=1
        super().__init__(nt, maxt, dt, ymin,ymax)
        
  
### ================================================================================
class realTimePlot_fct( realTimePlot):      # WORKS  
        
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

class realTimePlotSignal( baseRealTimePlot): #  WORKS for (nt>= 1 ) when Signal is properly provided
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

class realTimePlotSignal_EDF( realTimePlotSignal):        # WORKS WELL WHEN ONLY A SINGLE CHANNEL, DOES NOT REALLY SUPPORT MULTIPLE CHANNELS YET. 
    """
    """

    def __init__(self, nt=50, maxt=5000, iChan =0,fileEDF = 'Subject00_1.edf', datapath = 'C:/Users/Moi2/Documents/DONNEES/EEG_physionet/EEG_arithmeticTasks/'):
       
        ## EEG_EDFsignal instance : An attribute of the class realTimePlot_EDF
        Signal = EEG_EDFsignal(iChan,fileEDF ,datapath )    
        super().__init__(Signal, nt, maxt )  
  



### ==============================================================================================
### =============================================================================
class RTspectgram(baseRealTimePlot):                # en construction **************************************
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
        '''
    
        ## Signal can be anything as long as it has tSignal and ySignal 1-D arrays as attributes. 
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
        
        self.image = self.ax.pcolorfast(t,f , Sxx )       # ***  f and t ARE already in increasing order. (never reversed), and NOT transposed
            
        #nD =8  
        #self.ax.set_xlim( 0,  nD*self.ntIm*dt)   # pour tester 
        self.ax.set_xlim( 0,  self.maxt )   #        
        
        ## the extent of the data we want to display (apriori not the same as self.fmin and self.fmax)
        self.ax.set_ylim(self.ymin, self.ymax)          
        
        ##  extend of actual data.  fmin = min(f) , and fmax = max(f)  
        extent = ( 0, self.ntIm*dt, self.fmin, self.fmax )  
        self.image.set_extent(extent) # extent is data axes (left, right, bottom, top) for making image plots    
        '''
        print(f.shape)   # (513,) OK
        print(t.shape)   # (2,)   OK
        print(Sxx.shape) # (513, 2)         
        '''
        ##  Initialization after computation of first chunk spectrogram , ny et ntIm are previously defined, but NOT ymax and ymin
        ## We reverse the chunks only after taking the spectrogram: 
        t2= t[::-1]
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
        """Axes.specgram(x, NFFT=None, Fs=None, Fc=None, detrend=None, window=None, noverlap=None, cmap=None, xextent=None, pad_to=None, sides=None, scale_by_freq=None, mode=None, scale=None, vmin=None, vmax=None, *, data=None, **kwargs)
        https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.specgram.html
        """        
        pass
    
            
    def sendDataToUpdate(self):   
        '''  blablabla    
        
        tSignalChunk :  time of original signal chunk 
        tImChunk     :  time of spectrogram image chunk
        '''
        
        

        spectrChunkEnd = self.tdata[0][0]
        tSignalEnd = self.Signal.tSignal[-1]
        
        notEnd = spectrChunkEnd < tSignalEnd
        

        print(spectrChunkEnd)  # 413
        print(tSignalEnd)  #  90999 ms
        print(notEnd)
        
        if spectrChunkEnd> 90402.0 :
            print('stop')        
        if notEnd==False: 
            print('stop')
            
        while notEnd:            
            index = int(  self.getRightChunkLeftPoint(self.tdata, self.ntIm)/self.dt ) 
            print(index) #  must advance by step of ntIm=2       ***********
 
            tSignalChunk, ySignalChunk= self.Signal.returnSignalChunk(index, self.nt, 'non-reversed')
            tImChunk, fImChunk, SChunk = self.computeSpectgram_scipy( tSignalChunk, ySignalChunk)
            
            ## "yield" returns a generator for a 3-tuple to make them an image object 
            yield tImChunk, fImChunk, SChunk 
    

    def update(self, data):         
        ##  data is get from sendDataToUpdate throught the FuncAnimation function  , data = tChunk, fChunk, SChunk
        ## The deque tdata contains maxt/dt points, that means maxt/(dt*nt) chunks of nt points    
        ## We assume   Chunks to be in reversed-time order 
        
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
            
  
            if SChunk.shape[0]!=self.ntIm  : 
                print(len(tChunk))     #  1  *********
                print(tChunk[-1])   # 150.0
                print(tChunk[0])    # 150.0
                print(self.tdata[0][0])
                print(self.Signal.tSignal[-1])
                
                
                
            assert SChunk.shape[1]==self.ny
            #------------
            
        ## check the time -reversed order. 
        if tChunk[-1]> tChunk[0]: 
            tChunk = tChunk[::-1]
            SChunk = SChunk[::-1, :]  

        ## In the case of the spectrogram : each tChunk time begins to 0, because the spectrogram does not take tChunk into account but only fs, i.e. ds.    
        ## self.tdata[0][0] is the max value of the previous tChunk. The current tChunk is added after it. 
        
    
        ##  WEIRD BUG ????  !!!!!        
        '''
        #tChunk +=self.tdata[0][0]    # BIZRRE SEMBLE MODIFIER AUSSI tdata et pas que tChunk !!! 
        tChunk +=caca                # CETTE OPERATION MODIFIE self.tdata[0] !!!!    ?????   MAIS PKOI !!???????
        #tChunk = tChunk + caca       #  OK  CTTE OPERATION NE MODIFIE PAS  self.tdata[0]   
        print(self.tdata[0])         # 
        caca = self.tdata[0][0]       # 
        print(type(caca))               #   'numpy.float64'
        print(self.tdata[0])         # [0.413 0.15 ]   OK
        print(caca)                 # 0.413   OK
        print(tChunk )               #   [0.413 0.15 ]
        '''        
        
        tChunk = tChunk + self.tdata[0][0]         # OK NE MODIFIE PAS tdata  Mais utiliser '+=' ne fonctionne pas correctement !!  
        self.tdata, self.Sdata = self.updateDeque(self.tdata, self.Sdata, [tChunk,SChunk] )   # ICI on ajuste deja les axes...  ********

        
        ## concatenation         
        tConcat = concatenateDeque(self.tdata, 1)
        SConcat = concatenateDeque(self.Sdata, self.ny)
        
        ## reversed wrt time : 
        tConcat = tConcat[::-1]
        SConcat = SConcat[::-1, :]

        ntIm = self.ntIm 
        dt = self.dt
        assert len(tConcat)%ntIm == 0
        
        '''# pour tester  --- 
        nD =  int(len(tConcat)/ntIm)
        dt = self.dt 
        self.ax.set_xlim( tConcat[0], tConcat[0]+ nD*ntIm*dt) 
        '''

        
        self.ax.set_xlim( tConcat[0], tConcat[0]+ self.maxt) 
        #self.ax.set_xlim( self.tdata[-1][-1], self.tdata[-1][-1] + self.maxt)   #  self.tdata[-1] : the most left chunk :
        #self.ax.set_xlim( self.tdata[-1][-1], self.tdata[-1][-1] + nD*ntIm*dt)  
        
        self.ax.set_ylim(self.ymin, self.ymax)  
        
        '''
        #---- 
        #  si tConcat[-1]  continue de grimper: 
        
        if len(tConcat) ==8: 
            self.tcmax = tConcat[-1]
        elif  len(tConcat) >8:    
            if tConcat[-1] > self.tcmax : 
                tcmax = tConcat[-1]
                print('tConcat[-1] max =%s'%tcmax)
        #----
        '''
        
        print(tConcat[0])  # 0.15 au lieu de 0 ..... bof   ***???
        print( tConcat[-1])   # 0.826
        print(tConcat[0]+ self.maxt)
        #print( tConcat)      # OK croissant , sans se repeter. 
        print( tConcat.shape) 
        print(SConcat.shape) # (4, 513)      

        
        #-------------------
        '''
        if nD ==10:
            plt.clf
            plt.pcolormesh(tConcat, self.fdata, SConcat.T)
            plt.show()
        '''
        #------------------
        


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

    def displayFigure(self,tSignal, ySignal,t, f, Sxx):
        '''
        Only for testing 
        We simply display signal and spectrogram for a static first chunk of signal data. 
         Time axis of signal and spectre are not properly aligned: very ugly but not a problem for a simple test 
        '''

        ymin =self.ymin # 0.5 # Hz
        ymax = self.ymax #50  # Hz 
        #fs = 1000/dt # Sampling frequency of the t time series. Defaults to 1.0.   [Hz]  --->  1kHz 
        
        plt.close()
        
        plt.figure(1)
        plt.subplot(211)
        plt.plot(tSignal, ySignal)
        ax=plt.subplot(212)
        plt.pcolormesh(t, f, Sxx, shading='gouraud' )
        plt.ylim(ymin,ymax )
        plt.show()

        
            
### =========================================================================================
###  main 
### =========================================================================================
def returnOurSignal():
    '''return the channel #0 of an EEG extracted from an EDF file that has been picked on physionet.org 
       Return a Signal object, with tSignal and ySignal attributes.  
    
    '''
    
    fileEDF = 'Subject00_1.edf'
    datapath = 'C:/Users/Moi2/Documents/DONNEES/EEG_physionet/EEG_arithmeticTasks/'
    iChan =0
    Signal = EEG_EDFsignal(iChan,fileEDF, datapath )
    return Signal
    
def testMain():
    maxt = 6000 ; ntSignal = 600;  ny = 513
    Signal = returnOurSignal()
    spectro = RTspectgram(Signal, ntSignal, ny, maxt)
    spectro.test_update()
  
def plotEDFSignal_main():
    '''  Display an EDF signal in real-time, 
         Use double queue data structures '''
    maxt = 6000 # ms en principe,
    nt = 600 ; iChan =0
    
    RTplot = realTimePlotSignal_EDF( nt, maxt, iChan,fileEDF, datapath)
    RTplot.showAnimation()  

    
def main():
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


if __name__ == '__main__' :
    main()
    
'''
## Affichage RT du spectrogramme  
maxt = 6000 # ms en principe,
nt = 600;  Nwin=2 # 3  , # acceptable je pense.  Mais il faudrait refondre l'affichage apres chaque iteration
Prop = 0.5 

Signal = EEG_EDFsignal(iChan,fileEDF, datapath )

spectro = RTspectgram(Signal, nt, maxt)

## Affichage du spectrogramme de maniere statique 
maxt = 6000 # ms en principe,
nt = 600;  Nwin=2 # 3  , # acceptable je pense.  Mais il faudrait refondre l'affichage apres chaque iteration
Prop = 0.5 

Signal = EEG_EDFsignal(iChan,fileEDF, datapath )

spectro = RTspectgram(Signal, nt, maxt)
tSignalChunk, ySignalChunk= Signal.returnSignalChunk(index, nt, 'non-reversed')

spectro.computeSpectgram(tSignalChunk, ySignalChunk, Nwin, Prop)
spectro.displayFigure(tSignalChunk, ySignalChunk,t, f, Sxx)

'''


##=======-----  sinus:  for testing --------------------------------------------
'''
dt = 0.1
maxt = 4*pi
scopeSinus = realTimePlot_fct(maxt, dt)
scopeSinus.showAnimation()
'''

