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
        
        self.signalLabels = self.f.sendDataToUpdateLabels()
           
    
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
            
        if mode =='reversed':    
            if index > 0:
                tChunk = self.tSignal[index + nt -1 : index-1:-1 ] 
                yChunk = self.ySignal[index + nt -1 : index-1:-1 ]
            elif index==0: 
                tChunk = self.tSignal[index + nt -1 ::-1 ] 
                yChunk = self.ySignal[index + nt -1 ::-1 ]
                
        else:  
            print('mode error')
        return tChunk, yChunk  


### ================================================================================
###  A abstract class. Cannot be instanced directly. The sendDataToUpdate() method is abstract (empty) and is defined only in derived classes.      
class baseRealTimePlot(metaclass=ABCMeta):          #         FONCTIONNE BIEN  lorsque nt = firstData est un int ; Verifions aussi pour firstData un 2-tuple. ******
    """
    An abstract class.  
    Draw the animation of a signal  chunk by chunk, for chunk of arbitrary size
    And for signal chunk OR image chunk 
    
    data = (tchunk, ychunk)  
    ychunk size :  (nt, ny)
    
    When we plot a signal, ny = 1, ychunk is an 1D array  of float or a list of float. 
    Whe we plot/show an image :  ny >1 , ychunk is a 2D array 
    
    In any cases: tchunk is an 1D array  of float or a list of float
    
    Size of ychunk and tchunk are set when receiving firstData = (tchunk, ychunk) that initialize the deques. 
    
    **Particular case:
    However, when the value of the argument "firstData" is an 'int' (not a float or anything else), it is interpreted as being ny instead of a proper first chunk. 
    And then :  ny is set tothis value while nt=1
    
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
        getTchunkYchunk( firstData)
        setChunkSize( A)
        getRightChunkLeftPoint()  
        update( data)
        showAnimation()                  
        sendDataToUpdate()  :   An abstract method required to send data to update method via animation.FuncAnimation()
        (abstract methods in python : Equivalent to C++ virtual functions )
        About abstract methods :   https://stackoverflow.com/questions/5856963/abstract-methods-in-python
        
    """
    def __init__(self, firstData = ([0],[0]), maxt=100, dt=0.1, ymin = -1, ymax = 1):  

        self.maxt = maxt 
        self.dt = dt
        self.ymin = ymin
        self.ymax = ymax

        ## set and test the chunk size and return the first chunk 
        ## first chunk : firstData = (tChunk, yChunk),  
        ## UNLESS the value of the param firstData is a scalar. In this case it is directly nt for the initial zero array
        tChunk, yChunk = self.setChunkSize(firstData)
        
        #--- pour tester sans utiliser setChunkSize
        tChunk= np.arange(self.nt)
        yChunk =np.zeros(self.ny)         
        
        #---
        
        ## we require the chunks to enter a entire number of time in the total width of the plot. 
        assert maxt/dt%self.nt ==0 
        
        ## Two double queues,one for each dimensions: x, y.
        self.tdata = deque()    
        self.tdata.appendleft(tChunk)
        
        if self.ny==1:  
            self.ydata = deque()          
            self.ydata.appendleft(yChunk)
        
        if self.ny>1:  
            
            self.fdata= yChunk #np.zeros(self.ny)
            
            ## SChunk: data image chunk to draw
            Schunk = zeros(self.ny, self.nt)   # ny = NFFT, since frequencies in y
            self.Sdata = deque()
            self.Sdata.appendleft(SChunk)
        
        ## figure 
        self.fig, self.ax = plt.subplots()   
        
        self.ax.set_ylim(self.ymin , self.ymax)
        self.ax.set_xlim(0, self.maxt)
        
        
        ## DEVRAIT ALLER DANS LES Classes derivees !! 
        #self.createPlot()
        ## when yChunk is 1-dim:  line object: 
        if self.ny==1:
            self.line = Line2D(self.tdata, self.ydata)
            self.ax.add_line(self.line)
        elif self.ny > 1:
            self.image = ax.pcolorfast(tChunk,  fChunk, SChunk)     #AxesImage object              

        

        

    def getTchunkYchunk(self, firstData) :
        ''' return tchunk , ychunk 
            Both firstData[0] and firstData[1] must be either lists or arrays, with nt >=1
            For arrays, ny >=1 
            For lists : ny=1.
            In case a first data chunk is available from initialization. 
            set self.nt and self.ny
        '''
        ## Both firstData[0] and firstData[1] must be either lists or arrays
        tchunk = firstData[0] 
        ychunk = firstData[1]
        
        ## length of chunks along time axis 
        assert len(self.tchunk)==1    # True even when tchunk = [0], a list or an array. --->  nt=1             
        self.nt = tchunk.shape[0] 
        nt_2 = ychunk.shape[0]
        assert self.nt == nt_2
        
        ## reverse tchunk and ychunk along the time axis 
        tchunk =tchunk[::-1]
        ychunk = ychunk[::-1]   # reverse along the time axis (i.e. permute rows)
        
        if len(ychunk.shape)==1:   
            self.ny = 1        
        elif len(ychunk.shape)==2:  
            self.ny = ychunk.shape[1]   
        else: print('format error')             # SHOULD raise an error !
        return tchunk, ychunk          

        
    def setChunkSize(self, A):
        '''
        chunk size (nt, ny ) are set internally as attributes (self.nt, self.ny)
        
        *** A := first chunk :  firstData = (tchunk0, ychunk0)
        
        *** Except in the particular case when the value of the argument "firstData" is an 'int' (not a float or anything else),
        it is interpreted as being ny instead of a proper first chunk.  And then :  ny is set to this value while nt=1
        
        Time is in reversed order (i.e. decreasing ) in tChunk and yChunk compared with the order in firstData (which is increasing)
        '''
        ## particular case where A is scalar:
        if type(A) == int: 
            self.ny = 1
            self.nt = A
            tchunk = list(reversed(np.arange(A))  )   # because matplotlib does not support generators 
            ychunk = [0.0]*self.nt                    # a list [0.0, 0.0, ..... 0.0], with nt times 0.0
            
        else: 
            ##  first chunk :  firstData = (tchunk, ychunk)
            firstData =A # change the name only for clarity, but useless otherwise..
            
            if type(firstData)!=tuple:
                print("format error")                 # should raise an error
            else:
                ## Extract tchunk and ychunk from firstData 
                tchunk, ychunk = getTchunkYchunk(firstData)
        return tchunk, ychunk          
      
    def getRightChunkLeftPoint(self):
        '''  return the most left point  in the last chunk (i.e. the most right ), in case nt==1  or nt > 1
        '''
        if self.ny==1:
            if self.nt>1:
                ## tdata[0] : the last  ('most right') chunk that have been pushed in the deque. It is a list or a 1D array, that contains nt points.
                assert self.tdata[0][-1]  < self.tdata[0][0]    # Check that time is reversed in the chunk  
                return self.tdata[0][0]
                
            if self.nt==1:
                ## tdata[0] : the last point in the deque, a scalar. 
                ## index of the NEXT point to pick to push it in the deque. 
                return self.tdata[0]  

    @abstractmethod 
    def update(self, data):        #**********************

        ##  data is get from sendDataToUpdate throught the FuncAnimation function  , data = (t,y)
        ## The deque tdata contains maxt/dt points, that means maxt/(dt*nt) chunks of nt points       

        """
        assert self.ny ==1    # VERSION FOR 1D single, NOT images 
        if len(data) == 2 :
            tChunk, yChunk = data            
        else: print('format error !!  ')
           
        if len(self.tdata ) * self.nt >= self.maxt/self.dt:
            self.tdata.pop()            
            assert len(self.ydata)==len(self.tdata ) 
            self.ydata.pop()
        
        self.tdata.appendleft(tChunk)   
        self.ydata.appendleft(yChunk)              


        ## maximum is taken over ALL the elements of the lists or array that are piped in the deque ydata             
        ymaxNew = np.max(self.ydata)   
        yminNew = np.min(self.ydata)            

        if ymaxNew  > self.ymax:  self.ymax = ymaxNew   
        if yminNew  < self.ymin:  self.ymin = yminNew   
        
        ## REM:  self.tdata[-1]  corresponds to the first data  pushed in the deque (a chunk or a single point).
        self.ax.set_xlim(self.tdata[-1], self.tdata[-1] + self.maxt)         
        self.ax.set_ylim(self.ymin, self.ymax)        

        self.line.set_data(self.tdata, self.ydata)
        return self.line   
        """



        
    @abstractmethod 
    def sendDataToUpdate(self):
        ''' To be override '''
        pass
    
    ## ON POURRAIT DEVOIR AJOUTER UN DECORATEUR ICI POUR ADAPTER CTTE FCT ??    
    def showAnimation(self):
        
        ## a generator is passed from "sendDataToUpdate" to the update function. 
        ani = animation.FuncAnimation(self.fig, self.update, self.sendDataToUpdate, interval=self.dt )
        plt.show()  
    
 
### ================================================================================

##===========================================================================
class realTimePlot(baseRealTimePlot):       # WORKS   
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
        
        ## self.dt is declared in super().__init__ , but initialize at an arbitrary value
        ## Assuming that the time is the same between each points.         
        ## initialization of dt here: 
        self.dt = tSignal[2]-tSignal[1]

    def update(self, data):
        ## override an abstract method in the basic abstract class 
        
        ##  data is get from sendDataToUpdate throught the FuncAnimation function  , data = (t,y)
        ## The deque tdata contains maxt/dt points, that means maxt/(dt*nt) chunks of nt points       
    
        assert self.ny == 1   # case for 1D signal 
        
        if len(data) == 2 :
            tChunk, yChunk = data            
        else: print('format error !!  ')
           
        if len(self.tdata ) * self.nt >= self.maxt/self.dt:
            self.tdata.pop()            
            assert len(self.ydata)==len(self.tdata ) 
            self.ydata.pop()
        
        self.tdata.appendleft(tChunk)   
        self.ydata.appendleft(yChunk)              
    
    
        ## maximum is taken over ALL the elements of the lists or array that are piped in the deque ydata             
        ymaxNew = np.max(self.ydata)   
        yminNew = np.min(self.ydata)            
    
        if ymaxNew  > self.ymax:  self.ymax = ymaxNew   
        if yminNew  < self.ymin:  self.ymin = yminNew   
        
        ## REM:  self.tdata[-1]  corresponds to the first data  pushed in the deque (a chunk or a single point).
        self.ax.set_xlim(self.tdata[-1], self.tdata[-1] + self.maxt)         
        self.ax.set_ylim(self.ymin, self.ymax)        
    
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
            index = int(  self.getRightChunkLeftPoint()/self.dt ) 
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
    
    
    
    """
  
    def __init__(self,Signal, nt=50, maxt=6000):
        
        assert hasattr(Signal, 'tSignal') and hasattr(Signal, 'ySignal')
        
        self.Signal = Signal
        tSignal =Signal.tSignal  
        #ySignal =Signal.ySignal
        

        ## baseRealTimePlot initializer : initialize: 
        super().__init__(nt, maxt )   
        
     
         
    

    def computeSpectgram_scipy(self, tSignalChunk, ySignalChunk, Nwin, PropOverlap=0):
        """
        Using the scipy.signal library :
            from scipy.signal import spectrogram 
        Alternative possibility: 
            from matplotlib.Axes import specgram
        """
        SCALING = 'spectrum' #scaling : { 'density', 'spectrum' }
        
        dt = abs(tSignalChunk[2] - tSignalChunk[1]  )            # en supposant que c'est le temps en ms
        
        fmin = 0.5 # Hz
        fmax = 50  # Hz 
        fs = 1000/dt # Sampling frequency of the t time series. Defaults to 1.0.   [Hz]  --->  1kHz
    
        nt = len(tSignalChunk)  #   chunk length 
        print('nt=%s'%nt)
        PropOverlap = 0 #0.5
        
        print('Nwin=%s'%Nwin)
        
        NPERSEG= int(len(tSignalChunk)/Nwin )    #  number per segment  20#
        print('NPERSEG=%s'%NPERSEG) 
        
        exponent = int(ceil(log2(nt)) )
        
        NFFT = 2**exponent  # with zero padding         # number of samples in the FFT,i.e. in frequency.  --->  NFFT/2 +1 
        
        window_1=('tukey', PropOverlap)   # i.e. Tukey window with 1/4 of a window’s length overlap at each end, i.e. overall : 1/2 overlapped
        #  'when array-like : window must be 1-D 
        
        #  ( NE MARCHE PAS COMME CELA :   window_2 = 20e-3 # window width in time  [s]   20 ms ---> 1/20ms = 1/20 kHz = 0.05 kHz  )
        #framelen=int(window_2*fs)   #  frame length  ; --->  20 pts
        
    
        #expZeroPadding = 10 # default is 8, since 2**8=256
        #NFFT=2**expZeroPadding# Length of the FFT used, if a zero padded FFT is desired. If None, the FFT length is nperseg. otherwise: =256
        
        #NFFT = framelen
        #NOVERLAP=framelen/2  noverlap=NOVERLAP
       
        #f, t, Sxx  = spectrogram(ySignalChunk, fs) 
        ##  REM  tSignalChunk is not directly involved itself in spectrogram function. 
        f, t, Sxx  = spectrogram(ySignalChunk, fs, window=window_1, nperseg=NPERSEG,nfft = NFFT, detrend='constant', return_onesided=True, scaling=SCALING, axis=-1, mode='psd')
        
        #assert t==tSignalChunk 
        #print(t)   # entre 0.001 et 0.049    # tSignalChunk /fs   #  [s]
        #print(tSignalChunk)  # entre 1 et 50  :                    [ms]   Mais inverse 
        
        
        #print(Sxx.shape)   # (129, 1)    
        #print(Sxx[0,0])
        #print(ySignalChunk.shape)
        #print(t.shape)
        #print(f.shape)  
        #print(t[0:5])
        #print(tSignalChunk)
        
        #self.displayFigure(tSignalChunk, ySignalChunk,t, f, Sxx)
        
        return t, f, Sxx
        
        
    def computeSpectgram_mpl(self, tSignalChunk, ySignalChunk, Nwin, PropOverlap=0):
        """Axes.specgram(x, NFFT=None, Fs=None, Fc=None, detrend=None, window=None, noverlap=None, cmap=None, xextent=None, pad_to=None, sides=None, scale_by_freq=None, mode=None, scale=None, vmin=None, vmax=None, *, data=None, **kwargs)
        https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.specgram.html
        """        
        pass
    
    def displayFigure(self,tSignal, ySignal,t, f, Sxx):
        
        ## On affiche ici avec les signaux "a l'endroit" . avec plt.pcolormesh 
        fmin = 0.5 # Hz
        fmax = 50  # Hz 
        #fs = 1000/dt # Sampling frequency of the t time series. Defaults to 1.0.   [Hz]  --->  1kHz 
        
        
        # time is not exactly aligned . 
        plt.figure(1)
        plt.subplot(211)
        plt.plot(tSignal, ySignal)
        ax=plt.subplot(212)
        plt.pcolormesh(t, f, Sxx, shading='gouraud' )
        
        
        #plt.pcolormesh(t, y, f, Sxx)
        #plt.ylabel('Frequency [Hz]')  
        #plt.xlabel('Time [s]')
        plt.ylim(fmin,fmax )
        #plt.ylim(0, 60) # 1 minute
        plt.show()
        
        #### 
        #fig, axs = plt.subplots(2, 1)
        #fig.suptitle('Multiple images')
        
        ###   If X and Y are each equidistant, imshow can be a faster alternative.
        #ax.imshow(Sxx, interpolation = 'sinc', cmap='viridis')  # N'importe quoi ... PAS CORRECT


    def sendDataToUpdate(self):   
        '''  blablabla      
        '''
        while True:            
            # 
            index = int(  self.getRightChunkLeftPoint()/self.dt ) 
        
            tSignalChunk, ySignalChunk= Signal.returnSignalChunk(index, nt, 'non-reversed')
        
            tChunk, fChunk, SChunk=computeSpectgram(self, tSignalChunk, ySignalChunk, Nwin, PropOverlap=0)
            
            ## We reverse the chunk only after taking the spectrogram: 
            yield tChunk[::-1], fChunk[::-1], SChunk[::-1, :] # "yield" returns a generator for a 3-tuple to make them an image object 
    

    def update(self, data):              # ********   EN CONSTRUCTION     ***********
        ##  data is get from sendDataToUpdate throught the FuncAnimation function  , data = tChunk, fChunk, SChunk
        ## The deque tdata contains maxt/dt points, that means maxt/(dt*nt) chunks of nt points    
        
        if(len(data) == 3) : 
            tChunk, fChunk, SChunk = data
        else: 
            print('format error !! ')  # should raise an error *****
            
        ## tChunk, fChunk, SChunk  are to be pushed into the  deques:  tdata, Sdata. 
        ##  fChunk should always be the same . It is store into fdata and check at each iteration*********
        ## For frequencies : fdata is actually NOT a deque, but an 1D array, and there is NO need for ydata here. 
        
        if len(self.tdata ) * self.nt >= self.maxt/self.dt:
            assert len(self.tdata ) == self.Sdata.shape[0]
            self.Sdata.pop()
            self.tdata.pop()            
        

        self.tdata.appendleft(tChunk)   
        self.fdata = fChunk    # should almost always stay the same.
        self.Sdata.appendleft(SChunk)              
        
        
        self.ax.set_xlim(self.tdata[-1][-1], self.tdata[-1][-1] + self.maxt)   #  self.tdata[-1] : the most left chunk :
        self.ax.set_ylim(self.ymin, self.ymax)        #  for fdata , not ydata, actually  
        self.image.set_data(self.tdata, self.fdata, self.Sdata)    
        return self.image   

        
        
"""
        assert fChunk == self.fdata 
        
        if len(self.tdata ) * self.nt >= self.maxt/self.dt:
            assert self.Sdata.shape[0] == len(self.tdata)
            self.tdata.pop()
            self.Sdata.pop()

        if(len(data) == 3) : 
            self.tdata.appendleft(tChunk)   
            self.Sdata.appendleft(SChunk)              
        else: 
            print('format error !! ')  # should raise an error *****
        
        '''
        ## maximum is taken over ALL the elements of the lists or array that are piped in the deque ydata             
        ## In cases ny=1, nt ==1 or nt> 1 :  np.max() applies     
        ymaxNew = np.max(self.ydata)   
        yminNew = np.min(self.ydata)            

        if ymaxNew  > self.ymax:  self.ymax = ymaxNew   

        if yminNew  < self.ymin:  self.ymin = yminNew   
        '''
        
        ## REM:  self.tdata[-1]  corresponds to the first data  pushed in the deque (a chunk or a single point).
        if self.nt==1:
            self.ax.set_xlim(self.tdata[-1], self.tdata[-1] + self.maxt)         
        elif self.nt>1 :
            self.ax.set_xlim(self.tdata[-1][-1], self.tdata[-1][-1] + self.maxt)   #  self.tdata[-1] : the most left chunk :
        self.ax.set_ylim(self.ymin, self.ymax)        
        
        
        ##  COMPARE  image :   https://matplotlib.org/api/image_api.html 
        ##      and  Lines2D   https://matplotlib.org/api/_as_gen/matplotlib.lines.Line2D.html      
        #self.line.set_data(self.tdata, self.ydata)
        #return self.line     #Line2D object
    
        self.image.set_data(self.tdata, self.fdata, self.Sdata)  
        return self.image   # image object
"""



            
            
### =========================================================================================
###  main 
### =========================================================================================


#---  EEG from edf file, one point at a time. ------------------------------
fileEDF = 'Subject00_1.edf'
datapath = 'C:/Users/Moi2/Documents/DONNEES/EEG_physionet/EEG_arithmeticTasks/'
iChan =0

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

##=======  affichage du signal en RT avec deque . -----------------------------------
RTplot = realTimePlotSignal_EDF( nt, maxt, iChan,fileEDF, datapath)
RTplot.showAnimation()  



##=======-----  sinus:  for testing --------------------------------------------
'''
dt = 0.1
maxt = 4*pi
scopeSinus = realTimePlot_fct(maxt, dt)
scopeSinus.showAnimation()
'''