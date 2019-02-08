#-*- encoding: utf-8 -*-

"""
define a signal
Can return a signal chunk by chunk of an EEG signal that is extracted from an EDF file. 
Time: in ms

"""
# ****************      la DEFINITION DES SIGNAUX EDF et LSL est a faire d'une maniere coherente   ***********************************************
#                        qui permettent aux fonctions ...._EDF  et..._LSL de fonctionner de maniere apparemment similaire , en maximisant la reutilisation de code. 

import numpy as np
from numpy import pi, sin, zeros, ones


## ==========================================================================
###    signalClass class :  present what is common between diverse type of signals 
###    EEG_LSLsignal  and   EEG_EDFsignal  are subclasses 
###  attributes : 
###               tSignal ; ySignal ;   
## ==========================================================================
class signalClass():
    
    """
    Describes a basic signal that is independent of format or communication protocol  (ex:  EDF format and LSL protocol)
    In particular:  it wraps around pylsl, pyedf and other libraries if required. 
    
    Supports signal with several channels, and allows to independently work with each of them. 
    
    Attributes of signalClass: 
            nChannels
            Nsamples
            signalLabels
            tSignal 
            signalType :  a string; possible values are 
            for each channel:   an object of class singleChannel() 
                         
    """
    class singleChannel():    
        """
        Describre the signal for a single channel
        Cannot be instanced outside the signalClass class 

        Attributes: 
            ySignal
            iChan
        """
        def __init__(self, Nsamples,iChan =0, yChunk=[]):
            ''' singleChannel class initializer 
            If len(yChunk) is larger than Nsamples, then it is truncated. 
            
            '''
            try:  
                if yChunk==[]: yChunk = zeros(Nsamples)
                self.ySignal = yChunk[0:Nsamples-1]
            except: print('format error')    
            
            self.iChan =0

        
    def __init__(self,signalType, nChannels, Nsamples, maxt):
        '''  signalClass class initializer '''
        
        ## acceptable signal type list can be expanded if needed
        acceptableSignalTypes = ['EDF', 'BCI_LSL']
        assert signalType in acceptableSignalTypes
        
        ## initialization of a 'zero' signal, before importing any real data.         
        self.signalType = signalType
        self.nChannels =1                       # int
        self.Nsamples =1000                     # int
        self.signalLabels = ['zeroSignal']      # a list of string 
        self.maxt
        self.tSignal = np.arange(self.Nsamples) # an 1D array of length= Nsamples 
        self.s0 = singleChannel(self.Nsamples) 




## ================================================================================================================================
class EEG_LSLsignal(signalClass):                  #  EN CONSTRUCTION ********************   DOIT TENIR COMPTE DE signalClass !!! ***********
    """ 
    See also the personal mini-package of tools for pylsl: pylsl_persoUtils.py . 
    
    
    """
    
    def __init__(self):
        '''
        Must set 
        the attributes:
        tSignal and ySignal
        
        '''
        ## resolve an EEG stream from the network.   
        streams = resolve_stream('type', 'EEG')    
        self.inlet = StreamInlet(streams[0])

    

    def getLSLSignalChunk(self):  
        '''Resolve EEG stream with pyLSL  and return tChunk, yChunk in reversed order with respect to time axis'''

   
        # Where is the chunk size defined ?? *******  
  
        yChunk, tChunk = self.inlet.pull_chunk()
        if tChunk:
            return tChunk[::-1], yChunk[::-1]

    

##  ce qui precede:  A travailler.  
##================================================================================================
##================================================================================================
### EEG  signal from EDF file :        IS USED IN realTimePlot_singleChannel.py and is WORKING :  
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
