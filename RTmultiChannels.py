#-*- encoding: utf-8 -*- 

"""
RTmultiChannels.py

Use multiprocessing.Process to start the analysis and display of each channel in a dedicated process.  

Then:  N channels ---> N processes. 

"""

from time import sleep
import multiprocessing as mp 
import copy

import numpy as np

## personal "libraries"
import realTimePlot_singleChannel as RT 
import signalClass as sC
## WE DONT REALLY USE the class named SignalClass for now, but rather the EEG_EDFsignal class that is also on the misnamed signalClass file. *** WE BE TO CHANGE WHEN ADATING SignalClass. 

    
## 
class RTmultiChannelsMonitor():
    def __init__(self, Nchannels = 2, file=None, path=None):    
   
        ## list of channels we want the signals
        channelList =range(Nchannels)
        
        ## if OS == windows : 
        mp.freeze_support()   
                
        ## In general, SignalChList is a list of 'Signal' objects. 
        ## Each of these signal objects for each channel; Each signal object has attributes tSignal and ySignal. 
        if file is not None:
            ## returnEDFSignalCh is in the file named signalClass.py but is actually not in any class in this file. 
            ## In this case, SignalChList is a list of EEG_EDFsignal objects. 
            SignalChList = sC.returnEDFSignalCh(channelList, file, path)
            
            channelNameList = [ s.selectedChannel for s in SignalChList ] 
            
        else: 
            ## Rather than loading signals from a file, we use pyLSL communication with openBCI-python
            #Signal = sC.returnMyPersonalEEG()  # the purpose, instead of the edf EEG we test with for now.      :  TO DO  ************** 
            pass
        
  
        ## Signalch1 and Signalch2  are (for now) of EEG_EDFsignal class, and should be of a more general class in the future
        ## in ALL CASES, Signalch1 and Signalch2 have tSignal and ySignal as attributes, and that is the only requirement. 
        
        ## Each process calls 'inEachProcess' method, which instantiates a "realTimePlot_singleChannel"  object that animates a spectro.
        processList = list()
        for n in channelList: 
            p = mp.Process(target= self.inEachProcess , name=channelNameList[n], args=(SignalChList[n],) )
            processList.append(p)
        print('---Now we start the processes.---')    
        #  ****  IN CASE OF EDF at least, the processes seem to be not exactly simultaneous.. The two figures are NOT simultaneously update during the animation (Maybe the problem is )
        for p in processList:    p.start()
        
        print('---join---')
        for p in processList:    p.join() 
 
        print('All processes are closed')
        
 


    def inEachProcess(self, SignalCh):
        
        maxt = 6000 ; ntSignal = 600;  ny = 513
     
        ## SignalCh: signal objects: each of them represents a single-channel signal
        ## SignalCh attributes:  has tSignal, ySignal 
        
      
        self.spectro = RT.RTspectgram(SignalCh, ntSignal, ny, maxt)
        
        #  self.spectro.fig ,  self.spectro.axis  
        self.spectro.showAnimation()  
    
        
   
if __name__ == '__main__':
    

    Nchannels=2
    EDFfile = 'Subject00_1.edf'
    path = 'C:/Users/Moi2/Documents/DONNEES/EEG_physionet/EEG_arithmeticTasks/'    
    multiChannels = RTmultiChannelsMonitor(Nchannels,EDFfile , path )

    

    