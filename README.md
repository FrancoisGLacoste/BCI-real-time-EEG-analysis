# BCI-real-time-EEG-analysis


 EEG signal OR spectrogram shown in real-time via animation.

Show an animation chunk by chunk (and hopefully in real-time) of 
              1)   an EEG signal that is extracted from an EDF file. 
              2)   the spectrogram of the EEG, chunk by chunk in (hopefully) real-time   

FOR NOW:
THIS PROGRAM WORKS ONLY FOR EEG FILE in EDF FORMAT.  NOT WITH LSL yet. 
But the purpose in the near future is to make it work with signals acquired through pyLSL.

FOR NOW: THE REAL_TIME ANIMATION WORKS ONLY FOR SIGNAL WITH SINGLE CHANNEL. But it is to be extended to multichannel acquisition. 
FOR NOW: it only plots either the signal in real-time, OR the spectrogram. 

Data structure :  double queues (deque)
