C------------------------------------------------------------------------
C --- THE WAVEFORM COMMON BLOCK
C------------------------------------------------------------------------
      INTEGER*4 NWAVE
      REAL*8 TWAVE(N_ARR),AWAVE(N_ARR),SLOPE(N_ARR),D2WDT2(N_ARR)
C ----------
      COMMON /WAVEBL/ NWAVE,TWAVE,AWAVE,D2WDT2
C ----------
C   NWAVE:             Number of samples in the piecewise linear waveform
C   TWAVE(1:NWAVE):    Sample times of the piecewise linear waveform.
C   AWAVE(1:NWAVE):    Sample values of the piecewise linear waveform.
C   SLOPE (1:NWAVE-1): The 1st derivative of the waveform.
C   D2WDT2(1:NWAVE):   The 2nd derivative of the waveform.
C ----------
