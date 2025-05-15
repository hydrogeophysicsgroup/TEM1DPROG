C--------------------------------------------------------------------------
C --- THE RESPONSE COMMON BLOCK
C--------------------------------------------------------------------------
      INTEGER*4 IRESPTYPE,IDERIV,IRXFILT,IAMPFILT,IREP,IWCONV
      REAL*8 REPFREQ, RXFREQ, AMPFREQ
C ----------
      COMMON /RESPBL/ IRESPTYPE,IDERIV,IRXFILT,IAMPFILT,IREP,IWCONV,
     # REPFREQ, RXFREQ, AMPFREQ
C ----------
C   IRESPTYPE:   [0 | 1 | 2] :
C               [Step w/o filters, no rep & waveform convolution |
C               Impulse w/o filters, no rep, no waveform convolution |
C               Convolved responses]
C   IDERIV:     [0 | 1] : [Do not calculate derivatives | Calculate derivatives ]
C   ICEN:       [0 | 1] : [Offset loop | Central loop] configuration
C               The value of ICEN is automatically determined in the program.
C   IRXFILT:    [0 | 1 | 2] : [Do not model Rx filter | Model Rx filter as 1st or 2nd order filter]
C   IAMPFILT:   [0 | 1 | 2] : [Do not model Amp filter | Model Amp filter as 1st or 2nd order filter]
C   IREP:       [1 | 0] : [Do | Do not] model repetition
C   IWCONV:     [1 | 0] : [Do | Do not] convolve with waveform.
C ----------
