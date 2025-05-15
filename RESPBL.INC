C--------------------------------------------------------------------------
C --- THE RESPONSE COMMON BLOCK
C--------------------------------------------------------------------------
      INTEGER*4 IRESPTYPE,IDERIV,NFILT,IREP,IWCONV
      REAL*8    REPFREQ, FILTFREQ(16)
C ----------
      COMMON /RESPBL/
     #  IRESPTYPE,IDERIV,NFILT,IREP,IWCONV,REPFREQ,FILTFREQ
C ----------
C   IRESPTYPE:   [0 | 1 | 2] :
C               [Step w/o filters, no rep & waveform convolution |
C               Impulse w/o filters, no rep, no waveform convolution |
C               Convolved responses]
C   IDERIV:     [0 | 1] : [Do not calculate derivatives | Calculate derivatives ]
C   ICEN:       [0 | 1] : [Offset loop | Central loop] configuration
C               The value of ICEN is automatically determined in the program.
C   NFILT:      The number of 1st order filters to be applied.
C   FILTFREQ:   The cutoff frquencies of the filters.
C   IREP:       [1 | 0] : [Do | Do not] model repetition
C   IWCONV:     [1 | 0] : [Do | Do not] convolve with waveform.
C ----------
