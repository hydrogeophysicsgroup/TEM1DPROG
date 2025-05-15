C--------------------------------------------------------------------------
C --- THE IP COMMON BLOCK
C--------------------------------------------------------------------------
      INTEGER*4 IMODIP
      REAL*8 CHAIP(N_ARR),TAUIP(N_ARR),POWIP(N_ARR)
      COMMON /IPBL/ IMODIP, CHAIP, TAUIP, POWIP
C ----------
C   IMODIP:  IMODIP = [ 1 | 0 ] : [ Model IP effect | Do not model IP effect]
C   CHAIP:   The chargeability of the Cole-Cole model).
C   TAUIP:   The time constant of the Cole-Cole model.
C   POWIP:   The power of the Cole-Cole model.
C ----------
