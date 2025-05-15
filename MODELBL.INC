C--------------------------------------------------------------------------
C --- THE MODEL COMMON BLOCK
C--------------------------------------------------------------------------
      INTEGER*4 IMLM,NLAY
      REAL*8 RHON(N_ARR),DEPN(N_ARR),SIGN(N_ARR),THKN(N_ARR)
      COMMON /MODELBL/ IMLM, NLAY, RHON, DEPN, SIGN, THKN
C ----------
C   IMLM:   IMLM = [ 0 | 1 ] means [ Few-layer model | Multi-layer model]
C   NLAY:   Number of layers in the model.
C   RHON:   Layer resistivities.
C   SIGN:   Layer conductivities.
C   THKN:   Layer thicknesses.
C   DEPN:   Depths to top of layer boundaries, DEPN(1)=0.
C ----------
