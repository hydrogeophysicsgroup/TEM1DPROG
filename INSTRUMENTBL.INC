C-----------------------------------------------------------------------
C --- THE INSTRUMENT COMMON BLOCK
C-----------------------------------------------------------------------
      INTEGER*4 IZEROPOS,ICEN,ISHTX1,ISHTX2,ISHRX1,ISHRX2
      REAL*8 TXAREA,TXRAD,RTXRX,HTX1,HTX2,HRX1,HRX2,EQRXPOS
      COMMON /INSTRUMENTBL/ ICEN,IZEROPOS,ISHTX1,ISHTX2,ISHRX1,ISHRX2,
     #            TXAREA,TXRAD,RTXRX,HTX1,HTX2,HRX1,HRX2,EQRXPOS
C -----------
C   ICEN:     [0 | 1] : [Offset loop | Central loop] configuration
C   IZEROPOS: [1 | 0] : [Rx is zero coupled to Tx | Not zero coupled]
C   ISHTX1:   The polarity of Tx numer 1: [ -1 | 0 | 1 ].
C   ISHTX2:   The polarity of Tx numer 2: [ -1 | 0 | 1 ].
C   ISHRX1:   The polarity of Rx numer 1: [ -1 | 0 | 1 ].
C   ISHRX2:   The polarity of Rx numer 2: [ -1 | 0 | 1 ].
C   TXAREA:   Area of the Tx.
C   TXRAD:    Radius of a an equivalent circular Tx loop
C             (Calculated if IZEROPOS = 1).
C   RXSZ:     Horizontal distance between Tx centre and Rx dipole.
C   EQRXPOS:  Equivalent value of RTXRX for zero coupled Rx.
C ----------


