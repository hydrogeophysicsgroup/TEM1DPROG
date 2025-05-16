CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     S U B R O U T I N E   I N T E R P 1 D
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Performs a local spline interpolation. The interpolated
C   curve passes through the data points, and the first
C   derivative is continuous.
C   However, the second derivative is not nessessarily continuous.
C
C   Parameters in the call:
C  =========================
C    NA  : The number of elements in the arrays XA and YA.
C    XA  : Tabscissas of the data set - must be in increasing order.
C    YA  : The ordinates of the data set.
C    NB  : The number of elements in the arrays XB and YB.
C    XB  : The arguments of the interpolated function.
C    YB  : The interpolated function values.
C    IOP : [0|1|2|3]: Interpolated function value is: [0th|1st|2nd|3rd] derivative.
C    EXTVAL : If XB is outside the sample interval, YB i set to EXTVAL.
C
C   The principle of the interpolation is that a 3rd degree
C   polynomium is fitted between the data points by assuming
C   the polynomium to pass through the data points, and that
C   the slope of the curve in a data point is equal to the
C   slope of the straight line connecting the two adjacent
C   data points. At the two end points, the slope is assigned
C   assuming a 2nd degree polynomial (see below).
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE INTERP1D (NA,XA,YA,NB,XB,YB,IOP,EXTVAL)
C ----------
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C--------------------------------------------------------------------------
      INCLUDE 'ARRAYSDIMBL.INC'
C ----------
      INTEGER*4 NA,NB,IOP,IINT(N_ARR)
      REAL*8 XA(NA),YA(NA),XB(NB),YB(N_ARR)
      REAL*8 EXTVAL
      REAL*8 DX(N_ARR),DY(N_ARR),SLOPE(N_ARR)
C ----------

      IWRI = 0

C.....................................................................
      IF (IWRI.EQ.1) THEN
        WRITE (*,*) 'INTERP1D: ENTER SUBROUTINE'
        WRITE (*,*) 'NA,XA,YA,NB,XB,IOP'
        WRITE (*,*)  NA,(XA(J),J = 1,NA)
        WRITE (*,*)  NA,(YA(J),J = 1,NA)
        WRITE (*,*)  (XB(J),J = 1,NB)
        WRITE (*,*)  IOP
        WRITE (*,*)  EXTVAL
      ENDIF
C.....................................................................

C----------------------------------------------------------------------------
C --- FIND X- AND Y-DIFFERENCES AND SLOPES
C----------------------------------------------------------------------------
      DO I = 1,NA-1
        DX(I) = XA(I+1)-XA(I)
        DY(I) = YA(I+1)-YA(I)
      ENDDO      

      DO I = 2,NA-1
        SLOPE(I) = (YA(I+1)-YA(I-1))/(XA(I+1)-XA(I-1))
      ENDDO
      SLOPE(1) = DY(1)/DX(1)
      SLOPE(NA) = DY(NA-1)/DX(NA-1)

C----------------------------------------------------------------------------
C --- Improvement on the slope at the end points.
C --- The slope is now the slope of the 2nd degree polynomial
C --- determined by the points [XA(1),YA(1)],[XA(2),YA(2)],
C --- and the slope at XA(2) = SLOPE(2),
C --- and similarly at the right end point.
C----------------------------------------------------------------------------
      SLOPE(1)  = 2.D0*SLOPE(1)-SLOPE(2)
      SLOPE(NA) = 2.D0*SLOPE(NA)-SLOPE(NA-1)
      
C-----------------------------------------------------------------------------------
C --- FIND WHICH INTERVAL EACH OF THE INTERPOLATION POINT ABSCISSAS BELONGS TO
C-----------------------------------------------------------------------------------
      DO I = 1,NB
      IINT(I) = 0

      DO J = 1,NA-1
      HH = (XB(I)-XA(J))*(XB(I)-XA(J+1))
      IF (HH.LE.0) THEN
        IINT(I) = J
        GOTO 19

        ELSEIF (XB(I).LT.XA(1)) THEN

        IINT(I) = 0

        ELSEIF (XB(I).GT.XA(NA)) THEN

        IINT(I) = 0

      ENDIF
      ENDDO

   19 CONTINUE
      ENDDO

C.....................................................................
      IF (IWRI.EQ.1) THEN
      WRITE (*,*) NB,(IINT(J),J = 1,NB)
      ENDIF
C.....................................................................
      
C=================================================================
C --- CALCULATE VALUES OF THE INTERPOLATED FUNCTION
C=================================================================
      DO I = 1,NB

C------------------------------------------------------------------------
C --- EXTRAPOLATION IF XB IS OUTSIDE THE SAMPLE INTERVAL
C------------------------------------------------------------------------
      IF (IINT(I).EQ.0) THEN

      IF (IOP.EQ.0) THEN
        YB(I) = EXTVAL
      ELSEIF (IOP.EQ.1) THEN
        YB(I) = 0.D0
      ELSEIF (IOP.EQ.2) THEN
        YB(I) = 0.D0
      ELSEIF (IOP.EQ.3) THEN
        YB(I) = 0.D0
      ENDIF

      ELSE

      II = IINT(I)
      XD = XB(I)-XA(II)
      DXI =    DX(II)
      DYI =    DY(II)
      YML =   SLOPE(II)
      YMR =   SLOPE(II+1)
      DX2I = 1.D0/(DXI*DXI)
      DX3I = DX2I/DXI


      IF (IOP.EQ.0) THEN
       Y0 = (3.D0*DYI-(YMR+2.D0*YML)*DXI)*DX2I
     #       +XD*((YML+YMR)*DXI-2.D0*DYI)*DX3I
       Y0 = YA(II)+XD*(YML+XD*Y0)
       YB(I) = Y0

      ELSEIF (IOP.EQ.1) THEN
       Y1 = YML+2.D0*XD*(3.D0*DYI-(YMR+2.D0*YML)*DXI)*DX2I
     #      +3.D0*XD*XD*((YML+YMR)*DXI-2.D0*DYI)*DX3I
       YB(I) = Y1

      ELSEIF (IOP.EQ.2) THEN
       Y2 = 2.D0*(3.D0*DYI-(YMR+2.D0*YML)*DXI)*DX2I
     #      +6.D0*XD*((YML+YMR)*DXI-2.D0*DYI)*DX3I
       YB(I) = Y2

      ELSEIF (IOP.EQ.3) THEN
       Y3 = 6.D0*((YML+YMR)*DXI-2.D0*DYI)*DX3I
       YB(I) = Y3
      ENDIF

      ENDIF
C --- ENDIF: EXTRAPOLATION OR INTERPOLATION

      ENDDO
C --- ENDD0: LOOP OVER ALL INTERPOLATION POINTS

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     S U B R O U T I N E   I N T E R P 1 D I N T
C
C     29.01.2024 / NBC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Performs an integration over a function given by its sample points.
C   The integral is calculated from the first sample to the last sample.
C   The function integrated is supposed to be given by a local spline
C   interpolation between the input samples, i.e. as a 3rd degree polynomial
C   in each interval.
C
C   Parameters in the call:
C    NA  : The number of elements in the arrays XA and YA.
C    XA  : Tabscissas of the data set - must be in increasing order.
C    YA  : The ordinates of the data set.
C    YINT: The integral from XA(1) to XA(NA).
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE INTERP1DINT(NA,XA,YA,YINT)
C ----------
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C--------------------------------------------------------------------------
      INCLUDE 'ARRAYSDIMBL.INC'
C ----------
      INTEGER*4 NA
      REAL*8 XA(NA),YA(NA),YINT
      REAL*8 DX(N_ARR),DY(N_ARR),SLOPE(N_ARR)
C ----------

      IWRI = 0

C.....................................................................
      IF (IWRI.EQ.1) THEN
        WRITE (*,*) 'INTERP1DINT: ENTER SUBROUTINE'
        WRITE (*,*) 'NA,XA,YA'
        WRITE (*,*)  NA,(XA(J),J = 1,NA)
        WRITE (*,*)  NA,(YA(J),J = 1,NA)
      ENDIF
C.....................................................................

C----------------------------------------------------------------------------
C --- FIND X- AND Y-DIFFERENCES AND SLOPES
C----------------------------------------------------------------------------
      DO I = 1,NA-1
        DX(I) = XA(I+1)-XA(I)
        DY(I) = YA(I+1)-YA(I)
      ENDDO      

      DO I = 2,NA-1
        SLOPE(I) = (YA(I+1)-YA(I-1))/(XA(I+1)-XA(I-1))
      ENDDO
      SLOPE(1) = DY(1)/DX(1)
      SLOPE(NA) = DY(NA-1)/DX(NA-1)

C----------------------------------------------------------------------------
C --- Improvement on the slope at the end points.
C --- The slope is now the slope of the 2nd degree polynomial
C --- determined by the points [XA(1),YA(1)],[XA(2),YA(2)],
C --- and the slope at XA(2) = SLOPE(2),
C --- and similarly at the right end point.
C----------------------------------------------------------------------------
      SLOPE(1)  = 2.D0*SLOPE(1)-SLOPE(2)
      SLOPE(NA) = 2.D0*SLOPE(NA)-SLOPE(NA-1)

C=================================================================
C --- CALCULATE VALUES OF THE INTEGRAL IN EACH SUBINTERVAL
C=================================================================
      YINT = 0
      DO I = 1,NA-1

      DXI =    DX(I)
      DYI =    DY(I)
      YML =   SLOPE(I)
      YMR =   SLOPE(I+1)
      DX2I = 1.D0/(DXI*DXI)
      DX3I = 1.D0/(DXI*DXI*DXI)

      Y00 = (3.D0*DYI-(YMR+2.D0*YML)*DXI)*DX2I
      Y01 = ((YML+YMR)*DXI-2.D0*DYI)*DX3I

      YINTDEL = YA(I)*DXI + 0.5D0*YML*DXI*DXI
     #         + Y00*DXI*DXI*DXI/3.D0
     #         + 0.25D0*Y01*DXI*DXI*DXI*DXI

      YINT = YINT+YINTDEL      

      ENDDO
C --- ENDD0: LOOP OVER ALL INTERPOLATION POINTS

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   R E A L * 8   F U N C T I O N   E X P C
C
C   EXPC tests the argument to prevent under- and overflow.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL*8 FUNCTION EXPC(X)
C ----------
      REAL*8 X,XX,XMIN,XMAX
C ----------
      DATA XMIN,XMAX /-300.D0,300.D0/
C ----------

      XX = X
      XX = MIN(XX,XMAX)
      XX = MAX(XX,XMIN)

      EXPC = EXP(XX)

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   C O M P L E X * 1 6  F U N C T I O N   C E X P C
C
C   CEXPC tests the argument to prevent under- and overflow.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      COMPLEX*16 FUNCTION CEXPC(Z)
C ----------
      COMPLEX*16 Z,ZZ,JJ
      REAL*8 XMIN,XMAX,RR
C ----------
      DATA XMIN,XMAX /-300.D0,300.D0/
      DATA JJ /(0.D0,1.D0)/
C ----------
  
      RR = DREAL(Z)
      RR = MIN(RR,XMAX)
      RR = MAX(RR,XMIN)
      ZZ = RR+JJ*DIMAG(Z)

      CEXPC = EXP(ZZ)
  
      RETURN
      END 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     R E A L * 8   F U N C T I O N   J 1
C
C   Calculates the Bessel function J_1 using a series expansion
C   for small argumenst and the asymptotic formula for large argumnets.
C
C   Programmed after Abramowitz & Stegun:
C   Handbook of Mathematical Functions,
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL*8 FUNCTION J1(X)
C ----------
      IMPLICIT REAL*8 (A-H,O-Z)  
      IMPLICIT INTEGER*4 (I-N)
C ----------
      REAL*8 A(7),B(7),C(7)
      DATA A /+0.50000000D0,-0.56249985D0,+0.21093573D0,
     #        -0.03954289D0,+0.00443319D0,-0.00031761D0,
     #        +0.00001109D0/
      DATA B /+0.79788456D0,+0.00000156D0,+0.01659667D0,
     #        +0.00017105D0,-0.00249511D0,+0.00113653D0,
     #        -0.00020033D0/
      DATA C /-2.35619449D0,+0.12499612D0,+0.00005650D0,
     #        -0.00637879D0,+0.00074348D0,+0.00079824D0,
     #        -0.00029166D0/
C ----------

      XA = ABS(X)

      IF (XA.EQ.0.) THEN
      J1 = 0.D0
      RETURN

      ELSEIF (XA.LE.3.D0) THEN

      S = X/XA
      Y = XA/3.D0
      Y = Y*Y
      J1 = A(1)+Y*(A(2)+Y*(A(3)+Y*(A(4)+Y*(A(5)+Y*(A(6)+Y*A(7))))))
      J1 = X*J1
      RETURN

      ELSE

      S = X/XA
      Y = 3.D0/XA
      BJ = B(1)+Y*(B(2)+Y*(B(3)+Y*(B(4)+Y*(B(5)+Y*(B(6)+Y*B(7))))))
      CJ = XA+C(1)+Y*(C(2)+Y*(C(3)+Y*(C(4)+Y*(C(5)+Y*(C(6)+Y*C(7))))))
      J1 = S*BJ*COS(CJ)/SQRT(XA)
      RETURN

      ENDIF

      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   R E A L * 8   F U N C T I O N   Z E R O P O S
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   For a given area of a circular loop and a given height above the plane
C   of the loop, REAL*8 FUNCTION ZEROPOS finds the lateral position
C   where the vertical magnetic field is zero, i.e. where a vertical receiver
C   dipole would be zero-coupled to a circualar transmitter loop.
C
C    INPUT
C   =======
C   AREA:   The area of the circular loop
C   HEIGHT: The vertical height of the receiver dipole above the loop (HRx-HTx).
C
C   REFERENCE
C  ===========
C   The London, Edinburgh, and Dublin Philosophical Magazine and Journal of Science,
C   Vol XLL - Sixth Series, January-June 1921: XXXIII
C   Magnetic field of circular currents, by H. Nagaoka, Professor of Physics,
C   Imperial University, Tokyo.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL*8 FUNCTION ZEROPOS(AREA,HEIGHT)
C ----------
      IMPLICIT INTEGER*4 (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
C ----------
      REAL*8 R(1001),HZ(1001)
C ----------
      REAL*8 ELLK,ELLE
      EXTERNAL ELLK,ELLE
C ----------
      REAL*8 PI
      DATA PI /3.14159265358979D+0/
C ----------

C--------------------------------
C --- Radius of loop
C--------------------------------
      RAD = SQRT(AREA/PI)

C--------------------------------
C --- Height of observation
C--------------------------------
      Z = HEIGHT

C--------------------------------
C --- Test radii of observation
C--------------------------------
      DR = 0.001D0*Z
      NR = 1001
      R(1) = RAD
      DO I = 2,NR
        R(I) = R(I-1) + DR
      ENDDO

      DO I = 1,NR

      AP  = RAD + R(I)
      AM  = RAD - R(I)
      RM2 = AM*AM + Z*Z
      RP2 = AP*AP + Z*Z
      RP  = SQRT(RP2)
      ARG = 4.D0*RAD*R(I)/RP2

      ELE = ELLE(ARG)
      ELK = ELLK(ARG) 

C-----------------------------------------------------------
C --- THE RADIAL B-FIELD
C --- Turned off for now. Not used to find zero position.
C-----------------------------------------------------------
C      IF (R(I).EQ.0.D0) THEN
C        HR(I) = 0.D0
C      ELSE
C        HR(I) = 2*( 2*RAD*Z*ELE/RM2 - Z*EKE/R ) / RP
C        HR(I) = 0.25D0*HR(I)/PI
C      ENDIF
C-----------------------------------------------------------

C-----------------------------------------------------------
C --- THE VERTICAL B-FIELD
C-----------------------------------------------------------
      HZ(I) = 2.D0 * ( 2.D0*RAD*AM*ELE/RM2 + (ELK-ELE) )/RP
      HZ(I) = 0.25D0*HZ(I)/PI

C-----------------------------------------------------------
C --- THE AMPLITUDE OF THE B-FIELD
C --- Turned off for now. Not used to find zero position
C-----------------------------------------------------------
C      HA(I) = SQRT(HZ(I)*HZ(I)+HR(I)*HR(I))
C-----------------------------------------------------------

      ENDDO
C --- ENDDO: Loop over the test points

C-----------------------------------------------------------
C --- INTERPOLATION TO FIND ZERO CROSSING
C --- Turned off for now.
C --- Using a simple linear interpolation. Sampling is dense.
C-----------------------------------------------------------
C      CALL INTERP1D(1001,HZ,R, 1,0.D0,ZEROPOS, 0,0.D0)
C-----------------------------------------------------------

C-----------------------------------------------------------
C --- Simple linear interpolation, testing from one end
C-----------------------------------------------------------
      DO I = 1,NR-1
      IF ((HZ(I)*HZ(I+1)).LE.0.D0) THEN
        ZEROPOS = ( HZ(I+1)*R(I) - HZ(I)*R(I+1) ) / ( HZ(I+1)-HZ(I) )
        GOTO 9999
      ENDIF
      ENDDO

 9999 CONTINUE

      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   R E A L * 8   F U N C T I O N   E L L K
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   REAL*8 FUNCTION ELLK is the complete elliptic integral K(X).
C
C   Programmed after Abramowitz & Stegun:
C   Handbook of Mathematical Functions, 2nd edition, p165, formula 17.3.34.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL*8 FUNCTION ELLK(X)
C ----------
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C ----------
      REAL*8 A(0:4),B(0:4)
      DATA A /1.38629436112D0,0.09666344259D0,0.03590092383D0,
     #        0.03742563713D0,0.01451196212D0/
      DATA B /0.50000000000D0,0.12498593597D0,0.06880248576D0,
     #        0.03328355346D0,0.00441787012D0/
C ----------

      Y = 1.D0-X
      H1 = A(0)+Y*(A(1)+Y*(A(2)+Y*(A(3)+Y*A(4))))
      H2 = B(0)+Y*(B(1)+Y*(B(2)+Y*(B(3)+Y*B(4))))
      ELLK = H1 - H2*LOG(Y)

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   R E A L * 8   F U N C T I O N   E L L E
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   REAL*8 FUNCTION ELLE is the complete elliptic integral E(X)
C
C   Programmed after Abramowitz & Stegun:
C   Handbook of Mathematical Functions, 2nd edition, p165, formula 17.3.36.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL*8 FUNCTION ELLE(X)
C ----------
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C ----------
      REAL*8 A(4),B(4)
      DATA A /0.44325141463D0,0.06260601220D0,0.04757383546D0,
     #        0.01736506451D0/
      DATA B /0.24998368310D0,0.09200180037D0,0.04069697526D0,
     #        0.00526449639D0/
C ----------

      Y = 1.D0-X
      H1 = 1.D0+Y*(A(1)+Y*(A(2)+Y*(A(3)+Y*A(4))))
      H2 = Y*(B(1)+Y*(B(2)+Y*(B(3)+Y*B(4))))
      ELLE = H1 - H2*LOG(Y)

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     R E A L * 8   F U N C T I O N   E R F C
C
C   REAL*8 FUNCTION ERFC calculates the real complementary error
C   function using a continued fraction formula for the argument
C   greater than 1 and using the formula ERFC = 1-ERF for
C   arguments smaller than 1.
C
C   Programmed after Abramowitz & Stegun:
C   Handbook of Mathematical Functions, formula 7.1.14, using the
C   modified Lenz's method from Numerical Recipes in FORTRAN 2nd edition,
C   page 165, formulas 5.2.8 to 5.2.10.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL*8 FUNCTION ERFC(XX)
C ----------
      IMPLICIT REAL*8 (A-H,O-Z)  
      IMPLICIT INTEGER*4 (I-N)
C ----------
      REAL*8 X, ERF, EXPC
      EXTERNAL ERF, EXPC
C ----------
      REAL*8 SQPI,XLIM,TINY,EPS
      DATA SQPI,XLIM / 0.5641895835477563D0, 2.D0/
      DATA TINY,EPS  /1.D-50, 1.D-12/
C ----------

      X = ABS(XX)

      IF (X .LT. XLIM) THEN
        ERFC = 1.D0 - ERF(X)
      ELSEIF (X .GT. 10.D0) THEN
        ERFC = 0.D0
      ELSE

      F = TINY
      CO = F
      DO = 0.D0
      A = 1.D0

      DN = X+A*DO
      IF (DN.EQ.0.D0) DN = TINY
      CN = X+A/CO
      IF (CN.EQ.0.D0) CN = TINY
      DN = 1.D0/DN
      DEL = CN*DN
      F = F*DEL
      D0 = DN
      C0 = CN

      A = 0.D0
      DO I = 1,300
      A = A + 0.5D0
      DN = X + A*D0
      IF (DN.EQ.0.D0) DN = TINY
      CN = X + A/C0
      IF (CN.EQ.0.D0) CN = TINY
      DN = 1.D0/DN
      DEL = CN*DN
      F = F*DEL
      IF (ABS(DEL-1.D0).LT.EPS) GOTO 99
      D0 = DN
      C0 = CN
      ENDDO

C.................................................................
      WRITE (*,1001) X
 1001 FORMAT (' More than 300 iterations in ERFC: argument ',
     #          1PE12.5,' too small')
C.................................................................

   99 ERFC = F*SQPI*EXPC(-X*X)

      ENDIF

      IF (XX.LT.0.D0) ERFC = 2.D0 - ERFC

      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     R E A L * 8   F U N C T I O N   E R F
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   REAL*8 FUNCTION ERF cumputes the real error function using
C   a power series for arguments smaller than 1 and the
C   expression ERF = (1 - ERFC) for arguments greater than 1.
C
C   Programmed after Abramowitz & Stegun:
C   Handbook of Mathematical Functions, formula 7.1.6.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL*8 FUNCTION ERF(X)
C ----------
      IMPLICIT REAL*8 (A-H,O-Z)  
      IMPLICIT INTEGER*4 (I-N)
C ----------
      REAL*8 X,ERFC, EXPC
      EXTERNAL ERFC, EXPC
C ----------
      REAL*8 SQPI,YLIM,EPS
      DATA SQPI,YLIM,EPS /0.5641895835477563D0, 2.D0, 1.D-12/
C ----------

      Y = ABS(X)
      Y2 = Y*Y

      IF (Y.EQ.0.D0) THEN
        ERF = 0.D0
      ELSEIF (Y.LE.YLIM) THEN

      S = Y
      SDEL = S
      I = 0
   10 I = I+1
      SDEL = SDEL*Y2*2.D0 / (2*I+1)
      S = S+SDEL
      C = ABS(SDEL/S)
      IF (C.GT.EPS) GOTO 10
      ERF = S*2.0*SQPI*EXPC(-Y2)

      ELSE
        ERF = 1.D0 - ERFC(Y)
      ENDIF

      IF (X.LT.0.D0) ERF = -ERF

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     R E A L * 8   F U N C T I O N   D I S T M I N
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   REAL*8 FUNCTION DISTMIN calculates the shortest distance
c   between a point (XP,YP) and a line segment given by
C   its end point coordinates: (X1,Y1) and (X2,Y2).
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL*8 FUNCTION DISTMIN(X1,Y1,X2,Y2,XP,YP)
C ----------
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C ----------

      DX = X2-X1
      DY = Y2-Y1

      IF (DX.EQ.0.D0 .AND. DY.EQ.0.D0) THEN
        DISTMIN = SQRT( (XP-X1)**2 + (YP-Y1)**2 )
      RETURN
      ENDIF

C-----------------------------------------------------
C --- CALCULATE THE T THAT MINIMIZES THE DISTANCE.
C-----------------------------------------------------
      T = ((XP-X1)*DX + (YP-Y1)*DY) / (DX*DX + DY*DY)

C----------------------------------------------------------
C --- DISCERN THE CASES: P INSIDE OR OUTSIDE LINE SEGMENT
C----------------------------------------------------------
      IF (T.LT.0.D0) THEN
      DISTMIN = SQRT( (XP-X1)**2 + (YP-Y1)**2 )

      ELSEIF (T.GT.1.D0) THEN
      DISTMIN = SQRT( (XP-X2)**2 + (YP-Y2)**2 )

      ELSE

      D1 = (X1-XP)+T*DX
      D2 = (Y1-YP)+T*DY
      DISTMIN = SQRT( D1*D1 + D2*D2 )

      ENDIF

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      F U N C T I O N    S T I M E R
C
C     Function STIMER returns the monitor time value in seconds.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL*8 FUNCTION STIMER(DUMMY)
C ----------
      INTEGER*2 IH,IM,IS,IS100
C ----------

      CALL GETTIM(IH,IM,IS,IS100)
      STIMER=IH*3600.D0+IM*60.D0+IS+IS100/100.D0

      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      S U B R O U T I N E    N U T I D
C
C     Function NUTID returns the monitor date and time.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE NUTID(ST)
C ----------
      CHARACTER*20 ST
C ----------
      INTEGER*2 IH,IM,IS,IS100,IY,IMON,IDAY
C ----------

      ST='                    '

      CALL GETTIM(IH,IM,IS,IS100)
      CALL GETDAT(IY,IMON,IDAY)
      WRITE (ST(1:19),1001) IDAY,IMON,IY,IH,IM,IS
 1001 FORMAT(I2,'.',I2,'.',I4,'/',I2,':',I2,':',I2)

      DO I=1,19
      IF (ST(I:I).EQ.' ') ST(I:I)='0'
      ENDDO

      RETURN
      END
