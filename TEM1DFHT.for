CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     S U B R O U T I N E   F H T C O N V
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  Subroutine fhtconv calculates a J0 or J1 Hankel transform of a
C  REAL function. The integral is computed as a convolution between
C  sampled values of the kernel function and Fast Hankel Transform
C  filter coefficients. All calculations are in double precision
C  REAL*8.
C
C  What is calculated is the integral:
C  
C  G(r) = integral of [ FUNC(x) * (x**IPOW) * JC(x*r) * dx ]
C  
C  from zero to infinity, where JC is either J0 or J1.
C  
C   INPUT PARAMETERS:
C  ===================
C   NDEC specifies the number of samples per decade of the Hankel filters.
C   FTYPE (CHARACTER*2) determines the type of transform. FTYPE must be
C   given with capital letters.
C   FTYPE = 'J0' :  J0-transform with the calculated filter defined
C                 by NDEC
C   FTYPE = 'J1' :  J1-transform with the calculated filter defined
C                 by NDEC
C
C  The integral is calculated for logarithmically regularly spaced
C  values of R.
C
C  R = exp(N*DEL)  ,  min(N) = NGLO  ,   max(N) = NGHI
C
C  Calculation is performed for values of R between R1 and R2, where R1
C  is the second largest regular R-value smaller than RLO, and R2 is
C  the second smallest regular R-value greater than RHI.
C
C  If RHI is negative, computation is done for R = RLO only and RLO does
C  not have to be one of the regular R-values above, but can have any
C  value. This option makes it possible to avoid interpolation, if only
C  one R-value is needed, e.g. frequency soundings with only one
C  transmitter/receiver separation.
C
C  FUNC is the REAL kernel function.
C
C  TIME is an extra parameter included in the call to the kernel function FUNC.

C  IKEEP (INTEGER) determines if the kernel function values are to be kept from
C  the previous computation.
C  IKEEP = 1 :  keep the old values
C  IKEEP = 0 :  compute new values
C  If IKEEP = 1 then the kernel function FUNC must be the same as in the
C  previous call, but the input function FUNC is now multiplied with
C  a power of the integration variable.
C  
C  IPOW (INTEGER) the power of the integration variable with which the
C  input function is multiplied, when IKEEP = 1. Only the values
C  IPOW = -1, 1, 2 are allowed. If other powers or other functions are
C  needed, please use the convolution routine with the option of
C  including any user defined function.
C  This option has been implemented to save computation time for
C  related transforms often met in e.g. EM calculations.
C  
C  EPS is the desired relative accuracy of the calculation and
C  determines the truncation of the convolution.
C
C  If RLO, or RHI violate their restricted intervals determined by the
C  dimensioning of the subroutine, they are reset to their limiting
C  values.
C  
C   OUTPUT PARAMETERS:
C  ==================== 
C  NOUT is the number of radial values of ROUT, where the integral is calculated.
C  RCONV is a REAL array containing the values of the integral.
C  ROUT and RCONV must be dimensioned: REAL*8 ROUT(NG),RCONV(NG)
C  in the calling subroutine, where NG must be the same as in this
C  subroutine.
C  
C   DIMENSIONING OF THE ARRAYS:
C  =============================
C  The subroutine is dimensioned through the PARAMETER statements in
C  the first few lines of the code. The present dimensioning will accommodate
C  filters with up to 20 samples per decade (NDEC <= 20)
C
C  The filter coefficients used in the subroutine are from number NHLO 
C  to number NHHI.
C  
C  The maximum interval within which the output function can be
C  calculated is given by the numbers NGLO and NGHI (see above under
C  input parameters).
C  
C  NLIM is the number of terms in the discrete convolution that are
C  calculated initially without checking if the desired relative accuracy
C  has been reached. NLIM = 4*NDEC. This saves computation time.
C
C  J0 and J1 filters are calcualted the first time the
C  subroutine is called and if the parameters NDEC is changed.
C  This only takes 10-20ms.
C
C   18.12.2023:
C  =============
C  The subroutine is hardwired to produce filters with opening angle of
C  analyticity of 'pi' [ANGLE=1]. For a given sampling density NDEC, these will produce
C  the most accurate results. The parameter 'ANGLE' has been hardwired.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE FHTCONV(FTYPE,NDEC,RLO,RHI,FUNC,TIME,IKEEP,IPOW,EPS,
     #                   NOUT,ROUT,RCONV)
C ----------
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C ----------
      INCLUDE 'ARRAYSDIMBL.INC'
C ----------
      PARAMETER (NHLO = -600,NHHI = 300,NGLO = -60,NGHI = 150,
     #           NG = NGHI-NGLO+1)
      PARAMETER (NFLO = NGLO-NHHI,NFHI = NGHI-NHLO)
C ----------
      CHARACTER*2 FTYPE,FTYPEOLD
      INTEGER*4 NDEC,IKEEP,IPOW,NOUT,NDECOLD
C      REAL*8 RLO,RHI,EPS,ROUT(NG),RCONV(NG)
      REAL*8 RLO,RHI,EPS,ROUT(N_ARR),RCONV(N_ARR)
C ----------
      REAL*8 FILT(NHLO:NHHI),FC(NFLO:NFHI)
      COMMON /RCONVCOM/ FC,FILT,EPS1,KEEP,IPOT
C ----------
      REAL*8 TIME,FUNC
      EXTERNAL FUNC
C ----------
      REAL*8 J0FILT(NHLO:NHHI),J1FILT(NHLO:NHHI)
      REAL*8 HCFILT(NHLO:NHHI),HSFILT(NHLO:NHHI)
C ----------
      DATA ICALL /0/
      DATA SQPI2 /1.2533141373155D0/
C ----------
      REAL*8 EXPC
      EXTERNAL EXPC
C ----------
      SAVE ICALL,FTYPEOLD,NDECOLD,KMIN,KMAX
C ----------

      KEEP = IKEEP
      IPOT = IPOW
      NLIM = 4*NDEC
      NALO = NHLO
      NAHI = NHHI

C¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
C --- NDEC: THE NUMBER OF SAMPLES PER DECADE IS HARDWIRED HERE.
C --- ALL OTHER SETTINGS IN THE PROGRAM ARE HEREBY OVERRIDDEN.
C --- HIS MAKES IT EASY TO CHANGE THESAMPLING DENSITY.
C¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
      NDEC = 10
C¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
      
C----------------------------------------------------------------------
C --- CALCULATE MAX AND MIN RADIAL DISTANCE FROM INITIAL PARAMETERS
C----------------------------------------------------------------------
      IF (NDEC.NE.0) SC = NDEC/DLOG(100.D0)
      DEL = 0.5D0/SC
      RLOMIN = EXPC(NGLO*DEL)
      RHIMAX = EXPC(NGHI*DEL)

C----------------------------------------------------------------------
C --- CHECK FOR ERRONEUS PARAMETERS IN THE SUBROUTINE CALL.
C----------------------------------------------------------------------
      IF (RLO.LT.RLOMIN) THEN
      RLO = RLOMIN
      WRITE (*,1001) RLO,RLOMIN
 1001 FORMAT (1X,'RLO = ',1PE10.3,' IN "FHTCONVOL" OUT OF RANGE.',
     X/,' RLO HAS BEEN SET TO THE MINIMUM VALUE = ',1PE10.3/)
      END IF

      IF (RHI.GT.RHIMAX) THEN
      RHI = RHIMAX
      WRITE (*,1002) RHIMAX
 1002 FORMAT (1X,'RHI =  ',1PE10.3,' IN "FHTCONVOL" OUT OF RANGE.',
     #    /,' RHI HAS BEEN SET TO THE MAXIMUM VALUE = ',1PE10.3/)
      END IF

      IF ((IKEEP.LT.0).OR.(IKEEP.GT.1)) WRITE (*,1003)
 1003 FORMAT (1X,'PARAMETER IKEEP IN SUBROUTINE FHTCONV OUT OF RANGE')
 
      IF ((IKEEP.EQ.1).AND.(ICALL.EQ.0)) THEN
      WRITE (*,1004)
      STOP
      ENDIF
 1004 FORMAT (1X,'IKEEP = 1 AT FIRST CALL TO FHTCONV IS ILLEGAL',/,
     #          ' PROGRAM STOPPED!')

C-----------------------------------------------------------------------------
C --- CALCULATE FILTER COEFFICIENTS THE FIRST TIME THE SUBROUTINE IS CALLED
C-----------------------------------------------------------------------------
      IF (ICALL.EQ.0) THEN

      ANGLE = 1
      CALL HANKFILT(NDEC,ANGLE,  0.D0, 0.D0,NALO,NAHI,J0FILT)
      CALL HANKFILT(NDEC,ANGLE,  1.D0, 0.D0,NALO,NAHI,J1FILT)

      ICALL = 1
      FTYPEOLD = 'XX'
      NDECOLD = NDEC

      ENDIF  
C --- ENDIF: ENDIF this is the first time the routine is called
C--------------------------------------------------------------------

C------------------------------------------------------
C --- If NDEC has changed, calculate filters again
C------------------------------------------------------
      IF (NDEC.NE.NDECOLD) THEN
      CALL HANKFILT(NDEC,ANGLE,  0.D0, 0.D0,NALO,NAHI,J0FILT)
      CALL HANKFILT(NDEC,ANGLE,  1.D0, 0.D0,NALO,NAHI,J1FILT)
      NDECOLD = NDEC
      FTYPEOLD = 'XX'
      WRITE (21,'(A)') ' NEW FILTERS GENERATED FOR NDEC NEW'
      ENDIF

C----------------------------------------------------------------------
C --- THE CHOSEN FILTER COEFFICIENTS ARE PUT INTO THE ARRAY FILT
C----------------------------------------------------------------------
      IF (FTYPE.NE.FTYPEOLD) THEN

      IF (FTYPE.EQ.'J0') THEN
      DO I = NHLO,NHHI
      FILT(I) = J0FILT(I)
      ENDDO
      FTYPEOLD = 'J0'

      ELSEIF (FTYPE.EQ.'J1') THEN
      DO I = NHLO,NHHI
      FILT(I) = J1FILT(I)
      ENDDO
      FTYPEOLD = 'J1'

      ELSE
      WRITE (*,1005)
 1005 FORMAT
     #   (1X,' CHARACTER STRING "FTYPE" ILLEGAL IN CALL TO FHTCONV',/,
     #           ' PROGRAM STOPPED!')
      ENDIF

      ENDIF

C---------------------------------------------------------
C --- INITIALIZATIONS
C---------------------------------------------------------
C      WRITE (*,*) 'FHT: START INITIALISATIONS'

      EPS1 = EPS
      E    = EXPC(DEL)
      E1   = 1.D0/E
 
      IF (RHI.LT.0.D0) THEN

      R1   = RLO
      R2   = RLO
      NLO  = INT(LOG(RLO)/DEL+100)-100
      R10  = EXPC(NLO*DEL)
      XFAC = R10/R1
      NOUT = 1

C      WRITE (*,*) 'FHT: INIT IF RHI < 0'

      ELSE

      NLO = INT(LOG(RLO)/DEL+100)-101
      IF (NLO.LT.NGLO) NLO = NGLO
      NHI = INT(LOG(RHI)/DEL+100)-98
      IF (NHI.GT.NGHI) NHI = NGHI
      NOUT = NHI-NLO+1
      R1 = EXPC(NLO*DEL)
      R2 = EXPC(NHI*DEL)
      XFAC = 1.D0

      ENDIF

C-------------------------------------------------
C --- Define initial interval if FUNC is new
C-------------------------------------------------
      IF (IKEEP.EQ.0) THEN

      KMIN = NLO
      KMAX = NLO+NLIM
      
      X = E/R1
      DO K = KMIN,KMAX
      X = X*E1

C      WRITE (*,*) 'FHT: BEFORE CALLING FUNCTION',X,TIME,KMIN,KMAX,K

      FC(K) = FUNC(X,TIME)

      ENDDO

      ENDIF
      
C===========================================================
C --- CALCULATIONS BEGIN HERE
C===========================================================

C      WRITE (*,*) 'FHT: CALCULATIONS BEGIN HERE'
C---------------------------------------------------------
C --- Calculate for smallest R-value
C---------------------------------------------------------
      
      R = R1

      XFIRST = XFAC*E1**(KMIN-1)
      S = 0.D0
      CALL RCONVOF(NLO,KMIN,KMAX,1,S,XFIRST,E1)

      IF (S.EQ.0.D0) GOTO 9

      XFIRST = XFAC*E1**KMAX
      CALL RCONVON(NLO,KMAX+1,NLO-NHLO,1,S,XFIRST,E1,KLAST,FUNC,TIME)
      KMAX = KLAST

      XFIRST = XFAC*E1**KMIN
      CALL RCONVON(NLO,KMIN-1,NLO-NHHI,-1,S,XFIRST,E,KLAST,FUNC,TIME)
      KMIN = KLAST

    9 RCONV(1) = S/R
      ROUT(1) = R

      IF (NOUT.EQ.1) GOTO 8888

C---------------------------------------------
C --- Calculate for greatest R-value
C---------------------------------------------
      R = R2

      XFIRST = E1**(KMIN-1)
      S = 0.D0
      CALL RCONVOF(NHI,KMIN,KMAX,1,S,XFIRST,E1)
      IF (S.EQ.0.D0) GOTO 19

      XFIRST = E1**KMAX
      CALL RCONVON(NHI,KMAX+1,NHI-NHLO,1,S,XFIRST,E1,KLAST,FUNC,TIME)
      KMAX = KLAST

   19 RCONV(NOUT) = S/R
      ROUT(NOUT) = R
 
      IF (NOUT.EQ.2) GOTO 8888
       
C---------------------------------------------------------
C --- Calculate for all other R-values
C---------------------------------------------------------
      R = R1
      DO I = NLO+1,NHI-1
      S = 0.D0
      XFIRST = E1**(KMIN-1)
      CALL RCONVOF(I,KMIN,KMAX,1,S,XFIRST,E1)
      R = R*E
      RCONV(I-NLO+1) = S/R
      ROUT(I-NLO+1) = R
      ENDDO

 8888 CONTINUE

C --- END CALCULATIONS

      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                     C
C     S U B R O U T I N E   R C O N V O F                             C
C                                                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE RCONVOF(IR,K1,K2,KDEL,S,X,DX)
C ----------
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C ----------
      INTEGER*4 IR,K1,K2,KDEL
      REAL*8 S,X,DX
C ----------
      PARAMETER (NHLO = -600,NHHI = 300,NGLO = -60,NGHI = 150,
     #           NG = NGHI-NGLO+1)
      PARAMETER (NFLO = NGLO-NHHI,NFHI = NGHI-NHLO)
      REAL*8 FILT(NHLO:NHHI),FC(NFLO:NFHI)
      COMMON /RCONVCOM/ FC,FILT,EPS1,KEEP,IPOT
C ----------

      IF (IPOT.EQ.0) THEN
      DO K = K1,K2,KDEL
      S = S+FC(K)*FILT(IR-K)
      ENDDO
      ENDIF
      
      IF (IPOT.EQ.-1) THEN
      DO K = K1,K2,KDEL
      X = X*DX
      S = S+FC(K)*FILT(IR-K)/X
      ENDDO
      ENDIF

      IF (IPOT.EQ.1) THEN
      DO K = K1,K2,KDEL
      X = X*DX
      S = S+FC(K)*FILT(IR-K)*X
      ENDDO
      ENDIF

      IF (IPOT.EQ.2) THEN
      DO K = K1,K2,KDEL
      X = X*DX
      S = S+FC(K)*FILT(IR-K)*X*X
      ENDDO
      ENDIF

      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                     C
C     S U B R O U T I N E   R C O N V O N                             C
C                                                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE RCONVON(IR,K1,K2,KDEL,S,X,DX,KLAST,FUNC,TIME)
C ----------
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C ----------
      INTEGER*4 IR,K1,K2,KDEL,KLAST
      REAL*8 S,X,DX,FUNC
      EXTERNAL FUNC
C ----------
      PARAMETER (NHLO = -600,NHHI = 300,NGLO = -60,NGHI = 150,
     #           NG = NGHI-NGLO+1)
      PARAMETER (NFLO = NGLO-NHHI,NFHI = NGHI-NHLO)
      REAL*8 FILT(NHLO:NHHI),FC(NFLO:NFHI)
      COMMON /RCONVCOM/ FC,FILT,EPS1,KEEP,IPOT
C ----------

      IF (IPOT.EQ.0) THEN
      DO K = K1,K2,KDEL
      X = X*DX
      FC(K) = FUNC(X,TIME)
      SDEL = FC(K)*FILT(IR-K)
      S = S+SDEL
      IF (ABS(SDEL/S).LT.EPS1) GOTO 100
      ENDDO
      ENDIF

      IF (IPOT.EQ.-1) THEN
      DO K = K1,K2,KDEL
      X = X*DX
      FC(K) = FUNC(X,TIME)
      SDEL = FC(K)*FILT(IR-K)/X
      S = S+SDEL
      IF (ABS(SDEL/S).LT.EPS1) GOTO 100
      ENDDO
      ENDIF

      IF (IPOT.EQ.1) THEN
      DO K = K1,K2,KDEL
      X = X*DX
      FC(K) = FUNC(X,TIME)
      SDEL = FC(K)*FILT(IR-K)*X
      S = S+SDEL
      IF (ABS(SDEL/S).LT.EPS1) GOTO 100
      ENDDO
      ENDIF

      IF (IPOT.EQ.2) THEN
      DO K = K1,K2,KDEL
      X = X*DX
      FC(K) = FUNC(X,TIME)
      SDEL = FC(K)*FILT(IR-K)*X*X
      S = S+SDEL
      IF (ABS(SDEL/S).LT.EPS1) GOTO 100
      ENDDO
      ENDIF

  100 KLAST = K

      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                C
C     S U B R O U T I N E   H A N K F I L T                      C
C                                                                C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                C
C     Subroutine HANKFILT calculates filter coefficients         C
C     for Optimized Fast Hankel Transform filters for a given    C
C     set of parameters.                                         C
C                                                                C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE HANKFILT(NDEC,ANGLE,ANY,AMY,NLO,NHI,FILT)
C ----------
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C ----------
      PARAMETER (NHLO = -600,NHHI = 300,NGLO = -60,NGHI = 150,
     #           NG = NGHI-NGLO+1)
      REAL*8 FILT(NHLO:NHHI)
C ----------
      COMPLEX*16 J,GAMQLN,Z0
      EXTERNAL GAMQLN
C ----------
      REAL*8 RLO0(30),RMI0(800)
      COMPLEX*16 RLO1(30),RHI1(80),RMI1(800)
      COMMON/ARRAY/ RLO0,RMI0,RLO1,RHI1,RMI1
C ----------
      DATA J,PI /(0.D0,1.D0),3.141592653589793D0/
C ----------

      EPS = 1.D-16
      IOPT = 1
      OMEGA = ANGLE
      ISHIFT = 0
      DEL0 = 0.D0

      SC = FLOAT(NDEC)/LOG(100.D0)
      DEL = 0.5D0/SC
      A = 0.5D0/(SC*PI*OMEGA)

C      WRITE (21,*)
C     #        'ANY,AMY,NDEC,SC,A,EPS,IOPT,OMEGA,NLO,NHI,ISHIFT,DEL0'
C      WRITE (21,*) ANY,AMY,NDEC,SC,A,EPS,IOPT,OMEGA,NLO,NHI,ISHIFT,DEL0

      IF (ISHIFT.EQ.1) THEN
      Z0 = SC+J*A*SC
      DEL0 = -(LOG(2.D0)-DIMAG(GAMQLN(Z0,ANY,AMY))*DEL/PI)
      DEL0 = DEL0-DEL*INT(DEL0/DEL)
      IF (DEL0.LT.0.D0) DEL0 = DEL0+DEL
      END IF

C --- CALCULATION OF FILTER COEFFICIENTS

      CALL FILINIT(ANY,AMY,NDEC,SC,A,EPS,IOPT,OMEGA)

      DO I = NLO,NHI
      V = I*DEL-DEL0
      FILT(I) = FILCOA(V)
      ENDDO

C --- END CALCULATION OF FILTER COEFFICIENTS

C      ER = ERROR(A,SC,OMEGA,IOPT)

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                C
C     R E A L * 8   F U N C T I O N   E R R O R                  C
C                                                                C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                C
C     REAL*8 FUNCTION ERROR COMPUTES THE ERROR EXPRESSION IN     C
C     THE THEORY FOR FAST HANKEL TRANSFORMS.                     C
C                                                                C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL*8 FUNCTION ERROR(A,SC,OM,IOPT)
C ----------
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C ----------
      EXTERNAL SUM,SUM1
C ----------
      REAL*8 EXPC
      EXTERNAL EXPC
      DATA PI /3.141592653589793D0/
C ----------

      OMEGA = OM*PI

      IF (IOPT.EQ.1) THEN
        A = 1.D0/(2.D0*SC*OMEGA)
        PA = PI/A
        EPA = EXPC(-PA)
        ERROR = (SUM1(EPA)+EPA*(0.5D0*PA+0.25D0))/(PI*OMEGA)
      ELSE
        TA = 2.D0*A*SC*OMEGA
        TP = TA*PI/A
        F1 = PI*TA*EXPC(-TP)/SIN(PI*TA)
        F2 = 2.D0*TA*TA*SUM(TP,TA)
        ERROR = SC*(F1-F2)/TP
      ENDIF

      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                C
C     R E A L * 8   F U N C T I O N   S U M                      C
C                                                                C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL*8 FUNCTION SUM(TP,TA)
C ----------
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C ----------
      REAL*8 EXPC
      EXTERNAL EXPC
C ----------

      E = EXPC(-TP/TA)

      IF (E.LT.1.D-09) THEN
        SUM = E/(1.D0-TA*TA)
        RETURN
      ENDIF

      R = -1.D0
      SUM = 0.D0

      I = 0
   10 I = I+1
      R = -R*E
      SUMDEL = R/(I*I-TA*TA)
      SUM = SUM+SUMDEL
      C = ABS(SUMDEL/SUM)
      IF (C.GT.1.D-09) GOTO 10
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                C
C     R E A L *8   F U N C T I O N   S U M 1                     C
C                                                                C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL*8 FUNCTION SUM1(X)
C ----------
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C ----------
      R = X*X
      E = -X
      SUM1 = R/3.D0
      IF (X.LT.1.D-09) RETURN

      I = 2
   10 I = I+1
      R = R*E
      SUMDEL = R/(I*I-1)
      SUM1 = SUM1+SUMDEL
      C = ABS(SUMDEL/SUM1)
      IF (C.GT.1.D-09) GOTO 10
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                C
C     R E A L * 8   F U N C T I O N   F I L C O A                C
C                                                                C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                C
C     REAL*8 FUNCTION FILCOA COMPUTES FILTERCOEFFICIENTS         C
C     USED FOR THE FAST HANKEL TRANSFORM FORMULATED IN THE       C
C     FOLLOWING MANNER:                                          C
C                                                                C
C     INTEGRAL FROM ZERO TO INFINITY OF THE INTEGRAND            C
C                                                                C
C     F(X) * (X**AMY) * J ANY(XS) DX                             C
C                                                                C
C     THE FILTER COEFFICIENTS MAY BE COMPUTED TO A RELATIVE      C
C     ACCURACY OF APPROX. 1.E-12.                                C
C                                                                C
C     PARAMETERS:                                                C
C                                                                C
C     ANY1:   THE ORDER OF THE BESSEL FUNCTION.                  C
C     AMY1:   THE POWER TO WHICH X IS RAISED IN THE INTEGRAND.   C
C     NDEC1:  THE CUT-OFF FREQUENCY GIVEN AS SAMPLES PR DECADE.  C
C             IF NDEC1 = 0 SC1 IS USED AS CUT-OFF FREQUENCY.     C
C     SC1:    THE CUT-OFF FREQUENCY USED IF NDEC1 = 0.           C
C     A1:     THE SMOOTHNESS PARAMETER OF THE FILTER.            C
C     EPS1:   THE WANTED RELATIVE ERROR OF THE FILTER COEFF.     C
C     IOPT:   IF IOPT = 1 THEN A1 = 1./(2.*SC1*OMEGA*PI)         C
C             THAT IS, THE FUNCTION COMPUTES OPTIMIZED FILT.COEF.C
C     OMEGA:  THE OPENING ANGLE OF THE REGION OF ANALYTICITY.    C
C             OMEGA IS GIVEN IN FRACTION OF PI.                  C
C     V:      THE ARGUMENT OF THE FUNCTION.                      C
C                                                                C
C     RESTRICTIONS:                                              C
C                                                                C
C     ANY1+AMY1+1 > 1                                            C
C     AMY1 < 0.5                                                 C
C     5 < NDEC1 < 500                                            C
C     1.1 < SC < 108.5                                           C
C     1/(2.*SC*PI) < A < 8/(2.*SC*PI)                            C
C     EPS > 1.E-16                                               C
C     0.125 < OMEGA < 1.                                         C
C                                                                C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL*8 FUNCTION FILCOA(V)
C ----------
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C ----------
      COMPLEX*16 J
      COMMON/CONST/ ANY,AMY,SC,A,EPS,
     # PI,TPI,PIDA,TPIDA,TPISC,TPIASC,FPIASC,TASC,PIDASC,
     # DEL,DEL1,ANM1,ALN2,TIMY,TIMNY,NARR,VM,VP,J
C ----------

C----------------------------------------------------------------
C --- BRANCH INTO THREE DIFFERENT FUNCTIONS TO COMPUTE THE
C --- FILTER COEFFICIENTS ACCORDING TO V.
C----------------------------------------------------------------

      IF (V.LE.VM) THEN
        FILCOA = FILLOW(V)
      ELSEIF (V.GE.VP) THEN
        FILCOA = FILHIGH(V)
      ELSE
        FILCOA = FILMID(V)
      ENDIF

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                               C
C     S U B R O U T I N E   F I L I N I T                       C
C                                                               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                               C
C     SUBROUTINE FILINIT COMPUTES CONSTANTS AND PARAMETERS      C
C     USED IN THE FOLLOWING SUBROUTINES AND FUNCTIONS           C
C                                                               C
C     DEL:   THE SAMPLING DENSITY                               C
C     DEL1:  THE SAMPLING DENSITY IN THE FOURIER SUM IN THE     C
C            FUNCTION FILMID.                                   C
C     NARR:  THE NESCESSARY NUMBER OF TERMS IN THE FOURIER      C
C            SUM IN THE FUNCTION FILMID.                        C
C     VM,VP: IF V < VM THE FILT.COEF. IS COMPUTED IN THE LOWER  C
C            RANGE BY THE FUNCTION FILLOW.                      C
C            IF VM < V < VP THE FILT.COEF. IS COMPUTED IN THE   C
C            MIDDLE RANGE BY THE FUNCTION FILMID.               C
C            IF V > VP THE FILT.COEF. IS COMPUTED IN THE UPPER  C
C            RANGE BY THE FUNCTION FILHIGH                      C
C                                                               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE FILINIT(ANY1,AMY1,NDEC1,SC1,A1,EPS1,IOPT,OMEGA)
C ----------
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C ----------
      COMPLEX*16 J
      COMMON/CONST/ ANY,AMY,SC,A,EPS,
     # PI,TPI,PIDA,TPIDA,TPISC,TPIASC,FPIASC,TASC,PIDASC,
     # DEL,DEL1,ANM1,ALN2,TIMY,TIMNY,NARR,VM,VP,J
C ----------

      PI = 3.141592653589793D0
      J = (0.D0,1.D0)

C----------------------------------------------------------------
C --- THE PARAMETERS IN THE CALL OF REAL FUNCTION FILCOA
C --- ARE ENTERED IN THE COMMON BLOCK /CONST/
C----------------------------------------------------------------

      ANY = ANY1
      AMY = AMY1
      SC = SC1
      IF(NDEC1.NE.0) SC = NDEC1/LOG(100.D0)
      A = A1
      IF(IOPT.EQ.1) A = 1.D0/(2.D0*SC*OMEGA*PI)
      EPS = EPS1

C----------------------------------------------------------------
C --- CALCULATION OF CONSTANTS
C----------------------------------------------------------------

      TPI = 2.D0*PI
      PIDA = PI/A
      TPIDA = TPI/A
      TPISC = TPI*SC
      TPIASC = TPISC*A
      TASC = TPIASC/PI
      PIDASC = TPI/TASC
      FPIASC = 2.D0*TPIASC
      DEL = 0.5D0/SC
      ALN2 = LOG(2.D0)
      ANM1 = ANY+AMY+1.D0
      TIMY = EXPC(AMY*ALN2)
      TIMNY = EXPC(-ANY*ALN2)

      VM = 2.D0
      VP = 5.08D0*(SC**0.0625D0)

      M1 = INT(VM/DEL)
      M2 = INT(VP/DEL)+1
      NF3 = M2-M1-1+2
      DEL1 = 1.D0/(NF3*DEL)

C----------------------------------------------------------------
C --- 36.85 CORRESPONDS TO A RELATIVE ERROR ON THE FOURIER SUM OF 1E-16
C --- 46.05 CORRESPONDS TO A RELATIVE ERROR ON THE FOURIER SUM OF 1E-20
C --- 55.26 CORRESPONDS TO A RELATIVE ERROR ON THE FOURIER SUM OF 1E-24
C----------------------------------------------------------------
      NARR = INT((1.D0+50.D0*A/PI)*SC/DEL1)+1

C----------------------------------------------------------------
C --- COEFFICIENTS OF THE SUMMATIONS ARE CALCULATED
C----------------------------------------------------------------
      CALL ARRAYS

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                               C
C     S U B R O U T I N E   A R R A Y S                         C
C                                                               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                               C
C     SUBROUTINE ARRAYS COMPUTES ARRAYS USED IN THE SUMMATION   C
C     OF POLES IN THE LOWER, MIDDLE, AND UPPER RANGES           C
C                                                               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE ARRAYS
C ----------
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C ----------
      COMPLEX*16 J
      COMMON/CONST/ ANY,AMY,SC,A,EPS,
     # PI,TPI,PIDA,TPIDA,TPISC,TPIASC,FPIASC,TASC,PIDASC,
     # DEL,DEL1,ANM1,ALN2,TIMY,TIMNY,NARR,VM,VP,J
C ----------
      REAL*8 RLO0(30),RMI0(800)
      COMPLEX*16 RLO1(30),RHI1(80),RMI1(800)
      COMMON/ARRAY/ RLO0,RMI0,RLO1,RHI1,RMI1
C ----------
      COMPLEX*16 Z,GAMQLN
C ----------

C----------------------------------------------------------------
C --- THE ARRAY USED IN THE LOWER RANGE FOR SUMMATION OF
C --- POLES FROM THE GAMMA FUNCTION QUOTIENT
C----------------------------------------------------------------
      X = ANM1/TPI
      EX = 1.D0/PI
      R = 1.D0
      RLO0(1) = PIAX(X)

      DO I = 1,29
        X = X+EX
        R = -R/(I*(I+ANY))
        RLO0(I+1) = R*PIAX(X)
      ENDDO

C----------------------------------------------------------------
C --- THE ARRAY USED IN THE LOWER RANGE FOR SUMMATION OF
C --- POLES FROM THE FILTER.
C----------------------------------------------------------------
      DO I = 1,30
        Z = SC-J*TASC*(I-0.5D0)
        RLO1(I) = GAMQLN(Z,ANY,AMY)
      ENDDO

C----------------------------------------------------------------
C --- THE ARRAY USED IN THE UPPER RANGE FOR SUMMATION OF
C --- POLES FROM THE FILTER
C----------------------------------------------------------------
      DO I = 1,80
        Z = SC+J*TASC*(I-0.5D0)
        RHI1(I) = GAMQLN(Z,ANY,AMY)
      ENDDO

C----------------------------------------------------------------
C --- THE ARRAY USED IN THE MIDDLE RANGE IN THE FOURIER SUM
C----------------------------------------------------------------
      DO I = 1,NARR
        X = DEL1*(I-1)
        RMI0(I) = PRAX(X)
      ENDDO

C----------------------------------------------------------------
C --- THE ARRAY USED IN THE MIDDLE RANGE IN THE FOURIER SUM
C----------------------------------------------------------------
      DO I = 2,NARR
        Z = DEL1*(I-1)
        RMI1(I) = GAMQLN(Z,ANY,AMY)
      ENDDO

      IF((ANY-AMY+1.D0).LT.1.D-12) THEN
        RMI1(1) = -1.D35
      ELSE
        RMI1(1) = GAMQLN((0.D0,0.D0),ANY,AMY)
      ENDIF

C --------------------------------------------
C 1201 FORMAT (1X,1PE12.5)
C 1202 FORMAT (1X,1PE12.5,1X,1PE12.5)
C
C      WRITE (21,*) ' '
C      WRITE (21,*) ' RLO0(1:30)'
C      WRITE (21,1201) (RLO0(I),I = 1,30)
C
C      WRITE (21,*) ' '
C      WRITE (21,*) ' RMI0(1:800)'
C      WRITE (21,1201) (RMI0(I),I = 1,NARR+2)
C
C      WRITE (21,*) ' '
C      WRITE (21,*) ' RLO1(1:30)'
C      WRITE (21,1202) (RLO1(I),I = 1,30)
C
C      WRITE (21,*) ' '
C      WRITE (21,*) ' RHI1(1:80)'
C      WRITE (21,1202) (RHI1(I),I = 1,80)
C
C      WRITE (21,*) ' '
C      WRITE (21,*) ' RMI1(1:800)'
C      WRITE (21,1202) (RMI1(I),I = 1,NARR+2)
C --------------------------------------------

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                               C
C     R E A L * 8   F U N C T I O N   F I L L O W               C
C                                                               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                               C
C     REAL*8 FUNCTION FILLOW COMPUTES THE FILT.COEF. IN THE     C
C     LOWER RANGE; THAT IS FOR V < VM                           C
C                                                               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL*8 FUNCTION FILLOW(V)
C ----------
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C ----------
      COMPLEX*16 J
      COMMON/CONST/ ANY,AMY,SC,A,EPS,
     # PI,TPI,PIDA,TPIDA,TPISC,TPIASC,FPIASC,TASC,PIDASC,
     # DEL,DEL1,ANM1,ALN2,TIMY,TIMNY,NARR,VM,VP,J
C ----------
      REAL*8 RLO0(30),RMI0(800)
      COMPLEX*16 RLO1(30),RHI1(80),RMI1(800)
      COMMON/ARRAY/ RLO0,RMI0,RLO1,RHI1,RMI1
C ----------
      COMPLEX*16 SA2,SA2DEL
      COMPLEX*16 CEXPC
      EXTERNAL EXPC,CEXPC
C ----------

      VLN = V-ALN2

C----------------------------------------------------------------
C --- SUMMATION OF POLES FROM THE GAMMA FUNCTION QUOTIENT
C----------------------------------------------------------------
      R = 1.D0
      E = EXPC(2.D0*VLN)
      SA1 = RLO0(1)

      I = 1
   10 I = I+1
      R = R*E
      SA1DEL = R*RLO0(I)
      SA1 = SA1+SA1DEL
      C = ABS(SA1DEL/SA1)
      IF (C.GT.EPS) GO TO 10

      S1 = SA1*EXPC(V*ANM1)*TIMNY/RGAMMA(ANY+1.D0)
      N1 = I

C----------------------------------------------------------------
C --- SUMMATION OF POLES FROM THE FILTER
C----------------------------------------------------------------
      EA2 = FPIASC*VLN
      RA2 = 0.D0
      SA2 = CEXPC(RLO1(1))

      I = 1
   20 I = I+1
      RA2 = RA2+EA2
      SA2DEL = CEXPC(RA2+RLO1(I))
      SA2 = SA2+SA2DEL
      C = ABS(SA2DEL/SA2)
      IF(C.GT.EPS) GO TO 20
      SA2 = CEXPC(J*TPISC*VLN)*SA2
      S2 = -A*TIMY*EXPC(TPIASC*VLN)*2.D0*DIMAG(SA2)
      N2 = I

      FILLOW = S1+S2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                               C
C     R E A L * 8   F U N C T I O N   F I L M I D               C
C                                                               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                               C
C     REAL*8 FUNCTION FILMID COMPUTES THE FILT.COEF. IN THE     C
C     MIDDLE RANGE; THAT IS FOR VM < V < VP.                    C
C                                                               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL*8 FUNCTION FILMID(V)
C ----------
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C ----------
      COMPLEX*16 J
      COMMON/CONST/ ANY,AMY,SC,A,EPS,
     # PI,TPI,PIDA,TPIDA,TPISC,TPIASC,FPIASC,TASC,PIDASC,
     # DEL,DEL1,ANM1,ALN2,TIMY,TIMNY,NARR,VM,VP,J
C ----------
      REAL*8 RLO0(30),RMI0(800)
      COMPLEX*16 RLO1(30),RHI1(80),RMI1(800)
      COMMON/ARRAY/ RLO0,RMI0,RLO1,RHI1,RMI1
C ----------
      COMPLEX*16 E,R,SB1,SB1DEL,CA2,RA2,RA3,SA2,SA2DEL
      COMPLEX*16 CB2,RB2,RB3,SB2,SB2DEL
      COMPLEX*16 CEXPC
      EXTERNAL EXPC,CEXPC
C ----------

      I0 = INT(SC/DEL1)+2
      VLN = V-ALN2

C----------------------------------------
C --- THE FOURIER SUM
C----------------------------------------

      E = J*TPI*VLN*DEL1
      R = (0.D0,0.D0)
      SB1 = (0.D0,0.D0)

      I = 1
   10 I = I+1
      R = R+E
      SB1DEL = RMI0(I)*CEXPC(RMI1(I)+R)
      SB1 = SB1+SB1DEL
      IF(I.LE.I0) GO TO 10
      C = ABS(SB1DEL/SB1)
      IF(C.GT.EPS) GO TO 10
      S1 = 2.D0*DREAL(SB1)
      N0 = I

      S1 = TIMY*DEL1*(S1+RMI0(1)*DREAL(CEXPC(RMI1(1))))

C----------------------------------------------------------------
C --- THE SUM OF FILT.COEF. IN THE LOWER RANGE IS PERFORMED
C --- AS THE USUAL CALCULATION IN THE LOWER RANGE.
C----------------------------------------------------------------

C----------------------------------------------------------------
C --- THE SUM OF POLES FROM THE GAMMA FUNCTION QUOTIENT
C----------------------------------------------------------------
      R1 = 1.D0
      E1 = EXPC(2.D0*(VLN-1.D0/DEL1))
      CA1 = EXPC(-ANM1/DEL1)
      R2 = CA1
      E2 = EXPC(-2.D0/DEL1)
      R3 = 1.D0/(1.D0-R2)
      SA1 = RLO0(1)*R3

      I = 1
   20 I = I+1
      R1 = R1*E1
      R2 = R2*E2
      R3 = 1.D0/(1.D0-R2)
      SA1DEL = RLO0(I)*R1*R3
      SA1 = SA1+SA1DEL
      C = ABS(SA1DEL/SA1)
      IF(C.GT.EPS) GO TO 20
      N1 = I

      S2 = CA1*TIMNY*EXPC(V*ANM1)*SA1/RGAMMA(ANY+1.D0)

C----------------------------------------------------------------
C --- THE SUM OF POLES FROM THE FILTER.
C----------------------------------------------------------------

      R1 = 0.D0
      E1 = FPIASC*(VLN-1.D0/DEL1)
      CA2 = CEXPC(-(TPIASC+J*TPISC)/DEL1)
      RA2 = CA2
      EA2 = EXPC(-FPIASC/DEL1)
      RA3 = 1.D0/(1.D0-RA2)
      SA2 = CEXPC(RLO1(1))*RA3

      I = 1
   30 I = I+1
      R1 = R1+E1
      RA2 = RA2*EA2
      RA3 = 1.D0/(1.D0-RA2)
      SA2DEL = CEXPC(RLO1(I)+R1)*RA3
      SA2 = SA2+SA2DEL
      C = ABS(SA2DEL/SA2)
      IF(C.GT.EPS) GO TO 30
      SA2 = CA2*CEXPC(J*TPISC*VLN)*SA2
      N2 = I

      S3 = -A*TIMY*EXPC(TPIASC*VLN)*2.D0*DIMAG(SA2)

C----------------------------------------------------------------
C --- THE SUM OF FILT.COEF. IN THE UPPER RANGE IS PERFORMED
C --- AS THE USUSAL CALCULATION IN THE UPPER RANGE
C----------------------------------------------------------------

C----------------------------------------------------------------
C --- THE SUM OF POLES FROM THE FILTER.
C----------------------------------------------------------------
      R1 = 0.D0
      E1 = -FPIASC*(VLN+1.D0/DEL1)
      CB2 = CEXPC(-(TPIASC-J*TPISC)/DEL1)
      RB2 = CB2
      EB2 = EXPC(-FPIASC/DEL1)
      RB3 = 1.D0/(1.D0-RB2)
      SB2 = CEXPC(RHI1(1))*RB3

      I = 1
   40 I = I+1
      R1 = R1+E1
      RB2 = RB2*EB2
      RB3 = 1.D0/(1.D0-RB2)
      SB2DEL = CEXPC(RHI1(I)+R1)*RB3
      SB2 = SB2+SB2DEL
      C = ABS(SB2DEL/SB2)
      IF(C.GT.EPS) GO TO 40
      SB2 = CB2*CEXPC(J*TPISC*VLN)*SB2
      N3 = I

      S4 = A*TIMY*EXPC(-TPIASC*VLN)*2.D0*DIMAG(SB2)

      FILMID = S1-S2-S3-S4

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                               C
C     R E A L *8   F U N C T I O N   F I L H I G H              C
C                                                               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                               C
C     REAL*8 FUNTION FILHIGH COMPUTES THE FILT.COEF. IN THE     C
C     UPPER RANGE; THAT IS FOR V > VP.                          C
C                                                               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL*8 FUNCTION FILHIGH(V)
C ----------
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C ----------
      COMPLEX*16 J
      COMMON/CONST/ ANY,AMY,SC,A,EPS,
     # PI,TPI,PIDA,TPIDA,TPISC,TPIASC,FPIASC,TASC,PIDASC,
     # DEL,DEL1,ANM1,ALN2,TIMY,TIMNY,NARR,VM,VP,J
C ----------
      REAL*8 RLO0(30),RMI0(800)
      COMPLEX*16 RLO1(30),RHI1(80),RMI1(800)
      COMMON/ARRAY/ RLO0,RMI0,RLO1,RHI1,RMI1
C ----------
      COMPLEX*16 SB2,SB2DEL
      COMPLEX*16 CEXPC
      EXTERNAL EXPC,CEXPC
C ----------

      VLN = V-ALN2

C------------------------------------------------
C --- SUMMATION OF POLES FROM THE FILTER.
C------------------------------------------------
      EB2 = -FPIASC*VLN
      RB2 = 0.D0
      SB2 = CEXPC(RHI1(1))

      I = 1
   10 I = I+1
      RB2 = RB2+EB2
      SB2DEL = CEXPC(RB2+RHI1(I))
      SB2 = SB2+SB2DEL
      C = ABS(SB2DEL/SB2)
      IF(C.GT.EPS) GO TO 10
      SB2 = CEXPC(J*TPISC*VLN)*SB2
      N2 = I

      FILHIGH = A*TIMY*EXPC(-TPIASC*VLN)*2.D0*DIMAG(SB2)

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                               C
C     R E A L * 8   F U N C T I O N   P I A X                   C
C                                                               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                               C
C     REAL*8 FUNCTION PIAX COMPUTES VALUES OF THE               C
C     FILTER FUNCTION ON THE NEGATIVE IMAGINARY AXIS.           C
C                                                               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL*8 FUNCTION PIAX(Y)
C ----------
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C ----------
      COMPLEX*16 J
      COMMON/CONST/ ANY,AMY,SC,A,EPS,
     # PI,TPI,PIDA,TPIDA,TPISC,TPIASC,FPIASC,TASC,PIDASC,
     # DEL,DEL1,ANM1,ALN2,TIMY,TIMNY,NARR,VM,VP,J
C ----------
      EXTERNAL EXPC
C ----------

      E = EXPC(-PIDA)
      E2 = E*E
      AY = PIDASC*Y

      PIAX = DEL*(1.D0-E2)/(1.D0+E2+2.D0*E*COS(AY))
 
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                               C
C     R E A L * 8   F U N C T I O N   P R A X                   C
C                                                               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                               C
C     REAL*8 FUNCTION PRAX COMPUTES VALUES OF THE FILTER        C
C     FUNCTION ON THE REAL AXIS.                                C
C                                                               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL*8 FUNCTION PRAX(X)
C ----------
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C ----------
      COMPLEX*16 J
      COMMON/CONST/ ANY,AMY,SC,A,EPS,
     # PI,TPI,PIDA,TPIDA,TPISC,TPIASC,FPIASC,TASC,PIDASC,
     # DEL,DEL1,ANM1,ALN2,TIMY,TIMNY,NARR,VM,VP,J
C ----------

      A1 = PIDASC*(X-SC)
      A2 = -PIDASC*(X+SC)
      P0 = 1.D0 - EXPC(-TPIDA)
      P1 = 1.D0 + EXPC(A1)
      P2 = 1.D0 + EXPC(A2)

      PRAX = DEL*P0/(P1*P2)

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   C
C     C O M P L E X * 1 6   F U N C T I O N   G A M Q L N           C
C                                                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   C
C     COMPLEX*16 FUNCTION GAMQLN COMPUTES THE LOGARITHM OF THE      C
C     FOLLOWING FUNCTION:                                           C
C                                                                   C
C     GAMQ = G((ANY+AMY+1)/2-J*PI*Z) / G((ANY-AMY+1)/2+J*PI*Z)      C
C                                                                   C
C     WHERE G IS THE GAMMA FUNCTION                                 C
C           Z IS A COMPLEX NUMBER                                   C
C           ANY IS REAL                                             C
C           AMY IS REAL < 0.5                                       C
C           J = SQRT(-1)                                            C
C           PI = ARCCOS(-1)                                         C
C                                                                   C
C     RESTRICTIONS:  Z # -J*(N+(ANY+AMY+1)/2)/PI                    C
C                    Z # +J*(N+(ANY-AMY+1)/2)/PI                    C
C                    REAL PART OF (ANY+AMY+1) > 0                   C
C                                                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      COMPLEX*16 FUNCTION GAMQLN(Z,ANY,AMY)
C ----------
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C ----------
      COMPLEX*16 Z,J,Z1,Z2
      COMPLEX*16 CGAMMA,PRIVAL
      EXTERNAL CGAMMA,PRIVAL
C ----------
      DATA J,PI/(0.D0,1.D0),3.141592653589793D0/
C ----------

      ANY1 = 0.5D0*(ANY+1.D0)
      AMY1 = AMY*0.5D0
      Z1 = ANY1-J*PI*Z+AMY1
      Z2 = ANY1+J*PI*Z-AMY1

      GAMQLN = CGAMMA(Z1)-CGAMMA(Z2)
      GAMQLN = PRIVAL(GAMQLN)

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                     C
C     C O M P L E X * 1 6   F U N C T I O N   C G A M M A             C
C                                                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                     C
C     COMPLEX*16 function CGAMMA computes the logarithm               C
C     of the complex gamma function.                                  C
C                                                                     C
C     Restrictions :  Z not equal to -N, N = 0,1,2,3, ...             C
C                                                                     C
C    The value of CGAMMA does not nessessarily belong to the prin-    C
C    cipal branch of the logarithm.                                   C
C                                                                     C
C    COMPLEX*16 function CGAMMA calls the complex function CGAMASY,   C
C    which computes the logarithm of the complex gamma function       C
C    using the asymptotic expansion from "HANdbook of Mathematical    C
C    Functions" by Abramowitz and Stegun, formula 6.1.42.             C
C                                                                     C
C    The absolute error on CGAMMA may be set to different levels:     C
C    EPS< 1.e-09 :  set A = 4.3 and NCOF = 5                          C
C    EPS< 1.e-12 :  set A = 6.6 and NCOF = 6                          C
C    EPS< 1.e-16 :  set A = 7.9 and NCOF = 9                          C
C                                                                     C
C    04.05.94  NBC                                                    C
C                                                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      COMPLEX*16 FUNCTION CGAMMA(Z)
C ----------
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C ----------
      COMPLEX*16 Z,RECUR
C ----------
      COMPLEX*16 CGAMASY,PRIVAL
      EXTERNAL CGAMASY,PRIVAL
C ----------
      DATA A /7.9D0/
C ----------

      X = DREAL(Z)
      Y = DIMAG(Z)
      C = ABS(Z)

      IF (ABS(Y).GE.A) THEN
      N = 0
      IF (ABS(X).LT.1.D0) N = 2
      ELSEIF (X.GE.A) THEN
      N = 0
      ELSE
      N = INT(A-X)+1
      ENDIF

      RECUR = (0.D0,0.D0)
      IF (N.GT.0) THEN
        DO I = 1,N
          RECUR = RECUR+LOG(Z+I-1)
        ENDDO
      ENDIF

      CGAMMA = CGAMASY(Z+N)-RECUR

      CGAMMA = PRIVAL(CGAMMA)

      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                     C
C     C O M P L E X * 1 6   F U N C T I O N    C G A M A S Y          C
C                                                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                     C
C     04.05.94  NBC                                                   C
C                                                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      COMPLEX*16 FUNCTION CGAMASY(Z)
C ----------
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C ----------
      COMPLEX*16 Z,Z2,R
      REAL*8 C(11)
C ----------
      DATA NCOF /9/
      DATA C / 8.333333333333333D-02,-2.777777777777778D-03,
     #         7.936507936507937D-04,-5.952380952380952D-04,
     #         8.417508417508418D-04,-1.917526917526918D-03,
     #         6.410256410256410D-03,-2.955065359477124D-02,
     #         1.796443723688306D-01,-1.392432216905901D+00,
     #         1.340286404416839D+01/
C ----------

      CGAMASY = (Z-0.5D0)*LOG(Z)-Z+0.9189385332046727D0

      Z2 = Z*Z
      R = Z

      DO I = 1,NCOF
        R = R/Z2
        CGAMASY = CGAMASY+C(I)*R
      ENDDO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                               C
C     C O M P L E X * 1 6   F U N C T I O N   P R I V A L       C
C                                                               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                               C
C     COMPLEX*16 FUNCTION PRIVAL(W) RETURNS THE PRINCIPAL VALUE C
C     OF THE ARGUMENT W, E.I. THE REAL PART OF W IS UNCHANGED   C
C     WHILE THE IMAGINARY PART IS GREATER THAN -PI AND SMALLER  C
C     THAN OR EQUAL TO +PI.                                     C
C                                                               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      COMPLEX*16 FUNCTION PRIVAL(W)
C ----------
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C ----------
      COMPLEX*16 W,Z,J
C ----------
      DATA PI,TPI,J /3.14159265358979D0,6.28318530717959D0,(0.D0,1.D0)/
C ----------

      Z = W
      Z = Z-J*TPI*INT(DIMAG(Z)/TPI)

   10 Y = DIMAG(Z)
      IF (Y.LE.PI) GOTO 20
      Z = Z-J*TPI
      GOTO 10

   20 Y = DIMAG(Z)
      IF(Y.GT.-PI) GOTO 99
      Z = Z+J*TPI
      GOTO 20

   99 PRIVAL = Z

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                     C
C     R E A L * 8   F U N C T I O N   R G A M M A                     C
C                                                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                     C
C     Real function RGAMMA computes the real gamma function           C
C                                                                     C
C     Restrictions :  X not equal to -N, N = 0,1,2,3, ...             C
C                                                                     C
C    Real function CGAMMA calls the real function RGAMASY, which      C
C    computes the logarithm of the real gamma function using the      C
C    asymptotic expansion from "Handbook of Mathematical Functions    C
C    by Abramowitz and Stegun, formula 6.1.42.                        C
C                                                                     C
C    The relative error on RGAMMA may be set to different levels:     C
C                                                                     C
C    EPS < 1.e-09 :  set A = 3.7 and NCOF = 5                         C
C    EPS < 1.e-12 :  set A = 5.0 and NCOF = 7                         C
C    EPS < 1.e-16 :  set A = 7.1 and NCOF = 9                         C
C                                                                     C
C    21.05.94  NBC                                                    C
C                                                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL*8 FUNCTION RGAMMA(X)
C ----------
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C ----------
      EXTERNAL RGAMASY,EXPC
C ----------
      DATA A,PI /8.D0,3.141592653589793D0/
C ----------

      IF (X.LT.0.5D0) THEN
        X1 = 1.D0-X
        IREF = 1
      ELSE
        X1 = X
        IREF = 0
      ENDIF

      IF (X1.GE.A) THEN
        RGAMMA = EXPC(RGAMASY(X1))
      ELSE
        N = INT(A-X1)+1
        RECUR = 1.D0

        DO I = 1,N
          RECUR = RECUR*(X1+I-1)
        ENDDO

      RGAMMA = EXPC(RGAMASY(X1+N))/RECUR

      ENDIF

      IF (IREF.EQ.1) RGAMMA = PI/(RGAMMA*SIN(PI*X))

      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                     C
C     R E A L * 8   F U N C T I O N    R G A M A S Y                  C
C                                                                     C
C     21.05.94  NBC                                                   C
C                                                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL*8 FUNCTION RGAMASY(X)
C ----------
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C ----------
      REAL*8 C(11)
C ----------
      DATA NCOF /8/
      DATA C / 8.333333333333333D-02,-2.777777777777778D-03,
     #         7.936507936507937D-04,-5.952380952380952D-04,
     #         8.417508417508418D-04,-1.917526917526918D-03,
     #         6.410256410256410D-03,-2.955065359477124D-02,
     #         1.796443723688306D-01,-1.392432216905901D+00,
     #         1.340286404416839D+01/
C ----------

      RGAMASY = (X-0.5D0)*LOG(X)-X+0.9189385332046727D0

      X2 = X*X
      R = X

      DO I = 1,NCOF
        R = R/X2
        RGAMASY = RGAMASY+C(I)*R
      ENDDO

      RETURN
      END

