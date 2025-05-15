CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   R E A L * 8   S U B R O U T I N E   T E M R E S P I P
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   SUBROUTINE TEMRESPIP calculates the TEM response of a given model
C   including the IP effect, optionally including the system response
C   and the derivatives.
C
C   OUPUT
C  ======
C   NTOUT:     The number of samples of the output.
C   TIMESOUT:  The delay times of the samples.
C   RESPOUT:   The response samples of the samples.
C   DRESPOUT:  The derivatives.
C   
C   February 2024 / NBC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TEMRESPIP(NTOUT,TIMESOUT,RESPOUT,DRESPOUT)
C ----------
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C ----------
      INCLUDE 'ARRAYSDIMBL.INC'
      INCLUDE 'MODELBL.INC'
      INCLUDE 'INSTRUMENTBL.INC'
      INCLUDE 'POLYGONBL.INC'
      INCLUDE 'RESPBL.INC'
      INCLUDE 'IPBL.INC'
      INCLUDE 'WAVEBL.INC'
C ----------
      INTEGER*4 IPARM
      COMMON /DFTEMBL/ IPARM
C ----------
      CHARACTER FTYPE*2
C ----------
      REAL*8 TN(N_ARR),TNL(N_ARR),RESPSTP(N_ARR),RESPIMP(N_ARR)
      REAL*8 DSTPDH(N_ARR),DSTPDP(N_ARR,N_ARR)
      REAL*8 TTREP(N_ARR),RESPREP(N_ARR),HH(N_ARR)
      REAL*8 TTS(N_ARR),RESPCONV(N_ARR)
      REAL*8 ROUT(N_ARR),RCONV(N_ARR)
C ----------
      REAL*8 TIMESOUT(N_ARR),RESPOUT(N_ARR),DRESPOUT(N_ARR,N_ARR)
C ----------
      REAL*8 FTEMIP,DFTEMIP
      EXTERNAL FTEMIP,DFTEMIP
C ----------

C............................................................
C --- Option: write info to output.
C............................................................
      IWRI = 1

C............................................................
      IF (IWRI.EQ.1) THEN
        WRITE (21,'(A)') 'TEMRESPIP CALLED'
      ENDIF
C............................................................

C=============================================================================
C --- FIND BASIC TIME ARRAY AS A FUNCTION OF THE REPETITION FREQUENCY
C=============================================================================

C------------------------------------------------------------------
C --- Number of parameter derivatives according to IMLM
C------------------------------------------------------------------
C --- If IMLM = 1, then there are only derivatives wrt. layer conductivities,
C --- plus the derivative wrt. HTX.
C --- If IMLM = 0, then there are derivatives wrt. layer conductivities
C --- and layer thicknesses, plus the derivative wrt. HTX.
C------------------------------------------------------------------
      IF (IDERIV.EQ.0) THEN
        NPARM = 0
      ELSE
        NPARM = NLAY+1
        IF (IMLM.EQ.0) THEN
          NPARM = NLAY + (NLAY-1) + 1
        ENDIF
      ENDIF

C-----------------------------------------------------------------------------------------
C --- Trap detecting inconsistent input.
C-----------------------------------------------------------------------------------------
      HPER = 0.5D0/REPFREQ
      IF (ABS(TWAVE(1)).GT.HPER) THEN
        WRITE ( *,'(A)') 'ERROR: Input inconsistency: Program stopped'
        WRITE ( *,'(A)') 'ERROR: Waveform longer than halfperiod.'
        WRITE (21,'(A)') 'ERROR: Input inconsistency: Program stopped'
        WRITE (21,'(A)') 'ERROR: Waveform longer than halfperiod.'
        STOP
      ENDIF

C-----------------------------------------------------------------------------------------
C --- If no repetition, max time is two halfperiod.
C --- If repetition is modelled, max time is five halfperiods.
C-----------------------------------------------------------------------------------------
      IF (IREP.EQ.0) THEN
        TMAX = 2.D0*HPER
      ELSE
        TMAX = 5.D0*HPER
      ENDIF

C-----------------------------------------------------------------------------------------
C --- The earliest delay time in the primary array is always 10 ns.
C-----------------------------------------------------------------------------------------
      TMIN = 1.D-8
      TLO = TMIN
      THI = TMAX

C----------------------------------------------------------------------
C --- Set up basic time array for caluclations. Density: 10/decade.
C----------------------------------------------------------------------
      NUMDECADES = INT(LOG10(THI/TLO)) + 1
      FAC = 10.D0**0.1D0
      TN(1) = TLO*1.D0
      DO I = 2,10*NUMDECADES
        TN(I) = TN(I-1)*FAC
        TNL(I) = LOG(TN(I))
      ENDDO
      NTT = 10*NUMDECADES

C*******************************************************************************
C=============================================================================
C --- CALCULATION OF FORWARD RESPONSES AND DERIVATIVES.
C=============================================================================
C*******************************************************************************

C----------------------------------------------
C --- Choice: Central loop or offset loop.
C----------------------------------------------
      IF (ICEN.EQ.1) THEN
        FTYPE = 'J1'
        RLO   = TXRAD
      ELSE
        FTYPE = 'J0'
        RLO   = RTXRX
      ENDIF

      NDEC  = 15
      RHI   = -1.D0
      EPS = 0.0001D0

C--------------------------------------------------------------------------
C --- FOR ALL RESPONSE TYPES, CALCULATE STEP RESPONSE WITH OPTIONAL FILTERS
C--------------------------------------------------------------------------
C      IF (IRESPTYPE.EQ.0 .OR. IRESPTYPE.EQ.2) THEN

      DO IT = 1,NTT
        TIME = TN(IT)

      IKEEP = 0
      IPOW  = 0
        CALL FHTCONV (FTYPE,NDEC,RLO,RHI,FTEMIP,TIME,IKEEP,IPOW,EPS,
     #                NOUT,ROUT,RCONV)
        RESPSTP(IT) = RCONV(1)

C--------------------------------------------------------------------------
C --- CALCULATE DERIVATIVE OF STEP RESPONSE WRT. TX HEIGHT
C--------------------------------------------------------------------------
      IF (IDERIV.GT.0) THEN

      IKEEP = 1
      IPOW  = 1
        CALL FHTCONV (FTYPE,NDEC,RLO,RHI,FTEMIP,TIME,IKEEP,IPOW,EPS,
     #                NOUT,ROUT,RCONV)
        DSTPDH(IT) = -2.D0*RCONV(1)

      DSTPDP(IT,NPARM) = DSTPDH(IT)

      ENDIF

      ENDDO

C..........................................................................
      IF (IWRI.EQ.1) THEN
      WRITE (21,'(A)') 'FROM TEMRESP: STEP RESP W/FILTERS CALCULATED'
      WRITE (21,'(A)') 'FROM TEMRESP: TN, RESPSTP, DB/DH'
      WRITE (21,1012) (TN(I),RESPSTP(I),DSTPDH(I),I=1,NTT)
 1012 FORMAT (3(2X,1PE11.4))      
      ENDIF
C..........................................................................

C      ENDIF
C --- ENDIF: Step responses are calculated

C======================================================================
C --- ACCURATE DERIVATIVES.
C======================================================================
      IF (IDERIV.GT.0 .AND. IRESPTYPE.EQ.2) THEN

      NDEC  = 15
      RHI   = -1.D0
      IKEEP = 0
      IPOW  = 0
      EPS = 0.0001D0

      DO J = 1,NPARM-1
      IPARM = J

        DO I = 1,NTT
          TIME = TN(I)
          CALL FHTCONV (FTYPE,NDEC,RLO,RHI,DFTEMIP,TIME,IKEEP,IPOW,EPS,
     #                  NOUT,ROUT,RCONV)
          DSTPDP(I,J) = RCONV(1)
        ENDDO

      ENDDO

      ENDIF
C --- ENDIF: If accurate derivatives are calculated.

Cいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいい
C
C --- THIS CONCLUDES THE CALCULATIONS OF THE BASIC RESPONSE AND DERIVATIVE ARRAYS.
C --- THE FOLLOWING CODE CONVOLVES THE BASIC RESULTS WITH THE SYSTEM RESPONSE.
C
Cいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいい

C===============================================================
C --- PRINCIPAL STEP AND IMPULSE RESPONSES.
C --- VALUES ARE JUST TRANSFERRED TO OUTPUT ARRAYS.
C===============================================================
      IF (IRESPTYPE.EQ.0) THEN
        NTOUT = NTT
        DO I = 1,NTT
          TIMESOUT(I) = TN(I)
          RESPOUT(I) = RESPSTP(I)
        ENDDO
      RETURN
      ENDIF

C---------------------------------------------------------------
C --- IMPULSE RESPONSE IS NOW CALUCLATED
C --- AS THE TIME DERIVATIVE OF THE STEP RESPONSE.
C---------------------------------------------------------------
      IF (IRESPTYPE.EQ.1) THEN
      CALL DERIV1ST (NTT,TN,RESPSTP,RESPIMP)

        NTOUT = NTT
        DO I = 1,NTT
          TIMESOUT(I) = TN(I)
          RESPOUT(I)  = -RESPIMP(I)
        ENDDO
      RETURN
      ENDIF

C==========================================================================
C --- IN THE CASE OF IRESPTYPE = 2, TRANSFER RESPONSES TO OUTPUT ARRAYS
C==========================================================================
      IF (IRESPTYPE.EQ.2) THEN
        NTOUT = NTT
        DO I = 1,NTT
          TIMESOUT(I) = TN(I)
          RESPOUT(I) = RESPSTP(I)
        ENDDO
      ENDIF

C===============================================================
C --- MODELLING SYSTEM RESPONSE EFFECTS ON THE RESPONSE
C===============================================================

C---------------------------------------------------------------------
C --- TIME DOMAIN MODELING OF REPETITION.
C---------------------------------------------------------------------
      IF (IRESPTYPE.EQ.2 .AND. IREP.EQ.1) THEN

        RFREQ = REPFREQ
        CALL REPMOD (NTOUT,TIMESOUT,RESPOUT,RFREQ,NTREP,TTREP,RESPREP)

C------------------------------------------------------------
C --- TRANSFER TO OUTPUT ARRAYS
C------------------------------------------------------------
      NTOUT = NTREP
      DO I = 1,NTREP
        TIMESOUT(I) = TTREP(I)
        RESPOUT(I)  = RESPREP(I)
      ENDDO

      ENDIF
C --- ENDIF: END IF REPETITION IS MODELLED

C............................................................................
      IF (IWRI.EQ.1) THEN
      WRITE (21,'(A)')
     #      'FROM TEMRESP, RESPONSE: AFTER MODELLING OF REPETITION: '
      WRITE (21, 2001)
     #      (TIMESOUT(J),RESPOUT(J) , J=1,NTOUT)
      ENDIF
 2001 FORMAT (2(2X,1PE11.4))
C............................................................................

C---------------------------------------------------------------------
C --- CONVOLUTION WITH WAVEFORM
C---------------------------------------------------------------------
      IF (IRESPTYPE.EQ.2 .AND. IWCONV.EQ.1) THEN

      CALL WAVECONV(NTOUT,TIMESOUT,RESPOUT,NTTS,TTS,RESPCONV)

C------------------------------------------------------------
C --- TRANSFER TO OUTPUT ARRAYS
C------------------------------------------------------------
      NTOUT = NTTS
      DO I = 1,NTTS
        TIMESOUT(I) = TTS(I)
        RESPOUT(I) = RESPCONV(I)
      ENDDO

      ENDIF
C --- ENDIF: END IF WAVEFORM CONVOLUTION IS MODELLED.

C............................................................
      IF (IWRI.EQ.1) THEN
        WRITE (21,'(A)') 'FROM TEMRESP, RESPONSE: AFTER WAVE CONV'
        WRITE (21, 2002)
     #      (TIMESOUT(J),RESPOUT(J), J=1,NTOUT)
      ENDIF
 2002 FORMAT (2(2X,1PE11.4))
C............................................................

C***************************************************************************
C==========================================================================
C --- THIS CONCLUDES MODELLING SYSTEM RESPONSE EFFECTS ON THE RESPONSE
C==========================================================================
C***************************************************************************


C***********************************************************************
C=====================================================================
C --- SYSTEM RESPONSE MODELLING OF THE DERIVATIVES
C=====================================================================
C***********************************************************************

C-----------------------------------------------------------------------
C --- TRANSFER DERIVATIVES TO OUTPUT ARRAYS.
C --- THE DELAY TIMES AND THEIR NUMBER MIGHT HAVE BEEN MODIFIED
C --- IN THE SYSTEM RESPONSE MODELLING OF THE RESPONSE
C-----------------------------------------------------------------------
      IF (IRESPTYPE.EQ.2) THEN
        DO J = 1,NPARM
          DO I = 1,NTT
            DRESPOUT(I,J) = DSTPDP(I,J)
          ENDDO
        ENDDO
      ENDIF
C-----------------------------------------------------------------------

      IF (IDERIV.GT.0) THEN

C-------------------------------------------------------------------
C --- MODELLING OF REPETITION FOR ALL DERIVATIVES.
c --- THIS IS DONE FOR BOTH APPROXIMATE AND ACCURATE DERIVATIVES.
C-------------------------------------------------------------------
      IF (IRESPTYPE.EQ.2 .AND. IREP.EQ.1) THEN

      DO J = 1,NPARM

        DO I = 1,NTT
        HH(I) = DRESPOUT(I,J)
        ENDDO

      RFREQ = REPFREQ

      CALL REPMOD(NTT,TIMESOUT,HH,RFREQ,NTREP,TTREP,RESPREP)

C------------------------------------------------------------
C --- STORE RESULT
C------------------------------------------------------------
        DO I = 1,NTREP
          DRESPOUT(I,J) = RESPREP(I)
        ENDDO

      ENDDO

C------------------------------------------------------------
C --- TRANSFER TO OUTPUT ARRAYS
C------------------------------------------------------------
      NTOUT = NTREP
      DO I = 1,NTREP
        TIMESOUT(I) = TTREP(I)
      ENDDO

      ENDIF

C............................................................................
      IF (IWRI.EQ.1) THEN
      WRITE (21,'(A)')
     #      'FROM TEMRESP, DERIV: AFTER REPETITION: TIME,RESP,DERESP '
      DO I = 1,NTOUT
      WRITE (21, 2009)
     #      TIMESOUT(I),RESPOUT(I),(DRESPOUT(I,J) , J=1,NPARM)
      ENDDO
      ENDIF
 2009 FORMAT (20(2X,1PE11.4))
C............................................................................

C-----------------------------------------------------------------------
C --- CONVOLUTION WITH WAVEFORM FOR ALL DERIVATIVES.
c --- THIS IS DONE FOR BOTH APPROXIMATE AND ACCURATE DERIVATIVES.
C-----------------------------------------------------------------------
      IF (IRESPTYPE.EQ.2 .AND. IWCONV.EQ.1) THEN

      DO J = 1,NPARM

        DO I = 1,NTOUT
          HH(I) = DRESPOUT(I,J)
        ENDDO

      CALL WAVECONV(NTOUT,TIMESOUT,HH, NTTS,TTS,RESPCONV)

C------------------------------------------------------------
C --- STORE RESULT
C------------------------------------------------------------
      DO I = 1,NTTS
        DRESPOUT(I,J) = RESPCONV(I)
      ENDDO

      ENDDO
C --- ENDDO: LOOP OVER THE DERIVATIVES

C------------------------------------------------------------
C --- TRANSFER TO OUTPUT ARRAYS
C------------------------------------------------------------
      NTOUT = NTTS
      DO I = 1,NTTS
        TIMESOUT(I) = TTS(I)
      ENDDO

      ENDIF
C --- ENDIF: IF CONVOLUTION WITH WAVEFORM

C............................................................................
      IF (IWRI.EQ.1) THEN
      WRITE (21,'(A)')
     #  'FROM TEMRESP: DERIVATIVES: AFTER WAVE CONV: T,RESP,DERIV'
      DO I = 1,NTOUT
      WRITE (21, 2017) TIMESOUT(I),RESPOUT(I),(DRESPOUT(I,K), K=1,NPARM)
      ENDDO
      ENDIF
 2017 FORMAT(30(2X,1PE11.4))
C............................................................................

C=======================================================================
C --- THIS CONCLUDES MODELLING OF SYSTEM RESPONSE ON THE DERIVATIVES
C=======================================================================

      ENDIF
C --- ENDIF: IF DERIVATIVES ARE INCLUDED IN THE COMPUTATIONS

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     R E A L * 8   F U N C T I O N   F T E M I P
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  REAL*8 FUNCTION FTEMIP calculates the inverse Laplace transform of the
C  kernel function of the TEM response including the IP parameters.
C
C  FTEMIP accommodates both step and impulse responses.
C  It delivers the kernel function as a function of wavenumber, 'X' (scalar).
C  The input parameters 'TIME' is the delay time (scalar).
C
C  February 2024 / NBC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL*8 FUNCTION FTEMIP(X,TIME)
C----------
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C ----------
      INCLUDE 'ARRAYSDIMBL.INC'
      INCLUDE 'MODELBL.INC'
      INCLUDE 'IPBL.INC'
      INCLUDE 'INSTRUMENTBL.INC'
      INCLUDE 'RESPBL.INC'
      INCLUDE 'POLYGONBL.INC'
C ----------
      REAL*8 FTM(N_ARR),SIGIP(N_ARR),W(N_ARR)
      REAL*8 EXPC, J1
      EXTERNAL EXPC, J1
C ----------
      REAL*8 PI,MU0
      DATA PI,MU0 /3.14159265358979D+00,1.256637061435917D-06/
C ----------

C-------------------------------------------------------------
C --- The Gaver-Stehfest weight factors.
C-------------------------------------------------------------
      REAL*8 GSCOEF(16)
      DATA GSCOEF / 
     # -3.968253968253968D-04, 2.133730158730159D+00,
     # -5.510166666666667D+02, 3.350016111111111D+04,
     # -8.126651111111111D+05, 1.007618376666667D+07,
     # -7.324138297777778D+07, 3.390596320730159D+08,
     # -1.052539536278571D+09, 2.259013328583333D+09,
     # -3.399701984433333D+09, 3.582450461700000D+09,
     # -2.591494081366667D+09, 1.227049828766667D+09,
     # -3.427345554285714D+08, 4.284181942857143D+07/
C-------------------------------------------------------------

      Q = LOG(2.D0)/TIME
      X2 = X*X

C-------------------------------------------------
C --- Outer loop over the laplace variable s.
C-------------------------------------------------
      DO IS = 1,16
      S = Q * IS

C---------------------------------------------------------------
C --- Compute Cole-Cole conductivities to model IP effect.
C---------------------------------------------------------------
C --- The Cole-Cole parameters are:
C --- CHAIP, TAUIP, POWIP corresponding to "m" (chargeability),
C --- "tau" (time constant), and "c" (the power).
C---------------------------------------------------------------
      DO I = 1,NLAY
      HH = 1.D0/(1.D0+(S*TAUIP(I))**POWIP(I))
      RHOIP = RHON(I)*( 1.D0 - CHAIP(I)*(1.D0-HH) )
      SIGIP(I) = 1.D0/RHOIP
      ENDDO

C---------------------------------------------------
C --- Compute wavenumbers.
C---------------------------------------------------
      DO I = 1,NLAY
        W(I) = SQRT(X2 + MU0*S*SIGIP(I))
      ENDDO

C-------------------------------------------------
C --- RECURSION OF THE KERNEL FUNCTION
C-------------------------------------------------
      G = 0.D0

C---------------------------------------------------------
C --- Recursion through the layers if more than 1 layer.
C---------------------------------------------------------
      IF (NLAY.GT.1) THEN

        DO I = NLAY-1,1,-1
          H = W(I) + W(I+1)
          R = MU0*S*(SIGIP(I)-SIGIP(I+1))/(H*H)
          E = EXP(-2.D0*W(I)*THKN(I))
          G = E*(R+G)/(1.D0+R*G)
        ENDDO

      ENDIF
C --- END: IF THERE IS MORE THAN 1 LAYER

C---------------------------------------------------------
C --- Recursion step from 1st layer to the air halfspace.
C --- This is entry point for a halfspace model.
C---------------------------------------------------------
      H  = X + W(1)
      R  = -MU0*S*SIGIP(1)/(H*H)
      FT = (R+G)/(1.D0+R*G)

C-------------------------------------------------------------------------------------
C --- Multiply with filter functions for the recording system
C-------------------------------------------------------------------------------------
      IF (NFILT.GT.0) THEN
      DO I = 1,NFILT
        FT = FT / (1+S/(2*PI*FILTFREQ(I)))
      ENDDO
      ENDIF

C----------------------------------------------------------------
C --- The Step response is the basis for all response types.
C --- Therefore always divide by Laplace variable
C----------------------------------------------------------------
      FTM(IS) = -FT/S

      ENDDO
C --- ENDDO: LOOP OVER THE 16 S-VALUES.
C-------------------------------------------------------------------------------------
      
C-------------------------------------------------------------------------------
C --- MULTIPLY WITH THE GAVER-STEHFEST WEIGHT FACTORS.
C-------------------------------------------------------------------------------
      RESP = 0.D0
      DO IS = 1,16
      RESP = RESP + GSCOEF(IS)*FTM(IS)
      ENDDO

      RESP = Q*RESP

C----------------------------------------------------------
C --- Multiply with height terms
C----------------------------------------------------------
      FAC = 0.D0
      FAC = FAC + FLOAT(ISHTX1*ISHRX1) * EXPC(-X*(HTX1+HRX1))
      FAC = FAC + FLOAT(ISHTX1*ISHRX2) * EXPC(-X*(HTX1+HRX2))
      FAC = FAC + FLOAT(ISHTX2*ISHRX1) * EXPC(-X*(HTX2+HRX1))
      FAC = FAC + FLOAT(ISHTX2*ISHRX2) * EXPC(-X*(HTX2+HRX2))
      
      RESP = RESP * FAC

C----------------------------------------------------------------------------------
C --- MULTIPLY WITH WAVENUMBER TERMS AND NORMALISING FACTOR 1/(4*PI)
C----------------------------------------------------------------------------------
C --- For offset loop, the factor is:
C --- LAMBDA^2 * J1(LAMBDA*A) / (0.5*LAMBDA*A) = (2*LAMBDA/A)*J1(LAMBDA*A)
C --- where "a" is the equivalent radius of the Tx loop.
C --- The Bessel function involved in the FHT convolution is J0.
C----------------------------------------------------------------------------------
C --- For central loop, the factor is the same, except that the "J1(lambda*a)"
C --- term is the bessel function involved in the FHT convolution, meaning that
C --- the J1 term in the above expression is not part of the kernel function.
C --- Also notice that the factor MU0*0.25D0/PI is applied in this function.
C----------------------------------------------------------------------------------
C --- For a polygobal Tx loop, ICEN is set to ICEN=1. Then the kernel will be correct.
C --- The integration over the sides will be called from the TEMRESPIP function.
C --- For a polygobal loop, all that remains is the ds*Y/R factor in the integration.
C----------------------------------------------------------------------------------
      IF (NPOLY.GT.0) THEN
        RESP = RESP*X
      ENDIF

      IF (NPOLY.EQ.0) THEN
        RESP = RESP * (2.D0*X/TXRAD)
        IF (ICEN.NE.1) THEN
          RESP = RESP * J1(X*TXRAD)
        ENDIF
      ENDIF

      RESP = MU0*0.25D0*RESP/PI

      FTEMIP = RESP

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     R E A L * 8   F U N C T I O N   D F T E M I P
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   REAL*8 FUNCTION DFTEMIP calculates the derivative of the kernel
C   function wrt. the model parameters for a step response.
C   It delivers the derivatives of the kernel function as a function of
C   wavenumber, 'X', and the delay time'TIME', both scalars.
C   The derivative is identified by it number 'IPARM' which is
C   transferred through the COMMON block DFTEMBL.
C
C     February 2024 / NBC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL*8 FUNCTION DFTEMIP(X,TIME)
C ----------
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C ----------
      INCLUDE 'ARRAYSDIMBL.INC'
      INCLUDE 'MODELBL.INC'
      INCLUDE 'INSTRUMENTBL.INC'
      INCLUDE 'RESPBL.INC'
      INCLUDE 'IPBL.INC'
      INCLUDE 'POLYGONBL.INC'
C ----------
      INTEGER*4 IPARM
      COMMON /DFTEMBL/ IPARM
C ----------
      REAL*8 TIME,SIGS(0:N_ARR),THKS(0:N_ARR)
      REAL*8 W(0:N_ARR),G(0:N_ARR),E(0:N_ARR),R(0:N_ARR)
      REAL*8 CHAIPS(0:N_ARR),TAUIPS(0:N_ARR),POWIPS(0:N_ARR)
      REAL*8 FTM(N_ARR),SIGIPFAC(0:N_ARR),SIGIP(0:N_ARR)
C ----------
      REAL*8 EXPC,J1
      EXTERNAL EXPC,J1
C ----------
      REAL*8 PI,MU0
      DATA PI,MU0 /3.14159265358979D0,1.256637061435917D-06/
C-------------------------------------------------------------
C --- The Gaver-Stehfest weight factors.
C-------------------------------------------------------------
      REAL*8 GSCOEF(16)
      DATA GSCOEF / 
     # -3.968253968253968D-04, 2.133730158730159D+00,
     # -5.510166666666667D+02, 3.350016111111111D+04,
     # -8.126651111111111D+05, 1.007618376666667D+07,
     # -7.324138297777778D+07, 3.390596320730159D+08,
     # -1.052539536278571D+09, 2.259013328583333D+09,
     # -3.399701984433333D+09, 3.582450461700000D+09,
     # -2.591494081366667D+09, 1.227049828766667D+09,
     # -3.427345554285714D+08, 4.284181942857143D+07/
C-------------------------------------------------------------

C------------------------------------------------------
C --- Add the conductivity of the air halfspace
C------------------------------------------------------
      SIGS(0) = 0.D0
      THKS(0) = 0.D0
      CHAIPS(0) = 0.D0
      TAUIPS(0) = TAUIP(1)
      POWIPS(0) = POWIP(1)
      DO I = 1,NLAY-1
      SIGS(I) = SIGN(I)
      THKS(I) = THKN(I)
      CHAIPS(I) = CHAIP(I)
      TAUIPS(I) = TAUIP(I)
      POWIPS(I) = POWIP(I)
      ENDDO
      SIGS(NLAY) = SIGN(NLAY)
      CHAIPS(NLAY) = CHAIP(NLAY)
      TAUIPS(NLAY) = TAUIP(NLAY)
      POWIPS(NLAY) = POWIP(NLAY)

      Q = LOG(2.D0)/TIME
      X2 = X*X

C-------------------------------------------------
C --- Outer loop over the laplace variable s.
C-------------------------------------------------
      DO IS = 1,16
      S = Q * IS

C---------------------------------------------------
C --- Compute wavenumbers
C---------------------------------------------------
      DO I = 0,NLAY
      HH = 1.D0/(1.D0+(S*TAUIPS(I))**POWIPS(I))
      RHOFAC = ( 1.D0 - CHAIPS(I)*(1.D0-HH) )
      SIGIPFAC(I) = 1.D0/RHOFAC
      SIGIP(I) = SIGS(I) * SIGIPFAC(I)
      W(I) = SQRT(X2+MU0*S*SIGIP(I))      
      ENDDO

C---------------------------------------------------
C --- Computation of GAMMA(I)
C---------------------------------------------------
      G(NLAY) = 0.D0
      DO I = NLAY-1,0,-1
      H = W(I)+W(I+1)
      R(I+1) = MU0*S*(SIGIP(I)-SIGIP(I+1))/(H*H)
      E(I) = EXPC(-2.D0*W(I)*THKS(I))
      G(I) = E(I)*(R(I+1)+G(I+1))/(1.D0+R(I+1)*G(I+1))
      ENDDO

C---------------------------------------------------
C --- Derivatives wrt resistivities
C---------------------------------------------------
      IF (IPARM.LE.NLAY) THEN

      IPP = IPARM

C---------------------------------------------------
C --- Derivatives wrt lower halfspace resistivity
C---------------------------------------------------
      IF (IPP.EQ.NLAY) THEN
      H = W(NLAY) + W(NLAY-1)
      DT = -2.D0*W(NLAY-1)/(H*H)
      DG = E(NLAY-1)*DT
      IF (NLAY.EQ.1) GOTO 1000
      LM2 = NLAY-2
        DO I = LM2,0,-1
          H = 1.D0 + G(I+1)*R(I+1)
          DG = DG*E(I) * (1.D0-R(I+1)*R(I+1))/(H*H)
        ENDDO
      GOTO 1000
      ENDIF

C---------------------------------------------------
C --- Other resistivities
C---------------------------------------------------
      H = W(IPP)+W(IPP+1)
      DT = 2.0*W(IPP+1)/(H*H)
      H = 1.D0+G(IPP+1)*R(IPP+1)
      DG = -2.D0*THKS(IPP)*G(IPP)+DT*E(IPP)*(1.D0-G(IPP+1)
     #     *G(IPP+1))/(H*H)

      IPP1 = IPP-1
      H = W(IPP1)+W(IPP1+1)
      DT = -2.0*W(IPP1)/(H*H)
      H = 1.D0+G(IPP1+1)*R(IPP1+1)
      DG = (DG*(1.D0-R(IPP1+1)*R(IPP1+1))+DT*(1.D0-G(IPP1+1)
     #     *G(IPP1+1)))*E(IPP1)/(H*H)
      IF (IPP1.EQ.0) GOTO 1000

      IPP2 = IPP-2
      DO I = IPP2,0,-1
      H = 1.D0+G(I+1)*R(I+1)
      DG = DG*E(I)*(1.D0-R(I+1)*R(I+1)) / (H*H)
      ENDDO

      GOTO 1000

      ELSE

C---------------------------------------------------
C --- Derivatives wrt thicknesses
C---------------------------------------------------
      IPP = IPARM-NLAY
      DG = -2.D0*W(IPP)*G(IPP)
      IPP1 = IPP-1
      DO I = IPP1,0,-1
      H = 1.D0 + G(I+1)*R(I+1)
      DG = DG*E(I) * (1.D0-R(I+1)*R(I+1)) / (H*H)
      ENDDO
      GOTO 2000

      ENDIF

C--------------------------------------------------------------
C --- CHANGE FROM DERIVATIVES WITH RESPECT TO WAVENUMBER TO
C --- DERIVATIVES WITH RESPECT TO CONDUCTIVITIES.
C--------------------------------------------------------------
 1000 DG =  DG*0.5D0*MU0*S * SIGIPFAC(IPP)/W(IPP)

 2000 CONTINUE

C-------------------------------------------------------------------------------------
C --- Multiply with filter functions for the recording system
C-------------------------------------------------------------------------------------
      IF (NFILT.GT.0) THEN
      DO I = 1,NFILT
        DG = DG / (1+S/(2*PI*FILTFREQ(I)))
      ENDDO
      ENDIF

      FTM(IS) = -DG/S

      ENDDO
C --- ENDDO: LOOP OVER LAPLACE VARIABLE

C-------------------------------------------------------------------------------
C --- MULTIPLY WITH THE GAVER-STEHFEST WEIGHT FACTORS.
C-------------------------------------------------------------------------------
      RESP = 0.D0
      DO IS = 1,16
      RESP = RESP + GSCOEF(IS)*FTM(IS)
      ENDDO

      RESP = Q*RESP

C----------------------------------------------------------
C --- Multiply with height terms
C----------------------------------------------------------
      FAC = 0.D0
      FAC = FAC + FLOAT(ISHTX1*ISHRX1) * EXPC(-X*(HTX1+HRX1))
      FAC = FAC + FLOAT(ISHTX1*ISHRX2) * EXPC(-X*(HTX1+HRX2))
      FAC = FAC + FLOAT(ISHTX2*ISHRX1) * EXPC(-X*(HTX2+HRX1))
      FAC = FAC + FLOAT(ISHTX2*ISHRX2) * EXPC(-X*(HTX2+HRX2))
      
      RESP = RESP * FAC

C----------------------------------------------------------------------------------
C --- MULTIPLY WITH WAVENUMBER TERMS AND NORMALISING FACTOR 1/(4*PI)
C----------------------------------------------------------------------------------
C --- For offset loop, the factor is:
C --- LAMBDA^2 * J1(LAMBDA*A) / (0.5*LAMBDA*A) = (2*LAMBDA/A)*J1(LAMBDA*A)
C --- where "a" is the equivalent radius of the Tx loop.
C --- The Bessel function involved in the FHT convolution is J0.
C----------------------------------------------------------------------------------
C --- For central loop, the factor is the same, except that the "J1(lambda*a)"
C --- term is the bessel function involved in the FHT convolution, meaning that
C --- the J1 term in the above expression is not part of the kernel function.
C --- Also notice that the factor MU0*0.25D0/PI is applied in this function.
C----------------------------------------------------------------------------------
C --- For a polygobal Tx loop, ICEN is set to ICEN=1. Then the kernel will be correct.
C --- The integration over the sides will be called from the TEMRESPIP function.
C --- For a polygobal loop, all that remains is the ds*Y/R factor in the integration.
C----------------------------------------------------------------------------------
      IF (NPOLY.GT.0) THEN
        RESP = RESP*X
      ENDIF

      IF (NPOLY.EQ.0) THEN
        RESP = RESP * (2.D0*X/TXRAD)
        IF (ICEN.NE.1) THEN
          RESP = RESP * J1(X*TXRAD)
        ENDIF
      ENDIF

      RESP = MU0*0.25D0*RESP/PI

      DFTEMIP = RESP

      RETURN
      END
