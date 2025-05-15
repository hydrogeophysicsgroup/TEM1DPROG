CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   R E A L * 8   S U B R O U T I N E   T E M  R E S P P O L Y
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   SUBROUTINE TEMRESPPOLY calculates the TEM response for a polygonal
C   Tx loop of a given model, optionally including the effect of the system response.
C   Also derivatives are, optionally, delivered on output.
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
      SUBROUTINE TEMRESPPOLY(NTOUT,TIMESOUT,RESPOUT,DRESPOUT)
C ----------
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C ----------
      INCLUDE 'ARRAYSDIMBL.INC'
      INCLUDE 'MODELBL.INC'
      INCLUDE 'INSTRUMENTBL.INC'
      INCLUDE 'POLYGONBL.INC'
      INCLUDE 'RESPBL.INC'
      INCLUDE 'WAVEBL.INC'
C ----------
      INTEGER*4 IPARM
      COMMON /DFTEMBL/ IPARM
C ----------
      CHARACTER FTYPE*2
C ----------
      REAL*8 TN(N_ARR),TNL(N_ARR),RESPSTP(N_ARR),RESPIMP(N_ARR)
      REAL*8 DSTPDP(N_ARR,N_ARR)
      REAL*8 TTREP(N_ARR),RESPREP(N_ARR),HH(N_ARR)
      REAL*8 TTS(N_ARR),RESPCONV(N_ARR)
      REAL*8 ROUT(N_ARR),RCONV(N_ARR)
C ----------
      REAL*8 TIMESOUT(N_ARR),RESPOUT(N_ARR),DRESPOUT(N_ARR,N_ARR)
C ----------
      REAL*8 FTEM,DFTEM
      EXTERNAL FTEM,DFTEM
C ----------

C............................................................
C --- Option: write info to output.
C............................................................
      IWRI = 0

C............................................................
      IF (IWRI.EQ.1) THEN
        WRITE (21,'(A)') 'TEMRESPPOLY CALLED'
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

      FTYPE = 'J1'
      RLO   = RSAMPMIN
      RHI   = RSAMPMAX
      NDEC  = 15
      EPS = 0.0001D0

C--------------------------------------------------------------------------
C --- FOR ALL RESPONSE TYPES, CALCULATE STEP RESPONSE WITH OPTIONAL FILTERS
C--------------------------------------------------------------------------
C      IF (IRESPTYPE.EQ.0 .OR. IRESPTYPE.EQ.2) THEN

      DO IT = 1,NTT
        TIME = TN(IT)

      IKEEP = 0
      IPOW  = 0
        CALL FHTCONV (FTYPE,NDEC,RLO,RHI,FTEM,TIME,IKEEP,IPOW,EPS,
     #                NOUT,ROUT,RCONV)

C----------------------------------
C --- CALL INTEGRATION ROUTINE
C----------------------------------
      CALL POLYINT(NOUT,ROUT,RCONV,TXINT)
      RESPSTP(IT) = TXINT

C--------------------------------------------------------------------------
C --- CALCULATE DERIVATIVE OF STEP RESPONSE WRT. TX HEIGHT
C--------------------------------------------------------------------------
      IF (IDERIV.GT.0) THEN

      IKEEP = 1
      IPOW  = 1
        CALL FHTCONV (FTYPE,NDEC,RLO,RHI,FTEM,TIME,IKEEP,IPOW,EPS,
     #                NOUT,ROUT,RCONV)

C----------------------------------
C --- CALL INTEGRATION ROUTINE
C----------------------------------
      CALL POLYINT(NOUT,ROUT,RCONV,TXINT)
      DSTPDP(IT,NPARM) = -2*TXINT

      ENDIF

      ENDDO

C..........................................................................
      IF (IWRI.EQ.1) THEN
      WRITE (21,'(A)') 'FROM TEMRESPPOLY: STEP RESP W/FILTERS'
      WRITE (21,'(A)') 'FROM TEMRESPPOLY: TN, RESPSTP, DB/DH'
      WRITE (21,1012) (TN(I),RESPSTP(I),DSTPDP(I,NPARM),I=1,NTT)
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
      RLO   = RSAMPMIN
      RHI   = RSAMPMAX
      IKEEP = 0
      IPOW  = 0
      EPS = 0.0001D0

      DO J = 1,NPARM-1
      IPARM = J

        DO I = 1,NTT
          TIME = TN(I)
          CALL FHTCONV (FTYPE,NDEC,RLO,RHI,DFTEM,TIME,IKEEP,IPOW,EPS,
     #                  NOUT,ROUT,RCONV)

C----------------------------------
C --- CALL INTEGRATION ROUTINE
C----------------------------------
      CALL POLYINT(NOUT,ROUT,RCONV,TXINT)
      DSTPDP(I,J) = TXINT

      ENDDO

      ENDDO

      ENDIF
C --- ENDIF: If accurate derivatives are calculated.

Cいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいい
C==========================================================================
C --- THIS CONCLUDES THE CALCULATIONS OF THE BASIC RESPONSE AND DERIVATIVE ARRAYS.
C --- THE FOLLOWING CODE CONVOLVES THE BASIC RESULTS WITH THE SYSTEM RESPONSE.
C==========================================================================
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
     #     'FROM TEMRESPPOLY, RESPONSE: AFTER MODELLING OF REPETITION:'
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
        WRITE (21,'(A)') 'FROM TEMRESPPOLY, RESPONSE: AFTER WAVE CONV'
        WRITE (21, 2002)
     #      (TIMESOUT(J),RESPOUT(J), J=1,NTOUT)
      ENDIF
 2002 FORMAT (2(2X,1PE11.4))
C............................................................

C==========================================================================
C --- THIS CONCLUDES MODELLING SYSTEM RESPONSE EFFECTS ON THE RESPONSE
C==========================================================================


C=====================================================================
C --- SYSTEM RESPONSE MODELLING OF THE DERIVATIVES
C=====================================================================

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
C-------------------------------------------------------------------
      IF (IRESPTYPE.EQ.2 .AND. IREP.EQ.1) THEN

      DO J = 1,NPARM

C        DO I = 1,NTOUT
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
     #   'FROM TEMRESPPOLY, DERIV: AFTER REPETITION: TIME,RESP,DERESP'
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
     #  'FROM TEMRESPPOLY: DERIVATIVES: AFTER WAVE CONV: T,RESP,DERIV'
      WRITE (*,*) NTOUT
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

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   S U B R O U T I N E   P O L Y I N T
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  Subroutine POLYINT performs the integration along the loop perimeter
C  over a picewise linear Tx loop of the secondary magnetic field
C  response from an electric dipole
C
C  INPUT
C =======
C  NA:           Number of sides of the apices of the Tx loop.
C  XA(1:NSIDES): X-coordinates of the apices of the Tx loop.
C  YA(1:NSIDES): Y-coordinates of the apices of the Tx loop.
C
C  OUTPUT
C ========
C  TXINT:        The output integral over the Tx loop segments.
C
C  February 2024 / NBC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE POLYINT(NA,XA,YA,TXINT)
C ----------
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C ----------
      INCLUDE 'ARRAYSDIMBL.INC'
      INCLUDE 'POLYGONBL.INC'
C ----------
      REAL*8 XA(NA),YA(NA),XAL(N_ARR),TXINT
      REAL*8 VAL(N_ARR),PRESP(N_ARR),XI(N_ARR),YB(N_ARR)
C ----------
      REAL*8 EXPC
      EXTERNAL EXPC
C ----------

      IWRI = 0

C........................................................................
      IF (IWRI.EQ.1) THEN
      WRITE (*,*) 'FROM POLYINT: TRANSFERRED CONVOLUTION RESULTS'
      WRITE (*,*) NA
      DO J = 1,NA
      WRITE (*,*) XA(J),YA(J)
      ENDDO
      ENDIF
C........................................................................

      RINT = 0.D0

C------------------------------------------------------
C --- FIND LOGARITHM OF FHT CONVOLUTION ARRAYS
C------------------------------------------------------
      DO I = 1,NA
      XAL(I) = LOG(XA(I))
      ENDDO

C-----------------------------------------------------------------------
C --- INTERPOLATION IN LOG(X) - LIN(Y)
C-----------------------------------------------------------------------
      IOP = 0
      EXTVAL = 0.D0
      CALL INTERP1D (NA,XAL,YA,NRSAMP,RSAMPL,YB,IOP,EXTVAL)

C........................................................................
      IF (IWRI.EQ.1) THEN
      WRITE (*,*) 'FROM POLYINT: INTERPOLATED RESULTS'
      WRITE (*,*) NRSAMP
      DO J = 1,NRSAMP
      WRITE (*,*) YSAMP,RSAMP(J)
      ENDDO
      ENDIF
C........................................................................

C----------------------------------------------------------------
C --- MULTIPLY WITH Y/R
C----------------------------------------------------------------
      DO J = 1,NRSAMP
      PRESP(J) =  YB(J)*YSAMP(J)/RSAMP(J)
      ENDDO

C---------------------------------------------------------
C --- INTEGRATE OVER ALL SIDES, ONE BY ONE
C---------------------------------------------------------
      RINT = 0.D0
      KSUM = 0

      DO K = 1,NPOLY

C---------------------------------------------------------
C --- SELECT ARRAY ELEMENTS BEFORE 10-POINT INTEGRATION
C---------------------------------------------------------
      NPS1 = NPSAMP(K)+1

      DO J = KSUM+1,KSUM+NPS1
      VAL(J-KSUM) = PRESP(J)
      ENDDO
      KSUM = KSUM+NPS1

      DDS = DS(K)
      DO J = 1,NPS1
      XI(J) = (J-1)*DDS
      ENDDO

      CALL INTERP1DINT(NPS1,XI,VAL,YINT)
      RINT = RINT + YINT

      ENDDO
C --- ENDDO: LOOP OVER ALL POLYGON SIDES

      TXINT = RINT / AREAPOLY

      RETURN
      END
