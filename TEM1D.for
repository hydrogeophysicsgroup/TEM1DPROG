CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C   S U B R O U T I N E    T E M 1 D P R O G
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   SUBROUTINE TEM1DPROG is the entry point of the program TEM1DRESP.
C   The subroutine is called with all of the parameters that are necessary
C   to calculate the TEM forwards response, optionally the derivatives.
C
C   Informationon on the parameters of the call is written to unit 21
C   which must be defined/opened in the calling program.
C
C   February 2024 / NBC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TEM1D (
     # IMLMi,NLAYi,RHONi,DEPNi,
     # IMODIPi,CHAIPi,TAUIPi,POWIPi,
     # TXAREAi,RTXRXi,IZEROPOSi,
     # ISHTX1i,ISHTX2i,ISHRX1i,ISHRX2i,HTX1i,HTX2i,HRX1i,HRX2i,
     # NPOLYi,XPOLYi,YPOLYi,X0RXi,Y0RXi,
     # IRESPTYPEi,IDERIVi,IREPi,IWCONVi,NFILTi,REPFREQi,FILTFREQi,
     # NWAVEi,TWAVEi,AWAVEi)
C ----------
      IMPLICIT INTEGER*4 (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
C ----------
      INCLUDE 'ARRAYSDIMBL.INC'
      INCLUDE 'MODELBL.INC'
      INCLUDE 'IPBL.INC'
      INCLUDE 'INSTRUMENTBL.INC'
      INCLUDE 'POLYGONBL.INC'
      INCLUDE 'RESPBL.INC'
      INCLUDE 'WAVEBL.INC'
C ----------
      REAL*8 RHONi(N_ARR),DEPNi(N_ARR)
      REAL*8 CHAIPi(N_ARR),TAUIPi(N_ARR),POWIPi(N_ARR)
      REAL*8 XPOLYi(N_ARR),YPOLYi(N_ARR)
      REAL*8 TWAVEi(N_ARR),AWAVEi(N_ARR)
C ----------
      INTEGER*4 NFILTi
      REAL*8 FILTFREQi(16)
      REAL*8 DRHI(N_ARR),DRLO(N_ARR),RESPOUTANA(N_ARR)
      REAL*8 TIMESOUT(N_ARR),RESPOUT(N_ARR),DRESPOUT(N_ARR,N_ARR)
      REAL*8 DSTPSP(N_ARR,N_ARR),RESPOUT0(N_ARR)
C ----------
      REAL*8 ZEROPOS
      EXTERNAL ZEROPOS
C ----------
      REAL*8 PI
      DATA PI / 3.141592653589793D0 /
C ----------

C=====================================================================
C --- HARDWIRING OUTPUT OPTIONS INCLUDING NUMERICAL DERIVATIVES
C=====================================================================
C...............................................................................
      IWRI = 0
C --- IWRI = [1 | 0] : [Write to output | Do not write to output].
C...............................................................................

C...............................................................................
      IWRITE21 = 1
C --- IWRITE21 = [1 | 0] : [Write to fil, unit 21 | Do not write].
C...............................................................................

C...............................................................................
C --- INUMDERIV = [1 | 0] : [2-sided numerical derivatives | Analytic derivatives].
      INUMDERIV = 0
C...............................................................................

C...............................................................................
      IF (IWRI.EQ.1) THEN
      WRITE (*,*) 'TEM1DPROG ENTERED'
      ENDIF
C...............................................................................

C============================================================
C --- TRANSFER INPUT PARAMETERS TO COMMON PARAMETERS
C============================================================
      IMLM = IMLMi

      NLAY = NLAYi
      DO I = 1,NLAY
        RHON(I)   =  RHONi(I)
        DEPN(I)   =  DEPNi(I)
      ENDDO

      IMODIP    = IMODIPi
      IF (IMODIP.GT.0) THEN
        DO I = 1,NLAY
          CHAIP(I)  = CHAIPi(I)
          TAUIP(I)  = TAUIPi(I)
          POWIP(I)  = POWIPi(I)
        ENDDO
      ENDIF

C....................................................................
      IF (IWRI.EQ.1) THEN
      WRITE (*,*) 'TEM1DPROG: MODEL PARAMETER TRANSFER DONE'
      ENDIF
C....................................................................

      IRESPTYPE = IRESPTYPEi
      IDERIV    = IDERIVi
      IZEROPOS  = IZEROPOSi
      TXAREA    = TXAREAi

      IF (INUMDERIV.EQ.1) THEN
      IDERIVOLD = IDERIV
      IDERIV = 0
      ENDIF

C....................................................................
      IF (IWRI.EQ.1) THEN
      WRITE (*,*) 'TEM1DPROG: IRESPTYPE ETC PARAMETER TRANSFER DONE'
      ENDIF
C....................................................................

      ISHTX1 = ISHTX1i
      ISHTX2 = ISHTX2i
      ISHRX1 = ISHRX1i
      ISHRX2 = ISHRX2i
      HTX1   = ABS(HTX1i)
      HTX2   = ABS(HTX2i)
      HRX1   = ABS(HRX1i)
      HRX2   = ABS(HRX2i)
      RTXRX  = RTXRXi

C....................................................................
      IF (IWRI.EQ.1) THEN
      WRITE (*,*) 'TEM1DPROG: ISH & HTR TRANSFER DONE'
      ENDIF
C....................................................................

      NFILT = NFILTi
      IF (NFILT.GT.0) THEN
        DO I = 1,NFILT
        FILTFREQ(I) = FILTFREQi(I)
        ENDDO
      ENDIF

C....................................................................
      IF (IWRI.EQ.1) THEN
      WRITE (*,*) 'TEM1DPROG: FILFREQ LOOP DONE'
      ENDIF
C....................................................................

      IREP      = IREPi
      REPFREQ   = REPFREQi
      IWCONV    = IWCONVi

      IF (IRESPTYPE.LT.2) THEN
        IREP = 0
        IWCONV = 0
        NWAVE = 0
      ENDIF

C....................................................................
      IF (IWRI.EQ.1) THEN
      WRITE (*,*) 'TEM1DPROG: INSTRUMENT PARAMETER TRANSFER DONE'
      ENDIF
C....................................................................

      NPOLY     = NPOLYi
      
      IF (NPOLY.GT.0) THEN

C....................................................................
      IF (NPOLY.LT.3) THEN
      WRITE (*,*) 'POLYGON WITH LESS THAN 3 SIDES: FALSE INPUT'
      STOP
      ENDIF
C....................................................................

        DO I = 1,NPOLY
          XPOLY(I)  = XPOLYi(I)
          YPOLY(I)  = YPOLYi(I)
        ENDDO
      X0RX    = X0RXi
      Y0RX    = Y0RXi
      Z0RX    =          ABS(HRX1-HTX1)
      Z0RX    = MIN(Z0RX,ABS(HRX1-HTX2))
      Z0RX    = MIN(Z0RX,ABS(HRX2-HTX1))
      Z0RX    = MIN(Z0RX,ABS(HRX2-HTX2))
      ENDIF

C....................................................................
      IF (IWRI.EQ.1) THEN
      WRITE (*,*) 'TEM1DPROG: POLYGONAL PARAMETER TRANSFER DONE'
      ENDIF
C....................................................................

      NWAVE  = NWAVEi
      IF (NWAVE.GT.0) THEN
      DO I = 1,NWAVE
        TWAVE(I)     = TWAVEi(I)
        AWAVE(I)     = AWAVEi(I)
      ENDDO
      ENDIF

C....................................................................
      IF (IWRI.EQ.1) THEN
      WRITE (*,*) 'TEM1DPROG: WAVEFORM PART OF DATA TRANSFER DONE'
      ENDIF
C....................................................................

C-----------------------------------------------------------------------
C --- CALCULATIONS OF OTHER COMMON BLOCK PARAMETERS
C-----------------------------------------------------------------------
      DO I = 1,NLAY
      SIGN(I) = 1.D0/RHON(I)
      ENDDO

      IF (NLAY.GT.1) THEN
        DO I = 1,NLAY-1
        THKN(I) = DEPN(I+1)-DEPN(I)
        ENDDO
      ENDIF

C....................................................................
      IF (IWRI.EQ.1) THEN
      WRITE (*,*) 'TEM1DPROG: SIGN AND THK TRANSFERRED'
      ENDIF
C....................................................................

C--------------------------------------------------------------------------
C --- CALCULATIONS WHEN THE TX IS APPROXIMATED WITH A CIRCULAR COIL
C--------------------------------------------------------------------------
      TXRAD = SQRT(TXAREA / PI)

      IF (IZEROPOS.GT.0) THEN
        TXA = TXAREA
        HEIGHT = HRX1-HTX1
        EQRXPOS = ZEROPOS(TXA,HEIGHT)
      ENDIF

C--------------------------------------------------------------------------
C --- IDENTIFICATION OF A CENTRAL LOOP CONFIGURATION
C--------------------------------------------------------------------------
      ICEN = 0
      IF (RTXRX.LT.0.01D0 .AND. NPOLY.EQ.0) THEN
      ICEN = 1
      ENDIF

C--------------------------------------------------------------------------
C --- FIND SECOND DERIVATIVE OF THE WAVEFORM
C--------------------------------------------------------------------------
      DO I = 1,NWAVE-1
      SLOPE(I) = (AWAVE(I+1)-AWAVE(I))/(TWAVE(I+1)-TWAVE(I))
      ENDDO

C....................................................................
      IF (IWRI.EQ.1) THEN
      WRITE (*,*) 'TEM1DPROG: 1ST DERIVATIVE OF WAVEFORM DONE'
      ENDIF
C....................................................................

      D2WDT2(1) = SLOPE(1)
      DO I = 2,NWAVE-1
       D2WDT2(I) = SLOPE(I)-SLOPE(I-1)
      ENDDO
      D2WDT2(NWAVE) = -SLOPE(NWAVE-1)

C....................................................................
      IF (IWRI.EQ.1) THEN
      WRITE (*,*) 'TEM1DPROG: 2ND DERIVATIVE OF WAVEFORM DONE'
      ENDIF
C....................................................................

C....................................................................
      IF (IWRI.EQ.1) THEN
      WRITE (*,*) 'TEM1DPROG: ALL DATA PROCESSING DONE'
      ENDIF
C....................................................................

C**************************************************************************************
C=====================================================================================
C --- IF POLYGONAL LOOP, CALCULATE THE PARAMETERS NEEDED IN THE FURTHER PROCESSING.
C=====================================================================================
C**************************************************************************************
      IF (NPOLY .GT. 0) THEN

C--------------------------------------------------------------------
C --- This line ensures that the kernel function will be correct.
C--------------------------------------------------------------------
      ICEN = 1

C=====================================================================================
C --- CALCULATE TX SAMPLING POINTS Y-DISTANCES & RADIAL DISTANCES
C --- TO SAMPLING POINTS ON THE POLYGONAL SIDES.
C --- Y-DISTANCES ARE THE SAME FOR EVERY SAMPLING POINT ON A POLYGON SIDE SO: Y(1:NSIDES).
C --- RADIAL DISTANCES DIFFER FOR ALL NSIDES*(NPSAMP+1) POSIITONS,
C --- BUT THEY ARE THE SAME FOR EVERY MODEL & EVERY DELAY TIME.
C=====================================================================================
      II = 0
      RSAMPMIN = 1.D30
      RSAMPMAX = 0.D0
      PERIMETERPOLY = 0.D0
      NRSAMP = 0

C=====================================================================================
C --- LOOP OVER THE POLYGON SIDES (REPEAT THE LAST POINT)
C=====================================================================================
      XPOLY(NPOLY+1) = XPOLY(1)
      YPOLY(NPOLY+1) = YPOLY(1)

      DO K = 1,NPOLY

      X1 = XPOLY(K)
      X2 = XPOLY(K+1)
      Y1 = YPOLY(K)
      Y2 = YPOLY(K+1)
      SIDEL = SQRT((X2-X1)*(X2-X1)+(Y2-Y1)*(Y2-Y1))
      PERIMETERPOLY = PERIMETERPOLY + SIDEL

C-----------------------------------------------------------------------------
C --- For this side,find the minimum distance to the Rx.
C --- It is assumed that the Rx does not lie on one of tye polygon sides.
C-----------------------------------------------------------------------------
      DMIN = DISTMIN(X1,Y1,X2,Y2,X0RX,Y0RX)
      DMINA = SQRT(Z0RX*Z0RX + DMIN*DMIN)

C---------------------------------------------------------------
C --- CHOOSE THE SAMPLING DENSITY ON THE SIDES
C --- These lines select the local splineintegration formula
C---------------------------------------------------------------
      DD = DIPOLEFAC * DMINA
      NPSAMP(K) = INT(SIDEL/DD)+1
      DS(K) = SIDEL/FLOAT(NPSAMP(K)) 
      NRSAMP = NRSAMP + (NPSAMP(K)+1)

C---------------------------------------------------------
C --- COORDINATES OF THE UNIT VECTOR ALONG THE SIDE
C---------------------------------------------------------
      E1 = (X2-X1)/SIDEL
      E2 = (Y2-Y1)/SIDEL

      SP    =  E1*(X0RX-X1)+E2*(Y0RX-Y1)
      YSA   = -E2*(X0RX-X1)+E1*(Y0RX-Y1)

C=====================================================================================
C --- LOOP OVER THE SAMPLING POINTS ON THE PRESENT SIDE
C=====================================================================================
      DO I = 1,NPSAMP(K)+1
      II = II+1
      S = (I-1)*DS(K)
      XI = SP-S
      YSAMP(II) = YSA
      RSAMP(II) = SQRT(XI*XI+YSA*YSA)
      RSAMPL(II) =  LOG(RSAMP(II))

      RSAMPMIN = MIN(RSAMPMIN,RSAMP(II))
      RSAMPMAX = MAX(RSAMPMAX,RSAMP(II))

      ENDDO
C --- ENDDO: LOOP OVER ALL SAMPLING POINTS ON THE POLYGON SIDE

      ENDDO
C --- ENDDO: LOOP OVER ALL POLYGON SIDES

C------------------------------------------
C --- FIND AREA OF POLYGON
C------------------------------------------
      AREAPOLY = 0.D0
      DO I = 1,NPOLY
      AREAPOLY = AREAPOLY+(XPOLY(I)*YPOLY(I+1)-XPOLY(I+1)*YPOLY(I))
      ENDDO
      AREAPOLY = 0.5D0*AREAPOLY

C....................................................................
      IF (IWRI.EQ.1) THEN
      WRITE (*,*) 'TEM1DPROG: POLYGON PROCESSING DONE'
      ENDIF
C....................................................................

      ENDIF
C --- ENDIF: POLYGONAL TX LOOP
C=====================================================================================


C===================================================================
C --- CALL THE RESPONSE ROUTINE
C===================================================================
      T1 = STIMER(DUMMY)

      IF (IMODIP.EQ.0 .AND. NPOLY.EQ.0) THEN
      WRITE (*,*) 'TEM1DPROG: BEFORE CALLING TEMRESP'
        CALL TEMRESP (NTOUT,TIMESOUT,RESPOUT,DRESPOUT)
      WRITE (*,*) 'TEM1DPROG: AFTER CALLING TEMRESP'
      ENDIF

      IF (IMODIP.EQ.1 .AND. NPOLY.EQ.0) THEN
      WRITE (*,*) 'TEM1DPROG: BEFORE CALLING TEMRESPIP'
        CALL TEMRESPIP  (NTOUT,TIMESOUT,RESPOUT,DRESPOUT)
      WRITE (*,*) 'TEM1DPROG: AFTER CALLING TEMRESPIP'
      ENDIF

      IF (IMODIP.EQ.0 .AND. NPOLY.GT.0) THEN
      WRITE (*,*) 'TEM1DPROG: BEFORE CALLING TEMRESPPOLY'
        CALL TEMRESPPOLY (NTOUT,TIMESOUT,RESPOUT,DRESPOUT)
      WRITE (*,*) 'TEM1DPROG: AFTER CALLING TEMRESPPOLY'
      ENDIF

      IF (IMODIP.EQ.1 .AND. NPOLY.GT.0) THEN
      WRITE (*,*) 'TEM1DPROG: BEFORE CALLING TEMRESPPOLYIP'
        CALL TEMRESPPOLYIP (NTOUT,TIMESOUT,RESPOUT,DRESPOUT)
      WRITE (*,*) 'TEM1DPROG: AFTER CALLING TEMRESPPOLYIPP'
      ENDIF

C....................................................................
      T2 = STIMER(DUMMY)
      IF (IWRITE21.EQ.1) THEN
      WRITE (21,*) 'CALCULATION TIME: ',T2-T1
      ENDIF
C....................................................................

C--------------------------------------------------------------------------------
C --- SAVE RESPONSE (MIGHT BE MODIFIED BELOW)
C--------------------------------------------------------------------------------
      DO I = 1,NTOUT
      RESPOUT0(I) = RESPOUT(I)
      ENDDO

C********************************************************************
C===================================================================
C --- CALCULATE THE NUMERICAL DERIVAITVES
C===================================================================
C********************************************************************
      IF (INUMDERIV.EQ.1) THEN

C------------------------------------------
C --- DERIVATIVES WRT CONDUCTIVITIES
C------------------------------------------
      DO J=1,NLAY
      SIG0 = SIGN(J)
      SIGN(J) = SIG0 * 1.05D0
      RHON(J)=1.D0/SIGN(J)
      CALL TEMRESP (NTOUT,TIMESOUT,RESPOUT,DRESPOUT)

        DO I = 1,NTOUT
        DRHI(I) = RESPOUT(I)
        ENDDO

      SIGN(J) = SIG0 * 0.95D0
      RHON(J)=1.D0/SIGN(J)

      CALL TEMRESP (NTOUT,TIMESOUT,RESPOUT,DRESPOUT)

        DO I = 1,NTOUT
        DRLO(I) = RESPOUT(I)
        ENDDO

      SIGN(J) = SIG0

        DO I = 1,NTOUT
        DRESPOUT(I,J) = (DRHI(I)-DRLO(I)) / (0.1D0*SIG0)
        ENDDO

      ENDDO

      WRITE (*,*) 'NUMDERIV: COND DONE'

C------------------------------------------
C --- DERIVATIVES WRT THICKNESSES
C------------------------------------------
      IF (IMLM.EQ.0)  THEN

      DO J=1,NLAY-1

      THK0 = THKN(J)
      THKN(J) = THK0 * 1.05D0

        DEPN(1) = 0.D0
        DO K=1,NLAY-1
        DEPN(K+1) = DEPN(K)+THKN(K)
        ENDDO

      CALL TEMRESP(NTOUT,TIMESOUT,RESPOUT,DRESPOUT)

        DO I = 1,NTOUT
        DRHI(I) = RESPOUT(I)
        ENDDO

      THKN(J) = THK0 * 0.95D0

        DEPN(1) = 0.D0
        DO K=1,NLAY-1
        DEPN(K+1) = DEPN(K)+THKN(K)
        ENDDO

      CALL TEMRESP(NTOUT,TIMESOUT,RESPOUT,DRESPOUT)

        DO I = 1,NTOUT
        DRLO(I) = RESPOUT(I)
        ENDDO

      THKN(J) = THK0

        DEPN(1) = 0.D0
        DO K=1,NLAY-1
        DEPN(K+1) = DEPN(K)+THKN(K)
        ENDDO

      DO I = 1,NTOUT
      DRESPOUT(I,NLAY+J) = (DRHI(I)-DRLO(I)) / (0.1D0*THK0)
      ENDDO

      ENDDO

      ENDIF

      WRITE (*,*) 'NUMDERIV: THKS DONE'

C---------------------------------------------------------------------
C --- DERIVATIVES WRT HTX
C --- PERTURBATION MUST BE ADDITIVE TO MAINTAIN CONFIGURATION.
C --- ALL HTX/HRX PARAMETERSS MUST BE PERTURBED.
C --- THE ISHTX1 TYPE PARAMETERS WILL TAKE CARE OF THE REST.
C---------------------------------------------------------------------
      HTX10 = HTX1
      HRX10 = HRX1
      HTX20 = HTX2
      HRX20 = HRX2
      DELH = 0.1D0*HTX1

      HTX1 = HTX10 + DELH
      HRX1 = HRX10 + DELH
      HTX2 = HTX20 + DELH
      HRX2 = HRX20 + DELH

      CALL TEMRESP(NTOUT,TIMESOUT,RESPOUT,DRESPOUT)

        DO I = 1,NTOUT
        DRHI(I) = RESPOUT(I)
        ENDDO

      HTX1 = HTX10 - DELH
      HRX1 = HRX10 - DELH
      HTX2 = HTX20 - DELH
      HRX2 = HRX20 - DELH

      CALL TEMRESP(NTOUT,TIMESOUT,RESPOUT,DRESPOUT)

        DO I = 1,NTOUT
        DRLO(I) = RESPOUT(I)
        ENDDO

      HTX1 = HTX10
      HRX1 = HRX10
      HTX2 = HTX20
      HRX2 = HRX20

      WRITE (*,*) 'NUMDERIV: HTX DONE - NOT YET STORED'

      IF (IMLM.EQ.1) THEN

      DO I = 1,NTOUT
      DRESPOUT(I,NLAY+1) = (DRHI(I)-DRLO(I)) / (2.D0*DELH)
      ENDDO

      ELSEIF (IMLM.EQ.0) THEN

      DO I = 1,NTOUT
      DRESPOUT(I,2*NLAY) = (DRHI(I)-DRLO(I)) / (2.D0*DELH)
      ENDDO

      ENDIF

      WRITE (*,*) 'NUMDERIV: HTX DONE'

      IDERIV = IDERIVOLD

      ENDIF
C --- ENDIF: NUMERICAL DERIVATIVES ARE DONE 

C===========================================================================
C --- OUTPUT FROM THE PROGRAM
C --- THE ONLY OUTPUT FROM THE PROGRAM ITSELF IS RESPONSES AND DERIVATIVES
C===========================================================================

C***************************************************************
C============================================================
C --- WRITE SETTINGS TO OUTPUT
C============================================================
C***************************************************************
      IF (IWRITE21.EQ.1) THEN
      
      WRITE (21,'(A)') ' '
      WRITE (21,'(A)') ' '
      WRITE (21,'(A)') '=============================================='
      WRITE (21,'(A)') '  >> THE MODEL BLOCK <<'
      WRITE (21,'(A)') '=============================================='
      WRITE (21,'(A)') ' IMLM, NLAY:'
      WRITE (21,'(3X,I1,3X,I2)')    IMLM, NLAY
      WRITE (21,'(A)') ' RHON, DEPN:'
      WRITE (21,3002)  (RHON(J),J=1,NLAY)
      WRITE (21,3002)  (DEPN(J),J=1,NLAY)
      WRITE (21,'(A)') '=============================================='

C***************************************************************
      WRITE (21,'(A)') ' '
      WRITE (21,'(A)') '=============================================='
      WRITE (21,'(A)') '  >> THE INSTRUMENT BLOCK <<'
      WRITE (21,'(A)') '=============================================='
      WRITE (21,'(A)') ' ICEN, IZEROPOS:'
      WRITE (21,'(3X,I1,2X,4X,I1)')    ICEN, IZEROPOS
      WRITE (21,'(A)') ' TXAREA ,  TXRAD ,   RTXRX , EQRXPOS:'
      WRITE (21,'(4(1X,F7.3,1X))')
     #            TXAREA,TXRAD,RTXRX,EQRXPOS
      WRITE (21,'(A)') '  HTX1 , HTX2 , HRX1 , HRX2:'
      WRITE (21,'(4(1X,F5.2,1X))') HTX1,HTX2,HRX1,HRX2
      WRITE (21,'(A)') ' ISHTX1, ISHTX2, ISHRX1, ISHRX2:'
      WRITE (21,'(4(3X,I2,3X))') ISHTX1, ISHTX2, ISHRX1, ISHRX2
      WRITE (21,'(A)') '=============================================='

C***************************************************************
      IF (IMODIP.GT.0) THEN
      WRITE (21,'(A)') ' '
      WRITE (21,'(A)') '=============================================='
      WRITE (21,'(A)') '  >> THE IP BLOCK <<'
      WRITE (21,'(A)') '=============================================='
      WRITE (21,'(A)') ' CHARGEABILITIES:'
      WRITE (21,3002)    (CHAIP(J),J=1,NLAY)
      WRITE (21,'(A)') ' TIME CONSTANTS:'
      WRITE (21,3002)    (TAUIP(J),J=1,NLAY)
      WRITE (21,'(A)') ' POWERS :'
      WRITE (21,3002)    (POWIP(J),J=1,NLAY)
      WRITE (21,'(A)') '=============================================='
      ENDIF
C***************************************************************
      IF (NPOLY.GT.0) THEN
      WRITE (21,'(A)') ' '
      WRITE (21,'(A)') '=============================================='
      WRITE (21,'(A)') '  >> THE POLYGONAL TX LOOP BLOCK <<'
      WRITE (21,'(A)') '=============================================='
      WRITE (21,'(A)') ' NPOLY:'
      WRITE (21,3000)    NPOLY
      WRITE (21,'(A)') ' X-COORDINATES OF APICES:'
      WRITE (21,3002)    (XPOLY(J),J=1,NPOLY)
      WRITE (21,'(A)') ' Y-COORDINATES OF APICES:'
      WRITE (21,3002)    (YPOLY(J),J=1,NPOLY)
      WRITE (21,'(A)') ' DIPOLEFAC & X-, Y-, AND Z-COORDINATES OF RX:'
      WRITE (21,3002)    DIPOLEFAC,X0RX,Y0RX,Z0RX
      WRITE (21,'(A)') '=============================================='
      ENDIF
C***************************************************************
      WRITE (21,'(A)') ' '
      WRITE (21,'(A)') '=============================================='
      WRITE (21,'(A)') '  >> THE RESPONSE BLOCK <<'
      WRITE (21,'(A)') '=============================================='
      WRITE (21,'(A)') ' IRESPTYPE,IDERIV,IREP,IWCONV,NFILT:'
      WRITE (21,3000)    IRESPTYPE,IDERIV,IREP,IWCONV,NFILT
      WRITE (21,'(A)') ' REPFREQ,(FILTFREQ(J),J=1,NFILT):'
      WRITE (21,3002)    REPFREQ,(FILTFREQ(J),J=1,NFILT)
      WRITE (21,'(A)') '=============================================='

 3000 FORMAT (12(2X,I3))
 3001 FORMAT (30(2X,1PE11.4))
 3002 FORMAT (30(2X,F10.2))
C***************************************************************

      WRITE (21,'(A)') ' '
      WRITE (21,'(A)') ' '
      WRITE (21,'(A)')
     # '=============================================================='
      WRITE (21,'(A)') ' NUMBER   TIME       RESP         DERIV'
      WRITE (21,'(A)')
     # '=============================================================='

      IF (IDERIV.EQ.0) THEN
      DO I = 1,NTOUT
        WRITE (21,3003) I,TIMESOUT(I),RESPOUT0(I)
      ENDDO
      ENDIF

 3003 FORMAT (I3,20(2X,1PE11.4))

      IF (IDERIV.GT.0) THEN

        IF (IMLM.EQ.1) THEN
          NPARM = NLAY+1
        ELSEIF (IMLM.EQ.0) THEN
          NPARM = NLAY + (NLAY-1) + 1
        ENDIF

      DO I = 1,NTOUT
        WRITE (21,3003)
     #   I,TIMESOUT(I),RESPOUT0(I),(DRESPOUT(I,J),J=1,NPARM)
      ENDDO

      WRITE (21,'(A)')
     # '=============================================================='

      ENDIF

      ENDIF
C --- ENDIF: WRITING TO UNIT 21

C========================================================================
C --- WRITE TO OUTPUT FILE: FORWRITE, IN SIMPLE FORMAT, UNIT=8
C========================================================================
      WRITE (*,*) 'WRITING RESPONSES TO FORWRITE'

      WRITE (8,3004) (TIMESOUT(I),I=1,NTOUT)
      WRITE (8,3004) (RESPOUT0(I),I=1,NTOUT)

      IF (IDERIV.GT.0) THEN
        DO J = 1,NPARM
        WRITE (8,3004) (DRESPOUT(I,J),I=1,NTOUT)
        ENDDO
      ENDIF

 3004 FORMAT (128(2X,1PE11.4))
C========================================================================

 1010 FORMAT ( 2(2X,1PE11.4))
 1011 FORMAT (35(2X,1PE11.4))

      RETURN
      END















