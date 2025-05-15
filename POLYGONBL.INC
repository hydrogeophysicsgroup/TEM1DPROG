C--------------------------------------------------------------------------
C --- THE POLYGONAL TX LOOP COMMON BLOCK
C--------------------------------------------------------------------------
      INTEGER*4 NPOLY,NRSAMP,NPSAMP(N_ARR)
      REAL*8 XPOLY(N_ARR),YPOLY(N_ARR),X0RX,Y0RX,Z0RX,DIPOLEFAC
      REAL*8 YSAMP(N_ARR),RSAMP(N_ARR),RSAMPL(N_ARR),DS(N_ARR)
      REAL*8 RSAMPMIN,RSAMPMAX,AREAPOLY,PERIMETERPOLY
      COMMON /POLYGONBL/ NPOLY,NRSAMP,NPSAMP,XPOLY,YPOLY,X0RX,Y0RX,ZORX,
     #   DS,YSAMP,RSAMP,RSAMPL,RSAMPMIN,RSAMPMAX,AREAPOLY,PERIMETERPOLY
C ----------
      PARAMETER (DIPOLEFAC = 0.10D0)
C ----------
C  NPOLY:         Number of linear segments of the polygonal Tx coil. NPOLY=0: Circular Tx coil.
C                 For NPOLY=0: Circular Tx coil, and the parameters below are moot.
C  XPOLY:         X-coordinates of the polygonal Tx coil apices.
C  YPOLY:         Y-coordinates of the polygonal Tx coil apices.
C  X0RX:          For polygonal Tx coil: The X-coordinate of the Rx dipole.
C  Y0RX:          For polygonal Tx coil: The y-coordinate of the Rx dipole.
C  Z0RX:          For polygonal Tx coil: The Z-coordinate of the Rx dipole relative to the tx plane.
C                 Z0RX is not an input paranmeter. It is calculated as Z0RX = ABS(HRX - HTX).
C  NRSAMP:        Total number of sampling points on the Tx loop perimeter.
C  NPSAMP:        Number of subdivisions on the Tx loop sides. May change from side to side.
C  YSAMP:         Y-coordinate of the Rx position in the coordinate system of an X-directed El dipole.
C  RSAMP:         Radial distance to the Rx position in the coordinate system of an X-directed El dipole.
C  RSAMPL:        log(RSAMP).
C  DS:            Length of subdivision interval for the sides.
C  RSAMPMIN:      Minimum distance for the FHT routine.
C  RSAMPMAX:      Maximum distance for the FHT routine.
C  AREAPOLY:      Area of the Tx loop
C  PERIMETERPOLY: Length of the perimeter of the Tx loop.
C ----------
