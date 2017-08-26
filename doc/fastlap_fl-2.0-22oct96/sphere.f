C-----------------------------------------------------------------------
C
C     COPYRIGHT (C) 1992  MASSACHUSETTS INSTITUTE OF TECHNOLOGY
C
C-----------------------------------------------------------------------
C
C     Main program 'sphere'
C
C-----------------------------------------------------------------------
C
C This program generates a geometric data file for a sphere for 
C input to the FastLap driver: driver.c.  To modify the output 
C for the FastLap driver: driver.f, use c2fdata.
C
C
C Boundary conditions are for a sphere translating in an
C infinite fluid for which the solution is known in closed form.
C
C     Written by T. Korsmeyer, March 1992
C
C-------------------------------------------------------------------------
      PROGRAM SPHERE
      IMPLICIT NONE
      INTEGER MAX_NODE, MAX_PANEL, NODE_TOTAL, PANEL_TOTAL
      PARAMETER (MAX_NODE = 16384)
      PARAMETER (MAX_PANEL = 16384)
      CHARACTER*72 HEAD
      INTEGER NUM_SPHERES
      DOUBLE PRECISION X(3,MAX_NODE)
      DOUBLE PRECISION XC(3,4)
      DOUBLE PRECISION XVER(4,3), XVP(4),YVP(4),XCT(3),DIR(3,3),
     $     SS(15),SID(4),VEL(3),VROT(4,3),DIAG(2)
      INTEGER IDIAG
      DOUBLE PRECISION POTENTIAL, DBYDNPOTEN
      INTEGER FSCON(4,MAX_PANEL), MSCON(4,MAX_PANEL), FSPAN, MSPAN
      INTEGER NFS, NMS, NWL
      INTEGER LFS(MAX_NODE), LMS(MAX_NODE), LWL(MAX_NODE)
      INTEGER INDEX_FS, INDEX_MS, M_FIRST, COUNT_FS, COUNT_MS,
     $     COUNT_WL, NA, NM, OFFSET, I, J, K, N
      DOUBLE PRECISION THETA, PHI, D_THETA, D_PHI, R, D_R, PI, RAD,
     $     CLOSENESS
      DOUBLE PRECISION  PHI_NORM, PHIN_NORM, MAX_PHI, MAX_PHIN
      INTEGER TYPE
      DOUBLE PRECISION DUM
      DATA DUM / 1.0D00 /  
      PARAMETER(RAD = 1.0)
      DATA PI / 3.1415927 /
C
C Define output device number and name.
C
      CHARACTER*20 FASTIN
      DATA FASTIN /'fastlap.in'/
C
C Output heading, and query on output files.
C
      WRITE (6,935)
 935  FORMAT(1X,/,
     &' A new panelization will be generated for a sphere.'///)
 1    CONTINUE
C
C Query on panelization.
C
 2    CONTINUE
      WRITE (6,941)
 941  FORMAT(1X,//5X,71(1H-),//
     &  8X,' The number of panels is specified by two inputs:'//
     &  8X,' Number from pole to pole:   ',$)
      READ (5,*) NA
      NA = NA / 2
      WRITE (6,942)
 942  FORMAT(1X,/
     &  8X,' Number around any parallel of latitude:   ',$)
      READ (5,*) NM
      OPEN (UNIT=14,FILE=FASTIN,STATUS='UNKNOWN')
C
C Generate nodes on the upper and lower surfaces.  We do them both
C at once because they have common x,y coordinates, only their elevation
C is different. Bottom dead center and top dead center are common to 
C every meridian,  so don't repeat them.
C
      NFS = (NA-1) * NM + 1
      NMS = NFS
      NWL = NM
      FSPAN = NM * NA
      MSPAN = NM * NA
      NODE_TOTAL = NFS + NWL + NMS
      PANEL_TOTAL = FSPAN + MSPAN
      IF (NODE_TOTAL .GT. MAX_NODE) THEN
         WRITE(6,'('' Node dimension of '',I5, 
     $        '' is exceeded, try again.'')') MAX_NODE
         GO TO 2
      END IF
      IF (PANEL_TOTAL .GT. MAX_PANEL) THEN
         WRITE(6,'('' Panel dimension of '',I5,
     $        '' is exceeded, try again.'')') MAX_PANEL
         GO TO 2
      END IF
      D_R = RAD / NA
      D_THETA = PI / (2. * NA)
      D_PHI = (2.0 * PI) / NM
      M_FIRST = NFS + NWL + 1
      INDEX_FS = 1
      INDEX_MS = M_FIRST
      COUNT_FS = 1
      COUNT_MS = 1
      COUNT_WL = 1
      X(1,INDEX_FS) = 0.0
      X(2,INDEX_FS) = 0.0
      X(3,INDEX_FS) = RAD
      X(1,INDEX_MS) =  0.0
      X(2,INDEX_MS) =  0.0
      X(3,INDEX_MS) = -RAD
      LFS(COUNT_FS) = INDEX_FS
      LMS(COUNT_MS) = INDEX_MS
      COUNT_FS = COUNT_FS + 1
      COUNT_MS = COUNT_MS + 1
C
C Add these nodes to the node lists 
C
C
C Loop along the parallels.
C
      DO J = 0, NM-1
         PHI = J * D_PHI
C
C Loop along the meridians, but don't do the equator.
C
         DO I = 1, NA-1
            R = I * D_R
            THETA = PI - (I * D_THETA)
            INDEX_FS = INDEX_FS + 1
            INDEX_MS = INDEX_MS + 1
            X(1,INDEX_MS) = RAD * SIN(THETA) * COS(PHI)
            X(2,INDEX_MS) = RAD * SIN(THETA) * SIN(PHI)
            X(3,INDEX_MS) = RAD * COS(THETA)
            X(1,INDEX_FS) =  X(1,INDEX_MS)
            X(2,INDEX_FS) =  X(2,INDEX_MS)
            X(3,INDEX_FS) = -X(3,INDEX_MS)
            LFS(COUNT_FS) = INDEX_FS
            LMS(COUNT_MS) = INDEX_MS
            COUNT_FS = COUNT_FS + 1
            COUNT_MS = COUNT_MS + 1
         END DO
C
C Now do the waterline.
C
         THETA = PI - (I * D_THETA)
         INDEX_FS = INDEX_FS + 1
         X(1,INDEX_FS) = RAD * SIN(THETA) * COS(PHI)
         X(2,INDEX_FS) = RAD * SIN(THETA) * SIN(PHI)
         X(3,INDEX_FS) = 0.0
         LWL(COUNT_WL) = INDEX_FS
         COUNT_WL = COUNT_WL + 1
      END DO
C
C Generate the connectivity of the panels.  Vertex numbering is 
C clockwise from inside the fluid domain.
C
      INDEX_FS = 0
      INDEX_MS = 0
      DO J = 1, NM
C
C First do the top and bottom dead center panels which are special.
C
         INDEX_FS = INDEX_FS + 1
         INDEX_MS = INDEX_MS + 1
         FSCON(4,INDEX_FS) = 1
         FSCON(3,INDEX_FS) = 2 + NA * (J-1) 
         FSCON(2,INDEX_FS) = 2 + NA * J
         FSCON(1,INDEX_FS) = 1
         MSCON(1,INDEX_MS) = M_FIRST
         MSCON(2,INDEX_MS) = M_FIRST + 1 + (NA-1) * (J-1) 
         MSCON(3,INDEX_MS) = M_FIRST + 1 + (NA-1) * J
         MSCON(4,INDEX_MS) = M_FIRST
C
C Now do the rest EXCEPT for the panels sharing the equator.
C
         DO I = 2, NA-1
            INDEX_FS = INDEX_FS + 1
            INDEX_MS = INDEX_MS + 1
            FSCON(4,INDEX_FS) = NA * (J-1) + I
            FSCON(1,INDEX_FS) = NA * J + I
            FSCON(3,INDEX_FS) = FSCON(4,INDEX_MS) + 1
            FSCON(2,INDEX_FS) = FSCON(1,INDEX_MS) + 1
            MSCON(1,INDEX_MS) = (NA-1) * (J-1) + (M_FIRST-1) + I
            MSCON(4,INDEX_MS) = (NA-1) * J + (M_FIRST-1) + I
            MSCON(3,INDEX_MS) = MSCON(4,INDEX_MS) + 1
            MSCON(2,INDEX_MS) = MSCON(1,INDEX_MS) + 1
         END DO 
C
C Now do the equator panels.
C
         INDEX_FS = INDEX_FS + 1
         INDEX_MS = INDEX_MS + 1
         FSCON(4,INDEX_FS) = NA * (J-1) + NA
         FSCON(1,INDEX_FS) = NA * J + NA
         FSCON(3,INDEX_FS) = FSCON(4,INDEX_MS) + 1
         FSCON(2,INDEX_FS) = FSCON(1,INDEX_MS) + 1
         MSCON(1,INDEX_MS) = (NA-1) * (J-1) + (M_FIRST-1) + NA
         MSCON(4,INDEX_MS) = (NA-1) * J + (M_FIRST-1) + NA
         MSCON(2,INDEX_MS) = 1 + J * NA
         MSCON(3,INDEX_MS) = 1 + (J+1) * NA
      END DO
C
C Now correct the last strip of panels on which the node numbers will exceed
C the max node number unless they are wrapped into the first strip of nodes.
C
      OFFSET = NM*NA-NA
      DO I = 1,NA
         FSCON(2,OFFSET+I) = FSCON(3,I)
         FSCON(1,OFFSET+I) = FSCON(4,I)
         MSCON(4,OFFSET+I) = MSCON(1,I)
         MSCON(3,OFFSET+I) = MSCON(2,I)
      END DO
      WRITE (6,903)
 903  FORMAT (1X,/,' Enter header line for gdf file: ',$)
      READ (5,902) HEAD
 902  FORMAT (A72)
      WRITE (14,905) HEAD
 905  FORMAT ('0  ',A72)
      WRITE(*,'('' What type on trailing half? 0=Dirichlet''
     $     '' 1=Neumann.'')')
      READ(5,*) type
      MAX_PHI = 0.0
      MAX_PHIN = 0.0
      PHI_NORM = 0.0
      PHIN_NORM = 0.0
      DO I = 1, FSPAN
         DO J = 1,4
            XC(1,J) = X(1,FSCON(J,I))
            XC(2,J) = X(2,FSCON(J,I))
            XC(3,J) = X(3,FSCON(J,I))
         END DO
c     
c Find the centroid of the panel.
c     
         CALL PANEL(XC,XCT,DIR,XVP,YVP,SS,SID)
         POTENTIAL = -0.5 * XCT(3)
         DBYDNPOTEN = -XCT(3)
         IF (type .eq. 1) THEN 
C     
C This is a panel on which PHI is the unknown.
C
            PHI_NORM = PHI_NORM + ABS(POTENTIAL)
            MAX_PHI = MAX(ABS(POTENTIAL),MAX_PHI)
         ELSE 
C
C This is a panel on which PHIN is the unknown.
C
            PHIN_NORM = PHIN_NORM + ABS(DBYDNPOTEN)
            MAX_PHIN = MAX(ABS(DBYDNPOTEN),MAX_PHIN)
         END IF
         IF (FSCON(1,I) .EQ. FSCON(4,I)) THEN
            WRITE(14,951) ((XC(K,J),K=1,3),J=1,3),
     $           POTENTIAL,DBYDNPOTEN,TYPE
         ELSE
            WRITE(14,950) ((XC(K,J),K=1,3),J=1,4),
     $           POTENTIAL,DBYDNPOTEN,TYPE
         END IF
      END DO
      write(*,'('' What type on leading half? 0=Dirichlet''
     $     '' 1=Neumann.'')')
      read(5,*) type
      DO I = 1, MSPAN
         DO J = 1,4
            XC(1,J) = X(1,MSCON(J,I))
            XC(2,J) = X(2,MSCON(J,I))
            XC(3,J) = X(3,MSCON(J,I))
         END DO
c
c Find the centroid of the panel.
c
         CALL PANEL(XC,XCT,DIR,XVP,YVP,SS,SID)
         POTENTIAL = -0.5 * XCT(3)
         DBYDNPOTEN = -XCT(3)
         IF (MSCON(1,I) .EQ. MSCON(4,I)) THEN
            WRITE(14,951) ((XC(K,J),K=1,3),J=1,3),
     $           POTENTIAL,DBYDNPOTEN,TYPE
         ELSE
            WRITE(14,950) ((XC(K,J),K=1,3),J=1,4),
     $           POTENTIAL,DBYDNPOTEN,TYPE
         END IF
      END DO
      WRITE(*,'(''Number of panels in fastlap.in:'',I5)') 
     $     FSPAN + MSPAN
 950  FORMAT('Q  ',14f18.14,i3)
 951  FORMAT('T  ',11f18.14,i3)
 952  FORMAT(12E12.4)
      STOP
      END
C
C
C
      SUBROUTINE PANEL(V,XCTP,DIRP,XVP,YVP,S,SIDE)
C-----------------------------------------------------------------------
C
C     COPYRIGHT (C) 1988,1991 MASSACHUSETTS INSTITUTE OF TECHNOLOGY
C
C-----------------------------------------------------------------------
C
C     DESCRIPTION : Subroutine PANEL evaluates geometrical data for a
C     plane quadrilateral approximating in a mean sense the curvilinear
C     surface defined by four input vertices (not lying on a plane).
C     Adjacent vertices are connected by straight segments. The plane of
C     the quadrilateral is defined by the segment midpoints and its
C     vertices by the projection of the original vetrices on that plane.
C
C-----------------------------------------------------------------------
C
C     Input arguments (in alphabetical order) :
C
C     Name              Description 
C     ------------------------------------------------------------------
C     V(I,J)            Body-system coordinates (x,y,z)<=>J=1,2,3
C                       of I-th vertex, I=1,..,4 of current panel
C
C     Output Arguments (in alphabetical order) :
C
C     Symbolic       Description 
C     name
C     ------------------------------------------------------------------
C     DIRP(I,J)      Direction cosines of axes of panel coordinate
C                      system for current panel, relative to axes of
C                      body coordinate system
C                      I=1,2,3 <=> (x,y,z) axes of panel system,
C                      J=1,2,3 <=> (x,y,z) axes of body system
C     S(K)           Moments of inertia of current quadrilateral
C                      surface w.r.t. panel coordinate system, K=1,..,15
C                      Defined as in subroutine GEOM
C     XVP(I),YVP(I)    Panel coordinates of I-th vertex of current
C                      quadrilateral, I=1,..,4
C     SIDE(I)        Length of four sides (I=1,..,4) of current
C                      quadrilateral numbered couter-clockwise and
C                      defined so that SIDE(1) connects vertex 1 to
C                      vertex 2
C     XCTP(I)        Body-system coordinates of centroid of N-th
C                      quadrilateral, (x,y,z) <=> I=1,2,3
C
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V(3,4),XCTP(*),XVP(*),YVP(*),DIRP(3,*),S(*),SIDE(*)
C-----------------------------------------------------------------------
C     Allocate space for local arrays
C-----------------------------------------------------------------------
      DIMENSION XORP(3),UL(3),VL(3),N1(4),DX(4),DY(4),
     &          IX(4),IY(4),XP(6,4),YP(6,4),CS(15),K1(5)
      DATA  ONE /1.E+00/, HALF / 0.5E+00/, QUART / 0.25E+00/,
     &      ZERO/ 0.0E+0/, GND/ 0.5773503E+00 /,
     &      TOLS / 1.E-10/, TOLT / 1.E-12 /, THIRD / 0.3333333E+00/,
     &      CS / 1., 1., 1.5, 1.5, 3.75, 1., 3.0, 1.5, 7.5, 1.5, 1.5,
     &           3.75, 1.5, 7.5, 3.75 /,
     &      IX / 1, 1, -1, -1/, IY / 1, -1, -1, 1/,
     &      K1 / 1, 6, 10, 13, 15 /,  N1 / 2, 3, 4, 1 /
      DO 10 J=1,3
        XORP(J)=((V(J,1)+V(J,2))+(V(J,3)+V(J,4)))*QUART
        UL(J)=(V(J,1)+V(J,2))-(V(J,3)+V(J,4))
        VL(J)=(V(J,1)+V(J,4))-(V(J,2)+V(J,3))
10    CONTINUE
C-----------------------------------------------------------------------
C     Evaluation of direction cosines DIRP(3,3) of local coordinate
C     system
C-----------------------------------------------------------------------
      ULSQ=SQRT(UL(1)*UL(1)+UL(2)*UL(2)+UL(3)*UL(3))
      DIRP(1,1)=UL(1)/ULSQ
      DIRP(1,2)=UL(2)/ULSQ
      DIRP(1,3)=UL(3)/ULSQ
      UL(1)=( DIRP(1,2)*VL(3)-DIRP(1,3)*VL(2))
      UL(2)=(-DIRP(1,1)*VL(3)+DIRP(1,3)*VL(1))
      UL(3)=( DIRP(1,1)*VL(2)-DIRP(1,2)*VL(1))
      ULSQ=SQRT(UL(1)*UL(1)+UL(2)*UL(2)+UL(3)*UL(3))
      DIRP(3,1)=UL(1)/ULSQ
      DIRP(3,2)=UL(2)/ULSQ
      DIRP(3,3)=UL(3)/ULSQ
      DIRP(2,1)= DIRP(3,2)*DIRP(1,3)-DIRP(3,3)*DIRP(1,2)
      DIRP(2,2)=-DIRP(3,1)*DIRP(1,3)+DIRP(3,3)*DIRP(1,1)
      DIRP(2,3)= DIRP(3,1)*DIRP(1,2)-DIRP(3,2)*DIRP(1,1)
C-----------------------------------------------------------------------
C     Evaluation of vertex coordinates relative to local coordinate
C     system.
C-----------------------------------------------------------------------
      DO 30 I=1,4
        XVP(I)=DIRP(1,1)*(V(1,I)-XORP(1))+
     &          DIRP(1,2)*(V(2,I)-XORP(2)) + DIRP(1,3)*(V(3,I)-XORP(3))
        YVP(I)=DIRP(2,1)*(V(1,I)-XORP(1))+
     &          DIRP(2,2)*(V(2,I)-XORP(2)) + DIRP(2,3)*(V(3,I)-XORP(3))
30    CONTINUE
C-----------------------------------------------------------------------
C     Calculation of local coordinates of Gaussian integration
C     nodes and corresponding weights. Evaluation of quadrilateral
C     centroid. Transfer of local system origin at centroid
C-----------------------------------------------------------------------
      A1 = (XVP(1)+XVP(2))*HALF
      A2 = (XVP(1)+XVP(4))*HALF
      A3 = (XVP(1)+XVP(3))*HALF
      B1 = (YVP(1)+YVP(2))*HALF
      B2 = (YVP(1)+YVP(4))*HALF
      B3 = (YVP(1)+YVP(3))*HALF
      D0 = A1*B2-A2*B1
      D1 = A1*B3-A3*B1
      D2 = A3*B2-A2*B3
      D0I = THIRD / D0
      UL(1) = (D1*A1+D2*A2)*D0I
      UL(2) = (D1*B1+D2*B2)*D0I
      DO 40 J = 1,4
        XVP(J) = XVP(J) - UL(1)
        YVP(J) = YVP(J) - UL(2)
40    CONTINUE
C-----------------------------------------------------------------------
C    Evaluation of coordinates of centroid and Gaussian nodes relative
C    to reference system
C-----------------------------------------------------------------------
      DO 60 J=1,3
        DSUM = DIRP(1,J) * UL(1) + DIRP(2,J) * UL(2)
        XCTP(J) = DSUM + XORP(J)
60    CONTINUE
C-----------------------------------------------------------------------
C     Evaluation of moments of quadrilateral surface relative to
C     local system, array S(15).  First initialize array
C     Note that S(2)=S(6)=0 due to transfer above
C-----------------------------------------------------------------------
      DO 70 N=1,15
        S(N)=ZERO
70    CONTINUE
C-----------------------------------------------------------------------
C     Evaluation of higher powers of vertex coordinates.
C-----------------------------------------------------------------------
      DO 80 N=1,4
        XP(1,N)=XVP(N)
        YP(1,N)=YVP(N)
        DO 75 M=2,6
          XP(M,N)=XP(1,N)*XP(M-1,N)
          YP(M,N)=YP(1,N)*YP(M-1,N)
75      CONTINUE
80    CONTINUE
      DO 120 N=1,4
        NXT=N1(N)
        DX(N)=XVP(N1(N))-XVP(N)
        DY(N)=YVP(N1(N))-YVP(N)
C-----------------------------------------------------------------------
C   Evaluate length of each side and accumulate area of panel
C-----------------------------------------------------------------------
        S(1)=S(1)+HALF*DX(N)*(YVP(N)+YVP(N1(N)))
        DDX=DX(N)*DX(N)
        DDY=DY(N)*DY(N)
        DDS=DDX+DDY
        IF (DDS.LT.TOLT) THEN
          SIDE(N)=ZERO
          GO TO 120
        ENDIF
        SIDE(N)=SQRT(DDS)
        IF (DDY.GT.DDX) GOTO 100
        DYDX=DY(N)/DX(N)
        DO 90 I=2,4
          I1=I+1
          I2=I+2
          SI=(XP(I1,NXT)*YP(1,NXT)-XP(I1,N)*YP(1,N))/(I1)
     &          +DYDX*(XP(I2,N)-XP(I2,NXT))/(I1*I2)
          K=K1(I1)
          S(K)=S(K)+SI
          DO 85 J=1,I
            J1=J+1
            IJ1=I-J+1
            SI=(XP(IJ1,NXT)*YP(J1,NXT)-XP(IJ1,N)*YP(J1,N))/(IJ1*J1)-
     &                                          DYDX*J*SI/IJ1
            K=J+K1(IJ1)
            S(K)=S(K)+SI
85        CONTINUE
90      CONTINUE
        GO TO 120
100     DXDY=DX(N)/DY(N)
        DO 110 I=2,4
          I1=I+1
          I2=I+2
          SI=DXDY/(I1*I2)*(YP(I2,NXT)-YP(I2,N))
          S(I1)=S(I1)+SI
          DO 105 J=1,I
            IJ=I-J
            IJ1=IJ+1
            IJ2=IJ+2
            SI=DXDY*((XP(J,NXT)*YP(IJ2,NXT)-XP(J,N)*YP(IJ2,N))/
     &                                  (IJ1*IJ2)-J*SI/IJ1)
            K=IJ+K1(J+1)
            S(K)=S(K)+SI
105       CONTINUE
110     CONTINUE
120   CONTINUE
      DO 130 N=3,15
        S(N)=S(N)*CS(N)
130   CONTINUE
      RETURN
      END










