      SUBROUTINE centroid(NPANEL,V,NUMSHAPE,XVP,YVP,XCT,DIR,AREA,
     $     SIDE,VEL)
C-----------------------------------------------------------------------
C
C     COPYRIGHT (C) 1988-1993 MASSACHUSETTS INSTITUTE OF TECHNOLOGY
C
C-----------------------------------------------------------------------
C
C     DESCRIPTION:  This routine is a reduced form of PANEL which is a 
C     utility for computing all of the panel specific
C     quantities needed for integration of quantities over the
C     panels themselves and over the surfaces they represent.
C     The arg-list of panel has been shortened to include only what
C     we need here, and lines of code have been commented out but left 
C     here in case we would like to compute more panel attributes 
C     external to FastLap.
C     A group of four vertices (not necessarily lying in a plane) is
C     used to identify a panel.  If two veritces coincide, the
C     panel is a triangle.  Some panel quantities are returned in
C     local panel coordinates, and some in the global system of
C     XVER.  The panel coordinate system has x1 and x2 lying in a
C     plane which is defined as the plane which passes through the
C     four midpoints of the four straight lines which connect the
C     four vertices.  x3 points away from your nose if the ordering of
C     the vertices looks counterclockwise to you.
C
C-----------------------------------------------------------------------
C
C     Parameters:
C
C     name         description
C     ------------------------------------------------------------------
C     NPANEL       Number of panels.
C
C     V(J,K,I)     Cartesion coordinates of vertices of quadralateral
C                  (or triangular) panels.  [If two vertices coincide,
C                  the panel is a triangle.]
C                       I=1,2,..NPANEL is the panel index
C                       K=1,2,3,4 is the vertex number.
C                       J=1,2,3 corresponds to cartesian coordinates.
C                  These coordinates imply the global coordinate system.
C                  If the vertex ordering appears counterclockwise to
C                  you as you look at the panel, then the panel normal
C                  will point away from you.  The coordinates of a panel
C                  need not lie in a plane.  The veritces in the local
C                  panel coordinates will lie in a plane (see below).
C
C     Variables.
C
C     name           description
C     ------------------------------------------------------------------
C     XVP(I,N)       x-coordinate of I-th vertex of the N-th
C                    panel in panel coordinates.
C
C     YVP(I,N)       y-coordinate of I-th vertex of the N-th
C                    panel in panel coordinates.
C                    XVP and YVP lie in the plane formed by the 4
C                    midpoints of the 4 sides of the panel.
C
C     XCT(I,N)       Cartesian coordinates of the centroid of the N-th
C                    panel in the global coordinate system.
C
C
C     DIR(I,J,N)     Direction cosines between I-th axis of the panel
C                    system and J-th axis of global system for the
C                    N-th panel.  The z-axis of the panel system is
C                    coincident with the panel normal, and therefore
C                    points away from you if the panel vertices appear
C                    counterclockwise to you.
C
C     AREA(N)        Area of N-th panel.
C
C     SIDE(I,N)      Lengths of the 4 sides of the N-th
C                    panel numbered the in the same order as the 4
C                    vertices and defined so that SID(1,N) connects
C                    vertex 1 to vertex 2, etc.
C
C     VEL(I,N)       Normal velocities at the N'th panel centroid
C                    due to rigid-body modes of motion of the surface
C                    on which the panel lies.
C                      I=1,2,...6 corresponds to the 3 translations and
C                      3 rotations.
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
C Parameters.
C
      INTEGER NPANEL
      INTEGER NUMSHAPE(NPANEL)
      DOUBLE PRECISION V(3,4,NPANEL)
C
C Variables.
C
      DOUBLE PRECISION XVP(4,NPANEL), YVP(4,NPANEL), XCT(3,NPANEL),
     $     DIR(3,3,NPANEL), AREA(NPANEL), SIDE(4,NPANEL),VEL(6,NPANEL)
C
C Local variables.
C
      DOUBLE PRECISION XORP(3), UL(3), VL(3), DX(4), DY(4)
      DOUBLE PRECISION ULSQ, A1, A2, A3, B1, B2, B3, D0, D1, D2, 
     $     D0I, DDS
      INTEGER I, J, K, N, NEXTN, FLAG
      DOUBLE PRECISION ONE, HALF, QUART, ZERO, TOLS, TOLT, TOLU, THIRD
      INTEGER NINDX(4)
C
C TOLS is the minimum panel area, TOLT is the minimum square of a panel
C side, TOLU is used for underflow protection. 
C
      DATA ONE /1.0D+0/, HALF /0.5D+0/, QUART /0.25D+0/, ZERO/0.0D+0/,
     $     TOLS /1.D-10/, TOLT /1.D-12 /, TOLU /1.D-15 /,
     $     THIRD /0.3333333333333D+0/, NINDX / 2, 3, 4, 1 /
C
C Initialize the bad panel flag.
C
      FLAG = 0
C
C Loop through panels.
C
      DO 10 I=1,NPANEL
C
C This entire routine was written to automatically handle triangles by
C exploiting a repeated vertex.  On the other hand, FastLap uses shape
C indicators.  So here we duplicate a vertex if the shape indicator is
C a triangle.
C
         IF (NUMSHAPE(I) .EQ. 3) THEN
            DO 109 J=1,3
               V(J,4,I) = V(J,3,I)
 109        CONTINUE
         END IF
         DO 110 J=1,3
            XORP(J) = ((V(J,1,I)+V(J,2,I))
     $           + (V(J,3,I)+V(J,4,I))) * QUART
            UL(J) = (V(J,1,I)+V(J,2,I)) - (V(J,3,I)+V(J,4,I))
            VL(J) = (V(J,1,I)+V(J,4,I)) - (V(J,2,I)+V(J,3,I))
 110     CONTINUE
C
C Evaluation of direction cosines DIR(3,3) of local coordinate
C system
C
         ULSQ = ONE / SQRT(UL(1)*UL(1)+UL(2)*UL(2)+UL(3)*UL(3))
         DIR(1,1,I) = UL(1)*ULSQ
         DIR(1,2,I) = UL(2)*ULSQ
         DIR(1,3,I) = UL(3)*ULSQ
         UL(1) = DIR(1,2,I)*VL(3) - DIR(1,3,I)*VL(2)
         UL(2) = DIR(1,3,I)*VL(1) - DIR(1,1,I)*VL(3)
         UL(3) = DIR(1,1,I)*VL(2) - DIR(1,2,I)*VL(1)
         ULSQ = ONE / SQRT(UL(1)*UL(1)+UL(2)*UL(2)+UL(3)*UL(3))
         DIR(3,1,I) = UL(1)*ULSQ
         DIR(3,2,I) = UL(2)*ULSQ
         DIR(3,3,I) = UL(3)*ULSQ
         DIR(2,1,I) = DIR(3,2,I)*DIR(1,3,I) - DIR(3,3,I)*DIR(1,2,I)
         DIR(2,2,I) = DIR(3,3,I)*DIR(1,1,I) - DIR(3,1,I)*DIR(1,3,I)
         DIR(2,3,I) = DIR(3,1,I)*DIR(1,2,I) - DIR(3,2,I)*DIR(1,1,I)
C
C Evaluation of vertex coordinates relative to local coordinate system.
C
         DO 210 K=1,4
            VL(1) = V(1,K,I) - XORP(1)
            VL(2) = V(2,K,I) - XORP(2)
            VL(3) = V(3,K,I) - XORP(3)
            XVP(K,I) = DIR(1,1,I)*VL(1)
     $           + DIR(1,2,I)*VL(2) + DIR(1,3,I)*VL(3)
            YVP(K,I) = DIR(2,1,I)*VL(1)
     $           + DIR(2,2,I)*VL(2) + DIR(2,3,I)*VL(3)
 210     CONTINUE
C
C Transfer of local system origin at centroid.
C
         A1 = (XVP(1,I)+XVP(2,I)) * HALF
         A2 = (XVP(1,I)+XVP(4,I)) * HALF
         A3 = (XVP(1,I)+XVP(3,I)) * HALF
         B1 = (YVP(1,I)+YVP(2,I)) * HALF
         B2 = (YVP(1,I)+YVP(4,I)) * HALF
         B3 = (YVP(1,I)+YVP(3,I)) * HALF
         D0 = A1*B2-A2*B1
         D1 = A1*B3-A3*B1
         D2 = A3*B2-A2*B3
         D0I = THIRD / D0
         UL(1)=(D1*A1+D2*A2) * D0I
         UL(2)=(D1*B1+D2*B2) * D0I
C
C   Evaluation of coordinates of centroid in body coordinates.
C
         XCT(1,I) = DIR(1,1,I)*UL(1) + DIR(2,1,I)*UL(2) + XORP(1)
         XCT(2,I) = DIR(1,2,I)*UL(1) + DIR(2,2,I)*UL(2) + XORP(2)
         XCT(3,I) = DIR(1,3,I)*UL(1) + DIR(2,3,I)*UL(2) + XORP(3)
         DO 310 K=1,4
            XVP(K,I) = XVP(K,I) - UL(1)
            YVP(K,I) = YVP(K,I) - UL(2)
 310     CONTINUE
C
C Evaluation area of quadrilateral surface. 
C
         AREA(I) = ZERO
         DO 510 N = 1,4
            NEXTN = NINDX(N)
            DX(N) = XVP(NEXTN,I) - XVP(N,I)
            DY(N) = YVP(NEXTN,I) - YVP(N,I)
C
C Evaluate length of each side and accumulate area of panel.
C
            AREA(I) = AREA(I) + HALF*DX(N)*(YVP(N,I) + YVP(NEXTN,I))
            DDS = (DX(N)*DX(N)) + (DY(N)*DY(N))
            IF (DDS .GE. TOLT) THEN
               SIDE(N,I) = SQRT(DDS)
            ELSE 
               SIDE(N,I) = ZERO
            END IF
 510     CONTINUE
C
C We have enough info now to test for bad panels.
C
         IF (AREA(I) .LT. TOLS) THEN
            PRINT *,' FE-PANEL: panel',I,' area below tolerance.'
            STOP
         END IF
         DO 710 K=1,4
            IF (DX(NINDX(K))*DY(K) .LT. DX(K)*DY(NINDX(K))) THEN
               FLAG = FLAG+1
            END IF
 710     CONTINUE
         IF (FLAG .EQ. 1) THEN
            PRINT *,' FE-PANEL: panel',I,' side is convex '
            STOP
         ELSE IF (FLAG .EQ. 2) THEN
            PRINT *,' FE-PANEL: panel',I,' sides intersect '
            STOP
         END IF
C
C Evaluation of normal vector components at centroid of
C I-th panel relative to global system.
C
         VEL(1,I) = DIR(3,1,I)
         VEL(2,I) = DIR(3,2,I)
         VEL(3,I) = DIR(3,3,I)
         VEL(4,I) = XCT(2,I)*DIR(3,3,I) - XCT(3,I)*DIR(3,2,I)
         VEL(5,I) = XCT(3,I)*DIR(3,1,I) - XCT(1,I)*DIR(3,3,I)
         VEL(6,I) = XCT(1,I)*DIR(3,2,I) - XCT(2,I)*DIR(3,1,I)
 10   CONTINUE
      RETURN
      END


