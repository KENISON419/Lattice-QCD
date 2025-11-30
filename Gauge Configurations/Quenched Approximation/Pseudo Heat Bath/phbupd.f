C***********************************************************
C   PSEUDO-HEAT BATH UPDATE ALGORITHM
C      .Originally written by S. Hioki ( 1994 for Paragon )
C      .Modifyed by Matsufuru,H. ( 26 May 1995 )
C
C       .Modified               3 Sep 1998  Matsufuru,H.
C       .For vector proc.      14 Jul 1999  Matsufuru,H.
C***********************************************************
      SUBROUTINE PHBUPD( Ioe,MU,IT )
C
      INCLUDE 'lattsize.h'
      INCLUDE 'updcomm.h'
C------------Variables--------------------------------------
      COMMON /TEMP2/  Z(NS2,4),ID( NS2 )
C------------Main-------------------------------------------
C
      CALL  TMP1( Ioe,MU,IT )
      CALL SU2UP
      CALL SPRD1( Ioe,MU,IT )
C
      CALL  TMP2( Ioe,MU,IT )
      CALL SU2UP 
      CALL SPRD2( Ioe,MU,IT )
C
      CALL  TMP3( Ioe,MU,IT )
      CALL SU2UP 
      CALL SPRD3( Ioe,MU,IT )
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE SU2UP 
C                   ----- Kennedy-Pendleton -----
      INCLUDE 'lattsize.h'
      INCLUDE 'updcomm.h'
C------------Variables--------------------------------------
      PARAMETER(TWOPI=2.D0*PI )
      COMMON /TEMP2/  Z(NS2,4),ID( NS2 )
      DIMENSION A0(NS2),DK(NS2)
      REAL*4 RN(NS2*4)
C------------Main-------------------------------------------
      BETA15=3.D0/BETA
      NRN=NS2
C
      DO 1020 is2=1,NS2
        DK(is2) = 1.D0/DSQRT(  Z(is2,1)*Z(is2,1)
     &                       + Z(is2,2)*Z(is2,2)
     &                       + Z(is2,3)*Z(is2,3)
     &                       + Z(is2,4)*Z(is2,4) )
        ID(is2) = 0
 1020 CONTINUE
C
C --- Kennedy-Pendleton ---
1099  CONTINUE
      II=NS2
      DO 1040 i=1,NS2
         II=II-ID(i)
 1040 CONTINUE
C
      IF( II .EQ. 0 ) GOTO 1031
C
      CALL RNDS(RN,4*NRN,4*II)
      J=-3
      DO 1041 i=1,NS2
        IF( ID(i) .EQ. 1 ) GOTO 1041
        J=J+4
       TT=COS(TWOPI*RN(J+2))
       D2 = -(LOG(RN(J))+LOG(RN(J+1))*TT*TT)*DK(I)*BETA15
       IF( RN(J+3)*RN(J+3) .GT. (1.0-0.5*D2) ) GOTO 1041
         ID(i) = 1
         A0(i) = 1.0 - D2
 1041 CONTINUE
      GOTO 1099
C
 1031 CONTINUE
C
      CALL RNDS(RN,4*NRN,2*NRN)
C
      DO 1060 is2=1,NS2
        J=is2+is2-1
        RAD = 1. - A0(is2)*A0(is2)
        A3 = DSQRT(RAD)*( 2.*RN(J) -1. )
        RAD2 = DSQRT( DABS(RAD-A3*A3) )
        THETA = TWOPI*RN(J+1)
        A1 = RAD2*DCOS(THETA)
        A2 = RAD2*DSIN(THETA)
         V1 = Z(is2,1)
         V2 = Z(is2,2)
         V3 = Z(is2,3)
         V4 = Z(is2,4)
        Z(is2,1) = DK(is2)*( A0(is2)*V1+A1*V2+A2*V3+A3*V4)
        Z(is2,2) = DK(is2)*(-A0(is2)*V2+A1*V1+A2*V4-A3*V3)
        Z(is2,3) = DK(is2)*(-A0(is2)*V3-A1*V4+A2*V1+A3*V2)
        Z(is2,4) = DK(is2)*(-A0(is2)*V4+A1*V3-A2*V2+A3*V1)
C
 1060 CONTINUE
C 
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE SPRD1( Ioe,MU,IT )
C                      ----- Y[SU(3)]=C[SU(2)]*X[SU(3)]-----
      INCLUDE 'lattsize.h'
      INCLUDE 'updcomm.h'
C------------Variables--------------------------------------
      COMMON /TEMP2/  Z(NS2,4),ID( NS2 )
C------------Main-------------------------------------------
      DO 100 is2 = 1, NS2
C
      Y1  =  Z( is2, 1 ) * U( is2, 1,Ioe,MU,IT ) 
     &      -Z( is2, 4 ) * U( is2, 2,Ioe,MU,IT )
     &      +Z( is2, 3 ) * U( is2, 7,Ioe,MU,IT )
     &      -Z( is2, 2 ) * U( is2, 8,Ioe,MU,IT )
      Y2  =  Z( is2, 1 ) * U( is2, 2,Ioe,MU,IT ) 
     &      +Z( is2, 4 ) * U( is2, 1,Ioe,MU,IT )
     &      +Z( is2, 3 ) * U( is2, 8,Ioe,MU,IT )
     &      +Z( is2, 2 ) * U( is2, 7,Ioe,MU,IT )
      Y3  =  Z( is2, 1 ) * U( is2, 3,Ioe,MU,IT ) 
     &      -Z( is2, 4 ) * U( is2, 4,Ioe,MU,IT )
     &      +Z( is2, 3 ) * U( is2, 9,Ioe,MU,IT )
     &      -Z( is2, 2 ) * U( is2,10,Ioe,MU,IT )
      Y4  =  Z( is2, 1 ) * U( is2, 4,Ioe,MU,IT ) 
     &      +Z( is2, 4 ) * U( is2, 3,Ioe,MU,IT )
     &      +Z( is2, 3 ) * U( is2,10,Ioe,MU,IT )
     &      +Z( is2, 2 ) * U( is2, 9,Ioe,MU,IT )
      Y5  =  Z( is2, 1 ) * U( is2, 5,Ioe,MU,IT ) 
     &      -Z( is2, 4 ) * U( is2, 6,Ioe,MU,IT )
     &      +Z( is2, 3 ) * U( is2,11,Ioe,MU,IT )
     &      -Z( is2, 2 ) * U( is2,12,Ioe,MU,IT )
      Y6  =  Z( is2, 1 ) * U( is2, 6,Ioe,MU,IT ) 
     &      +Z( is2, 4 ) * U( is2, 5,Ioe,MU,IT )
     &      +Z( is2, 3 ) * U( is2,12,Ioe,MU,IT )
     &      +Z( is2, 2 ) * U( is2,11,Ioe,MU,IT )
      Y7  = -Z( is2, 3 ) * U( is2, 1,Ioe,MU,IT ) 
     &      -Z( is2, 2 ) * U( is2, 2,Ioe,MU,IT )
     &      +Z( is2, 1 ) * U( is2, 7,Ioe,MU,IT )
     &      +Z( is2, 4 ) * U( is2, 8,Ioe,MU,IT )
      Y8  = -Z( is2, 3 ) * U( is2, 2,Ioe,MU,IT ) 
     &      +Z( is2, 2 ) * U( is2, 1,Ioe,MU,IT )
     &      +Z( is2, 1 ) * U( is2, 8,Ioe,MU,IT )
     &      -Z( is2, 4 ) * U( is2, 7,Ioe,MU,IT )
      Y9  = -Z( is2, 3 ) * U( is2, 3,Ioe,MU,IT ) 
     &      -Z( is2, 2 ) * U( is2, 4,Ioe,MU,IT )
     &      +Z( is2, 1 ) * U( is2, 9,Ioe,MU,IT )
     &      +Z( is2, 4 ) * U( is2,10,Ioe,MU,IT )
      Y10 = -Z( is2, 3 ) * U( is2, 4,Ioe,MU,IT ) 
     &      +Z( is2, 2 ) * U( is2, 3,Ioe,MU,IT )
     &      +Z( is2, 1 ) * U( is2,10,Ioe,MU,IT )
     &      -Z( is2, 4 ) * U( is2, 9,Ioe,MU,IT )
      Y11 = -Z( is2, 3 ) * U( is2, 5,Ioe,MU,IT ) 
     &      -Z( is2, 2 ) * U( is2, 6,Ioe,MU,IT )
     &      +Z( is2, 1 ) * U( is2,11,Ioe,MU,IT )
     &      +Z( is2, 4 ) * U( is2,12,Ioe,MU,IT )
      Y12 = -Z( is2, 3 ) * U( is2, 6,Ioe,MU,IT ) 
     &      +Z( is2, 2 ) * U( is2, 5,Ioe,MU,IT )
     &      +Z( is2, 1 ) * U( is2,12,Ioe,MU,IT )
     &      -Z( is2, 4 ) * U( is2,11,Ioe,MU,IT )
C
      U( is2, 1,Ioe,MU,IT ) = Y1
      U( is2, 2,Ioe,MU,IT ) = Y2
      U( is2, 3,Ioe,MU,IT ) = Y3
      U( is2, 4,Ioe,MU,IT ) = Y4
      U( is2, 5,Ioe,MU,IT ) = Y5
      U( is2, 6,Ioe,MU,IT ) = Y6
      U( is2, 7,Ioe,MU,IT ) = Y7
      U( is2, 8,Ioe,MU,IT ) = Y8
      U( is2, 9,Ioe,MU,IT ) = Y9
      U( is2,10,Ioe,MU,IT ) = Y10
      U( is2,11,Ioe,MU,IT ) = Y11
      U( is2,12,Ioe,MU,IT ) = Y12
C
  100 CONTINUE
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE SPRD2( Ioe,MU,IT )
C
      INCLUDE 'lattsize.h'
      INCLUDE 'updcomm.h'
C------------Variables--------------------------------------
      COMMON /TEMP2/  Z(NS2,4),ID( NS2 )
C------------Main-------------------------------------------
      DO 100 is2 = 1, NS2
C
      Y7  = Z( is2,1 ) * U( is2, 7,Ioe,MU,IT ) 
     &     -Z( is2,4 ) * U( is2, 8,Ioe,MU,IT )
     &     -Z( is2,3 ) * U( is2,13,Ioe,MU,IT )
     &     -Z( is2,2 ) * U( is2,14,Ioe,MU,IT )
      Y8  = Z( is2,1 ) * U( is2, 8,Ioe,MU,IT ) 
     &     +Z( is2,4 ) * U( is2, 7,Ioe,MU,IT )
     &     -Z( is2,3 ) * U( is2,14,Ioe,MU,IT )
     &     +Z( is2,2 ) * U( is2,13,Ioe,MU,IT )
      Y9  = Z( is2,1 ) * U( is2, 9,Ioe,MU,IT ) 
     &     -Z( is2,4 ) * U( is2,10,Ioe,MU,IT )
     &     -Z( is2,3 ) * U( is2,15,Ioe,MU,IT )
     &     -Z( is2,2 ) * U( is2,16,Ioe,MU,IT )
      Y10 = Z( is2,1 ) * U( is2,10,Ioe,MU,IT ) 
     &     +Z( is2,4 ) * U( is2, 9,Ioe,MU,IT )
     &     -Z( is2,3 ) * U( is2,16,Ioe,MU,IT )
     &     +Z( is2,2 ) * U( is2,15,Ioe,MU,IT )
      Y11 = Z( is2,1 ) * U( is2,11,Ioe,MU,IT ) 
     &     -Z( is2,4 ) * U( is2,12,Ioe,MU,IT )
     &     -Z( is2,3 ) * U( is2,17,Ioe,MU,IT )
     &     -Z( is2,2 ) * U( is2,18,Ioe,MU,IT )
      Y12 = Z( is2,1 ) * U( is2,12,Ioe,MU,IT ) 
     &     +Z( is2,4 ) * U( is2,11,Ioe,MU,IT )
     &     -Z( is2,3 ) * U( is2,18,Ioe,MU,IT )
     &     +Z( is2,2 ) * U( is2,17,Ioe,MU,IT )
      Y13 = Z( is2,3 ) * U( is2, 7,Ioe,MU,IT ) 
     &     -Z( is2,2 ) * U( is2, 8,Ioe,MU,IT )
     &     +Z( is2,1 ) * U( is2,13,Ioe,MU,IT )
     &     +Z( is2,4 ) * U( is2,14,Ioe,MU,IT )
      Y14 = Z( is2,3 ) * U( is2, 8,Ioe,MU,IT ) 
     &     +Z( is2,2 ) * U( is2, 7,Ioe,MU,IT )
     &     +Z( is2,1 ) * U( is2,14,Ioe,MU,IT )
     &     -Z( is2,4 ) * U( is2,13,Ioe,MU,IT )
      Y15 = Z( is2,3 ) * U( is2, 9,Ioe,MU,IT ) 
     &     -Z( is2,2 ) * U( is2,10,Ioe,MU,IT )
     &     +Z( is2,1 ) * U( is2,15,Ioe,MU,IT )
     &     +Z( is2,4 ) * U( is2,16,Ioe,MU,IT )
      Y16 = Z( is2,3 ) * U( is2,10,Ioe,MU,IT ) 
     &     +Z( is2,2 ) * U( is2, 9,Ioe,MU,IT )
     &     +Z( is2,1 ) * U( is2,16,Ioe,MU,IT )
     &     -Z( is2,4 ) * U( is2,15,Ioe,MU,IT )
      Y17 = Z( is2,3 ) * U( is2,11,Ioe,MU,IT ) 
     &     -Z( is2,2 ) * U( is2,12,Ioe,MU,IT )
     &     +Z( is2,1 ) * U( is2,17,Ioe,MU,IT )
     &     +Z( is2,4 ) * U( is2,18,Ioe,MU,IT )
      Y18 = Z( is2,3 ) * U( is2,12,Ioe,MU,IT ) 
     &     +Z( is2,2 ) * U( is2,11,Ioe,MU,IT )
     &     +Z( is2,1 ) * U( is2,18,Ioe,MU,IT )
     &     -Z( is2,4 ) * U( is2,17,Ioe,MU,IT )
C
      U( is2, 7,Ioe,MU,IT ) = Y7
      U( is2, 8,Ioe,MU,IT ) = Y8
      U( is2, 9,Ioe,MU,IT ) = Y9
      U( is2,10,Ioe,MU,IT ) = Y10
      U( is2,11,Ioe,MU,IT ) = Y11
      U( is2,12,Ioe,MU,IT ) = Y12
      U( is2,13,Ioe,MU,IT ) = Y13
      U( is2,14,Ioe,MU,IT ) = Y14
      U( is2,15,Ioe,MU,IT ) = Y15
      U( is2,16,Ioe,MU,IT ) = Y16
      U( is2,17,Ioe,MU,IT ) = Y17
      U( is2,18,Ioe,MU,IT ) = Y18
C
  100 CONTINUE
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE SPRD3( Ioe,MU,IT )
C
      INCLUDE 'lattsize.h'
      INCLUDE 'updcomm.h'
C------------Variables--------------------------------------
      COMMON /TEMP2/  Z(NS2,4),ID( NS2 )
C------------Main-------------------------------------------
      DO 100 is2 = 1, NS2
C
      Y1  =  Z( is2,1 ) * U( is2, 1,Ioe,MU,IT ) 
     &      -Z( is2,4 ) * U( is2, 2,Ioe,MU,IT )
     &      +Z( is2,3 ) * U( is2,13,Ioe,MU,IT )
     &      -Z( is2,2 ) * U( is2,14,Ioe,MU,IT )
      Y2  =  Z( is2,1 ) * U( is2, 2,Ioe,MU,IT ) 
     &      +Z( is2,4 ) * U( is2, 1,Ioe,MU,IT )
     &      +Z( is2,3 ) * U( is2,14,Ioe,MU,IT )
     &      +Z( is2,2 ) * U( is2,13,Ioe,MU,IT )
      Y3  =  Z( is2,1 ) * U( is2, 3,Ioe,MU,IT ) 
     &      -Z( is2,4 ) * U( is2, 4,Ioe,MU,IT )
     &      +Z( is2,3 ) * U( is2,15,Ioe,MU,IT )
     &      -Z( is2,2 ) * U( is2,16,Ioe,MU,IT )
      Y4  =  Z( is2,1 ) * U( is2, 4,Ioe,MU,IT ) 
     &      +Z( is2,4 ) * U( is2, 3,Ioe,MU,IT )
     &      +Z( is2,3 ) * U( is2,16,Ioe,MU,IT )
     &      +Z( is2,2 ) * U( is2,15,Ioe,MU,IT )
      Y5  =  Z( is2,1 ) * U( is2, 5,Ioe,MU,IT ) 
     &      -Z( is2,4 ) * U( is2, 6,Ioe,MU,IT )
     &      +Z( is2,3 ) * U( is2,17,Ioe,MU,IT )
     &      -Z( is2,2 ) * U( is2,18,Ioe,MU,IT )
      Y6  =  Z( is2,1 ) * U( is2, 6,Ioe,MU,IT ) 
     &      +Z( is2,4 ) * U( is2, 5,Ioe,MU,IT )
     &      +Z( is2,3 ) * U( is2,18,Ioe,MU,IT )
     &      +Z( is2,2 ) * U( is2,17,Ioe,MU,IT )
      Y13 = -Z( is2,3 ) * U( is2, 1,Ioe,MU,IT ) 
     &      -Z( is2,2 ) * U( is2, 2,Ioe,MU,IT )
     &      +Z( is2,1 ) * U( is2,13,Ioe,MU,IT )
     &      +Z( is2,4 ) * U( is2,14,Ioe,MU,IT )
      Y14 = -Z( is2,3 ) * U( is2, 2,Ioe,MU,IT ) 
     &      +Z( is2,2 ) * U( is2, 1,Ioe,MU,IT )
     &      +Z( is2,1 ) * U( is2,14,Ioe,MU,IT )
     &      -Z( is2,4 ) * U( is2,13,Ioe,MU,IT )
      Y15 = -Z( is2,3 ) * U( is2, 3,Ioe,MU,IT ) 
     &      -Z( is2,2 ) * U( is2, 4,Ioe,MU,IT )
     &      +Z( is2,1 ) * U( is2,15,Ioe,MU,IT )
     &      +Z( is2,4 ) * U( is2,16,Ioe,MU,IT )
      Y16 = -Z( is2,3 ) * U( is2, 4,Ioe,MU,IT ) 
     &      +Z( is2,2 ) * U( is2, 3,Ioe,MU,IT )
     &      +Z( is2,1 ) * U( is2,16,Ioe,MU,IT )
     &      -Z( is2,4 ) * U( is2,15,Ioe,MU,IT )
      Y17 = -Z( is2,3 ) * U( is2, 5,Ioe,MU,IT ) 
     &      -Z( is2,2 ) * U( is2, 6,Ioe,MU,IT )
     &      +Z( is2,1 ) * U( is2,17,Ioe,MU,IT )
     &      +Z( is2,4 ) * U( is2,18,Ioe,MU,IT )
      Y18 = -Z( is2,3 ) * U( is2, 6,Ioe,MU,IT ) 
     &      +Z( is2,2 ) * U( is2, 5,Ioe,MU,IT )
     &      +Z( is2,1 ) * U( is2,18,Ioe,MU,IT )
     &      -Z( is2,4 ) * U( is2,17,Ioe,MU,IT )
C
          U( is2, 1,Ioe,MU,IT ) = Y1
          U( is2, 2,Ioe,MU,IT ) = Y2
          U( is2, 3,Ioe,MU,IT ) = Y3
          U( is2, 4,Ioe,MU,IT ) = Y4
          U( is2, 5,Ioe,MU,IT ) = Y5
          U( is2, 6,Ioe,MU,IT ) = Y6
          U( is2,13,Ioe,MU,IT ) = Y13
          U( is2,14,Ioe,MU,IT ) = Y14
          U( is2,15,Ioe,MU,IT ) = Y15
          U( is2,16,Ioe,MU,IT ) = Y16
          U( is2,17,Ioe,MU,IT ) = Y17
          U( is2,18,Ioe,MU,IT ) = Y18
C
  100 CONTINUE
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE TMP1( Ioe,MU,IT )
C
      INCLUDE 'lattsize.h'
      INCLUDE 'updcomm.h'
C------------Variables--------------------------------------
      COMMON /TEMP2/  Z(NS2,4),ID( NS2 )
C------------Main-------------------------------------------
      DO 100 is2 = 1, NS2
C
      Z( is2,1 ) =   U( is2, 1,Ioe,MU,IT ) * C( is2, 1 )
     &             + U( is2, 2,Ioe,MU,IT ) * C( is2, 2 )
     &             + U( is2, 3,Ioe,MU,IT ) * C( is2, 3 )
     &             + U( is2, 4,Ioe,MU,IT ) * C( is2, 4 )
     &             + U( is2, 5,Ioe,MU,IT ) * C( is2, 5 )
     &             + U( is2, 6,Ioe,MU,IT ) * C( is2, 6 )
     &             + U( is2, 7,Ioe,MU,IT ) * C( is2, 7 )
     &             + U( is2, 8,Ioe,MU,IT ) * C( is2, 8 )
     &             + U( is2, 9,Ioe,MU,IT ) * C( is2, 9 )
     &             + U( is2,10,Ioe,MU,IT ) * C( is2,10 )
     &             + U( is2,11,Ioe,MU,IT ) * C( is2,11 )
     &             + U( is2,12,Ioe,MU,IT ) * C( is2,12 )
C
      Z( is2,2 ) = - U( is2, 1,Ioe,MU,IT ) * C( is2, 8 )
     &             + U( is2, 2,Ioe,MU,IT ) * C( is2, 7 )
     &             - U( is2, 3,Ioe,MU,IT ) * C( is2,10 )
     &             + U( is2, 4,Ioe,MU,IT ) * C( is2, 9 )
     &             - U( is2, 5,Ioe,MU,IT ) * C( is2,12 )
     &             + U( is2, 6,Ioe,MU,IT ) * C( is2,11 )
     &             - U( is2, 7,Ioe,MU,IT ) * C( is2, 2 )
     &             + U( is2, 8,Ioe,MU,IT ) * C( is2, 1 )
     &             - U( is2, 9,Ioe,MU,IT ) * C( is2, 4 )
     &             + U( is2,10,Ioe,MU,IT ) * C( is2, 3 )
     &             - U( is2,11,Ioe,MU,IT ) * C( is2, 6 )
     &             + U( is2,12,Ioe,MU,IT ) * C( is2, 5 )
C
      Z( is2,3 ) =   U( is2, 1,Ioe,MU,IT ) * C( is2, 7 )
     &             + U( is2, 2,Ioe,MU,IT ) * C( is2, 8 )
     &             + U( is2, 3,Ioe,MU,IT ) * C( is2, 9 )
     &             + U( is2, 4,Ioe,MU,IT ) * C( is2,10 )
     &             + U( is2, 5,Ioe,MU,IT ) * C( is2,11 )
     &             + U( is2, 6,Ioe,MU,IT ) * C( is2,12 )
     &             - U( is2, 7,Ioe,MU,IT ) * C( is2, 1 )
     &             - U( is2, 8,Ioe,MU,IT ) * C( is2, 2 )
     &             - U( is2, 9,Ioe,MU,IT ) * C( is2, 3 )
     &             - U( is2,10,Ioe,MU,IT ) * C( is2, 4 )
     &             - U( is2,11,Ioe,MU,IT ) * C( is2, 5 )
     &             - U( is2,12,Ioe,MU,IT ) * C( is2, 6 )
C
      Z( is2,4 ) = - U( is2, 1,Ioe,MU,IT ) * C( is2, 2 )
     &             + U( is2, 2,Ioe,MU,IT ) * C( is2, 1 )
     &             - U( is2, 3,Ioe,MU,IT ) * C( is2, 4 )
     &             + U( is2, 4,Ioe,MU,IT ) * C( is2, 3 )
     &             - U( is2, 5,Ioe,MU,IT ) * C( is2, 6 )
     &             + U( is2, 6,Ioe,MU,IT ) * C( is2, 5 )
     &             + U( is2, 7,Ioe,MU,IT ) * C( is2, 8 )
     &             - U( is2, 8,Ioe,MU,IT ) * C( is2, 7 )
     &             + U( is2, 9,Ioe,MU,IT ) * C( is2,10 )
     &             - U( is2,10,Ioe,MU,IT ) * C( is2, 9 )
     &             + U( is2,11,Ioe,MU,IT ) * C( is2,12 )
     &             - U( is2,12,Ioe,MU,IT ) * C( is2,11 )
C
 100   CONTINUE
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE TMP2( Ioe,MU,IT )
C
      INCLUDE 'lattsize.h'
      INCLUDE 'updcomm.h'
C------------Variables--------------------------------------
      COMMON /TEMP2/  Z(NS2,4),ID( NS2 )
C------------Main-------------------------------------------
      DO 100 is2 = 1, NS2
C
      Z( is2,1 ) =   U( is2, 7,Ioe,MU,IT ) * C( is2, 7 )
     &             + U( is2, 8,Ioe,MU,IT ) * C( is2, 8 )
     &             + U( is2, 9,Ioe,MU,IT ) * C( is2, 9 )
     &             + U( is2,10,Ioe,MU,IT ) * C( is2,10 )
     &             + U( is2,11,Ioe,MU,IT ) * C( is2,11 )
     &             + U( is2,12,Ioe,MU,IT ) * C( is2,12 )
     &             + U( is2,13,Ioe,MU,IT ) * C( is2,13 )
     &             + U( is2,14,Ioe,MU,IT ) * C( is2,14 )
     &             + U( is2,15,Ioe,MU,IT ) * C( is2,15 )
     &             + U( is2,16,Ioe,MU,IT ) * C( is2,16 )
     &             + U( is2,17,Ioe,MU,IT ) * C( is2,17 )
     &             + U( is2,18,Ioe,MU,IT ) * C( is2,18 )
C
      Z( is2,2 ) = - U( is2, 7,Ioe,MU,IT ) * C( is2,14 )
     &             + U( is2, 8,Ioe,MU,IT ) * C( is2,13 )
     &             - U( is2, 9,Ioe,MU,IT ) * C( is2,16 )
     &             + U( is2,10,Ioe,MU,IT ) * C( is2,15 )
     &             - U( is2,11,Ioe,MU,IT ) * C( is2,18 )
     &             + U( is2,12,Ioe,MU,IT ) * C( is2,17 )
     &             - U( is2,13,Ioe,MU,IT ) * C( is2, 8 )
     &             + U( is2,14,Ioe,MU,IT ) * C( is2, 7 )
     &             - U( is2,15,Ioe,MU,IT ) * C( is2,10 )
     &             + U( is2,16,Ioe,MU,IT ) * C( is2, 9 )
     &             - U( is2,17,Ioe,MU,IT ) * C( is2,12 )
     &             + U( is2,18,Ioe,MU,IT ) * C( is2,11 )
C
      Z( is2,3 ) =   U( is2,13,Ioe,MU,IT ) * C( is2, 7 )
     &             + U( is2,14,Ioe,MU,IT ) * C( is2, 8 )
     &             + U( is2,15,Ioe,MU,IT ) * C( is2, 9 )
     &             + U( is2,16,Ioe,MU,IT ) * C( is2,10 )
     &             + U( is2,17,Ioe,MU,IT ) * C( is2,11 )
     &             + U( is2,18,Ioe,MU,IT ) * C( is2,12 )
     &             - U( is2, 7,Ioe,MU,IT ) * C( is2,13 )
     &             - U( is2, 8,Ioe,MU,IT ) * C( is2,14 )
     &             - U( is2, 9,Ioe,MU,IT ) * C( is2,15 )
     &             - U( is2,10,Ioe,MU,IT ) * C( is2,16 )
     &             - U( is2,11,Ioe,MU,IT ) * C( is2,17 )
     &             - U( is2,12,Ioe,MU,IT ) * C( is2,18 )
C
      Z( is2,4 ) = - U( is2, 7,Ioe,MU,IT ) * C( is2, 8 )
     &             + U( is2, 8,Ioe,MU,IT ) * C( is2, 7 )
     &             - U( is2, 9,Ioe,MU,IT ) * C( is2,10 )
     &             + U( is2,10,Ioe,MU,IT ) * C( is2, 9 )
     &             - U( is2,11,Ioe,MU,IT ) * C( is2,12 )
     &             + U( is2,12,Ioe,MU,IT ) * C( is2,11 )
     &             + U( is2,13,Ioe,MU,IT ) * C( is2,14 )
     &             - U( is2,14,Ioe,MU,IT ) * C( is2,13 )
     &             + U( is2,15,Ioe,MU,IT ) * C( is2,16 )
     &             - U( is2,16,Ioe,MU,IT ) * C( is2,15 )
     &             + U( is2,17,Ioe,MU,IT ) * C( is2,18 )
     &             - U( is2,18,Ioe,MU,IT ) * C( is2,17 )
C
 100  CONTINUE
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE TMP3( Ioe,MU,IT )
C
      INCLUDE 'lattsize.h'
      INCLUDE 'updcomm.h'
C------------Variables--------------------------------------
      COMMON /TEMP2/  Z(NS2,4),ID( NS2 )
C------------Main-------------------------------------------
      DO 100 is2 = 1, NS2
C
      Z( is2,1 ) =   U( is2, 1,Ioe,MU,IT ) * C( is2, 1 )
     &             + U( is2, 2,Ioe,MU,IT ) * C( is2, 2 )
     &             + U( is2, 3,Ioe,MU,IT ) * C( is2, 3 )
     &             + U( is2, 4,Ioe,MU,IT ) * C( is2, 4 )
     &             + U( is2, 5,Ioe,MU,IT ) * C( is2, 5 )
     &             + U( is2, 6,Ioe,MU,IT ) * C( is2, 6 )
     &             + U( is2,13,Ioe,MU,IT ) * C( is2,13 )
     &             + U( is2,14,Ioe,MU,IT ) * C( is2,14 )
     &             + U( is2,15,Ioe,MU,IT ) * C( is2,15 )
     &             + U( is2,16,Ioe,MU,IT ) * C( is2,16 )
     &             + U( is2,17,Ioe,MU,IT ) * C( is2,17 )
     &             + U( is2,18,Ioe,MU,IT ) * C( is2,18 )
C
      Z( is2,2 ) = - U( is2, 1,Ioe,MU,IT ) * C( is2,14 )
     &             + U( is2, 2,Ioe,MU,IT ) * C( is2,13 )
     &             - U( is2, 3,Ioe,MU,IT ) * C( is2,16 )
     &             + U( is2, 4,Ioe,MU,IT ) * C( is2,15 )
     &             - U( is2, 5,Ioe,MU,IT ) * C( is2,18 )
     &             + U( is2, 6,Ioe,MU,IT ) * C( is2,17 )
     &             - U( is2,13,Ioe,MU,IT ) * C( is2, 2 )
     &             + U( is2,14,Ioe,MU,IT ) * C( is2, 1 )
     &             - U( is2,15,Ioe,MU,IT ) * C( is2, 4 )
     &             + U( is2,16,Ioe,MU,IT ) * C( is2, 3 )
     &             - U( is2,17,Ioe,MU,IT ) * C( is2, 6 )
     &             + U( is2,18,Ioe,MU,IT ) * C( is2, 5 )
C
      Z( is2,3 ) =   U( is2, 1,Ioe,MU,IT ) * C( is2,13 )
     &             + U( is2, 2,Ioe,MU,IT ) * C( is2,14 )
     &             + U( is2, 3,Ioe,MU,IT ) * C( is2,15 )
     &             + U( is2, 4,Ioe,MU,IT ) * C( is2,16 )
     &             + U( is2, 5,Ioe,MU,IT ) * C( is2,17 )
     &             + U( is2, 6,Ioe,MU,IT ) * C( is2,18 )
     &             - U( is2,13,Ioe,MU,IT ) * C( is2, 1 )
     &             - U( is2,14,Ioe,MU,IT ) * C( is2, 2 )
     &             - U( is2,15,Ioe,MU,IT ) * C( is2, 3 )
     &             - U( is2,16,Ioe,MU,IT ) * C( is2, 4 )
     &             - U( is2,17,Ioe,MU,IT ) * C( is2, 5 )
     &             - U( is2,18,Ioe,MU,IT ) * C( is2, 6 )
C
      Z( is2,4 ) = - U( is2, 1,Ioe,MU,IT ) * C( is2, 2 )
     &             + U( is2, 2,Ioe,MU,IT ) * C( is2, 1 )
     &             - U( is2, 3,Ioe,MU,IT ) * C( is2, 4 )
     &             + U( is2, 4,Ioe,MU,IT ) * C( is2, 3 )
     &             - U( is2, 5,Ioe,MU,IT ) * C( is2, 6 )
     &             + U( is2, 6,Ioe,MU,IT ) * C( is2, 5 )
     &             + U( is2,13,Ioe,MU,IT ) * C( is2,14 )
     &             - U( is2,14,Ioe,MU,IT ) * C( is2,13 )
     &             + U( is2,15,Ioe,MU,IT ) * C( is2,16 )
     &             - U( is2,16,Ioe,MU,IT ) * C( is2,15 )
     &             + U( is2,17,Ioe,MU,IT ) * C( is2,18 )
     &             - U( is2,18,Ioe,MU,IT ) * C( is2,17 )
C
 100   CONTINUE
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
C***************************************************END*****
