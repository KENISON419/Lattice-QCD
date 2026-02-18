C***********************************************************
C   SU(3) matrix * spinor field 
C                  and \gamma matrices multiplication
C
C  .For vector processor version.
C  .Need to include parameter file 'lattsize.h'
C       in which the precision is specified.
C
C       UxPrdv_N(UX,UP,WP)       --  UX = UP   * WP 
C       UxPrdv_D(UX,UP,WP)       --  UX = UP^d * WP 
C
C                               07 Jun 1998  Matsufuru,H.
C        For vector processor   12 Jun 1999  M.,H.
C        \gamma matrices        14 Jun 1999  M.,H.
C***********************************************************
      SUBROUTINE UxPrdv_N(Nv,UX,UP,WP)
C                            ---  Ux = Up * wp  ---
      INCLUDE 'lattsize.h'
C-------------------Variables-------------------------------
      DIMENSION UX(Nv,Nvc,ND),UP(Nv,Ndf),WP(Nv,Nvc,ND)
C---------------------Main----------------------------------
C
      DO 100 id = 1, ND
      DO 100 iv = 1, Nv
C
        UX(iv,1,id)
     &     =  UP(iv, 1)*WP(iv,1,id) - UP(iv, 2)*WP(iv,2,id)
     &      + UP(iv, 3)*WP(iv,3,id) - UP(iv, 4)*WP(iv,4,id)
     &      + UP(iv, 5)*WP(iv,5,id) - UP(iv, 6)*WP(iv,6,id)
C
        UX(iv,2,id)
     &     =  UP(iv, 1)*WP(iv,2,id) + UP(iv, 2)*WP(iv,1,id)
     &      + UP(iv, 3)*WP(iv,4,id) + UP(iv, 4)*WP(iv,3,id)
     &      + UP(iv, 5)*WP(iv,6,id) + UP(iv, 6)*WP(iv,5,id)
C
        UX(iv,3,id)
     &     =  UP(iv, 7)*WP(iv,1,id) - UP(iv, 8)*WP(iv,2,id)
     &      + UP(iv, 9)*WP(iv,3,id) - UP(iv,10)*WP(iv,4,id)
     &      + UP(iv,11)*WP(iv,5,id) - UP(iv,12)*WP(iv,6,id)
C
        UX(iv,4,id)
     &     =  UP(iv, 7)*WP(iv,2,id) + UP(iv, 8)*WP(iv,1,id)
     &      + UP(iv, 9)*WP(iv,4,id) + UP(iv,10)*WP(iv,3,id)
     &      + UP(iv,11)*WP(iv,6,id) + UP(iv,12)*WP(iv,5,id)
C
        UX(iv,5,id)
     &     =  UP(iv,13)*WP(iv,1,id) - UP(iv,14)*WP(iv,2,id)
     &      + UP(iv,15)*WP(iv,3,id) - UP(iv,16)*WP(iv,4,id)
     &      + UP(iv,17)*WP(iv,5,id) - UP(iv,18)*WP(iv,6,id)
C
        UX(iv,6,id)
     &     =  UP(iv,13)*WP(iv,2,id) + UP(iv,14)*WP(iv,1,id)
     &      + UP(iv,15)*WP(iv,4,id) + UP(iv,16)*WP(iv,3,id)
     &      + UP(iv,17)*WP(iv,6,id) + UP(iv,18)*WP(iv,5,id)
C
 100  CONTINUE
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE UxPrdv_D(Nv,UX,UP,WP)
C                           ---  Ux = Up^d * wp  ---
      INCLUDE 'lattsize.h'
C-------------------Variables-------------------------------
      DIMENSION UX(Nv,Nvc,ND),UP(Nv,Ndf),WP(Nv,Nvc,ND)
C---------------------Main----------------------------------
C
      DO 100 id=1,ND
      DO 100 iv = 1, Nv
C
        UX(iv,1,id)
     &     =  UP(iv, 1)*WP(iv,1,id) + UP(iv, 2)*WP(iv,2,id)
     &      + UP(iv, 7)*WP(iv,3,id) + UP(iv, 8)*WP(iv,4,id)
     &      + UP(iv,13)*WP(iv,5,id) + UP(iv,14)*WP(iv,6,id)
C
        UX(iv,2,id)
     &     =  UP(iv, 1)*WP(iv,2,id) - UP(iv, 2)*WP(iv,1,id)
     &      + UP(iv, 7)*WP(iv,4,id) - UP(iv, 8)*WP(iv,3,id)
     &      + UP(iv,13)*WP(iv,6,id) - UP(iv,14)*WP(iv,5,id)
C
        UX(iv,3,id)
     &     =  UP(iv, 3)*WP(iv,1,id) + UP(iv, 4)*WP(iv,2,id)
     &      + UP(iv, 9)*WP(iv,3,id) + UP(iv,10)*WP(iv,4,id)
     &      + UP(iv,15)*WP(iv,5,id) + UP(iv,16)*WP(iv,6,id)
C
        UX(iv,4,id)
     &     =  UP(iv, 3)*WP(iv,2,id) - UP(iv, 4)*WP(iv,1,id)
     &      + UP(iv, 9)*WP(iv,4,id) - UP(iv,10)*WP(iv,3,id)
     &      + UP(iv,15)*WP(iv,6,id) - UP(iv,16)*WP(iv,5,id)
C
        UX(iv,5,id)
     &     =  UP(iv, 5)*WP(iv,1,id) + UP(iv, 6)*WP(iv,2,id)
     &      + UP(iv,11)*WP(iv,3,id) + UP(iv,12)*WP(iv,4,id)
     &      + UP(iv,17)*WP(iv,5,id) + UP(iv,18)*WP(iv,6,id)
C
        UX(iv,6,id)
     &     =  UP(iv, 5)*WP(iv,2,id) - UP(iv, 6)*WP(iv,1,id)
     &      + UP(iv,11)*WP(iv,4,id) - UP(iv,12)*WP(iv,3,id)
     &      + UP(iv,17)*WP(iv,6,id) - UP(iv,18)*WP(iv,5,id)
C
 100  CONTINUE
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE G1Prd(Nv,Ipm,V2,V1,BF)
C          ---  V2 = BF * ( 1 +/- \gamma_1 ) * V1 ---
C           Ipm = 1: +,  -1: -
C           BF: Boundary condition
C
      INCLUDE 'lattsize.h'
C------------Variables--------------------------------------
      DIMENSION V2(Nv,2,Ncol,ND),V1(Nv,2,Ncol,ND),BF(Nv)
C------------Main-------------------------------------------
      IF(Ipm.EQ.1) THEN
       DO 100 ic=1,Ncol
       DO 100 iv=1,Nv
        V2(iv,1,ic,1)=BF(iv)*(V1(iv,1,ic,1)+V1(iv,2,ic,4))
        V2(iv,2,ic,1)=BF(iv)*(V1(iv,2,ic,1)-V1(iv,1,ic,4))
        V2(iv,1,ic,2)=BF(iv)*(V1(iv,1,ic,2)+V1(iv,2,ic,3))
        V2(iv,2,ic,2)=BF(iv)*(V1(iv,2,ic,2)-V1(iv,1,ic,3))
        V2(iv,1,ic,3)=BF(iv)*(V1(iv,1,ic,3)-V1(iv,2,ic,2))
        V2(iv,2,ic,3)=BF(iv)*(V1(iv,2,ic,3)+V1(iv,1,ic,2))
        V2(iv,1,ic,4)=BF(iv)*(V1(iv,1,ic,4)-V1(iv,2,ic,1))
        V2(iv,2,ic,4)=BF(iv)*(V1(iv,2,ic,4)+V1(iv,1,ic,1))
 100   CONTINUE
      ELSE
       DO 200 ic=1,Ncol
       DO 200 iv=1,Nv
        V2(iv,1,ic,1)=BF(iv)*(V1(iv,1,ic,1)-V1(iv,2,ic,4))
        V2(iv,2,ic,1)=BF(iv)*(V1(iv,2,ic,1)+V1(iv,1,ic,4))
        V2(iv,1,ic,2)=BF(iv)*(V1(iv,1,ic,2)-V1(iv,2,ic,3))
        V2(iv,2,ic,2)=BF(iv)*(V1(iv,2,ic,2)+V1(iv,1,ic,3))
        V2(iv,1,ic,3)=BF(iv)*(V1(iv,1,ic,3)+V1(iv,2,ic,2))
        V2(iv,2,ic,3)=BF(iv)*(V1(iv,2,ic,3)-V1(iv,1,ic,2))
        V2(iv,1,ic,4)=BF(iv)*(V1(iv,1,ic,4)+V1(iv,2,ic,1))
        V2(iv,2,ic,4)=BF(iv)*(V1(iv,2,ic,4)-V1(iv,1,ic,1))
 200   CONTINUE
      ENDIF
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE G2Prd(Nv,Ipm,V2,V1,BF)
C          ---  V2 = BF * ( 1 +/- \gamma_2 ) * V1 ---
C
      INCLUDE 'lattsize.h'
C------------Variables--------------------------------------
      DIMENSION V2(Nv,2,Ncol,ND),V1(Nv,2,Ncol,ND),BF(Nv)
C------------Main-------------------------------------------
      IF(Ipm.EQ.1) THEN
       DO 100 ic=1,Ncol
       DO 100 iv=1,Nv
        V2(iv,1,ic,1)=BF(iv)*(V1(iv,1,ic,1)-V1(iv,1,ic,4))
        V2(iv,2,ic,1)=BF(iv)*(V1(iv,2,ic,1)-V1(iv,2,ic,4))
        V2(iv,1,ic,2)=BF(iv)*(V1(iv,1,ic,2)+V1(iv,1,ic,3))
        V2(iv,2,ic,2)=BF(iv)*(V1(iv,2,ic,2)+V1(iv,2,ic,3))
        V2(iv,1,ic,3)=BF(iv)*(V1(iv,1,ic,3)+V1(iv,1,ic,2))
        V2(iv,2,ic,3)=BF(iv)*(V1(iv,2,ic,3)+V1(iv,2,ic,2))
        V2(iv,1,ic,4)=BF(iv)*(V1(iv,1,ic,4)-V1(iv,1,ic,1))
        V2(iv,2,ic,4)=BF(iv)*(V1(iv,2,ic,4)-V1(iv,2,ic,1))
 100   CONTINUE
      ELSE
       DO 200 ic=1,Ncol
       DO 200 iv=1,Nv
        V2(iv,1,ic,1)=BF(iv)*(V1(iv,1,ic,1)+V1(iv,1,ic,4))
        V2(iv,2,ic,1)=BF(iv)*(V1(iv,2,ic,1)+V1(iv,2,ic,4))
        V2(iv,1,ic,2)=BF(iv)*(V1(iv,1,ic,2)-V1(iv,1,ic,3))
        V2(iv,2,ic,2)=BF(iv)*(V1(iv,2,ic,2)-V1(iv,2,ic,3))
        V2(iv,1,ic,3)=BF(iv)*(V1(iv,1,ic,3)-V1(iv,1,ic,2))
        V2(iv,2,ic,3)=BF(iv)*(V1(iv,2,ic,3)-V1(iv,2,ic,2))
        V2(iv,1,ic,4)=BF(iv)*(V1(iv,1,ic,4)+V1(iv,1,ic,1))
        V2(iv,2,ic,4)=BF(iv)*(V1(iv,2,ic,4)+V1(iv,2,ic,1))
 200   CONTINUE
      ENDIF
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE G3Prd(Nv,Ipm,V2,V1,BF)
C          ---  V2 = BF * ( 1 +/- \gamma_3 ) * V1 ---
C
      INCLUDE 'lattsize.h'
C------------Variables--------------------------------------
      DIMENSION V2(Nv,2,Ncol,ND),V1(Nv,2,Ncol,ND),BF(Nv)
C------------Main-------------------------------------------
      IF(Ipm.EQ.1) THEN
       DO 100 ic=1,Ncol
       DO 100 iv=1,Nv
        V2(iv,1,ic,1)=BF(iv)*(V1(iv,1,ic,1)+V1(iv,2,ic,3))
        V2(iv,2,ic,1)=BF(iv)*(V1(iv,2,ic,1)-V1(iv,1,ic,3))
        V2(iv,1,ic,2)=BF(iv)*(V1(iv,1,ic,2)-V1(iv,2,ic,4))
        V2(iv,2,ic,2)=BF(iv)*(V1(iv,2,ic,2)+V1(iv,1,ic,4))
        V2(iv,1,ic,3)=BF(iv)*(V1(iv,1,ic,3)-V1(iv,2,ic,1))
        V2(iv,2,ic,3)=BF(iv)*(V1(iv,2,ic,3)+V1(iv,1,ic,1))
        V2(iv,1,ic,4)=BF(iv)*(V1(iv,1,ic,4)+V1(iv,2,ic,2))
        V2(iv,2,ic,4)=BF(iv)*(V1(iv,2,ic,4)-V1(iv,1,ic,2))
 100   CONTINUE
      ELSE
       DO 200 ic=1,Ncol
       DO 200 iv=1,Nv
        V2(iv,1,ic,1)=BF(iv)*(V1(iv,1,ic,1)-V1(iv,2,ic,3))
        V2(iv,2,ic,1)=BF(iv)*(V1(iv,2,ic,1)+V1(iv,1,ic,3))
        V2(iv,1,ic,2)=BF(iv)*(V1(iv,1,ic,2)+V1(iv,2,ic,4))
        V2(iv,2,ic,2)=BF(iv)*(V1(iv,2,ic,2)-V1(iv,1,ic,4))
        V2(iv,1,ic,3)=BF(iv)*(V1(iv,1,ic,3)+V1(iv,2,ic,1))
        V2(iv,2,ic,3)=BF(iv)*(V1(iv,2,ic,3)-V1(iv,1,ic,1))
        V2(iv,1,ic,4)=BF(iv)*(V1(iv,1,ic,4)-V1(iv,2,ic,2))
        V2(iv,2,ic,4)=BF(iv)*(V1(iv,2,ic,4)+V1(iv,1,ic,2))
 200   CONTINUE
      ENDIF
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE G4Prd(Nv,Ipm,V2,V1)
C          ---  V2 = BF * ( 1 +/- \gamma_4 ) * V1 ---
C            ( No boundary factor )
C
      INCLUDE 'lattsize.h'
C------------Variables--------------------------------------
      DIMENSION V2(Nv,2,Ncol,ND),V1(Nv,2,Ncol,ND)
C------------Main-------------------------------------------
      IF(Ipm.EQ.1) THEN
       DO 100 ic=1,Ncol
       DO 100 iv=1,Nv
        V2(iv,1,ic,1) = V1(iv,1,ic,1) + V1(iv,1,ic,1)
        V2(iv,2,ic,1) = V1(iv,2,ic,1) + V1(iv,2,ic,1)
        V2(iv,1,ic,2) = V1(iv,1,ic,2) + V1(iv,1,ic,2)
        V2(iv,2,ic,2) = V1(iv,2,ic,2) + V1(iv,2,ic,2)
        V2(iv,1,ic,3) = V1(iv,1,ic,3) - V1(iv,1,ic,3)
        V2(iv,2,ic,3) = V1(iv,2,ic,3) - V1(iv,2,ic,3)
        V2(iv,1,ic,4) = V1(iv,1,ic,4) - V1(iv,1,ic,4)
        V2(iv,2,ic,4) = V1(iv,2,ic,4) - V1(iv,2,ic,4)
 100   CONTINUE
      ELSE
       DO 200 ic=1,Ncol
       DO 200 iv=1,Nv
        V2(iv,1,ic,1) = V1(iv,1,ic,1) - V1(iv,1,ic,1)
        V2(iv,2,ic,1) = V1(iv,2,ic,1) - V1(iv,2,ic,1)
        V2(iv,1,ic,2) = V1(iv,1,ic,2) - V1(iv,1,ic,2)
        V2(iv,2,ic,2) = V1(iv,2,ic,2) - V1(iv,2,ic,2)
        V2(iv,1,ic,3) = V1(iv,1,ic,3) + V1(iv,1,ic,3)
        V2(iv,2,ic,3) = V1(iv,2,ic,3) + V1(iv,2,ic,3)
        V2(iv,1,ic,4) = V1(iv,1,ic,4) + V1(iv,1,ic,4)
        V2(iv,2,ic,4) = V1(iv,2,ic,4) + V1(iv,2,ic,4)
 200   CONTINUE
      ENDIF
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
C***************************************************END*****
