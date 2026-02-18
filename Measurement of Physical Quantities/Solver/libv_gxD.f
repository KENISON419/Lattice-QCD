C***********************************************************
C   \gamma matrices multiplication
C               ---  Dirac representation  ---
C
C      .For Vector processor
C                               2 Nov 1999  Matsufuru,H.
C***********************************************************
      SUBROUTINE AddG1xPrd(Nv,V2,V1,RW,FJD)
C          ---  V2 = \gamma_1  * V1 ---
C
      INCLUDE 'lattsize.h'
C------------Variables--------------------------------------
      DIMENSION V2(Nv,2,Ncol,ND),V1(Nv,2,Ncol,ND)
C------------Main-------------------------------------------
      DO 100 ic=1,Ncol
      DO 100 iv=1,Nv
       V2(iv,1,ic,1) = RW*V1(iv,1,ic,1) + FJD*V1(iv,2,ic,4)
       V2(iv,2,ic,1) = RW*V1(iv,2,ic,1) - FJD*V1(iv,1,ic,4)
       V2(iv,1,ic,2) = RW*V1(iv,1,ic,2) + FJD*V1(iv,2,ic,3)
       V2(iv,2,ic,2) = RW*V1(iv,2,ic,2) - FJD*V1(iv,1,ic,3)
       V2(iv,1,ic,3) = RW*V1(iv,1,ic,3) - FJD*V1(iv,2,ic,2)
       V2(iv,2,ic,3) = RW*V1(iv,2,ic,3) + FJD*V1(iv,1,ic,2)
       V2(iv,1,ic,4) = RW*V1(iv,1,ic,4) - FJD*V1(iv,2,ic,1)
       V2(iv,2,ic,4) = RW*V1(iv,2,ic,4) + FJD*V1(iv,1,ic,1)
 100  CONTINUE
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE AddG2xPrd(Nv,V2,V1,RW,FJD)
C          ---  V2 = \gamma_2  * V1 ---
C
      INCLUDE 'lattsize.h'
C------------Variables--------------------------------------
      DIMENSION V2(Nv,2,Ncol,ND),V1(Nv,2,Ncol,ND)
C------------Main-------------------------------------------
      DO 100 ic=1,Ncol
      DO 100 iv=1,Nv
       V2(iv,1,ic,1) = RW*V1(iv,1,ic,1) - FJD*V1(iv,1,ic,4)
       V2(iv,2,ic,1) = RW*V1(iv,2,ic,1) - FJD*V1(iv,2,ic,4)
       V2(iv,1,ic,2) = RW*V1(iv,1,ic,2) + FJD*V1(iv,1,ic,3)
       V2(iv,2,ic,2) = RW*V1(iv,2,ic,2) + FJD*V1(iv,2,ic,3)
       V2(iv,1,ic,3) = RW*V1(iv,1,ic,3) + FJD*V1(iv,1,ic,2)
       V2(iv,2,ic,3) = RW*V1(iv,2,ic,3) + FJD*V1(iv,2,ic,2)
       V2(iv,1,ic,4) = RW*V1(iv,1,ic,4) - FJD*V1(iv,1,ic,1)
       V2(iv,2,ic,4) = RW*V1(iv,2,ic,4) - FJD*V1(iv,2,ic,1)
 100  CONTINUE
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE AddG3xPrd(Nv,V2,V1,RW,FJD)
C          ---  V2 = \gamma_3  * V1 ---
C
      INCLUDE 'lattsize.h'
C------------Variables--------------------------------------
      DIMENSION V2(Nv,2,Ncol,ND),V1(Nv,2,Ncol,ND)
C------------Main-------------------------------------------
      DO 100 ic=1,Ncol
      DO 100 iv=1,Nv
       V2(iv,1,ic,1) = RW*V1(iv,1,ic,1) + FJD*V1(iv,2,ic,3)
       V2(iv,2,ic,1) = RW*V1(iv,2,ic,1) - FJD*V1(iv,1,ic,3)
       V2(iv,1,ic,2) = RW*V1(iv,1,ic,2) - FJD*V1(iv,2,ic,4)
       V2(iv,2,ic,2) = RW*V1(iv,2,ic,2) + FJD*V1(iv,1,ic,4)
       V2(iv,1,ic,3) = RW*V1(iv,1,ic,3) - FJD*V1(iv,2,ic,1)
       V2(iv,2,ic,3) = RW*V1(iv,2,ic,3) + FJD*V1(iv,1,ic,1)
       V2(iv,1,ic,4) = RW*V1(iv,1,ic,4) + FJD*V1(iv,2,ic,2)
       V2(iv,2,ic,4) = RW*V1(iv,2,ic,4) - FJD*V1(iv,1,ic,2)
 100  CONTINUE
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE AddG4xPrd(Nv,V2,V1,RW,FJD)
C          ---  V2 = \gamma_4  * V1 ---
C
      INCLUDE 'lattsize.h'
C------------Variables--------------------------------------
      DIMENSION V2(Nv,2,Ncol,ND),V1(Nv,2,Ncol,ND)
C------------Main-------------------------------------------
      DO 100 ic=1,Ncol
      DO 100 iv=1,Nv
       V2(iv,1,ic,1) = RW*V1(iv,1,ic,1) + FJD*V1(iv,1,ic,1)
       V2(iv,2,ic,1) = RW*V1(iv,2,ic,1) + FJD*V1(iv,2,ic,1)
       V2(iv,1,ic,2) = RW*V1(iv,1,ic,2) + FJD*V1(iv,1,ic,2)
       V2(iv,2,ic,2) = RW*V1(iv,2,ic,2) + FJD*V1(iv,2,ic,2)
       V2(iv,1,ic,3) = RW*V1(iv,1,ic,3) - FJD*V1(iv,1,ic,3)
       V2(iv,2,ic,3) = RW*V1(iv,2,ic,3) - FJD*V1(iv,2,ic,3)
       V2(iv,1,ic,4) = RW*V1(iv,1,ic,4) - FJD*V1(iv,1,ic,4)
       V2(iv,2,ic,4) = RW*V1(iv,2,ic,4) - FJD*V1(iv,2,ic,4)
 100  CONTINUE
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE SG23xPrd(Nv,V2,V1)
C          ---  V2 = \sigma_23 * V1 ---
C
      INCLUDE 'lattsize.h'
C------------Variables--------------------------------------
      DIMENSION V2(Nv,2,Ncol,ND),V1(Nv,2,Ncol,ND)
C------------Main-------------------------------------------
      DO 100 ic=1,Ncol
      DO 100 iv=1,Nv
        V2(iv,1,ic,1) =   V1(iv,1,ic,2)
        V2(iv,2,ic,1) =   V1(iv,2,ic,2)
        V2(iv,1,ic,2) =   V1(iv,1,ic,1)
        V2(iv,2,ic,2) =   V1(iv,2,ic,1)
        V2(iv,1,ic,3) =   V1(iv,1,ic,4)
        V2(iv,2,ic,3) =   V1(iv,2,ic,4)
        V2(iv,1,ic,4) =   V1(iv,1,ic,3)
        V2(iv,2,ic,4) =   V1(iv,2,ic,3)
 100  CONTINUE
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE SG31xPrd(Nv,V2,V1)
C          ---  V2 = \sigma_31 * V1 ---
C
      INCLUDE 'lattsize.h'
C------------Variables--------------------------------------
      DIMENSION V2(Nv,2,Ncol,ND),V1(Nv,2,Ncol,ND)
C------------Main-------------------------------------------
      DO 100 ic=1,Ncol
      DO 100 iv=1,Nv
        V2(iv,1,ic,1) =   V1(iv,2,ic,2)
        V2(iv,2,ic,1) = - V1(iv,1,ic,2)
        V2(iv,1,ic,2) = - V1(iv,2,ic,1)
        V2(iv,2,ic,2) =   V1(iv,1,ic,1)
        V2(iv,1,ic,3) =   V1(iv,2,ic,4)
        V2(iv,2,ic,3) = - V1(iv,1,ic,4)
        V2(iv,1,ic,4) = - V1(iv,2,ic,3)
        V2(iv,2,ic,4) =   V1(iv,1,ic,3)
 100  CONTINUE
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE SG12xPrd(Nv,V2,V1)
C          ---  V2 = \sigma_12 * V1 ---
C
      INCLUDE 'lattsize.h'
C------------Variables--------------------------------------
      DIMENSION V2(Nv,2,Ncol,ND),V1(Nv,2,Ncol,ND)
C------------Main-------------------------------------------
      DO 100 ic=1,Ncol
      DO 100 iv=1,Nv
        V2(iv,1,ic,1) =   V1(iv,1,ic,1)
        V2(iv,2,ic,1) =   V1(iv,2,ic,1)
        V2(iv,1,ic,2) = - V1(iv,1,ic,2)
        V2(iv,2,ic,2) = - V1(iv,2,ic,2)
        V2(iv,1,ic,3) =   V1(iv,1,ic,3)
        V2(iv,2,ic,3) =   V1(iv,2,ic,3)
        V2(iv,1,ic,4) = - V1(iv,1,ic,4)
        V2(iv,2,ic,4) = - V1(iv,2,ic,4)
 100  CONTINUE
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE SG14xPrd(Nv,V2,V1)
C          ---  V2 = \sigma_14 * V1 ---
C
      INCLUDE 'lattsize.h'
C------------Variables--------------------------------------
      DIMENSION V2(Nv,2,Ncol,ND),V1(Nv,2,Ncol,ND)
C------------Main-------------------------------------------
      DO 100 ic=1,Ncol
      DO 100 iv=1,Nv
        V2(iv,1,ic,1) = - V1(iv,1,ic,4)
        V2(iv,2,ic,1) = - V1(iv,2,ic,4)
        V2(iv,1,ic,2) = - V1(iv,1,ic,3)
        V2(iv,2,ic,2) = - V1(iv,2,ic,3)
        V2(iv,1,ic,3) = - V1(iv,1,ic,2)
        V2(iv,2,ic,3) = - V1(iv,2,ic,2)
        V2(iv,1,ic,4) = - V1(iv,1,ic,1)
        V2(iv,2,ic,4) = - V1(iv,2,ic,1)
 100  CONTINUE
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE SG24xPrd(Nv,V2,V1)
C          ---  V2 = \sigma_24 * V1 ---
C
      INCLUDE 'lattsize.h'
C------------Variables--------------------------------------
      DIMENSION V2(Nv,2,Ncol,ND),V1(Nv,2,Ncol,ND)
C------------Main-------------------------------------------
      DO 100 ic=1,Ncol
      DO 100 iv=1,Nv
        V2(iv,1,ic,1)= -V1(iv,2,ic,4)
        V2(iv,2,ic,1)=  V1(iv,1,ic,4)
        V2(iv,1,ic,2)=  V1(iv,2,ic,3)
        V2(iv,2,ic,2)= -V1(iv,1,ic,3)
        V2(iv,1,ic,3)= -V1(iv,2,ic,2)
        V2(iv,2,ic,3)=  V1(iv,1,ic,2)
        V2(iv,1,ic,4)=  V1(iv,2,ic,1)
        V2(iv,2,ic,4)= -V1(iv,1,ic,1)
 100  CONTINUE
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE SG34xPrd(Nv,V2,V1)
C          ---  V2 = \sigma_34 * V1 ---
C
      INCLUDE 'lattsize.h'
C------------Variables--------------------------------------
      DIMENSION V2(Nv,2,Ncol,ND),V1(Nv,2,Ncol,ND)
C------------Main-------------------------------------------
      DO 100 ic=1,Ncol
      DO 100 iv=1,Nv
        V2(iv,1,ic,1)= - V1(iv,1,ic,3)
        V2(iv,2,ic,1)= - V1(iv,2,ic,3)
        V2(iv,1,ic,2)=   V1(iv,1,ic,4)
        V2(iv,2,ic,2)=   V1(iv,2,ic,4)
        V2(iv,1,ic,3)= - V1(iv,1,ic,1)
        V2(iv,2,ic,3)= - V1(iv,2,ic,1)
        V2(iv,1,ic,4)=   V1(iv,1,ic,2)
        V2(iv,2,ic,4)=   V1(iv,2,ic,2)
 100  CONTINUE
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
C***************************************************END*****










