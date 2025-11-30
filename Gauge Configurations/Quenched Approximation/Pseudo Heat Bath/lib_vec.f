C***********************************************************
C   Vector operators
C
C  .Need to include parameter file 'lattsize.h'
C       in which the precision is specified.
C  .The data size(Nv) must be giben. 
C  .Followngs are entry;
C       SetConst(Nv,V,A)         --  V = a
C       Equate(Nv,Y,W)           --  Y = W
C       SetVec(Nv,V,W,PRF)      --   V = prf * W
C       SelfAdd(Nv,V,W,PRF)      --  V = V + prf * W 
C       SelfSPrd(Nv,V,PRF)       --  V = prf * V
C       NormV(Nv,AV,V)           --  av = |V|^2
C       VecPrd(Nv,AV,V,W)        --  av = V * W
C
C                               07 Jun 1998  Matsufuru, H.
C              SetVec added     15 Jun 1999  M.,H.
C***********************************************************
      SUBROUTINE SetConst(Nv,V,A)
C                                     ---  V = a  ---
      INCLUDE 'lattsize.h'
C------------------Variables--------------------------------
      DIMENSION V(Nv) 
C--------------------Main-----------------------------------
C
      DO 100 i = 1,Nv
        V(i) = A
 100  CONTINUE
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE Equate(Nv,Y,W)
C                                   ---  Y = W  ---
      INCLUDE 'lattsize.h'
C------------------Variables--------------------------------
      DIMENSION Y(Nv),W(Nv)
C--------------------Main-----------------------------------
C
      DO 100 i = 1, Nv
           Y(i) = W(i)
 100  CONTINUE
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE SetVec(Nv,V,W,PRF)
C                             ---  V = prf * W  ---
      INCLUDE 'lattsize.h'
C-------------------Variables-------------------------------
      DIMENSION V(Nv),W(Nv) 
C---------------------Main----------------------------------
C
      DO 100 i = 1, Nv
        V(i) = PRF * W(i)
 100  CONTINUE
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE SelfAdd(Nv,V,W,PRF)
C                             ---  V = V + prf * W  ---
      INCLUDE 'lattsize.h'
C-------------------Variables-------------------------------
      DIMENSION V(Nv),W(Nv) 
C---------------------Main----------------------------------
C
      DO 100 i = 1, Nv
        V(i) = V(i) + PRF * W(i)
 100  CONTINUE
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE SelfSPrd(Nv,V,PRF)
C                                ---  V = prf * V  ---
      INCLUDE 'lattsize.h'
C------------------Variables--------------------------------
      DIMENSION V(Nv) 
C--------------------Main-----------------------------------
C
      DO 100 i = 1, Nv
        V(i) = PRF *V(i)
 100  CONTINUE
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE NormV(Nv,AV,V)
C                                  ---  av = |V|^2  ---
      INCLUDE 'lattsize.h'
C------------------Variables--------------------------------
      DIMENSION V(Nv) 
C--------------------Main-----------------------------------
C
      AV = 0.D0
C
      DO 100 i = 1, Nv
        AV = AV + V(i) * V(i)
 100   CONTINUE
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE VecPrd(Nv,AV,V,W)
C                                 ---  av = V * W  ---
      INCLUDE 'lattsize.h'
C-----------------Variables---------------------------------
      DIMENSION V(Nv),W(Nv) 
C-------------------Main------------------------------------
C
      AV = 0.D0
C
      DO 100 i=1, Nv
        AV = AV + V(i)*W(i)
 100   CONTINUE
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
C***************************************************END*****


