C***********************************************************
C   Vector operation with Boundary factor
C
C  .Need to include parameter file 'lattsize.h'
C       in which the precision is specified.
C  .The data size(Nv1,Nv2) must be giben. 
C  .Followngs are entry;
C
C       VVPrd(Nv,V,W1,W2,PRF)
C          --  V(i1,i2) = V(i1,i2) + Prf * BF(i)*W2(i1,i2)
C       AddVVPrd(Nv,V,W1,W2,PRF)
C          --  V(i1,i2) =  Prf * BF(i)*W2(i1,i2)
C                               2 Nov 1999  Matsufuru,H.
C***********************************************************
      SUBROUTINE AddBFV(Nv1,Nv2,V,W,BF,PRF)
C          --  V(i1,i2) = V(i1,i2) + Prf * BF(i)*W(i1,i2)
      INCLUDE 'lattsize.h'
C------------------Variables--------------------------------
      DIMENSION V(Nv1,Nv2),W(Nv1,Nv2),BF(Nv1)
C--------------------Main-----------------------------------
C
      DO 100 i2 = 1, Nv2
C
      DO 100 i1 = 1, Nv1
C
        V(i1,i2) = V(i1,i2) + PRF * BF(i1) * W(i1,i2)
C
 100  CONTINUE
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE SetBFV(Nv1,Nv2,V,W,BF,PRF)
C                  --  V(i1,i2) =  Prf * BF(i)*W(i1,i2)
      INCLUDE 'lattsize.h'
C------------------Variables--------------------------------
      DIMENSION V(Nv1,Nv2),W(Nv1,Nv2),BF(Nv1)
C--------------------Main-----------------------------------
C
      DO 100 i2 = 1, Nv2
C
      DO 100 i1 = 1, Nv1
C
        V(i1,i2) = PRF * BF(i1) * W(i1,i2)
C
 100  CONTINUE
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
C***************************************************END*****


