C***********************************************************
C   SU(3) matrix production operators
C         --- For vector processor ---
C
C   .SetUnitv(Nv,A)         - Set A(Ndf,Nv) to unity
C   .TracevDN(Nv,TrB,A1,A2) - Tr[ A1^d * A2 ]
C   .UPrdvNN(Nv,B,A1,A2)       - B = A1   * A2
C   .UPrdvND(Nv,B,A1,A2)       - B = A1   * A2^d
C   .UPrdvDN(Nv,B,A1,A2)       - B = A1^d * A2
C   .Daggerv(Nv,A,B)           - B  = A^d
C
C                                24 Sep 1998  Matsufuru,H.
C        For vector prosessor    09 Jun 1999  M.,H.
C        `REUNITv' is moved to other dedicated file.
C                                12 Jun 1999  M.,H.
C***********************************************************
      SUBROUTINE SetUnitv(Nv,A)
C                                      ---  A = 1  ---
      INCLUDE 'lattsize.h'
C------------Variables--------------------------------------
      DIMENSION A(Nv,Ndf)
C------------Main-------------------------------------------
      DO 100 iv = 1, Nv
C
      A(iv, 1) = 1.0
      A(iv, 2) =  0.0
      A(iv, 3) =  0.0
      A(iv, 4) =  0.0
      A(iv, 5) =  0.0
      A(iv, 6) =  0.0
C
      A(iv, 7) =  0.0
      A(iv, 8) =  0.0
      A(iv, 9) = 1.0
      A(iv,10) =  0.0
      A(iv,11) =  0.0
      A(iv,12) =  0.0
C
      A(iv,13) =  0.0
      A(iv,14) =  0.0
      A(iv,15) =  0.0
      A(iv,16) =  0.0
      A(iv,17) = 1.0
      A(iv,18) =  0.0
C
 100  CONTINUE
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE TracevDN( Nv,TrB, A1,A2 ) 
C                  --- Tr B where B = A1^d * A2
      INCLUDE 'lattsize.h'
C------------Variables--------------------------------------
      DIMENSION A1(Nv,Ndf),A2(Nv,Ndf)
C------------Main-------------------------------------------
C
      TrB = 0.D0
C
      DO 100 idf=1,Ndf
      DO 100 iv =1,Nv
        TrB = TrB + A1(iv,idf) * A2(iv,idf)
 100  CONTINUE
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE UPrdvNN(Nv,B,A1,A2)
C                                   ---  B = A1 * A2  ---
      INCLUDE 'lattsize.h'
C------------Variables--------------------------------------
      DIMENSION B(Nv,Ndf),A1(Nv,Ndf),A2(Nv,Ndf)
C------------Main-------------------------------------------
      DO 100 iv=1,Nv
C
      B(iv, 1)= A1(iv, 1)*A2(iv, 1) -A1(iv, 2)*A2(iv, 2) 
     &         +A1(iv, 3)*A2(iv, 7) -A1(iv, 4)*A2(iv, 8) 
     &         +A1(iv, 5)*A2(iv,13) -A1(iv, 6)*A2(iv,14) 
C
      B(iv, 2)= A1(iv, 2)*A2(iv, 1) +A1(iv, 1)*A2(iv, 2) 
     &         +A1(iv, 4)*A2(iv, 7) +A1(iv, 3)*A2(iv, 8) 
     &         +A1(iv, 6)*A2(iv,13) +A1(iv, 5)*A2(iv,14) 
C
      B(iv, 3)= A1(iv, 1)*A2(iv, 3) -A1(iv, 2)*A2(iv, 4) 
     &         +A1(iv, 3)*A2(iv, 9) -A1(iv, 4)*A2(iv,10) 
     &         +A1(iv, 5)*A2(iv,15) -A1(iv, 6)*A2(iv,16) 
C
      B(iv, 4)= A1(iv, 2)*A2(iv, 3) +A1(iv, 1)*A2(iv, 4) 
     &         +A1(iv, 4)*A2(iv, 9) +A1(iv, 3)*A2(iv,10) 
     &         +A1(iv, 6)*A2(iv,15) +A1(iv, 5)*A2(iv,16) 
C
      B(iv, 5)= A1(iv, 1)*A2(iv, 5) -A1(iv, 2)*A2(iv, 6) 
     &         +A1(iv, 3)*A2(iv,11) -A1(iv, 4)*A2(iv,12) 
     &         +A1(iv, 5)*A2(iv,17) -A1(iv, 6)*A2(iv,18) 
C
      B(iv, 6)= A1(iv, 2)*A2(iv, 5) +A1(iv, 1)*A2(iv, 6) 
     &         +A1(iv, 4)*A2(iv,11) +A1(iv, 3)*A2(iv,12) 
     &         +A1(iv, 6)*A2(iv,17) +A1(iv, 5)*A2(iv,18) 
C---
      B(iv, 7)= A1(iv, 7)*A2(iv, 1) -A1(iv, 8)*A2(iv, 2) 
     &         +A1(iv, 9)*A2(iv, 7) -A1(iv,10)*A2(iv, 8) 
     &         +A1(iv,11)*A2(iv,13) -A1(iv,12)*A2(iv,14) 
C
      B(iv, 8)= A1(iv, 8)*A2(iv, 1) +A1(iv, 7)*A2(iv, 2) 
     &         +A1(iv,10)*A2(iv, 7) +A1(iv, 9)*A2(iv, 8) 
     &         +A1(iv,12)*A2(iv,13) +A1(iv,11)*A2(iv,14) 
C
      B(iv, 9)= A1(iv, 7)*A2(iv, 3) -A1(iv, 8)*A2(iv, 4) 
     &         +A1(iv, 9)*A2(iv, 9) -A1(iv,10)*A2(iv,10) 
     &         +A1(iv,11)*A2(iv,15) -A1(iv,12)*A2(iv,16) 
C
      B(iv,10)= A1(iv, 8)*A2(iv, 3) +A1(iv, 7)*A2(iv, 4) 
     &         +A1(iv,10)*A2(iv, 9) +A1(iv, 9)*A2(iv,10) 
     &         +A1(iv,12)*A2(iv,15) +A1(iv,11)*A2(iv,16) 
C
      B(iv,11)= A1(iv, 7)*A2(iv, 5) -A1(iv, 8)*A2(iv, 6) 
     &         +A1(iv, 9)*A2(iv,11) -A1(iv,10)*A2(iv,12) 
     &         +A1(iv,11)*A2(iv,17) -A1(iv,12)*A2(iv,18) 
C
      B(iv,12)= A1(iv, 8)*A2(iv, 5) +A1(iv, 7)*A2(iv, 6) 
     &         +A1(iv,10)*A2(iv,11) +A1(iv, 9)*A2(iv,12) 
     &         +A1(iv,12)*A2(iv,17) +A1(iv,11)*A2(iv,18) 
C---
      B(iv,13)= A1(iv,13)*A2(iv, 1) -A1(iv,14)*A2(iv, 2) 
     &         +A1(iv,15)*A2(iv, 7) -A1(iv,16)*A2(iv, 8) 
     &         +A1(iv,17)*A2(iv,13) -A1(iv,18)*A2(iv,14) 
C
      B(iv,14)= A1(iv,14)*A2(iv, 1) +A1(iv,13)*A2(iv, 2) 
     &         +A1(iv,16)*A2(iv, 7) +A1(iv,15)*A2(iv, 8) 
     &         +A1(iv,18)*A2(iv,13) +A1(iv,17)*A2(iv,14) 
C
      B(iv,15)= A1(iv,13)*A2(iv, 3) -A1(iv,14)*A2(iv, 4) 
     &         +A1(iv,15)*A2(iv, 9) -A1(iv,16)*A2(iv,10) 
     &         +A1(iv,17)*A2(iv,15) -A1(iv,18)*A2(iv,16) 
C
      B(iv,16)= A1(iv,14)*A2(iv, 3) +A1(iv,13)*A2(iv, 4) 
     &         +A1(iv,16)*A2(iv, 9) +A1(iv,15)*A2(iv,10) 
     &         +A1(iv,18)*A2(iv,15) +A1(iv,17)*A2(iv,16) 
C
      B(iv,17)= A1(iv,13)*A2(iv, 5) -A1(iv,14)*A2(iv, 6) 
     &         +A1(iv,15)*A2(iv,11) -A1(iv,16)*A2(iv,12) 
     &         +A1(iv,17)*A2(iv,17) -A1(iv,18)*A2(iv,18) 
C
      B(iv,18)= A1(iv,14)*A2(iv, 5) +A1(iv,13)*A2(iv, 6) 
     &         +A1(iv,16)*A2(iv,11) +A1(iv,15)*A2(iv,12) 
     &         +A1(iv,18)*A2(iv,17) +A1(iv,17)*A2(iv,18) 
C
 100  CONTINUE
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE UPrdvND(Nv,B,A1,A2)
C                                 ---  B = A1 * A2^d  ---
      INCLUDE 'lattsize.h'
C------------Variables--------------------------------------
      DIMENSION B(Nv,Ndf), A1(Nv,Ndf), A2(Nv,Ndf)
C------------Main-------------------------------------------
      DO 100 iv=1,Nv
C
      B(iv, 1)= A1(iv, 1)*A2(iv, 1) +A1(iv, 2)*A2(iv, 2) 
     &         +A1(iv, 3)*A2(iv, 3) +A1(iv, 4)*A2(iv, 4) 
     &         +A1(iv, 5)*A2(iv, 5) +A1(iv, 6)*A2(iv, 6) 
C
      B(iv, 2)= A1(iv, 2)*A2(iv, 1) -A1(iv, 1)*A2(iv, 2) 
     &         +A1(iv, 4)*A2(iv, 3) -A1(iv, 3)*A2(iv, 4) 
     &         +A1(iv, 6)*A2(iv, 5) -A1(iv, 5)*A2(iv, 6) 
C
      B(iv, 3)= A1(iv, 1)*A2(iv, 7) +A1(iv, 2)*A2(iv, 8) 
     &         +A1(iv, 3)*A2(iv, 9) +A1(iv, 4)*A2(iv,10) 
     &         +A1(iv, 5)*A2(iv,11) +A1(iv, 6)*A2(iv,12) 
C
      B(iv, 4)= A1(iv, 2)*A2(iv, 7) -A1(iv, 1)*A2(iv, 8) 
     &         +A1(iv, 4)*A2(iv, 9) -A1(iv, 3)*A2(iv,10) 
     &         +A1(iv, 6)*A2(iv,11) -A1(iv, 5)*A2(iv,12) 
C
      B(iv, 5)= A1(iv, 1)*A2(iv,13) +A1(iv, 2)*A2(iv,14) 
     &         +A1(iv, 3)*A2(iv,15) +A1(iv, 4)*A2(iv,16) 
     &         +A1(iv, 5)*A2(iv,17) +A1(iv, 6)*A2(iv,18) 
C
      B(iv, 6)= A1(iv, 2)*A2(iv,13) -A1(iv, 1)*A2(iv,14) 
     &         +A1(iv, 4)*A2(iv,15) -A1(iv, 3)*A2(iv,16) 
     &         +A1(iv, 6)*A2(iv,17) -A1(iv, 5)*A2(iv,18) 
C---
      B(iv, 7)= A1(iv, 7)*A2(iv, 1) +A1(iv, 8)*A2(iv, 2) 
     &         +A1(iv, 9)*A2(iv, 3) +A1(iv,10)*A2(iv, 4) 
     &         +A1(iv,11)*A2(iv, 5) +A1(iv,12)*A2(iv, 6) 
C
      B(iv, 8)= A1(iv, 8)*A2(iv, 1) -A1(iv, 7)*A2(iv, 2) 
     &         +A1(iv,10)*A2(iv, 3) -A1(iv, 9)*A2(iv, 4) 
     &         +A1(iv,12)*A2(iv, 5) -A1(iv,11)*A2(iv, 6) 
C
      B(iv, 9)= A1(iv, 7)*A2(iv, 7) +A1(iv, 8)*A2(iv, 8) 
     &         +A1(iv, 9)*A2(iv, 9) +A1(iv,10)*A2(iv,10) 
     &         +A1(iv,11)*A2(iv,11) +A1(iv,12)*A2(iv,12) 
C
      B(iv,10)= A1(iv, 8)*A2(iv, 7) -A1(iv, 7)*A2(iv, 8) 
     &         +A1(iv,10)*A2(iv, 9) -A1(iv, 9)*A2(iv,10) 
     &         +A1(iv,12)*A2(iv,11) -A1(iv,11)*A2(iv,12) 
C
      B(iv,11)= A1(iv, 7)*A2(iv,13) +A1(iv, 8)*A2(iv,14) 
     &         +A1(iv, 9)*A2(iv,15) +A1(iv,10)*A2(iv,16) 
     &         +A1(iv,11)*A2(iv,17) +A1(iv,12)*A2(iv,18) 
C
      B(iv,12)= A1(iv, 8)*A2(iv,13) -A1(iv, 7)*A2(iv,14) 
     &         +A1(iv,10)*A2(iv,15) -A1(iv, 9)*A2(iv,16) 
     &         +A1(iv,12)*A2(iv,17) -A1(iv,11)*A2(iv,18) 
C---
      B(iv,13)= A1(iv,13)*A2(iv, 1) +A1(iv,14)*A2(iv, 2) 
     &         +A1(iv,15)*A2(iv, 3) +A1(iv,16)*A2(iv, 4) 
     &         +A1(iv,17)*A2(iv, 5) +A1(iv,18)*A2(iv, 6) 
C
      B(iv,14)= A1(iv,14)*A2(iv, 1) -A1(iv,13)*A2(iv, 2) 
     &         +A1(iv,16)*A2(iv, 3) -A1(iv,15)*A2(iv, 4) 
     &         +A1(iv,18)*A2(iv, 5) -A1(iv,17)*A2(iv, 6) 
C
      B(iv,15)= A1(iv,13)*A2(iv, 7) +A1(iv,14)*A2(iv, 8) 
     &         +A1(iv,15)*A2(iv, 9) +A1(iv,16)*A2(iv,10) 
     &         +A1(iv,17)*A2(iv,11) +A1(iv,18)*A2(iv,12) 
C
      B(iv,16)= A1(iv,14)*A2(iv, 7) -A1(iv,13)*A2(iv, 8) 
     &         +A1(iv,16)*A2(iv, 9) -A1(iv,15)*A2(iv,10) 
     &         +A1(iv,18)*A2(iv,11) -A1(iv,17)*A2(iv,12) 
C
      B(iv,17)= A1(iv,13)*A2(iv,13) +A1(iv,14)*A2(iv,14)
     &         +A1(iv,15)*A2(iv,15) +A1(iv,16)*A2(iv,16)
     &         +A1(iv,17)*A2(iv,17) +A1(iv,18)*A2(iv,18)
C
      B(iv,18)= A1(iv,14)*A2(iv,13) -A1(iv,13)*A2(iv,14) 
     &         +A1(iv,16)*A2(iv,15) -A1(iv,15)*A2(iv,16) 
     &         +A1(iv,18)*A2(iv,17) -A1(iv,17)*A2(iv,18) 
C
 100  CONTINUE
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE UPrdvDN(Nv,B,A1,A2)
C                                 ---  B = A1^d * A2  ---
      INCLUDE 'lattsize.h'
C------------Variables--------------------------------------
      DIMENSION B(Nv,Ndf),A1(Nv,Ndf),A2(Nv,Ndf)
C------------Main-------------------------------------------
      DO 100 iv = 1, Nv
C
      B(iv, 1)= A1(iv, 1)*A2(iv, 1) +A1(iv, 2)*A2(iv, 2) 
     &         +A1(iv, 7)*A2(iv, 7) +A1(iv, 8)*A2(iv, 8) 
     &         +A1(iv,13)*A2(iv,13) +A1(iv,14)*A2(iv,14) 
C
      B(iv, 2)=-A1(iv, 2)*A2(iv, 1) +A1(iv, 1)*A2(iv, 2) 
     &         -A1(iv, 8)*A2(iv, 7) +A1(iv, 7)*A2(iv, 8) 
     &         -A1(iv,14)*A2(iv,13) +A1(iv,13)*A2(iv,14) 
C
      B(iv, 3)= A1(iv, 1)*A2(iv, 3) +A1(iv, 2)*A2(iv, 4) 
     &         +A1(iv, 7)*A2(iv, 9) +A1(iv, 8)*A2(iv,10) 
     &         +A1(iv,13)*A2(iv,15) +A1(iv,14)*A2(iv,16) 
C
      B(iv, 4)=-A1(iv, 2)*A2(iv, 3) +A1(iv, 1)*A2(iv, 4) 
     &         -A1(iv, 8)*A2(iv, 9) +A1(iv, 7)*A2(iv,10) 
     &         -A1(iv,14)*A2(iv,15) +A1(iv,13)*A2(iv,16) 
C
      B(iv, 5)= A1(iv, 1)*A2(iv, 5) +A1(iv, 2)*A2(iv, 6) 
     &         +A1(iv, 7)*A2(iv,11) +A1(iv, 8)*A2(iv,12) 
     &         +A1(iv,13)*A2(iv,17) +A1(iv,14)*A2(iv,18) 
C
      B(iv, 6)=-A1(iv, 2)*A2(iv, 5) +A1(iv, 1)*A2(iv, 6) 
     &         -A1(iv, 8)*A2(iv,11) +A1(iv, 7)*A2(iv,12) 
     &         -A1(iv,14)*A2(iv,17) +A1(iv,13)*A2(iv,18) 
C---
      B(iv, 7)= A1(iv, 3)*A2(iv, 1) +A1(iv, 4)*A2(iv, 2) 
     &         +A1(iv, 9)*A2(iv, 7) +A1(iv,10)*A2(iv, 8) 
     &         +A1(iv,15)*A2(iv,13) +A1(iv,16)*A2(iv,14) 
C
      B(iv, 8)=-A1(iv, 4)*A2(iv, 1) +A1(iv, 3)*A2(iv, 2) 
     &         -A1(iv,10)*A2(iv, 7) +A1(iv, 9)*A2(iv, 8) 
     &         -A1(iv,16)*A2(iv,13) +A1(iv,15)*A2(iv,14) 
C
      B(iv, 9)= A1(iv, 3)*A2(iv, 3) +A1(iv, 4)*A2(iv, 4) 
     &         +A1(iv, 9)*A2(iv, 9) +A1(iv,10)*A2(iv,10) 
     &         +A1(iv,15)*A2(iv,15) +A1(iv,16)*A2(iv,16) 
C
      B(iv,10)=-A1(iv, 4)*A2(iv, 3) +A1(iv, 3)*A2(iv, 4) 
     &         -A1(iv,10)*A2(iv, 9) +A1(iv, 9)*A2(iv,10) 
     &         -A1(iv,16)*A2(iv,15) +A1(iv,15)*A2(iv,16) 
C
      B(iv,11)= A1(iv, 3)*A2(iv, 5) +A1(iv, 4)*A2(iv, 6) 
     &         +A1(iv, 9)*A2(iv,11) +A1(iv,10)*A2(iv,12) 
     &         +A1(iv,15)*A2(iv,17) +A1(iv,16)*A2(iv,18) 
C
      B(iv,12)=-A1(iv, 4)*A2(iv, 5) +A1(iv, 3)*A2(iv, 6) 
     &         -A1(iv,10)*A2(iv,11) +A1(iv, 9)*A2(iv,12) 
     &         -A1(iv,16)*A2(iv,17) +A1(iv,15)*A2(iv,18) 
C---
      B(iv,13)= A1(iv, 5)*A2(iv, 1) +A1(iv, 6)*A2(iv, 2) 
     &         +A1(iv,11)*A2(iv, 7) +A1(iv,12)*A2(iv, 8) 
     &         +A1(iv,17)*A2(iv,13) +A1(iv,18)*A2(iv,14) 
C
      B(iv,14)=-A1(iv, 6)*A2(iv, 1) +A1(iv, 5)*A2(iv, 2) 
     &         -A1(iv,12)*A2(iv, 7) +A1(iv,11)*A2(iv, 8) 
     &         -A1(iv,18)*A2(iv,13) +A1(iv,17)*A2(iv,14) 
C
      B(iv,15)= A1(iv, 5)*A2(iv, 3) +A1(iv, 6)*A2(iv, 4)
     &         +A1(iv,11)*A2(iv, 9) +A1(iv,12)*A2(iv,10) 
     &         +A1(iv,17)*A2(iv,15) +A1(iv,18)*A2(iv,16) 
C
      B(iv,16)=-A1(iv, 6)*A2(iv, 3) +A1(iv, 5)*A2(iv, 4) 
     &         -A1(iv,12)*A2(iv, 9) +A1(iv,11)*A2(iv,10) 
     &         -A1(iv,18)*A2(iv,15) +A1(iv,17)*A2(iv,16) 
C
      B(iv,17)= A1(iv, 5)*A2(iv, 5) +A1(iv, 6)*A2(iv, 6) 
     &         +A1(iv,11)*A2(iv,11) +A1(iv,12)*A2(iv,12) 
     &         +A1(iv,17)*A2(iv,17) +A1(iv,18)*A2(iv,18) 
C
      B(iv,18)=-A1(iv, 6)*A2(iv, 5) +A1(iv, 5)*A2(iv, 6) 
     &         -A1(iv,12)*A2(iv,11) +A1(iv,11)*A2(iv,12) 
     &         -A1(iv,18)*A2(iv,17) +A1(iv,17)*A2(iv,18) 
C
 100  CONTINUE
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE Daggerv( Nv, B, A )
C                             ------    B = A^d    ------
C
      INCLUDE 'lattsize.h'
C------------Variables--------------------------------------
      DIMENSION A(Nv,Ndf), B(Nv,Ndf)
C------------Main-------------------------------------------
      DO 100 iv=1,Nv
C
      B( iv,  1 ) =  A( iv,  1 )
      B( iv,  2 ) = -A( iv,  2 )
      B( iv,  3 ) =  A( iv,  7 )
      B( iv,  4 ) = -A( iv,  8 )
      B( iv,  5 ) =  A( iv, 13 )
      B( iv,  6 ) = -A( iv, 14 )
      B( iv,  7 ) =  A( iv,  3 )
      B( iv,  8 ) = -A( iv,  4 )
      B( iv,  9 ) =  A( iv,  9 )
      B( iv, 10 ) = -A( iv, 10 )
      B( iv, 11 ) =  A( iv, 15 )
      B( iv, 12 ) = -A( iv, 16 )
      B( iv, 13 ) =  A( iv,  5 )
      B( iv, 14 ) = -A( iv,  6 )
      B( iv, 15 ) =  A( iv, 11 )
      B( iv, 16 ) = -A( iv, 12 )
      B( iv, 17 ) =  A( iv, 17 )
      B( iv, 18 ) = -A( iv, 18 )
C
 100  CONTINUE
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
C***************************************************END*****


