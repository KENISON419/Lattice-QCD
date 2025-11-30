C***********************************************************
C   SU(3) matrix reunitarization
C         --- (e,o) version, for vector processor ---
C
C                                12 Jun 1999  Matsufuru,H.
C***********************************************************
      SUBROUTINE REUNIT2v(Nv,G)
C  ----- Schmidt unitarization of SU(3) matrices -----
C                                             Matsufuru,H.
C      Modified (real variable)   04 Sep 1998 Matsufuru,H.
C      Vectorized                 08 Jun 1999 Matsufuru,H.
C      Further vectorized         12 Jun 1999 Matsufuru,H.
C
      INCLUDE 'lattsize.h'
C------------Variables--------------------------------------
      DIMENSION G(NS2,Ndf,Nv),F(NS2,Ndf)
C------------Main-------------------------------------------
C
      DO 5000 iv = 1, Nv
C
      DO 1000 is = 1, NS2
C
C --- Normalize 1st row vector ---
      SN1 = 1.D0/DSQRT( G(is,1,iv)**2 +G(is,2,iv)**2
     &                 +G(is,3,iv)**2 +G(is,4,iv)**2
     &                 +G(is,5,iv)**2 +G(is,6,iv)**2 )
      DO 100 idf = 1,6
        F(is,idf)=G(is,idf,iv)*SN1
 100  CONTINUE
C
C --- Orthogonalize 2nd vector to 1st vector ---
      SP1r =   F(is, 1)*G(is, 7,iv) + F(is, 2)*G(is, 8,iv)
     &       + F(is, 3)*G(is, 9,iv) + F(is, 4)*G(is,10,iv)
     &       + F(is, 5)*G(is,11,iv) + F(is, 6)*G(is,12,iv)
      SP1i =   F(is, 1)*G(is, 8,iv) - F(is, 2)*G(is, 7,iv)
     &       + F(is, 3)*G(is,10,iv) - F(is, 4)*G(is, 9,iv)
     &       + F(is, 5)*G(is,12,iv) - F(is, 6)*G(is,11,iv)
C
      F(is, 7) = G(is, 7,iv) - SP1r*F(is, 1) + SP1i*F(is, 2)
      F(is, 8) = G(is, 8,iv) - SP1r*F(is, 2) - SP1i*F(is, 1)
      F(is, 9) = G(is, 9,iv) - SP1r*F(is, 3) + SP1i*F(is, 4)
      F(is,10) = G(is,10,iv) - SP1r*F(is, 4) - SP1i*F(is, 3)
      F(is,11) = G(is,11,iv) - SP1r*F(is, 5) + SP1i*F(is, 6)
      F(is,12) = G(is,12,iv) - SP1r*F(is, 6) - SP1i*F(is, 5)

C --- Normalize the 2nd row vector ---
      SN2 = 1.0/DSQRT(   F(is, 7)**2 + F(is, 8)**2 + F(is, 9)**2
     &                 + F(is,10)**2 + F(is,11)**2 + F(is,12)**2 )
      DO 200 idf = 7,12
        F(is,idf)=F(is,idf)*SN2
 200  CONTINUE
C
C --- Orthogonalize 3rd vector to 1st and 2nd vector ---
      SP2r =   F(is, 1)*G(is,13,iv) + F(is, 2)*G(is,14,iv)
     &       + F(is, 3)*G(is,15,iv) + F(is, 4)*G(is,16,iv)
     &       + F(is, 5)*G(is,17,iv) + F(is, 6)*G(is,18,iv)
      SP2i =   F(is, 1)*G(is,14,iv) - F(is, 2)*G(is,13,iv)
     &       + F(is, 3)*G(is,16,iv) - F(is, 4)*G(is,15,iv)
     &       + F(is, 5)*G(is,18,iv) - F(is, 6)*G(is,17,iv)
C
      SP3r =   F(is, 7)*G(is,13,iv) + F(is, 8)*G(is,14,iv)
     &       + F(is, 9)*G(is,15,iv) + F(is,10)*G(is,16,iv)
     &       + F(is,11)*G(is,17,iv) + F(is,12)*G(is,18,iv)
      SP3i =   F(is, 7)*G(is,14,iv) - F(is, 8)*G(is,13,iv)
     &       + F(is, 9)*G(is,16,iv) - F(is,10)*G(is,15,iv)
     &       + F(is,11)*G(is,18,iv) - F(is,12)*G(is,17,iv)
C
      F(is,13) = G(is,13,iv) - SP2r*F(is, 1) + SP2i*F(is, 2)
     &                       - SP3r*F(is, 7) + SP3i*F(is, 8)
      F(is,14) = G(is,14,iv) - SP2r*F(is, 2) - SP2i*F(is, 1)
     &                       - SP3r*F(is, 8) - SP3i*F(is, 7)
      F(is,15) = G(is,15,iv) - SP2r*F(is, 3) + SP2i*F(is, 4)
     &                       - SP3r*F(is, 9) + SP3i*F(is,10)
      F(is,16) = G(is,16,iv) - SP2r*F(is, 4) - SP2i*F(is, 3)
     &                       - SP3r*F(is,10) - SP3i*F(is, 9)
      F(is,17) = G(is,17,iv) - SP2r*F(is, 5) + SP2i*F(is, 6)
     &                       - SP3r*F(is,11) + SP3i*F(is,12)
      F(is,18) = G(is,18,iv) - SP2r*F(is, 6) - SP2i*F(is, 5)
     &                       - SP3r*F(is,12) - SP3i*F(is,11)
C
C --- Normalize 3rd row vector ---
      SN3 = 1.D0/DSQRT(   F(is,13)**2 + F(is,14)**2 + F(is,15)**2
     &                  + F(is,16)**2 + F(is,17)**2 + F(is,18)**2 )
      DO 300 idf = 13,18
        F(is,idf)=F(is,idf)*SN3
 300  CONTINUE
C
C --- Determinant and phase factor ---
C
      DETr= F(is, 1)*( F(is, 9)*F(is,17) - F(is,10)*F(is,18) )
     &     -F(is, 2)*( F(is, 9)*F(is,18) + F(is,10)*F(is,17) )
     &    + F(is, 3)*( F(is,11)*F(is,13) - F(is,12)*F(is,14) )
     &     -F(is, 4)*( F(is,11)*F(is,14) + F(is,12)*F(is,13) )
     &    + F(is, 7)*( F(is,15)*F(is, 5) - F(is,16)*F(is, 6) )
     &     -F(is, 8)*( F(is,15)*F(is, 6) + F(is,16)*F(is, 5) )
     &    - F(is, 5)*( F(is, 9)*F(is,13) - F(is,10)*F(is,14) )
     &     +F(is, 6)*( F(is, 9)*F(is,14) + F(is,10)*F(is,13) )
     &    - F(is, 3)*( F(is, 7)*F(is,17) - F(is, 8)*F(is,18) )
     &     +F(is, 4)*( F(is, 7)*F(is,18) + F(is, 8)*F(is,17) )
     &    - F(is,11)*( F(is,15)*F(is, 1) - F(is,16)*F(is, 2) )
     &     +F(is,12)*( F(is,15)*F(is, 2) + F(is,16)*F(is, 1) )
C
      DETi= F(is, 1)*( F(is, 9)*F(is,18) + F(is,10)*F(is,17) )
     &     +F(is, 2)*( F(is, 9)*F(is,17) - F(is,10)*F(is,18) )
     &    + F(is, 3)*( F(is,11)*F(is,14) + F(is,12)*F(is,13) )
     &     +F(is, 4)*( F(is,11)*F(is,13) - F(is,12)*F(is,14) )
     &    + F(is, 7)*( F(is,15)*F(is, 6) + F(is,16)*F(is, 5) )
     &     +F(is, 8)*( F(is,15)*F(is, 5) - F(is,16)*F(is, 6) )
     &    - F(is, 5)*( F(is, 9)*F(is,14) + F(is,10)*F(is,13) )
     &     -F(is, 6)*( F(is, 9)*F(is,13) - F(is,10)*F(is,14) )
     &    - F(is, 3)*( F(is, 7)*F(is,18) + F(is, 8)*F(is,17) )
     &     -F(is, 4)*( F(is, 7)*F(is,17) - F(is, 8)*F(is,18) )
     &    - F(is,11)*( F(is,15)*F(is, 2) + F(is,16)*F(is, 1) )
     &     -F(is,12)*( F(is,15)*F(is, 1) - F(is,16)*F(is, 2) )
C
      Gnorm = 1.0/DSQRT( DETr**2 + DETi**2 )
      Arg = DATAN2( Gnorm*DETi, Gnorm*DETr )
      PHr = DCOS( -Arg/3.D0 )
      PHi = DSIN( -Arg/3.D0 )
C
C --- Set the determinant to unity ---
      DO 500 icp = 1,Ncp
       G(is,2*icp-1,iv) = PHr*F(is,2*icp-1) -PHi*F(is,2*icp  )
       G(is,2*icp  ,iv) = PHr*F(is,2*icp  ) +PHi*F(is,2*icp-1)
  500 CONTINUE
C
 1000 CONTINUE
C
 5000 CONTINUE
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
C***************************************************END*****


