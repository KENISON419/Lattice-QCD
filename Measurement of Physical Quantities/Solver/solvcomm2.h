C***********************************************************
C  This is the common variables for quark solver (type-2)
C        on the anisotropic lattice.
C
C                                28 Aug 1998  Matsufuru,H.
C          For vector processor  12 Jun 1999  M.,H.
C------------Variables--------------------------------------
      COMMON /START/  Iconf,Ikapp,NBC(4),Isrc,NtimeS,
     &                Niter,Istart
      COMMON /SPARAM/ CKs,gammaF,Enorm 
C
      COMMON /LINKV/  U( NS,Ndf,4,NT)
      COMMON /LISTV/  LV( NS,3,2 )
      COMMON /LISTB/  BF( NS,3,2 )
C
      COMMON /SRCFNC/ Bsrc(NS)
      COMMON /FVECT1/ X( NS,Nvc,ND,NT) 
      COMMON /FVECT2/ R( NS,Nvc,ND,NT) 
      COMMON /FVECT3/ S( NS,Nvc,ND,NT) 
C
C***********************************************************
