C***********************************************************
C  This is the common variables for pseudo-heat bath
C      quenched lattice update program.
C
C                                12 Mar 1998  Matsufuru,H.
C      .Modified for vectoriz.   14 Jul 1999  Matsufuru,H.
C------------Variables--------------------------------------
      COMMON /START/  Ilatt, Iconf, Istart,
     &                Nsweep, Nprev, Nhit, Iplaq
      COMMON /SPARAM/ Beta, Gamma
C                                      (Parameters)
      COMMON /LINKV/  U(NS2,Ndf,2,4,NT)
C                                      (Link var.)
      COMMON /LISTV/  LV(NS2,3,2,2)
C                                      (List vector)
      COMMON /STAPLE/ C(NS2,Ndf)
C                                      (Staple)
C***********************************************************
