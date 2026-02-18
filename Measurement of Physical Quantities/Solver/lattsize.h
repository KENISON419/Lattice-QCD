C***********************************************************
C  This is the lattice size parameters
C               and definitions for convenience.
C                                 1998.03.12  Matsufuru,H.
C-----------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
C-----------------------------------------------------------
      PARAMETER ( Nsize=12, Ntime=24 )
      PARAMETER ( Nx=Nsize, Ny=Nsize, Nz=Nsize )
      PARAMETER ( NPE=1, NT=Ntime/NPE )
      PARAMETER ( Nsite=Nx*Ny*Nz, Nlink=Nsite*4 )
      PARAMETER ( NS=Nsite,NS2=Nsite/2,Nx2=Nx/2 )
      PARAMETER ( Ncol=3,Nvc=Ncol*2 )
      PARAMETER ( Ndf=Ncol*Ncol*2, Ncp=Ncol*Ncol )
      PARAMETER ( ND = 4, ND2 = 2 )
      PARAMETER ( PI = 3. 14159 26535 89793 )
      PARAMETER ( NBint  = 4, NBreal = 8 )
C              (NBreal=4 for single & 8 for double prec.)
C***********************************************************

