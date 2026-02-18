C***********************************************************
C  This is the common variables for light meson correlators
C      in quark solver program.
C
C                                 1998.06.07  Matsufuru,H.
C                      Modified   1998.08.28  Matsufuru,H.
C------------Variables--------------------------------------
      PARAMETER(Nqn=8)
      COMPLEX*16 Z1(NS,Ncol,ND,ND,NT)
       COMMON /LQP1/ Z1
      INTEGER    IGM(ND,0:9)
      COMPLEX*16 ZGM(ND,0:9)
         COMMON /GAMMA/ IGM,ZGM
      COMPLEX*16 ZWT(0:Nqn,NT),ZWZ(0:Nqn,0:Nz-1,NT)
      COMPLEX*16 ZWW(0:Nqn,Nz/2,NT)
         COMMON /MESON/ ZWT,ZWZ,ZWW
C***********************************************************
