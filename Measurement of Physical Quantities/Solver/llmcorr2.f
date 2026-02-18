C                   < llmcorr2.f >
C***********************************************************
      SUBROUTINE LLMCORR
C         ----- Meson correlator and wave function -----
C
C    .Version for use in quark solver          12 May 1998
C    .Ver.2.0: with wave function              28 Mar 1999
C
C                      Copyright (c)  Hideo Matsufuru 1999
C***********************************************************
      INCLUDE 'lattsize.h'
      INCLUDE 'llmcomm2.h'
C----------------------Variables----------------------------
      COMPLEX*16 ZWZT(0:Nz-1,Nt)
C-----------------Meson correlator--------------------------
C
      DO 5000 iqn = 0, Nqn
C
      IF(IQN.EQ.0) THEN
       IGM0 = 0
       IGMx = 0
      ELSE IF(IQN.EQ.1) THEN
       IGM0 = 6
       IGMx = 6
      ELSE IF(IQN.EQ.2) THEN
       IGM0 = 7
       IGMx = 7
      ELSE IF(IQN.EQ.3) THEN
       IGM0 = 8
       IGMx = 8
      ELSE IF(IQN.EQ.4) THEN
       IGM0 = 0
       IGMx = 4
      ELSE IF(IQN.EQ.5) THEN
       IGM0 = 5
       IGMx = 5
      ELSE IF(IQN.EQ.6) THEN
       IGM0 = 1
       IGMx = 1
      ELSE IF(IQN.EQ.7) THEN
       IGM0 = 2
       IGMx = 2
      ELSE IF(IQN.EQ.8) THEN
       IGM0 = 3
       IGMx = 3
      END IF
C
C ---Correlator in z-direction---
      DO 1000 it = 1,Nt
      DO 1000 iz = 0,Nz-1
        ZWZT( iz,it ) = DCMPLX( 0.D0, 0.D0 )
 1000 CONTINUE
C
      DO 1500 it = 1, Nt
      DO 1500 iz = 0, Nz-1
C
       DO 1600 id0 = 1,ND
C      ( Souce-color(IC0) loop is not here )
       DO 1600 idx = 1,ND
       DO 1600 icx = 1,Ncol
       DO 1600 ixy = 1, Nx*Ny
          is = ixy + iz *Nx*Ny
         ZWZT(iz,it) = ZWZT(iz,it)
     &           + ZGM(idx,IGMx)*ZGM(id0,IGM0)
     &            *Z1( is,icx,IGM(idx,IGMx),id0,it)
     &            *DCONJG( Z1(is,icx,idx,IGM(id0,IGM0),it) )
 1600  CONTINUE
C
 1500 CONTINUE
C
      DO 1800 it = 1,Nt
      DO 1800 iz = 0,Nz-1
         ZWT(iqn,   it) = ZWT(iqn,   it) + ZWZT(iz,it)
         ZWZ(iqn,iz,it) = ZWZ(iqn,iz,it) + ZWZT(iz,it)
 1800 CONTINUE
C
C ---Wave function---
      DO 2500 it = 1, Nt
      DO 2500 ir = 1, Nz/2
C
       DO 2600 id0 = 1,ND
C      ( Souce-color(IC0) loop is not here )
       DO 2600 idx = 1,ND
       DO 2600 icx = 1,Ncol
       DO 2600 iz = 0, Nz-1
          izr = MOD(iz+ir,Nz)
       DO 2600 ixy = 1, Nx*Ny
          is  = ixy +  iz  *Nx*Ny
          isr = ixy +  izr *Nx*Ny
         ZWW(iqn,ir,it) = ZWW(iqn,ir,it)
     &          + ZGM(idx,IGMx)*ZGM(id0,IGM0)
     &           *Z1( is,icx,IGM(idx,IGMx),id0,it)
     &           *CONJG( Z1(isr,icx,idx,IGM(id0,IGM0),it) )
 2600 CONTINUE
C
 2500 CONTINUE
C
C ---iqn loop end---
 5000 CONTINUE
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE SetZ1(IDS,XN)
C
      INCLUDE 'lattsize.h'
      INCLUDE 'llmcomm2.h'
C------------Variables--------------------------------------
      DIMENSION XN(NS,Nvc,ND,NT) 
C------------Main-------------------------------------------
C
      DO 2800 it  = 1, Nt
      DO 2800 id  = 1, ND
      DO 2800 ic  = 1, Ncol
      DO 2800 is  = 1, NS
        Z1(is,ic,id,IDS,it)
     &        = CMPLX( XN(is,2*ic-1,id,it),
     &                 XN(is, 2*ic ,id,it) )
 2800  CONTINUE
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE LLMinit
C
      INCLUDE 'lattsize.h'
      INCLUDE 'llmcomm2.h'
C------------Main-------------------------------------------
C ---Gamma matrices initial set---
      CALL SETGM(IGM,ZGM)
C
      DO 100 iqn = 0,Nqn
      DO 100 it = 1,NT
C
          ZWT(iqn,it) = DCMPLX( 0.D0, 0.D0 )
C
       DO 120 iz = 0,Nz-1
          ZWZ(iqn,iz,it) = DCMPLX( 0.D0, 0.D0 )
 120   CONTINUE
C
       DO 140 ir = 1,Nz/2
          ZWW(iqn,ir,it) = DCMPLX( 0.D0, 0.D0 )
 140   CONTINUE
C
 100  CONTINUE
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE WRITEW(Iconf, Ikapp, NODEID)
C
      INCLUDE 'lattsize.h'
      INCLUDE 'llmcomm2.h'
C-------------------Variables-------------------------------
      CHARACTER*120 Fname
      CHARACTER*4 CHconf
C----------------------Main---------------------------------
C
      CHconf = '0' //CHAR( ICHAR('0')+    Iconf/100 )
     &             //CHAR( ICHAR('0')+MOD(Iconf/10,10) )
     &             //CHAR( ICHAR('0')+MOD(Iconf,10) )
C
      DO 2000 itime = 0,Ntime-1
        it    = itime + 1
C
C ---z-correlator---
      Fname = './LLM/' // CHconf
     &               //'/M' // CHconf //'_'
     &               //CHAR(ICHAR('0')+Ikapp/10)
     &               //CHAR(ICHAR('0')+MOD(Ikapp,10))
     &               //'.'
     &               //CHAR(ICHAR('0')+itime/10)
     &               //CHAR(ICHAR('0')+MOD(itime,10))
      OPEN(UNIT=20,FILE=Fname,STATUS='UNKNOWN',
     &                          FORM='FORMATTED')
        DO 100 iqn= 0, Nqn
             WRITE(20,*) DBLE(  ZWT(iqn,it) ),
     &                   DIMAG( ZWT(iqn,it) )
 100     CONTINUE
        DO 200 iz = 0, Nz-1
        DO 200 iqn  = 0, Nqn
             WRITE(20,*) DBLE(  ZWZ(iqn,iz,it) ),
     &                   DIMAG( ZWZ(iqn,iz,it) )
 200     CONTINUE
      CLOSE(20)
C
C ---Wave function---
      Fname = './LLM/' // CHconf
     &               //'/MW' // CHconf //'_'
     &               //CHAR(ICHAR('0')+Ikapp/10)
     &               //CHAR(ICHAR('0')+MOD(Ikapp,10))
     &               //'.'
     &               //CHAR(ICHAR('0')+itime/10)
     &               //CHAR(ICHAR('0')+MOD(itime,10))
      OPEN(UNIT=21,FILE=Fname,STATUS='UNKNOWN',
     &                          FORM='FORMATTED')
        DO 500 iqn= 0, Nqn
             WRITE(21,*) DBLE(  ZWT(iqn,it) ),
     &                   DIMAG( ZWT(iqn,it) )
 500     CONTINUE
        DO 600 ir = 1, Nz/2
        DO 600 iqn  = 0, Nqn
             WRITE(21,*) DBLE(  ZWW(iqn,ir,it) ),
     &                   DIMAG( ZWW(iqn,ir,it) )
 600     CONTINUE
      CLOSE(21)
C
 2000 CONTINUE
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE SETGM(IGM,ZGM)
C          -----  Definition of Gamma matrices  ------
C
      INTEGER    IGM(4,0:9)
      COMPLEX*16 ZGM(4,0:9)
C
C-------------------------Main------------------------------
C
C----- GM0 = 1 ( UNIT MATRIX ) -----
          IGM(1,0) = 1
          IGM(2,0) = 2
          IGM(3,0) = 3
          IGM(4,0) = 4
       ZGM(1,0) = DCMPLX(  1.D0 ,  0.D0 )
       ZGM(2,0) = DCMPLX(  1.D0 ,  0.D0 )
       ZGM(3,0) = DCMPLX(  1.D0 ,  0.D0 )
       ZGM(4,0) = DCMPLX(  1.D0 ,  0.D0 )
C
C----- GM1 = GAMMA(1) -----
          IGM(1,1) = 4
          IGM(2,1) = 3
          IGM(3,1) = 2
          IGM(4,1) = 1
       ZGM(1,1) = DCMPLX(  0.D0 , -1.D0 )
       ZGM(2,1) = DCMPLX(  0.D0 , -1.D0 )
       ZGM(3,1) = DCMPLX(  0.D0 ,  1.D0 )
       ZGM(4,1) = DCMPLX(  0.D0 ,  1.D0 )
C
C----- GM2 = GAMMA(2) -----
          IGM(1,2) = 4
          IGM(2,2) = 3
          IGM(3,2) = 2
          IGM(4,2) = 1
       ZGM(1,2) = DCMPLX( -1.D0 ,  0.D0 )
       ZGM(2,2) = DCMPLX(  1.D0 ,  0.D0 )
       ZGM(3,2) = DCMPLX(  1.D0 ,  0.D0 )
       ZGM(4,2) = DCMPLX( -1.D0 ,  0.D0 )
C
C----- GM3 = GAMMA(3) -----
          IGM(1,3) = 3
          IGM(2,3) = 4
          IGM(3,3) = 1
          IGM(4,3) = 2
       ZGM(1,3) = DCMPLX(  0.D0 , -1.D0 )
       ZGM(2,3) = DCMPLX(  0.D0 ,  1.D0 )
       ZGM(3,3) = DCMPLX(  0.D0 ,  1.D0 )
       ZGM(4,3) = DCMPLX(  0.D0 , -1.D0 )
C
C----- GM4 = GAMMA(4) -----
          IGM(1,4) = 1
          IGM(2,4) = 2
          IGM(3,4) = 3
          IGM(4,4) = 4
       ZGM(1,4) = DCMPLX(  1.D0 ,  0.D0 )
       ZGM(2,4) = DCMPLX(  1.D0 ,  0.D0 )
       ZGM(3,4) = DCMPLX( -1.D0 ,  0.D0 )
       ZGM(4,4) = DCMPLX( -1.D0 ,  0.D0 )
C
C----- GM5 = GAMMA(5) -----
          IGM(1,5) = 3
          IGM(2,5) = 4
          IGM(3,5) = 1
          IGM(4,5) = 2
       ZGM(1,5) = DCMPLX(  1.D0 ,  0.D0 )
       ZGM(2,5) = DCMPLX(  1.D0 ,  0.D0 )
       ZGM(3,5) = DCMPLX(  1.D0 ,  0.D0 )
       ZGM(4,5) = DCMPLX(  1.D0 ,  0.D0 )
C
C----- GM6 = GAMMA(1)*GAMMA(5) -----
          IGM(1,6) = 2
          IGM(2,6) = 1
          IGM(3,6) = 4
          IGM(4,6) = 3
       ZGM(1,6) = DCMPLX(  0.D0 , -1.D0 )
       ZGM(2,6) = DCMPLX(  0.D0 , -1.D0 )
       ZGM(3,6) = DCMPLX(  0.D0 ,  1.D0 )
       ZGM(4,6) = DCMPLX(  0.D0 ,  1.D0 )
C
C----- GM7 = GAMMA(2)*GAMMA(5) -----
          IGM(1,7) = 2
          IGM(2,7) = 1
          IGM(3,7) = 4
          IGM(4,7) = 3
       ZGM(1,7) = DCMPLX( -1.D0 ,  0.D0 )
       ZGM(2,7) = DCMPLX(  1.D0 ,  0.D0 )
       ZGM(3,7) = DCMPLX(  1.D0 ,  0.D0 )
       ZGM(4,7) = DCMPLX( -1.D0 ,  0.D0 )
C
C----- GM8 = GAMMA(3)*GAMMA(5) -----
          IGM(1,8) = 1
          IGM(2,8) = 2
          IGM(3,8) = 3
          IGM(4,8) = 4
       ZGM(1,8) = DCMPLX(  0.D0 , -1.D0 )
       ZGM(2,8) = DCMPLX(  0.D0 ,  1.D0 )
       ZGM(3,8) = DCMPLX(  0.D0 ,  1.D0 )
       ZGM(4,8) = DCMPLX(  0.D0 , -1.D0 )
C
C----- GM9 = GAMMA(4)*GAMMA(5) -----
          IGM(1,9) = 3
          IGM(2,9) = 4
          IGM(3,9) = 1
          IGM(4,9) = 2
       ZGM(1,9) = DCMPLX(  1.D0 ,  0.D0 )
       ZGM(2,9) = DCMPLX(  1.D0 ,  0.D0 )
       ZGM(3,9) = DCMPLX( -1.D0 ,  0.D0 )
       ZGM(4,9) = DCMPLX( -1.D0 ,  0.D0 )
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
C***************************************************END*****

