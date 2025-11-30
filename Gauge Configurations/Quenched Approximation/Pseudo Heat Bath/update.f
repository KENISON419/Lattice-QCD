C               <  update.f  >
C***********************************************************
C   Quenched SU(3) Anisotropic Lattice Update Program
C                                   Ver. 1.0   14 Jul 1999
C
C    .Anisotropic lattice
C    .(e,o) site index version
C    .For single vector processor
C
C  Notice:
C   Pseudoheat bath update routine (in phbupd.f) is written
C     by S.Hioki.  c.f. QCDMPI
C             (http://kurobe.sci.hiroshima-u.ac.jp/QCDMPI/)
C   Random number generator (in lib_ran.f) is written
C     by J.Makino. c.f. Parallel Computing 20 (1994) 1357.
C
C                      Copyright (c)  Hideo Matsufuru 1999
C***********************************************************
      PROGRAM UPDATE
C
      INCLUDE 'lattsize.h'
      INCLUDE 'updcomm.h'
C------------Main-------------------------------------------
C
      Rnorm=1.0/DBLE( NS *Ntime *6 *3 )
C
      CALL SETPRM
C       ( Read Parameters )
      CALL SITEINDEX2(LV)
C       ( Set list vector: (e,o) version )
      CALL Uinit
C       ( Set link variables )
      Ndelay = Iconf + 1
c      Ndelay = Iconf - 1
      CALL CINIT3(Ndelay)
C       ( Random number initialization )
C
C ---Anouncement---
      WRITE(*,*) '--- SU(3) Anisotropic lattice update ---'
      WRITE(*,*) '              (vector processor version)'
      WRITE(*,13) ' Size =',Nx,' *',Ny,' *',Nz,' *',Ntime
      WRITE(*,12) ' Nt    =',NT
      WRITE(*,11) ' Iconf  =',Iconf
      WRITE(*,14) ' beta   =',Beta,'  gamma =',Gamma
      WRITE(*,15) ' Nprev  =',Nprev
      WRITE(*,15) ' Nsweep =',Nsweep
      WRITE(*,11) ' Nhit   =',Nhit
      WRITE(*,12) ' Istart =',Istart,'  Iplaq =',Iplaq
      WRITE(*,*)
      WRITE(*,*) '  sweep    Uplaq(s)    Uplaq(t)      time'
      WRITE(*,*)
     &        ' -------------------------------------------'
 11   FORMAT( A,I4  )
 12   FORMAT( 2(A,I4) )
 13   FORMAT( 4(A,I4) )
 14   FORMAT( 2(A,F7.4) )
 15   FORMAT( A,I10 )
C
      DO 1000 isweep=1,Nsweep
C
       Plaqs = 0.D0
       Plaqt = 0.D0
c       Etime1 = DCLOCK()
C
       DO 1200  MU = 1, 4
C         ( MU : direction --> 1:x, 2:y, 3:z, 4:t )
C 
       DO 1200 Iturn = 1,2
C         ( Iturn : (odd,even) turn )
C
       DO 1200  IT = 1, NT
C         ( IT : time coordinate )
C
        Ioe = MOD( Iturn +IT , 2 ) + 1
        Joe = 3-Ioe
C         ( (Ioe,Joe) = (1,2) : odd,  (2,1) : even )
C
       CALL STPL(Ioe,Joe,MU,IT)
C                ( Staple construction )
C
       IF(Iplaq.EQ.0) CALL PHBUPD( Ioe,MU,IT )
C                ( Pseudo-heat bath update )
C
       CALL TracevDN( NS2,Plaq,U(1,1,Ioe,MU,IT),C(1,1) )
       IF(MU.EQ.4) Plaqt = Plaqt + Plaq
       IF(MU.LT.4) Plaqs = Plaqs + Plaq
C                ( Plaquette calculation )
C
 1200  CONTINUE
C
       Plaqs = (Plaqs-Plaqt) *0.5D0 * Rnorm *Gamma
       Plaqt = Plaqt * Rnorm /Gamma
c       Etime2=DCLOCK()
C
c       WRITE(6,1201) isweep, Plaqs, Plaqt, Etime2-Etime1
       WRITE(6,1201) isweep, Plaqs, Plaqt
 1201  FORMAT(I8,3F12.6)
C
C ---Reunitarization---
       IF( MOD(ISWEEP,10) .EQ. 0 ) THEN
c         Etime1 = DCLOCK()
         Nvt = 2 * 4 * NT
         CALL REUNIT2v(Nvt,U(1,1,1,1,1))
c         Etime2 = DCLOCK()
c          WRITE(*,1001) '    Reunitalized',Etime2-Etime1
       ENDIF
 1001  FORMAT(A,16X,F12.6 )
C
 1000 CONTINUE
C
      CALL FIN
C
C-----------------------------------------------------------
      STOP
      END
C***********************************************************
      SUBROUTINE SETPRM
C
      INCLUDE 'lattsize.h'
      INCLUDE 'updcomm.h'
C------------Notation---------------------------------------
C    Ilatt  :  Lattice identification
C    Iconf  :  Configuration number to be generated
C    Istart :  (  0: cold / no O,   1: cold / O, 
C                 2: read / no O,   3: read / O  )
C    Beta   :  \beta
C    Gamma  :  \gamma
C    Ntherm :  Thermalization
C    Nsep   :  Separation between confs
C    Iplaq  :  (  0: update,  1: calculate plaquette  )
C------------Variables--------------------------------------
      CHARACTER*40 FNprm
C------------Main-------------------------------------------
C
C ---Read from standard input---
      READ(*,*)  Ilatt
      READ(*,*)  Iconf
      READ(*,*)  Istart
C
C ---Read from parameter file---
      FNprm = './gaugeprm.'
     &             //CHAR( ICHAR('0')+MOD(Ilatt/10,10) )
     &             //CHAR( ICHAR('0')+MOD(Ilatt,10) )
C
      OPEN(18,FILE=FNprm,STATUS='OLD',FORM='FORMATTED')
        READ(18,*) Beta, Gamma
        READ(18,*) Ntherm, Nsep
        READ(18,*) Nhit
        READ(18,*) Iplaq
      CLOSE(18)
C
      IF(Iconf.EQ.1) THEN
        Nsweep = Ntherm
        Nprev = 0
      ELSE
        Nsweep = Nsep
        Nprev = Ntherm + (Iconf-2)*Nsep
      ENDIF
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE Uinit
C
      INCLUDE 'lattsize.h'
      INCLUDE 'updcomm.h'
C------------Variables--------------------------------------
      CHARACTER*100 RFILE
      CHARACTER*4 CHconf
c      REAL*4 UR(NS,Ndf,4)
      DIMENSION UR(NS,Ndf,4)
C------------Main-------------------------------------------
C
      IF( Istart.LE.1 ) THEN
C           ( Istart 0,1: cold start ; 2,3: read file )
C
C ---Cold start---
      DO 500 it  = 1, NT
      DO 500 mu  = 1, 4
      DO 500 ioe = 1, 2
        CALL SetUnitv( NS2,U(1,1,ioe,mu,it) )
 500  CONTINUE
C
      ELSE 
C           ( Read configuration from file )
C
      Iprev = Iconf - 1
      CHconf = '0' //CHAR( ICHAR('0')+    Iprev/100 )
     &             //CHAR( ICHAR('0')+MOD(Iprev/10,10) )
     &             //CHAR( ICHAR('0')+MOD(Iprev,10) )
C
      DO 1000 it = 1,Ntime
       itime = it - 1
C
      RFILE = './ConfU/' // CHconf //'/G' // CHconf //'.'
     &          //CHAR( ICHAR('0')+    itime/10  )
     &          //CHAR( ICHAR('0')+MOD(itime,10) )
      OPEN(UNIT=12, FILE=RFILE,
     &                 STATUS='OLD', FORM='UNFORMATTED')
        READ(12) UR
      CLOSE(12)
C
C ---(e,o) <--> normal convert--- 
      DO 1700 mu  = 1, 4
      DO 1700 iz  = 0, Nz-1
      DO 1700 iy  = 0, Ny-1
      DO 1700 ix2 = 0, Nx2-1
        Ioe = MOD(iz+iy,2)+1
        Joe = 3-Ioe
          is2 = 1 + ix2 + iy *Nx2 + iz *Nx2*Ny
          isA = 2*is2 -1
          isB = 2*is2
        DO 1800 idf = 1, Ndf
          U( is2,idf,Ioe,mu,it) = DBLE( UR(isA,idf,mu) )
          U( is2,idf,Joe,mu,it) = DBLE( UR(isB,idf,mu) )
 1800   CONTINUE
 1700 CONTINUE
C
 1000 CONTINUE
C
      ENDIF
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE FIN
C
      INCLUDE 'lattsize.h'
      INCLUDE 'updcomm.h'
C------------Variables--------------------------------------
      CHARACTER*100 WFILE
      CHARACTER*4 CHconf
c      REAL*4 UW(NS,Ndf,4)
      DIMENSION UW(NS,Ndf,4)
C------------Main-------------------------------------------
C
      IF( MOD(Istart,2) .EQ. 1 ) THEN
C         ( Istart =0,2: No output ; 1,3: output )
C
      CHconf = '0' //CHAR( ICHAR('0')+    Iconf/100 )
     &             //CHAR( ICHAR('0')+MOD(Iconf/10,10) )
     &             //CHAR( ICHAR('0')+MOD(Iconf,10) )
C
      DO 1000 it = 1,Ntime
       itime = it - 1
C
C --- Convert to normal order ---
      DO 1700 mu  = 1, 4
      DO 1700 iz  = 0, Nz-1
      DO 1700 iy  = 0, Ny-1
      DO 1700 ix2 = 0, Nx2-1
        Ioe = MOD(iz+iy,2)+1
        Joe = 3-Ioe
          is2 = 1 + ix2 + iy *Nx2 + iz *Nx2*Ny
          isA = 2*is2 -1
          isB = 2*is2
        DO 1800 idf=1,Ndf
          UW(isA,idf,mu) = U(is2,idf,Ioe,mu,it)
          UW(isB,idf,mu) = U(is2,idf,Joe,mu,it)
 1800   CONTINUE
 1700 CONTINUE
C
      WFILE = './ConfU/' // CHconf //'/G' // CHconf //'.'
     &          //CHAR( ICHAR('0')+    itime/10  )
     &          //CHAR( ICHAR('0')+MOD(itime,10) )
      OPEN(UNIT=13, FILE=WFILE,
     &                 STATUS='NEW', FORM='UNFORMATTED')
        WRITE(13) UW
      CLOSE(13)
C
 1000 CONTINUE
C
      ENDIF
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
C***************************************************END*****

