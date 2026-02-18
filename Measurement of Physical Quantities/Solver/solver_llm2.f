C               <  solver_llm2.f  >
C***********************************************************
C   SU(3) WILSON FERMION SOLVER (Type-2)
C                                   Ver. 1.2   28 Jan 2000
C  .Anisotropic lattice
C  .Type-2 (monotonous site index)
C  .For single vector processor
C  .Meson correlators are calculated
C    If you don't need them, comment out the subroutines
C        LLMinit, SetZ1, LLMCorr, WRITEW,
C    and compile without 'llmcorr.f'
C
C                      Copyright (c)  Hideo Matsufuru 2000
C***********************************************************
      INCLUDE 'lattsize.h'
      INCLUDE 'solvcomm2.h'
C------------Variables--------------------------------------
      DIMENSION CORR(NT)
C------------Main-------------------------------------------
      CALL SETPRM
C       ( Read Parameters )
      CALL SITEINDEX(LV)
C       ( Set list vector )
      CALL BCINDEX
C       ( Set boundary conditions )
      CALL Uinit
C       ( Set link variables )
      CALL SETBsrc
C       ( Set source function: normal ordered vector )
      CALL LLMinit
C       ( Light meson correlators )
C
C ---Anouncement---
      WRITE(*,*) 'SU(3) Anisotropic Quark Solver (type-2)'
      WRITE(*,13) ' size =',NX,'*',NY,'*',NZ,'*',Ntime
      WRITE(*,13) ' b.c. :', NBC(1),' ',NBC(2),' ',NBC(3),
     &                   ' ',NBC(4)
      WRITE(*,12) ' NPE  =',NPE,   '  Nt    =',NT
      WRITE(*,11) ' Iconf  =',Iconf
      WRITE(*,11) ' Ikapp  =',Ikapp
      WRITE(*,14) ' CKs    =',CKs,'  CKt =',CKs*gammaF,
     &           '  gamma_F =',gammaF
      WRITE(*,12) ' Isrc   =',Isrc, '  NtimeS = ',NtimeS
      WRITE(*,15) ' Niter  =',Niter,'  Enorm =',Enorm
      WRITE(*,11) ' Istart =',Istart 
 11   FORMAT(   A,I4  )
 12   FORMAT( 2(A,I4) )
 13   FORMAT( 4(A,' ',I3) )
 14   FORMAT( 3(A,F7.4) )
 15   FORMAT( A,I6,A,E10.2 )
C
C ---Correlator---
      DO 100 it = 1,NT
         CORR(it) = 0.D0
 100  CONTINUE
C
C ---Solver start---
      DO 3000 ICS=1,Ncol
C
      DO 2000 IDS=1,ND
C
      WRITE(*,*)
      WRITE(*,*) 'IDS =',IDS,'  ICS=',ICS 
      WRITE(*,*)
C
      CALL SETSRC(IDS,ICS)
C       ( Set source )
      WRITE(*,*) 'Source set o.k.'
C
      CALL SOLVE
C       ( Main solver routine )
      WRITE(*,*) 'Iteration finished.'
C
      CALL CheckX( IDS, ICS, Rdiff )
C       ( Calculate difference: b-Dx )
C
       WRITE(*,*) 'Diff. =', Rdiff
      CALL FIN(IDS,ICS)
C       ( Write propagator on disk )
C
      Nvs = NS * Nvc * ND
      DO 2200 it = 1, NT
        CALL NormV(Nvs,CORRF,X(1,1,1,it))
        CORR(it) = CORR(it) + CORRF
 2200  CONTINUE
C
      CALL SetZ1(IDS,X)
C       ( Complex vector with Dirac index at the source )
C
 2000 CONTINUE
C      --- ( Here souce-spinor(IDS) do loop end )
C
      CALL LLMCORR
C       ( Light meson correlators )
C
 3000 CONTINUE
C      --- ( Here souce-color(ICS) do loop end )
C
      CALL WRITEW(Iconf, Ikapp)
C       ( Write meson corr. to disk )
C
      WRITE(*,*)
      DO 4000 itime = 0, Ntime-1
        WRITE(*,*) itime, CORR(MOD(itime,NT)+1)  
 4000 CONTINUE
C
C-----------------------------------------------------------
      STOP
      END
C***********************************************************
      SUBROUTINE SETPRM
C
      INCLUDE 'lattsize.h'
      INCLUDE 'solvcomm2.h'
C------------Notation---------------------------------------
C    Iconf   : Conf. number
C    Ikapp   : Kappa number
C    CKs     : Kappa in spacial direction
C    gammaF  : Anisotropic param. for quark
C    NBC(1-4): Boundary conditions
C             (1: periodic, -1: anti-periodic, 0: Diriclet)
C    Isrc    : Label of source (0:point, etc...)
C    NtimeS  : Time slice on which source is put [0,Ntime-1]
C    Niter   : Maximum iteration number
C    Istart  : ( 0: I(free) / no O,   1: I(free)/ O, 
C                2: I(read) / no O,   3: I(read)/ O  )
C    Enorm   : Convergence criterion
C------------Variables--------------------------------------
      CHARACTER*20 FNPRM
C------------Main-------------------------------------------
C ---Read from standard input---
      READ(*,*)  Iconf
      READ(*,*)  Ikapp
C
C ---Read from parameter file---
      FNPRM = 'quarkprm.'
     &          //CHAR( ICHAR('0')+MOD(Ikapp/10,10) )
     &          //CHAR( ICHAR('0')+MOD(Ikapp,10) )
      OPEN(18,FILE=FNPRM,STATUS='OLD',FORM='FORMATTED')
       READ(18,*) IkappD
       READ(18,*) CKs,gammaF
       READ(18,*) NBC(1),NBC(2),NBC(3),NBC(4)
       READ(18,*) Isrc, NtimeS
       READ(18,*) Niter,Istart
       READ(18,*) Enorm
      CLOSE(18)
      IF(IkappD.NE.Ikapp) WRITE(*,*) 'Ikappa is wrong !'
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE Uinit
C
      INCLUDE 'lattsize.h'
      INCLUDE 'solvcomm2.h'
C------------Variables--------------------------------------
      CHARACTER*120 RFILE
      CHARACTER*4 CHconf
c      REAL*4 UR(NS,Ndf,4)
      DIMENSION UR(NS,Ndf,4)
C------------Main-------------------------------------------
C
      IF( Istart .LE. 1 ) THEN
C                                     ( Free field )
C
       DO 500 it = 1, NT
       DO 500 idir = 1,4
         CALL SetUnitv( NS,U(1,1,idir,it) )
 500   CONTINUE
C
      ELSE 
C           --- Read data from file ---
C
      CHconf = '0' //CHAR( ICHAR('0')+    Iconf/100 )
     &             //CHAR( ICHAR('0')+MOD(Iconf/10,10) )
     &             //CHAR( ICHAR('0')+MOD(Iconf,10) )
C
      DO 1000 it = 1, NT
        itime = it - 1
C
      RFILE = './ConfU/' // CHconf
     &          //'/G' // CHconf //'.'
     &          //CHAR( ICHAR('0')+    itime/10  )
     &          //CHAR( ICHAR('0')+MOD(itime,10) )
      OPEN(UNIT=12, FILE=RFILE,
     &                 STATUS='OLD', FORM='UNFORMATTED')
       READ(12) UR
      CLOSE(12)
C
      DO 1700 idir = 1, 4
      DO 1700 idf  = 1, Ndf
      DO 1700 is   = 1, NS
        U(is,idf,idir,it) = DBLE( UR(is,idf,idir) )
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
      SUBROUTINE FIN( IDS, ICS )
C
      INCLUDE 'lattsize.h'
      INCLUDE 'solvcomm2.h'
C------------Variables--------------------------------------
      CHARACTER*120 WFILE
      CHARACTER*4 CHconf
c      REAL*4 XW(NS,Nvc,ND)
      DIMENSION XW(NS,Nvc,ND)
C------------Main-------------------------------------------
C
      IF( MOD(ISTART,2) .EQ. 1 ) THEN
C
      CHconf = '0' //CHAR( ICHAR('0')+    Iconf/100 )
     &             //CHAR( ICHAR('0')+MOD(Iconf/10,10) )
     &             //CHAR( ICHAR('0')+MOD(Iconf,10) )
C
      DO 1000 itime = 0,Ntime-1
C
       DO 1200 id  = 1, ND
       DO 1200 ivc = 1, Nvc
       DO 1200 is  = 1, NS
         XW(is,ivc,id) = X(is,ivc,id,itime+1)
 1200  CONTINUE
C
      WFILE = './ConfX/' // CHconf
     &           //'/G' // CHconf //'_'
     &           //CHAR( ICHAR('0')+MOD(Ikapp/10,10) )
     &           //CHAR( ICHAR('0')+MOD(Ikapp,10) )
     &           //CHAR( ICHAR('0')+MOD(IDS,10) )
     &           //CHAR( ICHAR('0')+MOD(ICS,10) )
     &           //'.'
     &           //CHAR(ICHAR('0')+ itime/10)
     &           //CHAR(ICHAR('0')+MOD(itime,10))
       OPEN(UNIT=13, FILE=WFILE,
     &           STATUS='UNKNOWN', FORM='UNFORMATTED')
       WRITE(13) XW
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

