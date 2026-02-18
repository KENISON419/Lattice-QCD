C                     <  opr.f  >  
C***********************************************************
C  Operator for Anisotropic Wilson quark (Type-2)
C
C    For Intel Paragon XP/S     14 Sep 1995  Matsufuru, H.
C    Modified:                  13 Feb 1998, 28 Aug 1998
C    For WS (single PE)         25 May 1999  Matsufuru, H.
C    Modified, for Vector processor  15 Jun 1999 M.,H. 
C***********************************************************
      SUBROUTINE Dopr(V,W,PRF,JD,Iflg)
C             ----- variation of K[U] Matrix Operation -----
C
C       .Iflg=0: V =     PRF *[ 1 - kM ] W
C        Iflg=1: V = V + PRF *[ 1 - kM ] W
C       .PRF is a real number 
C       .JD=1 : K , JD=-1 : K^dagger
C-----------------------------------------------------------
      INCLUDE 'lattsize.h'
      INCLUDE 'solvcomm2.h'
C------------Variables--------------------------------------
      DIMENSION  V(NS,Nvc,ND,NT), W(NS,Nvc,ND,NT)
      DIMENSION W2(NS,Nvc,ND,NT) 
C------------Main-------------------------------------------
      Nvst = NS * Nvc * ND * NT
C
      CALL Mopr( W2, W, JD )
C
      IF(Iflg.EQ.0) THEN
        CALL SetVec( Nvst,V,W , PRF)
        CALL SelfAdd(Nvst,V,W2,-PRF)
      ELSE
        CALL SelfAdd(Nvst,V,W , PRF)
        CALL SelfAdd(Nvst,V,W2,-PRF)
      ENDIF
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE Mopr(V,W,JD)
C
      INCLUDE 'lattsize.h'
      INCLUDE 'solvcomm2.h'
C------------Variables--------------------------------------
      DIMENSION V(NS,Nvc,ND,NT),W(NS,Nvc,ND,NT) 
      DIMENSION UX1(NS,Nvc,ND),UX2(NS,Nvc,ND)
C------------Main-------------------------------------------
      RW = 1.D0
      Cunit = 1.D0
      CKt = CKs * gammaF
      FJD = DBLE(JD)
C
      Ncd = Nvc * ND
      Nvs = NS * Nvc * ND
C
      DO 1000 IT = 1,NT
C
C ---< mu = + 1 >---
      CALL Sortv(NS,Ncd,UX1,W(1,1,1,IT),LV(1,1,1))
      CALL UxPrdv_N(NS,UX2,U(1,1,1,IT),UX1)
      CALL AddG1xPrd(NS,UX1,UX2,RW,-FJD)
      CALL SetBFV(NS,Ncd,V(1,1,1,IT),UX1,BF(1,1,1),CKs)
C
C ---< mu = - 1 >---
      CALL UxPrdv_D(NS,UX1,U(1,1,1,IT),W(1,1,1,IT))
      CALL Sortv(NS,Ncd,UX2,UX1,LV(1,1,2))
      CALL AddG1xPrd(NS,UX1,UX2,RW, FJD)
      CALL AddBFV(NS,Ncd,V(1,1,1,IT),UX1,BF(1,1,2),CKs)
C
C ---< mu = + 2 >---
      CALL Sortv(NS,Ncd,UX1,W(1,1,1,IT),LV(1,2,1))
      CALL UxPrdv_N(NS,UX2,U(1,1,2,IT),UX1)
      CALL AddG2xPrd(NS,UX1,UX2,RW,-FJD)
      CALL AddBFV(NS,Ncd,V(1,1,1,IT),UX1,BF(1,2,1),CKs)
C
C ---< mu = - 2 >---
      CALL UxPrdv_D(NS,UX1,U(1,1,2,IT),W(1,1,1,IT))
      CALL Sortv(NS,Ncd,UX2,UX1,LV(1,2,2))
      CALL AddG2xPrd(NS,UX1,UX2,RW, FJD)
      CALL AddBFV(NS,Ncd,V(1,1,1,IT),UX1,BF(1,2,2),CKs)
C
C ---< mu = + 3 >---
      CALL Sortv(NS,Ncd,UX1,W(1,1,1,IT),LV(1,3,1))
      CALL UxPrdv_N(NS,UX2,U(1,1,3,IT),UX1)
      CALL AddG3xPrd(NS,UX1,UX2,RW,-FJD)
      CALL AddBFV(NS,Ncd,V(1,1,1,IT),UX1,BF(1,3,1),CKs)
C
C ---< mu = - 3 >---
      CALL UxPrdv_D(NS,UX1,U(1,1,3,IT),W(1,1,1,IT))
      CALL Sortv(NS,Ncd,UX2,UX1,LV(1,3,2))
      CALL AddG3xPrd(NS,UX1,UX2,RW, FJD)
      CALL AddBFV(NS,Ncd,V(1,1,1,IT),UX1,BF(1,3,2),CKs)
C
C ---< mu = + 4 >---
        ITP = MOD(IT,NT)+1
        BF4 = 1.D0
        IF( IT .EQ. Ntime ) BF4 = DBLE( NBC(4) )
        CKtBF4 = CKt * BF4
      CALL UxPrdv_N(NS,UX2,U(1,1,4,IT),W(1,1,1,ITP))
      CALL AddG4xPrd(NS,UX1,UX2,Cunit,-FJD)
      CALL SelfAdd(Nvs,V(1,1,1,IT),UX1,CKtBF4)
C
C ---< mu = - 4 >---
        ITM = MOD(IT-2+NT,NT) + 1
        BF4 = 1.D0
        IF( IT .EQ. 1 ) BF4 = DBLE( NBC(4) )
        CKtBF4 = CKt * BF4
      CALL UxPrdv_D(NS,UX2,U(1,1,4,ITM),W(1,1,1,ITM))
      CALL AddG4xPrd(NS,UX1,UX2,Cunit, FJD)
      CALL SelfAdd(Nvs,V(1,1,1,IT),UX1,CKtBF4)
C
 1000 CONTINUE
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE SETBsrc
C
      INCLUDE 'lattsize.h'
      INCLUDE 'solvcomm2.h'
C----------------Variables----------------------------------
      CHARACTER*20 FNsrc
C------------Main-------------------------------------------
C
      Czero = 0.D0
C
      NPEsrc = NtimeS/NT
      NTsrc  = MOD(NtimeS,NT)+1
C
C ---Source set on Normal Ordered Vector---
      IF( Isrc.EQ.0 ) THEN
C                          ------------ Point source ---
        NSsrc  = 1
C
        Nv = NS
        CALL SETCONST( Nv, Bsrc, Czero )
        Bsrc( NSsrc ) = 1.D0
      END IF
C
      IF( Isrc.GT.0 ) THEN
C                        ------------ Smeared source ---
        FNsrc = 'srcfunc.'
     &           //CHAR( ICHAR('0')+MOD(Isrc/10,10) )
     &           //CHAR( ICHAR('0')+MOD(Isrc,10) )
        OPEN(18,FILE=FNsrc,STATUS='OLD',FORM='FORMATTED')
          DO 1000 is=1,NS
             READ(18,*) isD, Bsrcf
             Bsrc(is) = DBLE( Bsrcf )
 1000     CONTINUE
        CLOSE(18)
C
      ENDIF
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE SETSRC(IDS,ICS)
C
      INCLUDE 'lattsize.h'
      INCLUDE 'solvcomm2.h'
C------------Main-------------------------------------------
C
      Czero = 0.D0
      Cunit = 1.D0
      Nvst = Nvc *ND *NS *NT
C
      NTsrc  = MOD(NtimeS,NT)+1
C
      CALL SetConst(Nvst,S,Czero)
C
      IVCS   = ICS*2-1
C
      CALL Equate(NS,S(1,IVCS,IDS,NTsrc),Bsrc(1))
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE CheckX(IDS,ICS,Rdiff)
C
      INCLUDE 'lattsize.h'
      INCLUDE 'solvcomm2.h'
C------------Main-------------------------------------------
      Cunit = 1.D0
      Nvst = Nvc *ND *NS *NT
      Snorm = 1.0
c      CALL NormV( Nvst,SR, S )
c      Snorm = 1.0 / SR
cc      Snorm = 1.0 / DBLE(NS*Ntime)
C
      CALL SETSRC(IDS,ICS)
      CALL Dopr(S,X,-Cunit,1,1)
      CALL NormV( Nvst,Rdiff,S )
      Rdiff = Rdiff * Snorm
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE  BCINDEX
C   --- Set the spacial boundary conditions ---
C
      INCLUDE 'lattsize.h'
      INCLUDE 'solvcomm2.h'
C-----------------Main--------------------------------------
C
      DO 1000  is = 1,NS
        BF(is,1,1) = 1.D0
        BF(is,1,2) = 1.D0
        BF(is,2,1) = 1.D0
        BF(is,2,2) = 1.D0
        BF(is,3,1) = 1.D0
        BF(is,3,2) = 1.D0
 1000 CONTINUE
C
      DO 1400  iz = 0, NZ-1
      DO 1400  iy = 0, NY-1
      DO 1400  ix = 0, NX-1
C
          is = 1 + ix + iy * NX + iz * NX*NY
C
       IF( ix .EQ. NX-1 )  BF( is,1,1 ) = DBLE( NBC(1) )
       IF( ix .EQ.    0 )  BF( is,1,2 ) = DBLE( NBC(1) )
C
       IF( iy .EQ. NY-1 )  BF( is,2,1 ) = DBLE( NBC(2) )
       IF( iy .EQ.    0 )  BF( is,2,2 ) = DBLE( NBC(2) )
C
       IF( iz .EQ. NZ-1 )  BF( is,3,1 ) = DBLE( NBC(3) )
       IF( iz .EQ.    0 )  BF( is,3,2 ) = DBLE( NBC(3) )
C
 1400 CONTINUE
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
C***************************************************END*****


