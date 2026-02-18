C***********************************************************
C   BCGStab Solver with ILU preconditioner
C    .Notation 
C         X : Solution
C         R : Residual vector
C         S : Working vector
C
C                                20 Dec 1995  H.Matsufuru
C***********************************************************
      SUBROUTINE SOLVE
C
      INCLUDE 'lattsize.h'
      INCLUDE 'solvcomm2.h'
C------------Main-------------------------------------------
C ---Preperarion---
C
      Nvst = NS *Nvc *ND *NT
      CALL NormV( Nvst,SR, S )
C
      Snorm = 1.D0 / SR
C      Snorm = 1.D0 / DBLE(NS*Ntime)
C
C ---Initial step---
      CALL BSinit(RN)
      Rnorm = RN * Snorm
      WRITE(*,*) 'init:',Rnorm
C ---Iteration---
      DO 1000 iter=1,Niter
C
         CALL BSstep(RN)
         Rnorm = RN*Snorm
         WRITE(*,*) iter,Rnorm
         IF( Rnorm.LT.Enorm ) GOTO 1001
C
 1000  CONTINUE
 1001  CONTINUE
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE BSinit(RR)
C
      INCLUDE 'lattsize.h'
      INCLUDE 'solvcomm2.h'
C------------Variables--------------------------------------
      COMMON /FVECT4/ P( NS,Nvc,ND,NT) 
      COMMON /FVECT5/ RH(NS,Nvc,ND,NT) 
      COMMON /FVECT6/ V( NS,Nvc,ND,NT) 
      COMMON /FVECT7/ T( NS,Nvc,ND,NT) 
      COMMON /CONSTS/ RHOP,ALPP,OMGP
C------------Main-------------------------------------------
      Cunit = 1.D0
      Czero = 0.D0
      Nvst = NS *Nvc *ND *NT
C
      CALL EQUATE( Nvst,R,S )
      CALL EQUATE( Nvst,X,S )
C
      CALL Dopr( R, S, -Cunit, 1, 1 )
C
      CALL EQUATE( Nvst,RH,R )
      CALL NORMV( Nvst,RR,R )
        RHOP = 1.D0
        ALPP = 1.D0
        OMGP = 1.D0
C
      CALL SETCONST(Nvst,P,Czero)
      CALL SETCONST(Nvst,V,Czero)
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE BSstep(RRN)
C
      INCLUDE 'lattsize.h'
      INCLUDE 'solvcomm2.h'
C------------Variables--------------------------------------
      COMMON /FVECT4/ P( NS,Nvc,ND,NT) 
      COMMON /FVECT5/ RH(NS,Nvc,ND,NT) 
      COMMON /FVECT6/ V( NS,Nvc,ND,NT) 
      COMMON /FVECT7/ T( NS,Nvc,ND,NT) 
      COMMON /CONSTS/ RHOP,ALPP,OMGP
C------------Main-------------------------------------------
      Cunit = 1.D0
      Nvst = NS *Nvc *ND *NT
C
      CALL VecPrd(Nvst,RHO,RH,R)
C
      BET = RHO*ALPP/(RHOP*OMGP)
C
      CALL SelfSPrd(Nvst,P,BET)
      CALL SelfAdd( Nvst,P,R,Cunit )
      CALL SelfAdd( Nvst,P,V,-OMGP*BET )
C
      CALL Dopr( V, P, Cunit, 1, 0 )
C
      CALL VecPrd(Nvst,ADEN,RH,V)
      ALP = RHO/ADEN
C
      CALL Equate( Nvst,S,R )
      CALL SelfAdd( Nvst,S,V,-ALP )
C
      CALL Dopr( T, S, Cunit, 1, 0 )
C
      CALL VecPrd( Nvst,OMGN,T,S)
      CALL NormV( Nvst,OMGD,T)
      OMG=OMGN/OMGD
C
      CALL SelfAdd( Nvst,X,S,OMG )
      CALL SelfAdd( Nvst,X,P,ALP )
C
      CALL Equate( Nvst,R,S )
      CALL SelfAdd( Nvst,R,T,-OMG )
C
      RHOP = RHO
      ALPP = ALP
      OMGP = OMG
C
      CALL NormV( Nvst,RRN,R)
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
C***************************************************END*****
