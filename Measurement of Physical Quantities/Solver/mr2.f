C***********************************************************
C   Minimal Residual Solver with ILU preconditioner
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
      Nvst = Nvc *ND *NS *NT
      CALL NormV( Nvst,SR, S )
C
      Snorm = 1.D0 / SR
C      Snorm = 1.D0 / DBLE(NS*Ntime)
C
C ---Initial step---
C
      CALL MRinit( RN )
      Rnorm = RN * Snorm
      WRITE(*,*) 'init:',Rnorm
C ---Iteration---
      DO 1000 iter = 1, Niter
C
         CALL MRstep( RN )
         Rnorm = RN * Snorm
C
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
      SUBROUTINE MRinit( RR )
C
      INCLUDE 'lattsize.h'
      INCLUDE 'solvcomm2.h'
C------------Main-------------------------------------------
      Cunit = 1.D0
      Nvst = Nvc *ND *NS *NT
C
      CALL Equate( Nvst,R, S )
      CALL Equate( Nvst,X, S )
C
      CALL Dopr( R, S, -Cunit, 1, 1 )
C
      CALL NormV( Nvst,RR, R )
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE MRstep( RRN )
C
      INCLUDE 'lattsize.h'
      INCLUDE 'solvcomm2.h'
C------------Main-------------------------------------------
      Cunit = 1.D0
      Nvst = Nvc *ND *NS *NT
C
      CALL Dopr( S, R, Cunit, 1, 0 )
C
      CALL NormV( Nvst,ADEN, S )
      CALL VecPrd( Nvst,ANUM, R, S )
      ALPH = ANUM/ADEN
      CALL SelfAdd( Nvst,X, R,  ALPH ) 
      CALL SelfAdd( Nvst,R, S, -ALPH ) 
C
      CALL NormV( Nvst,RRN, R )
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
C***************************************************END*****
