C***********************************************************
C   Staple Construction Program - (e,o) version
C
C         .Vectorized version      14 Jul 1999 H.Matsufuru
C***********************************************************
      SUBROUTINE STPL(Ioe,Joe,MU,IT)
C
      INCLUDE 'lattsize.h'
      INCLUDE 'updcomm.h'
C------------Local vaiables---------------------------------
      DIMENSION U1(NS2,Ndf),U2(NS2,Ndf),U3(NS2,Ndf)
C------------Main-------------------------------------------
      Czero = 0.D0
      GammaR = 1.D0/Gamma
      Nvs2 = NS2 * Ndf
C
         CALL SetConst(Nvs2,C(1,1),Czero)
C
      ITP = MOD(IT,NT)+1
      ITM = MOD(IT-2+NT,NT)+1
C
      IF( MU .EQ. 4 ) THEN
C                           ( Mu = 4 case )
      DO 1000 nu=1,3
C
C --- Upper mu-nu staples ---         
C
      CALL Sortv( NS2,Ndf,
     &            U1(1,1),U(1,1,Joe,MU,IT),LV(1,nu,1,Ioe) )
      CALL UPrdvND( NS2, U2,U1,U(1,1,Ioe,nu,ITP) )
      CALL UPrdvNN( NS2, U3,U(1,1,Ioe,nu,IT),U2 )
      CALL SelfAdd(Nvs2,C,U3,Gamma)
C
C --- Lower mu-nu staples ---         
      CALL UPrdvNN( NS2, U1,
     &                U(1,1,Joe,MU,IT),U(1,1,Joe,nu,ITP) )
      CALL UPrdvDN( NS2, U2,U(1,1,Joe,nu,IT),U1 )
      CALL Sortv( NS2,Ndf, U3,U2,LV(1,nu,2,Ioe) )
      CALL SelfAdd(Nvs2,C,U3,Gamma)
C
 1000 CONTINUE
C
      ELSE
C          ( Mu = i case )
C
C----- Nu = 4 case -----
C
C --- Upper mu-nu staples ---
      CALL Sortv( NS2,Ndf,
     &                U1,U(1,1,Joe,4,IT),LV(1,MU,1,Ioe) )
      CALL UPrdvND( NS2, U2,U(1,1,Ioe,MU,ITP),U1 )
      CALL UPrdvNN( NS2, U3,U(1,1,Ioe,4,IT),U2 )
      CALL SelfAdd(Nvs2,C,U3,Gamma)
C
C --- Lower mu-nu staples ---         
      CALL Sortv( NS2,Ndf,
     &                U1,U(1,1,Joe,4,ITM),LV(1,MU,1,Ioe) )
      CALL UPrdvNN( NS2, U2,U(1,1,Ioe,MU,ITM),U1 )
      CALL UPrdvDN( NS2, U3,U(1,1,Ioe,4,ITM),U2 )
      CALL SelfAdd(Nvs2,C,U3,Gamma)
C
C----- Nu = i case -----
      DO 3000 nu = 1, 3
      IF( nu .EQ. MU ) GOTO 3000
C
C --- Upper mu-nu staples ---         
      CALL Sortv( NS2,Ndf,
     &                U1,U(1,1,Joe,MU,IT),LV(1,nu,1,Ioe) )
      CALL Sortv( NS2,Ndf,
     &                U2,U(1,1,Joe,nu,IT),LV(1,MU,1,Ioe) )
      CALL UPrdvND( NS2, U3,U1,U2 )
      CALL UPrdvNN( NS2, U1,U(1,1,Ioe,nu,IT),U3 )
      CALL SelfAdd(Nvs2,C,U1,GammaR)
C
C --- Lower mu-nu staples ---         
      CALL Sortv( NS2,Ndf,
     &                U1,U(1,1,Ioe,nu,IT),LV(1,MU,1,Joe) )
      CALL UPrdvNN( NS2, U2,U(1,1,Joe,MU,IT),U1 )
      CALL UPrdvDN( NS2, U3,U(1,1,Joe,nu,IT),U2 )
      CALL Sortv( NS2,Ndf, U1,U3,LV(1,nu,2,Ioe) )
      CALL SelfAdd(Nvs2,C,U1,GammaR)
C         
 3000 CONTINUE
C
      ENDIF
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
C***************************************************END*****
