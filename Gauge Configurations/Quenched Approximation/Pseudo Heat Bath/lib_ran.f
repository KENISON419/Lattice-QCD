      SUBROUTINE RNDS(RN,N,M)
C
C-----------------------------------------------------------
      PARAMETER( IP= 521, IQ= 32 )
      PARAMETER( MACRM= 40, MACRI= 1 )
      INTEGER P,Q
      PARAMETER(P=IP,Q=IQ)
      PARAMETER(Fnorm =0.465661E-9, Fadd = 1.E-12)
      INTEGER IR(0:P-1)
      COMMON /RANDM/ IR,J,K
      DIMENSION RN(N)
C
C------------Main-------------------------------------------
C
      DO 100 i=1,m 
        IR(J)=XOR(IR(J),IR(K))
C        RN(i)=IR(J)*Fnorm
        RN(i)=IR(J)*Fnorm + Fadd
c        IF (RN(i).eq.1.0) STOP
        J=J+1
        IF(J.GE.P)J=J-P
        K=K+1
        IF(K.GE.P)K=K-P
  100 CONTINUE
C    
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      FUNCTION RANF()
C-----------------------------------------------------------
      PARAMETER( IP= 521, IQ= 32 )
      PARAMETER( MACRM= 40, MACRI= 1 )
      INTEGER P,Q
      PARAMETER(P=IP,Q=IQ)
      PARAMETER(Fnorm =0.465661E-9, Fadd = 1.E-12)
      INTEGER IR(0:P-1)
      COMMON /RANDM/ IR,J,K
C
C------------Main-------------------------------------------
C
      IR(J)=XOR(IR(J),IR(K))
      IRND=IR(J)
c      FRND=IRND*Fnorm
c      RANF=FRND
c      IF (ranf.eq.1.0) STOP
      RANF = IRND * Fnorm + Fadd
      J=J+1
      IF(J.GE.P)J=J-P
      K=K+1
      IF(K.GE.P)K=K-P
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE CINIT3(Ndelay)
C
C------------Main-------------------------------------------
       CALL INIT3
       CALL DELAY3(Ndelay)
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE INIT3
C
C------------Variables--------------------------------------
      PARAMETER( IP= 521, IQ= 32 )
      PARAMETER( MACRM= 40, MACRI= 1 )
      COMMON /RANDM/ IW(0:IP-1),JR,KR
      INTEGER IB(0:IP-1)
C------------Main-------------------------------------------
      IX=MACRI
      DO 10 I=0,IP-1
        IX=IX*69069
        IB(I)=ISHFT(IX,-31)
   10 CONTINUE
      JR=0
      KR=IP-IQ
      DO 30 J=0,IP-1
        IWORK=0
        DO 20 I=0,31
          IWORK=IWORK*2+IB(JR)
          IB(JR)=XOR(IB(JR),IB(KR))
          JR=JR+1
          IF (JR.EQ.IP) JR=0
          KR=KR+1
          IF (KR.EQ.IP) KR=0
   20   CONTINUE
        IW(J)=ISHFT(IWORK,-1)
   30 CONTINUE
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE DELAY3(LAMBDA)
C
C------------Variables--------------------------------------
      PARAMETER( IP= 521, IQ= 32 )
      PARAMETER( MACRM= 40, MACRI= 1 )
      COMMON /RANDM/ IW(0:IP-1),JR,KR
      INTEGER IWK(0:2*IP-2),C(0:2*IP-1),IB(0:IP+31)
C------------Main-------------------------------------------
       MU=MACRM
      DO 110 I=0,IP-1
        IWK(I)=IW(I)
  110 CONTINUE
      DO 120 I=IP,2*IP-2
        IWK(I)=XOR(IWK(I-IP),IWK(I-IQ))
  120 CONTINUE
      DO 210 I=0,MU-1
        IB(I)=0
  210 CONTINUE
      M=LAMBDA
      NB=MU-1
  220 CONTINUE
        IF(M.LE.IP-1) GOTO 300
        NB=NB+1
        IB(NB)=MOD(M,2)
        M=M/2
        GOTO 220
  300 DO 310 I=0,IP-1
        C(I)=0
  310 CONTINUE
      C(M)=1
      DO 340 J=NB,0,-1
        DO 320 I=IP-1,0,-1
          C(2*I+IB(J))=C(I)
          C(2*I+1-IB(J))=0
  320   CONTINUE
        DO 330 I=2*IP-1,IP,-1
          C(I-IP)=XOR(C(I-IP),C(I))
          C(I-IQ)=XOR(C(I-IQ),C(I))
  330   CONTINUE
  340 CONTINUE
      DO 420 J=0,IP-1
        IWORK=0
        DO 410 I=0,IP-1
          IWORK=XOR(IWORK,C(I)*IWK(J+I))
  410   CONTINUE
        IW(J)=IWORK
  420 CONTINUE
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
C***************************************************END*****
