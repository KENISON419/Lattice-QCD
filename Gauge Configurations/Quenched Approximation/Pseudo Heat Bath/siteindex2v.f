C***********************************************************
      SUBROUTINE  SITEINDEX2(LV)
C
      INCLUDE 'lattsize.h'
C---------------Variable------------------------------------
      DIMENSION LV(NS2,3,2,2)
C-----------------Main--------------------------------------
      DO 1000  iz  = 0,Nz-1
      DO 1000  iy  = 0,Ny-1
      DO 1000  ix2 = 0,Nx2-1
C
        is2 = 1 + ix2 + iy * Nx2  + iz * Nx2*Ny
C
       LV(is2,2,1,1)= 1 +ix2 +  MOD(iy+1,Ny)*Nx2  +iz*Nx2*Ny
       LV(is2,2,1,2)= 1 +ix2 +  MOD(iy+1,Ny)*Nx2  +iz*Nx2*Ny
C 
       LV(is2,2,2,1)= 1 +ix2 +MOD(iy+Ny-1,Ny)*Nx2 +iz*Nx2*Ny
       LV(is2,2,2,2)= 1 +ix2 +MOD(iy+Ny-1,Ny)*Nx2 +iz*Nx2*Ny
C 
       LV(is2,3,1,1)= 1 +ix2 +iy*Nx2  +MOD(iz+1,Nz)*Nx2*Ny
       LV(is2,3,1,2)= 1 +ix2 +iy*Nx2  +MOD(iz+1,Nz)*Nx2*Ny
C 
       LV(is2,3,2,1)= 1 +ix2 +iy*Nx2 +MOD(iz+Nz-1,Nz)*Nx2*Ny
       LV(is2,3,2,2)= 1 +ix2 +iy*Nx2 +MOD(iz+Nz-1,Nz)*Nx2*Ny
C 
      IF( MOD(iy+iz,2) .EQ. 0 ) THEN
C
       LV(is2,1,1,1)= is2
       LV(is2,1,2,1)= 1+MOD(ix2+Nx2-1,Nx2) +iy*Nx2+iz*Nx2*Ny
C
       LV(is2,1,1,2)= 1  +MOD(ix2+1,Nx2)   +iy*Nx2+iz*Nx2*Ny
       LV(is2,1,2,2)= is2
C
      ELSE
C
       LV(is2,1,1,1)= 1  +MOD(ix2+1,Nx2)   +iy*Nx2+iz*Nx2*Ny
       LV(is2,1,2,1)= is2
C
       LV(is2,1,1,2)= is2
       LV(is2,1,2,2)= 1+MOD(ix2+Nx2-1,Nx2) +iy*Nx2+iz*Nx2*Ny
C
      ENDIF
C
 1000 CONTINUE
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
      SUBROUTINE ConvSitev(Nv1,Nv2,A,B,Iflg)
C        (e,o) <--> normal oreder converter:
C          A: (e,o) order var. ,  B: normal order
C         Iflg=1:  A <-- B
C         Iflg=2:  A --> B
C     .Vector processor version
C
      INCLUDE 'lattsize.h'
C-----------------Variables---------------------------------
      DIMENSION A(NS2,Nv1,2,Nv2), B(NS,Nv1,Nv2)
C-------------------Main------------------------------------
      IF(Iflg.EQ.1) THEN
C        ------  A <-- B
       DO 1000 i2  = 1, Nv2
        DO 1000 iz  = 0, Nz-1
        DO 1000 iy  = 0, Ny-1
        DO 1000 ix2 = 0, Nx2-1
          Ioe = MOD(iz+iy,2)+1
          Joe = 3-Ioe
            is2 = 1 + ix2 + iy *Nx2 + iz *Nx2*Ny
            isA = 2*is2 -1
            isB = 2*is2
          DO 1200 i1 = 1, Nv1
            A( is2,i1,Ioe,i2 ) = B( isA,i1,i2 )
            A( is2,i1,Joe,i2 ) = B( isB,i1,i2 )
 1200     CONTINUE
 1000  CONTINUE
C
      ELSE IF(Iflg.EQ.2) THEN
C        ------  A --> B
       DO 2000 i2  = 1, Nv2
        DO 2000 iz  = 0, Nz-1
        DO 2000 iy  = 0, Ny-1
        DO 2000 ix2 = 0, Nx2-1
          Ioe = MOD(iz+iy,2)+1
          Joe = 3-Ioe
            is2 = 1 + ix2 + iy *Nx2 + iz *Nx2*Ny
            isA = 2*is2 -1
            isB = 2*is2
         DO 2200 i1 = 1, Nv1
           B( isA,i1,i2 ) = A( is2,i1,Ioe,i2 )
           B( isB,i1,i2 ) = A( is2,i1,Joe,i2 )
 2200    CONTINUE
 2000  CONTINUE
C
      ELSE
        WRITE(*,*) 'Converter: Iflg is wrong.'
      END IF
C
C-----------------------------------------------------------
      RETURN
      END
C***************************************************END*****



