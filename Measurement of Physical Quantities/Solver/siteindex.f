C***********************************************************
      SUBROUTINE  SITEINDEX(LV)
C
      INCLUDE 'lattsize.h'
C-----------------Variables---------------------------------
      DIMENSION LV(NS,3,2)
C-----------------Main--------------------------------------
C
      DO 1000  iz = 0, Nz-1
      DO 1000  iy = 0, Ny-1
      DO 1000  ix = 0, Nx-1
C
         is = 1 + ix + iy * Nx + iz * Nx * Ny
C
       LV(is,1,1) = 1 +   MOD(ix+1,Nx) +iy*Nx +iz*Nx*Ny
       LV(is,1,2) = 1 +MOD(ix+Nx-1,Nx) +iy*Nx +iz*Nx*Ny
C
       LV(is,2,1) = 1 +ix +   MOD(iy+1,Ny)*NX +iz*Nx*Ny
       LV(is,2,2) = 1 +ix +MOD(iy+Ny-1,Ny)*NX +iz*Nx*Ny
C 
       LV(is,3,1) = 1 +ix +iy*Nx +   MOD(iz+1,Nz)*Nx*Ny
       LV(is,3,2) = 1 +ix +iy*Nx +MOD(iz+Nz-1,Nz)*Nx*Ny
C 
 1000  CONTINUE
C
C-----------------------------------------------------------
      RETURN
      END
C***********************************************************
C***************************************************END*****


