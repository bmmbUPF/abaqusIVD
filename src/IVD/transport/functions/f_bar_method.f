C--------------------------------------------------------------------------
C Take this opportunity to perform calculations at the element centroid. 
C F-bar method
C--------------------------------------------------------------------------
C Obtain shape functions and their local gradients at the element
C centriod, that means ξ=η=ζ=0.0, and nIntPt=1
C shapeFunc_C3D20(numInt,Cip,i_intPT,N,dN,NNODE,NDIM)

      IF(NNODE .EQ. 20) THEN
C Coordinates of the centroid (ξ=η=ζ=0.0)
      Cip0(NDIM) = ZERO
      CALL shapeFunc_C3D20(1,Cip0,1,N0,dNxi,d2Nxi,NNODE,NDIM)

      ELSE
      OPEN(15,FILE=infoFilePath_T_D, status='old', position='append')
      WRITE(15,*) ' '
      WRITE(15,*) '//////////////////////////////////////////////'
      WRITE(15,*) 'ERROR: Number of nodes: NNODE .NE. 20'
      WRITE(15,*) '//////////////////////////////////////////////'
      WRITE(15,*) ' '
      CLOSE(15)
      PRINT *, ' '
      PRINT *, '//////////////////////////////////////////////'
      PRINT *, 'ERROR: Number of nodes: NNODE .NE. 20'
      PRINT *, '//////////////////////////////////////////////'
      PRINT *, ' '
C Terminate the analysis
      CALL XIT
      END IF

C Map shape functions from local to global reference coordinate system
C mapShape3D(NNODE,dNxi,COORDS,dN,DETMapJ,stat)
      CALL mapShape3D(NNODE,dNxi,d2Nxi,COORDS,dN0,d2N0,
     +                DETMapJ0,stat)
      IF (stat .EQ. 0) THEN
            PNEWDT = HALF
            RETURN
      END IF

C Map shape functions from local to global current coordinate system
C this work to check if the determinant of the deformation gradient is zero
C mapShape3D(NNODE,dNxi,COORDS,dN,DETMapJ,stat)
      CALL mapShape3D(NNODE,dNxi,d2Nxi,COORDSC,dN0C,d2N0C,
     +                DETMapJ0C,stat)
      IF (stat .EQ. 0) THEN
            PNEWDT = HALF
            RETURN
      END IF