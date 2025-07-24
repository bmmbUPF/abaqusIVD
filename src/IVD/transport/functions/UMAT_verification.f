C
C--------------------------------------------------------------------
C    STATEV(1)=Water Content (WATER CONTENT) (1/100)
C    STATEV(2)=Cell viability (% wet weight)
C--------------------------------------------------------------------
C
C********************************************************************
C
C             Calculate the actual porosity (Water content)
C                 based on the deformation gradient
C
C********************************************************************
C

C--------------------------------------------------------------------------
C Definitions of the mechanical behaviour of the IVD
C to calculate the water content (NF)
C--------------------------------------------------------------------------

C defining a mechanical behaviour (dummy one in this case)
C Calculation of the D_LAMDA parameter for mechanical behavior.
      D_LAMDA = NU*E/ (ONE+NU)/(ONE-TWO*NU)

C Calculation of the D_MU parameter for mechanical behavior.
      D_MU = E/TWO/(ONE+NU)

C
C--------------------------------------------------------------------
C     defining a mechanical behaviour (dummy one in this case)
C--------------------------------------------------------------------
C
C DDSDDE: Tangent stiffness matrix
C NTENS: Number of components of the stress and strain vectors

      DO I = 1,NTENS
            DO J = 1,NTENS
                  DDSDDE(I,J) = ZERO
            ENDDO
      ENDDO

      DDSDDE(1,1)=D_LAMDA+TWO*D_MU
      DDSDDE(2,2)=DDSDDE(1,1)
      DDSDDE(3,3)=DDSDDE(1,1)
      DDSDDE(4,4)=D_MU
      DDSDDE(5,5)=D_MU
      DDSDDE(6,6)=D_MU
      DDSDDE(1,2)=D_LAMDA
      DDSDDE(1,3)=D_LAMDA
      DDSDDE(2,3)=D_LAMDA      
      DDSDDE(2,1)=DDSDDE(1,2)
      DDSDDE(3,1)=DDSDDE(1,3)
      DDSDDE(3,2)=DDSDDE(2,3)

C
C STRESS: Stress vector
C DSTRAN: Strain increment vector
      DO I = 1,NTENS
            DO J = 1,NTENS
                  STRESS(I)=STRESS(I)+DDSDDE(I,J)*DSTRAN(J)
            END DO
      END DO
C
C DTIME: Time increment
C DFGRD0: Deformation gradient at the beginning of the increment
C DFGRD1: Deformation gradient at the end of the increment
C DFGR: Deformation gradient increment
C
      IF (DTIME .GT. ZERO) THEN
            DO I = 1,3
                  DO J = 1,3
                        DFGR(I,J)=(DFGRD1(I,J)-DFGRD0(I,J))
                  END DO
            END DO
      END IF
C
C  --------------------------------------------
C    Determinant of the deformation gradient
C       
C  is determined from the boundary displcement
C  that for each increment correspond to those
C  of the "global" model (poromechanical one)
C
C  --------------------------------------------
C
C DETDFG: Determinant of the deformation gradient

C Uncomment to use the determinant function developed by Malandrino
C      CALL mDETDFG(DFGRD1,DETDFG,NSHR)

C matInvDet3D(A,A_inv,A_invT,det_A,istat)
      CALL matInvDet3D(DFGRD1,DFGRD1_inv,DFGRD1_invT,DETDFG,stat)
C
C--------------------------------------------------------------------
C       Initial values of water content (porosity) by tissue
C--------------------------------------------------------------------
C
C NF0
c      NF0 = SDV_NF0(NOEL,NPT)
C NF0h
      NF0 = NF0h

      IF (DETDFG .EQ. ZERO) THEN
            NF = NF0
            DETDFG = ONE
      ELSE
            NF =  (NF0-ONE+DETDFG)/DETDFG
      END IF

C Calculate the gradient of the diffusion: TRACE_D = TRACE(F-I)
      TRACE_D = DFGRD1(1,1) + DFGRD1(2,2) + DFGRD1(3,3) - THREE

C      dD_SOLdX = FOUR*NF*(ONE-NF0)*TRACE_D/(DETDFG*(TWO-NF)**THREE)
      dD_SOLdX = (4.0D0*NF*(1.0D0-NF0))*TRACE_D/(DETDFG*(TWO-NF)**3)

C Uncomment the following lines to set ∇D = 0 for testing
C      dD_SOLdX = ZERO
C      CELL_rho_0 = ZERO

C--------------------------------------------------------------------
C
C                 VOLUMETRIC HEAT GENERATION
C   it's the way to create the concentration source or dissipation
C   for the subsequent diffusion-convection analysis
C
C--------------------------------------------------------------------
C
C Initial Values
      glu = 5.0D0 ! VALOR INICIAL GRANDE PARA QUE NO AFECTE
      pH_val = 7.5D0
C
C--------------------------------------------------------------------
C
C         Parameters for cell viability calculation
C
C--------------------------------------------------------------------
C
      diff = glu - 0.5D0
      adiff = ABS(diff)
      summ1 = glu + kk1
      summ2 = glu + kk2
      ALPHA_GLUC = ALPHA_GCF*((diff/summ1)-(adiff/summ2)) ! glucose death rate
      CELL_viab = ONE ! Cell viability for optimal nutrition conditions
      ti = TIME(2) ! Current time
      deti = DTIME ! Current time increment
C
C--------------------------------------------------------------------
C
C          Parameters for reaction calculations
C
C--------------------------------------------------------------------
C
      CC = 1.46D0
      B1 = 4.95D0
      D1 = 4.03D0
C
      BB = pH_val-B1
      DD = D1*BB

      ! IF (KSTEP.eq.1.and.KINC.eq.1) then
      !       pH_val = 7.5D0
      !       BB=pH_val-B1
      !       DD=D1*BB
      !       glu=5.0 ! VALOR INICIAL GRANDE PARA QUE NO AFECTE
      !       CELL_viab = ONE
      ! ELSE
      !       BB=pH_val-B1
      !       DD=D1*BB
      !       glu=5.0
      ! end if
c
C      if ((glu.gt.0.5) .and. (pH_val.gt.6.78)) then
c
C            CELL_viab=STATEV(2)
c
C      else if ((glu.le.0.5) .and. (pH_val.gt.6.78)) then
c
C            CELL_viab = STATEV(2)*EXP(ALPHA_GLUC*deti)
c
C      else if ((glu.le.0.5) .and. (pH_val.le.6.78)) then
c
C            CELL_viab = STATEV(2)*EXP(ALPHA_GLUC*deti)
C            CELL_viab = CELL_viab*EXP(ALPHA_PH*deti)
c
C      else
c
C            CELL_viab = STATEV(2)*EXP(ALPHA_PH*deti)
c
c      end if


C
C--------------------------------------------------------------------
C            Initial cell density by tissue
C--------------------------------------------------------------------
C
C CELL_rho_0
C
      CELL_rho = CELL_rho_0*CELL_viab/DETDFG ! ACTUAL CELL

C********************************************************************
C
C             Metabolic reactions calculation
C              Based on Bibby aat al. (2005)
C  
C********************************************************************

C--------------------------------------------------------------------
C
C                 Reaction for Oxygen by tissue
C
C--------------------------------------------------------------------
C

      AA = (7.28D0*CELL_rho*NF)/(3600.D0*SO2)
      RPL = -(AA*BB*TEMP)/(CC+DD+TEMP)
      DRPLDT = -(AA*BB*(CC+DD))/((CC+DD+TEMP)**TWO)


C--------------------------------------------------------------------------
C Save the state variables at this integ point
      STATEV(1) = NF
      STATEV(2) = CELL_viab
      STATEV(3) = dD_SOLdX
C

C Uncoment the following lines to print the information
!       IF (run_count_share .EQ. numElemInt 
!      +.OR. run_count_share .EQ. 1) THEN !write info
!             PRINT *, ''
!             PRINT *, '--------------------------------------------'
!             PRINT *, 'UMAT'
!             PRINT *, 'run_count = ', run_count_share
!             PRINT *, 'Current step number: ', KSTEP
!             PRINT *, 'Current increment number: ', KINC
!             PRINT *, 'MATNAME = ', MATNAME
!             PRINT *, 'NOEL = ', NOEL
!             PRINT *, 'NPT = ', NPT
!             PRINT *, 'DTIME = ', DTIME
!             PRINT *, 'E',E
!             PRINT *, 'NU',NU
!             PRINT *, 'ALPHA_GCF',ALPHA_GCF
!             PRINT *, 'kk1',kk1
!             PRINT *, 'kk2',kk2
!             PRINT *, 'NF0h',NF0h
!             PRINT *, 'CELL_rho_0',CELL_rho_0
!             PRINT *, 'D_LAMDA',D_LAMDA
!             PRINT *, 'D_MU',D_MU
!             PRINT *, 'NTENS',NTENS
!             PRINT *, 'STRESS(I)=',STRESS(1),STRESS(2),STRESS(3)
!             PRINT *, 'STRESS(I)=',STRESS(4),STRESS(5),STRESS(6)
!             PRINT *, 'DSTRAN(I)=',DSTRAN(1),DSTRAN(2),DSTRAN(3)
!             PRINT *, 'DSTRAN(I)=',DSTRAN(4),DSTRAN(5),DSTRAN(6)
!             PRINT *, 'DETDFG=',DETDFG
!             PRINT *, 'NF=',NF
!             PRINT *, 'STATEV(1)=',STATEV(1)
!             PRINT *, 'COORDS = ', COORDS
! C      CALL XIT
!       END IF
