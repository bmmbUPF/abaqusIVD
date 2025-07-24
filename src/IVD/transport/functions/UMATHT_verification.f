C
C--------------------------------------------------------------------
C    PROPS(1) - SOLUTE DIFFUSITIVTY IN WATER
C    PROPS(2) - SPECHT
C--------------------------------------------------------------------
C
C--------------------------------------------------------------------
C     Diffusivity calculation - (Mackie and Meares 1955)
C--------------------------------------------------------------------
C
      COND = PROPS(1)*((STATEV(1)/(TWO-STATEV(1)))**TWO)
      COND2 = PROPS(1)*STATEV(3)
      SPECHT = PROPS(2)
C
      DUDT = SPECHT
      DU = DUDT*DTEMP
      U = U+DU

C
C--------------------------------------------------------------------
C
C           Solute Concentration Calculation (FLUX)
C              (loop over the integration points)
C
C--------------------------------------------------------------------
C
      DO I=1, NTGRD
            FLUX(I) = -COND*DTEMDX(I) - COND2*TEMP(I)
            DFDG(I,I) = -COND
            DFDT(I) = -COND2      
      END DO