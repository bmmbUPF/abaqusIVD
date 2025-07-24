      PROGRAM PNeq_Main
C-----------------------------------------------------------------
C Main Program to Run the PNeq Subroutine
C-----------------------------------------------------------------
      IMPLICIT NONE

C Declare variables
      DOUBLE PRECISION DTIME, GLUC_int, pH_val, MAG, FREQ
      INTEGER PN_TIME
      INTEGER KSTEP, NDIM, stat

C Declare outputs from PNeq
      DOUBLE PRECISION TNF, IL1B, IL1B_PROT, TNF_PROT
      DOUBLE PRECISION AGG_IMN_ACC, AGG_IL1B_ACC
      DOUBLE PRECISION AGG_TNF_ACC, AGG_IL1B_TNF_ACC
      DOUBLE PRECISION COLI_IMN_ACC, COLI_IL1B_ACC
      DOUBLE PRECISION COLI_TNF_ACC, COLI_IL1B_TNF_ACC
      DOUBLE PRECISION COLII_IMN_ACC, COLII_IL1B_ACC
      DOUBLE PRECISION COLII_TNF_ACC, COLII_IL1B_TNF_ACC
      DOUBLE PRECISION MMP3_IMN_ACC, MMP3_IL1B_ACC
      DOUBLE PRECISION MMP3_TNF_ACC, MMP3_IL1B_TNF_ACC
      DOUBLE PRECISION ADM4_IMN_ACC, ADM4_IL1B_ACC
      DOUBLE PRECISION ADM4_TNF_ACC, ADM4_IL1B_TNF_ACC

C Initialize input variables
      AGG_IMN_ACC = 0.0D0
      AGG_IL1B _ACC= 0.0D0
      AGG_TNF_ACC = 0.0D0
      AGG_IL1B_TNF_ACC = 0.0D0
      
      COLI_IMN_ACC = 0.0D0
      COLI_IL1B_ACC = 0.0D0
      COLI_TNF_ACC = 0.0D0
      COLI_IL1B_TNF_ACC = 0.0D0

      COLII_IMN_ACC = 0.0D0
      COLII_IL1B_ACC = 0.0D0
      COLII_TNF_ACC = 0.0D0
      COLII_IL1B_TNF_ACC = 0.0D0

      MMP3_IMN_ACC = 0.0D0
      MMP3_IL1B_ACC = 0.0D0
      MMP3_TNF_ACC = 0.0D0
      MMP3_IL1B_TNF_ACC = 0.0D0

      ADM4_IMN_ACC = 0.0D0
      ADM4_IL1B_ACC = 0.0D0
      ADM4_TNF_ACC = 0.0D0
      ADM4_IL1B_TNF_ACC = 0.0D0

      KSTEP = 1            ! Current step
      NDIM = 3             ! Number of dimensions
      stat = 0             ! Status variable
      DTIME = 0.01D0       ! Time increment (hours)

      GLUC_int = 1.6138D0     ! Glucose concentration
      pH_val = 7.0240D0       ! pH value
      MAG = 0.25D0        ! Pressure magnitude
      FREQ = 0.0D0         ! Pressure frequency
      PN_TIME = 1          ! Total simulation time (hours)

C Call PNeq to calculate AGG_IMN
      CALL PNeq(
     +           GLUC_int, pH_val,
     +           MAG, FREQ,
     +           IL1B, TNF, IL1B_PROT, TNF_PROT,
     +           AGG_IMN_ACC, AGG_IL1B_ACC,
     +           AGG_TNF_ACC, AGG_IL1B_TNF_ACC,
     +           COLI_IMN_ACC, COLI_IL1B_ACC,
     +           COLI_TNF_ACC, COLI_IL1B_TNF_ACC,
     +           COLII_IMN_ACC, COLII_IL1B_ACC,
     +           COLII_TNF_ACC, COLII_IL1B_TNF_ACC,
     +           MMP3_IMN_ACC, MMP3_IL1B_ACC,
     +           MMP3_TNF_ACC, MMP3_IL1B_TNF_ACC,
     +           ADM4_IMN_ACC, ADM4_IL1B_ACC,
     +           ADM4_TNF_ACC, ADM4_IL1B_TNF_ACC,
     +           PN_TIME)

C Print results
      PRINT *, " "
      PRINT *, "PNeq Results"
      PRINT *, " "
      PRINT *, "Time: ", PN_TIME, " hours"
      PRINT *, " "
      PRINT *, "IL−1β & TNF Expression Levels"
      PRINT *, "IL−1β: ", IL1B
      PRINT *, "TNF: ", TNF
      PRINT *, " "
      PRINT *, "IL−1β & TNF Protein Levels"
      PRINT *, "IL−1β: ", IL1B_PROT
      PRINT *, "TNF: ", TNF_PROT
      PRINT *, " "
      PRINT *, "Aggrecan Expression"
      PRINT *, AGG_IMN_ACC, AGG_IL1B_ACC,
     +         AGG_TNF_ACC, AGG_IL1B_TNF_ACC
      PRINT *, " "
      PRINT *, "Collagen I Expression"
      PRINT *, COLI_IMN_ACC, COLI_IL1B_ACC,
     +         COLI_TNF_ACC, COLI_IL1B_TNF_ACC
      PRINT *, " "
      PRINT *, "Collagen II Expression"
      PRINT *, COLII_IMN_ACC, COLII_IL1B_ACC,
     +         COLII_TNF_ACC, COLII_IL1B_TNF_ACC
      PRINT *, " "
      PRINT *, "MMP-3 Expression"
      PRINT *, MMP3_IMN_ACC, MMP3_IL1B_ACC,
     +         MMP3_TNF_ACC, MMP3_IL1B_TNF_ACC
      PRINT *, " "
      PRINT *, "ADAMTS-4 Expression"
      PRINT *, ADM4_IMN_ACC, ADM4_IL1B_ACC, ADM4_TNF_ACC,
     +         ADM4_IL1B_TNF_ACC
      PRINT *, " "

      STOP
      END PROGRAM PNeq_Main


      INCLUDE './UEL_PNeq.f'
      INCLUDE './UEL_PNeq_Functions.f'