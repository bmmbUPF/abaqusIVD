C
C ___  ___  ________     ___    ___ ___  ___  ________     ___    ___ ___    ___  _______  ________
C |\  \|\  \|\   __  \   |\  \  /  /|\  \|\  \|\   __  \   |\  \  /  /|\  \  /  /|/  ___  \|\_____  \
C \ \  \\\  \ \  \|\  \  \ \  \/  / | \  \\\  \ \  \|\  \  \ \  \/  / | \  \/  / /__/|_/  /\|____|\ /_
C  \ \   __  \ \   __  \  \ \    / / \ \   __  \ \   __  \  \ \    / / \ \    / /|__|//  / /     \|\  \
C   \ \  \ \  \ \  \ \  \  /     \/   \ \  \ \  \ \  \ \  \  /     \/   /     \/     /  /_/__   __\_\  \
C    \ \__\ \__\ \__\ \__\/  /\   \    \ \__\ \__\ \__\ \__\/  /\   \  /  /\   \    |\________\|\_______\
C     \|__|\|__|\|__|\|__/__/ /\ __\    \|__|\|__|\|__|\|__/__/ /\ __\/__/ /\ __\    \|_______|\|_______|
C                        |__|/ \|__|                       |__|/ \|__||__|/ \|__|
C          01001000 01100001 01111000 01001000 01100001 01111000 01111000 00110010 00110011
C

C **
C ****************************************************************************
C ******    ABAQUS SUBROUTINE                                           ******
C ******                                                                ******
C ******    Diffusion-Reaction Transport IVD Simulation                 ******
C ******    Obs: Coupled with a Mechanical IVD Simulation               ******
C ******    Type: User Element (UEL)                                    ******
C ******    Function:                                                   ******
C ******    - Network Modeling: the Parallel Networks (PN) equations    ******
C ******      to comprehensively estimate dynamics of multicellular     ******
C ******      environments whilst allowing for changes of the nature of ******
C ******      of the links  77 according to the dose and exposure       ******
C ******      time of the cell stimuli.                                 ******
C ******                                                                ******
C ******    Auth: Estefano Muñoz-Moya                                   ******
C ******    Linktree: https://linktr.ee/estefano23                      ******
C ******    Web Page: https://estefano23.github.io/                     ******
C ******    GitHub: estefano23                                          ******
C ******    Email: estefano.munoz.moya@gmail.com                        ******
C ******    HaxHaxx23                                                   ******
C ******    01001000 01100001 01111000                                  ******
C ******    01001000 01100001 01111000 01111000                         ******
C ******    00110010 00110011                                           ******
C ****************************************************************************
C **

C----------------------------------------------------------------------------------------------------------------------------------------------
C Introduction to PN-Methodology
C ------------------------------

C This subroutine will investigate the behaviour of NP cells under global and local stimuli,
C including nutritional conditions and proin-flammatory environments. The model specifically
C incorporates glucose (GLUC_int) and pH (pH_val) as primary (global) stimuli (first-order stimuli),
C and proinflammatory cytokines TNF and IL−1β as secondary (local) stimuli (second-orderstimuli),
C which can result in four distinct proinflammatory cell states (CS): not exposed to proinflammatory cytokines,
C exposed to IL−1β, exposed to TNF and exposed to both, IL−1β and TNF proinflammatory cytokines.
C

C
C-------------------------------------------------------------------------------------------------------------------------------------------------------------------
C Code
C ------------

C--------------------------------------------------------------------
C Subroutine to Compute Theta Values and PN-Denominator
C--------------------------------------------------------------------
      SUBROUTINE Compute_Theta_Values(
     1           aTHETA_IL1B_GLUC, aTHETA_IL1B_LACT,
     2           aTHETA_IL1B_MAG, iTHETA_IL1B_MAG,
     3           aTHETA_IL1B_FREQ, iTHETA_IL1B_FREQ,
     4           aTHETA_TNF_GLUC, aTHETA_TNF_LACT,
     5           aTHETA_TNF_MAG, iTHETA_TNF_MAG,
     6           aTHETA_TNF_FREQ, iTHETA_TNF_FREQ,
     7           aTHETA_AGG_GLUC, aTHETA_AGG_LACT,
     8           iTHETA_AGG_IL1B, iTHETA_AGG_TNF,
     9           aTHETA_AGG_MAG, iTHETA_AGG_MAG,
     1           aTHETA_AGG_FREQ, iTHETA_AGG_FREQ,
     2           aTHETA_COLI_GLUC, aTHETA_COLI_LACT,
     3           aTHETA_COLI_IL1B, iTHETA_COLI_TNF,
     4           aTHETA_COLI_MAG, iTHETA_COLI_MAG,
     5           aTHETA_COLI_FREQ, iTHETA_COLI_FREQ,
     6           aTHETA_COLII_GLUC, aTHETA_COLII_LACT,
     7           iTHETA_COLII_IL1B, iTHETA_COLII_TNF,
     8           aTHETA_COLII_MAG, iTHETA_COLII_MAG,
     9           aTHETA_COLII_FREQ, iTHETA_COLII_FREQ,
     1           aTHETA_MMP3_GLUC, aTHETA_MMP3_LACT,
     2           aTHETA_MMP3_IL1B, aTHETA_MMP3_TNF,
     3           aTHETA_MMP3_MAG, iTHETA_MMP3_MAG,
     4           aTHETA_MMP3_FREQ, iTHETA_MMP3_FREQ,
     5           aTHETA_ADM4_GLUC, aTHETA_ADM4_LACT,
     6           iTHETA_ADM4_IL1B, aTHETA_ADM4_TNF,
     7           aTHETA_ADM4_MAG, iTHETA_ADM4_MAG,
     8           aTHETA_ADM4_FREQ, iTHETA_ADM4_FREQ,
     9           aTHETA_SUM)

C INPUTS AND OUTPUTS:
C All variables are DOUBLE PRECISION representing theta values.
C `aTHETA_SUM` is the computed PN-denominator.

      IMPLICIT NONE

C Input/Output Variables
      DOUBLE PRECISION aTHETA_IL1B_GLUC, aTHETA_IL1B_LACT
      DOUBLE PRECISION aTHETA_IL1B_MAG, iTHETA_IL1B_MAG
      DOUBLE PRECISION aTHETA_IL1B_FREQ, iTHETA_IL1B_FREQ
      DOUBLE PRECISION aTHETA_TNF_GLUC, aTHETA_TNF_LACT
      DOUBLE PRECISION aTHETA_TNF_MAG, iTHETA_TNF_MAG
      DOUBLE PRECISION aTHETA_TNF_FREQ, iTHETA_TNF_FREQ
      DOUBLE PRECISION aTHETA_AGG_GLUC, aTHETA_AGG_LACT
      DOUBLE PRECISION iTHETA_AGG_IL1B, iTHETA_AGG_TNF
      DOUBLE PRECISION aTHETA_AGG_MAG, iTHETA_AGG_MAG
      DOUBLE PRECISION aTHETA_AGG_FREQ, iTHETA_AGG_FREQ
      DOUBLE PRECISION aTHETA_COLI_GLUC, aTHETA_COLI_LACT
      DOUBLE PRECISION aTHETA_COLI_IL1B, iTHETA_COLI_TNF
      DOUBLE PRECISION aTHETA_COLI_MAG, iTHETA_COLI_MAG
      DOUBLE PRECISION aTHETA_COLI_FREQ, iTHETA_COLI_FREQ
      DOUBLE PRECISION aTHETA_COLII_GLUC, aTHETA_COLII_LACT
      DOUBLE PRECISION iTHETA_COLII_IL1B, iTHETA_COLII_TNF
      DOUBLE PRECISION aTHETA_COLII_MAG, iTHETA_COLII_MAG
      DOUBLE PRECISION aTHETA_COLII_FREQ, iTHETA_COLII_FREQ
      DOUBLE PRECISION aTHETA_MMP3_GLUC, aTHETA_MMP3_LACT
      DOUBLE PRECISION aTHETA_MMP3_IL1B, aTHETA_MMP3_TNF
      DOUBLE PRECISION aTHETA_MMP3_MAG, iTHETA_MMP3_MAG
      DOUBLE PRECISION aTHETA_MMP3_FREQ, iTHETA_MMP3_FREQ
      DOUBLE PRECISION aTHETA_ADM4_GLUC, aTHETA_ADM4_LACT
      DOUBLE PRECISION iTHETA_ADM4_IL1B, aTHETA_ADM4_TNF
      DOUBLE PRECISION aTHETA_ADM4_MAG, iTHETA_ADM4_MAG
      DOUBLE PRECISION aTHETA_ADM4_FREQ, iTHETA_ADM4_FREQ
      DOUBLE PRECISION aTHETA_SUM

C--------------------------------------------------------------------
C Initialize Theta Values
C--------------------------------------------------------------------
C IL1B
      aTHETA_IL1B_GLUC = 0.01D0
      aTHETA_IL1B_LACT = 2.8223D0
      aTHETA_IL1B_MAG  = 0.01D0
      iTHETA_IL1B_MAG  = 0.01D0
      aTHETA_IL1B_FREQ = 0.01D0
      iTHETA_IL1B_FREQ = 0.01D0

C TNF
      aTHETA_TNF_GLUC = 0.01D0
      aTHETA_TNF_LACT = 0.01D0
      aTHETA_TNF_MAG  = 0.01D0
      iTHETA_TNF_MAG  = 0.01D0
      aTHETA_TNF_FREQ = 0.01D0
      iTHETA_TNF_FREQ = 0.01D0

C AGG
      aTHETA_AGG_GLUC  = 0.01D0
      aTHETA_AGG_LACT  = 0.0942D0
      iTHETA_AGG_IL1B  = 0.0774D0
      iTHETA_AGG_TNF   = 0.01D0
      aTHETA_AGG_MAG   = 0.3484D0
      iTHETA_AGG_MAG   = 0.3484D0
      aTHETA_AGG_FREQ  = 0.3318D0
      iTHETA_AGG_FREQ  = 0.3318D0

C COLI
      aTHETA_COLI_GLUC = 0.01D0
      aTHETA_COLI_LACT = 0.01D0
      aTHETA_COLI_IL1B = 0.01D0
      iTHETA_COLI_TNF  = 0.1124D0
      aTHETA_COLI_MAG  = 0.8711D0
      iTHETA_COLI_MAG  = 0.8711D0
      aTHETA_COLI_FREQ = 0.9955D0
      iTHETA_COLI_FREQ = 0.9955D0

C COLII
      aTHETA_COLII_GLUC = 0.01D0
      aTHETA_COLII_LACT = 0.0553D0
      iTHETA_COLII_IL1B = 0.01D0
      iTHETA_COLII_TNF  = 0.5807D0
      aTHETA_COLII_MAG  = 0.1394D0
      iTHETA_COLII_MAG  = 0.1394D0
      aTHETA_COLII_FREQ = 0.01D0
      iTHETA_COLII_FREQ = 0.01D0

C MMP3
      aTHETA_MMP3_GLUC = 0.01D0
      aTHETA_MMP3_LACT = 1.0D0
      aTHETA_MMP3_IL1B = 0.3763D0
      aTHETA_MMP3_TNF  = 0.9355D0
      aTHETA_MMP3_MAG  = 0.01D0
      iTHETA_MMP3_MAG  = 0.01D0
      aTHETA_MMP3_FREQ = 0.5124D0
      iTHETA_MMP3_FREQ = 0.5124D0

C ADM4
      aTHETA_ADM4_GLUC = 0.01D0
      aTHETA_ADM4_LACT = 0.1986D0
      iTHETA_ADM4_IL1B = 0.01D0
      aTHETA_ADM4_TNF  = 0.2010D0
      aTHETA_ADM4_MAG  = 0.1394D0
      iTHETA_ADM4_MAG  = 0.1394D0
      aTHETA_ADM4_FREQ = 0.1048D0
      iTHETA_ADM4_FREQ = 0.1048D0

C--------------------------------------------------------------------
C Compute PN-denominator (aTHETA_SUM)
C--------------------------------------------------------------------
      aTHETA_SUM = aTHETA_AGG_GLUC + aTHETA_AGG_LACT +
     1           aTHETA_COLI_GLUC + aTHETA_COLI_LACT +
     2           aTHETA_COLII_GLUC + aTHETA_COLII_LACT +
     3           aTHETA_MMP3_GLUC + aTHETA_MMP3_LACT +
     4           aTHETA_ADM4_GLUC + aTHETA_ADM4_LACT +
     5           aTHETA_COLI_IL1B + aTHETA_MMP3_IL1B +
     6           aTHETA_MMP3_TNF + aTHETA_ADM4_TNF +
     7           aTHETA_AGG_MAG + aTHETA_AGG_FREQ +
     8           aTHETA_COLI_MAG + aTHETA_COLI_FREQ +
     9           aTHETA_COLII_MAG + aTHETA_COLII_FREQ +
     1           aTHETA_MMP3_MAG + aTHETA_MMP3_FREQ +
     2           aTHETA_ADM4_MAG + aTHETA_ADM4_FREQ

C --------------------------------------------------------------------------
      RETURN
      END SUBROUTINE Compute_Theta_Values

C---------------------------------------------------------------------------------------------------------------------------------------------
C Compute the normalized activity of a cell
C ------------------------------------------------------------

      SUBROUTINE Compute_X_S(GLUC_int, pH_val, MAG, FREQ,
     1                      X_S_IL1B_GLUC, X_S_IL1B_LACT, X_S_IL1B_MAG,
     2                      X_S_IL1B_FREQ, X_S_TNF_GLUC, X_S_TNF_LACT,
     3                      X_S_TNF_MAG, X_S_TNF_FREQ, X_S_AGG_GLUC,
     4                      X_S_AGG_LACT, X_S_AGG_MAG, X_S_AGG_FREQ,
     5                      X_S_COLI_GLUC, X_S_COLI_LACT, X_S_COLI_MAG,
     6                      X_S_COLI_FREQ, X_S_COLII_GLUC,
     7                      X_S_COLII_LACT, X_S_COLII_MAG,
     8                      X_S_COLII_FREQ, X_S_MMP3_GLUC,
     9                      X_S_MMP3_LACT, X_S_MMP3_MAG, X_S_MMP3_FREQ,
     1                      X_S_ADM4_GLUC, X_S_ADM4_LACT,
     2                      X_S_ADM4_MAG, X_S_ADM4_FREQ,
     3                      I_A_MAG_gFun, A_I_MAG_gFun,
     4                      I_A_FREQ_gFun, A_I_FREQ_gFun)
C
C ********************************************************************
C
      IMPLICIT NONE
C
C Intermediate Variables
      DOUBLE PRECISION X_S_curve

C INPUT VARIABLES
      DOUBLE PRECISION GLUC_int, pH_val, MAG, FREQ

C OUTPUT VARIABLES
      DOUBLE PRECISION X_S_IL1B_GLUC, X_S_IL1B_LACT
      DOUBLE PRECISION X_S_IL1B_MAG, X_S_IL1B_FREQ
      DOUBLE PRECISION X_S_TNF_GLUC, X_S_TNF_LACT
      DOUBLE PRECISION X_S_TNF_MAG, X_S_TNF_FREQ
      DOUBLE PRECISION X_S_AGG_GLUC, X_S_AGG_LACT
      DOUBLE PRECISION X_S_AGG_MAG, X_S_AGG_FREQ
      DOUBLE PRECISION X_S_COLI_GLUC, X_S_COLI_LACT
      DOUBLE PRECISION X_S_COLI_MAG, X_S_COLI_FREQ
      DOUBLE PRECISION X_S_COLII_GLUC, X_S_COLII_LACT
      DOUBLE PRECISION X_S_COLII_MAG, X_S_COLII_FREQ
      DOUBLE PRECISION X_S_MMP3_GLUC, X_S_MMP3_LACT
      DOUBLE PRECISION X_S_MMP3_MAG, X_S_MMP3_FREQ
      DOUBLE PRECISION X_S_ADM4_GLUC, X_S_ADM4_LACT
      DOUBLE PRECISION X_S_ADM4_MAG, X_S_ADM4_FREQ

C LOCAL VARIABLES
      DOUBLE PRECISION I_A_MAG_gFun, A_I_MAG_gFun
      DOUBLE PRECISION I_A_FREQ_gFun, A_I_FREQ_gFun

C CONSTANT PARAMETERS
      DOUBLE PRECISION, PARAMETER :: ONE = 1.0D0, ZERO = 0.0D0

C--------------------------------------------------------------------------
C Initializing variables
      A_I_MAG_gFun = ZERO
      I_A_MAG_gFun = ZERO
      A_I_FREQ_gFun = ZERO
      I_A_FREQ_gFun = ZERO
      X_S_IL1B_GLUC = ZERO
      X_S_IL1B_LACT = ZERO
      X_S_IL1B_MAG = ZERO
      X_S_IL1B_FREQ = ZERO
      X_S_TNF_GLUC = ZERO
      X_S_TNF_LACT = ZERO
      X_S_TNF_MAG = ZERO
      X_S_TNF_FREQ = ZERO
      X_S_AGG_GLUC = ZERO
      X_S_AGG_LACT = ZERO
      X_S_AGG_MAG = ZERO
      X_S_AGG_FREQ = ZERO
      X_S_COLI_GLUC = ZERO
      X_S_COLI_LACT = ZERO
      X_S_COLI_MAG = ZERO
      X_S_COLI_FREQ = ZERO
      X_S_COLII_GLUC = ZERO
      X_S_COLII_LACT = ZERO
      X_S_COLII_MAG = ZERO
      X_S_COLII_FREQ = ZERO
      X_S_MMP3_GLUC = ZERO
      X_S_MMP3_LACT = ZERO
      X_S_MMP3_MAG = ZERO
      X_S_MMP3_FREQ = ZERO
      X_S_ADM4_GLUC = ZERO
      X_S_ADM4_LACT = ZERO
      X_S_ADM4_MAG = ZERO
      X_S_ADM4_FREQ = ZERO

C
C--------------------------------------------------------------------
C Setup functions (x values)
C--------------------------------------------------------------------
C
C knowledge-based generic functions to approxmiate the effect of magnitude and frequency, respectively

C Magnitude
      A_I_MAG_gFun = (((-1.99999D0 * EXP(14.0D0 * MAG)) /
     1               (EXP(14.0D0 * MAG) + 1.2D0 * (10.0D0**6.0D0))) +
     2               (1.99999D0 / 2.0D0))

      I_A_MAG_gFun = (((1.99999D0 * EXP(14.0D0 * MAG)) /
     1               (EXP(14.0D0 * MAG) + 1.2D0 * (10.0D0**6.0D0))) -
     2               (1.99999D0 / 2.0D0))

C Frequency
      A_I_FREQ_gFun = (((-1.99999D0 * EXP(6.0D0 * FREQ)) /
     1                (EXP(6.0D0 * FREQ) + 6.6D0 * (10.0D0**7.0D0))) +
     2                (1.99999D0 / 2.0D0))

      I_A_FREQ_gFun = (((1.99999D0 * EXP(6.0D0 * FREQ)) /
     1                (EXP(6.0D0 * FREQ) + 6.6D0 * (10.0D0**7.0D0))) -
     2                (1.99999D0 / 2.0D0))

C
C -------- X_S_curve --------
C
C      EXPONENT = C1 * S
C      X_S_curve = C0 * EXP(EXPONENT) /
C     1            (EXP(EXPONENT) + C2 * 10.0D0**C3) + C4

C
C -------- for IL1B --------
C
C GLUCOSE
      IF (GLUC_int .GE. ZERO .AND. GLUC_int .LE. 1.4578D0) THEN
            X_S_IL1B_GLUC = X_S_curve(ONE, 55.0D0, 1.62507D0,
     1                      19.0D0, ZERO, GLUC_int)

      ELSE IF (GLUC_int .GT. 1.4578D0 .AND. GLUC_int .LE. 5.0D0) THEN
            X_S_IL1B_GLUC = X_S_curve(-0.648D0, 4.0D0, 1.314D0,
     1                      6.0D0, 1.000168D0, GLUC_int)

      END IF

C pH
      IF (pH_val .GE. 6.5D0 .AND. pH_val .LE. 7.0372D0) THEN
            X_S_IL1B_LACT = X_S_curve(-ONE / 1.036D0, 55.0D0, 1.94D0,
     1                      4.0D0, ONE, pH_val - 6.5D0)

      ELSE IF (pH_val .GT. 7.0372D0 .AND. pH_val .LE. 7.4D0) THEN
            X_S_IL1B_LACT = X_S_curve(-ONE / 28.7D0, 35.0D0, 5.4D0,
     1                      10.0D0, ONE / 28.7D0, pH_val - 6.5D0)

      END IF

C MAG
      X_S_IL1B_MAG = I_A_MAG_gFun 

C FREQ
      X_S_IL1B_FREQ = I_A_FREQ_gFun;

C
C -------- for TNF ---------
C
C GLUCOSE
      IF (GLUC_int .GE. ZERO .AND. GLUC_int .LE. 0.8134D0) THEN
            X_S_TNF_GLUC = X_S_curve(ONE, 36.0D0, 0.135196D0,
     1                     9.0D0, ZERO, GLUC_int)
      ELSE IF (GLUC_int .GT. 0.8134D0 .AND. GLUC_int .LE. 5.0D0) THEN
            X_S_TNF_GLUC = X_S_curve(-0.5389D0, 24.0D0, 0.53304D0,
     1                     11.0D0, 1.003D0, GLUC_int)

      END IF

C pH
      IF (pH_val .GE. 6.5D0 .AND. pH_val .LE. 6.8042D0) THEN
            X_S_TNF_LACT = X_S_curve(ONE, 32.0D0, 0.6474D0,
     1                     ZERO, ZERO, pH_val - 6.5D0)

      ELSE IF (pH_val .GT. 6.8042D0 .AND. pH_val .LE. 7.1115D0) THEN
            X_S_TNF_LACT = X_S_curve(-ONE, 68.0D0, 2.5D0,
     1                     13.0D0, ONE, pH_val - 6.5D0)

      ELSE IF (pH_val .GT. 7.1115D0 .AND. pH_val .LE. 7.4D0) THEN
            X_S_TNF_LACT = X_S_curve(ONE, 42.0D0, 6.5D0,
     1                     15.0D0, ZERO, pH_val - 6.5D0)

      END IF

C MAG
      X_S_TNF_MAG = I_A_MAG_gFun

C FREQ
      X_S_TNF_FREQ = I_A_FREQ_gFun

C
C -------- for Agg ---------
C
C GLUCOSE
      X_S_AGG_GLUC = X_S_curve(ONE, 18.0D0, 3.0D0,
     1                         4.0D0, ZERO, GLUC_int)

C pH
      IF (6.5D0 .LE. pH_val .AND. pH_val .LE. 6.892D0) THEN
           X_S_AGG_LACT = X_S_curve(ONE, 68.0D0, 2.0D0,
     1                               4.0D0, ZERO, pH_val - 6.5D0)

      ELSE IF (pH_val .GT. 6.892D0 .AND. pH_val .LE. 7.4D0) THEN
           X_S_AGG_LACT = X_S_curve(-ONE, 32.0D0, 5.035D0,
     1                               12.0D0, ONE, pH_val - 6.5D0)

      END IF

C MAG
      X_S_AGG_MAG = A_I_MAG_gFun

C FREQ
      X_S_AGG_FREQ = A_I_FREQ_gFun

C
C ------- for Col-I --------
C
C GLUCOSE
      IF (GLUC_int .GE. ZERO .AND. GLUC_int .LE. 0.5043D0) THEN
           X_S_COLI_GLUC = X_S_curve(ONE, 20.0D0, 1.083D0,
     1                               ZERO, ZERO, GLUC_int)

      ELSE IF (GLUC_int .GT. 0.5043D0 .AND. GLUC_int .LE. 5.0D0) THEN
           X_S_COLI_GLUC = X_S_curve(-ONE, 28.0D0, 3.0D0,
     1                               10.0D0, ONE, GLUC_int)
      END IF

C pH
      IF (pH_val .GE. 6.5D0 .AND. pH_val .LT. 7.1207D0) THEN
           X_S_COLI_LACT = X_S_curve(-ONE, 35.0D0, 6.002D0,
     1                     4.0D0, ONE, pH_val - 6.5D0)

      ELSE IF (pH_val .GE. 7.1207D0 .AND. pH_val .LE. 7.4D0) THEN
           X_S_COLI_LACT = X_S_curve(ONE, 35.0D0, 1.2376D0,
     1                     14.0D0, ZERO, pH_val - 6.5D0)

      END IF

C MAG
      X_S_COLI_MAG = I_A_MAG_gFun

C FREQ
      X_S_COLI_FREQ = I_A_FREQ_gFun

C
C ------- for Col-II -------
C
C GLUCOSE
      X_S_COLII_GLUC = X_S_curve(ONE, 18.0D0, 9.4D0,
     1                           4.0D0, ZERO, GLUC_int)

C pH
      X_S_COLII_LACT = X_S_curve(ONE, 34.0D0, 2.845D0,
     1                 4.0D0, ZERO, pH_val - 6.5D0)

C MAG
      X_S_COLII_MAG = A_I_MAG_gFun

C FREQ
      X_S_COLII_FREQ = A_I_FREQ_gFun

C
C ------- for MMP-3 --------
C
C GLUCOSE
      X_S_MMP3_GLUC = X_S_curve(-ONE, 18.0D0, 2.283D0,
     1                4.0D0, ONE, GLUC_int)

C pH
      X_S_MMP3_LACT = X_S_curve(ONE, 20.0D0, 4.62D0,
     1                9.0D0, ZERO, pH_val - 6.5D0)

C MAG
      X_S_MMP3_MAG = I_A_MAG_gFun

C FREQ
      X_S_MMP3_FREQ = I_A_FREQ_gFun

C
C ------ for ADAMTS-4 ------
C
C GLUCOSE
      IF (GLUC_int .GE. ZERO .AND. GLUC_int .LE. 1.197D0) THEN
             X_S_ADM4_GLUC = X_S_curve(0.4817D0, 38.0D0, ONE,
     1                          4.0D0, ZERO, GLUC_int)

      ELSE IF (GLUC_int .GT. 1.197D0 .AND. GLUC_int .LE. 5.0D0) THEN
             X_S_ADM4_GLUC = X_S_curve(0.520D0, 6.0D0, 4.0D0,
     1                          5.0D0, 0.48D0, GLUC_int)

      END IF

C pH
      X_S_ADM4_LACT = X_S_curve(-ONE, 44.0D0, ONE,
     1                   3.0D0, ONE, pH_val - 6.5D0)

C MAG
      X_S_ADM4_MAG = I_A_MAG_gFun

C FREQ
      X_S_ADM4_FREQ = I_A_FREQ_gFun
         
C --------------------------------------------------------------------------
      RETURN
      END SUBROUTINE Compute_X_S

C --------------------------------------------------------------------------
C Subroutine to Compute Sensitivity Factors (Gamma Values)
C --------------------------------------------------------------------------

      SUBROUTINE timeSenFact(
     1           IL1B_MAG_FACT, IL1B_FREQ_FACT, 
     2           TNF_MAG_FACT, TNF_FREQ_FACT,
     3           AGG_MAG_FACT, AGG_MAG_FACT_2,
     4           AGG_FREQ_FACT, AGG_FREQ_FACT_2,
     5           COLII_MAG_FACT, COLII_MAG_FACT_2,
     6           COLII_FREQ_FACT, COLII_FREQ_FACT_2,
     7           COLI_MAG_FACT, COLI_MAG_FACT_2,
     8           COLI_FREQ_FACT, COLI_FREQ_FACT_2,
     9           MMP3_MAG_FACT, MMP3_MAG_FACT_2,
     1           MMP3_FREQ_FACT, MMP3_FREQ_FACT_2,
     2           ADM4_MAG_FACT, ADM4_MAG_FACT_2,
     3           ADM4_FREQ_FACT, ADM4_FREQ_FACT_2)

C INPUTS:
C   sigma_mag   - Sensitivity scaling factor for magnitude
C   sigma_freq  - Sensitivity scaling factor for frequency

C OUTPUTS:
C   Sensitivity factors for each CA and stimulus

      IMPLICIT NONE

C Output Variables
      DOUBLE PRECISION IL1B_MAG_FACT, IL1B_FREQ_FACT
      DOUBLE PRECISION TNF_MAG_FACT, TNF_FREQ_FACT
      DOUBLE PRECISION AGG_MAG_FACT, AGG_MAG_FACT_2
      DOUBLE PRECISION AGG_FREQ_FACT, AGG_FREQ_FACT_2
      DOUBLE PRECISION COLII_MAG_FACT, COLII_MAG_FACT_2
      DOUBLE PRECISION COLII_FREQ_FACT, COLII_FREQ_FACT_2
      DOUBLE PRECISION COLI_MAG_FACT, COLI_MAG_FACT_2
      DOUBLE PRECISION COLI_FREQ_FACT, COLI_FREQ_FACT_2
      DOUBLE PRECISION MMP3_MAG_FACT, MMP3_MAG_FACT_2
      DOUBLE PRECISION MMP3_FREQ_FACT, MMP3_FREQ_FACT_2
      DOUBLE PRECISION ADM4_MAG_FACT, ADM4_MAG_FACT_2
      DOUBLE PRECISION ADM4_FREQ_FACT, ADM4_FREQ_FACT_2

C CONSTANT PARAMETERS
      DOUBLE PRECISION, PARAMETER :: ONE = 1.0D0, ZERO = 0.0D0

C ------ sigma value (give biological sense to purely mathematically obtained generic functions) ------
      DOUBLE PRECISION, PARAMETER :: sigma_mag = 1.5D0
      DOUBLE PRECISION, PARAMETER :: sigma_freq = 1.5D0

C--------------------------------------------------------------------------
C Initializing variables
      IL1B_MAG_FACT = ZERO
      IL1B_FREQ_FACT = ZERO
      TNF_MAG_FACT = ZERO
      TNF_FREQ_FACT = ZERO
      AGG_MAG_FACT = ZERO
      AGG_MAG_FACT_2 = ZERO
      AGG_FREQ_FACT = ZERO
      AGG_FREQ_FACT_2 = ZERO
      COLII_MAG_FACT = ZERO
      COLII_MAG_FACT_2 = ZERO
      COLII_FREQ_FACT = ZERO
      COLII_FREQ_FACT_2 = ZERO
      COLI_MAG_FACT = ZERO
      COLI_MAG_FACT_2 = ZERO
      COLI_FREQ_FACT = ZERO
      COLI_FREQ_FACT_2 = ZERO
      MMP3_MAG_FACT = ZERO
      MMP3_MAG_FACT_2 = ZERO
      MMP3_FREQ_FACT = ZERO
      MMP3_FREQ_FACT_2 = ZERO
      ADM4_MAG_FACT = ZERO
      ADM4_MAG_FACT_2 = ZERO
      ADM4_FREQ_FACT = ZERO
      ADM4_FREQ_FACT_2 = ZERO

C --------------------------------------------------------------------------
C Compute Gamma Values (Sensitivity Factors)
C --------------------------------------------------------------------------

C IL1B
      IL1B_MAG_FACT = 0.3143D0 * sigma_mag
      IL1B_FREQ_FACT = 0.3143D0 * sigma_freq

C TNF
      TNF_MAG_FACT = 0.3143D0 * sigma_mag
      TNF_FREQ_FACT = 0.3143D0 * sigma_freq

C Agg
      AGG_MAG_FACT = 3.1429D0 * sigma_mag
      AGG_MAG_FACT_2 = 1.0000D0 * sigma_mag
      AGG_FREQ_FACT = 3.1429D0 * sigma_freq
      AGG_FREQ_FACT_2 = 1.0000D0 * sigma_freq

C ColII
      COLII_MAG_FACT = 1.7500D0 * sigma_mag
      COLII_MAG_FACT_2 = 0.5568D0 * sigma_mag
      COLII_FREQ_FACT = 1.7500D0 * sigma_freq
      COLII_FREQ_FACT_2 = 0.5568D0 * sigma_freq

C ColI
      COLI_MAG_FACT = 2.7647D0 * sigma_mag
      COLI_MAG_FACT_2 = 0.8797D0 * sigma_mag
      COLI_FREQ_FACT = 2.7647D0 * sigma_freq
      COLI_FREQ_FACT_2 = 0.8797D0 * sigma_freq

C MMP3
      MMP3_MAG_FACT = 2.1571D0 * sigma_mag
      MMP3_MAG_FACT_2 = 0.6864D0 * sigma_mag
      MMP3_FREQ_FACT = 2.1571D0 * sigma_freq
      MMP3_FREQ_FACT_2 = 0.6864D0 * sigma_freq

C ADM4
      ADM4_MAG_FACT = 3.1481D0 * sigma_mag
      ADM4_MAG_FACT_2 = 1.0017D0 * sigma_mag
      ADM4_FREQ_FACT = 3.1481D0 * sigma_freq
      ADM4_FREQ_FACT_2 = 1.0017D0 * sigma_freq

C --------------------------------------------------------------------------
      RETURN
      END SUBROUTINE timeSenFact

C--------------------------------------------------------------------------------------------------------------------------------------------- 
C Evaluate the IL−1β expression
C -----------------------------

      SUBROUTINE PNeq_IL1B(IL1B, TICKS,
     1           FREQ_0, FREQ_0_time,
     2           MAG_0, MAG_0_time,
     3           GLUC_int, pH_val, MAG, FREQ,
     4           aTHETA_IL1B_GLUC, aTHETA_IL1B_LACT,
     5           aTHETA_IL1B_MAG, iTHETA_IL1B_MAG,
     6           aTHETA_IL1B_FREQ, iTHETA_IL1B_FREQ,
     7           AC_int_IL1B_dir, IC_int_IL1B_dir,
     8           X_S_IL1B_GLUC, X_S_IL1B_LACT,
     9           X_S_IL1B_MAG, X_S_IL1B_FREQ,
     1           A_I_MAG_ANA, A_I_MAG_CAT_1, A_I_MAG_CAT_2,
     2           I_A_MAG_ANA, I_A_MAG_CAT_1, I_A_MAG_CAT_2,
     3           A_I_FREQ_ANA, A_I_FREQ_CAT_1, A_I_FREQ_CAT_2,
     4           I_A_FREQ_ANA, I_A_FREQ_CAT_1, I_A_FREQ_CAT_2,
     5           IL1B_MAG_FACT, IL1B_FREQ_FACT,
     6           X_S_IL1B_MAG_COMPL, X_S_IL1B_FREQ_COMPL)

C VARIABLES:
C X_S: SENSITIVITY OF A CELL ACTIVITY (CA) TO A STIMULUS CONCENTRATION
C THETA: SENSITIVITY OF A CELL ACTIVITY (CA) TO A STIMULUS TYPE

C NECESSARY VARIABLES OF THE PROBLEM:
C GLUCOSE CONCENTRATIONS: GLUC_int
C pH: pH_val
C MAGNITUDE: MAG
C FREQUENCY: FREQ 

C OUTPUTS:
C IL1B: EXPRESSION OF IL−1β

C ********************************************************************
C
      IMPLICIT NONE
C
C FUNCTION
      DOUBLE PRECISION ROUND

C INTEGER VARIABLES
      INTEGER TICKS

C SPECIAL CASE FOR FREQ = 0
      DOUBLE PRECISION FREQ_0, FREQ_0_time
      DOUBLE PRECISION MAG_0, MAG_0_time

C Intermediate Variables
      DOUBLE PRECISION X_S_IL1B_GLUC, X_S_IL1B_LACT
      DOUBLE PRECISION X_S_IL1B_MAG, X_S_IL1B_FREQ

C Input and Output Variables
      DOUBLE PRECISION GLUC_int, pH_val, MAG, FREQ
      DOUBLE PRECISION IL1B

C Constants for Sensitivity Factors
      DOUBLE PRECISION aTHETA_IL1B_GLUC, aTHETA_IL1B_LACT
      DOUBLE PRECISION aTHETA_IL1B_MAG, iTHETA_IL1B_MAG
      DOUBLE PRECISION aTHETA_IL1B_FREQ, iTHETA_IL1B_FREQ
      DOUBLE PRECISION AC_int_IL1B_dir, IC_int_IL1B_dir

C Time dependent funtions
      DOUBLE PRECISION A_I_MAG_ANA, A_I_MAG_CAT_1, A_I_MAG_CAT_2
      DOUBLE PRECISION I_A_MAG_ANA, I_A_MAG_CAT_1, I_A_MAG_CAT_2
      DOUBLE PRECISION A_I_FREQ_ANA, A_I_FREQ_CAT_1, A_I_FREQ_CAT_2
      DOUBLE PRECISION I_A_FREQ_ANA, I_A_FREQ_CAT_1, I_A_FREQ_CAT_2

C Time sensitivity factors
      DOUBLE PRECISION IL1B_MAG_FACT, IL1B_FREQ_FACT

C Time sensitivity range
      DOUBLE PRECISION IL1B_MAG_ANA_RANGE
      DOUBLE PRECISION IL1B_FREQ_ANA_RANGE

C Completive factors
      DOUBLE PRECISION X_S_IL1B_MAG_COMPL
      DOUBLE PRECISION X_S_IL1B_MAG_DEPt
      DOUBLE PRECISION X_S_IL1B_FREQ_COMPL
      DOUBLE PRECISION X_S_IL1B_FREQ_DEPt

C CONSTANT PARAMETERS:
C SIMPLE DOUBLE PRECISION NUMBERS
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.D0, ONE=1.0D0

C Initialize the variables
      IL1B_MAG_ANA_RANGE = ZERO
      X_S_IL1B_MAG_DEPt = ZERO
      IL1B_FREQ_ANA_RANGE = ZERO
      X_S_IL1B_FREQ_DEPt = ZERO

C Evaluate the IL−1β expression

C---------------------------------------------------------------------------
C Magnitude <= 0.998 and Frequency <= 3.0 (Magnitude and Frequency anabolic)
C---------------------------------------------------------------------------

      IF (MAG .LE. 0.998D0 .AND. FREQ .LE. 3.0D0 .AND.
     1    FREQ .NE. ZERO) THEN

C Magnitude Anabolic
            IL1B_MAG_ANA_RANGE = I_A_MAG_ANA * IL1B_MAG_FACT

            IF (IL1B_MAG_ANA_RANGE .GE. ABS(X_S_IL1B_MAG)) THEN
                  IL1B_MAG_ANA_RANGE = ABS(X_S_IL1B_MAG)
                  X_S_IL1B_MAG_COMPL = X_S_IL1B_MAG_COMPL +
     1                                  (I_A_MAG_CAT_1 * IL1B_MAG_FACT)

                  IF (X_S_IL1B_MAG_COMPL .GT. ONE) THEN
                        X_S_IL1B_MAG_COMPL = ONE
                  END IF
            ELSE
                  X_S_IL1B_MAG_COMPL = ZERO
            END IF
            X_S_IL1B_MAG_DEPt = ABS(X_S_IL1B_MAG + IL1B_MAG_ANA_RANGE)

C Frequency Anabolic
            IL1B_FREQ_ANA_RANGE = I_A_FREQ_ANA * IL1B_FREQ_FACT

            IF (IL1B_FREQ_ANA_RANGE .GE. ABS(X_S_IL1B_FREQ)) THEN
                  IL1B_FREQ_ANA_RANGE = ABS(X_S_IL1B_FREQ)
                  X_S_IL1B_FREQ_COMPL = X_S_IL1B_FREQ_COMPL +
     1            (I_A_FREQ_CAT_1 * IL1B_FREQ_FACT)

                  IF (X_S_IL1B_FREQ_COMPL .GT. ONE) THEN
                        X_S_IL1B_FREQ_COMPL = ONE
                  END IF
            ELSE
                  X_S_IL1B_FREQ_COMPL = ZERO
            END IF

            X_S_IL1B_FREQ_DEPt = ABS(X_S_IL1B_FREQ +
     1                           IL1B_FREQ_ANA_RANGE)

            ! IL1B mRNA Expression
            IL1B = ROUND(
     1             (((ONE + AC_int_IL1B_dir) / (AC_int_IL1B_dir)) *
     2             (((aTHETA_IL1B_LACT * X_S_IL1B_LACT) +
     3             (aTHETA_IL1B_GLUC * X_S_IL1B_GLUC) +
     4             (aTHETA_IL1B_MAG * X_S_IL1B_MAG_COMPL) +
     5             (aTHETA_IL1B_FREQ * X_S_IL1B_FREQ_COMPL)) / (ONE +
     6             (aTHETA_IL1B_LACT * X_S_IL1B_LACT) +
     7             (aTHETA_IL1B_GLUC * X_S_IL1B_GLUC) +
     8             (aTHETA_IL1B_MAG * X_S_IL1B_MAG_COMPL) +
     9             (aTHETA_IL1B_FREQ * X_S_IL1B_FREQ_COMPL)))) * (ONE -
     1             (((iTHETA_IL1B_MAG + iTHETA_IL1B_FREQ) /
     2             (aTHETA_IL1B_GLUC + aTHETA_IL1B_LACT +
     3             iTHETA_IL1B_MAG + iTHETA_IL1B_FREQ)) * ((ONE +
     4             IC_int_IL1B_dir) / (IC_int_IL1B_dir)) *
     5             (((iTHETA_IL1B_MAG * X_S_IL1B_MAG_DEPt) +
     6             (iTHETA_IL1B_FREQ * X_S_IL1B_FREQ_DEPt)) / (ONE +
     7             (iTHETA_IL1B_MAG * X_S_IL1B_MAG_DEPt) +
     8             (iTHETA_IL1B_FREQ * X_S_IL1B_FREQ_DEPt)))))
     9             , 4)

C---------------------------------------------------------------------------
C Magnitude <= 0.998 and Frequency == 0 (Magnitude anabolic, Frequency 0)
C---------------------------------------------------------------------------

      ELSE IF (MAG .LE. 0.998D0 .AND. FREQ .EQ. ZERO) THEN

C Magnitude Anabolic Range
            IL1B_MAG_ANA_RANGE = I_A_MAG_ANA * IL1B_MAG_FACT

            IF (IL1B_MAG_ANA_RANGE .GE. ABS(X_S_IL1B_MAG)) THEN
                  IL1B_MAG_ANA_RANGE = ABS(X_S_IL1B_MAG)
                  X_S_IL1B_MAG_COMPL = X_S_IL1B_MAG_COMPL +
     1                                  (I_A_MAG_CAT_1 * IL1B_MAG_FACT)

            IF (X_S_IL1B_MAG_COMPL .GT. ONE) THEN
                  X_S_IL1B_MAG_COMPL = ONE
            END IF
            
            ELSE
                  X_S_IL1B_MAG_COMPL = ZERO
            END IF

            X_S_IL1B_MAG_DEPt = ABS(X_S_IL1B_MAG + IL1B_MAG_ANA_RANGE)

C Frequency Anabolic Range
            IL1B_FREQ_ANA_RANGE = FREQ_0_TIME *
     1                            DBLE(TICKS) * IL1B_FREQ_FACT

            IF (IL1B_FREQ_ANA_RANGE .GE. FREQ_0) THEN
                  IL1B_FREQ_ANA_RANGE = FREQ_0
                  X_S_IL1B_FREQ_COMPL = X_S_IL1B_FREQ_COMPL +
     1                                   (FREQ_0_TIME * IL1B_FREQ_FACT)

                  IF (X_S_IL1B_FREQ_COMPL .GT. ONE) THEN
                  X_S_IL1B_FREQ_COMPL = ONE
                  END IF
            ELSE
                  X_S_IL1B_FREQ_COMPL = ZERO
            END IF

            X_S_IL1B_FREQ_DEPt = FREQ_0 - IL1B_FREQ_ANA_RANGE

C IL1b mRNA Expression

            IL1B = ROUND(
     1             (((ONE + AC_int_IL1B_dir) / (AC_int_IL1B_dir)) *
     2             (((aTHETA_IL1B_LACT * X_S_IL1B_LACT) +
     3             (aTHETA_IL1B_GLUC * X_S_IL1B_GLUC) +
     4             (aTHETA_IL1B_MAG * X_S_IL1B_MAG_COMPL) +
     5             (aTHETA_IL1B_FREQ * X_S_IL1B_FREQ_COMPL)) / (ONE +
     6             (aTHETA_IL1B_LACT * X_S_IL1B_LACT) +
     7             (aTHETA_IL1B_GLUC * X_S_IL1B_GLUC) +
     8             (aTHETA_IL1B_MAG * X_S_IL1B_MAG_COMPL) +
     9             (aTHETA_IL1B_FREQ * X_S_IL1B_FREQ_COMPL)))) * (ONE -
     1             (((iTHETA_IL1B_MAG + iTHETA_IL1B_FREQ) /
     2             (aTHETA_IL1B_GLUC + aTHETA_IL1B_LACT +
     3             iTHETA_IL1B_MAG + iTHETA_IL1B_FREQ)) * ((ONE +
     4             IC_int_IL1B_dir) / (IC_int_IL1B_dir)) *
     5             (((iTHETA_IL1B_MAG * X_S_IL1B_MAG_DEPt) +
     6             (iTHETA_IL1B_FREQ * X_S_IL1B_FREQ_DEPt)) / (ONE +
     7             (iTHETA_IL1B_MAG * X_S_IL1B_MAG_DEPt) +
     8             (iTHETA_IL1B_FREQ * X_S_IL1B_FREQ_DEPt)))))
     9             , 4)

C---------------------------------------------------------------------------
C Magnitude <= 0.998 and Frequency > 3.0 (Magnitude anabolic, Frequency catabolic)
C---------------------------------------------------------------------------

      ELSE IF (MAG .LE. 0.998D0 .AND. FREQ .GT. 3.0D0) THEN

C Magnitude Anabolic Range
            IL1B_MAG_ANA_RANGE = I_A_MAG_ANA * IL1B_MAG_FACT

            IF (IL1B_MAG_ANA_RANGE .GE. ABS(X_S_IL1B_MAG)) THEN
                  IL1B_MAG_ANA_RANGE = ABS(X_S_IL1B_MAG)
                  X_S_IL1B_MAG_COMPL = X_S_IL1B_MAG_COMPL +
     1                                  (I_A_MAG_CAT_1 * IL1B_MAG_FACT)

                  IF (X_S_IL1B_MAG_COMPL .GT. ONE) THEN
                        X_S_IL1B_MAG_COMPL = ONE
                  END IF
            ELSE
                  X_S_IL1B_MAG_COMPL = ZERO
            END IF

            X_S_IL1B_MAG_DEPt = ABS(X_S_IL1B_MAG + IL1B_MAG_ANA_RANGE)

C Frequency Catabolic Range
            X_S_IL1B_FREQ_DEPt = X_S_IL1B_FREQ ! no latency time of proinflammatory cytokines, value of generic function is taken

C IL1b mRNA Expression
            IL1B = ROUND(
     1             (((ONE + AC_int_IL1B_dir) / (AC_int_IL1B_dir)) *
     2             (((aTHETA_IL1B_LACT * X_S_IL1B_LACT) +
     3             (aTHETA_IL1B_GLUC * X_S_IL1B_GLUC) +
     4             (aTHETA_IL1B_MAG * X_S_IL1B_MAG_COMPL) +
     5             (aTHETA_IL1B_FREQ * X_S_IL1B_FREQ_DEPt)) / (ONE +
     6             (aTHETA_IL1B_LACT * X_S_IL1B_LACT) +
     7             (aTHETA_IL1B_GLUC * X_S_IL1B_GLUC) +
     8             (aTHETA_IL1B_MAG * X_S_IL1B_MAG_COMPL) +
     9             (aTHETA_IL1B_FREQ * X_S_IL1B_FREQ_DEPt)))) * (ONE -
     1             (((iTHETA_IL1B_MAG) /
     2             (aTHETA_IL1B_GLUC + aTHETA_IL1B_LACT +
     3             iTHETA_IL1B_MAG + aTHETA_IL1B_FREQ)) * ((ONE +
     4             IC_int_IL1B_dir) / (IC_int_IL1B_dir)) *
     5             (((iTHETA_IL1B_MAG * X_S_IL1B_MAG_DEPt)) / (ONE +
     6             (iTHETA_IL1B_MAG * X_S_IL1B_MAG_DEPt))))), 4)

C--------------------------------------------------------------------------------
C Magnitude >= 1 and Frequency <= 3 (Magnituide catabolic and Frequency anabolic)
C--------------------------------------------------------------------------------

      ELSE IF (MAG .GE. ONE .AND. FREQ .LE. 3.0D0 .AND.
     1         FREQ .NE. ZERO) THEN

C Magnitude Catabolic Range
            X_S_IL1B_MAG_DEPt = I_A_MAG_ANA

C Frequency Anabolic Range
            IL1B_FREQ_ANA_RANGE = I_A_FREQ_ANA * IL1B_FREQ_FACT

            IF (IL1B_FREQ_ANA_RANGE .GE. ABS(X_S_IL1B_FREQ)) THEN
                  IL1B_FREQ_ANA_RANGE = ABS(X_S_IL1B_FREQ)
                  X_S_IL1B_FREQ_COMPL = X_S_IL1B_FREQ_COMPL +
     1                                 (I_A_FREQ_CAT_1 * IL1B_FREQ_FACT)

                  IF (X_S_IL1B_FREQ_COMPL .GT. ONE) THEN
                        X_S_IL1B_FREQ_COMPL = ONE
                  END IF
            ELSE
                  X_S_IL1B_FREQ_COMPL = ZERO
            END IF

            X_S_IL1B_FREQ_DEPt = ABS(X_S_IL1B_FREQ +
     1                           IL1B_FREQ_ANA_RANGE)

C IL1b mRNA Expression
            IL1B = ROUND(
     1             (((ONE + AC_int_IL1B_dir) / (AC_int_IL1B_dir)) *
     2             (((aTHETA_IL1B_LACT * X_S_IL1B_LACT) +
     3             (aTHETA_IL1B_GLUC * X_S_IL1B_GLUC) +
     4             (aTHETA_IL1B_MAG * X_S_IL1B_MAG_DEPt) +
     5             (aTHETA_IL1B_FREQ * X_S_IL1B_FREQ_COMPL)) / (ONE +
     6             (aTHETA_IL1B_LACT * X_S_IL1B_LACT) +
     7             (aTHETA_IL1B_GLUC * X_S_IL1B_GLUC) +
     8             (aTHETA_IL1B_MAG * X_S_IL1B_MAG_DEPt) +
     9             (aTHETA_IL1B_FREQ * X_S_IL1B_FREQ_COMPL)))) * (ONE -
     1             (((iTHETA_IL1B_FREQ) /
     2             (aTHETA_IL1B_GLUC + aTHETA_IL1B_LACT +
     3             aTHETA_IL1B_MAG + iTHETA_IL1B_FREQ)) * ((ONE +
     4             IC_int_IL1B_dir) / (IC_int_IL1B_dir)) *
     5             (((iTHETA_IL1B_FREQ * X_S_IL1B_FREQ_DEPt)) / (ONE +
     6             (iTHETA_IL1B_FREQ * X_S_IL1B_FREQ_DEPt))))), 4)

C---------------------------------------------------------------------------
C Magnitude >= 1 and Frequency == 0 (Magnitude canabolic, Frequency 0)
C---------------------------------------------------------------------------

      ELSE IF (MAG .GE. ONE .AND. FREQ .EQ. ZERO) THEN

C Magnitude Catabolic Range
            X_S_IL1B_MAG_DEPt = X_S_IL1B_MAG  ! no latency time of proinflammatory cytokines, value of generic function is taken

C Frequency Anabolic Range
            IL1B_FREQ_ANA_RANGE = FREQ_0_TIME *
     1                            DBLE(TICKS) * IL1B_FREQ_FACT

            IF (IL1B_FREQ_ANA_RANGE .GE. FREQ_0) THEN
                  IL1B_FREQ_ANA_RANGE = FREQ_0
                  X_S_IL1B_FREQ_COMPL = X_S_IL1B_FREQ_COMPL +
     1                                   (FREQ_0_TIME * IL1B_FREQ_FACT)

                  IF (X_S_IL1B_FREQ_COMPL .GT. ONE) THEN
                        X_S_IL1B_FREQ_COMPL = ONE
                  END IF
            ELSE
                  X_S_IL1B_FREQ_COMPL = ZERO
            END IF

            X_S_IL1B_FREQ_DEPt = FREQ_0 - IL1B_FREQ_ANA_RANGE

C IL1b mRNA Expression
            IL1B = ROUND(
     1             (((ONE + AC_int_IL1B_dir) / (AC_int_IL1B_dir)) *
     2             (((aTHETA_IL1B_LACT * X_S_IL1B_LACT) +
     3             (aTHETA_IL1B_GLUC * X_S_IL1B_GLUC) +
     4             (aTHETA_IL1B_MAG * X_S_IL1B_MAG_DEPt) +
     5             (aTHETA_IL1B_FREQ * X_S_IL1B_FREQ_COMPL)) / (ONE + 
     6             (aTHETA_IL1B_LACT * X_S_IL1B_LACT) +
     7             (aTHETA_IL1B_GLUC * X_S_IL1B_GLUC) +
     8             (aTHETA_IL1B_MAG * X_S_IL1B_MAG_DEPt) +
     9             (aTHETA_IL1B_FREQ * X_S_IL1B_FREQ_COMPL)))) * (ONE - 
     1             (((iTHETA_IL1B_FREQ) /
     2             (aTHETA_IL1B_GLUC + aTHETA_IL1B_LACT +
     3             aTHETA_IL1B_MAG + iTHETA_IL1B_FREQ)) * ((ONE +
     4             IC_int_IL1B_dir) / (IC_int_IL1B_dir)) *
     5             (((iTHETA_IL1B_FREQ * X_S_IL1B_FREQ_DEPt)) / (ONE + 
     6             (iTHETA_IL1B_FREQ * X_S_IL1B_FREQ_DEPt))))), 4)

C-----------------------------------------------------------------------------
C Magnitude >= 1 and Frequency >= 3 (Magnitude catabolic, Frequency catabolic)
C-----------------------------------------------------------------------------

      ELSE IF (MAG .GE. ONE .AND. FREQ .GT. 3.0D0) THEN

C Magnitude Catabolic Range
            X_S_IL1B_MAG_DEPt = X_S_IL1B_MAG

C Frequency Catabolic Range
            X_S_IL1B_FREQ_DEPt = X_S_IL1B_FREQ

C IL1b mRNA Expression
            IL1B = ROUND(
     1             (((ONE + AC_int_IL1B_dir) / (AC_int_IL1B_dir)) *
     2             (((aTHETA_IL1B_LACT * X_S_IL1B_LACT) +
     3             (aTHETA_IL1B_GLUC * X_S_IL1B_GLUC) +
     4             (aTHETA_IL1B_MAG * X_S_IL1B_MAG_DEPt) +
     5             (aTHETA_IL1B_FREQ * X_S_IL1B_FREQ_DEPt)) / (ONE +
     6             (aTHETA_IL1B_LACT * X_S_IL1B_LACT) +
     7             (aTHETA_IL1B_GLUC * X_S_IL1B_GLUC) +
     8             (aTHETA_IL1B_MAG * X_S_IL1B_MAG_DEPt) +
     9             (aTHETA_IL1B_FREQ * X_S_IL1B_FREQ_DEPt)))), 4)

      END IF

C--------------------------------------------------------------------------
      RETURN
      END SUBROUTINE PNeq_IL1B

C--------------------------------------------------------------------------------------------------------------------------------------------- 
C Evaluate the TNF expression
C -----------------------------

      SUBROUTINE PNeq_TNF(TNF, TICKS,
     1           FREQ_0, FREQ_0_time,
     2           MAG_0, MAG_0_time,
     3           GLUC_int, pH_val, MAG, FREQ,
     4           aTHETA_TNF_GLUC, aTHETA_TNF_LACT,
     5           aTHETA_TNF_MAG, iTHETA_TNF_MAG,
     6           aTHETA_TNF_FREQ, iTHETA_TNF_FREQ,
     7           AC_int_TNF_dir, IC_int_TNF_dir,
     8           X_S_TNF_GLUC, X_S_TNF_LACT,
     9           X_S_TNF_MAG, X_S_TNF_FREQ,
     1           A_I_MAG_ANA, A_I_MAG_CAT_1, A_I_MAG_CAT_2,
     2           I_A_MAG_ANA, I_A_MAG_CAT_1, I_A_MAG_CAT_2,
     3           A_I_FREQ_ANA, A_I_FREQ_CAT_1, A_I_FREQ_CAT_2,
     4           I_A_FREQ_ANA, I_A_FREQ_CAT_1, I_A_FREQ_CAT_2,
     5           TNF_MAG_FACT, TNF_FREQ_FACT,
     6           X_S_TNF_MAG_COMPL, X_S_TNF_FREQ_COMPL)

C VARIABLES:
C X_S: SENSITIVITY OF A CELL ACTIVITY (CA) TO A STIMULUS CONCENTRATION
C THETA: SENSITIVITY OF A CELL ACTIVITY (CA) TO A STIMULUS TYPE

C NECESSARY VARIABLES OF THE PROBLEM:
C GLUCOSE CONCENTRATIONS: GLUC_int
C pH: pH_val
C MAGNITUDE: MAG
C FREQUENCY: FREQ 

C OUTPUTS:
C TNF: EXPRESSION OF TNF

C ********************************************************************
C
      IMPLICIT NONE
C
C FUNCTION
      DOUBLE PRECISION ROUND

C INTEGER VARIABLES
      INTEGER TICKS

C SPECIAL CASE FOR FREQ = 0
      DOUBLE PRECISION FREQ_0, FREQ_0_time
      DOUBLE PRECISION MAG_0, MAG_0_time

C Input and Output Variables
      DOUBLE PRECISION TNF, GLUC_int, pH_val, MAG, FREQ

C Constants for Sensitivity Factors
      DOUBLE PRECISION aTHETA_TNF_GLUC, aTHETA_TNF_LACT
      DOUBLE PRECISION aTHETA_TNF_MAG, iTHETA_TNF_MAG
      DOUBLE PRECISION aTHETA_TNF_FREQ, iTHETA_TNF_FREQ
      DOUBLE PRECISION AC_int_TNF_dir, IC_int_TNF_dir

C Intermediate Variables
      DOUBLE PRECISION X_S_TNF_GLUC, X_S_TNF_LACT
      DOUBLE PRECISION X_S_TNF_MAG, X_S_TNF_FREQ
      DOUBLE PRECISION TNF_MAG_ANA_RANGE, TNF_FREQ_ANA_RANGE
      DOUBLE PRECISION X_S_TNF_MAG_COMPL, X_S_TNF_FREQ_COMPL
      DOUBLE PRECISION X_S_TNF_MAG_DEPt, X_S_TNF_FREQ_DEPt

C Time dependent funtions
      DOUBLE PRECISION A_I_MAG_ANA, A_I_MAG_CAT_1, A_I_MAG_CAT_2
      DOUBLE PRECISION I_A_MAG_ANA, I_A_MAG_CAT_1, I_A_MAG_CAT_2
      DOUBLE PRECISION A_I_FREQ_ANA, A_I_FREQ_CAT_1, A_I_FREQ_CAT_2
      DOUBLE PRECISION I_A_FREQ_ANA, I_A_FREQ_CAT_1, I_A_FREQ_CAT_2

C Time sensitivity factors
      DOUBLE PRECISION TNF_MAG_FACT, TNF_FREQ_FACT

C CONSTANT PARAMETERS:
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.D0, ONE = 1.D0

C Initialize variables
      TNF_MAG_ANA_RANGE = ZERO
      TNF_FREQ_ANA_RANGE = ZERO
      X_S_TNF_MAG_DEPt = ZERO
      X_S_TNF_FREQ_DEPt = ZERO

C---------------------------------------------------------------------------
C Magnitude <= 0.998 and Frequency <= 3.0 (Magnitude and Frequency anabolic)
C---------------------------------------------------------------------------

      IF (MAG .LE. 0.998D0 .AND. FREQ .LE. 3.0D0 .AND.
     1    FREQ .NE. ZERO) THEN

C Magnitude Anabolic
            TNF_MAG_ANA_RANGE = I_A_MAG_ANA * TNF_MAG_FACT

            IF (TNF_MAG_ANA_RANGE .GE. ABS(X_S_TNF_MAG)) THEN
                  TNF_MAG_ANA_RANGE = ABS(X_S_TNF_MAG)
                  X_S_TNF_MAG_COMPL = X_S_TNF_MAG_COMPL +
     1                                (I_A_MAG_CAT_1 * TNF_MAG_FACT)

                  IF (X_S_TNF_MAG_COMPL .GT. ONE) THEN
                        X_S_TNF_MAG_COMPL = ONE
                  END IF
            ELSE
                  X_S_TNF_MAG_COMPL = ZERO
            END IF

            X_S_TNF_MAG_DEPt = ABS(X_S_TNF_MAG + TNF_MAG_ANA_RANGE)

C Frequency Anabolic
            TNF_FREQ_ANA_RANGE = I_A_FREQ_ANA * TNF_FREQ_FACT

            IF (TNF_FREQ_ANA_RANGE .GE. ABS(X_S_TNF_FREQ)) THEN
                  TNF_FREQ_ANA_RANGE = ABS(X_S_TNF_FREQ)
                  X_S_TNF_FREQ_COMPL = X_S_TNF_FREQ_COMPL +
     1                                 (I_A_FREQ_CAT_1 * TNF_FREQ_FACT)

                  IF (X_S_TNF_FREQ_COMPL .GT. ONE) THEN
                        X_S_TNF_FREQ_COMPL = ONE
                  END IF
            ELSE
                  X_S_TNF_FREQ_COMPL = ZERO
            END IF

            X_S_TNF_FREQ_DEPt = ABS(X_S_TNF_FREQ + TNF_FREQ_ANA_RANGE)

C TNF mRNA Expression
            TNF = ROUND(
     1             (((ONE + AC_int_TNF_dir) / AC_int_TNF_dir) *
     2             (((aTHETA_TNF_LACT * X_S_TNF_LACT) +
     3             (aTHETA_TNF_GLUC * X_S_TNF_GLUC) +
     4             (aTHETA_TNF_MAG * X_S_TNF_MAG_COMPL) +
     5             (aTHETA_TNF_FREQ * X_S_TNF_FREQ_COMPL)) / (ONE +
     6             (aTHETA_TNF_LACT * X_S_TNF_LACT) +
     7             (aTHETA_TNF_GLUC * X_S_TNF_GLUC) +
     8             (aTHETA_TNF_MAG * X_S_TNF_MAG_COMPL) +
     9             (aTHETA_TNF_FREQ * X_S_TNF_FREQ_COMPL)))) * (ONE -
     1             (((iTHETA_TNF_MAG + iTHETA_TNF_FREQ) /
     2             (aTHETA_TNF_GLUC + aTHETA_TNF_LACT +
     3             iTHETA_TNF_MAG + iTHETA_TNF_FREQ)) * ((ONE +
     4             IC_int_TNF_dir) / IC_int_TNF_dir) *
     5             (((iTHETA_TNF_MAG * X_S_TNF_MAG_DEPt) +
     6             (iTHETA_TNF_FREQ * X_S_TNF_FREQ_DEPt)) / (ONE +
     7             (iTHETA_TNF_MAG * X_S_TNF_MAG_DEPt) +
     8             (iTHETA_TNF_FREQ * X_S_TNF_FREQ_DEPt))))), 4)

C---------------------------------------------------------------------------
C Magnitude <= 0.998 and Frequency == 0 (Magnitude anabolic, Frequency 0)
C---------------------------------------------------------------------------

      ELSE IF (MAG .LE. 0.998D0 .AND. FREQ .EQ. ZERO) THEN

C Magnitude Anabolic Range
            TNF_MAG_ANA_RANGE = I_A_MAG_ANA * TNF_MAG_FACT

            IF (TNF_MAG_ANA_RANGE .GE. ABS(X_S_TNF_MAG)) THEN
                  TNF_MAG_ANA_RANGE = ABS(X_S_TNF_MAG)
                  X_S_TNF_MAG_COMPL = X_S_TNF_MAG_COMPL +
     1                                (I_A_MAG_CAT_1 * TNF_MAG_FACT)

                  IF (X_S_TNF_MAG_COMPL .GT. ONE) THEN
                        X_S_TNF_MAG_COMPL = ONE
                  END IF
            ELSE
                  X_S_TNF_MAG_COMPL = ZERO
            END IF

            X_S_TNF_MAG_DEPt = ABS(X_S_TNF_MAG + TNF_MAG_ANA_RANGE)

C Frequency Anabolic Range
            TNF_FREQ_ANA_RANGE = FREQ_0_TIME *
     1                           DBLE(TICKS) * TNF_FREQ_FACT

            IF (TNF_FREQ_ANA_RANGE .GE. FREQ_0) THEN
                  TNF_FREQ_ANA_RANGE = FREQ_0
                  X_S_TNF_FREQ_COMPL = X_S_TNF_FREQ_COMPL +
     1                                 (FREQ_0_TIME * TNF_FREQ_FACT)

                  IF (X_S_TNF_FREQ_COMPL .GT. ONE) THEN
                        X_S_TNF_FREQ_COMPL = ONE
                  END IF
            ELSE
                  X_S_TNF_FREQ_COMPL = ZERO
            END IF

            X_S_TNF_FREQ_DEPt = FREQ_0 - TNF_FREQ_ANA_RANGE

C TNF mRNA Expression
            TNF = ROUND(
     1             (((ONE + AC_int_TNF_dir) / (AC_int_TNF_dir)) *
     2             (((aTHETA_TNF_LACT * X_S_TNF_LACT) +
     3             (aTHETA_TNF_GLUC * X_S_TNF_GLUC) +
     4             (aTHETA_TNF_MAG * X_S_TNF_MAG_COMPL) +
     5             (aTHETA_TNF_FREQ * X_S_TNF_FREQ_COMPL)) / (ONE +
     6             (aTHETA_TNF_LACT * X_S_TNF_LACT) +
     7             (aTHETA_TNF_GLUC * X_S_TNF_GLUC) +
     8             (aTHETA_TNF_MAG * X_S_TNF_MAG_COMPL) +
     9             (aTHETA_TNF_FREQ * X_S_TNF_FREQ_COMPL)))) * (ONE -
     1             (((iTHETA_TNF_MAG + iTHETA_TNF_FREQ) /
     2             (aTHETA_TNF_GLUC + aTHETA_TNF_LACT +
     3             iTHETA_TNF_MAG + iTHETA_TNF_FREQ)) * ((ONE +
     4             IC_int_TNF_dir) / (IC_int_TNF_dir)) *
     5             (((iTHETA_TNF_MAG * X_S_TNF_MAG_DEPt) +
     6             (iTHETA_TNF_FREQ * X_S_TNF_FREQ_DEPt)) / (ONE +
     7             (iTHETA_TNF_MAG * X_S_TNF_MAG_DEPt) +
     8             (iTHETA_TNF_FREQ * X_S_TNF_FREQ_DEPt))))), 4)

C---------------------------------------------------------------------------------
C Magnitude <= 0.998 and Frequency > 3.0 (Magnitude anabolic, Frequency catabolic)
C---------------------------------------------------------------------------------
C protease => mag inhibits within the first period of time, freq activates
C
      ELSE IF (MAG .LE. 0.998D0 .AND. FREQ .GT. 3.0D0) THEN

C Magnitude Anabolic Range
            TNF_MAG_ANA_RANGE = I_A_MAG_ANA * TNF_MAG_FACT

            IF (TNF_MAG_ANA_RANGE .GE. ABS(X_S_TNF_MAG)) THEN
                  TNF_MAG_ANA_RANGE = ABS(X_S_TNF_MAG)
                  X_S_TNF_MAG_COMPL = X_S_TNF_MAG_COMPL +
     1                                (I_A_MAG_CAT_1 * TNF_MAG_FACT)

                  IF (X_S_TNF_MAG_COMPL .GT. ONE) THEN
                        X_S_TNF_MAG_COMPL = ONE
                  END IF
            ELSE
                  X_S_TNF_MAG_COMPL = ZERO
            END IF

            X_S_TNF_MAG_DEPt = ABS(X_S_TNF_MAG + TNF_MAG_ANA_RANGE)

C Frequency Catabolic Range
            X_S_TNF_FREQ_DEPt = X_S_TNF_FREQ     ! no latency time of proinflammatory cytokines, value of generic function is taken

C TNF mRNA Expression
            TNF = ROUND(
     1             (((ONE + AC_int_TNF_dir) / (AC_int_TNF_dir)) *
     2             (((aTHETA_TNF_LACT * X_S_TNF_LACT) +
     3             (aTHETA_TNF_GLUC * X_S_TNF_GLUC) +
     4             (aTHETA_TNF_MAG * X_S_TNF_MAG_COMPL) +
     5             (aTHETA_TNF_FREQ * X_S_TNF_FREQ_DEPt)) / (ONE +
     6             (aTHETA_TNF_LACT * X_S_TNF_LACT) +
     7             (aTHETA_TNF_GLUC * X_S_TNF_GLUC) +
     8             (aTHETA_TNF_MAG * X_S_TNF_MAG_COMPL) +
     9             (aTHETA_TNF_FREQ * X_S_TNF_FREQ_DEPt)))) * (ONE -
     1             (((iTHETA_TNF_MAG) /
     2             (aTHETA_TNF_GLUC + aTHETA_TNF_LACT +
     3             iTHETA_TNF_MAG + aTHETA_TNF_FREQ)) * ((ONE +
     4             IC_int_TNF_dir) / (IC_int_TNF_dir)) *
     5             (((iTHETA_TNF_MAG * X_S_TNF_MAG_DEPt)) / (ONE +
     6             (iTHETA_TNF_MAG * X_S_TNF_MAG_DEPt))))), 4)

C--------------------------------------------------------------------------------
C Magnitude >= 1.0 and Frequency <= 3.0 (Magnitude catabolic, Frequency anabolic)
C--------------------------------------------------------------------------------
C Protease => mag activates, freq inhibits within the first period of time
C
      ELSE IF (MAG .GE. ONE .AND. FREQ .LE. 3.0D0 .AND.
     1         FREQ .NE. ZERO) THEN

C Magnitude Catabolic Range
            X_S_TNF_MAG_DEPt = X_S_TNF_MAG   ! no latency time of proinflammatory cytokines, value of generic function is taken

C Frequency Anabolic Range
            TNF_FREQ_ANA_RANGE = I_A_FREQ_ANA * TNF_FREQ_FACT

            IF (TNF_FREQ_ANA_RANGE .GE. ABS(X_S_TNF_FREQ)) THEN
                  TNF_FREQ_ANA_RANGE = ABS(X_S_TNF_FREQ)
                  X_S_TNF_FREQ_COMPL = X_S_TNF_FREQ_COMPL +
     1                                 (I_A_FREQ_CAT_1 * TNF_FREQ_FACT)

                  IF (X_S_TNF_FREQ_COMPL .GT. ONE) THEN
                        X_S_TNF_FREQ_COMPL = ONE
                  END IF
            ELSE
                  X_S_TNF_FREQ_COMPL = ZERO
            END IF

            X_S_TNF_FREQ_DEPt = ABS(X_S_TNF_FREQ + TNF_FREQ_ANA_RANGE)

C TNF mRNA Expression
            TNF = ROUND(
     1             (((ONE + AC_int_TNF_dir) / (AC_int_TNF_dir)) *
     2             (((aTHETA_TNF_LACT * X_S_TNF_LACT) +
     3             (aTHETA_TNF_GLUC * X_S_TNF_GLUC) +
     4             (aTHETA_TNF_MAG * X_S_TNF_MAG_DEPt) +
     5             (aTHETA_TNF_FREQ * X_S_TNF_FREQ_COMPL)) / (ONE +
     6             (aTHETA_TNF_LACT * X_S_TNF_LACT) +
     7             (aTHETA_TNF_GLUC * X_S_TNF_GLUC) +
     8             (aTHETA_TNF_MAG * X_S_TNF_MAG_DEPt) +
     9             (aTHETA_TNF_FREQ * X_S_TNF_FREQ_COMPL)))) * (ONE -
     1             (((iTHETA_TNF_FREQ) /
     2             (aTHETA_TNF_GLUC + aTHETA_TNF_LACT +
     3             aTHETA_TNF_MAG + iTHETA_TNF_FREQ)) * ((ONE +
     4             IC_int_TNF_dir) / (IC_int_TNF_dir)) *
     5             (((iTHETA_TNF_FREQ * X_S_TNF_FREQ_DEPt)) / (ONE +
     6             (iTHETA_TNF_FREQ * X_S_TNF_FREQ_DEPt))))), 4)


C---------------------------------------------------------------------------
C Magnitude >= 1.0 and Frequency == 0 (Magnitude catabolic, Frequency 0)
C---------------------------------------------------------------------------

      ELSE IF (MAG .GE. ONE .AND. FREQ .EQ. ZERO) THEN

C Magnitude Catabolic Range
            X_S_TNF_MAG_DEPt = X_S_TNF_MAG     ! no latency time of proinflammatory cytokines, value of generic function is taken

C Frequency Anabolic Range
            TNF_FREQ_ANA_RANGE = FREQ_0_TIME *
     1                           DBLE(TICKS) * TNF_FREQ_FACT

            IF (TNF_FREQ_ANA_RANGE .GE. FREQ_0) THEN
                  TNF_FREQ_ANA_RANGE = FREQ_0
                  X_S_TNF_FREQ_COMPL = X_S_TNF_FREQ_COMPL +
     1                                 (FREQ_0_TIME * TNF_FREQ_FACT)

                  IF (X_S_TNF_FREQ_COMPL .GT. ONE) THEN
                        X_S_TNF_FREQ_COMPL = ONE
                  END IF
            ELSE
                  X_S_TNF_FREQ_COMPL = ZERO
            END IF

            X_S_TNF_FREQ_DEPt = FREQ_0 - TNF_FREQ_ANA_RANGE

C TNF mRNA Expression
            TNF = ROUND(
     1             (((ONE + AC_int_TNF_dir) / (AC_int_TNF_dir)) *
     2             (((aTHETA_TNF_LACT * X_S_TNF_LACT) +
     3             (aTHETA_TNF_GLUC * X_S_TNF_GLUC) +
     4             (aTHETA_TNF_MAG * X_S_TNF_MAG_DEPt) +
     5             (aTHETA_TNF_FREQ * X_S_TNF_FREQ_COMPL)) / (ONE +
     6             (aTHETA_TNF_LACT * X_S_TNF_LACT) +
     7             (aTHETA_TNF_GLUC * X_S_TNF_GLUC) +
     8             (aTHETA_TNF_MAG * X_S_TNF_MAG_DEPt) +
     9             (aTHETA_TNF_FREQ * X_S_TNF_FREQ_COMPL)))) * (ONE -
     1             (((iTHETA_TNF_FREQ) /
     2             (aTHETA_TNF_GLUC + aTHETA_TNF_LACT +
     3             aTHETA_TNF_MAG + iTHETA_TNF_FREQ)) * ((ONE +
     4             IC_int_TNF_dir) / (IC_int_TNF_dir)) *
     5             (((iTHETA_TNF_FREQ * X_S_TNF_FREQ_DEPt)) / (ONE +
     6             (iTHETA_TNF_FREQ * X_S_TNF_FREQ_DEPt))))), 4)

C---------------------------------------------------------------------------
C Magnitude >= 1.0 and Frequency > 3.0 (Magnitude and Frequency catabolic)
C---------------------------------------------------------------------------
C protease => mag and freq activate
C
      ELSE IF (MAG .GE. ONE .AND. FREQ .GT. 3.0D0) THEN

C Magnitude Catabolic Range
            X_S_TNF_MAG_DEPt = X_S_TNF_MAG    ! no latency time of proinflammatory cytokines, value of generic function is taken

C Frequency Catabolic Range
            X_S_TNF_FREQ_DEPt = X_S_TNF_FREQ  ! no latency time of proinflammatory cytokines, value of generic function is taken

C TNF mRNA Expression
            TNF = ROUND(
     1            (((ONE + AC_int_TNF_dir) / (AC_int_TNF_dir)) *
     2            (((aTHETA_TNF_LACT * X_S_TNF_LACT) +
     3            (aTHETA_TNF_GLUC * X_S_TNF_GLUC) +
     4            (aTHETA_TNF_MAG * X_S_TNF_MAG_DEPt) +
     5            (aTHETA_TNF_FREQ * X_S_TNF_FREQ_DEPt)) / (ONE +
     6            (aTHETA_TNF_LACT * X_S_TNF_LACT) +
     7            (aTHETA_TNF_GLUC * X_S_TNF_GLUC) +
     8            (aTHETA_TNF_MAG * X_S_TNF_MAG_DEPt) +
     9            (aTHETA_TNF_FREQ * X_S_TNF_FREQ_DEPt))))
     1            , 4 )

      END IF

C--------------------------------------------------------------------------
      RETURN
      END SUBROUTINE PNeq_TNF

C---------------------------------------------------------------------------------------------------------------------------------------------
C Evaluate the Aggrecan expression
C ------------------------------------------------------------

      SUBROUTINE PNeq_AGG(AGG_IMN, AGG_IL1B, AGG_TNF, AGG_IL1B_TNF,
     1           TICKS,
     2           IL1B_PROT, TNF_PROT,
     3           GLUC_int, pH_val, MAG, FREQ,
     4           FREQ_0, FREQ_0_time,
     5           MAG_0, MAG_0_time,
     6           aTHETA_AGG_GLUC, aTHETA_AGG_LACT,
     7           iTHETA_AGG_IL1B, iTHETA_AGG_TNF,
     8           aTHETA_AGG_MAG, iTHETA_AGG_MAG,
     9           aTHETA_AGG_FREQ, iTHETA_AGG_FREQ,
     1           AC_int_IL1B_dir, AC_int_TNF_dir, IC_int_AGG_dir,
     2           X_S_AGG_GLUC, X_S_AGG_LACT,
     3           X_S_AGG_MAG, X_S_AGG_FREQ,
     4           A_I_MAG_ANA, A_I_MAG_CAT_1, A_I_MAG_CAT_2,
     5           I_A_MAG_ANA, I_A_MAG_CAT_1, I_A_MAG_CAT_2,
     6           A_I_FREQ_ANA, A_I_FREQ_CAT_1, A_I_FREQ_CAT_2,
     7           I_A_FREQ_ANA, I_A_FREQ_CAT_1, I_A_FREQ_CAT_2,
     8           AGG_MAG_FACT, AGG_MAG_FACT_2,
     9           AGG_FREQ_FACT, AGG_FREQ_FACT_2,
     1           X_S_AGG_MAG_COMPL, X_S_AGG_FREQ_COMPL,
     2           aTHETA_SUM)

C VARIABLES:
C X_S: SENSITIVITY OF A CELL ACTIVITY (CA) TO A STIMULUS CONCENTRATION
C THETA: SENSITIVITY OF A CELL ACTIVITY (CA) TO A STIMULUS TYPE
C OUTPUTS:
C AGG_IMN: AGGRECAN EXPRESSION FOR IMMUNONEGATIVE CELLS
C AGG_IL1B: CONTRIBUTION FROM IL1B
C AGG_TNF: CONTRIBUTION FROM TNF
C AGG_IL1B_TNF: COMBINED CONTRIBUTION OF IL1B AND TNF

C ********************************************************************
C
      IMPLICIT NONE
C
C FUNCTION
      DOUBLE PRECISION ROUND

C INTEGER VARIABLES
      INTEGER TICKS

C SPECIAL CASE FOR FREQ = 0
      DOUBLE PRECISION FREQ_0, FREQ_0_time
      DOUBLE PRECISION MAG_0, MAG_0_time

C Input and Output Variables
      DOUBLE PRECISION AGG_IMN, AGG_IL1B, AGG_TNF, AGG_IL1B_TNF
      DOUBLE PRECISION IL1B_PROT, TNF_PROT
      DOUBLE PRECISION GLUC_int, pH_val, MAG, FREQ

C Constants for Sensitivity Factors
      DOUBLE PRECISION aTHETA_AGG_GLUC, aTHETA_AGG_LACT
      DOUBLE PRECISION iTHETA_AGG_IL1B, iTHETA_AGG_TNF
      DOUBLE PRECISION aTHETA_AGG_MAG, iTHETA_AGG_MAG
      DOUBLE PRECISION aTHETA_AGG_FREQ, iTHETA_AGG_FREQ
      DOUBLE PRECISION AC_int_IL1B_dir, AC_int_TNF_dir
      DOUBLE PRECISION IC_int_AGG_dir
      DOUBLE PRECISION aTHETA_SUM

C Intermediate Variables
      DOUBLE PRECISION X_S_AGG_GLUC, X_S_AGG_LACT
      DOUBLE PRECISION X_S_AGG_MAG, X_S_AGG_FREQ
      DOUBLE PRECISION AGG_MAG_FACT, AGG_MAG_FACT_2
      DOUBLE PRECISION AGG_FREQ_FACT, AGG_FREQ_FACT_2

C Completive factors
      DOUBLE PRECISION X_S_AGG_MAG_COMPL
      DOUBLE PRECISION X_S_AGG_MAG_DEPt
      DOUBLE PRECISION X_S_AGG_FREQ_COMPL
      DOUBLE PRECISION X_S_AGG_FREQ_DEPt

C Time dependent funtions
      DOUBLE PRECISION A_I_MAG_ANA, A_I_MAG_CAT_1, A_I_MAG_CAT_2
      DOUBLE PRECISION I_A_MAG_ANA, I_A_MAG_CAT_1, I_A_MAG_CAT_2
      DOUBLE PRECISION A_I_FREQ_ANA, A_I_FREQ_CAT_1, A_I_FREQ_CAT_2
      DOUBLE PRECISION I_A_FREQ_ANA, I_A_FREQ_CAT_1, I_A_FREQ_CAT_2

C Time sensitivity range
      DOUBLE PRECISION AGG_MAG_ANA_RANGE
      DOUBLE PRECISION AGG_FREQ_ANA_RANGE
      DOUBLE PRECISION AGG_MAG_CAT_RANGE
      DOUBLE PRECISION AGG_FREQ_CAT_RANGE

C CONSTANT PARAMETERS
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.D0, ONE = 1.D0

C---------------------------------------------------------------------------
C Initialize Variables
      AGG_MAG_ANA_RANGE = ZERO
      AGG_FREQ_ANA_RANGE = ZERO
      AGG_MAG_CAT_RANGE = ZERO
      AGG_FREQ_CAT_RANGE = ZERO
      X_S_AGG_MAG_DEPt = ZERO
      X_S_AGG_FREQ_DEPt = ZERO

C---------------------------------------------------------------------------
C Magnitude <= 0.998 and Frequency <= 3.0 (Magnitude and Frequency anabolic)
C---------------------------------------------------------------------------

      IF (MAG .LE. 0.998D0 .AND. FREQ .LE. 3.0D0 .AND.
     1    FREQ .NE. ZERO) THEN

C Magnitude Anabolic
            AGG_MAG_ANA_RANGE = A_I_MAG_ANA * AGG_MAG_FACT

            IF (AGG_MAG_ANA_RANGE .GE. ABS(X_S_AGG_MAG)) THEN
                  AGG_MAG_ANA_RANGE = ABS(X_S_AGG_MAG)
                  X_S_AGG_MAG_COMPL = X_S_AGG_MAG_COMPL +
     1                                  (A_I_MAG_CAT_1 * AGG_MAG_FACT)

                  IF (X_S_AGG_MAG_COMPL .GT. ONE) THEN
                        X_S_AGG_MAG_COMPL = ONE
                  END IF
            ELSE
                  X_S_AGG_MAG_COMPL = ZERO
            END IF

            X_S_AGG_MAG_DEPt = X_S_AGG_MAG - AGG_MAG_ANA_RANGE  ! becomes 0 if no anabolism due to time-sensitivity

C Frequency Anabolic
            AGG_FREQ_ANA_RANGE = A_I_FREQ_ANA * AGG_FREQ_FACT

            IF (AGG_FREQ_ANA_RANGE .GE. ABS(X_S_AGG_FREQ)) THEN
                  AGG_FREQ_ANA_RANGE = ABS(X_S_AGG_FREQ)
                  X_S_AGG_FREQ_COMPL = X_S_AGG_FREQ_COMPL +
     1                                  (A_I_FREQ_CAT_1 * AGG_FREQ_FACT)

                  IF (X_S_AGG_FREQ_COMPL .GT. ONE) THEN
                        X_S_AGG_FREQ_COMPL = ONE
                  END IF
            ELSE
                  X_S_AGG_FREQ_COMPL = ZERO
            END IF

            X_S_AGG_FREQ_DEPt = X_S_AGG_FREQ - AGG_FREQ_ANA_RANGE

C Aggrecan Expression for Immunonegative Cells
            AGG_IMN = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_AGG_LACT * X_S_AGG_LACT) +
     2             (aTHETA_AGG_GLUC * X_S_AGG_GLUC) +
     3             (aTHETA_AGG_MAG * X_S_AGG_MAG_DEPt) +
     4             (aTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt)) / (ONE +
     5             (aTHETA_AGG_LACT * X_S_AGG_LACT) +
     6             (aTHETA_AGG_GLUC * X_S_AGG_GLUC) +
     7             (aTHETA_AGG_MAG * X_S_AGG_MAG_DEPt) +
     8             (aTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt)))) * (ONE -
     9             (((iTHETA_AGG_MAG + iTHETA_AGG_FREQ) /
     1             (aTHETA_AGG_GLUC + aTHETA_AGG_LACT +
     2             iTHETA_AGG_MAG + iTHETA_AGG_FREQ)) *
     3             ((ONE + IC_int_AGG_dir) / (IC_int_AGG_dir)) *
     4             (((iTHETA_AGG_MAG * X_S_AGG_MAG_COMPL) +
     5             (iTHETA_AGG_FREQ * X_S_AGG_FREQ_COMPL)) / (ONE +
     6             (iTHETA_AGG_MAG * X_S_AGG_MAG_COMPL) +
     7             (iTHETA_AGG_FREQ * X_S_AGG_FREQ_COMPL)))))
     8             , 4)
            
C Aggrecan Expression for IL1B Immunopositive Cells
      AGG_IL1B = ROUND(
     1             (((ONE + aTHETA_SUM) / aTHETA_SUM) *
     2             (((aTHETA_AGG_LACT * X_S_AGG_LACT) +
     3             (aTHETA_AGG_GLUC * X_S_AGG_GLUC) +
     4             (aTHETA_AGG_MAG * X_S_AGG_MAG_DEPt) +
     5             (aTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt)) / (ONE +
     6             (aTHETA_AGG_LACT * X_S_AGG_LACT) +
     7             (aTHETA_AGG_GLUC * X_S_AGG_GLUC) +
     8             (aTHETA_AGG_MAG * X_S_AGG_MAG_DEPt) +
     9             (aTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt)))) * (ONE -
     1             (((iTHETA_AGG_MAG + iTHETA_AGG_FREQ +
     2             iTHETA_AGG_IL1B) /
     3             (aTHETA_AGG_GLUC + aTHETA_AGG_LACT +
     4             iTHETA_AGG_MAG + iTHETA_AGG_FREQ +
     5             iTHETA_AGG_IL1B)) *
     6             ((ONE + IC_int_AGG_dir) / IC_int_AGG_dir) *
     7             (((iTHETA_AGG_MAG * X_S_AGG_MAG_COMPL) +
     8             (iTHETA_AGG_FREQ * X_S_AGG_FREQ_COMPL) +
     9             (iTHETA_AGG_IL1B * IL1B_PROT)) / (ONE +
     1             (iTHETA_AGG_MAG * X_S_AGG_MAG_COMPL) +
     2             (iTHETA_AGG_FREQ * X_S_AGG_FREQ_COMPL) +
     3             (iTHETA_AGG_IL1B * IL1B_PROT)))))
     4             , 4)

C Aggrecan Expression for TNF Immunopositive Cells
      AGG_TNF = ROUND(
     1             (((ONE + aTHETA_SUM) / aTHETA_SUM) *
     2             (((aTHETA_AGG_LACT * X_S_AGG_LACT) +
     3             (aTHETA_AGG_GLUC * X_S_AGG_GLUC) +
     4             (aTHETA_AGG_MAG * X_S_AGG_MAG_DEPt) +
     5             (aTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt)) / (ONE +
     6             (aTHETA_AGG_LACT * X_S_AGG_LACT) +
     7             (aTHETA_AGG_GLUC * X_S_AGG_GLUC) +
     8             (aTHETA_AGG_MAG * X_S_AGG_MAG_DEPt) +
     9             (aTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt)))) * (ONE -
     1             (((iTHETA_AGG_MAG + iTHETA_AGG_FREQ +
     2             iTHETA_AGG_TNF) /
     3             (aTHETA_AGG_GLUC + aTHETA_AGG_LACT +
     4             iTHETA_AGG_MAG + iTHETA_AGG_FREQ +
     5             iTHETA_AGG_TNF)) *
     6             ((ONE + IC_int_AGG_dir) / IC_int_AGG_dir) *
     7             (((iTHETA_AGG_MAG * X_S_AGG_MAG_COMPL) +
     8             (iTHETA_AGG_FREQ * X_S_AGG_FREQ_COMPL) +
     9             (iTHETA_AGG_TNF * TNF_PROT)) / (ONE +
     1             (iTHETA_AGG_MAG * X_S_AGG_MAG_COMPL) +
     2             (iTHETA_AGG_FREQ * X_S_AGG_FREQ_COMPL) +
     3             (iTHETA_AGG_TNF * TNF_PROT)))))
     4             , 4)

C Aggrecan Expression for IL1B & TNF Immunopositive Cells
      AGG_IL1B_TNF = ROUND(
     1             (((ONE + aTHETA_SUM) / aTHETA_SUM) *
     2             (((aTHETA_AGG_LACT * X_S_AGG_LACT) +
     3             (aTHETA_AGG_GLUC * X_S_AGG_GLUC) +
     4             (aTHETA_AGG_MAG * X_S_AGG_MAG_DEPt) +
     5             (aTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt)) / (ONE +
     6             (aTHETA_AGG_LACT * X_S_AGG_LACT) +
     7             (aTHETA_AGG_GLUC * X_S_AGG_GLUC) +
     8             (aTHETA_AGG_MAG * X_S_AGG_MAG_DEPt) +
     9             (aTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt)))) * (ONE -
     1             (((iTHETA_AGG_MAG + iTHETA_AGG_FREQ +
     2             iTHETA_AGG_IL1B + iTHETA_AGG_TNF) /
     3             (aTHETA_AGG_GLUC + aTHETA_AGG_LACT +
     4             iTHETA_AGG_MAG + iTHETA_AGG_FREQ +
     5             iTHETA_AGG_IL1B + iTHETA_AGG_TNF)) *
     6             ((ONE + IC_int_AGG_dir) / IC_int_AGG_dir) *
     7             (((iTHETA_AGG_MAG * X_S_AGG_MAG_COMPL) +
     8             (iTHETA_AGG_FREQ * X_S_AGG_FREQ_COMPL) +
     9             (iTHETA_AGG_IL1B * IL1B_PROT) +
     1             (iTHETA_AGG_TNF * TNF_PROT)) / (ONE +
     2             (iTHETA_AGG_MAG * X_S_AGG_MAG_COMPL) +
     3             (iTHETA_AGG_FREQ * X_S_AGG_FREQ_COMPL) +
     4             (iTHETA_AGG_IL1B * IL1B_PROT) +
     5             (iTHETA_AGG_TNF * TNF_PROT)))))
     6             , 4)

C---------------------------------------------------------------------------
C Magnitude <= 0.998 and Frequency == 0 (Magnitude anabolic, Frequency 0)
C---------------------------------------------------------------------------

      ELSE IF (MAG .LE. 0.998D0 .AND. FREQ .EQ. ZERO) THEN

C Magnitude Anabolic
            AGG_MAG_ANA_RANGE = A_I_MAG_ANA * AGG_MAG_FACT

            IF (AGG_MAG_ANA_RANGE .GE. X_S_AGG_MAG) THEN
                  AGG_MAG_ANA_RANGE = X_S_AGG_MAG
                  X_S_AGG_MAG_COMPL = X_S_AGG_MAG_COMPL +
     1                                  (A_I_MAG_CAT_1 * AGG_MAG_FACT)

                  IF (X_S_AGG_MAG_COMPL .GT. ONE) THEN
                        X_S_AGG_MAG_COMPL = ONE
                  END IF
            ELSE
                  X_S_AGG_MAG_COMPL = ZERO
            END IF

            X_S_AGG_MAG_DEPt = X_S_AGG_MAG - AGG_MAG_ANA_RANGE  ! becomes 0 if no anabolism due to time-sensitivity

C Frequency Anabolic (Frequency == 0 Case)
            AGG_FREQ_ANA_RANGE = FREQ_0_TIME *
     1                           DBLE(TICKS) * AGG_FREQ_FACT

            IF (AGG_FREQ_ANA_RANGE .GE. FREQ_0) THEN
                  AGG_FREQ_ANA_RANGE = FREQ_0
                  X_S_AGG_FREQ_COMPL = X_S_AGG_FREQ_COMPL +
     1                                  (FREQ_0_TIME * AGG_FREQ_FACT)

                  IF (X_S_AGG_FREQ_COMPL .GT. ONE) THEN
                        X_S_AGG_FREQ_COMPL = ONE  ! Maximum value is 1
                  END IF
            ELSE
                  X_S_AGG_FREQ_COMPL = ZERO
            END IF

            X_S_AGG_FREQ_DEPt = FREQ_0 - AGG_FREQ_ANA_RANGE

C Aggrecan Expression for Immunonegative Cells
            AGG_IMN = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1                (((aTHETA_AGG_LACT * X_S_AGG_LACT) +
     2                (aTHETA_AGG_GLUC * X_S_AGG_GLUC) +
     3                (aTHETA_AGG_MAG * X_S_AGG_MAG_DEPt) +
     4                (aTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt)) / (ONE +
     5                (aTHETA_AGG_LACT * X_S_AGG_LACT) +
     6                (aTHETA_AGG_GLUC * X_S_AGG_GLUC) +
     7                (aTHETA_AGG_MAG * X_S_AGG_MAG_DEPt) +
     8                (aTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt)))) * (ONE -
     9                (((iTHETA_AGG_MAG + iTHETA_AGG_FREQ) /
     1                (aTHETA_AGG_GLUC + aTHETA_AGG_LACT +
     2                iTHETA_AGG_MAG + iTHETA_AGG_FREQ)) *
     3                ((ONE + IC_int_AGG_dir) / IC_int_AGG_dir) *
     4                (((iTHETA_AGG_MAG * X_S_AGG_MAG_COMPL) +
     5                (iTHETA_AGG_FREQ * X_S_AGG_FREQ_COMPL)) / (ONE +
     6                (iTHETA_AGG_MAG * X_S_AGG_MAG_COMPL) +
     7                (iTHETA_AGG_FREQ * X_S_AGG_FREQ_COMPL)))))
     8                , 4)


C Aggrecan Expression for IL1B Immunopositive Cells
            AGG_IL1B = ROUND((((ONE + aTHETA_SUM) / aTHETA_SUM) *
     1                (((aTHETA_AGG_LACT * X_S_AGG_LACT) +
     2                (aTHETA_AGG_GLUC * X_S_AGG_GLUC) +
     3                (aTHETA_AGG_MAG * X_S_AGG_MAG_DEPt) +
     4                (aTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt)) / (ONE +
     5                (aTHETA_AGG_LACT * X_S_AGG_LACT) +
     6                (aTHETA_AGG_GLUC * X_S_AGG_GLUC) +
     7                (aTHETA_AGG_MAG * X_S_AGG_MAG_DEPt) +
     8                (aTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt)))) * (ONE -
     9                (((iTHETA_AGG_MAG + iTHETA_AGG_FREQ +
     1                iTHETA_AGG_IL1B) /
     2                (aTHETA_AGG_GLUC + aTHETA_AGG_LACT +
     3                iTHETA_AGG_MAG + iTHETA_AGG_FREQ +
     4                iTHETA_AGG_IL1B)) *
     5                ((ONE + IC_int_AGG_dir) / IC_int_AGG_dir) *
     6                (((iTHETA_AGG_MAG * X_S_AGG_MAG_COMPL) +
     7                (iTHETA_AGG_FREQ * X_S_AGG_FREQ_COMPL) +
     8                (iTHETA_AGG_IL1B * IL1B_PROT)) / (ONE +
     9                (iTHETA_AGG_MAG * X_S_AGG_MAG_COMPL) +
     1                (iTHETA_AGG_FREQ * X_S_AGG_FREQ_COMPL) +
     2                (iTHETA_AGG_IL1B * IL1B_PROT)))))
     3                , 4)

C Aggrecan Expression for TNF Immunopositive Cells
            AGG_TNF = ROUND((((ONE + aTHETA_SUM) / aTHETA_SUM) *
     1                (((aTHETA_AGG_LACT * X_S_AGG_LACT) +
     2                (aTHETA_AGG_GLUC * X_S_AGG_GLUC) +
     3                (aTHETA_AGG_MAG * X_S_AGG_MAG_DEPt) +
     4                (aTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt)) / (ONE +
     5                (aTHETA_AGG_LACT * X_S_AGG_LACT) +
     6                (aTHETA_AGG_GLUC * X_S_AGG_GLUC) +
     7                (aTHETA_AGG_MAG * X_S_AGG_MAG_DEPt) +
     8                (aTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt)))) * (ONE -
     9                (((iTHETA_AGG_MAG + iTHETA_AGG_FREQ +
     1                iTHETA_AGG_TNF) /
     2                (aTHETA_AGG_GLUC + aTHETA_AGG_LACT +
     3                iTHETA_AGG_MAG + iTHETA_AGG_FREQ +
     4                iTHETA_AGG_TNF)) *
     5                ((ONE + IC_int_AGG_dir) / IC_int_AGG_dir) *
     6                (((iTHETA_AGG_MAG * X_S_AGG_MAG_COMPL) +
     7                (iTHETA_AGG_FREQ * X_S_AGG_FREQ_COMPL) +
     8                (iTHETA_AGG_TNF * TNF_PROT)) / (ONE +
     9                (iTHETA_AGG_MAG * X_S_AGG_MAG_COMPL) +
     1                (iTHETA_AGG_FREQ * X_S_AGG_FREQ_COMPL) +
     2                (iTHETA_AGG_TNF * TNF_PROT)))))
     3                , 4)

C Aggrecan Expression for IL1B & TNF Immunopositive Cells
            AGG_IL1B_TNF = ROUND((((ONE + aTHETA_SUM) / aTHETA_SUM) *
     1                (((aTHETA_AGG_LACT * X_S_AGG_LACT) +
     2                (aTHETA_AGG_GLUC * X_S_AGG_GLUC) +
     3                (aTHETA_AGG_MAG * X_S_AGG_MAG_DEPt) +
     4                (aTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt)) / (ONE +
     5                (aTHETA_AGG_LACT * X_S_AGG_LACT) +
     6                (aTHETA_AGG_GLUC * X_S_AGG_GLUC) +
     7                (aTHETA_AGG_MAG * X_S_AGG_MAG_DEPt) +
     8                (aTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt)))) * (ONE -
     9                (((iTHETA_AGG_MAG + iTHETA_AGG_FREQ +
     1                iTHETA_AGG_IL1B + iTHETA_AGG_TNF) /
     2                (aTHETA_AGG_GLUC + aTHETA_AGG_LACT +
     3                iTHETA_AGG_MAG + iTHETA_AGG_FREQ +
     4                iTHETA_AGG_IL1B + iTHETA_AGG_TNF)) *
     5                ((ONE + IC_int_AGG_dir) / IC_int_AGG_dir) *
     6                (((iTHETA_AGG_MAG * X_S_AGG_MAG_COMPL) +
     7                (iTHETA_AGG_FREQ * X_S_AGG_FREQ_COMPL)
     8                + (iTHETA_AGG_IL1B * IL1B_PROT) + 
     9                iTHETA_AGG_TNF * TNF_PROT)) / (ONE +
     1                (iTHETA_AGG_MAG * X_S_AGG_MAG_COMPL) +
     2                (iTHETA_AGG_FREQ * X_S_AGG_FREQ_COMPL) +
     3                (iTHETA_AGG_IL1B * IL1B_PROT) +
     4                (iTHETA_AGG_TNF * TNF_PROT))))
     5                , 4)

C---------------------------------------------------------------------------
C Magnitude <= 0.998 and Frequency > 3 (Magnitude anabolic, Frequency catabolic)
C---------------------------------------------------------------------------

      ELSE IF (MAG .LE. 0.998D0 .AND. FREQ .GT. 3.0D0) THEN

C Magnitude Anabolic
            AGG_MAG_ANA_RANGE = A_I_MAG_ANA * AGG_MAG_FACT

            IF (AGG_MAG_ANA_RANGE .GE. X_S_AGG_MAG) THEN
                  AGG_MAG_ANA_RANGE = X_S_AGG_MAG
                  X_S_AGG_MAG_COMPL = X_S_AGG_MAG_COMPL +
     1                                  (A_I_MAG_CAT_1 * AGG_MAG_FACT)

                  IF (X_S_AGG_MAG_COMPL .GT. ONE) THEN
                        X_S_AGG_MAG_COMPL = ONE
                  END IF
            ELSE
                  X_S_AGG_MAG_COMPL = ZERO
            END IF

            X_S_AGG_MAG_DEPt = X_S_AGG_MAG - AGG_MAG_ANA_RANGE

C Frequency Catabolic
            AGG_FREQ_CAT_RANGE = A_I_FREQ_CAT_2 * AGG_FREQ_FACT_2

            IF (AGG_FREQ_CAT_RANGE .GE. ONE) THEN
                  AGG_FREQ_CAT_RANGE = ONE
            END IF

            X_S_AGG_FREQ_DEPt = ABS(ZERO - AGG_FREQ_CAT_RANGE)

C Aggrecan Expression for Immunonegative Cells
            AGG_IMN = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_AGG_LACT * X_S_AGG_LACT) +
     2             (aTHETA_AGG_GLUC * X_S_AGG_GLUC) +
     3             (aTHETA_AGG_MAG * X_S_AGG_MAG_DEPt)) / (ONE +
     4             (aTHETA_AGG_LACT * X_S_AGG_LACT) +
     5             (aTHETA_AGG_GLUC * X_S_AGG_GLUC) +
     6             (aTHETA_AGG_MAG * X_S_AGG_MAG_DEPt)))) * (ONE -
     7             (((iTHETA_AGG_MAG + iTHETA_AGG_FREQ) /
     8             (aTHETA_AGG_GLUC + aTHETA_AGG_LACT +
     9             iTHETA_AGG_MAG + iTHETA_AGG_FREQ)) *
     1             ((ONE + IC_int_AGG_dir) / (IC_int_AGG_dir)) *
     2             (((iTHETA_AGG_MAG * X_S_AGG_MAG_COMPL) +
     3             (iTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt)) / (ONE +
     4             (iTHETA_AGG_MAG * X_S_AGG_MAG_COMPL) +
     5             (iTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt)))))
     6             , 4)

C Aggrecan Expression for IL1B Immunopositive Cells
      AGG_IL1B = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_AGG_LACT * X_S_AGG_LACT) +
     2             (aTHETA_AGG_GLUC * X_S_AGG_GLUC) +
     3             (aTHETA_AGG_MAG * X_S_AGG_MAG_DEPt)) / (ONE +
     4             (aTHETA_AGG_LACT * X_S_AGG_LACT) +
     5             (aTHETA_AGG_GLUC * X_S_AGG_GLUC) +
     6             (aTHETA_AGG_MAG * X_S_AGG_MAG_DEPt))) ) * (ONE -
     7             (((iTHETA_AGG_MAG + iTHETA_AGG_FREQ +
     8             iTHETA_AGG_IL1B) /
     9             (aTHETA_AGG_GLUC + aTHETA_AGG_LACT +
     1             iTHETA_AGG_MAG + iTHETA_AGG_FREQ +
     2             iTHETA_AGG_IL1B)) *
     3             ((ONE+ IC_int_AGG_dir) / (IC_int_AGG_dir)) *
     4             (((iTHETA_AGG_MAG * X_S_AGG_MAG_COMPL) +
     5             (iTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt) +
     6             (iTHETA_AGG_IL1B * IL1B_PROT)) / (ONE +
     7             (iTHETA_AGG_MAG * X_S_AGG_MAG_COMPL) +
     8             (iTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt) +
     9             (iTHETA_AGG_IL1B * IL1B_PROT)))))
     1             , 4)

C Aggrecan Expression for TNF Immunopositive Cells
      AGG_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_AGG_LACT * X_S_AGG_LACT) +
     2             (aTHETA_AGG_GLUC * X_S_AGG_GLUC) +
     3             (aTHETA_AGG_MAG * X_S_AGG_MAG_DEPt)) / (ONE +
     4             (aTHETA_AGG_LACT * X_S_AGG_LACT) +
     5             (aTHETA_AGG_GLUC * X_S_AGG_GLUC) +
     6             (aTHETA_AGG_MAG * X_S_AGG_MAG_DEPt))) ) * (ONE -
     7             (((iTHETA_AGG_MAG + iTHETA_AGG_FREQ +
     8             iTHETA_AGG_TNF) /
     9             (aTHETA_AGG_GLUC + aTHETA_AGG_LACT +
     1             iTHETA_AGG_MAG + iTHETA_AGG_FREQ +
     2             iTHETA_AGG_TNF)) *
     3             ((ONE + IC_int_AGG_dir) / (IC_int_AGG_dir)) *
     4             (((iTHETA_AGG_MAG * X_S_AGG_MAG_COMPL) +
     5             (iTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt) +
     6             (iTHETA_AGG_TNF * TNF_PROT)) / (ONE +
     7             (iTHETA_AGG_MAG * X_S_AGG_MAG_COMPL) +
     8             (iTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt) +
     9             (iTHETA_AGG_TNF * TNF_PROT)))))
     1             , 4)

C Aggrecan Expression for IL1B & TNF Immunopositive Cells
      AGG_IL1B_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_AGG_LACT * X_S_AGG_LACT) +
     2             (aTHETA_AGG_GLUC * X_S_AGG_GLUC) +
     3             (aTHETA_AGG_MAG * X_S_AGG_MAG_DEPt)) / (ONE +
     4             (aTHETA_AGG_LACT * X_S_AGG_LACT) +
     5             (aTHETA_AGG_GLUC * X_S_AGG_GLUC) +
     6             (aTHETA_AGG_MAG * X_S_AGG_MAG_DEPt))) ) * (ONE -
     7             (((iTHETA_AGG_MAG + iTHETA_AGG_FREQ +
     8             iTHETA_AGG_IL1B + iTHETA_AGG_TNF) /
     9             (aTHETA_AGG_GLUC + aTHETA_AGG_LACT +
     1             iTHETA_AGG_MAG + iTHETA_AGG_FREQ +
     2             iTHETA_AGG_IL1B + iTHETA_AGG_TNF)) *
     3             ((ONE + IC_int_AGG_dir) / (IC_int_AGG_dir)) *
     4             (((iTHETA_AGG_MAG * X_S_AGG_MAG_COMPL) +
     5             (iTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt) +
     6             (iTHETA_AGG_IL1B * IL1B_PROT) +
     7             (iTHETA_AGG_TNF * TNF_PROT)) / (ONE +
     8             (iTHETA_AGG_MAG * X_S_AGG_MAG_COMPL) +
     9             (iTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt) +
     1             (iTHETA_AGG_IL1B * IL1B_PROT) +
     2             (iTHETA_AGG_TNF * TNF_PROT)))))
     3             , 4)


C---------------------------------------------------------------------------
C Magnitude >= 1 and Frequency <= 3 (Magnitude catabolic, Frequency anabolic)
C---------------------------------------------------------------------------

C Magnitude >= 1 and Frequency <= 3 (Magnitude catabolic, Frequency anabolic)
      ELSE IF (MAG .GE. ONE .AND. FREQ .LE. 3.D0 .AND.
     1         FREQ .NE. ZERO) THEN

C Magnitude Catabolic
            AGG_MAG_CAT_RANGE = A_I_MAG_CAT_2 * AGG_MAG_FACT_2

            IF (AGG_MAG_CAT_RANGE .GT. ONE) THEN
                  AGG_MAG_CAT_RANGE = ONE
            END IF

            X_S_AGG_MAG_DEPt = ABS(ZERO - AGG_MAG_CAT_RANGE)  ! Must be positive for PN-Equation

C Frequency Anabolic
            AGG_FREQ_ANA_RANGE = A_I_FREQ_ANA * AGG_FREQ_FACT

            IF (AGG_FREQ_ANA_RANGE .GE. X_S_AGG_FREQ) THEN
                  AGG_FREQ_ANA_RANGE = X_S_AGG_FREQ
                  X_S_AGG_FREQ_COMPL = X_S_AGG_FREQ_COMPL +
     1                                  (A_I_FREQ_CAT_1 * AGG_FREQ_FACT)

                  IF (X_S_AGG_FREQ_COMPL .GT. ONE) THEN
                        X_S_AGG_FREQ_COMPL = ONE
                  END IF
            ELSE
                  X_S_AGG_FREQ_COMPL = ZERO
            END IF

            X_S_AGG_FREQ_DEPt = X_S_AGG_FREQ - AGG_FREQ_ANA_RANGE  ! Becomes 0 if no anabolism due to time-sensitivity

C Aggrecan Expression for Immunonegative Cells
      AGG_IMN = ROUND((((ONE + aTHETA_SUM) / aTHETA_SUM) *
     1             (((aTHETA_AGG_LACT * X_S_AGG_LACT) +
     2             (aTHETA_AGG_GLUC * X_S_AGG_GLUC) +
     3             (aTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt))) / (ONE +
     4             (aTHETA_AGG_LACT * X_S_AGG_LACT) +
     5             (aTHETA_AGG_GLUC * X_S_AGG_GLUC) +
     6             (aTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt))) * (ONE -
     7             (((iTHETA_AGG_MAG + iTHETA_AGG_FREQ) /
     8             (aTHETA_AGG_GLUC + aTHETA_AGG_LACT +
     9             iTHETA_AGG_MAG + iTHETA_AGG_FREQ)) *
     1             ((ONE + IC_int_AGG_dir) / IC_int_AGG_dir) *
     2             (((iTHETA_AGG_MAG * X_S_AGG_MAG_DEPt) +
     3             (iTHETA_AGG_FREQ * X_S_AGG_FREQ_COMPL)) / (ONE +
     4             (iTHETA_AGG_MAG * X_S_AGG_MAG_DEPt) +
     5             (iTHETA_AGG_FREQ * X_S_AGG_FREQ_COMPL)))))
     6             , 4)

C Aggrecan Expression for IL1B Immunopositive Cells
      AGG_IL1B = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_AGG_LACT * X_S_AGG_LACT) +
     2             (aTHETA_AGG_GLUC * X_S_AGG_GLUC) +
     3             (aTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt)) /
     4             (ONE + (aTHETA_AGG_LACT * X_S_AGG_LACT) +
     5             (aTHETA_AGG_GLUC * X_S_AGG_GLUC) +
     6             (aTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt))) ) * (ONE -
     7             (((iTHETA_AGG_MAG + iTHETA_AGG_FREQ +
     8             iTHETA_AGG_IL1B) /
     9             (aTHETA_AGG_GLUC + aTHETA_AGG_LACT +
     1             iTHETA_AGG_MAG + iTHETA_AGG_FREQ +
     2             iTHETA_AGG_IL1B)) *
     3             ((ONE + IC_int_AGG_dir) / (IC_int_AGG_dir)) *
     4             (((iTHETA_AGG_MAG * X_S_AGG_MAG_DEPt) +
     5             (iTHETA_AGG_FREQ * X_S_AGG_FREQ_COMPL) +
     6             (iTHETA_AGG_IL1B * IL1B_PROT)) / (ONE +
     7             (iTHETA_AGG_MAG * X_S_AGG_MAG_DEPt) +
     8             (iTHETA_AGG_FREQ * X_S_AGG_FREQ_COMPL) +
     9             (iTHETA_AGG_IL1B * IL1B_PROT)))))
     1             , 4 )

C Aggrecan Expression for TNF Immunopositive Cells
      AGG_TNF = ROUND((((ONE + aTHETA_SUM) / aTHETA_SUM) *
     1             (((aTHETA_AGG_LACT * X_S_AGG_LACT) +
     2             (aTHETA_AGG_GLUC * X_S_AGG_GLUC) +
     3             (aTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt))) / (ONE +
     4             (aTHETA_AGG_LACT * X_S_AGG_LACT) +
     5             (aTHETA_AGG_GLUC * X_S_AGG_GLUC) +
     6             (aTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt))) * (ONE -
     7             (((iTHETA_AGG_MAG + iTHETA_AGG_FREQ +
     8             iTHETA_AGG_TNF) /
     9             (aTHETA_AGG_GLUC + aTHETA_AGG_LACT +
     1             iTHETA_AGG_MAG + iTHETA_AGG_FREQ +
     2             iTHETA_AGG_TNF)) *
     3             ((ONE + IC_int_AGG_dir) / IC_int_AGG_dir) *
     4             (((iTHETA_AGG_MAG * X_S_AGG_MAG_DEPt) +
     5             (iTHETA_AGG_FREQ * X_S_AGG_FREQ_COMPL) +
     6             (iTHETA_AGG_TNF * TNF_PROT)) / (ONE +
     7             (iTHETA_AGG_MAG * X_S_AGG_MAG_DEPt) +
     8             (iTHETA_AGG_FREQ * X_S_AGG_FREQ_COMPL) +
     9             (iTHETA_AGG_TNF * TNF_PROT)))))
     1             , 4)

C Aggrecan Expression for IL1B & TNF Immunopositive Cells
      AGG_IL1B_TNF = ROUND((((ONE + aTHETA_SUM) / aTHETA_SUM) *
     1             (((aTHETA_AGG_LACT * X_S_AGG_LACT) +
     2             (aTHETA_AGG_GLUC * X_S_AGG_GLUC) +
     3             (aTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt))) / (ONE +
     4             (aTHETA_AGG_LACT * X_S_AGG_LACT) +
     5             (aTHETA_AGG_GLUC * X_S_AGG_GLUC) +
     6             (aTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt))) * (ONE -
     7             (((iTHETA_AGG_MAG + iTHETA_AGG_FREQ +
     8             iTHETA_AGG_IL1B + iTHETA_AGG_TNF) /
     9             (aTHETA_AGG_GLUC + aTHETA_AGG_LACT +
     1             iTHETA_AGG_MAG + iTHETA_AGG_FREQ +
     2             iTHETA_AGG_IL1B + iTHETA_AGG_TNF)) *
     3             ((ONE + IC_int_AGG_dir) / IC_int_AGG_dir) *
     4             (((iTHETA_AGG_MAG * X_S_AGG_MAG_DEPt) +
     5             (iTHETA_AGG_FREQ * X_S_AGG_FREQ_COMPL) +
     6             (iTHETA_AGG_IL1B * IL1B_PROT) +
     7             (iTHETA_AGG_TNF * TNF_PROT)) / (ONE +
     8             (iTHETA_AGG_MAG * X_S_AGG_MAG_DEPt) +
     9             (iTHETA_AGG_FREQ * X_S_AGG_FREQ_COMPL) +
     1             (iTHETA_AGG_IL1B * IL1B_PROT) +
     2             (iTHETA_AGG_TNF * TNF_PROT)))))
     3             , 4)

C---------------------------------------------------------------------------
C Magnitude >= 1 and Frequency > 3 (Magnitude catabolic, Frequency catabolic)
C---------------------------------------------------------------------------

      ELSE IF (MAG .GE. ONE .AND. FREQ .EQ. ZERO) THEN

C Magnitude Catabolic
         AGG_MAG_CAT_RANGE = A_I_MAG_CAT_2 * AGG_MAG_FACT_2
         IF (AGG_MAG_CAT_RANGE .GT. ONE) THEN
            AGG_MAG_CAT_RANGE = ONE
         END IF
         X_S_AGG_MAG_DEPt = ABS(ZERO - AGG_MAG_CAT_RANGE)

C Frequency Anabolic (Frequency == 0 Case)
         AGG_FREQ_ANA_RANGE = FREQ_0_TIME *
     1                        DBLE(TICKS) * AGG_FREQ_FACT
         IF (AGG_FREQ_ANA_RANGE .GE. FREQ_0) THEN
            AGG_FREQ_ANA_RANGE = FREQ_0
            X_S_AGG_FREQ_COMPL = X_S_AGG_FREQ_COMPL +
     1                           (FREQ_0_TIME * AGG_FREQ_FACT)
            IF (X_S_AGG_FREQ_COMPL .GT. ONE) THEN
               X_S_AGG_FREQ_COMPL = ONE
            END IF
         ELSE
            X_S_AGG_FREQ_COMPL = ZERO
         END IF
         X_S_AGG_FREQ_DEPt = FREQ_0 - AGG_FREQ_ANA_RANGE

C Aggrecan Expression for Immunonegative Cells
         AGG_IMN = ROUND((((ONE + aTHETA_SUM) / aTHETA_SUM) *
     1             (((aTHETA_AGG_LACT * X_S_AGG_LACT) +
     2             (aTHETA_AGG_GLUC * X_S_AGG_GLUC) +
     3             (aTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt))) / (ONE +
     4             (aTHETA_AGG_LACT * X_S_AGG_LACT) +
     5             (aTHETA_AGG_GLUC * X_S_AGG_GLUC) +
     6             (aTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt))) * (ONE -
     7             (((iTHETA_AGG_MAG + iTHETA_AGG_FREQ) /
     8             (aTHETA_AGG_GLUC + aTHETA_AGG_LACT +
     9             iTHETA_AGG_MAG + iTHETA_AGG_FREQ)) *
     1             ((ONE + IC_int_AGG_dir) / IC_int_AGG_dir) *
     2             (((iTHETA_AGG_MAG * X_S_AGG_MAG_DEPt) +
     3             (iTHETA_AGG_FREQ * X_S_AGG_FREQ_COMPL)) / (ONE +
     4             (iTHETA_AGG_MAG * X_S_AGG_MAG_DEPt) +
     5             (iTHETA_AGG_FREQ * X_S_AGG_FREQ_COMPL))))), 4)

C Aggrecan Expression for IL1B Immunopositive Cells
         AGG_IL1B = ROUND((((ONE + aTHETA_SUM) / aTHETA_SUM) *
     1             (((aTHETA_AGG_LACT * X_S_AGG_LACT) +
     2             (aTHETA_AGG_GLUC * X_S_AGG_GLUC) +
     3             (aTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt))) / (ONE +
     4             (aTHETA_AGG_LACT * X_S_AGG_LACT) +
     5             (aTHETA_AGG_GLUC * X_S_AGG_GLUC) +
     6             (aTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt))) * (ONE -
     7             (((iTHETA_AGG_MAG + iTHETA_AGG_FREQ +
     8             iTHETA_AGG_IL1B) / (aTHETA_AGG_GLUC +
     9             aTHETA_AGG_LACT + iTHETA_AGG_MAG +
     1             iTHETA_AGG_FREQ + iTHETA_AGG_IL1B)) *
     2             ((ONE + IC_int_AGG_dir) / IC_int_AGG_dir) *
     3             (((iTHETA_AGG_MAG * X_S_AGG_MAG_DEPt) +
     4             (iTHETA_AGG_FREQ * X_S_AGG_FREQ_COMPL) +
     5             (iTHETA_AGG_IL1B * IL1B_PROT)) / (ONE +
     6             (iTHETA_AGG_MAG * X_S_AGG_MAG_DEPt) +
     7             (iTHETA_AGG_FREQ * X_S_AGG_FREQ_COMPL) +
     8             (iTHETA_AGG_IL1B * IL1B_PROT))))), 4)

C Aggrecan Expression for TNF Immunopositive Cells
         AGG_TNF = ROUND((((ONE + aTHETA_SUM) / aTHETA_SUM) *
     1             (((aTHETA_AGG_LACT * X_S_AGG_LACT) +
     2             (aTHETA_AGG_GLUC * X_S_AGG_GLUC) +
     3             (aTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt))) / (ONE +
     4             (aTHETA_AGG_LACT * X_S_AGG_LACT) +
     5             (aTHETA_AGG_GLUC * X_S_AGG_GLUC) +
     6             (aTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt))) * (ONE -
     7             (((iTHETA_AGG_MAG + iTHETA_AGG_FREQ +
     8             iTHETA_AGG_TNF) / (aTHETA_AGG_GLUC +
     9             aTHETA_AGG_LACT + iTHETA_AGG_MAG +
     1             iTHETA_AGG_FREQ + iTHETA_AGG_TNF)) *
     2             ((ONE + IC_int_AGG_dir) / IC_int_AGG_dir) *
     3             (((iTHETA_AGG_MAG * X_S_AGG_MAG_DEPt) +
     4             (iTHETA_AGG_FREQ * X_S_AGG_FREQ_COMPL) +
     5             (iTHETA_AGG_TNF * TNF_PROT)) / (ONE +
     6             (iTHETA_AGG_MAG * X_S_AGG_MAG_DEPt) +
     7             (iTHETA_AGG_FREQ * X_S_AGG_FREQ_COMPL) +
     8             (iTHETA_AGG_TNF * TNF_PROT))))), 4)

C Aggrecan Expression for IL1B & TNF Immunopositive Cells
         AGG_IL1B_TNF = ROUND((((ONE + aTHETA_SUM) / aTHETA_SUM) *
     1             (((aTHETA_AGG_LACT * X_S_AGG_LACT) +
     2             (aTHETA_AGG_GLUC * X_S_AGG_GLUC) +
     3             (aTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt))) / (ONE +
     4             (aTHETA_AGG_LACT * X_S_AGG_LACT) +
     5             (aTHETA_AGG_GLUC * X_S_AGG_GLUC) +
     6             (aTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt))) * (ONE -
     7             (((iTHETA_AGG_MAG + iTHETA_AGG_FREQ +
     8             iTHETA_AGG_IL1B + iTHETA_AGG_TNF) /
     9             (aTHETA_AGG_GLUC + aTHETA_AGG_LACT +
     1             iTHETA_AGG_MAG + iTHETA_AGG_FREQ +
     2             iTHETA_AGG_IL1B + iTHETA_AGG_TNF)) *
     3             ((ONE + IC_int_AGG_dir) / IC_int_AGG_dir) *
     4             (((iTHETA_AGG_MAG * X_S_AGG_MAG_DEPt) +
     5             (iTHETA_AGG_FREQ * X_S_AGG_FREQ_COMPL) +
     6             (iTHETA_AGG_IL1B * IL1B_PROT) +
     7             (iTHETA_AGG_TNF * TNF_PROT)) / (ONE +
     8             (iTHETA_AGG_MAG * X_S_AGG_MAG_DEPt) +
     9             (iTHETA_AGG_FREQ * X_S_AGG_FREQ_COMPL) +
     1             (iTHETA_AGG_IL1B * IL1B_PROT) +
     2             (iTHETA_AGG_TNF * TNF_PROT))))), 4)

C---------------------------------------------------------------------------
C Magnitude >= 1 and Frequency >= 3 (Magnitude and Frequency catabolic)
C---------------------------------------------------------------------------

      ELSE IF (MAG .GE. ONE .AND. FREQ .GT. 3.0D0) THEN

C Magnitude Catabolic
         AGG_MAG_CAT_RANGE = A_I_MAG_CAT_2 * AGG_MAG_FACT_2
         IF (AGG_MAG_CAT_RANGE .GT. ONE) THEN
            AGG_MAG_CAT_RANGE = ONE
         END IF
         X_S_AGG_MAG_DEPt = ABS(ZERO - AGG_MAG_CAT_RANGE)

C Frequency Catabolic
         AGG_FREQ_CAT_RANGE = A_I_FREQ_CAT_2 * AGG_FREQ_FACT_2
         IF (AGG_FREQ_CAT_RANGE .GT. ONE) THEN
            AGG_FREQ_CAT_RANGE = ONE
         END IF
         X_S_AGG_FREQ_DEPt = ABS(ZERO - AGG_FREQ_CAT_RANGE)

C Aggrecan Expression for Immunonegative Cells
         AGG_IMN = ROUND((((ONE + aTHETA_SUM) / aTHETA_SUM) *
     1             (((aTHETA_AGG_LACT * X_S_AGG_LACT) +
     2             (aTHETA_AGG_GLUC * X_S_AGG_GLUC))) / (ONE +
     3             (aTHETA_AGG_LACT * X_S_AGG_LACT) +
     4             (aTHETA_AGG_GLUC * X_S_AGG_GLUC))) * (ONE -
     5             (((iTHETA_AGG_MAG + iTHETA_AGG_FREQ) /
     6             (aTHETA_AGG_GLUC + aTHETA_AGG_LACT +
     7             iTHETA_AGG_MAG + iTHETA_AGG_FREQ)) *
     8             ((ONE + IC_int_AGG_dir) / IC_int_AGG_dir) *
     9             (((iTHETA_AGG_MAG * X_S_AGG_MAG_DEPt) +
     1             (iTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt)) / (ONE +
     2             (iTHETA_AGG_MAG * X_S_AGG_MAG_DEPt) +
     3             (iTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt))))), 4)

C Aggrecan Expression for IL1B Immunopositive Cells
         AGG_IL1B = ROUND((((ONE + aTHETA_SUM) / aTHETA_SUM) *
     1             (((aTHETA_AGG_LACT * X_S_AGG_LACT) +
     2             (aTHETA_AGG_GLUC * X_S_AGG_GLUC))) / (ONE +
     3             (aTHETA_AGG_LACT * X_S_AGG_LACT) +
     4             (aTHETA_AGG_GLUC * X_S_AGG_GLUC))) * (ONE -
     5             (((iTHETA_AGG_MAG + iTHETA_AGG_FREQ +
     6             iTHETA_AGG_IL1B) / (aTHETA_AGG_GLUC +
     7             aTHETA_AGG_LACT + iTHETA_AGG_MAG +
     8             iTHETA_AGG_FREQ + iTHETA_AGG_IL1B)) *
     9             ((ONE + IC_int_AGG_dir) / IC_int_AGG_dir) *
     1             (((iTHETA_AGG_MAG * X_S_AGG_MAG_DEPt) +
     2             (iTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt) +
     3             (iTHETA_AGG_IL1B * IL1B_PROT)) / (ONE +
     4             (iTHETA_AGG_MAG * X_S_AGG_MAG_DEPt) +
     5             (iTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt) +
     6             (iTHETA_AGG_IL1B * IL1B_PROT))))), 4)

C Aggrecan Expression for TNF Immunopositive Cells
         AGG_TNF = ROUND((((ONE + aTHETA_SUM) / aTHETA_SUM) *
     1             (((aTHETA_AGG_LACT * X_S_AGG_LACT) +
     2             (aTHETA_AGG_GLUC * X_S_AGG_GLUC))) / (ONE +
     3             (aTHETA_AGG_LACT * X_S_AGG_LACT) +
     4             (aTHETA_AGG_GLUC * X_S_AGG_GLUC))) * (ONE -
     5             (((iTHETA_AGG_MAG + iTHETA_AGG_FREQ +
     6             iTHETA_AGG_TNF) / (aTHETA_AGG_GLUC +
     7             aTHETA_AGG_LACT + iTHETA_AGG_MAG +
     8             iTHETA_AGG_FREQ + iTHETA_AGG_TNF)) *
     9             ((ONE + IC_int_AGG_dir) / IC_int_AGG_dir) *
     1             (((iTHETA_AGG_MAG * X_S_AGG_MAG_DEPt) +
     2             (iTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt) +
     3             (iTHETA_AGG_TNF * TNF_PROT)) / (ONE +
     4             (iTHETA_AGG_MAG * X_S_AGG_MAG_DEPt) +
     5             (iTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt) +
     6             (iTHETA_AGG_TNF * TNF_PROT))))), 4)

C Aggrecan Expression for IL1B & TNF Immunopositive Cells
         AGG_IL1B_TNF = ROUND((((ONE + aTHETA_SUM) / aTHETA_SUM) *
     1             (((aTHETA_AGG_LACT * X_S_AGG_LACT) +
     2             (aTHETA_AGG_GLUC * X_S_AGG_GLUC))) / (ONE +
     3             (aTHETA_AGG_LACT * X_S_AGG_LACT) +
     4             (aTHETA_AGG_GLUC * X_S_AGG_GLUC))) * (ONE -
     5             (((iTHETA_AGG_MAG + iTHETA_AGG_FREQ +
     6             iTHETA_AGG_IL1B + iTHETA_AGG_TNF) /
     7             (aTHETA_AGG_GLUC + aTHETA_AGG_LACT +
     8             iTHETA_AGG_MAG + iTHETA_AGG_FREQ +
     9             iTHETA_AGG_IL1B + iTHETA_AGG_TNF)) *
     1             ((ONE + IC_int_AGG_dir) / IC_int_AGG_dir) *
     2             (((iTHETA_AGG_MAG * X_S_AGG_MAG_DEPt) +
     3             (iTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt) +
     4             (iTHETA_AGG_IL1B * IL1B_PROT) +
     5             (iTHETA_AGG_TNF * TNF_PROT)) / (ONE +
     6             (iTHETA_AGG_MAG * X_S_AGG_MAG_DEPt) +
     7             (iTHETA_AGG_FREQ * X_S_AGG_FREQ_DEPt) +
     8             (iTHETA_AGG_IL1B * IL1B_PROT) +
     9             (iTHETA_AGG_TNF * TNF_PROT))))), 4)

      END IF

C--------------------------------------------------------------------------
      RETURN
      END SUBROUTINE PNeq_AGG

C---------------------------------------------------------------------------------------------------------------------------------------------
C Evaluate the collagen I expression
C --------------------------------------------------------------

      SUBROUTINE PNeq_COLI(COLI_IMN, COLI_IL1B, COLI_TNF, COLI_IL1B_TNF,
     1           TICKS,
     2           IL1B_PROT, TNF_PROT,
     3           GLUC_int, pH_val, MAG, FREQ,
     4           FREQ_0, FREQ_0_time,
     5           MAG_0, MAG_0_time,
     6           aTHETA_COLI_GLUC, aTHETA_COLI_LACT,
     7           aTHETA_COLI_IL1B, iTHETA_COLI_TNF,
     8           aTHETA_COLI_MAG, iTHETA_COLI_MAG,
     9           aTHETA_COLI_FREQ, iTHETA_COLI_FREQ,
     1           AC_int_IL1B_dir, AC_int_TNF_dir,
     2           IC_int_IL1B_dir, IC_int_TNF_dir,
     3           IC_int_COLI_dir,
     4           X_S_COLI_GLUC, X_S_COLI_LACT,
     5           X_S_COLI_MAG, X_S_COLI_FREQ,
     6           A_I_MAG_ANA, A_I_MAG_CAT_1, A_I_MAG_CAT_2,
     7           I_A_MAG_ANA, I_A_MAG_CAT_1, I_A_MAG_CAT_2,
     8           A_I_FREQ_ANA, A_I_FREQ_CAT_1, A_I_FREQ_CAT_2,
     9           I_A_FREQ_ANA, I_A_FREQ_CAT_1, I_A_FREQ_CAT_2,
     1           COLI_MAG_FACT, COLI_MAG_FACT_2,
     2           COLI_FREQ_FACT, COLI_FREQ_FACT_2,
     3           X_S_COLI_MAG_COMPL, X_S_COLI_FREQ_COMPL,
     4           aTHETA_SUM)

C VARIABLES:
C X_S: SENSITIVITY OF A CELL ACTIVITY (CA) TO A STIMULUS CONCENTRATION
C THETA: SENSITIVITY OF A CELL ACTIVITY (CA) TO A STIMULUS TYPE
C OUTPUTS:
C COLI_IMN: COLIRECAN EXPRESSION FOR IMMUNONEGATIVE CELLS
C COLI_IL1B: CONTRIBUTION FROM IL1B
C COLI_TNF: CONTRIBUTION FROM TNF
C COLI_IL1B_TNF: COMBINED CONTRIBUTION OF IL1B AND TNF

C ********************************************************************
C
      IMPLICIT NONE
C
C FUNCTION
      DOUBLE PRECISION ROUND

C INTEGER VARIABLES
      INTEGER TICKS

C SPECIAL CASE FOR FREQ = 0
      DOUBLE PRECISION FREQ_0, FREQ_0_time
      DOUBLE PRECISION MAG_0, MAG_0_time

C Input and Output Variables
      DOUBLE PRECISION COLI_IMN, COLI_IL1B, COLI_TNF, COLI_IL1B_TNF
      DOUBLE PRECISION IL1B_PROT, TNF_PROT
      DOUBLE PRECISION GLUC_int, pH_val, MAG, FREQ

C Constants for Sensitivity Factors
      DOUBLE PRECISION aTHETA_COLI_GLUC, aTHETA_COLI_LACT
      DOUBLE PRECISION aTHETA_COLI_IL1B, iTHETA_COLI_TNF
      DOUBLE PRECISION aTHETA_COLI_MAG, iTHETA_COLI_MAG
      DOUBLE PRECISION aTHETA_COLI_FREQ, iTHETA_COLI_FREQ
      DOUBLE PRECISION AC_int_IL1B_dir, AC_int_TNF_dir
      DOUBLE PRECISION IC_int_IL1B_dir, IC_int_TNF_dir
      DOUBLE PRECISION IC_int_COLI_dir
      DOUBLE PRECISION aTHETA_SUM

C Intermediate Variables
      DOUBLE PRECISION X_S_COLI_GLUC, X_S_COLI_LACT
      DOUBLE PRECISION X_S_COLI_MAG, X_S_COLI_FREQ
      DOUBLE PRECISION COLI_MAG_FACT, COLI_MAG_FACT_2
      DOUBLE PRECISION COLI_FREQ_FACT, COLI_FREQ_FACT_2

C Completive factors
      DOUBLE PRECISION X_S_COLI_MAG_COMPL
      DOUBLE PRECISION X_S_COLI_MAG_DEPt
      DOUBLE PRECISION X_S_COLI_FREQ_COMPL
      DOUBLE PRECISION X_S_COLI_FREQ_DEPt

C Time dependent funtions
      DOUBLE PRECISION A_I_MAG_ANA, A_I_MAG_CAT_1, A_I_MAG_CAT_2
      DOUBLE PRECISION I_A_MAG_ANA, I_A_MAG_CAT_1, I_A_MAG_CAT_2
      DOUBLE PRECISION A_I_FREQ_ANA, A_I_FREQ_CAT_1, A_I_FREQ_CAT_2
      DOUBLE PRECISION I_A_FREQ_ANA, I_A_FREQ_CAT_1, I_A_FREQ_CAT_2

C Time sensitivity range
      DOUBLE PRECISION COLI_MAG_ANA_RANGE
      DOUBLE PRECISION COLI_FREQ_ANA_RANGE
      DOUBLE PRECISION COLI_MAG_CAT_RANGE
      DOUBLE PRECISION COLI_FREQ_CAT_RANGE

C CONSTANT PARAMETERS
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.D0, ONE = 1.D0

C---------------------------------------------------------------------------
C Initialize Variables
      COLI_MAG_ANA_RANGE = ZERO
      COLI_FREQ_ANA_RANGE = ZERO
      COLI_MAG_CAT_RANGE = ZERO
      COLI_FREQ_CAT_RANGE = ZERO
      X_S_COLI_MAG_DEPt = ZERO
      X_S_COLI_FREQ_DEPt = ZERO

C---------------------------------------------------------------------------
C Magnitude <= 0.998 and Frequency <= 3 (Magnitude anabolic, Frequency anabolic)
C---------------------------------------------------------------------------

      IF (MAG .LE. 0.998D0 .AND. FREQ .LE. 3.0D0 .AND.
     1    FREQ .NE. ZERO) THEN

C Magnitude Anabolic
         COLI_MAG_ANA_RANGE = I_A_MAG_ANA * COLI_MAG_FACT

         IF (COLI_MAG_ANA_RANGE .GE. ABS(X_S_COLI_MAG)) THEN
               COLI_MAG_ANA_RANGE = ABS(X_S_COLI_MAG)
               X_S_COLI_MAG_COMPL = X_S_COLI_MAG_COMPL +
     1                               (I_A_MAG_CAT_1 * COLI_MAG_FACT)

               IF (X_S_COLI_MAG_COMPL .GT. ONE) THEN
                     X_S_COLI_MAG_COMPL = ONE  ! Maximum value is 1
               END IF
         ELSE
               X_S_COLI_MAG_COMPL = ZERO
         END IF

         X_S_COLI_MAG_DEPt = ABS(X_S_COLI_MAG + COLI_MAG_ANA_RANGE)

C Frequency Anabolic
         COLI_FREQ_ANA_RANGE = I_A_FREQ_ANA * COLI_FREQ_FACT

         IF (COLI_FREQ_ANA_RANGE .GE. ABS(X_S_COLI_FREQ)) THEN
               COLI_FREQ_ANA_RANGE = ABS(X_S_COLI_FREQ)
               X_S_COLI_FREQ_COMPL = X_S_COLI_FREQ_COMPL +
     1                               (I_A_FREQ_CAT_1 * COLI_FREQ_FACT)

               IF (X_S_COLI_FREQ_COMPL .GT. ONE) THEN
                     X_S_COLI_FREQ_COMPL = ONE
               END IF
         ELSE
               X_S_COLI_FREQ_COMPL = ZERO
         END IF

         X_S_COLI_FREQ_DEPt = ABS(X_S_COLI_FREQ + COLI_FREQ_ANA_RANGE)

C ColI mRNA expression for immunonegative cells
            COLI_IMN = ROUND((((ONE + aTHETA_SUM) / aTHETA_SUM) *
     1             (((aTHETA_COLI_LACT * X_S_COLI_LACT) +
     2             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     3             (aTHETA_COLI_MAG * X_S_COLI_MAG_COMPL) +
     4             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_COMPL)) / 
     5             (ONE + (aTHETA_COLI_LACT * X_S_COLI_LACT) +
     6             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     7             (aTHETA_COLI_MAG * X_S_COLI_MAG_COMPL) +
     8             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_COMPL)))) *
     9             (ONE - (((iTHETA_COLI_MAG + iTHETA_COLI_FREQ) /
     1             (aTHETA_COLI_GLUC + aTHETA_COLI_LACT +
     2             iTHETA_COLI_MAG + iTHETA_COLI_FREQ)) *
     3             ((ONE + IC_int_COLI_dir) / IC_int_COLI_dir) *
     4             (((iTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     5             (iTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt)) / 
     6             (ONE + (iTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     7             (iTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt)))))
     8             , 4)

C ColI mRNA expression for IL1B immunopositive cells
            COLI_IL1B = ROUND((((ONE + aTHETA_SUM) / aTHETA_SUM) *
     1             (((aTHETA_COLI_LACT * X_S_COLI_LACT) +
     2             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     3             (aTHETA_COLI_MAG * X_S_COLI_MAG_COMPL) +
     4             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_COMPL) +
     5             (aTHETA_COLI_IL1B * IL1B_PROT)) / 
     6             (ONE + (aTHETA_COLI_LACT * X_S_COLI_LACT) +
     7             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     8             (aTHETA_COLI_MAG * X_S_COLI_MAG_COMPL) +
     9             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_COMPL) +
     1             (aTHETA_COLI_IL1B * IL1B_PROT)))) *
     2             (ONE - (((iTHETA_COLI_MAG + iTHETA_COLI_FREQ) /
     3             (aTHETA_COLI_GLUC + aTHETA_COLI_LACT +
     4             aTHETA_COLI_IL1B +
     5             iTHETA_COLI_MAG + iTHETA_COLI_FREQ)) *
     6             ((ONE + IC_int_COLI_dir) / IC_int_COLI_dir) *
     7             (((iTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     8             (iTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt)) / 
     9             (ONE + (iTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     1             (iTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt)))))
     2             , 4)

C ColI mRNA expression for TNF immunopositive cells
            COLI_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_COLI_LACT * X_S_COLI_LACT) +
     2             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     3             (aTHETA_COLI_MAG * X_S_COLI_MAG_COMPL) +
     4             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_COMPL)) / (ONE +
     5             (aTHETA_COLI_LACT * X_S_COLI_LACT) +
     6             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     7             (aTHETA_COLI_MAG * X_S_COLI_MAG_COMPL) +
     8             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_COMPL))) ) * (ONE -
     9             (((iTHETA_COLI_MAG + iTHETA_COLI_FREQ +
     1             iTHETA_COLI_TNF) /
     2             (aTHETA_COLI_GLUC + aTHETA_COLI_LACT +
     3             iTHETA_COLI_MAG + iTHETA_COLI_FREQ +
     4             iTHETA_COLI_TNF)) *
     5             ((ONE + IC_int_COLI_dir) / (IC_int_COLI_dir)) *
     6             (((iTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     7             (iTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt) +
     8             (iTHETA_COLI_TNF * TNF_PROT)) / (ONE +
     9             (iTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     1             (iTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt) +
     2             (iTHETA_COLI_TNF * TNF_PROT))) ))
     3             , 4)

C ColI mRNA expression for IL1B and TNF immunopositive cells
            COLI_IL1B_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_COLI_LACT * X_S_COLI_LACT) +
     2             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     3             (aTHETA_COLI_MAG * X_S_COLI_MAG_COMPL) +
     4             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_COMPL) +
     5             (aTHETA_COLI_IL1B * IL1B_PROT)) / (ONE +
     6             (aTHETA_COLI_LACT * X_S_COLI_LACT) +
     7             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     8             (aTHETA_COLI_MAG * X_S_COLI_MAG_COMPL) +
     9             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_COMPL) +
     1             (aTHETA_COLI_IL1B * IL1B_PROT))) ) * (ONE -
     2             (((iTHETA_COLI_MAG + iTHETA_COLI_FREQ +
     3             iTHETA_COLI_TNF) /
     4             (aTHETA_COLI_GLUC + aTHETA_COLI_LACT +
     5             aTHETA_COLI_IL1B +
     6             iTHETA_COLI_MAG + iTHETA_COLI_FREQ +
     7             iTHETA_COLI_TNF)) *
     8             ((ONE + IC_int_COLI_dir) / (IC_int_COLI_dir)) *
     9             (((iTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     1             (iTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt) +
     2             (iTHETA_COLI_TNF * TNF_PROT)) / (ONE +
     3             (iTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     4             (iTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt) +
     5             (iTHETA_COLI_TNF * TNF_PROT)))))
     6             , 4 )

C---------------------------------------------------------------------------
C Magnitude <= 0.998 and Frequency = 0 (Magnitude anabolic, Frequency catabolic)
C---------------------------------------------------------------------------

      ELSE IF (MAG .LE. 0.998D0 .AND. FREQ .EQ. ZERO) THEN

C Magnitude Anabolic
      COLI_MAG_ANA_RANGE = I_A_MAG_ANA * COLI_MAG_FACT

      IF (COLI_MAG_ANA_RANGE .GE. ABS(X_S_COLI_MAG)) THEN
            COLI_MAG_ANA_RANGE = ABS(X_S_COLI_MAG)
            X_S_COLI_MAG_COMPL = X_S_COLI_MAG_COMPL +
     1                               (I_A_MAG_CAT_1 * COLI_MAG_FACT)

            IF (X_S_COLI_MAG_COMPL .GE. ONE) THEN
                  X_S_COLI_MAG_COMPL = ONE
            END IF
      ELSE
            X_S_COLI_MAG_COMPL = ZERO
      END IF

      X_S_COLI_MAG_DEPt = ABS(X_S_COLI_MAG + COLI_MAG_ANA_RANGE)

C Frequency Anabolic (Frequency == 0 Case)
      COLI_FREQ_ANA_RANGE = FREQ_0_TIME *
     1                        DBLE(TICKS) * COLI_FREQ_FACT

      IF (COLI_FREQ_ANA_RANGE .GE. FREQ_0) THEN
            COLI_FREQ_ANA_RANGE = FREQ_0
            X_S_COLI_FREQ_COMPL = X_S_COLI_FREQ_COMPL +
     1                               (FREQ_0_TIME * COLI_FREQ_FACT)

            IF (X_S_COLI_FREQ_COMPL .GT. ONE) THEN
                  X_S_COLI_FREQ_COMPL = ONE
            END IF
      ELSE
            X_S_COLI_FREQ_COMPL = ZERO
      END IF

      X_S_COLI_FREQ_DEPt = FREQ_0 - COLI_FREQ_ANA_RANGE

C ColI mRNA expression for immunonegative cells
            COLI_IMN = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_COLI_LACT * X_S_COLI_LACT) +
     2             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     3             (aTHETA_COLI_MAG * X_S_COLI_MAG_COMPL) +
     4             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_COMPL)) / (ONE +
     5             (aTHETA_COLI_LACT * X_S_COLI_LACT) +
     6             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     7             (aTHETA_COLI_MAG * X_S_COLI_MAG_COMPL) +
     8             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_COMPL)))) * (ONE -
     9             (((iTHETA_COLI_MAG + iTHETA_COLI_FREQ) /
     1             (aTHETA_COLI_GLUC + aTHETA_COLI_LACT +
     2             iTHETA_COLI_MAG + iTHETA_COLI_FREQ)) *
     3             ((ONE + IC_int_COLI_dir) / (IC_int_COLI_dir)) *
     4             (((iTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     5             (iTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt)) / (ONE +
     6             (iTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     7             (iTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt)))))
     8             , 4)

C ColI mRNA expression for IL1B immunopositive cells
            COLI_IL1B = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_COLI_LACT * X_S_COLI_LACT) +
     2             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     3             (aTHETA_COLI_MAG * X_S_COLI_MAG_COMPL) +
     4             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_COMPL) +
     5             (aTHETA_COLI_IL1B * IL1B_PROT)) / (ONE +
     6             (aTHETA_COLI_LACT * X_S_COLI_LACT) +
     7             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     8             (aTHETA_COLI_MAG * X_S_COLI_MAG_COMPL) +
     9             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_COMPL) +
     1             (aTHETA_COLI_IL1B * IL1B_PROT)))) * (ONE -
     2             (((iTHETA_COLI_MAG + iTHETA_COLI_FREQ) /
     3             (aTHETA_COLI_GLUC + aTHETA_COLI_LACT +
     4             aTHETA_COLI_IL1B +
     5             iTHETA_COLI_MAG + iTHETA_COLI_FREQ)) *
     6             ((ONE + IC_int_COLI_dir) / (IC_int_COLI_dir)) *
     7             (((iTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     9             (iTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt)) / (ONE +
     1             (iTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     2             (iTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt)))))
     3             , 4)

C ColI mRNA expression for TNF immunopositive cells
            COLI_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_COLI_LACT * X_S_COLI_LACT) +
     2             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     3             (aTHETA_COLI_MAG * X_S_COLI_MAG_COMPL) +
     4             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_COMPL)) / (ONE +
     5             (aTHETA_COLI_LACT * X_S_COLI_LACT) +
     6             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     7             (aTHETA_COLI_MAG * X_S_COLI_MAG_COMPL) +
     8             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_COMPL)))) * (ONE -
     9             (((iTHETA_COLI_MAG + iTHETA_COLI_FREQ +
     1             iTHETA_COLI_TNF) /
     2             (aTHETA_COLI_GLUC + aTHETA_COLI_LACT +
     3             iTHETA_COLI_MAG + iTHETA_COLI_FREQ +
     4             iTHETA_COLI_TNF)) *
     5             ((ONE + IC_int_COLI_dir) / (IC_int_COLI_dir)) *
     6             (((iTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     7             (iTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt) +
     8             (iTHETA_COLI_TNF * TNF_PROT)) / (ONE +
     9             (iTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     1             (iTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt) +
     2             (iTHETA_COLI_TNF * TNF_PROT)))))
     3             , 4)

C ColI mRNA expression for IL1B and TNF immunopositive cells
            COLI_IL1B_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_COLI_LACT * X_S_COLI_LACT) +
     2             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     3             (aTHETA_COLI_MAG * X_S_COLI_MAG_COMPL) +
     4             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_COMPL) +
     5             (aTHETA_COLI_IL1B * IL1B_PROT)) / (ONE +
     6             (aTHETA_COLI_LACT * X_S_COLI_LACT) +
     7             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     8             (aTHETA_COLI_MAG * X_S_COLI_MAG_COMPL) +
     9             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_COMPL) +
     1             (aTHETA_COLI_IL1B * IL1B_PROT))) ) * (ONE -
     2             (((iTHETA_COLI_MAG + iTHETA_COLI_FREQ +
     3             iTHETA_COLI_TNF) /
     4             (aTHETA_COLI_GLUC + aTHETA_COLI_LACT +
     5             aTHETA_COLI_IL1B +
     6             iTHETA_COLI_MAG + iTHETA_COLI_FREQ +
     7             iTHETA_COLI_TNF)) *
     8             ((ONE + IC_int_COLI_dir) / (IC_int_COLI_dir)) *
     9             (((iTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     1             (iTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt) +
     2             (iTHETA_COLI_TNF * TNF_PROT)) / (ONE +
     3             (iTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     4             (iTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt) +
     5             (iTHETA_COLI_TNF * TNF_PROT))) ))
     6             , 4)

C---------------------------------------------------------------------------
C Magnitude <= 0.998 and Frequency > 3 (Magnitude anabolic, Frequency catabolic)
C---------------------------------------------------------------------------

      ELSE IF (MAG .LE. 0.998D0 .AND. FREQ .GT. 3.0D0) THEN

C Magnitude Anabolic
            COLI_MAG_ANA_RANGE = I_A_MAG_ANA * COLI_MAG_FACT

            IF (COLI_MAG_ANA_RANGE .GE. ABS(X_S_COLI_MAG)) THEN
                  COLI_MAG_ANA_RANGE = ABS(X_S_COLI_MAG)
                  X_S_COLI_MAG_COMPL = X_S_COLI_MAG_COMPL +
     1                               (I_A_MAG_CAT_1 * COLI_MAG_FACT)

                  IF (X_S_COLI_MAG_COMPL .GT. ONE) THEN
                        X_S_COLI_MAG_COMPL = ONE
                  END IF
            ELSE
                  X_S_COLI_MAG_COMPL = ZERO
            END IF

            X_S_COLI_MAG_DEPt = ABS(X_S_COLI_MAG + COLI_MAG_ANA_RANGE)

C Frequency Catabolic
            COLI_FREQ_CAT_RANGE = I_A_FREQ_CAT_2 * COLI_FREQ_FACT_2

            IF (COLI_FREQ_CAT_RANGE .GT. ONE) THEN
                  COLI_FREQ_CAT_RANGE = ONE
            END IF

            X_S_COLI_FREQ_DEPt = ZERO + COLI_FREQ_CAT_RANGE

C ColI mRNA expression for immunonegative cells
            COLI_IMN = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_COLI_LACT * X_S_COLI_LACT) +
     2             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     3             (aTHETA_COLI_MAG * X_S_COLI_MAG_COMPL) +
     4             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt)) / (ONE +
     5             (aTHETA_COLI_LACT * X_S_COLI_LACT) +
     6             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     7             (aTHETA_COLI_MAG * X_S_COLI_MAG_COMPL) +
     8             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt))) ) * (ONE -
     9             ((iTHETA_COLI_MAG /
     1             (aTHETA_COLI_GLUC + aTHETA_COLI_LACT +
     2             iTHETA_COLI_MAG + aTHETA_COLI_FREQ)) *
     3             ((ONE + IC_int_COLI_dir) / (IC_int_COLI_dir)) *
     4             ((iTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) / (ONE +
     5             (iTHETA_COLI_MAG * X_S_COLI_MAG_DEPt)))))
     6             , 4)

C ColI mRNA expression for IL1B immunopositive cells
            COLI_IL1B = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_COLI_LACT * X_S_COLI_LACT) +
     2             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     3             (aTHETA_COLI_MAG * X_S_COLI_MAG_COMPL) +
     4             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt) +
     5             (aTHETA_COLI_IL1B * IL1B_PROT)) / (ONE +
     6             (aTHETA_COLI_LACT * X_S_COLI_LACT) +
     7             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     8             (aTHETA_COLI_MAG * X_S_COLI_MAG_COMPL) +
     9             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt) +
     1             (aTHETA_COLI_IL1B * IL1B_PROT))) ) * (ONE -
     2             ((iTHETA_COLI_MAG /
     3             (aTHETA_COLI_GLUC + aTHETA_COLI_LACT +
     4             aTHETA_COLI_IL1B +
     5             iTHETA_COLI_MAG + aTHETA_COLI_FREQ)) *
     6             ((ONE + IC_int_COLI_dir) / (IC_int_COLI_dir)) *
     7             ((iTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) / (ONE +
     8             (iTHETA_COLI_MAG * X_S_COLI_MAG_DEPt)))))
     9             , 4)

C ColI mRNA expression for TNF immunopositive cells
            COLI_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_COLI_LACT * X_S_COLI_LACT) +
     2             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     3             (aTHETA_COLI_MAG * X_S_COLI_MAG_COMPL) +
     4             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt)) / (ONE +
     5             (aTHETA_COLI_LACT * X_S_COLI_LACT) +
     6             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     7             (aTHETA_COLI_MAG * X_S_COLI_MAG_COMPL) +
     8             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt)))) * (ONE -
     9             ((iTHETA_COLI_MAG +
     1             iTHETA_COLI_TNF /
     2             (aTHETA_COLI_GLUC + aTHETA_COLI_LACT +
     3             iTHETA_COLI_MAG + aTHETA_COLI_FREQ +
     4             iTHETA_COLI_TNF)) *
     5             ((ONE + IC_int_COLI_dir) / (IC_int_COLI_dir)) *
     6             (((iTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     7             (iTHETA_COLI_TNF * TNF_PROT)) / (ONE +
     8             (iTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     9             (iTHETA_COLI_TNF * TNF_PROT)))))
     1             , 4)

C ColI mRNA expression for IL1B and TNF immunopositive cells
            COLI_IL1B_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_COLI_LACT * X_S_COLI_LACT) +
     2             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     3             (aTHETA_COLI_MAG * X_S_COLI_MAG_COMPL) +
     4             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt) +
     5             (aTHETA_COLI_IL1B * IL1B_PROT)) / (ONE +
     6             (aTHETA_COLI_LACT * X_S_COLI_LACT) +
     7             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     8             (aTHETA_COLI_MAG * X_S_COLI_MAG_COMPL) +
     9             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt) +
     1             (aTHETA_COLI_IL1B * IL1B_PROT))) ) * (ONE -
     2             ((iTHETA_COLI_MAG +
     3             iTHETA_COLI_TNF /
     4             (aTHETA_COLI_GLUC + aTHETA_COLI_LACT +
     5             aTHETA_COLI_FREQ +
     6             aTHETA_COLI_IL1B +
     7             iTHETA_COLI_MAG +
     8             iTHETA_COLI_TNF)) *
     9             ((IC_int_COLI_dir) / (IC_int_COLI_dir)) *
     1             (((iTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     2             (iTHETA_COLI_TNF * TNF_PROT)) / (ONE +
     3             (iTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     4             (iTHETA_COLI_TNF * TNF_PROT)))))
     5             , 4)

C---------------------------------------------------------------------------
C Magnitude >= 1 and Frequency <= 3 (Magnitude catabolic, Frequency anabolic)
C---------------------------------------------------------------------------

      ELSE IF (MAG .GE. ONE .AND. FREQ .LE. 3.0D0 .AND.
     1         FREQ .NE. ZERO) THEN

C Magnitude Catabolic
            COLI_MAG_CAT_RANGE = I_A_MAG_CAT_2 * COLI_MAG_FACT_2

            IF (COLI_MAG_CAT_RANGE .GT. ONE) THEN
                  COLI_MAG_CAT_RANGE = ONE
            END IF

            X_S_COLI_MAG_DEPt = ZERO + COLI_MAG_CAT_RANGE

C Frequency Anabolic
            COLI_FREQ_ANA_RANGE = I_A_FREQ_ANA * COLI_FREQ_FACT

            IF (COLI_FREQ_ANA_RANGE .GE. ABS(X_S_COLI_FREQ)) THEN
                  COLI_FREQ_ANA_RANGE = ABS(X_S_COLI_FREQ)
                  X_S_COLI_FREQ_COMPL = X_S_COLI_FREQ_COMPL +
     1                                 (I_A_FREQ_CAT_1 * COLI_FREQ_FACT)

                  IF (X_S_COLI_FREQ_COMPL .GT. ONE) THEN
                        X_S_COLI_FREQ_COMPL = ONE
                  END IF
            ELSE
                  X_S_COLI_FREQ_COMPL = ZERO
            END IF

            X_S_COLI_FREQ_DEPt = ABS(X_S_COLI_FREQ +
     1                           COLI_FREQ_ANA_RANGE)

C ColI mRNA expression for immunonegative cells
            COLI_IMN = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_COLI_LACT * X_S_COLI_LACT) +
     2             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     3             (aTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     4             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_COMPL)) / (ONE +
     5             (aTHETA_COLI_LACT * X_S_COLI_LACT) +
     6             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     7             (aTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     8             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_COMPL)))) * (ONE -
     9             (((iTHETA_COLI_FREQ) /
     1             (aTHETA_COLI_GLUC + aTHETA_COLI_LACT +
     2             aTHETA_COLI_MAG + iTHETA_COLI_FREQ)) *
     3             ((ONE + IC_int_COLI_dir) / (IC_int_COLI_dir)) *
     4             (((iTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt)) / (ONE +
     5             (iTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt)))))
     6             , 4)

C ColI mRNA expression for IL1B immunopositive cells
            COLI_IL1B = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_COLI_LACT * X_S_COLI_LACT) +
     2             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     3             (aTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     4             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_COMPL) +
     5             (aTHETA_COLI_IL1B * IL1B_PROT)) / (ONE +
     6             (aTHETA_COLI_LACT * X_S_COLI_LACT) +
     7             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     8             (aTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     9             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_COMPL) +
     1             (aTHETA_COLI_IL1B * IL1B_PROT)))) * (ONE -
     2             (((iTHETA_COLI_FREQ) /
     3             (aTHETA_COLI_GLUC + aTHETA_COLI_LACT +
     4             aTHETA_COLI_IL1B +
     5             aTHETA_COLI_MAG + iTHETA_COLI_FREQ)) *
     6             ((ONE + IC_int_COLI_dir) / (IC_int_COLI_dir)) *
     7             (((iTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt)) / (ONE +
     8             (iTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt)))))
     9             , 4)

C ColI mRNA expression for TNF immunopositive cells
            COLI_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_COLI_LACT * X_S_COLI_LACT) +
     2             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     3             (aTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     4             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_COMPL)) / (ONE +
     5             (aTHETA_COLI_LACT * X_S_COLI_LACT) +
     6             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     7             (aTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     8             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_COMPL)))) * (ONE -
     9             (((iTHETA_COLI_FREQ + iTHETA_COLI_TNF) /
     1             (aTHETA_COLI_GLUC + aTHETA_COLI_LACT +
     2             aTHETA_COLI_MAG + iTHETA_COLI_FREQ +
     3             iTHETA_COLI_TNF)) *
     4             ((ONE + IC_int_COLI_dir) / (IC_int_COLI_dir)) *
     5             (((iTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt) +
     6             (iTHETA_COLI_TNF * TNF_PROT)) / (ONE +
     7             (iTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt) +
     8             (iTHETA_COLI_TNF * TNF_PROT)))))
     9             , 4)

C ColI mRNA expression for IL1B and TNF immunopositive cells
            COLI_IL1B_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_COLI_LACT * X_S_COLI_LACT) +
     2             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     3             (aTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     4             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_COMPL) +
     5             (aTHETA_COLI_IL1B * IL1B_PROT)) / (ONE +
     6             (aTHETA_COLI_LACT * X_S_COLI_LACT) +
     7             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     8             (aTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     9             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_COMPL) +
     1             (aTHETA_COLI_IL1B * IL1B_PROT)))) * (ONE -
     2             (((iTHETA_COLI_FREQ +
     3             iTHETA_COLI_TNF) /
     4             (aTHETA_COLI_GLUC + aTHETA_COLI_LACT +
     5             aTHETA_COLI_IL1B +
     6             aTHETA_COLI_MAG + iTHETA_COLI_FREQ +
     7             iTHETA_COLI_TNF)) *
     8             ((ONE + IC_int_COLI_dir) / (IC_int_COLI_dir)) *
     9             (((iTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt) +
     1             (iTHETA_COLI_TNF * TNF_PROT)) / (ONE +
     2             (iTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt) +
     3             (iTHETA_COLI_TNF * TNF_PROT)))))
     4             , 4)


C---------------------------------------------------------------------------
C Magnitude >= 1 and Frequency =  0 (Magnitude catabolic, Frequency anabolic)
C---------------------------------------------------------------------------

      ELSE IF (MAG .GE. ONE .AND. FREQ .EQ. ZERO) THEN

C Magnitude Catabolic
            COLI_MAG_CAT_RANGE = I_A_MAG_CAT_2 * COLI_MAG_FACT_2
            IF (COLI_MAG_CAT_RANGE .GT. ONE) THEN
                  COLI_MAG_CAT_RANGE = ONE
            END IF
            X_S_COLI_MAG_DEPt = ZERO + COLI_MAG_CAT_RANGE

C Frequency Anabolic (Frequency == 0 Case)
            COLI_FREQ_ANA_RANGE = FREQ_0_TIME *
     1                        DBLE(TICKS) * COLI_FREQ_FACT
            IF (COLI_FREQ_ANA_RANGE .GE. FREQ_0) THEN
                  COLI_FREQ_ANA_RANGE = FREQ_0
                  X_S_COLI_FREQ_COMPL = X_S_COLI_FREQ_COMPL +
     1                                  (FREQ_0_TIME * COLI_FREQ_FACT)
                  IF (X_S_COLI_FREQ_COMPL .GT. ONE) THEN
                        X_S_COLI_FREQ_COMPL = ONE
                  END IF
            ELSE
                  X_S_COLI_FREQ_COMPL = ZERO
            END IF
            
            X_S_COLI_FREQ_DEPt = FREQ_0 - COLI_FREQ_ANA_RANGE

C ColI mRNA expression for immunonegative cells
            COLI_IMN = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_COLI_LACT * X_S_COLI_LACT) +
     2             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     3             (aTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     4             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_COMPL)) / (ONE +
     5             (aTHETA_COLI_LACT * X_S_COLI_LACT) +
     6             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     7             (aTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     8             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_COMPL)))) * (ONE -
     9             (((iTHETA_COLI_FREQ) /
     1             (aTHETA_COLI_GLUC + aTHETA_COLI_LACT +
     2             aTHETA_COLI_MAG + iTHETA_COLI_FREQ)) *
     3             ((ONE + IC_int_COLI_dir) / (IC_int_COLI_dir)) *
     4             (((iTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt)) / (ONE +
     5             (iTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt))) ))
     6             , 4)

C ColI mRNA expression for IL1B immunopositive cells
            COLI_IL1B = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_COLI_LACT * X_S_COLI_LACT) +
     2             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     3             (aTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     4             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_COMPL) +
     5             (aTHETA_COLI_IL1B * IL1B_PROT)) / (ONE +
     6             (aTHETA_COLI_LACT * X_S_COLI_LACT) +
     7             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     8             (aTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     9             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_COMPL) +
     1             (aTHETA_COLI_IL1B * IL1B_PROT)))) * (ONE -
     2             (((iTHETA_COLI_FREQ) /
     3             (aTHETA_COLI_GLUC + aTHETA_COLI_LACT +
     4             aTHETA_COLI_IL1B +
     5             aTHETA_COLI_MAG + iTHETA_COLI_FREQ)) *
     6             ((ONE + IC_int_COLI_dir) / (IC_int_COLI_dir)) *
     7             (((iTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt)) / (ONE +
     8             (iTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt)))))
     9             , 4)

C ColI mRNA expression for TNF immunopositive cells
            COLI_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_COLI_LACT * X_S_COLI_LACT) +
     2             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     3             (aTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     4             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_COMPL)) / (ONE +
     5             (aTHETA_COLI_LACT * X_S_COLI_LACT) +
     6             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     7             (aTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     8             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_COMPL)))) * (ONE -
     9             (((iTHETA_COLI_FREQ +
     1             iTHETA_COLI_TNF) /
     2             (aTHETA_COLI_GLUC + aTHETA_COLI_LACT +
     3             aTHETA_COLI_MAG + iTHETA_COLI_FREQ +
     4             iTHETA_COLI_TNF)) *
     5             ((ONE + IC_int_COLI_dir) / (IC_int_COLI_dir)) *
     6             (((iTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt) +
     7             (iTHETA_COLI_TNF * TNF_PROT)) / (ONE +
     8             (iTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt) +
     9             (iTHETA_COLI_TNF * TNF_PROT)))))
     1             , 4)

C ColI mRNA expression for IL1B and TNF immunopositive cells
            COLI_IL1B_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_COLI_LACT * X_S_COLI_LACT) +
     2             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     3             (aTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     4             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_COMPL) +
     5             (aTHETA_COLI_IL1B * IL1B_PROT)) / (ONE +
     6             (aTHETA_COLI_LACT * X_S_COLI_LACT) +
     7             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     8             (aTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     9             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_COMPL) +
     1             (aTHETA_COLI_IL1B * IL1B_PROT)))) * (ONE -
     2             (((iTHETA_COLI_FREQ +
     3             iTHETA_COLI_TNF) /
     4             (aTHETA_COLI_GLUC + aTHETA_COLI_LACT +
     5             aTHETA_COLI_IL1B + aTHETA_COLI_MAG +
     6             iTHETA_COLI_FREQ +
     7             iTHETA_COLI_TNF)) *
     8             ((ONE + IC_int_COLI_dir) / (IC_int_COLI_dir)) *
     9             (((iTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt) +
     1             (iTHETA_COLI_TNF * TNF_PROT)) / (ONE +
     2             (iTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt) +
     3             (iTHETA_COLI_TNF * TNF_PROT)))))
     4             , 4)

C---------------------------------------------------------------------------
C Magnitude >= 1.0 and Frequency > 3.0 (Magnitude & Frequency catabolic)
C---------------------------------------------------------------------------

      ELSE IF (MAG .GE. ONE .AND. FREQ .GT. 3.0D0) THEN

C Magnitude Catabolic
            COLI_MAG_CAT_RANGE = I_A_MAG_CAT_2 * COLI_MAG_FACT_2
            IF (COLI_MAG_CAT_RANGE .GT. ONE) THEN
                  COLI_MAG_CAT_RANGE = ONE
            END IF
            X_S_COLI_MAG_DEPt = ZERO + COLI_MAG_CAT_RANGE

C Frequency Catabolic
            COLI_FREQ_CAT_RANGE = I_A_FREQ_CAT_2 * COLI_FREQ_FACT_2
            IF (COLI_FREQ_CAT_RANGE .GT. ONE) THEN
                  COLI_FREQ_CAT_RANGE = ONE
            END IF
            X_S_COLI_FREQ_DEPt = ZERO + COLI_FREQ_CAT_RANGE

C ColI mRNA expression for immunonegative cells
            COLI_IMN = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_COLI_LACT * X_S_COLI_LACT) +
     2             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     3             (aTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     4             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt)) / (ONE +
     5             (aTHETA_COLI_LACT * X_S_COLI_LACT) +
     6             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     7             (aTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     8             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt))))
     9             , 4)

C ColI mRNA expression for IL1B immunopositive cells
            COLI_IL1B = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_COLI_LACT * X_S_COLI_LACT) +
     2             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     3             (aTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     4             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt) +
     5             (aTHETA_COLI_IL1B * IL1B_PROT)) / (ONE +
     6             (aTHETA_COLI_LACT * X_S_COLI_LACT) +
     7             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     8             (aTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     9             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt) +
     1             (aTHETA_COLI_IL1B * IL1B_PROT))))
     2             , 4)

C ColI mRNA expression for TNF immunopositive cells
            COLI_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_COLI_LACT * X_S_COLI_LACT) +
     2             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     3             (aTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     4             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt)) / (ONE +
     5             (aTHETA_COLI_LACT * X_S_COLI_LACT) +
     6             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     7             (aTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     8             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt)))) * (ONE -
     9             (((iTHETA_COLI_TNF) /
     1             (aTHETA_COLI_GLUC + aTHETA_COLI_LACT +
     2             aTHETA_COLI_MAG + aTHETA_COLI_FREQ +
     3             iTHETA_COLI_TNF)) *
     4             ((ONE + IC_int_COLI_dir) / (IC_int_COLI_dir)) *
     5             (((iTHETA_COLI_TNF * TNF_PROT)) / (ONE +
     6             (iTHETA_COLI_TNF * TNF_PROT)))))
     7             , 4)

C ColI mRNA expression for IL1B and TNF immunopositive cells
            COLI_IL1B_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_COLI_LACT * X_S_COLI_LACT) +
     2             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     3             (aTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     4             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt) +
     5             (aTHETA_COLI_IL1B * IL1B_PROT)) / (ONE +
     6             (aTHETA_COLI_LACT * X_S_COLI_LACT) +
     7             (aTHETA_COLI_GLUC * X_S_COLI_GLUC) +
     8             (aTHETA_COLI_MAG * X_S_COLI_MAG_DEPt) +
     9             (aTHETA_COLI_FREQ * X_S_COLI_FREQ_DEPt) +
     1             (aTHETA_COLI_IL1B * IL1B_PROT)))) * (ONE -
     2             (((iTHETA_COLI_TNF) /
     3             (aTHETA_COLI_GLUC + aTHETA_COLI_LACT +
     4             aTHETA_COLI_FREQ +
     5             aTHETA_COLI_IL1B +
     6             aTHETA_COLI_MAG +
     7             iTHETA_COLI_TNF)) *
     8             ((ONE + IC_int_COLI_dir) / (IC_int_COLI_dir)) *
     9             (((iTHETA_COLI_TNF * TNF_PROT)) / (ONE +
     1             (iTHETA_COLI_TNF * TNF_PROT)))))
     2             , 4)


      END IF

C--------------------------------------------------------------------------
      RETURN
      END SUBROUTINE PNeq_COLI

C---------------------------------------------------------------------------------------------------------------------------------------------
C Evaluate the collagen II expression
C --------------------------------------------------------------

      SUBROUTINE PNeq_COLII(COLII_IMN, COLII_IL1B,
     1           COLII_TNF, COLII_IL1B_TNF,
     2           TICKS,
     3           IL1B_PROT, TNF_PROT,
     4           GLUC_int, pH_val, MAG, FREQ,
     5           FREQ_0, FREQ_0_time,
     6           MAG_0, MAG_0_time,
     7           aTHETA_COLII_GLUC, aTHETA_COLII_LACT,
     8           iTHETA_COLII_IL1B, iTHETA_COLII_TNF,
     9           aTHETA_COLII_MAG, iTHETA_COLII_MAG,
     1           aTHETA_COLII_FREQ, iTHETA_COLII_FREQ,
     2           AC_int_IL1B_dir, AC_int_TNF_dir,
     3           IC_int_IL1B_dir, IC_int_TNF_dir,
     4           IC_int_COLII_dir,
     5           X_S_COLII_GLUC, X_S_COLII_LACT,
     6           X_S_COLII_MAG, X_S_COLII_FREQ,
     7           A_I_MAG_ANA, A_I_MAG_CAT_1, A_I_MAG_CAT_2,
     8           I_A_MAG_ANA, I_A_MAG_CAT_1, I_A_MAG_CAT_2,
     9           A_I_FREQ_ANA, A_I_FREQ_CAT_1, A_I_FREQ_CAT_2,
     1           I_A_FREQ_ANA, I_A_FREQ_CAT_1, I_A_FREQ_CAT_2,
     2           COLII_MAG_FACT, COLII_MAG_FACT_2,
     3           COLII_FREQ_FACT, COLII_FREQ_FACT_2,
     4           X_S_COLII_MAG_COMPL, X_S_COLII_FREQ_COMPL,
     5           aTHETA_SUM)

C VARIABLES:
C X_S: SENSITIVITY OF A CELL ACTIVITY (CA) TO A STIMULUS CONCENTRATION
C THETA: SENSITIVITY OF A CELL ACTIVITY (CA) TO A STIMULUS TYPE
C OUTPUTS:
C COLII_IMN: COLIIRECAN EXPRESSION FOR IMMUNONEGATIVE CELLS
C COLII_IL1B: CONTRIBUTION FROM IL1B
C COLII_TNF: CONTRIBUTION FROM TNF
C COLII_IL1B_TNF: COMBINED CONTRIBUTION OF IL1B AND TNF

C ********************************************************************
C
      IMPLICIT NONE
C
C FUNCTION
      DOUBLE PRECISION ROUND

C INTEGER VARIABLES
      INTEGER TICKS

C SPECIAL CASE FOR FREQ = 0
      DOUBLE PRECISION FREQ_0, FREQ_0_time
      DOUBLE PRECISION MAG_0, MAG_0_time

C Input and Output Variables
      DOUBLE PRECISION COLII_IMN, COLII_IL1B, COLII_TNF, COLII_IL1B_TNF
      DOUBLE PRECISION IL1B_PROT, TNF_PROT
      DOUBLE PRECISION GLUC_int, pH_val, MAG, FREQ

C Constants for Sensitivity Factors
      DOUBLE PRECISION aTHETA_COLII_GLUC, aTHETA_COLII_LACT
      DOUBLE PRECISION iTHETA_COLII_IL1B, iTHETA_COLII_TNF
      DOUBLE PRECISION aTHETA_COLII_MAG, iTHETA_COLII_MAG
      DOUBLE PRECISION aTHETA_COLII_FREQ, iTHETA_COLII_FREQ
      DOUBLE PRECISION AC_int_IL1B_dir, AC_int_TNF_dir
      DOUBLE PRECISION IC_int_IL1B_dir, IC_int_TNF_dir
      DOUBLE PRECISION IC_int_COLII_dir
      DOUBLE PRECISION aTHETA_SUM

C Intermediate Variables
      DOUBLE PRECISION X_S_COLII_GLUC, X_S_COLII_LACT
      DOUBLE PRECISION X_S_COLII_MAG, X_S_COLII_FREQ
      DOUBLE PRECISION COLII_MAG_FACT, COLII_MAG_FACT_2
      DOUBLE PRECISION COLII_FREQ_FACT, COLII_FREQ_FACT_2

C Completive factors
      DOUBLE PRECISION X_S_COLII_MAG_COMPL
      DOUBLE PRECISION X_S_COLII_MAG_DEPt
      DOUBLE PRECISION X_S_COLII_FREQ_COMPL
      DOUBLE PRECISION X_S_COLII_FREQ_DEPt

C Time dependent funtions
      DOUBLE PRECISION A_I_MAG_ANA, A_I_MAG_CAT_1, A_I_MAG_CAT_2
      DOUBLE PRECISION I_A_MAG_ANA, I_A_MAG_CAT_1, I_A_MAG_CAT_2
      DOUBLE PRECISION A_I_FREQ_ANA, A_I_FREQ_CAT_1, A_I_FREQ_CAT_2
      DOUBLE PRECISION I_A_FREQ_ANA, I_A_FREQ_CAT_1, I_A_FREQ_CAT_2

C Time sensitivity range
      DOUBLE PRECISION COLII_MAG_ANA_RANGE
      DOUBLE PRECISION COLII_FREQ_ANA_RANGE
      DOUBLE PRECISION COLII_MAG_CAT_RANGE
      DOUBLE PRECISION COLII_FREQ_CAT_RANGE

C CONSTANT PARAMETERS
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.D0, ONE = 1.D0

C---------------------------------------------------------------------------
C Initialize Variables
      COLII_MAG_ANA_RANGE = ZERO
      COLII_FREQ_ANA_RANGE = ZERO
      COLII_MAG_CAT_RANGE = ZERO
      COLII_FREQ_CAT_RANGE = ZERO
      X_S_COLII_MAG_DEPt = ZERO
      X_S_COLII_FREQ_DEPt = ZERO

C---------------------------------------------------------------------------
C Magnitude <= 0.998 and Frequency <= 3 (Magnitude anabolic, Frequency anabolic)
C---------------------------------------------------------------------------

      IF (MAG .LE. 0.998D0 .AND. FREQ .LE. 3.0D0 .AND.
     1    FREQ .NE. ZERO) THEN

C Magnitude Anabolic
            COLII_MAG_ANA_RANGE = A_I_MAG_ANA * COLII_MAG_FACT

            IF (COLII_MAG_ANA_RANGE .GE. ABS(X_S_COLII_MAG)) THEN
                  COLII_MAG_ANA_RANGE = ABS(X_S_COLII_MAG)
                  X_S_COLII_MAG_COMPL = X_S_COLII_MAG_COMPL +
     1                                  (A_I_MAG_CAT_1 * COLII_MAG_FACT)

                  IF (X_S_COLII_MAG_COMPL .GT. ONE) THEN
                        X_S_COLII_MAG_COMPL = ONE
                  END IF
            ELSE
                  X_S_COLII_MAG_COMPL = ZERO
            END IF

            X_S_COLII_MAG_DEPt = X_S_COLII_MAG - COLII_MAG_ANA_RANGE

C Frequency Anabolic
            COLII_FREQ_ANA_RANGE = A_I_FREQ_ANA * COLII_FREQ_FACT

            IF (COLII_FREQ_ANA_RANGE .GE. ABS(X_S_COLII_FREQ)) THEN
                  COLII_FREQ_ANA_RANGE = ABS(X_S_COLII_FREQ)
                  X_S_COLII_FREQ_COMPL = X_S_COLII_FREQ_COMPL +
     1                                (A_I_FREQ_CAT_1 * COLII_FREQ_FACT)

                  IF (X_S_COLII_FREQ_COMPL .GT. ONE) THEN
                        X_S_COLII_FREQ_COMPL = ONE
                  END IF
            ELSE
                  X_S_COLII_FREQ_COMPL = ZERO
            END IF

            X_S_COLII_FREQ_DEPt = X_S_COLII_FREQ - COLII_FREQ_ANA_RANGE

C ColII mRNA expression for immunonegative cells
            COLII_IMN = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_COLII_LACT * X_S_COLII_LACT) +
     2             (aTHETA_COLII_GLUC * X_S_COLII_GLUC) +
     3             (aTHETA_COLII_MAG * X_S_COLII_MAG_DEPt) +
     4             (aTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt)) / (ONE +
     5             (aTHETA_COLII_LACT * X_S_COLII_LACT) +
     6             (aTHETA_COLII_GLUC * X_S_COLII_GLUC) +
     7             (aTHETA_COLII_MAG * X_S_COLII_MAG_DEPt) +
     8             (aTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt)))) * (ONE -
     9             (((iTHETA_COLII_MAG + iTHETA_COLII_FREQ) /
     1             (aTHETA_COLII_GLUC + aTHETA_COLII_LACT +
     2             iTHETA_COLII_MAG + iTHETA_COLII_FREQ)) *
     3             ((ONE + IC_int_COLII_dir) / (IC_int_COLII_dir)) *
     4             (((iTHETA_COLII_MAG * X_S_COLII_MAG_COMPL) +
     5             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_COMPL)) / (ONE +
     6             (iTHETA_COLII_MAG * X_S_COLII_MAG_COMPL) +
     7             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_COMPL)))))
     8             , 4)

C ColII mRNA expression for IL1B immunopositive cells
            COLII_IL1B = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_COLII_LACT * X_S_COLII_LACT) +
     2             (aTHETA_COLII_GLUC * X_S_COLII_GLUC) +
     3             (aTHETA_COLII_MAG * X_S_COLII_MAG_DEPt) +
     4             (aTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt)) / (ONE +
     5             (aTHETA_COLII_LACT * X_S_COLII_LACT) +
     6             (aTHETA_COLII_GLUC * X_S_COLII_GLUC) +
     7             (aTHETA_COLII_MAG * X_S_COLII_MAG_DEPt) +
     8             (aTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt)))) * (ONE -
     9             (((iTHETA_COLII_MAG + iTHETA_COLII_FREQ +
     1             iTHETA_COLII_IL1B) /
     2             (aTHETA_COLII_GLUC + aTHETA_COLII_LACT +
     3             iTHETA_COLII_MAG + iTHETA_COLII_FREQ +
     4             iTHETA_COLII_IL1B)) *
     5             ((ONE + IC_int_COLII_dir) / (IC_int_COLII_dir)) *
     6             (((iTHETA_COLII_MAG * X_S_COLII_MAG_COMPL) +
     7             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_COMPL) +
     8             (iTHETA_COLII_IL1B * IL1B_PROT)) / (ONE +
     9             (iTHETA_COLII_MAG * X_S_COLII_MAG_COMPL) +
     1             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_COMPL) +
     2             (iTHETA_COLII_IL1B * IL1B_PROT)))))
     3             , 4)

C ColII mRNA expression for TNF immunopositive cells
            COLII_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_COLII_LACT * X_S_COLII_LACT) +
     2             (aTHETA_COLII_GLUC * X_S_COLII_GLUC) +
     3             (aTHETA_COLII_MAG * X_S_COLII_MAG_DEPt) +
     4             (aTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt)) / (ONE +
     5             (aTHETA_COLII_LACT * X_S_COLII_LACT) +
     6             (aTHETA_COLII_GLUC * X_S_COLII_GLUC) +
     7             (aTHETA_COLII_MAG * X_S_COLII_MAG_DEPt) +
     8             (aTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt)))) * (ONE -
     9             (((iTHETA_COLII_MAG + iTHETA_COLII_FREQ +
     1             iTHETA_COLII_TNF) /
     2             (aTHETA_COLII_GLUC + aTHETA_COLII_LACT +
     3             iTHETA_COLII_MAG + iTHETA_COLII_FREQ +
     4             iTHETA_COLII_TNF)) *
     5             ((ONE + IC_int_COLII_dir) / (IC_int_COLII_dir)) *
     6             (((iTHETA_COLII_MAG * X_S_COLII_MAG_COMPL) +
     7             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_COMPL) +
     8             (iTHETA_COLII_TNF * TNF_PROT)) / (ONE +
     9             (iTHETA_COLII_MAG * X_S_COLII_MAG_COMPL) +
     1             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_COMPL) +
     2             (iTHETA_COLII_TNF * TNF_PROT)))))
     3             , 4)

C ColII mRNA expression for IL1B and TNF immunopositive cells
            COLII_IL1B_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM))*
     1             (((aTHETA_COLII_LACT * X_S_COLII_LACT) +
     2             (aTHETA_COLII_GLUC * X_S_COLII_GLUC) +
     3             (aTHETA_COLII_MAG * X_S_COLII_MAG_DEPt) +
     4             (aTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt)) / (ONE +
     5             (aTHETA_COLII_LACT * X_S_COLII_LACT) +
     6             (aTHETA_COLII_GLUC * X_S_COLII_GLUC) +
     7             (aTHETA_COLII_MAG * X_S_COLII_MAG_DEPt) +
     8             (aTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt)))) * (ONE -
     9             (((iTHETA_COLII_MAG + iTHETA_COLII_FREQ +
     1             iTHETA_COLII_IL1B + iTHETA_COLII_TNF) /
     2             (aTHETA_COLII_GLUC + aTHETA_COLII_LACT +
     3             iTHETA_COLII_MAG + iTHETA_COLII_FREQ +
     4             iTHETA_COLII_IL1B + iTHETA_COLII_TNF)) *
     5             ((ONE + IC_int_COLII_dir) / (IC_int_COLII_dir)) *
     6             (((iTHETA_COLII_MAG * X_S_COLII_MAG_COMPL) +
     7             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_COMPL) +
     8             (iTHETA_COLII_IL1B * IL1B_PROT) +
     9             (iTHETA_COLII_TNF * TNF_PROT)) / (ONE +
     1             (iTHETA_COLII_MAG * X_S_COLII_MAG_COMPL) +
     2             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_COMPL) +
     3             (iTHETA_COLII_IL1B * IL1B_PROT) +
     4             (iTHETA_COLII_TNF * TNF_PROT)))))
     5             , 4)

C---------------------------------------------------------------------------
C Magnitude <= 0.998 and Frequency = 0 (Magnitude & Frequency anabolic)
C---------------------------------------------------------------------------

      ELSE IF (MAG .LE. 0.998D0 .AND. FREQ .EQ. ZERO) THEN

C Magnitude Anabolic
            COLII_MAG_ANA_RANGE = A_I_MAG_ANA * COLII_MAG_FACT

            IF (COLII_MAG_ANA_RANGE .GE. ABS(X_S_COLII_MAG)) THEN
                  COLII_MAG_ANA_RANGE = ABS(X_S_COLII_MAG)
                  X_S_COLII_MAG_COMPL = X_S_COLII_MAG_COMPL +
     1                                  (A_I_MAG_CAT_1 * COLII_MAG_FACT)

                  IF (X_S_COLII_MAG_COMPL .GT. ONE) THEN
                        X_S_COLII_MAG_COMPL = ONE
                  END IF
            ELSE
                  X_S_COLII_MAG_COMPL = ZERO
            END IF

            X_S_COLII_MAG_DEPt = X_S_COLII_MAG - COLII_MAG_ANA_RANGE

C Frequency Anabolic (Frequency == 0 Case)
            COLII_FREQ_ANA_RANGE = FREQ_0_TIME *
     1                              DBLE(TICKS) * COLII_FREQ_FACT

            IF (COLII_FREQ_ANA_RANGE .GE. FREQ_0) THEN
                  COLII_FREQ_ANA_RANGE = FREQ_0
                  X_S_COLII_FREQ_COMPL = X_S_COLII_FREQ_COMPL +
     1                                   (FREQ_0_TIME * COLII_FREQ_FACT)

                  IF (X_S_COLII_FREQ_COMPL .GT. ONE) THEN
                        X_S_COLII_FREQ_COMPL = ONE
                  END IF
            ELSE
                  X_S_COLII_FREQ_COMPL = ZERO
            END IF

            X_S_COLII_FREQ_DEPt = FREQ_0 - COLII_FREQ_ANA_RANGE

C ColII mRNA expression for immunonegative cells
            COLII_IMN = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_COLII_LACT * X_S_COLII_LACT) +
     2             (aTHETA_COLII_GLUC * X_S_COLII_GLUC) +
     3             (aTHETA_COLII_MAG * X_S_COLII_MAG_DEPt) +
     4             (aTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt)) / (ONE +
     5             (aTHETA_COLII_LACT * X_S_COLII_LACT) +
     6             (aTHETA_COLII_GLUC * X_S_COLII_GLUC) +
     7             (aTHETA_COLII_MAG * X_S_COLII_MAG_DEPt) +
     8             (aTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt)))) * (ONE -
     9             (((iTHETA_COLII_MAG + iTHETA_COLII_FREQ) /
     1             (aTHETA_COLII_GLUC + aTHETA_COLII_LACT +
     2             iTHETA_COLII_MAG + iTHETA_COLII_FREQ)) *
     3             ((ONE + IC_int_COLII_dir) / (IC_int_COLII_dir)) *
     4             (((iTHETA_COLII_MAG * X_S_COLII_MAG_COMPL) +
     5             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_COMPL)) / (ONE +
     6             (iTHETA_COLII_MAG * X_S_COLII_MAG_COMPL) +
     7             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_COMPL)))))
     8             , 4)

C ColII mRNA expression for IL1B immunopositive cells
            COLII_IL1B = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_COLII_LACT * X_S_COLII_LACT) +
     2             (aTHETA_COLII_GLUC * X_S_COLII_GLUC) +
     3             (aTHETA_COLII_MAG * X_S_COLII_MAG_DEPt) +
     4             (aTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt)) / (ONE +
     5             (aTHETA_COLII_LACT * X_S_COLII_LACT) +
     6             (aTHETA_COLII_GLUC * X_S_COLII_GLUC) +
     7             (aTHETA_COLII_MAG * X_S_COLII_MAG_DEPt) +
     8             (aTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt)))) * (ONE -
     9             (((iTHETA_COLII_MAG + iTHETA_COLII_FREQ +
     1             iTHETA_COLII_IL1B) /
     2             (aTHETA_COLII_GLUC + aTHETA_COLII_LACT +
     3             iTHETA_COLII_MAG + iTHETA_COLII_FREQ +
     4             iTHETA_COLII_IL1B)) *
     5             ((ONE + IC_int_COLII_dir) / (IC_int_COLII_dir)) *
     6             (((iTHETA_COLII_MAG * X_S_COLII_MAG_COMPL) +
     7             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_COMPL) +
     8             (iTHETA_COLII_IL1B * IL1B_PROT)) / (ONE +
     9             (iTHETA_COLII_MAG * X_S_COLII_MAG_COMPL) +
     1             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_COMPL) +
     2             (iTHETA_COLII_IL1B * IL1B_PROT)))))
     3             , 4)

C ColII mRNA expression for TNF immunopositive cells
            COLII_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_COLII_LACT * X_S_COLII_LACT) +
     2             (aTHETA_COLII_GLUC * X_S_COLII_GLUC) +
     3             (aTHETA_COLII_MAG * X_S_COLII_MAG_DEPt) +
     4             (aTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt)) / (ONE +
     5             (aTHETA_COLII_LACT * X_S_COLII_LACT) +
     6             (aTHETA_COLII_GLUC * X_S_COLII_GLUC) +
     7             (aTHETA_COLII_MAG * X_S_COLII_MAG_DEPt) +
     8             (aTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt)))) * (ONE -
     9             (((iTHETA_COLII_MAG + iTHETA_COLII_FREQ +
     1             iTHETA_COLII_TNF) /
     2             (aTHETA_COLII_GLUC + aTHETA_COLII_LACT +
     3             iTHETA_COLII_MAG + iTHETA_COLII_FREQ +
     4             iTHETA_COLII_TNF)) *
     5             ((ONE + IC_int_COLII_dir) / (IC_int_COLII_dir)) *
     6             (((iTHETA_COLII_MAG * X_S_COLII_MAG_COMPL) +
     7             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_COMPL) +
     8             (iTHETA_COLII_TNF * TNF_PROT)) / (ONE +
     9             (iTHETA_COLII_MAG * X_S_COLII_MAG_COMPL) +
     1             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_COMPL) +
     2             (iTHETA_COLII_TNF * TNF_PROT)))))
     3             , 4)

C ColII mRNA expression for IL1B and TNF immunopositive cells
            COLII_IL1B_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM))*
     1             (((aTHETA_COLII_LACT * X_S_COLII_LACT) +
     2             (aTHETA_COLII_GLUC * X_S_COLII_GLUC) +
     3             (aTHETA_COLII_MAG * X_S_COLII_MAG_DEPt) +
     4             (aTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt)) / (ONE +
     5             (aTHETA_COLII_LACT * X_S_COLII_LACT) +
     6             (aTHETA_COLII_GLUC * X_S_COLII_GLUC) +
     7             (aTHETA_COLII_MAG * X_S_COLII_MAG_DEPt) +
     8             (aTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt)))) * (ONE -
     9             (((iTHETA_COLII_MAG + iTHETA_COLII_FREQ +
     1             iTHETA_COLII_IL1B + iTHETA_COLII_TNF) /
     2             (aTHETA_COLII_GLUC + aTHETA_COLII_LACT +
     3             iTHETA_COLII_MAG + iTHETA_COLII_FREQ +
     4             iTHETA_COLII_IL1B + iTHETA_COLII_TNF)) *
     5             ((ONE + IC_int_COLII_dir) / (IC_int_COLII_dir)) *
     6             (((iTHETA_COLII_MAG * X_S_COLII_MAG_COMPL) +
     7             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_COMPL) +
     8             (iTHETA_COLII_IL1B * IL1B_PROT) +
     9             (iTHETA_COLII_TNF * TNF_PROT)) / (ONE +
     1             (iTHETA_COLII_MAG * X_S_COLII_MAG_COMPL) +
     2             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_COMPL) +
     3             (iTHETA_COLII_IL1B * IL1B_PROT) +
     4             (iTHETA_COLII_TNF * TNF_PROT)))))
     5             , 4)

C---------------------------------------------------------------------------
C Magnitude <= 0.998 and Frequency > 3 (Magnitude anabolic, Frequency catabolic)
C---------------------------------------------------------------------------

      ELSE IF (MAG .LE. 0.998D0 .AND. FREQ .GT. 3.0D0) THEN

C Magnitude Anabolic
            COLII_MAG_ANA_RANGE = A_I_MAG_ANA * COLII_MAG_FACT

            IF (COLII_MAG_ANA_RANGE .GE. ABS(X_S_COLII_MAG)) THEN
                  COLII_MAG_ANA_RANGE = ABS(X_S_COLII_MAG)
                  X_S_COLII_MAG_COMPL = X_S_COLII_MAG_COMPL +
     1                                  (A_I_MAG_CAT_1 * COLII_MAG_FACT)

                  IF (X_S_COLII_MAG_COMPL .GT. ONE) THEN
                        X_S_COLII_MAG_COMPL = ONE
                  END IF
            ELSE
                  X_S_COLII_MAG_COMPL = ZERO
            END IF

            X_S_COLII_MAG_DEPt = X_S_COLII_MAG - COLII_MAG_ANA_RANGE

C Frequency Catabolic
            COLII_FREQ_CAT_RANGE = A_I_FREQ_CAT_2 * COLII_FREQ_FACT_2

            IF (COLII_FREQ_CAT_RANGE .GT. ONE) THEN
                  COLII_FREQ_CAT_RANGE = ONE
            END IF

            X_S_COLII_FREQ_DEPt = ABS(ZERO - COLII_FREQ_CAT_RANGE)

C ColII mRNA expression for immunonegative cells
            COLII_IMN = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_COLII_LACT * X_S_COLII_LACT) +
     2             (aTHETA_COLII_GLUC * X_S_COLII_GLUC) +
     3             (aTHETA_COLII_MAG * X_S_COLII_MAG_DEPt)) / (ONE +
     4             (aTHETA_COLII_LACT * X_S_COLII_LACT) +
     5             (aTHETA_COLII_GLUC * X_S_COLII_GLUC) +
     6             (aTHETA_COLII_MAG * X_S_COLII_MAG_DEPt)))) * (ONE -
     7             (((iTHETA_COLII_MAG + iTHETA_COLII_FREQ) /
     8             (aTHETA_COLII_GLUC + aTHETA_COLII_LACT +
     9             iTHETA_COLII_MAG + iTHETA_COLII_FREQ)) *
     1             ((ONE + IC_int_COLII_dir) / (IC_int_COLII_dir)) *
     2             (((iTHETA_COLII_MAG * X_S_COLII_MAG_COMPL) +
     3             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt)) / (ONE +
     4             (iTHETA_COLII_MAG * X_S_COLII_MAG_COMPL) +
     5             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt)))))
     6             , 4)

C ColII mRNA expression for IL1B immunopositive cells
            COLII_IL1B = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_COLII_LACT * X_S_COLII_LACT) +
     2             (aTHETA_COLII_GLUC * X_S_COLII_GLUC) +
     3             (aTHETA_COLII_MAG * X_S_COLII_MAG_DEPt)) / (ONE +
     4             (aTHETA_COLII_LACT * X_S_COLII_LACT) +
     5             (aTHETA_COLII_GLUC * X_S_COLII_GLUC) +
     6             (aTHETA_COLII_MAG * X_S_COLII_MAG_DEPt)))) * (ONE -
     7             (((iTHETA_COLII_MAG + iTHETA_COLII_FREQ +
     8             iTHETA_COLII_IL1B) /
     9             (aTHETA_COLII_GLUC + aTHETA_COLII_LACT +
     1             iTHETA_COLII_MAG + iTHETA_COLII_FREQ +
     2             iTHETA_COLII_IL1B)) *
     3             ((ONE + IC_int_COLII_dir) / (IC_int_COLII_dir)) *
     4             (((iTHETA_COLII_MAG * X_S_COLII_MAG_COMPL) +
     5             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt) +
     6             (iTHETA_COLII_IL1B * IL1B_PROT)) / (ONE +
     7             (iTHETA_COLII_MAG * X_S_COLII_MAG_COMPL) +
     8             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt) +
     9             (iTHETA_COLII_IL1B * IL1B_PROT)))))
     1             , 4)

C ColII mRNA expression for TNF immunopositive cells
            COLII_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_COLII_LACT * X_S_COLII_LACT) +
     2             (aTHETA_COLII_GLUC * X_S_COLII_GLUC) +
     3             (aTHETA_COLII_MAG * X_S_COLII_MAG_DEPt)) / (ONE +
     4             (aTHETA_COLII_LACT * X_S_COLII_LACT) +
     5             (aTHETA_COLII_GLUC * X_S_COLII_GLUC) +
     6             (aTHETA_COLII_MAG * X_S_COLII_MAG_DEPt)))) * (ONE -
     7             (((iTHETA_COLII_MAG + iTHETA_COLII_FREQ +
     8             iTHETA_COLII_TNF) /
     9             (aTHETA_COLII_GLUC + aTHETA_COLII_LACT +
     1             iTHETA_COLII_MAG + iTHETA_COLII_FREQ +
     2             iTHETA_COLII_TNF)) *
     3             ((ONE + IC_int_COLII_dir) / (IC_int_COLII_dir)) *
     4             (((iTHETA_COLII_MAG * X_S_COLII_MAG_COMPL) +
     5             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt) +
     6             (iTHETA_COLII_TNF * TNF_PROT)) / (ONE +
     7             (iTHETA_COLII_MAG * X_S_COLII_MAG_COMPL) +
     8             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt) +
     9             (iTHETA_COLII_TNF * TNF_PROT)))))
     1             , 4)

C ColII mRNA expression for IL1B and TNF immunopositive cells
            COLII_IL1B_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM))*
     1             (((aTHETA_COLII_LACT * X_S_COLII_LACT) +
     2             (aTHETA_COLII_GLUC * X_S_COLII_GLUC) +
     3             (aTHETA_COLII_MAG * X_S_COLII_MAG_DEPt)) / (ONE +
     4             (aTHETA_COLII_LACT * X_S_COLII_LACT) +
     5             (aTHETA_COLII_GLUC * X_S_COLII_GLUC) +
     6             (aTHETA_COLII_MAG * X_S_COLII_MAG_DEPt)))) * (ONE -
     7             (((iTHETA_COLII_MAG + iTHETA_COLII_FREQ +
     8             iTHETA_COLII_IL1B + iTHETA_COLII_TNF) /
     9             (aTHETA_COLII_GLUC + aTHETA_COLII_LACT +
     1             iTHETA_COLII_MAG + iTHETA_COLII_FREQ +
     2             iTHETA_COLII_IL1B + iTHETA_COLII_TNF)) *
     3             ((ONE + IC_int_COLII_dir) / (IC_int_COLII_dir)) *
     4             (((iTHETA_COLII_MAG * X_S_COLII_MAG_COMPL) +
     5             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt) +
     6             (iTHETA_COLII_IL1B * IL1B_PROT) +
     7             (iTHETA_COLII_TNF * TNF_PROT)) / (ONE +
     8             (iTHETA_COLII_MAG * X_S_COLII_MAG_COMPL) +
     9             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt) +
     1             (iTHETA_COLII_IL1B * IL1B_PROT) +
     2             (iTHETA_COLII_TNF * TNF_PROT)))))
     3             , 4)

C--------------------------------------------------------------------------
C Magnitude >= 1.0 and Frequency <= 3 (Magnitude catabolic, Frequency anabolic)
C--------------------------------------------------------------------------

      ELSE IF (MAG .GE. ONE .AND. FREQ .LE. 3.0D0 .AND.
     1         FREQ .NE. ZERO) THEN

C Magnitude Catabolic
            COLII_MAG_CAT_RANGE = A_I_MAG_CAT_2 * COLII_MAG_FACT_2

            IF (COLII_MAG_CAT_RANGE .GT. ONE) THEN
                  COLII_MAG_CAT_RANGE = ONE
            END IF

            X_S_COLII_MAG_DEPt = ABS(ZERO - COLII_MAG_CAT_RANGE)

C Frequency Anabolic
            COLII_FREQ_ANA_RANGE = A_I_FREQ_ANA * COLII_FREQ_FACT

            IF (COLII_FREQ_ANA_RANGE .GE. ABS(X_S_COLII_FREQ)) THEN
                  COLII_FREQ_ANA_RANGE = ABS(X_S_COLII_FREQ)
                  X_S_COLII_FREQ_COMPL = X_S_COLII_FREQ_COMPL +
     1                                (A_I_FREQ_CAT_1 * COLII_FREQ_FACT)

                  IF (X_S_COLII_FREQ_COMPL .GT. ONE) THEN
                        X_S_COLII_FREQ_COMPL = ONE
                  END IF
            ELSE
                  X_S_COLII_FREQ_COMPL = ZERO
            END IF

            X_S_COLII_FREQ_DEPt = ABS(X_S_COLII_FREQ -
     1                            COLII_FREQ_ANA_RANGE)

C ColII mRNA expression for immunonegative cells
            COLII_IMN = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_COLII_LACT * X_S_COLII_LACT) +
     2             (aTHETA_COLII_GLUC * X_S_COLII_GLUC) +
     3             (aTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt)) / (ONE +
     4             (aTHETA_COLII_LACT * X_S_COLII_LACT) +
     5             (aTHETA_COLII_GLUC * X_S_COLII_GLUC) +
     6             (aTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt)))) * (ONE -
     7             (((iTHETA_COLII_MAG + iTHETA_COLII_FREQ) /
     8             (aTHETA_COLII_GLUC + aTHETA_COLII_LACT +
     9             iTHETA_COLII_MAG + iTHETA_COLII_FREQ)) *
     1             ((ONE + IC_int_COLII_dir) / (IC_int_COLII_dir)) *
     2             (((iTHETA_COLII_MAG * X_S_COLII_MAG_DEPt) +
     3             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_COMPL)) / (ONE +
     4             (iTHETA_COLII_MAG * X_S_COLII_MAG_DEPt) +
     5             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_COMPL)))))
     6             , 4)

C ColII mRNA expression for IL1B immunopositive cells
            COLII_IL1B = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_COLII_LACT * X_S_COLII_LACT) +
     2             (aTHETA_COLII_GLUC * X_S_COLII_GLUC) +
     3             (aTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt)) / (ONE +
     4             (aTHETA_COLII_LACT * X_S_COLII_LACT) +
     5             (aTHETA_COLII_GLUC * X_S_COLII_GLUC) +
     6             (aTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt)))) * (ONE -
     7             (((iTHETA_COLII_MAG + iTHETA_COLII_FREQ +
     8             iTHETA_COLII_IL1B) /
     9             (aTHETA_COLII_GLUC + aTHETA_COLII_LACT +
     1             iTHETA_COLII_MAG + iTHETA_COLII_FREQ +
     2             iTHETA_COLII_IL1B)) *
     3             ((ONE + IC_int_COLII_dir) / (IC_int_COLII_dir)) *
     4             (((iTHETA_COLII_MAG * X_S_COLII_MAG_DEPt) +
     5             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_COMPL) +
     6             (iTHETA_COLII_IL1B * IL1B_PROT)) / (ONE +
     7             (iTHETA_COLII_MAG * X_S_COLII_MAG_DEPt) +
     8             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_COMPL) +
     9             (iTHETA_COLII_IL1B * IL1B_PROT)))))
     1             , 4)

C ColII mRNA expression for TNF immunopositive cells
            COLII_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_COLII_LACT * X_S_COLII_LACT) +
     2             (aTHETA_COLII_GLUC * X_S_COLII_GLUC) +
     3             (aTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt)) / (ONE +
     4             (aTHETA_COLII_LACT * X_S_COLII_LACT) +
     5             (aTHETA_COLII_GLUC * X_S_COLII_GLUC) +
     6             (aTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt)))) * (ONE -
     7             (((iTHETA_COLII_MAG + iTHETA_COLII_FREQ +
     8             iTHETA_COLII_TNF) /
     9             (aTHETA_COLII_GLUC + aTHETA_COLII_LACT +
     1             iTHETA_COLII_MAG + iTHETA_COLII_FREQ +
     2             iTHETA_COLII_TNF)) *
     3             ((ONE + IC_int_COLII_dir) / (IC_int_COLII_dir)) *
     4             (((iTHETA_COLII_MAG * X_S_COLII_MAG_DEPt) +
     5             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_COMPL) +
     6             (iTHETA_COLII_TNF * TNF_PROT)) / (ONE +
     7             (iTHETA_COLII_MAG * X_S_COLII_MAG_DEPt) +
     8             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_COMPL) +
     9             (iTHETA_COLII_TNF * TNF_PROT)))))
     1             , 4)

C ColII mRNA expression for IL1B and TNF immunopositive cells
            COLII_IL1B_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM))*
     1             (((aTHETA_COLII_LACT * X_S_COLII_LACT) +
     2             (aTHETA_COLII_GLUC * X_S_COLII_GLUC) +
     3             (aTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt)) / (ONE +
     4             (aTHETA_COLII_LACT * X_S_COLII_LACT) +
     5             (aTHETA_COLII_GLUC * X_S_COLII_GLUC) +
     6             (aTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt)))) * (ONE -
     7             (((iTHETA_COLII_MAG + iTHETA_COLII_FREQ +
     8             iTHETA_COLII_IL1B + iTHETA_COLII_TNF) /
     9             (aTHETA_COLII_GLUC + aTHETA_COLII_LACT +
     1             iTHETA_COLII_MAG + iTHETA_COLII_FREQ +
     2             iTHETA_COLII_IL1B + iTHETA_COLII_TNF)) *
     3             ((ONE + IC_int_COLII_dir) / (IC_int_COLII_dir)) *
     4             (((iTHETA_COLII_MAG * X_S_COLII_MAG_DEPt) +
     5             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_COMPL) +
     6             (iTHETA_COLII_IL1B * IL1B_PROT) +
     7             (iTHETA_COLII_TNF * TNF_PROT)) / (ONE +
     8             (iTHETA_COLII_MAG * X_S_COLII_MAG_DEPt) +
     9             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_COMPL) +
     1             (iTHETA_COLII_IL1B * IL1B_PROT) +
     2             (iTHETA_COLII_TNF * TNF_PROT)))))
     3             , 4)

C--------------------------------------------------------------------------
C Magnitude >= 1.0 and Frequency = 0 (Magnitude & Frequency catabolic)
C--------------------------------------------------------------------------

      ELSE IF (MAG .GE. ONE .AND. FREQ .EQ. ZERO) THEN

C Magnitude Catabolic
            COLII_MAG_CAT_RANGE = A_I_MAG_CAT_2 * COLII_MAG_FACT_2
            IF (COLII_MAG_CAT_RANGE .GT. ONE) THEN
                  COLII_MAG_CAT_RANGE = ONE
            END IF
            X_S_COLII_MAG_DEPt = ABS(ZERO - COLII_MAG_CAT_RANGE)

C Frequency Anabolic (Frequency == 0 Case)
            COLII_FREQ_ANA_RANGE = FREQ_0_TIME *
     1                             DBLE(TICKS) * COLII_FREQ_FACT
            IF (COLII_FREQ_ANA_RANGE .GE. FREQ_0) THEN
                  COLII_FREQ_ANA_RANGE = FREQ_0
                  X_S_COLII_FREQ_COMPL = X_S_COLII_FREQ_COMPL +
     1                                  (FREQ_0_TIME * COLII_FREQ_FACT)
                  IF (X_S_COLII_FREQ_COMPL .GT. ONE) THEN
                        X_S_COLII_FREQ_COMPL = ONE
                  END IF
            ELSE
                  X_S_COLII_FREQ_COMPL = ZERO
            END IF
            X_S_COLII_FREQ_DEPt = FREQ_0 - COLII_FREQ_ANA_RANGE

C ColII mRNA expression for immunonegative cells
            COLII_IMN = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_COLII_LACT * X_S_COLII_LACT) +
     2             (aTHETA_COLII_GLUC * X_S_COLII_GLUC) +
     3             (aTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt)) / (ONE +
     4             (aTHETA_COLII_LACT * X_S_COLII_LACT) +
     5             (aTHETA_COLII_GLUC * X_S_COLII_GLUC) +
     6             (aTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt)))) * (ONE -
     7             (((iTHETA_COLII_MAG + iTHETA_COLII_FREQ) /
     8             (aTHETA_COLII_GLUC + aTHETA_COLII_LACT +
     9             iTHETA_COLII_MAG + iTHETA_COLII_FREQ)) *
     1             ((ONE + IC_int_COLII_dir) / (IC_int_COLII_dir)) *
     2             (((iTHETA_COLII_MAG * X_S_COLII_MAG_DEPt) +
     3             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_COMPL)) / (ONE +
     4             (iTHETA_COLII_MAG * X_S_COLII_MAG_DEPt) +
     5             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_COMPL)))))
     6             , 4)

C ColII mRNA expression for IL1B immunopositive cells
            COLII_IL1B = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_COLII_LACT * X_S_COLII_LACT) +
     2             (aTHETA_COLII_GLUC * X_S_COLII_GLUC) +
     3             (aTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt)) / (ONE +
     4             (aTHETA_COLII_LACT * X_S_COLII_LACT) +
     5             (aTHETA_COLII_GLUC * X_S_COLII_GLUC) +
     6             (aTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt)))) * (ONE -
     7             (((iTHETA_COLII_MAG + iTHETA_COLII_FREQ +
     8             iTHETA_COLII_IL1B) /
     9             (aTHETA_COLII_GLUC + aTHETA_COLII_LACT +
     1             iTHETA_COLII_MAG + iTHETA_COLII_FREQ +
     2             iTHETA_COLII_IL1B)) *
     3             ((ONE + IC_int_COLII_dir) / (IC_int_COLII_dir)) *
     4             (((iTHETA_COLII_MAG * X_S_COLII_MAG_DEPt) +
     5             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_COMPL) +
     6             (iTHETA_COLII_IL1B * IL1B_PROT)) / (ONE +
     7             (iTHETA_COLII_MAG * X_S_COLII_MAG_DEPt) +
     8             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_COMPL) +
     9             (iTHETA_COLII_IL1B * IL1B_PROT)))))
     1             , 4)

C ColII mRNA expression for TNF immunopositive cells
            COLII_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_COLII_LACT * X_S_COLII_LACT) +
     2             (aTHETA_COLII_GLUC * X_S_COLII_GLUC) +
     3             (aTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt)) / (ONE +
     4             (aTHETA_COLII_LACT * X_S_COLII_LACT) +
     5             (aTHETA_COLII_GLUC * X_S_COLII_GLUC) +
     6             (aTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt)))) * (ONE -
     7             (((iTHETA_COLII_MAG + iTHETA_COLII_FREQ +
     8             iTHETA_COLII_TNF) /
     9             (aTHETA_COLII_GLUC + aTHETA_COLII_LACT +
     1             iTHETA_COLII_MAG + iTHETA_COLII_FREQ +
     2             iTHETA_COLII_TNF)) *
     3             ((ONE + IC_int_COLII_dir) / (IC_int_COLII_dir)) *
     4             (((iTHETA_COLII_MAG * X_S_COLII_MAG_DEPt) +
     5             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_COMPL) +
     6             (iTHETA_COLII_TNF * TNF_PROT)) / (ONE +
     7             (iTHETA_COLII_MAG * X_S_COLII_MAG_DEPt) +
     8             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_COMPL) +
     9             (iTHETA_COLII_TNF * TNF_PROT)))))
     1             , 4)

C ColII mRNA expression for IL1B and TNF immunopositive cells
            COLII_IL1B_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM))*
     1             (((aTHETA_COLII_LACT * X_S_COLII_LACT) +
     2             (aTHETA_COLII_GLUC * X_S_COLII_GLUC) +
     3             (aTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt)) / (ONE +
     4             (aTHETA_COLII_LACT * X_S_COLII_LACT) +
     5             (aTHETA_COLII_GLUC * X_S_COLII_GLUC) +
     6             (aTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt)))) * (ONE -
     7             (((iTHETA_COLII_MAG + iTHETA_COLII_FREQ +
     8             iTHETA_COLII_IL1B + iTHETA_COLII_TNF) /
     9             (aTHETA_COLII_GLUC + aTHETA_COLII_LACT +
     1             iTHETA_COLII_MAG + iTHETA_COLII_FREQ +
     2             iTHETA_COLII_IL1B + iTHETA_COLII_TNF)) *
     3             ((ONE + IC_int_COLII_dir) / (IC_int_COLII_dir)) *
     4             (((iTHETA_COLII_MAG * X_S_COLII_MAG_DEPt) +
     5             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_COMPL) +
     6             (iTHETA_COLII_IL1B * IL1B_PROT) +
     7             (iTHETA_COLII_TNF * TNF_PROT)) / (ONE +
     8             (iTHETA_COLII_MAG * X_S_COLII_MAG_DEPt) +
     9             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_COMPL) +
     1             (iTHETA_COLII_IL1B * IL1B_PROT) +
     2             (iTHETA_COLII_TNF * TNF_PROT)))))
     3             , 4)

C--------------------------------------------------------------------------
C Magnitude >= 1.0 and Frequency > 3 (Magnitude catabolic, Frequency catabolic)
C--------------------------------------------------------------------------

      ELSE IF (MAG .GE. ONE .AND. FREQ .GT. 3.0D0) THEN

C Magnitude Catabolic
         COLII_MAG_CAT_RANGE = A_I_MAG_CAT_2 * COLII_MAG_FACT_2
         IF (COLII_MAG_CAT_RANGE .GT. ONE) THEN
            COLII_MAG_CAT_RANGE = ONE
         END IF
         X_S_COLII_MAG_DEPt = ABS(ZERO - COLII_MAG_CAT_RANGE)

C Frequency Catabolic
         COLII_FREQ_CAT_RANGE = A_I_FREQ_CAT_2 * COLII_FREQ_FACT_2
         IF (COLII_FREQ_CAT_RANGE .GT. ONE) THEN
            COLII_FREQ_CAT_RANGE = ONE
         END IF
         X_S_COLII_FREQ_DEPt = ABS(ZERO - COLII_FREQ_CAT_RANGE)

C ColII mRNA expression for immunonegative cells
         COLII_IMN = ROUND((((ONE + aTHETA_SUM) / aTHETA_SUM) *
     1             (((aTHETA_COLII_LACT * X_S_COLII_LACT) +
     2             (aTHETA_COLII_GLUC * X_S_COLII_GLUC))) / (ONE +
     3             (aTHETA_COLII_LACT * X_S_COLII_LACT) +
     4             (aTHETA_COLII_GLUC * X_S_COLII_GLUC))) * (ONE -
     5             (((iTHETA_COLII_MAG + iTHETA_COLII_FREQ) /
     6             (aTHETA_COLII_GLUC + aTHETA_COLII_LACT +
     7             iTHETA_COLII_MAG + iTHETA_COLII_FREQ)) *
     8             ((ONE + IC_int_COLII_dir) / IC_int_COLII_dir) *
     9             (((iTHETA_COLII_MAG * X_S_COLII_MAG_DEPt) +
     1             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt)) / (ONE +
     2             (iTHETA_COLII_MAG * X_S_COLII_MAG_DEPt) +
     3             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt))))), 4)

C ColII mRNA expression for IL1B immunopositive cells
         COLII_IL1B = ROUND((((ONE + aTHETA_SUM) / aTHETA_SUM) *
     1             (((aTHETA_COLII_LACT * X_S_COLII_LACT) +
     2             (aTHETA_COLII_GLUC * X_S_COLII_GLUC))) / (ONE +
     3             (aTHETA_COLII_LACT * X_S_COLII_LACT) +
     4             (aTHETA_COLII_GLUC * X_S_COLII_GLUC))) * (ONE -
     5             (((iTHETA_COLII_MAG + iTHETA_COLII_FREQ +
     6             iTHETA_COLII_IL1B) /
     7             (aTHETA_COLII_GLUC + aTHETA_COLII_LACT +
     8             iTHETA_COLII_MAG + iTHETA_COLII_FREQ +
     9             iTHETA_COLII_IL1B)) *
     1             ((ONE + IC_int_COLII_dir) / IC_int_COLII_dir) *
     2             (((iTHETA_COLII_MAG * X_S_COLII_MAG_DEPt) +
     3             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt) +
     4             (iTHETA_COLII_IL1B * IL1B_PROT)) / (ONE +
     5             (iTHETA_COLII_MAG * X_S_COLII_MAG_DEPt) +
     6             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt) +
     7             (iTHETA_COLII_IL1B * IL1B_PROT))))), 4)

C ColII mRNA expression for TNF immunopositive cells
         COLII_TNF = ROUND((((ONE + aTHETA_SUM) / aTHETA_SUM) *
     1             (((aTHETA_COLII_LACT * X_S_COLII_LACT) +
     2             (aTHETA_COLII_GLUC * X_S_COLII_GLUC))) / (ONE +
     3             (aTHETA_COLII_LACT * X_S_COLII_LACT) +
     4             (aTHETA_COLII_GLUC * X_S_COLII_GLUC))) * (ONE -
     5             (((iTHETA_COLII_MAG + iTHETA_COLII_FREQ +
     6             iTHETA_COLII_TNF) /
     7             (aTHETA_COLII_GLUC + aTHETA_COLII_LACT +
     8             iTHETA_COLII_MAG + iTHETA_COLII_FREQ +
     9             iTHETA_COLII_TNF)) *
     1             ((ONE + IC_int_COLII_dir) / IC_int_COLII_dir) *
     2             (((iTHETA_COLII_MAG * X_S_COLII_MAG_DEPt) +
     3             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt) +
     4             (iTHETA_COLII_TNF * TNF_PROT)) / (ONE +
     5             (iTHETA_COLII_MAG * X_S_COLII_MAG_DEPt) +
     6             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt) +
     7             (iTHETA_COLII_TNF * TNF_PROT))))), 4)

C ColII mRNA expression for IL1B and TNF immunopositive cells
         COLII_IL1B_TNF = ROUND((((ONE + aTHETA_SUM) / aTHETA_SUM) *
     1             (((aTHETA_COLII_LACT * X_S_COLII_LACT) +
     2             (aTHETA_COLII_GLUC * X_S_COLII_GLUC))) / (ONE +
     3             (aTHETA_COLII_LACT * X_S_COLII_LACT) +
     4             (aTHETA_COLII_GLUC * X_S_COLII_GLUC))) * (ONE -
     5             (((iTHETA_COLII_MAG + iTHETA_COLII_FREQ +
     6             iTHETA_COLII_IL1B + iTHETA_COLII_TNF) /
     7             (aTHETA_COLII_GLUC + aTHETA_COLII_LACT +
     8             iTHETA_COLII_MAG + iTHETA_COLII_FREQ +
     9             iTHETA_COLII_IL1B + iTHETA_COLII_TNF)) *
     1             ((ONE + IC_int_COLII_dir) / IC_int_COLII_dir) *
     2             (((iTHETA_COLII_MAG * X_S_COLII_MAG_DEPt) +
     3             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt) +
     4             (iTHETA_COLII_IL1B * IL1B_PROT) +
     5             (iTHETA_COLII_TNF * TNF_PROT)) / (ONE +
     6             (iTHETA_COLII_MAG * X_S_COLII_MAG_DEPt) +
     7             (iTHETA_COLII_FREQ * X_S_COLII_FREQ_DEPt) +
     8             (iTHETA_COLII_IL1B * IL1B_PROT) +
     9             (iTHETA_COLII_TNF * TNF_PROT))))), 4)

      END IF

C--------------------------------------------------------------------------
      RETURN
      END SUBROUTINE PNeq_COLII

C--------------------------------------------------------------------------------------------------------------------------------------------- 
C Evaluate the MMP-3 expression
C ---------------------------------------------------------

      SUBROUTINE PNeq_MMP3(MMP3_IMN, MMP3_IL1B,
     1           MMP3_TNF, MMP3_IL1B_TNF,
     2           TICKS,
     3           IL1B_PROT, TNF_PROT,
     4           GLUC_int, pH_val, MAG, FREQ,
     5           FREQ_0, FREQ_0_time,
     6           MAG_0, MAG_0_time,
     7           aTHETA_MMP3_GLUC, aTHETA_MMP3_LACT,
     8           aTHETA_MMP3_IL1B, aTHETA_MMP3_TNF,
     9           aTHETA_MMP3_MAG, iTHETA_MMP3_MAG,
     1           aTHETA_MMP3_FREQ, iTHETA_MMP3_FREQ,
     2           AC_int_IL1B_dir, AC_int_TNF_dir,
     3           IC_int_IL1B_dir, IC_int_TNF_dir,
     4           IC_int_MMP3_dir,
     5           X_S_MMP3_GLUC, X_S_MMP3_LACT,
     6           X_S_MMP3_MAG, X_S_MMP3_FREQ,
     7           A_I_MAG_ANA, A_I_MAG_CAT_1, A_I_MAG_CAT_2,
     8           I_A_MAG_ANA, I_A_MAG_CAT_1, I_A_MAG_CAT_2,
     9           A_I_FREQ_ANA, A_I_FREQ_CAT_1, A_I_FREQ_CAT_2,
     1           I_A_FREQ_ANA, I_A_FREQ_CAT_1, I_A_FREQ_CAT_2,
     2           MMP3_MAG_FACT, MMP3_MAG_FACT_2,
     3           MMP3_FREQ_FACT, MMP3_FREQ_FACT_2,
     4           X_S_MMP3_MAG_COMPL, X_S_MMP3_FREQ_COMPL,
     5           aTHETA_SUM)

C VARIABLES:
C X_S: SENSITIVITY OF A CELL ACTIVITY (CA) TO A STIMULUS CONCENTRATION
C THETA: SENSITIVITY OF A CELL ACTIVITY (CA) TO A STIMULUS TYPE
C OUTPUTS:
C MMP3_IMN: MMP3RECAN EXPRESSION FOR IMMUNONEGATIVE CELLS
C MMP3_IL1B: CONTRIBUTION FROM IL1B
C MMP3_TNF: CONTRIBUTION FROM TNF
C MMP3_IL1B_TNF: COMBINED CONTRIBUTION OF IL1B AND TNF

C ********************************************************************
C
      IMPLICIT NONE
C
C FUNCTION
      DOUBLE PRECISION ROUND

C INTEGER VARIABLES
      INTEGER TICKS

C SPECIAL CASE FOR FREQ = 0
      DOUBLE PRECISION FREQ_0, FREQ_0_time
      DOUBLE PRECISION MAG_0, MAG_0_time

C Input and Output Variables
      DOUBLE PRECISION MMP3_IMN, MMP3_IL1B, MMP3_TNF, MMP3_IL1B_TNF
      DOUBLE PRECISION IL1B_PROT, TNF_PROT
      DOUBLE PRECISION GLUC_int, pH_val, MAG, FREQ

C Constants for Sensitivity Factors
      DOUBLE PRECISION aTHETA_MMP3_GLUC, aTHETA_MMP3_LACT
      DOUBLE PRECISION aTHETA_MMP3_IL1B, aTHETA_MMP3_TNF
      DOUBLE PRECISION aTHETA_MMP3_MAG, iTHETA_MMP3_MAG
      DOUBLE PRECISION aTHETA_MMP3_FREQ, iTHETA_MMP3_FREQ
      DOUBLE PRECISION AC_int_IL1B_dir, AC_int_TNF_dir
      DOUBLE PRECISION IC_int_IL1B_dir, IC_int_TNF_dir
      DOUBLE PRECISION IC_int_MMP3_dir
      DOUBLE PRECISION aTHETA_SUM

C Intermediate Variables
      DOUBLE PRECISION X_S_MMP3_GLUC, X_S_MMP3_LACT
      DOUBLE PRECISION X_S_MMP3_MAG, X_S_MMP3_FREQ
      DOUBLE PRECISION MMP3_MAG_FACT, MMP3_MAG_FACT_2
      DOUBLE PRECISION MMP3_FREQ_FACT, MMP3_FREQ_FACT_2

C Completive factors
      DOUBLE PRECISION X_S_MMP3_MAG_COMPL
      DOUBLE PRECISION X_S_MMP3_MAG_DEPt
      DOUBLE PRECISION X_S_MMP3_FREQ_COMPL
      DOUBLE PRECISION X_S_MMP3_FREQ_DEPt

C Time dependent funtions
      DOUBLE PRECISION A_I_MAG_ANA, A_I_MAG_CAT_1, A_I_MAG_CAT_2
      DOUBLE PRECISION I_A_MAG_ANA, I_A_MAG_CAT_1, I_A_MAG_CAT_2
      DOUBLE PRECISION A_I_FREQ_ANA, A_I_FREQ_CAT_1, A_I_FREQ_CAT_2
      DOUBLE PRECISION I_A_FREQ_ANA, I_A_FREQ_CAT_1, I_A_FREQ_CAT_2

C Time sensitivity range
      DOUBLE PRECISION MMP3_MAG_ANA_RANGE
      DOUBLE PRECISION MMP3_FREQ_ANA_RANGE
      DOUBLE PRECISION MMP3_MAG_CAT_RANGE
      DOUBLE PRECISION MMP3_FREQ_CAT_RANGE

C CONSTANT PARAMETERS
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.D0, ONE = 1.D0

C---------------------------------------------------------------------------
C Initialize Variables
      MMP3_MAG_ANA_RANGE = ZERO
      MMP3_FREQ_ANA_RANGE = ZERO
      MMP3_MAG_CAT_RANGE = ZERO
      MMP3_FREQ_CAT_RANGE = ZERO
      X_S_MMP3_MAG_DEPt = ZERO
      X_S_MMP3_FREQ_DEPt = ZERO

C---------------------------------------------------------------------------
C Magnitude <= 0.998 and Frequency <= 3.0 (Magnitude anabolic, Frequency anabolic)
C---------------------------------------------------------------------------

      IF (MAG .LE. 0.998D0 .AND. FREQ .LE. 3.0D0 .AND.
     1    FREQ .NE. ZERO) THEN

C Magnitude Anabolic
         MMP3_MAG_ANA_RANGE = I_A_MAG_ANA * MMP3_MAG_FACT
         IF (MMP3_MAG_ANA_RANGE .GE. ABS(X_S_MMP3_MAG)) THEN
            MMP3_MAG_ANA_RANGE = ABS(X_S_MMP3_MAG)
            X_S_MMP3_MAG_COMPL = X_S_MMP3_MAG_COMPL +
     1                          (I_A_MAG_CAT_1 * MMP3_MAG_FACT)
            IF (X_S_MMP3_MAG_COMPL .GT. ONE) THEN
               X_S_MMP3_MAG_COMPL = ONE
            END IF
         ELSE
            X_S_MMP3_MAG_COMPL = ZERO
         END IF
         X_S_MMP3_MAG_DEPt = ABS(X_S_MMP3_MAG + MMP3_MAG_ANA_RANGE)

C Frequency Anabolic
         MMP3_FREQ_ANA_RANGE = I_A_FREQ_ANA * MMP3_FREQ_FACT
         IF (MMP3_FREQ_ANA_RANGE .GE. ABS(X_S_MMP3_FREQ)) THEN
            MMP3_FREQ_ANA_RANGE = ABS(X_S_MMP3_FREQ)
            X_S_MMP3_FREQ_COMPL = X_S_MMP3_FREQ_COMPL +
     1                            (I_A_FREQ_CAT_1 * MMP3_FREQ_FACT)
            IF (X_S_MMP3_FREQ_COMPL .GT. ONE) THEN
               X_S_MMP3_FREQ_COMPL = ONE
            END IF
         ELSE
            X_S_MMP3_FREQ_COMPL = ZERO
         END IF
         X_S_MMP3_FREQ_DEPt = ABS(X_S_MMP3_FREQ + MMP3_FREQ_ANA_RANGE)

C MMP3 Expression for Immunonegative Cells
         MMP3_IMN = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     2             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     3             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_COMPL) +
     4             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_COMPL)) / (ONE +
     5             (aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     6             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     7             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_COMPL) +
     8             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_COMPL)))) * (ONE -
     9             (((iTHETA_MMP3_MAG + iTHETA_MMP3_FREQ) /
     1             (aTHETA_MMP3_GLUC + aTHETA_MMP3_LACT +
     2             iTHETA_MMP3_MAG + iTHETA_MMP3_FREQ)) *
     3             ((ONE + IC_int_MMP3_dir) / (IC_int_MMP3_dir)) *
     4             (((iTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt) +
     5             (iTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt)) / (ONE +
     6             (iTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt) +
     7             (iTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt)))))
     8             , 4)

C MMP3 Expression for IL1B Immunopositive Cells
         MMP3_IL1B = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     2             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     3             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_COMPL) +
     4             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_COMPL) +
     5             (aTHETA_MMP3_IL1B * IL1B_PROT)) / (ONE +
     6             (aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     7             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     8             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_COMPL) +
     9             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_COMPL) +
     1             (aTHETA_MMP3_IL1B * IL1B_PROT)))) * (ONE -
     2             (((iTHETA_MMP3_MAG + iTHETA_MMP3_FREQ) /
     3             (aTHETA_MMP3_GLUC + aTHETA_MMP3_LACT +
     4             aTHETA_MMP3_IL1B +
     5             iTHETA_MMP3_MAG + iTHETA_MMP3_FREQ)) *
     6             ((ONE + IC_int_MMP3_dir) / (IC_int_MMP3_dir)) *
     7             (((iTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt) +
     8             (iTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt)) / (ONE +
     9             (iTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt) +
     1             (iTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt)))))
     2             , 4)

C MMP3 Expression for TNF Immunopositive Cells
         MMP3_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     2             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     3             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_COMPL) +
     4             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_COMPL) +
     5             (aTHETA_MMP3_TNF * TNF_PROT)) / (ONE +
     6             (aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     7             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     8             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_COMPL) +
     9             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_COMPL) +
     1             (aTHETA_MMP3_TNF * TNF_PROT)))) * (ONE -
     2             (((iTHETA_MMP3_MAG + iTHETA_MMP3_FREQ) /
     3             (aTHETA_MMP3_GLUC + aTHETA_MMP3_LACT +
     4             aTHETA_MMP3_TNF +
     5             iTHETA_MMP3_MAG + iTHETA_MMP3_FREQ)) *
     6             ((ONE + IC_int_MMP3_dir) / (IC_int_MMP3_dir)) *
     7             (((iTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt) +
     8             (iTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt)) / (ONE +
     9             (iTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt) +
     1             (iTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt)))))
     2             , 4)

C MMP3 Expression for IL1B & TNF Immunopositive Cells
         MMP3_IL1B_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     2             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     3             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_COMPL) +
     4             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_COMPL) +
     5             (aTHETA_MMP3_IL1B * IL1B_PROT) +
     6             (aTHETA_MMP3_TNF * TNF_PROT)) / (ONE +
     7             (aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     8             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     9             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_COMPL) +
     1             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_COMPL) +
     2             (aTHETA_MMP3_IL1B * IL1B_PROT) +
     3             (aTHETA_MMP3_TNF * TNF_PROT)))) * (ONE -
     4             (((iTHETA_MMP3_MAG + iTHETA_MMP3_FREQ) /
     5             (aTHETA_MMP3_GLUC + aTHETA_MMP3_LACT +
     6             aTHETA_MMP3_IL1B + aTHETA_MMP3_TNF +
     7             iTHETA_MMP3_MAG + iTHETA_MMP3_FREQ)) *
     8             ((ONE + IC_int_MMP3_dir) / (IC_int_MMP3_dir)) *
     9             (((iTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt) +
     1             (iTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt)) / (ONE +
     2             (iTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt) +
     3             (iTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt)))))
     4             , 4)

C--------------------------------------------------------------------------
C Magnitude <= 0.998 and Frequency == 0 (Magnitude & Frequency anabolic)
C--------------------------------------------------------------------------

      ELSE IF (MAG .LE. 0.998D0 .AND. FREQ .EQ. ZERO) THEN

C Magnitude Anabolic
         MMP3_MAG_ANA_RANGE = I_A_MAG_ANA * MMP3_MAG_FACT
         IF (MMP3_MAG_ANA_RANGE .GE. ABS(X_S_MMP3_MAG)) THEN
            MMP3_MAG_ANA_RANGE = ABS(X_S_MMP3_MAG)
            X_S_MMP3_MAG_COMPL = X_S_MMP3_MAG_COMPL +
     1                          (I_A_MAG_CAT_1 * MMP3_MAG_FACT)
            IF (X_S_MMP3_MAG_COMPL .GT. ONE) THEN
               X_S_MMP3_MAG_COMPL = ONE
            END IF
         ELSE
            X_S_MMP3_MAG_COMPL = ZERO
         END IF
         X_S_MMP3_MAG_DEPt = ABS(X_S_MMP3_MAG + MMP3_MAG_ANA_RANGE)

C Frequency Anabolic (Frequency == 0 Case)
         MMP3_FREQ_ANA_RANGE = FREQ_0_TIME *
     1                         DBLE(TICKS) * MMP3_FREQ_FACT
         IF (MMP3_FREQ_ANA_RANGE .GE. FREQ_0) THEN
            MMP3_FREQ_ANA_RANGE = FREQ_0
            X_S_MMP3_FREQ_COMPL = X_S_MMP3_FREQ_COMPL +
     1                           (FREQ_0_TIME * MMP3_FREQ_FACT)
            IF (X_S_MMP3_FREQ_COMPL .GT. ONE) THEN
               X_S_MMP3_FREQ_COMPL = ONE
            END IF
         ELSE
            X_S_MMP3_FREQ_COMPL = ZERO
         END IF
         X_S_MMP3_FREQ_DEPt = FREQ_0 - MMP3_FREQ_ANA_RANGE

C MMP3 Expression for Immunonegative Cells
         MMP3_IMN = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     2             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     3             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_COMPL) +
     4             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_COMPL)) / (ONE +
     5             (aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     6             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     7             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_COMPL) +
     8             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_COMPL)))) * (ONE -
     9             (((iTHETA_MMP3_MAG + iTHETA_MMP3_FREQ) /
     1             (aTHETA_MMP3_GLUC + aTHETA_MMP3_LACT +
     2             iTHETA_MMP3_MAG + iTHETA_MMP3_FREQ)) *
     3             ((ONE + IC_int_MMP3_dir) / (IC_int_MMP3_dir)) *
     4             (((iTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt) +
     5             (iTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt)) / (ONE +
     6             (iTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt) +
     7             (iTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt)))))
     8             , 4)

C MMP3 Expression for IL1B Immunopositive Cells
         MMP3_IL1B = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     2             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     3             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_COMPL) +
     4             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_COMPL) +
     5             (aTHETA_MMP3_IL1B * IL1B_PROT)) / (ONE +
     6             (aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     7             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     8             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_COMPL) +
     9             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_COMPL) +
     1             (aTHETA_MMP3_IL1B * IL1B_PROT)))) * (ONE -
     2             (((iTHETA_MMP3_MAG + iTHETA_MMP3_FREQ) /
     3             (aTHETA_MMP3_GLUC + aTHETA_MMP3_LACT +
     4             aTHETA_MMP3_IL1B +
     5             iTHETA_MMP3_MAG + iTHETA_MMP3_FREQ)) *
     6             ((ONE + IC_int_MMP3_dir) / (IC_int_MMP3_dir)) *
     7             (((iTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt) +
     8             (iTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt)) / (ONE +
     9             (iTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt) +
     1             (iTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt)))))
     2             , 4)

C MMP3 Expression for TNF Immunopositive Cells
         MMP3_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     2             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     3             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_COMPL) +
     4             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_COMPL) +
     5             (aTHETA_MMP3_TNF * TNF_PROT)) / (ONE +
     6             (aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     7             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     8             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_COMPL) +
     9             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_COMPL) +
     1             (aTHETA_MMP3_TNF * TNF_PROT)))) * (ONE -
     2             (((iTHETA_MMP3_MAG + iTHETA_MMP3_FREQ) /
     3             (aTHETA_MMP3_GLUC + aTHETA_MMP3_LACT +
     4             aTHETA_MMP3_TNF +
     5             iTHETA_MMP3_MAG + iTHETA_MMP3_FREQ)) *
     6             ((ONE + IC_int_MMP3_dir) / (IC_int_MMP3_dir)) *
     7             (((iTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt) +
     8             (iTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt)) / (ONE +
     9             (iTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt) +
     1             (iTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt)))))
     2             , 4)

C MMP3 Expression for IL1B & TNF Immunopositive Cells
         MMP3_IL1B_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     2             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     3             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_COMPL) +
     4             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_COMPL) +
     5             (aTHETA_MMP3_IL1B * IL1B_PROT) +
     6             (aTHETA_MMP3_TNF * TNF_PROT)) / (ONE +
     7             (aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     8             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     9             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_COMPL) +
     1             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_COMPL) +
     2             (aTHETA_MMP3_IL1B * IL1B_PROT) +
     3             (aTHETA_MMP3_TNF * TNF_PROT)))) * (ONE -
     4             (((iTHETA_MMP3_MAG + iTHETA_MMP3_FREQ) /
     5             (aTHETA_MMP3_GLUC + aTHETA_MMP3_LACT +
     6             aTHETA_MMP3_IL1B + aTHETA_MMP3_TNF +
     7             iTHETA_MMP3_MAG + iTHETA_MMP3_FREQ)) *
     8             ((ONE + IC_int_MMP3_dir) / (IC_int_MMP3_dir)) *
     9             (((iTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt) +
     1             (iTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt)) / (ONE +
     2             (iTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt) +
     3             (iTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt)))))
     4             , 4)

C--------------------------------------------------------------------------
C Magnitude <= 0.998 and Frequency > 3.0 (Magnitude anabolic, Frequency catabolic)
C--------------------------------------------------------------------------

      ELSE IF (MAG .LE. 0.998D0 .AND. FREQ .GT. 3.0D0) THEN

C Magnitude Anabolic
         MMP3_MAG_ANA_RANGE = I_A_MAG_ANA * MMP3_MAG_FACT
         IF (MMP3_MAG_ANA_RANGE .GE. ABS(X_S_MMP3_MAG)) THEN
            MMP3_MAG_ANA_RANGE = ABS(X_S_MMP3_MAG)
            X_S_MMP3_MAG_COMPL = X_S_MMP3_MAG_COMPL +
     1                          (I_A_MAG_CAT_1 * MMP3_MAG_FACT)
            IF (X_S_MMP3_MAG_COMPL .GT. ONE) THEN
               X_S_MMP3_MAG_COMPL = ONE
            END IF
         ELSE
            X_S_MMP3_MAG_COMPL = ZERO
         END IF
         X_S_MMP3_MAG_DEPt = ABS(X_S_MMP3_MAG + MMP3_MAG_ANA_RANGE)

C Frequency Catabolic
         MMP3_FREQ_CAT_RANGE = I_A_FREQ_CAT_2 * MMP3_FREQ_FACT_2
         IF (MMP3_FREQ_CAT_RANGE .GT. ONE) THEN
            MMP3_FREQ_CAT_RANGE = ONE
         END IF
         X_S_MMP3_FREQ_DEPt = ZERO + MMP3_FREQ_CAT_RANGE

C MMP3 Expression for Immunonegative Cells
         MMP3_IMN = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     2             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     3             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_COMPL) +
     4             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt)) / (ONE +
     5             (aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     6             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     7             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_COMPL) +
     8             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt)))) * (ONE -
     9             (((iTHETA_MMP3_MAG) /
     1             (aTHETA_MMP3_GLUC + aTHETA_MMP3_LACT +
     2             iTHETA_MMP3_MAG + aTHETA_MMP3_FREQ)) *
     3             ((ONE + IC_int_MMP3_dir) / (IC_int_MMP3_dir)) *
     4             (((iTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt)) / (ONE +
     5             (iTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt)))))
     6             , 4)

C MMP3 Expression for IL1B Immunopositive Cells
         MMP3_IL1B = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     2             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     3             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_COMPL) +
     4             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt) +
     5             (aTHETA_MMP3_IL1B * IL1B_PROT)) / (ONE +
     6             (aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     7             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     8             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_COMPL) +
     9             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt) +
     1             (aTHETA_MMP3_IL1B * IL1B_PROT)))) * (ONE -
     2             (((iTHETA_MMP3_MAG) /
     3             (aTHETA_MMP3_GLUC + aTHETA_MMP3_LACT +
     4             aTHETA_MMP3_IL1B +
     5             iTHETA_MMP3_MAG + aTHETA_MMP3_FREQ)) *
     6             ((ONE + IC_int_MMP3_dir) / (IC_int_MMP3_dir)) *
     7             (((iTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt)) / (ONE +
     8             (iTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt)))))
     9             , 4)

C MMP3 Expression for TNF Immunopositive Cells
         MMP3_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     2             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     3             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_COMPL) +
     4             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt) +
     5             (aTHETA_MMP3_TNF * TNF_PROT)) / (ONE +
     6             (aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     7             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     8             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_COMPL) +
     9             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt) +
     1             (aTHETA_MMP3_TNF * TNF_PROT)))) * (ONE -
     2             (((iTHETA_MMP3_MAG) /
     3             (aTHETA_MMP3_GLUC + aTHETA_MMP3_LACT +
     4             aTHETA_MMP3_TNF +
     5             iTHETA_MMP3_MAG + aTHETA_MMP3_FREQ)) *
     6             ((ONE + IC_int_MMP3_dir) / (IC_int_MMP3_dir)) *
     7             (((iTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt)) / (ONE +
     8             (iTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt)))))
     9             , 4)

C MMP3 Expression for IL1B & TNF Immunopositive Cells
         MMP3_IL1B_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     2             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     3             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_COMPL) +
     4             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt) +
     5             (aTHETA_MMP3_IL1B * IL1B_PROT) +
     6             (aTHETA_MMP3_TNF * TNF_PROT)) / (ONE +
     7             (aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     8             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     9             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_COMPL) +
     1             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt) +
     2             (aTHETA_MMP3_IL1B * IL1B_PROT) +
     3             (aTHETA_MMP3_TNF * TNF_PROT)))) * (ONE -
     4             (((iTHETA_MMP3_MAG) /
     5             (aTHETA_MMP3_GLUC + aTHETA_MMP3_LACT +
     6             aTHETA_MMP3_IL1B + aTHETA_MMP3_TNF +
     7             iTHETA_MMP3_MAG + aTHETA_MMP3_FREQ)) *
     8             ((ONE + IC_int_MMP3_dir) / (IC_int_MMP3_dir)) *
     9             (((iTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt)) / (ONE +
     1             (iTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt)))))
     2             , 4)

C--------------------------------------------------------------------------
C Magnitude >= 1.0 and Frequency <= 3.0 (Magnitude catabolic, Frequency anabolic)
C--------------------------------------------------------------------------

      ELSE IF (MAG .GE. ONE .AND. FREQ .LE. 3.0D0 .AND.
     1         FREQ .NE. ZERO) THEN

C Magnitude Catabolic
            MMP3_MAG_CAT_RANGE = I_A_MAG_CAT_2 * MMP3_MAG_FACT_2
            IF (MMP3_MAG_CAT_RANGE .GT. ONE) THEN
                  MMP3_MAG_CAT_RANGE = ONE
            END IF
            X_S_MMP3_MAG_DEPt = ZERO + MMP3_MAG_CAT_RANGE

C Frequency Anabolic
            MMP3_FREQ_ANA_RANGE = I_A_FREQ_ANA * MMP3_FREQ_FACT
            IF (MMP3_FREQ_ANA_RANGE .GE. ABS(X_S_MMP3_FREQ)) THEN
                  MMP3_FREQ_ANA_RANGE = ABS(X_S_MMP3_FREQ)
                  X_S_MMP3_FREQ_COMPL = X_S_MMP3_FREQ_COMPL +
     1                                 (I_A_FREQ_CAT_1 * MMP3_FREQ_FACT)
                  IF (X_S_MMP3_FREQ_COMPL .GT. ONE) THEN
                        X_S_MMP3_FREQ_COMPL = ONE
                  END IF
            ELSE
                  X_S_MMP3_FREQ_COMPL = ZERO
            END IF
            X_S_MMP3_FREQ_DEPt = ABS(X_S_MMP3_FREQ +
     1                           MMP3_FREQ_ANA_RANGE)

C MMP3 Expression for Immunonegative Cells
            MMP3_IMN = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     2             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     3             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt) +
     4             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_COMPL)) / (ONE +
     5             (aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     6             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     7             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt) +
     8             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_COMPL)))) * (ONE -
     9             (((iTHETA_MMP3_FREQ) /
     1             (aTHETA_MMP3_GLUC + aTHETA_MMP3_LACT +
     2             aTHETA_MMP3_MAG + iTHETA_MMP3_FREQ)) *
     3             ((ONE + IC_int_MMP3_dir) / (IC_int_MMP3_dir)) *
     4             (((iTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt)) / (ONE +
     5             (iTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt)))))
     6             , 4)

C MMP3 Expression for IL1B Immunopositive Cells
            MMP3_IL1B = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     2             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     3             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt) +
     4             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_COMPL) +
     5             (aTHETA_MMP3_IL1B * IL1B_PROT)) / (ONE +
     6             (aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     7             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     8             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt) +
     9             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_COMPL) +
     1             (aTHETA_MMP3_IL1B * IL1B_PROT)))) * (ONE -
     2             (((iTHETA_MMP3_FREQ) /
     3             (aTHETA_MMP3_GLUC + aTHETA_MMP3_LACT +
     4             aTHETA_MMP3_MAG + aTHETA_MMP3_IL1B +
     5             iTHETA_MMP3_FREQ)) *
     6             ((ONE + IC_int_MMP3_dir) / (IC_int_MMP3_dir)) *
     7             (((iTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt)) / (ONE +
     8             (iTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt)))))
     9             , 4)

C MMP3 Expression for TNF Immunopositive Cells
            MMP3_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     2             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     3             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt) +
     4             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_COMPL) +
     5             (aTHETA_MMP3_TNF * TNF_PROT)) / (ONE +
     6             (aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     7             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     8             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt) +
     9             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_COMPL) +
     1             (aTHETA_MMP3_TNF * TNF_PROT)))) * (ONE -
     2             (((iTHETA_MMP3_FREQ) /
     3             (aTHETA_MMP3_GLUC + aTHETA_MMP3_LACT +
     4             aTHETA_MMP3_MAG + aTHETA_MMP3_TNF +
     5             iTHETA_MMP3_FREQ)) *
     6             ((ONE + IC_int_MMP3_dir) / (IC_int_MMP3_dir)) *
     7             (((iTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt)) / (ONE +
     8             (iTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt)))))
     9             , 4)

C MMP3 Expression for IL1B & TNF Immunopositive Cells
            MMP3_IL1B_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     2             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     3             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt) +
     4             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_COMPL) +
     5             (aTHETA_MMP3_IL1B * IL1B_PROT) +
     6             (aTHETA_MMP3_TNF * TNF_PROT)) / (ONE +
     7             (aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     8             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     9             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt) +
     1             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_COMPL) +
     2             (aTHETA_MMP3_IL1B * IL1B_PROT) +
     3             (aTHETA_MMP3_TNF * TNF_PROT)))) * (ONE -
     4             (((iTHETA_MMP3_FREQ) /
     5             (aTHETA_MMP3_GLUC + aTHETA_MMP3_LACT +
     6             aTHETA_MMP3_MAG + aTHETA_MMP3_IL1B +
     7             aTHETA_MMP3_TNF + iTHETA_MMP3_FREQ)) *
     8             ((ONE + IC_int_MMP3_dir) / (IC_int_MMP3_dir)) *
     9             (((iTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt)) / (ONE +
     1             (iTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt)))))
     2             , 4)

C--------------------------------------------------------------------------
C Magnitude >= 1.0 and Frequency == 0 (Magnitude catabolic, Frequency anabolic)
C--------------------------------------------------------------------------

      ELSE IF (MAG .GE. ONE .AND. FREQ .EQ. ZERO) THEN

C Magnitude Catabolic
            MMP3_MAG_CAT_RANGE = I_A_MAG_CAT_2 * MMP3_MAG_FACT_2
            IF (MMP3_MAG_CAT_RANGE .GT. ONE) THEN
                  MMP3_MAG_CAT_RANGE = ONE
            END IF
            X_S_MMP3_MAG_DEPt = ZERO + MMP3_MAG_CAT_RANGE

C Frequency Anabolic (Frequency == 0 Case)
            MMP3_FREQ_ANA_RANGE = FREQ_0_TIME *
     1                            DBLE(TICKS) * MMP3_FREQ_FACT
            IF (MMP3_FREQ_ANA_RANGE .GE. FREQ_0) THEN
                  MMP3_FREQ_ANA_RANGE = FREQ_0
                  X_S_MMP3_FREQ_COMPL = X_S_MMP3_FREQ_COMPL +
     1                                  (FREQ_0_TIME * MMP3_FREQ_FACT)
                  IF (X_S_MMP3_FREQ_COMPL .GT. ONE) THEN
                        X_S_MMP3_FREQ_COMPL = ONE
                  END IF
            ELSE
                  X_S_MMP3_FREQ_COMPL = ZERO
            END IF
            X_S_MMP3_FREQ_DEPt = FREQ_0 - MMP3_FREQ_ANA_RANGE

C MMP3 Expression for Immunonegative Cells
            MMP3_IMN = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     2             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     3             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt) +
     4             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_COMPL)) / (ONE +
     5             (aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     6             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     7             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt) +
     8             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_COMPL)))) * (ONE -
     9             (((iTHETA_MMP3_FREQ) /
     1             (aTHETA_MMP3_GLUC + aTHETA_MMP3_LACT +
     2             aTHETA_MMP3_MAG + iTHETA_MMP3_FREQ)) *
     3             ((ONE + IC_int_MMP3_dir) / (IC_int_MMP3_dir)) *
     4             (((iTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt)) / (ONE +
     5             (iTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt)))))
     6             , 4)

C MMP3 Expression for IL1B Immunopositive Cells
            MMP3_IL1B = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     2             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     3             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt) +
     4             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_COMPL) +
     5             (aTHETA_MMP3_IL1B * IL1B_PROT)) / (ONE +
     6             (aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     7             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     8             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt) +
     9             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_COMPL) +
     1             (aTHETA_MMP3_IL1B * IL1B_PROT)))) * (ONE -
     2             (((iTHETA_MMP3_FREQ) /
     3             (aTHETA_MMP3_GLUC + aTHETA_MMP3_LACT +
     4             aTHETA_MMP3_MAG + aTHETA_MMP3_IL1B +
     5             iTHETA_MMP3_FREQ)) *
     6             ((ONE + IC_int_MMP3_dir) / (IC_int_MMP3_dir)) *
     7             (((iTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt)) / (ONE +
     8             (iTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt)))))
     9             , 4)

C MMP3 Expression for TNF Immunopositive Cells
            MMP3_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     2             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     3             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt) +
     4             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_COMPL) +
     5             (aTHETA_MMP3_TNF * TNF_PROT)) / (ONE +
     6             (aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     7             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     8             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt) +
     9             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_COMPL) +
     1             (aTHETA_MMP3_TNF * TNF_PROT)))) * (ONE -
     2             (((iTHETA_MMP3_FREQ) /
     3             (aTHETA_MMP3_GLUC + aTHETA_MMP3_LACT +
     4             aTHETA_MMP3_MAG + aTHETA_MMP3_TNF +
     5             iTHETA_MMP3_FREQ)) *
     6             ((ONE + IC_int_MMP3_dir) / (IC_int_MMP3_dir)) *
     7             (((iTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt)) / (ONE +
     8             (iTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt)))))
     9             , 4)

C MMP3 Expression for IL1B & TNF Immunopositive Cells
            MMP3_IL1B_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM))*
     1             (((aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     2             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     3             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt) +
     4             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_COMPL) +
     5             (aTHETA_MMP3_IL1B * IL1B_PROT) +
     6             (aTHETA_MMP3_TNF * TNF_PROT)) / (ONE +
     7             (aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     8             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     9             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt) +
     1             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_COMPL) +
     2             (aTHETA_MMP3_IL1B * IL1B_PROT) +
     3             (aTHETA_MMP3_TNF * TNF_PROT)))) * (ONE -
     4             (((iTHETA_MMP3_FREQ) /
     5             (aTHETA_MMP3_GLUC + aTHETA_MMP3_LACT +
     6             aTHETA_MMP3_MAG + aTHETA_MMP3_IL1B +
     7             aTHETA_MMP3_TNF + iTHETA_MMP3_FREQ)) *
     8             ((ONE + IC_int_MMP3_dir) / (IC_int_MMP3_dir)) *
     9             (((iTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt)) / (ONE +
     1             (iTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt)))))
     2             , 4)

C--------------------------------------------------------------------------
C Magnitude >= 1.0 and Frequency > 3.0 (Magnitude catabolic, Frequency catabolic)
C--------------------------------------------------------------------------
      ELSE IF (MAG .GE. ONE .AND. FREQ .GT. 3.0D0) THEN

C Magnitude Catabolic
            MMP3_MAG_CAT_RANGE = I_A_MAG_CAT_2 * MMP3_MAG_FACT_2
            IF (MMP3_MAG_CAT_RANGE .GT. ONE) THEN
                  MMP3_MAG_CAT_RANGE = ONE
            END IF
            X_S_MMP3_MAG_DEPt = MMP3_MAG_CAT_RANGE

C Frequency Catabolic
            MMP3_FREQ_CAT_RANGE = I_A_FREQ_CAT_2 * MMP3_FREQ_FACT_2
            IF (MMP3_FREQ_CAT_RANGE .GT. ONE) THEN
                  MMP3_FREQ_CAT_RANGE = ONE
            END IF
            X_S_MMP3_FREQ_DEPt = MMP3_FREQ_CAT_RANGE

C MMP3 Expression for Immunonegative Cells
            MMP3_IMN = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     2             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     3             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt) +
     4             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt)) / (ONE +
     5             (aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     6             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     7             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt) +
     8             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt))))
     9             , 4)

C MMP3 Expression for IL1B Immunopositive Cells
            MMP3_IL1B = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     2             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     3             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt) +
     4             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt) +
     5             (aTHETA_MMP3_IL1B * IL1B_PROT)) / (ONE +
     6             (aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     7             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     8             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt) +
     9             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt) +
     1             (aTHETA_MMP3_IL1B * IL1B_PROT))))
     2             , 4)

C MMP3 Expression for TNF Immunopositive Cells
            MMP3_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     2             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     3             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt) +
     4             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt) +
     5             (aTHETA_MMP3_TNF * TNF_PROT)) / (ONE +
     6             (aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     7             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     8             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt) +
     9             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt) +
     1             (aTHETA_MMP3_TNF * TNF_PROT))))
     2             , 4)

C MMP3 Expression for IL1B & TNF Immunopositive Cells
            MMP3_IL1B_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     2             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     3             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt) +
     4             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt) +
     5             (aTHETA_MMP3_IL1B * IL1B_PROT) +
     6             (aTHETA_MMP3_TNF * TNF_PROT)) / (ONE +
     7             (aTHETA_MMP3_LACT * X_S_MMP3_LACT) +
     8             (aTHETA_MMP3_GLUC * X_S_MMP3_GLUC) +
     9             (aTHETA_MMP3_MAG * X_S_MMP3_MAG_DEPt) +
     1             (aTHETA_MMP3_FREQ * X_S_MMP3_FREQ_DEPt) +
     2             (aTHETA_MMP3_IL1B * IL1B_PROT) +
     3             (aTHETA_MMP3_TNF * TNF_PROT))))
     4             , 4)
      
      END IF
      
C---------------------------------------------------------------------------
      RETURN
      END SUBROUTINE PNeq_MMP3

C--------------------------------------------------------------------------------------------------------------------------------------------- 
C Evaluate the ADAMTS-4 expression
C ---------------------------------------------------------

      SUBROUTINE PNeq_ADM4(ADM4_IMN, ADM4_IL1B,
     1           ADM4_TNF, ADM4_IL1B_TNF,
     2           TICKS,
     3           IL1B_PROT, TNF_PROT,
     4           GLUC_int, pH_val, MAG, FREQ,
     5           FREQ_0, FREQ_0_time,
     6           MAG_0, MAG_0_time,
     7           aTHETA_ADM4_GLUC, aTHETA_ADM4_LACT,
     8           iTHETA_ADM4_IL1B, aTHETA_ADM4_TNF,
     9           aTHETA_ADM4_MAG, iTHETA_ADM4_MAG,
     1           aTHETA_ADM4_FREQ, iTHETA_ADM4_FREQ,
     2           AC_int_IL1B_dir, AC_int_TNF_dir,
     3           IC_int_IL1B_dir, IC_int_TNF_dir,
     4           IC_int_ADM4_dir,
     5           X_S_ADM4_GLUC, X_S_ADM4_LACT,
     6           X_S_ADM4_MAG, X_S_ADM4_FREQ,
     7           A_I_MAG_ANA, A_I_MAG_CAT_1, A_I_MAG_CAT_2,
     8           I_A_MAG_ANA, I_A_MAG_CAT_1, I_A_MAG_CAT_2,
     9           A_I_FREQ_ANA, A_I_FREQ_CAT_1, A_I_FREQ_CAT_2,
     1           I_A_FREQ_ANA, I_A_FREQ_CAT_1, I_A_FREQ_CAT_2,
     2           ADM4_MAG_FACT, ADM4_MAG_FACT_2,
     3           ADM4_FREQ_FACT, ADM4_FREQ_FACT_2,
     4           X_S_ADM4_MAG_COMPL, X_S_ADM4_FREQ_COMPL,
     5           aTHETA_SUM)

C VARIABLES:
C X_S: SENSITIVITY OF A CELL ACTIVITY (CA) TO A STIMULUS CONCENTRATION
C THETA: SENSITIVITY OF A CELL ACTIVITY (CA) TO A STIMULUS TYPE
C OUTPUTS:
C ADM4_IMN: ADM4RECAN EXPRESSION FOR IMMUNONEGATIVE CELLS
C ADM4_IL1B: CONTRIBUTION FROM IL1B
C ADM4_TNF: CONTRIBUTION FROM TNF
C ADM4_IL1B_TNF: COMBINED CONTRIBUTION OF IL1B AND TNF

C ********************************************************************
C
      IMPLICIT NONE
C
C FUNCTION
      DOUBLE PRECISION ROUND

C INTEGER VARIABLES
      INTEGER TICKS

C SPECIAL CASE FOR FREQ = 0
      DOUBLE PRECISION FREQ_0, FREQ_0_time
      DOUBLE PRECISION MAG_0, MAG_0_time

C Input and Output Variables
      DOUBLE PRECISION ADM4_IMN, ADM4_IL1B, ADM4_TNF, ADM4_IL1B_TNF
      DOUBLE PRECISION IL1B_PROT, TNF_PROT
      DOUBLE PRECISION GLUC_int, pH_val, MAG, FREQ

C Constants for Sensitivity Factors
      DOUBLE PRECISION aTHETA_ADM4_GLUC, aTHETA_ADM4_LACT
      DOUBLE PRECISION iTHETA_ADM4_IL1B, aTHETA_ADM4_TNF
      DOUBLE PRECISION aTHETA_ADM4_MAG, iTHETA_ADM4_MAG
      DOUBLE PRECISION aTHETA_ADM4_FREQ, iTHETA_ADM4_FREQ
      DOUBLE PRECISION AC_int_IL1B_dir, AC_int_TNF_dir
      DOUBLE PRECISION IC_int_IL1B_dir, IC_int_TNF_dir
      DOUBLE PRECISION IC_int_ADM4_dir
      DOUBLE PRECISION aTHETA_SUM

C Intermediate Variables
      DOUBLE PRECISION X_S_ADM4_GLUC, X_S_ADM4_LACT
      DOUBLE PRECISION X_S_ADM4_MAG, X_S_ADM4_FREQ
      DOUBLE PRECISION ADM4_MAG_FACT, ADM4_MAG_FACT_2
      DOUBLE PRECISION ADM4_FREQ_FACT, ADM4_FREQ_FACT_2

C Completive factors
      DOUBLE PRECISION X_S_ADM4_MAG_COMPL
      DOUBLE PRECISION X_S_ADM4_MAG_DEPt
      DOUBLE PRECISION X_S_ADM4_FREQ_COMPL
      DOUBLE PRECISION X_S_ADM4_FREQ_DEPt

C Time dependent funtions
      DOUBLE PRECISION A_I_MAG_ANA, A_I_MAG_CAT_1, A_I_MAG_CAT_2
      DOUBLE PRECISION I_A_MAG_ANA, I_A_MAG_CAT_1, I_A_MAG_CAT_2
      DOUBLE PRECISION A_I_FREQ_ANA, A_I_FREQ_CAT_1, A_I_FREQ_CAT_2
      DOUBLE PRECISION I_A_FREQ_ANA, I_A_FREQ_CAT_1, I_A_FREQ_CAT_2

C Time sensitivity range
      DOUBLE PRECISION ADM4_MAG_ANA_RANGE
      DOUBLE PRECISION ADM4_FREQ_ANA_RANGE
      DOUBLE PRECISION ADM4_MAG_CAT_RANGE
      DOUBLE PRECISION ADM4_FREQ_CAT_RANGE

C CONSTANT PARAMETERS
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.D0, ONE = 1.D0

C---------------------------------------------------------------------------
C Initialize Variables
      ADM4_MAG_ANA_RANGE = ZERO
      ADM4_FREQ_ANA_RANGE = ZERO
      ADM4_MAG_CAT_RANGE = ZERO
      ADM4_FREQ_CAT_RANGE = ZERO
      X_S_ADM4_MAG_DEPt = ZERO
      X_S_ADM4_FREQ_DEPt = ZERO

C---------------------------------------------------------------------------
C Magnitude <= 0.998 and Frequency <= 3.0 (Magnitude & Frequency anabolic)
C---------------------------------------------------------------------------

      IF (MAG .LE. 0.998D0 .AND. FREQ .LE. 3.0D0 .AND.
     1    FREQ .NE. ZERO) THEN

C Magnitude Anabolic
            ADM4_MAG_ANA_RANGE = I_A_MAG_ANA * ADM4_MAG_FACT
            IF (ADM4_MAG_ANA_RANGE .GE. ABS(X_S_ADM4_MAG)) THEN
                  ADM4_MAG_ANA_RANGE = ABS(X_S_ADM4_MAG)
                  X_S_ADM4_MAG_COMPL = X_S_ADM4_MAG_COMPL +
     1                                 (I_A_MAG_CAT_1 * ADM4_MAG_FACT)
                  IF (X_S_ADM4_MAG_COMPL .GT. ONE) THEN
                        X_S_ADM4_MAG_COMPL = ONE
                  END IF
            ELSE
            X_S_ADM4_MAG_COMPL = ZERO
            END IF
            X_S_ADM4_MAG_DEPt = ABS(X_S_ADM4_MAG +
     1                             ADM4_MAG_ANA_RANGE)

C Frequency Anabolic
            ADM4_FREQ_ANA_RANGE = I_A_FREQ_ANA * ADM4_FREQ_FACT
            IF (ADM4_FREQ_ANA_RANGE .GE. ABS(X_S_ADM4_FREQ)) THEN
            ADM4_FREQ_ANA_RANGE = ABS(X_S_ADM4_FREQ)
            X_S_ADM4_FREQ_COMPL = X_S_ADM4_FREQ_COMPL +
     1                            (I_A_FREQ_CAT_1 * ADM4_FREQ_FACT)
                  IF (X_S_ADM4_FREQ_COMPL .GT. ONE) THEN
                        X_S_ADM4_FREQ_COMPL = ONE
                  END IF
            ELSE
            X_S_ADM4_FREQ_COMPL = ZERO
            END IF
            X_S_ADM4_FREQ_DEPt = ABS(X_S_ADM4_FREQ +
     1                              ADM4_FREQ_ANA_RANGE)

C ADM4 mRNA expression for immunonegative cells
            ADM4_IMN = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     2             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     3             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_COMPL) +
     4             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_COMPL)) /
     5             (ONE + (aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     6             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     7             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_COMPL) +
     8             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_COMPL)))) * (ONE -
     9             (((iTHETA_ADM4_MAG + iTHETA_ADM4_FREQ) /
     1             (aTHETA_ADM4_GLUC + aTHETA_ADM4_LACT +
     2             iTHETA_ADM4_MAG + iTHETA_ADM4_FREQ)) *
     3             ((ONE + IC_int_ADM4_dir) / (IC_int_ADM4_dir)) *
     4             (((iTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     5             (iTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt)) / (ONE +
     6             (iTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     7             (iTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt)))))
     8             , 4)

C ADM4 mRNA expression for IL1B immunopositive cells
            ADM4_IL1B = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     2             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     3             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_COMPL) +
     4             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_COMPL)) / (ONE +
     5             (aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     6             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     7             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_COMPL) +
     8             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_COMPL)))) * (ONE -
     9             (((iTHETA_ADM4_MAG + iTHETA_ADM4_FREQ +
     1             iTHETA_ADM4_IL1B) /
     2             (aTHETA_ADM4_GLUC + aTHETA_ADM4_LACT +
     3             iTHETA_ADM4_MAG + iTHETA_ADM4_FREQ +
     4             iTHETA_ADM4_IL1B)) *
     5             ((ONE + IC_int_ADM4_dir) / (IC_int_ADM4_dir)) *
     6             (((iTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     7             (iTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt) +
     8             (iTHETA_ADM4_IL1B * IL1B_PROT)) / (ONE +
     9             (iTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     1             (iTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt) +
     2             (iTHETA_ADM4_IL1B * IL1B_PROT)))))
     3             , 4)

C ADM4 mRNA expression for TNF immunopositive cells
            ADM4_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     2             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     3             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_COMPL) +
     4             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_COMPL) +
     5             (aTHETA_ADM4_TNF * TNF_PROT)) / (ONE +
     6             (aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     7             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     8             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_COMPL) +
     9             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_COMPL) +
     1             (aTHETA_ADM4_TNF * TNF_PROT)))) * (ONE -
     2             (((iTHETA_ADM4_MAG + iTHETA_ADM4_FREQ) /
     3             (aTHETA_ADM4_GLUC + aTHETA_ADM4_LACT +
     4             aTHETA_ADM4_TNF +
     5             iTHETA_ADM4_MAG + iTHETA_ADM4_FREQ)) *
     6             ((ONE + IC_int_ADM4_dir) / (IC_int_ADM4_dir)) *
     7             (((iTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     8             (iTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt)) / (ONE +
     9             (iTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     1             (iTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt)))))
     2             , 4)

C ADM4 mRNA expression for IL1B and TNF immunopositive cells
            ADM4_IL1B_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     2             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     3             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_COMPL) +
     4             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_COMPL) +
     5             (aTHETA_ADM4_TNF * TNF_PROT)) / (ONE +
     6             (aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     7             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     8             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_COMPL) +
     9             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_COMPL) +
     1             (aTHETA_ADM4_TNF * TNF_PROT)))) * (ONE -
     2             (((iTHETA_ADM4_MAG + iTHETA_ADM4_FREQ +
     3             iTHETA_ADM4_IL1B) /
     4             (aTHETA_ADM4_GLUC + aTHETA_ADM4_LACT +
     5             aTHETA_ADM4_TNF +
     6             iTHETA_ADM4_MAG + iTHETA_ADM4_FREQ +
     7             iTHETA_ADM4_IL1B)) *
     8             ((ONE + IC_int_ADM4_dir) / (IC_int_ADM4_dir)) *
     9             (((iTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     1             (iTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt) +
     2             (iTHETA_ADM4_IL1B * IL1B_PROT)) / (ONE +
     3             (iTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     4             (iTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt) +
     5             (iTHETA_ADM4_IL1B * IL1B_PROT)))))
     6             , 4)

C---------------------------------------------------------------------------
C Magnitude <= 0.998 and Frequency = 0 (Magnitude anabolic, Frequency catabolic)
C---------------------------------------------------------------------------

      ELSE IF (MAG .LE. 0.998D0 .AND. FREQ .EQ. ZERO) THEN

C Magnitude Anabolic
            ADM4_MAG_ANA_RANGE = I_A_MAG_ANA * ADM4_MAG_FACT
            IF (ADM4_MAG_ANA_RANGE .GE. ABS(X_S_ADM4_MAG)) THEN
                  ADM4_MAG_ANA_RANGE = ABS(X_S_ADM4_MAG)
                  X_S_ADM4_MAG_COMPL = X_S_ADM4_MAG_COMPL +
     1                                 (I_A_MAG_CAT_1 * ADM4_MAG_FACT)
                  IF (X_S_ADM4_MAG_COMPL .GE. ONE) THEN
                        X_S_ADM4_MAG_COMPL = ONE
                  END IF
            ELSE
                  X_S_ADM4_MAG_COMPL = ZERO
            END IF
            X_S_ADM4_MAG_DEPt = ABS(X_S_ADM4_MAG + ADM4_MAG_ANA_RANGE)

C Frequency Anabolic (Frequency == 0 Case)
            ADM4_FREQ_ANA_RANGE = FREQ_0_TIME *
     1                            DBLE(TICKS) * ADM4_FREQ_FACT
            IF (ADM4_FREQ_ANA_RANGE .GE. FREQ_0) THEN
                  ADM4_FREQ_ANA_RANGE = FREQ_0
                  X_S_ADM4_FREQ_COMPL = X_S_ADM4_FREQ_COMPL +
     1                                  (FREQ_0_TIME * ADM4_FREQ_FACT)
                  IF (X_S_ADM4_FREQ_COMPL .GE. ONE) THEN
                        X_S_ADM4_FREQ_COMPL = ONE
                  END IF
            ELSE
                  X_S_ADM4_FREQ_COMPL = ZERO
            END IF
            X_S_ADM4_FREQ_DEPt = FREQ_0 - ADM4_FREQ_ANA_RANGE

C ADM4 mRNA expression for immunonegative cells
            ADM4_IMN = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     2             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     3             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_COMPL) +
     4             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_COMPL)) / (ONE +
     5             (aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     6             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     7             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_COMPL) +
     8             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_COMPL)))) * (ONE -
     9             (((iTHETA_ADM4_MAG + iTHETA_ADM4_FREQ) /
     1             (aTHETA_ADM4_GLUC + aTHETA_ADM4_LACT +
     2             iTHETA_ADM4_MAG + iTHETA_ADM4_FREQ)) *
     3             ((ONE + IC_int_ADM4_dir) / (IC_int_ADM4_dir)) *
     4             (((iTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     5             (iTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt)) / (ONE +
     6             (iTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     7             (iTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt)))))
     8             , 4)

C ADM4 mRNA expression for IL1B immunopositive cells
            ADM4_IL1B = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     2             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     3             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_COMPL) +
     4             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_COMPL)) / (ONE +
     5             (aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     6             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     7             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_COMPL) +
     8             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_COMPL)))) * (ONE -
     9             (((iTHETA_ADM4_MAG + iTHETA_ADM4_FREQ +
     1             iTHETA_ADM4_IL1B) /
     2             (aTHETA_ADM4_GLUC + aTHETA_ADM4_LACT +
     3             iTHETA_ADM4_MAG + iTHETA_ADM4_FREQ +
     4             iTHETA_ADM4_IL1B)) *
     5             ((ONE + IC_int_ADM4_dir) / (IC_int_ADM4_dir)) *
     6             (((iTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     7             (iTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt) +
     8             (iTHETA_ADM4_IL1B * IL1B_PROT)) / (ONE +
     9             (iTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     1             (iTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt) +
     2             (iTHETA_ADM4_IL1B * IL1B_PROT)))))
     3             , 4)

C ADM4 mRNA expression for TNF immunopositive cells
            ADM4_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     2             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     3             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_COMPL) +
     4             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_COMPL) +
     5             (aTHETA_ADM4_TNF * TNF_PROT)) / (ONE +
     6             (aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     7             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     8             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_COMPL) +
     9             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_COMPL) +
     1             (aTHETA_ADM4_TNF * TNF_PROT)))) * (ONE -
     2             (((iTHETA_ADM4_MAG + iTHETA_ADM4_FREQ) /
     3             (aTHETA_ADM4_GLUC + aTHETA_ADM4_LACT +
     4             aTHETA_ADM4_TNF +
     5             iTHETA_ADM4_MAG + iTHETA_ADM4_FREQ)) *
     6             ((ONE + IC_int_ADM4_dir) / (IC_int_ADM4_dir)) *
     7             (((iTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     8             (iTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt)) / (ONE +
     9             (iTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     1             (iTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt)))))
     2             , 4)

C ADM4 mRNA expression for IL1B and TNF immunopositive cells
            ADM4_IL1B_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     2             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     3             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_COMPL) +
     4             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_COMPL) +
     5             (aTHETA_ADM4_TNF * TNF_PROT)) / (ONE +
     6             (aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     7             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     8             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_COMPL) +
     9             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_COMPL) +
     1             (aTHETA_ADM4_TNF * TNF_PROT)))) * (ONE -
     2             (((iTHETA_ADM4_MAG + iTHETA_ADM4_FREQ +
     3             iTHETA_ADM4_IL1B) /
     4             (aTHETA_ADM4_GLUC + aTHETA_ADM4_LACT +
     5             aTHETA_ADM4_TNF +
     6             iTHETA_ADM4_MAG + iTHETA_ADM4_FREQ +
     7             iTHETA_ADM4_IL1B)) *
     8             ((ONE + IC_int_ADM4_dir) / (IC_int_ADM4_dir)) *
     9             (((iTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     1             (iTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt) +
     2             (iTHETA_ADM4_IL1B * IL1B_PROT)) / (ONE +
     3             (iTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     4             (iTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt) +
     5             (iTHETA_ADM4_IL1B * IL1B_PROT)))))
     6             , 4)

C---------------------------------------------------------------------------
C Magnitude <= 0.998 and Frequency > 3.0 (Magnitude anabolic, Frequency catabolic)
C---------------------------------------------------------------------------

      ELSE IF (MAG .LE. 0.998D0 .AND. FREQ .GT. 3.0D0) THEN

C Magnitude Anabolic
            ADM4_MAG_ANA_RANGE = I_A_MAG_ANA * ADM4_MAG_FACT
            IF (ADM4_MAG_ANA_RANGE .GE. ABS(X_S_ADM4_MAG)) THEN
                  ADM4_MAG_ANA_RANGE = ABS(X_S_ADM4_MAG)
                  X_S_ADM4_MAG_COMPL = X_S_ADM4_MAG_COMPL +
     1                                 (I_A_MAG_CAT_1 * ADM4_MAG_FACT)
                  IF (X_S_ADM4_MAG_COMPL .GE. ONE) THEN
                        X_S_ADM4_MAG_COMPL = ONE
                  END IF
            ELSE
                  X_S_ADM4_MAG_COMPL = ZERO
            END IF
            X_S_ADM4_MAG_DEPt = ABS(X_S_ADM4_MAG + ADM4_MAG_ANA_RANGE)

C Frequency Catabolic
            ADM4_FREQ_CAT_RANGE = I_A_FREQ_CAT_2 * ADM4_FREQ_FACT_2
            IF (ADM4_FREQ_CAT_RANGE .GT. ONE) THEN
                  ADM4_FREQ_CAT_RANGE = ONE
            END IF
            X_S_ADM4_FREQ_DEPt = ZERO + ADM4_FREQ_CAT_RANGE

C ADM4 mRNA expression for immunonegative cells
            ADM4_IMN = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     2             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     3             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_COMPL) +
     4             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt)) / (ONE +
     5             (aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     6             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     7             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_COMPL) +
     8             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt)))) * (ONE -
     9             (((iTHETA_ADM4_MAG) /
     1             (aTHETA_ADM4_GLUC + aTHETA_ADM4_LACT +
     2             iTHETA_ADM4_MAG + aTHETA_ADM4_FREQ)) *
     3             ((ONE + IC_int_ADM4_dir) / (IC_int_ADM4_dir)) *
     4             (((iTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt)) / (ONE +
     5             (iTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt)))))
     6             , 4)

C ADM4 mRNA expression for IL1B immunopositive cells
            ADM4_IL1B = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     2             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     3             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_COMPL) +
     4             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt)) / (ONE +
     5             (aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     6             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     7             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_COMPL) +
     8             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt)))) * (ONE -
     9             (((iTHETA_ADM4_MAG + iTHETA_ADM4_IL1B) /
     1             (aTHETA_ADM4_GLUC + aTHETA_ADM4_LACT +
     2             iTHETA_ADM4_MAG + aTHETA_ADM4_FREQ +
     3             iTHETA_ADM4_IL1B)) *
     4             ((ONE + IC_int_ADM4_dir) / (IC_int_ADM4_dir)) *
     5             (((iTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     6             (iTHETA_ADM4_IL1B * IL1B_PROT)) / (ONE +
     7             (iTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     8             (iTHETA_ADM4_IL1B * IL1B_PROT)))))
     9             , 4)

C ADM4 mRNA expression for TNF immunopositive cells
            ADM4_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     2             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     3             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_COMPL) +
     4             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt) +
     5             (aTHETA_ADM4_TNF * TNF_PROT)) / (ONE +
     6             (aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     7             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     8             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_COMPL) +
     9             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt) +
     1             (aTHETA_ADM4_TNF * TNF_PROT)))) * (ONE -
     2             (((iTHETA_ADM4_MAG) /
     3             (aTHETA_ADM4_GLUC + aTHETA_ADM4_LACT +
     4             aTHETA_ADM4_TNF +
     5             iTHETA_ADM4_MAG + aTHETA_ADM4_FREQ)) *
     6             ((ONE + IC_int_ADM4_dir) / (IC_int_ADM4_dir)) *
     7             (((iTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt)) / (ONE +
     8             (iTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt)))))
     9             , 4)

C ADM4 mRNA expression for IL1B and TNF immunopositive cells
            ADM4_IL1B_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     2             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     3             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_COMPL) +
     4             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt) +
     5             (aTHETA_ADM4_TNF * TNF_PROT)) / (ONE +
     6             (aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     7             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     8             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_COMPL) +
     9             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt) +
     1             (aTHETA_ADM4_TNF * TNF_PROT)))) * (ONE -
     2             (((iTHETA_ADM4_MAG + iTHETA_ADM4_IL1B) /
     3             (aTHETA_ADM4_GLUC + aTHETA_ADM4_LACT +
     4             aTHETA_ADM4_TNF +
     5             iTHETA_ADM4_MAG + aTHETA_ADM4_FREQ +
     6             iTHETA_ADM4_IL1B)) *
     7             ((ONE + IC_int_ADM4_dir) / (IC_int_ADM4_dir)) *
     8             (((iTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     9             (iTHETA_ADM4_IL1B * IL1B_PROT)) / (ONE +
     1             (iTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     2             (iTHETA_ADM4_IL1B * IL1B_PROT)))))
     3             , 4)

C---------------------------------------------------------------------------
C Magnitude >= 1.0 and Frequency <= 3.0 (Magnitude catabolic, Frequency anabolic)
C---------------------------------------------------------------------------

      ELSE IF (MAG .GE. ONE .AND. FREQ .LE. 3.0D0 .AND.
     1         FREQ .NE. ZERO) THEN

C Magnitude Catabolic
            ADM4_MAG_CAT_RANGE = I_A_MAG_CAT_2 * ADM4_MAG_FACT_2
            IF (ADM4_MAG_CAT_RANGE .GT. ONE) THEN
                  ADM4_MAG_CAT_RANGE = ONE
            END IF
            X_S_ADM4_MAG_DEPt = ZERO + ADM4_MAG_CAT_RANGE

C Frequency Anabolic
            ADM4_FREQ_ANA_RANGE = I_A_FREQ_ANA * ADM4_FREQ_FACT
            IF (ADM4_FREQ_ANA_RANGE .GE. ABS(X_S_ADM4_FREQ)) THEN
                  ADM4_FREQ_ANA_RANGE = ABS(X_S_ADM4_FREQ)
                  X_S_ADM4_FREQ_COMPL = X_S_ADM4_FREQ_COMPL +
     1                                 (I_A_FREQ_CAT_1 * ADM4_FREQ_FACT)
                  IF (X_S_ADM4_FREQ_COMPL .GT. ONE) THEN
                        X_S_ADM4_FREQ_COMPL = ONE
                  END IF
            ELSE
                  X_S_ADM4_FREQ_COMPL = ZERO
            END IF
            X_S_ADM4_FREQ_DEPt = ABS(X_S_ADM4_FREQ +
     1                           ADM4_FREQ_ANA_RANGE)

C ADM4 mRNA expression for immunonegative cells
            ADM4_IMN = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     2             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     3             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     4             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_COMPL)) / (ONE +
     5             (aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     6             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     7             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     8             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_COMPL)))) * (ONE -
     9             (((iTHETA_ADM4_FREQ) /
     1             (aTHETA_ADM4_GLUC + aTHETA_ADM4_LACT +
     2             aTHETA_ADM4_MAG + iTHETA_ADM4_FREQ)) *
     3             ((ONE + IC_int_ADM4_dir) / (IC_int_ADM4_dir)) *
     4             (((iTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt)) / (ONE +
     5             (iTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt)))))
     6             , 4)

C ADM4 mRNA expression for IL1B immunopositive cells
            ADM4_IL1B = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     2             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     3             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     4             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_COMPL)) / (ONE +
     5             (aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     6             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     7             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     8             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_COMPL)))) * (ONE -
     9             (((iTHETA_ADM4_FREQ + iTHETA_ADM4_IL1B) /
     1             (aTHETA_ADM4_GLUC + aTHETA_ADM4_LACT +
     2             aTHETA_ADM4_MAG + iTHETA_ADM4_FREQ +
     3             iTHETA_ADM4_IL1B)) *
     4             ((ONE + IC_int_ADM4_dir) / (IC_int_ADM4_dir)) *
     5             (((iTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt) +
     6             (iTHETA_ADM4_IL1B * IL1B_PROT)) / (ONE +
     7             (iTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt) +
     8             (iTHETA_ADM4_IL1B * IL1B_PROT)))))
     9             , 4)

C ADM4 mRNA expression for TNF immunopositive cells
            ADM4_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     2             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     3             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     4             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_COMPL) +
     5             (aTHETA_ADM4_TNF * TNF_PROT)) / (ONE +
     6             (aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     7             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     8             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     9             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_COMPL) +
     1             (aTHETA_ADM4_TNF * TNF_PROT)))) * (ONE -
     2             (((iTHETA_ADM4_FREQ) /
     3             (aTHETA_ADM4_GLUC + aTHETA_ADM4_LACT +
     4             aTHETA_ADM4_MAG + aTHETA_ADM4_TNF +
     5             iTHETA_ADM4_FREQ)) *
     6             ((ONE + IC_int_ADM4_dir) / (IC_int_ADM4_dir)) *
     7             (((iTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt)) / (ONE +
     8             (iTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt)))))
     9             , 4)

C ADM4 mRNA expression for IL1B and TNF immunopositive cells
            ADM4_IL1B_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     2             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     3             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     4             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_COMPL) +
     5             (aTHETA_ADM4_TNF * TNF_PROT)) / (ONE +
     6             (aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     7             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     8             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     9             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_COMPL) +
     1             (aTHETA_ADM4_TNF * TNF_PROT)))) * (ONE -
     2             (((iTHETA_ADM4_FREQ + iTHETA_ADM4_IL1B) /
     3             (aTHETA_ADM4_GLUC + aTHETA_ADM4_LACT +
     4             aTHETA_ADM4_MAG + aTHETA_ADM4_TNF +
     5             iTHETA_ADM4_FREQ + iTHETA_ADM4_IL1B)) *
     6             ((ONE + IC_int_ADM4_dir) / (IC_int_ADM4_dir)) *
     7             (((iTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt) +
     8             (iTHETA_ADM4_IL1B * IL1B_PROT)) / (ONE +
     9             (iTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt) +
     1             (iTHETA_ADM4_IL1B * IL1B_PROT)))))
     2             , 4)

C---------------------------------------------------------------------------
C Magnitude >= 1.0 and Frequency = 0 (Magnitude catabolic, Frequency 0)
C---------------------------------------------------------------------------

      ELSE IF (MAG .GE. ONE .AND. FREQ .EQ. ZERO) THEN

C Magnitude Catabolic
            ADM4_MAG_CAT_RANGE = I_A_MAG_CAT_2 * ADM4_MAG_FACT_2
            IF (ADM4_MAG_CAT_RANGE .GE. ONE) THEN
                  ADM4_MAG_CAT_RANGE = ONE
            END IF
            X_S_ADM4_MAG_DEPt = ZERO + ADM4_MAG_CAT_RANGE

C Frequency Anabolic (Frequency == 0 Case)
            ADM4_FREQ_ANA_RANGE = FREQ_0_TIME *
     1                            DBLE(TICKS) * ADM4_FREQ_FACT
            IF (ADM4_FREQ_ANA_RANGE .GE. FREQ_0) THEN
                  ADM4_FREQ_ANA_RANGE = FREQ_0
                  X_S_ADM4_FREQ_COMPL = X_S_ADM4_FREQ_COMPL +
     1                                  (FREQ_0_TIME * ADM4_FREQ_FACT)
                  IF (X_S_ADM4_FREQ_COMPL .GE. ONE) THEN
                        X_S_ADM4_FREQ_COMPL = ONE
                  END IF
            ELSE
                  X_S_ADM4_FREQ_COMPL = ZERO
            END IF
            X_S_ADM4_FREQ_DEPt = FREQ_0 - ADM4_FREQ_ANA_RANGE

C ADM4 mRNA expression for immunonegative cells
            ADM4_IMN = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     2             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     3             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     4             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_COMPL)) / (ONE +
     5             (aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     6             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     7             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     8             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_COMPL)))) * (ONE -
     9             (((iTHETA_ADM4_FREQ) /
     1             (aTHETA_ADM4_GLUC + aTHETA_ADM4_LACT +
     2             aTHETA_ADM4_MAG + iTHETA_ADM4_FREQ)) *
     3             ((ONE + IC_int_ADM4_dir) / (IC_int_ADM4_dir)) *
     4             (((iTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt)) / (ONE +
     5             (iTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt)))))
     6             , 4)

C ADM4 mRNA expression for IL1B immunopositive cells
            ADM4_IL1B = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     2             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     3             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     4             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_COMPL)) / (ONE +
     5             (aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     6             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     7             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     8             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_COMPL)))) * (ONE -
     9             (((iTHETA_ADM4_FREQ + iTHETA_ADM4_IL1B) /
     1             (aTHETA_ADM4_GLUC + aTHETA_ADM4_LACT +
     2             aTHETA_ADM4_MAG + iTHETA_ADM4_FREQ +
     3             iTHETA_ADM4_IL1B)) *
     4             ((ONE + IC_int_ADM4_dir) / (IC_int_ADM4_dir)) *
     5             (((iTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt) +
     6             (iTHETA_ADM4_IL1B * IL1B_PROT)) / (ONE +
     7             (iTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt) +
     8             (iTHETA_ADM4_IL1B * IL1B_PROT)))))
     9             , 4)

C ADM4 mRNA expression for TNF immunopositive cells
            ADM4_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     2             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     3             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     4             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_COMPL) +
     5             (aTHETA_ADM4_TNF * TNF_PROT)) / (ONE +
     6             (aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     7             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     8             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     9             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_COMPL) +
     1             (aTHETA_ADM4_TNF * TNF_PROT)))) * (ONE -
     2             (((iTHETA_ADM4_FREQ) /
     3             (aTHETA_ADM4_GLUC + aTHETA_ADM4_LACT +
     4             aTHETA_ADM4_MAG + aTHETA_ADM4_TNF +
     5             iTHETA_ADM4_FREQ)) *
     6             ((ONE + IC_int_ADM4_dir) / (IC_int_ADM4_dir)) *
     7             (((iTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt)) / (ONE +
     8             (iTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt)))))
     9             , 4)

C ADM4 mRNA expression for IL1B and TNF immunopositive cells
            ADM4_IL1B_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     2             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     3             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     4             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_COMPL) +
     5             (aTHETA_ADM4_TNF * TNF_PROT)) / (ONE +
     6             (aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     7             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     8             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     9             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_COMPL) +
     1             (aTHETA_ADM4_TNF * TNF_PROT)))) * (ONE -
     2             (((iTHETA_ADM4_FREQ + iTHETA_ADM4_IL1B) /
     3             (aTHETA_ADM4_GLUC + aTHETA_ADM4_LACT +
     4             aTHETA_ADM4_MAG + aTHETA_ADM4_TNF +
     5             iTHETA_ADM4_FREQ + iTHETA_ADM4_IL1B)) *
     6             ((ONE + IC_int_ADM4_dir) / (IC_int_ADM4_dir)) *
     7             (((iTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt) +
     8             (iTHETA_ADM4_IL1B * IL1B_PROT)) / (ONE +
     9             (iTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt) +
     1             (iTHETA_ADM4_IL1B * IL1B_PROT)))))
     2             , 4)

C---------------------------------------------------------------------------
C Magnitude >= 1.0 and Frequency > 3.0 (Magnitude & Frequency catabolic)
C---------------------------------------------------------------------------

      ELSE IF (MAG .GE. ONE .AND. FREQ .GT. 3.0D0) THEN

C Magnitude Catabolic
            ADM4_MAG_CAT_RANGE = I_A_MAG_CAT_2 * ADM4_MAG_FACT_2
            IF (ADM4_MAG_CAT_RANGE .GE. ONE) THEN
                  ADM4_MAG_CAT_RANGE = ONE
            END IF
            X_S_ADM4_MAG_DEPt = ZERO + ADM4_MAG_CAT_RANGE

C Frequency Catabolic
            ADM4_FREQ_CAT_RANGE = I_A_FREQ_CAT_2 * ADM4_FREQ_FACT_2
            IF (ADM4_FREQ_CAT_RANGE .GE. ONE) THEN
                  ADM4_FREQ_CAT_RANGE = ONE
            END IF
            X_S_ADM4_FREQ_DEPt = ZERO + ADM4_FREQ_CAT_RANGE

C ADM4 mRNA expression for immunonegative cells
            ADM4_IMN = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     2             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     3             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     4             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt)) / (ONE +
     5             (aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     6             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     7             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     8             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt)))) , 4)

C ADM4 mRNA expression for IL1B immunopositive cells
            ADM4_IL1B = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     2             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     3             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     4             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt)) / (ONE +
     5             (aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     6             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     7             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     8             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt)))) * (ONE -
     9             ((iTHETA_ADM4_IL1B / (aTHETA_ADM4_LACT +
     1             aTHETA_ADM4_GLUC + aTHETA_ADM4_MAG +
     2             aTHETA_ADM4_FREQ + iTHETA_ADM4_IL1B)) *
     3             ((ONE + IC_int_ADM4_dir) / (IC_int_ADM4_dir)) *
     4             ((iTHETA_ADM4_IL1B * IL1B_PROT) / (ONE +
     5             (iTHETA_ADM4_IL1B * IL1B_PROT))))) , 4)

C ADM4 mRNA expression for TNF immunopositive cells
            ADM4_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     2             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     3             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     4             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt) +
     5             (aTHETA_ADM4_TNF * TNF_PROT)) / (ONE +
     6             (aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     7             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     8             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     9             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt) +
     1             (aTHETA_ADM4_TNF * TNF_PROT)))) , 4)

C ADM4 mRNA expression for IL1B and TNF immunopositive cells
            ADM4_IL1B_TNF = ROUND((((ONE + aTHETA_SUM) / (aTHETA_SUM)) *
     1             (((aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     2             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     3             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     4             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt) +
     5             (aTHETA_ADM4_TNF * TNF_PROT)) / (ONE +
     6             (aTHETA_ADM4_LACT * X_S_ADM4_LACT) +
     7             (aTHETA_ADM4_GLUC * X_S_ADM4_GLUC) +
     8             (aTHETA_ADM4_MAG * X_S_ADM4_MAG_DEPt) +
     9             (aTHETA_ADM4_FREQ * X_S_ADM4_FREQ_DEPt) +
     1             (aTHETA_ADM4_TNF * TNF_PROT)))) * (ONE -
     2             ((iTHETA_ADM4_IL1B / (aTHETA_ADM4_LACT +
     3             aTHETA_ADM4_GLUC + aTHETA_ADM4_MAG +
     4             aTHETA_ADM4_FREQ + aTHETA_ADM4_TNF +
     5             iTHETA_ADM4_IL1B)) *
     6             ((ONE + IC_int_ADM4_dir) / (IC_int_ADM4_dir)) *
     7             ((iTHETA_ADM4_IL1B * IL1B_PROT) / (ONE +
     8             (iTHETA_ADM4_IL1B * IL1B_PROT))))) , 4)

      END IF

C---------------------------------------------------------------------------
      RETURN
      END SUBROUTINE PNeq_ADM4

C---------------------------------------------------------------------------------------------------------------------------------------------
C Evaluate Sensitivity (X_S) of a Cell activity (CA) to a stimulus (S) concentration
C ----------------------------------------------------------------------------------
C This function calculates the value of the curve using the formula:
C X_S_curve = c0 / (1.0 + exp(-c1 * (S - c2))) - c3
C Inputs:
C   c0    - Amplitude of the curve (scalar)
C   c1    - Steepness of the curve (scalar)
C   c2    - Midpoint of the curve (scalar)
C   c3    - Offset of the curve (scalar)
C   S     - Stimulus concentration (scalar)
C Outputs:
C   X_S_curve - Evaluated result of the curve (scalar)
C---------------------------------------------------------------------------------------------------------------------------------------------
      FUNCTION X_S_curve(C0, C1, C2, C3, C4, S)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION X_S_curve
      DOUBLE PRECISION C0, C1, C2, C3, C4, S
      DOUBLE PRECISION EXPONENT

C Initialize the variables
      X_S_curve = 0.0D0
      EXPONENT = 0.0D0

C Evaluate the curve
      EXPONENT = C1 * S
      X_S_curve = C0 * EXP(EXPONENT) /
     1            (EXP(EXPONENT) + C2 * 10.0D0**C3) + C4

C--------------------------------------------------------------------------      
      RETURN
      END FUNCTION X_S_curve

C---------------------------------------------------------------------------
C Function to round a value to a specified number of decimal places
C---------------------------------------------------------------------------
C This function rounds a value to a specified number of decimal places using:
C ROUND(VAL, N) = INT(VAL * 10^N + 0.5) / 10^N
C Inputs:
C   VAL - Value to be rounded (scalar, DOUBLE PRECISION)
C   N   - Number of decimal places to round to (scalar, INTEGER)
C Output:
C   The rounded value (DOUBLE PRECISION)
C---------------------------------------------------------------------------
      FUNCTION ROUND(VAL, N)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION :: VAL, ROUND, SCALE
      INTEGER :: N

      ! Scale the value to the desired decimal precision
      SCALE = 10.0D0 ** N

      ! Perform rounding using INT for Fortran 77 compatibility
      IF (VAL .GE. 0.0D0) THEN
            ROUND = DBLE(INT(VAL * SCALE + 0.5D0)) / SCALE
      ELSE
            ROUND = DBLE(INT(VAL * SCALE - 0.5D0)) / SCALE
      END IF

C--------------------------------------------------------------------------      
      RETURN
      END FUNCTION ROUND
