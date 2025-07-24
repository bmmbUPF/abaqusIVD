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
      SUBROUTINE PNeq(
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
     +           PN_TIME, PNeq_STEP, PN_FREQ_file)

C VARIABLES:
C I,J,K: Loop integers
C DTIME: Time increment
C KSTEP: Current step
C NDIM: Number of dimensions
C stat: Status of the subroutine

C NECESSARY VARIABLES OF THE PROBLEM:
C GLUC_int: GLUCOSE CONCENTRATIONS 0 - 5.0 mM
C CELL_rho: ACTUAL CELL DENSITY
C pH_val: ACTUAL pH 6.5 - 7.4
C MAG: PRESSURE MAGNITUDE 0.1 - 3.5 MPa
C FREQ: PRESSURE FREQUENCY 0 - 40 Hz
C WEIGHTHING FACTORS (THETA VALUES)

C INPUTS:
C PN_TIME: TOTAL TIME OF THE PN EQ. SIMULATION


C OUTPUTS:
C FOUR AGGRECANS STATES:
C AGG_IMN, AGG_IL1B, AGG_TNF, AGG_IL1B_TNF
C FOUR COLLAGEN I STATES:
C COLI_IMN, COLI_IL1B, COLI_TNF, COLI_IL1B_TNF
C FOUR COLLAGEN II STATES:
C COLII_IMN, COLII_IL1B, COLII_TNF, COLII_IL1B_TNF
C FOUR MMP-3 STATES:
C MMP3_IMN, MMP3_IL1B, MMP3_TNF, MMP3_IL1B_TNF
C FOUR ADAMTS-4 STATES:
C ADM4_IMN, ADM4_IL1B, ADM4_TNF, ADM4_IL1B_TNF
C TNF AND IL−1β EXPRESSION LEVELS:
C TNF, IL1B

C
C ********************************************************************
C
      IMPLICIT NONE
C
C VARIABLES:
C LOOP INTEGERS
      INTEGER I, TICKS

C NECESSARY VARIABLES OF THE PROBLEM:
C DEFORMATION GRADIENT
c GLUCOSE CONCENTRATION, CELL DENSITY, pH, PRESSURE MAGNITUDE AND FREQUENCY
      DOUBLE PRECISION GLUC_int, pH_val, MAG, FREQ

C INPUTS:
C TOTAL TIME (IN HOURS) OF THE PN EQ. SIMULATION
      INTEGER PN_TIME
C Are we in the PN eq. step?
      INTEGER PNeq_STEP
C If yes, then the frequency for every hour
      DOUBLE PRECISION PN_FREQ_file(*)

C OUTPUTS:
C FOUR AGGRECANS STATES
      DOUBLE PRECISION AGG_IMN, AGG_IL1B
      DOUBLE PRECISION AGG_TNF, AGG_IL1B_TNF
C ACCUMULATED PROTEIN LEVELS
      DOUBLE PRECISION AGG_IMN_ACC, AGG_IL1B_ACC
      DOUBLE PRECISION AGG_TNF_ACC, AGG_IL1B_TNF_ACC

C FOUR COLLAGEN I STATES
      DOUBLE PRECISION COLI_IMN, COLI_IL1B
      DOUBLE PRECISION COLI_TNF, COLI_IL1B_TNF
C ACCUMULATED PROTEIN LEVELS
      DOUBLE PRECISION COLI_IMN_ACC, COLI_IL1B_ACC
      DOUBLE PRECISION COLI_TNF_ACC, COLI_IL1B_TNF_ACC

C FOUR COLLAGEN II STATES
      DOUBLE PRECISION COLII_IMN, COLII_IL1B
      DOUBLE PRECISION COLII_TNF, COLII_IL1B_TNF
C ACCUMULATED PROTEIN LEVELS
      DOUBLE PRECISION COLII_IMN_ACC, COLII_IL1B_ACC
      DOUBLE PRECISION COLII_TNF_ACC, COLII_IL1B_TNF_ACC

C FOUR MMP-3 STATES
      DOUBLE PRECISION MMP3_IMN, MMP3_IL1B
      DOUBLE PRECISION MMP3_TNF, MMP3_IL1B_TNF
C ACCUMULATED PROTEIN LEVELS
      DOUBLE PRECISION MMP3_IMN_ACC, MMP3_IL1B_ACC
      DOUBLE PRECISION MMP3_TNF_ACC, MMP3_IL1B_TNF_ACC

C FOUR ADAMTS-4 STATES
      DOUBLE PRECISION ADM4_IMN, ADM4_IL1B
      DOUBLE PRECISION ADM4_TNF, ADM4_IL1B_TNF
C ACCUMULATED PROTEIN LEVELS
      DOUBLE PRECISION ADM4_IMN_ACC, ADM4_IL1B_ACC
      DOUBLE PRECISION ADM4_TNF_ACC, ADM4_IL1B_TNF_ACC

C TNF AND IL−1β EXPRESSION LEVELS
      DOUBLE PRECISION IL1B, TNF
      DOUBLE PRECISION IL1B_previous, TNF_previous
      DOUBLE PRECISION IL1B_PROT, TNF_PROT

C WEIGHTHING FACTORS (THETA VALUES)
C IL−1β
      DOUBLE PRECISION aTHETA_IL1B_GLUC, aTHETA_IL1B_LACT
      DOUBLE PRECISION aTHETA_IL1B_MAG, aTHETA_IL1B_FREQ
      DOUBLE PRECISION iTHETA_IL1B_MAG, iTHETA_IL1B_FREQ
C TNF
      DOUBLE PRECISION aTHETA_TNF_GLUC, aTHETA_TNF_LACT
      DOUBLE PRECISION aTHETA_TNF_MAG, aTHETA_TNF_FREQ
      DOUBLE PRECISION iTHETA_TNF_MAG, iTHETA_TNF_FREQ
C AGG
      DOUBLE PRECISION aTHETA_AGG_GLUC, aTHETA_AGG_LACT
      DOUBLE PRECISION iTHETA_AGG_IL1B, iTHETA_AGG_TNF
      DOUBLE PRECISION aTHETA_AGG_MAG, aTHETA_AGG_FREQ
      DOUBLE PRECISION iTHETA_AGG_MAG, iTHETA_AGG_FREQ
C COLI
      DOUBLE PRECISION aTHETA_COLI_GLUC, aTHETA_COLI_LACT
      DOUBLE PRECISION aTHETA_COLI_IL1B, iTHETA_COLI_TNF
      DOUBLE PRECISION aTHETA_COLI_MAG, aTHETA_COLI_FREQ
      DOUBLE PRECISION iTHETA_COLI_MAG, iTHETA_COLI_FREQ
C COLII
      DOUBLE PRECISION aTHETA_COLII_GLUC, aTHETA_COLII_LACT
      DOUBLE PRECISION iTHETA_COLII_IL1B, iTHETA_COLII_TNF
      DOUBLE PRECISION aTHETA_COLII_MAG, aTHETA_COLII_FREQ
      DOUBLE PRECISION iTHETA_COLII_MAG, iTHETA_COLII_FREQ
C MMP-3
      DOUBLE PRECISION aTHETA_MMP3_GLUC, aTHETA_MMP3_LACT
      DOUBLE PRECISION aTHETA_MMP3_IL1B, aTHETA_MMP3_TNF
      DOUBLE PRECISION aTHETA_MMP3_MAG, aTHETA_MMP3_FREQ
      DOUBLE PRECISION iTHETA_MMP3_MAG, iTHETA_MMP3_FREQ
C ADAMTS-4
      DOUBLE PRECISION aTHETA_ADM4_GLUC, aTHETA_ADM4_LACT
      DOUBLE PRECISION iTHETA_ADM4_IL1B, aTHETA_ADM4_TNF
      DOUBLE PRECISION aTHETA_ADM4_MAG, aTHETA_ADM4_FREQ
      DOUBLE PRECISION iTHETA_ADM4_MAG, iTHETA_ADM4_FREQ
C PN-denominator
      DOUBLE PRECISION aTHETA_SUM

C INTERNAL COMPARISON VARIABLES
C IL−1β
      DOUBLE PRECISION AC_int_IL1B_dir
      DOUBLE PRECISION IC_int_IL1B_dir
C TNF
      DOUBLE PRECISION AC_int_TNF_dir
      DOUBLE PRECISION IC_int_TNF_dir
C AGG
      DOUBLE PRECISION IC_int_AGG_dir
C COLI      
      DOUBLE PRECISION IC_int_COLI_dir
C COLII
      DOUBLE PRECISION IC_int_COLII_dir
C MMP-3
      DOUBLE PRECISION IC_int_MMP3_dir
C ADAMTS-4
      DOUBLE PRECISION IC_int_ADM4_dir

C SETUP FUNCTION (X_S VALUES)
C IL−1β
      DOUBLE PRECISION X_S_IL1B_GLUC, X_S_IL1B_LACT
      DOUBLE PRECISION X_S_IL1B_MAG, X_S_IL1B_FREQ
C TNF
      DOUBLE PRECISION X_S_TNF_GLUC, X_S_TNF_LACT
      DOUBLE PRECISION X_S_TNF_MAG, X_S_TNF_FREQ
C AGG
      DOUBLE PRECISION X_S_AGG_GLUC, X_S_AGG_LACT
      DOUBLE PRECISION X_S_AGG_MAG, X_S_AGG_FREQ
C COLI
      DOUBLE PRECISION X_S_COLI_GLUC, X_S_COLI_LACT
      DOUBLE PRECISION X_S_COLI_MAG, X_S_COLI_FREQ
C COLII
      DOUBLE PRECISION X_S_COLII_GLUC, X_S_COLII_LACT
      DOUBLE PRECISION X_S_COLII_MAG, X_S_COLII_FREQ
C MMP-3
      DOUBLE PRECISION X_S_MMP3_GLUC, X_S_MMP3_LACT
      DOUBLE PRECISION X_S_MMP3_MAG, X_S_MMP3_FREQ
C ADAMTS-4
      DOUBLE PRECISION X_S_ADM4_GLUC, X_S_ADM4_LACT
      DOUBLE PRECISION X_S_ADM4_MAG, X_S_ADM4_FREQ

C TIME SENSITIVITY FACTORS
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

C SPECIAL CASE FOR FREQ = 0
      DOUBLE PRECISION FREQ_0, FREQ_0_time
      DOUBLE PRECISION MAG_0, MAG_0_time

C GENERIC FUNCTIONS TIME DEPENDENT
      DOUBLE PRECISION A_I_MAG_gFun, I_A_MAG_gFun
      DOUBLE PRECISION A_I_FREQ_gFun, I_A_FREQ_gFun
      DOUBLE PRECISION A_I_MAG_ANA
      DOUBLE PRECISION A_I_MAG_CAT_1
      DOUBLE PRECISION A_I_MAG_CAT_2
      DOUBLE PRECISION I_A_MAG_ANA
      DOUBLE PRECISION I_A_MAG_CAT_1
      DOUBLE PRECISION I_A_MAG_CAT_2
      DOUBLE PRECISION A_I_FREQ_ANA
      DOUBLE PRECISION A_I_FREQ_CAT_1
      DOUBLE PRECISION A_I_FREQ_CAT_2
      DOUBLE PRECISION I_A_FREQ_ANA
      DOUBLE PRECISION I_A_FREQ_CAT_1
      DOUBLE PRECISION I_A_FREQ_CAT_2

C TEMPORAL VARIABLES
      DOUBLE PRECISION X_S_IL1B_MAG_COMPL, X_S_IL1B_FREQ_COMPL
      DOUBLE PRECISION X_S_TNF_MAG_COMPL, X_S_TNF_FREQ_COMPL
      DOUBLE PRECISION X_S_AGG_MAG_COMPL, X_S_AGG_FREQ_COMPL
      DOUBLE PRECISION X_S_COLI_MAG_COMPL, X_S_COLI_FREQ_COMPL
      DOUBLE PRECISION X_S_COLII_MAG_COMPL, X_S_COLII_FREQ_COMPL
      DOUBLE PRECISION X_S_MMP3_MAG_COMPL, X_S_MMP3_FREQ_COMPL
      DOUBLE PRECISION X_S_ADM4_MAG_COMPL, X_S_ADM4_FREQ_COMPL

C CONSTANT PARAMETERS:
C SIMPLE DOUBLE PRECISION NUMBERS
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.D0, ONE=1.0D0

C--------------------------------------------------------------------------
C Initializing variables
C LOOP INTEGERS
      I = 0
C PROTEINS
      IL1B_previous = ZERO
      TNF_previous = ZERO

C AGGRECAN
      AGG_IMN = ZERO
      AGG_IL1B = ZERO
      AGG_TNF = ZERO
      AGG_IL1B_TNF = ZERO

C COLLAGEN I
      COLI_IMN = ZERO
      COLI_IL1B = ZERO
      COLI_TNF = ZERO
      COLI_IL1B_TNF = ZERO

C COLLAGEN II
      COLII_IMN = ZERO
      COLII_IL1B = ZERO
      COLII_TNF = ZERO
      COLII_IL1B_TNF = ZERO

C MMP-3
      MMP3_IMN = ZERO
      MMP3_IL1B = ZERO
      MMP3_TNF = ZERO
      MMP3_IL1B_TNF = ZERO

C ADAMTS-4
      ADM4_IMN = ZERO
      ADM4_IL1B = ZERO
      ADM4_TNF = ZERO
      ADM4_IL1B_TNF = ZERO

C PN variables
      aTHETA_SUM = ZERO
      AC_int_IL1B_dir = ZERO
      IC_int_IL1B_dir = ZERO
      AC_int_TNF_dir = ZERO
      IC_int_TNF_dir = ZERO
      IC_int_AGG_dir = ZERO
      IC_int_COLI_dir = ZERO
      IC_int_COLII_dir = ZERO
      IC_int_MMP3_dir = ZERO
      IC_int_ADM4_dir = ZERO

C Temporal variables
      X_S_IL1B_MAG_COMPL = ZERO
      X_S_IL1B_FREQ_COMPL = ZERO
      X_S_TNF_MAG_COMPL = ZERO
      X_S_TNF_FREQ_COMPL = ZERO
      X_S_AGG_MAG_COMPL = ZERO
      X_S_AGG_FREQ_COMPL = ZERO
      X_S_COLI_MAG_COMPL = ZERO
      X_S_COLI_FREQ_COMPL = ZERO
      X_S_COLII_MAG_COMPL = ZERO
      X_S_COLII_FREQ_COMPL = ZERO
      X_S_MMP3_MAG_COMPL = ZERO
      X_S_MMP3_FREQ_COMPL = ZERO
      X_S_ADM4_MAG_COMPL = ZERO
      X_S_ADM4_FREQ_COMPL = ZERO

C--------------------------------------------------------------------------
C
C--------------------------------------------------------------------
C Weighting factors (theta values):
C--------------------------------------------------------------------
C
C  Theta-a: weighting factors for all activating S-CA relaitonships 
C  Theta-b: all of the inhibiting weighting factors within the same CA(independent of the CS)

C PN-activity is defined by values to the 4th decimal place in order to
C reflect the sensitivity in changes within the CA 

      CALL Compute_Theta_Values(
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
C
C--------------------------------------------------------------------
C Internal comparison within different CS
C--------------------------------------------------------------------
C
C Comparison within different CS (c-int; stands for "internal comparison"); to build the second inhibiting term of the PN-Equation (and in terms of IL1B and TNF, also the first activating ter of their individual PN-Equation
C
C -------- for IL1B --------
      AC_int_IL1B_dir = aTHETA_IL1B_GLUC + aTHETA_IL1B_LACT +
     1                  aTHETA_IL1B_MAG + aTHETA_IL1B_FREQ
      IC_int_IL1B_dir = iTHETA_IL1B_MAG + iTHETA_IL1B_FREQ

C -------- for TNF ---------
      AC_int_TNF_dir = aTHETA_TNF_GLUC + aTHETA_TNF_LACT +
     1                 aTHETA_TNF_MAG + aTHETA_TNF_FREQ
      IC_int_TNF_dir = iTHETA_TNF_MAG + iTHETA_TNF_FREQ

C -------- for Agg ---------
      IC_int_AGG_dir = iTHETA_AGG_IL1B + iTHETA_AGG_TNF +
     1                 iTHETA_AGG_MAG + iTHETA_AGG_FREQ

C ------- for Col-I --------
      IC_int_COLI_dir = iTHETA_COLI_TNF +
     1                  iTHETA_COLI_MAG + iTHETA_COLI_FREQ

C ------- for Col-II -------
      IC_int_COLII_dir = iTHETA_COLII_IL1B + iTHETA_COLII_TNF +
     1                   iTHETA_COLII_MAG + iTHETA_COLII_FREQ

C ------- for MMP-3 --------
      IC_int_MMP3_dir = iTHETA_MMP3_MAG + iTHETA_MMP3_FREQ

C ------ for ADAMTS-4 ------
      IC_int_ADM4_dir = iTHETA_ADM4_IL1B + 
     1                     iTHETA_ADM4_MAG + iTHETA_ADM4_FREQ

C
C--------------------------------------------------------------------
C Time sensitiviy and individualization of generic functions
C--------------------------------------------------------------------
C
      CALL timeSenFact(
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

C--------------------------------------------------------------------------
C
C--------------------------------------------------------------------
C Evaluate the protein expressions
C--------------------------------------------------------------------
C

      DO TICKS = 0, PN_TIME - 1

C Check if we are in the PN eq. step
      IF (PNeq_STEP .EQ. 1) THEN
            FREQ = PN_FREQ_file(TICKS + 1)
      END IF

C
C--------------------------------------------------------------------
C Setup functions (x values)
C--------------------------------------------------------------------
C
      CALL Compute_X_S(GLUC_int, pH_val, MAG, FREQ,
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
C--------------------------------------------------------------------
C Special case for FREQ = 0
C--------------------------------------------------------------------
C
C Special case for FREQ = 0 (to account for discontiuity of the generic function between 0 Hz 0.1 Hz, check manusrcipt for further explanation) 

      FREQ_0 = 0.75D0
      FREQ_0_time = 1.0D0 - A_I_FREQ_gFun       ! A_I_FREQ_gFun_generic_fct
      MAG_0 = 0.75D0                            ! Same values as for freq assumed
      MAG_0_time = 1.0D0 - A_I_MAG_gFun         ! A_I_MAG_gFun_generic_fct

C
C--------------------------------------------------------------------------
C ----- Generic Functions Time Dependent -----
C Magnitude
      A_I_MAG_ANA = (ONE - A_I_MAG_gFun) * DBLE(TICKS)
      A_I_MAG_CAT_1 = (ONE - A_I_MAG_gFun)
      A_I_MAG_CAT_2 = ABS(A_I_MAG_gFun) * (DBLE(TICKS) + ONE)

      I_A_MAG_ANA = ABS((-ONE) + ABS(I_A_MAG_gFun)) * DBLE(TICKS)
      I_A_MAG_CAT_1 = ABS((-ONE) + ABS(I_A_MAG_gFun))
      I_A_MAG_CAT_2 = I_A_MAG_gFun * (DBLE(TICKS) + ONE)

C Frequency
      A_I_FREQ_ANA = (ONE - A_I_FREQ_gFun) * DBLE(TICKS)
      A_I_FREQ_CAT_1 = (ONE - A_I_FREQ_gFun)
      A_I_FREQ_CAT_2 = ABS(A_I_FREQ_gFun) * (DBLE(TICKS) + ONE)

      I_A_FREQ_ANA = ABS((-ONE) + ABS(I_A_FREQ_gFun)) * DBLE(TICKS)
      I_A_FREQ_CAT_1 = ABS((-ONE) + ABS(I_A_FREQ_gFun))
      I_A_FREQ_CAT_2 = I_A_FREQ_gFun * (DBLE(TICKS) + ONE)

C--------------------------------------------------------------------------
C Evaluate the TNF and IL−1β expression levels
      CALL PNeq_IL1B(IL1B, TICKS,
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

      CALL PNeq_TNF(TNF, TICKS,
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

C--------------------------------------------------------------------------
C ----- estimation proinflammatory cytokines (at protein level, to determine 2nd and 3rd level CA)
      IF (TICKS .EQ. 0) THEN
            IL1B_previous = IL1B
            TNF_previous = TNF
      END IF

C considering half-life in terms of adaptation of IL1B due to a change in IL1B mRNA expression
      IL1B_PROT = (IL1B + IL1B_previous) / 2.0D0
      TNF_PROT = TNF

C Update the previous values for the next iteration
      IL1B_previous = IL1B
      TNF_previous = TNF

C--------------------------------------------------------------------------
C Evaluate the Aggrecan expression
      CALL PNeq_AGG(AGG_IMN, AGG_IL1B, AGG_TNF, AGG_IL1B_TNF,
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

C--------------------------------------------------------------------------
C Evaluate the ColI expression
      CALL PNeq_COLI(COLI_IMN, COLI_IL1B, COLI_TNF, COLI_IL1B_TNF,
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

C--------------------------------------------------------------------------
C Evaluate the ColII expression
      CALL PNeq_COLII(COLII_IMN, COLII_IL1B,
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

C--------------------------------------------------------------------------
C Evaluate the MMP-3 expression
      CALL PNeq_MMP3(MMP3_IMN, MMP3_IL1B,
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

C--------------------------------------------------------------------------
C Evaluate the ADAMTS-4 expression
      CALL PNeq_ADM4(ADM4_IMN, ADM4_IL1B,
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

C--------------------------------------------------------------------------
C Calculate the accumulated protein levels
      IF (TICKS .EQ. 0) THEN
            AGG_IMN_ACC = AGG_IMN
            AGG_IL1B_ACC = AGG_IL1B
            AGG_TNF_ACC = AGG_TNF
            AGG_IL1B_TNF_ACC = AGG_IL1B_TNF

            COLI_IMN_ACC = COLI_IMN
            COLI_IL1B_ACC = COLI_IL1B
            COLI_TNF_ACC = COLI_TNF
            COLI_IL1B_TNF_ACC = COLI_IL1B_TNF

            COLII_IMN_ACC = COLII_IMN
            COLII_IL1B_ACC = COLII_IL1B
            COLII_TNF_ACC = COLII_TNF
            COLII_IL1B_TNF_ACC = COLII_IL1B_TNF

            MMP3_IMN_ACC = MMP3_IMN
            MMP3_IL1B_ACC = MMP3_IL1B
            MMP3_TNF_ACC = MMP3_TNF
            MMP3_IL1B_TNF_ACC = MMP3_IL1B_TNF

            ADM4_IMN_ACC = ADM4_IMN
            ADM4_IL1B_ACC = ADM4_IL1B
            ADM4_TNF_ACC = ADM4_TNF
            ADM4_IL1B_TNF_ACC = ADM4_IL1B_TNF
      ELSE
            AGG_IMN_ACC = AGG_IMN_ACC + AGG_IMN
            AGG_IL1B_ACC = AGG_IL1B_ACC + AGG_IL1B
            AGG_TNF_ACC = AGG_TNF_ACC + AGG_TNF
            AGG_IL1B_TNF_ACC = AGG_IL1B_TNF_ACC + AGG_IL1B_TNF

            COLI_IMN_ACC = COLI_IMN_ACC + COLI_IMN
            COLI_IL1B_ACC = COLI_IL1B_ACC + COLI_IL1B
            COLI_TNF_ACC = COLI_TNF_ACC + COLI_TNF
            COLI_IL1B_TNF_ACC = COLI_IL1B_TNF_ACC + COLI_IL1B_TNF

            COLII_IMN_ACC = COLII_IMN_ACC + COLII_IMN
            COLII_IL1B_ACC = COLII_IL1B_ACC + COLII_IL1B
            COLII_TNF_ACC = COLII_TNF_ACC + COLII_TNF
            COLII_IL1B_TNF_ACC = COLII_IL1B_TNF_ACC + COLII_IL1B_TNF

            MMP3_IMN_ACC = MMP3_IMN_ACC + MMP3_IMN
            MMP3_IL1B_ACC = MMP3_IL1B_ACC + MMP3_IL1B
            MMP3_TNF_ACC = MMP3_TNF_ACC + MMP3_TNF
            MMP3_IL1B_TNF_ACC = MMP3_IL1B_TNF_ACC + MMP3_IL1B_TNF

            ADM4_IMN_ACC = ADM4_IMN_ACC + ADM4_IMN
            ADM4_IL1B_ACC = ADM4_IL1B_ACC + ADM4_IL1B
            ADM4_TNF_ACC = ADM4_TNF_ACC + ADM4_TNF
            ADM4_IL1B_TNF_ACC = ADM4_IL1B_TNF_ACC + ADM4_IL1B_TNF

      END IF

      END DO

C--------------------------------------------------------------------------

      RETURN
      END SUBROUTINE PNeq
C--------------------------------------------------------------------------