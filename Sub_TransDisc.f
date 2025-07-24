C
C  ___  ___  ________     ___    ___ ___  ___  ________     ___    ___ ___    ___  _______  ________
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
C Introduction
C ------------
C

C all the equation can be viewed in https://latex.codecogs.com/eqneditor/editor.php

C Mechanical simulation of the IVD in Abaqus is requiered, which has been run with C3D20P type elements. Now that the simulation is finished
C the idea is to couple the results to a solute transport model, which considers the diffusion-reaction of oxygen, lactate, and glucose.
C This is the equation written in latex for a better understanding:

C \frac{\partial }{\partial t}
C \begin{pmatrix}
C \text{C}_{\text{O}_2}\\
C \text{C}_{\text{lact}}\\
C \text{C}_{\text{gluc}}
C \end{pmatrix}
C -
C \begin{pmatrix}
C \text{D}_{\text{O}_2} & 0 & 0\\
C 0 & \text{D}_{lact} & 0\\
C 0 & 0 & \text{D}_{gluc}
C \end{pmatrix}
C \nabla^{2}
C \begin{pmatrix}
C \text{C}_{\text{O}_2}\\
C \text{C}_{\text{lact}}\\
C \text{C}_{\text{gluc}}
C \end{pmatrix}
C =
C \begin{pmatrix}
C \text{R}_{\text{O}_2}\\
C \text{R}_{\text{lact}}\\
C \text{R}_{\text{gluc}}
C \end{pmatrix}
C

C where $\text{C}_i$, $\text{D}_i$ and $\text{R}_i$ are the concentrations, tissue diffusion coefficients, and the reactions of oxygen ($i=\text{O}_2$),
C lactate ($i=\text{lact}$) and glucose ($i=\text{gluc}$), respectively.
C ###

C----------------------------------------------------------------------------------------------------------------------------------------------
C MODULES
C -------------------------

C**********************************************************************************************************************************************
C Module to define the global variables of the simulation
C**********************************************************************************************************************************************

C
      MODULE GLOBAL

C  This module is used to transfer SDV's from the UEL
C  to the UVARM so that SDV's can be visualized on a
C  dummy mesh

C  globalSDVs(X,Y,Z)
C  X - element pointer
C  Y - integration point pointer
C  Z - SDV pointer

C  numElem
C  Total number of elements in the real mesh, the dummy
C  mesh needs to have the same number of elements, and 
C  the dummy mesh needs to have the same number of integ
C  points.  You must set that parameter value here.

C  ElemOffset
C  Offset between element numbers on the real mesh and
C  dummy mesh.  That is set in the input file, and 
C  that value must be set here the same.

C Declare variables to store date and time components:
C From https://gcc.gnu.org/onlinedocs/gfortran/DATE_005fAND_005fTIME.html
C date_ini: Initial date and time values
C date_fin: Current date and time values
C value(1): Year
C value(2): Month
C value(3): Day
C value(4): Time difference from UTC in minutes
C value(5): Hour
C value(6): Minute
C value(7): Second
C value(8): Millisecond
      INTEGER :: date_ini(8), date_fin(8)
C date_diff: Difference between the initial and current date and time values
C date_diff(1): Days
C date_diff(2): Hours
C date_diff(3): Minutes
C date_diff(4): Seconds
C date_diff(5): Total seconds
      INTEGER :: date_diff(5)

C error variable
      INTEGER err

C String variable to store necessary information
      CHARACTER*256 userText

C Declare the ITERATION array to access the iteration number
C This is for printin and writing the iteration number
      INTEGER :: ITERATION = 0
      INTEGER :: LastKINC = -1
      INTEGER :: LastStep = -1
      INTEGER :: LastElement = -1

C -----------------------------------------------------------
C Number of dimensions
      INTEGER, PARAMETER :: NDIM = 3

C -----------------------------------------------------------
C Set the number of UEL elements used here
C IVD model: 19392 elements
      INTEGER, PARAMETER :: numElem = 19392
C Diffusion Chamber model: 52944 elements
C      INTEGER, PARAMETER :: numElem = 5148

C -----------------------------------------------------------
C Set the larget id number for the UEL elements
C IVD model: 52944 largest id element
      INTEGER, PARAMETER :: numIdElem = 52944
C Diffusion Chamber model: 35148 largest id element
C      INTEGER, PARAMETER :: numIdElem = 35148

C -----------------------------------------------------------
C Set the offset here for UVARM plotting, must match input file
C IVD model: 30000
      INTEGER, PARAMETER :: ElemOffset = 30000
C Diffusion Chamber model: 30000
C      INTEGER, PARAMETER :: ElemOffset = 30000

C -----------------------------------------------------------
C Set the number of UEL integration points (same as dummy-DECOY)
      INTEGER, PARAMETER :: numInt = 27

C -----------------------------------------------------------
C Set number of elements * number of integration points (previously called sou)
      INTEGER, PARAMETER :: numElemInt = numElem*numInt

C Set the number of SDV's per each integration point considering numIdElem
      INTEGER, PARAMETER :: numSDVxIdElem = numIdElem-ElemOffset

C -----------------------------------------------------------
C SDV's form previous simulation. These should have the dimensions as the integration points
C (Elem, IP = SDV)
      DOUBLE PRECISION, DIMENSION(numSDVxIdElem,numInt) :: SDV_NF0
      DOUBLE PRECISION, DIMENSION(numSDVxIdElem,numInt) :: SDV_SP

C -----------------------------------------------------------
C Set the number of SDV's here
C this variable is used to store the SDV's from the UEL
      DOUBLE PRECISION, ALLOCATABLE :: globalSDVs(:,:,:)
C      DOUBLE PRECISION, ALLOCATABLE :: localSDVs(:,:,:)

C -----------------------------------------------------------
C Set the step names here
      INTEGER NSTEPS
      CHARACTER*256, ALLOCATABLE :: STEPNAMES(:)

C -----------------------------------------------------------
C Defining the SDVs
      INTEGER, PARAMETER :: UVAR_NF0 = 1                                !Initial water content
      INTEGER, PARAMETER :: UVAR_NF = 2                                 !Water content
      INTEGER, PARAMETER :: UVAR_O2 = 3                                 !O2 concentration
      INTEGER, PARAMETER :: UVAR_Lact = 4                               !Lactate concentration
      INTEGER, PARAMETER :: UVAR_Gluc = 5                               !Glucose concentration
      INTEGER, PARAMETER :: UVAR_pH_val = 6                             !pH value
      INTEGER, PARAMETER :: UVAR_CELL_rho = 7                           !Cell density
      INTEGER, PARAMETER :: UVAR_CELL_viab = 8                          !Cell Viability
      INTEGER, PARAMETER :: UVAR_NabD = 9                               !Solubility gradient term
      INTEGER, PARAMETER :: UVAR_SP = 10                                !Total Stress Pressure
      INTEGER, PARAMETER :: UVAR_IL1B = 11                              !IL1B
      INTEGER, PARAMETER :: UVAR_IL1B_PROT = 12                         !IL1B Protein
      INTEGER, PARAMETER :: UVAR_TNF = 13                               !TNF
      INTEGER, PARAMETER :: UVAR_TNF_PROT = 14                          !TNF Protein
      INTEGER, PARAMETER :: UVAR_AGG_IMN = 15                           !Agg IMN
      INTEGER, PARAMETER :: UVAR_AGG_IL1B = 16                          !Agg IL1B
      INTEGER, PARAMETER :: UVAR_AGG_TNF = 17                           !Agg TNF
      INTEGER, PARAMETER :: UVAR_AGG_IL1B_TNF = 18                      !Agg IL1B TNF
      INTEGER, PARAMETER :: UVAR_COLI_IMN = 19                          !Col I IMN
      INTEGER, PARAMETER :: UVAR_COLI_IL1B = 20                         !Col I IL1B
      INTEGER, PARAMETER :: UVAR_COLI_TNF = 21                          !Col I TNF
      INTEGER, PARAMETER :: UVAR_COLI_IL1B_TNF = 22                     !Col I IL1B TNF
      INTEGER, PARAMETER :: UVAR_COLII_IMN = 23                         !Col II IMN
      INTEGER, PARAMETER :: UVAR_COLII_IL1B = 24                        !Col II IL1B
      INTEGER, PARAMETER :: UVAR_COLII_TNF = 25                         !Col II TNF
      INTEGER, PARAMETER :: UVAR_COLII_IL1B_TNF = 26                    !Col II IL1B TNF
      INTEGER, PARAMETER :: UVAR_MMP3_IMN = 27                          !MMP3 IMN
      INTEGER, PARAMETER :: UVAR_MMP3_IL1B = 28                         !MMP3 IL1B
      INTEGER, PARAMETER :: UVAR_MMP3_TNF = 29                          !MMP3 TNF
      INTEGER, PARAMETER :: UVAR_MMP3_IL1B_TNF = 30                     !MMP3 IL1B TNF
      INTEGER, PARAMETER :: UVAR_ADM4_IMN = 31                          !ADAM4 IMN
      INTEGER, PARAMETER :: UVAR_ADM4_IL1B = 32                         !ADAM4 IL1B
      INTEGER, PARAMETER :: UVAR_ADM4_TNF = 33                          !ADAM4 TNF
      INTEGER, PARAMETER :: UVAR_ADM4_IL1B_TNF = 34                     !ADAM4 IL1B TNF
      INTEGER, PARAMETER :: UVAR_FREQ = 35                              !Frequency

C -----------------------------------------------------------
C Defining vector with UVARM numbers and names
      CHARACTER*256 :: UVARM_NAMES(35)

C -----------------------------------------------------------
C Defining variables for the frequency file
      INTEGER NActivities                                               ! max number of activity columns in the CSV
      INTEGER PN_NActivities                                            ! max number of activity columns in the CSV for PN step
      CHARACTER*256, ALLOCATABLE :: STEP_LIST_STR(:)                    ! e.g., "5-9-13": STEP_LIST_STR(NActivities)
      CHARACTER*256, ALLOCATABLE :: HOUR_RANGE_STR(:)                   ! e.g., "9-12": HOUR_RANGE_STR(NActivities)
      DOUBLE PRECISION, ALLOCATABLE :: FREQ_VAL(:)                      ! e.g., 0.1: FREQ_VAL(NActivities)
      CHARACTER*256, ALLOCATABLE :: PN_HOUR_RANGE_STR(:)                ! PN eq. e.g., "9-12": HOUR_RANGE_STR(NActivities)
      DOUBLE PRECISION, ALLOCATABLE :: PN_FREQ_VAL(:)                   ! PN EQ. e.g., 0.1: FREQ_VAL(NActivities)
      DOUBLE PRECISION FREQUENCY_CURRENT                                ! current frequency value
      INTEGER FREQfileSwitch                                            ! 0: not file, 1: file
      INTEGER FREQ_activation                                           ! 0: not activated, 1: activated
      INTEGER PN_TIME_file                                              ! PN eq. time defined in the file
      DOUBLE PRECISION, ALLOCATABLE :: PN_FREQ_file(:)                  ! PN eq. frequency defined in the file

C -----------------------------------------------------------
C Some constant parameters
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0, ONE=1.0D0, TWO=2.0D0
      DOUBLE PRECISION, PARAMETER :: THREE=3.0D0, FOUR=4.0D0, SIX=6.0D0
      DOUBLE PRECISION, PARAMETER :: HALF=0.5D0, FOURTH=0.25D0
      DOUBLE PRECISION, PARAMETER :: EIGHTH=0.125D0

C--------------------------------------------------------------------------
      END MODULE GLOBAL

C-------------------------------------------------------------------------------------------------------------------------------------------------------------------
C SUBROUTINES
C------------

C**********************************************************************************************************************************************
C UVARM NAMES subroutine
C**********************************************************************************************************************************************

C This subroutine is used to initialize the UVARM_NAMES array
C
      SUBROUTINE Initialize_UVARM_NAMES()
C
C ********************************************************************
C
      USE GLOBAL
C

C Initialize the UVARM_NAMES array
      UVARM_NAMES(UVAR_NF0) = 'Initial Water Content'
      UVARM_NAMES(UVAR_NF) = 'Water Content'
      UVARM_NAMES(UVAR_O2) = 'O2 Concentration'
      UVARM_NAMES(UVAR_Lact) = 'Lactate Concentration'
      UVARM_NAMES(UVAR_Gluc) = 'Glucose Concentration'
      UVARM_NAMES(UVAR_pH_val) = 'pH Value'
      UVARM_NAMES(UVAR_CELL_rho) = 'Cell Density'
      UVARM_NAMES(UVAR_CELL_viab) = 'Cell Viability'
      UVARM_NAMES(UVAR_NabD) = 'Solubility Gradient term'
      UVARM_NAMES(UVAR_SP) = 'Total Stress Pressure'
      UVARM_NAMES(UVAR_IL1B) = 'IL1B'
      UVARM_NAMES(UVAR_IL1B_PROT) = 'IL1B Protein'
      UVARM_NAMES(UVAR_TNF) = 'TNF'
      UVARM_NAMES(UVAR_TNF_PROT) = 'TNF Protein'
      UVARM_NAMES(UVAR_AGG_IMN) = 'Agg IMN'
      UVARM_NAMES(UVAR_AGG_IL1B) = 'Agg IL1B'
      UVARM_NAMES(UVAR_AGG_TNF) = 'Agg TNF'
      UVARM_NAMES(UVAR_AGG_IL1B_TNF) = 'Agg IL1B TNF'
      UVARM_NAMES(UVAR_COLI_IMN) = 'Col I IMN'
      UVARM_NAMES(UVAR_COLI_IL1B) = 'Col I IL1B'
      UVARM_NAMES(UVAR_COLI_TNF) = 'Col I TNF'
      UVARM_NAMES(UVAR_COLI_IL1B_TNF) = 'Col I IL1B TNF'
      UVARM_NAMES(UVAR_COLII_IMN) = 'Col II IMN'
      UVARM_NAMES(UVAR_COLII_IL1B) = 'Col II IL1B'
      UVARM_NAMES(UVAR_COLII_TNF) = 'Col II TNF'
      UVARM_NAMES(UVAR_COLII_IL1B_TNF) = 'Col II IL1B TNF'
      UVARM_NAMES(UVAR_MMP3_IMN) = 'MMP3 IMN'
      UVARM_NAMES(UVAR_MMP3_IL1B) = 'MMP3 IL1B'
      UVARM_NAMES(UVAR_MMP3_TNF) = 'MMP3 TNF'
      UVARM_NAMES(UVAR_MMP3_IL1B_TNF) = 'MMP3 IL1B TNF'
      UVARM_NAMES(UVAR_ADM4_IMN) = 'ADAM4 IMN'
      UVARM_NAMES(UVAR_ADM4_IL1B) = 'ADAM4 IL1B'
      UVARM_NAMES(UVAR_ADM4_TNF) = 'ADAM4 TNF'
      UVARM_NAMES(UVAR_ADM4_IL1B_TNF) = 'ADAM4 IL1B TNF'
      UVARM_NAMES(UVAR_FREQ) = 'Frequency'

C--------------------------------------------------------------------------
      RETURN
      END SUBROUTINE Initialize_UVARM_NAMES

C**********************************************************************************************************************************************
C UVARM subroutine input and workspace definitions
C**********************************************************************************************************************************************

C This subroutine is used to transfer SDV's from the UEL
C onto the dummy mesh for viewing.  Note that an offset of
C ElemOffset is used between the real mesh and the dummy mesh.
C If your model has more than ElemOffset UEL elements, then
C this will need to be modified.
C
      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     1 NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     2 JMAC,JMATYP,MATLAYO,LACCFLA)
C
C ********************************************************************
C
      USE GLOBAL
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(15)
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)
      INTEGER ELEM_Curr

C current element number
      ELEM_Curr = 0
      ELEM_Curr = NOEL + ElemOffset

C Transfer the SDV's from the UEL to the UVARM
      UVAR(UVAR_NF0) = globalSDVs(ELEM_Curr,NPT,UVAR_NF0)               !Initial water content
      UVAR(UVAR_NF) = globalSDVs(ELEM_Curr,NPT,UVAR_NF)                 !Water content
      UVAR(UVAR_O2) = globalSDVs(ELEM_Curr,NPT,UVAR_O2)                 !O2 concentration
      UVAR(UVAR_Lact) = globalSDVs(ELEM_Curr,NPT,UVAR_Lact)             !Lactate concentration
      UVAR(UVAR_Gluc) = globalSDVs(ELEM_Curr,NPT,UVAR_Gluc)             !Glucose concentration
      UVAR(UVAR_pH_val) = globalSDVs(ELEM_Curr,NPT,UVAR_pH_val)         !pH value
      UVAR(UVAR_CELL_rho) = globalSDVs(ELEM_Curr,NPT,UVAR_CELL_rho)     !Cell density
      UVAR(UVAR_CELL_viab) = globalSDVs(ELEM_Curr,NPT,UVAR_CELL_viab)   !Cell viability
      UVAR(UVAR_NabD) = globalSDVs(ELEM_Curr,NPT,UVAR_NabD)             !solubility gradient term
      UVAR(UVAR_SP) = globalSDVs(ELEM_Curr,NPT,UVAR_SP)                 !Total Stress Pressure

C PN_Eq SDV's
C TNF & IL1B
      UVAR(UVAR_IL1B) = globalSDVs(ELEM_Curr,NPT,UVAR_IL1B)             !IL1B
      UVAR(UVAR_IL1B_PROT) = globalSDVs(ELEM_Curr,NPT,UVAR_IL1B_PROT)   !IL1B Protein
      UVAR(UVAR_TNF) = globalSDVs(ELEM_Curr,NPT,UVAR_TNF)               !TNF
      UVAR(UVAR_TNF_PROT) = globalSDVs(ELEM_Curr,NPT,UVAR_TNF_PROT)     !TNF Protein

C Agg, Col I, Col II, MMP3, ADAM4 SDV's
      UVAR(UVAR_AGG_IMN) = globalSDVs(ELEM_Curr,NPT,UVAR_AGG_IMN)       !Agg IMN
      UVAR(UVAR_AGG_IL1B) = globalSDVs(ELEM_Curr,NPT,UVAR_AGG_IL1B)     !Agg IL1B
      UVAR(UVAR_AGG_TNF) = globalSDVs(ELEM_Curr,NPT,UVAR_AGG_TNF)       !Agg TNF
      UVAR(UVAR_AGG_IL1B_TNF) = globalSDVs(ELEM_Curr,NPT,               !Agg IL1B TNF
     1                          UVAR_AGG_IL1B_TNF)

      UVAR(UVAR_COLI_IMN) = globalSDVs(ELEM_Curr,NPT,UVAR_COLI_IMN)     !Col I IMN
      UVAR(UVAR_COLI_IL1B) = globalSDVs(ELEM_Curr,NPT,UVAR_COLI_IL1B)   !Col I IL1B
      UVAR(UVAR_COLI_TNF) = globalSDVs(ELEM_Curr,NPT,UVAR_COLI_TNF)     !Col I TNF
      UVAR(UVAR_COLI_IL1B_TNF) = globalSDVs(ELEM_Curr,NPT,              !Col I IL1B TNF
     1                           UVAR_COLI_IL1B_TNF)

      UVAR(UVAR_COLII_IMN) = globalSDVs(ELEM_Curr,NPT,UVAR_COLII_IMN)   !Col II IMN
      UVAR(UVAR_COLII_IL1B) = globalSDVs(ELEM_Curr,NPT,UVAR_COLII_IL1B) !Col II IL1B
      UVAR(UVAR_COLII_TNF) = globalSDVs(ELEM_Curr,NPT,UVAR_COLII_TNF)   !Col II TNF
      UVAR(UVAR_COLII_IL1B_TNF) = globalSDVs(ELEM_Curr,NPT,             !Col II IL1B TNF
     1                            UVAR_COLII_IL1B_TNF)

      UVAR(UVAR_MMP3_IMN) = globalSDVs(ELEM_Curr,NPT,UVAR_MMP3_IMN)     !MMP3 IMN
      UVAR(UVAR_MMP3_IL1B) = globalSDVs(ELEM_Curr,NPT,UVAR_MMP3_IL1B)   !MMP3 IL1B
      UVAR(UVAR_MMP3_TNF) = globalSDVs(ELEM_Curr,NPT,UVAR_MMP3_TNF)     !MMP3 TNF
      UVAR(UVAR_MMP3_IL1B_TNF) = globalSDVs(ELEM_Curr,NPT,              !MMP3 IL1B TNF
     1                           UVAR_MMP3_IL1B_TNF)

      UVAR(UVAR_ADM4_IMN) = globalSDVs(ELEM_Curr,NPT,UVAR_ADM4_IMN)     !ADAM4 IMN
      UVAR(UVAR_ADM4_IL1B) = globalSDVs(ELEM_Curr,NPT,UVAR_ADM4_IL1B)   !ADAM4 IL1B
      UVAR(UVAR_ADM4_TNF) = globalSDVs(ELEM_Curr,NPT,UVAR_ADM4_TNF)     !ADAM4 TNF
      UVAR(UVAR_ADM4_IL1B_TNF) = globalSDVs(ELEM_Curr,NPT,              !ADAM4 IL1B TNF
     1                            UVAR_ADM4_IL1B_TNF)
C Frequency
      UVAR(UVAR_FREQ) = globalSDVs(ELEM_Curr,NPT,UVAR_FREQ)             !Frequency

C--------------------------------------------------------------------------
      RETURN
      END SUBROUTINE UVARM

C**********************************************************************************************************************************************
C UEXTERNALDB input and workspace definitions
C**********************************************************************************************************************************************

C Subroutine to create a information file
C
      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
C
C ********************************************************************
C
      USE GLOBAL
C
      INCLUDE 'ABA_PARAM.INC'
C
C Define the COMMON block for the info file
      COMMON /comm_info/ infoFilePath_T_D,cwd,outdir
C Define info file name
      CHARACTER*256 cwd,outdir,infoFile_T_D,infoFilePath_T_D,JOBNAME
C Define input file name
      CHARACTER*256 inpFile, inpFilePath
C SRC: Path to the transport source folder
      CHARACTER*256 srcFolderTransp
C SDV file name
C NF0: Initial water content
      CHARACTER*256 sdvNF0File,sdvNF0FilePath
C SP: Total Stress Pressure
      CHARACTER*256 sdvSPFile, sdvSPFolder, sdvSPFilePath
      CHARACTER*10  kstep_str, kinc_str
C
      DIMENSION TIME(2)
      INTEGER error, LENJOBNAME
C For SDVs
      INTEGER iElem,iIP
      DOUBLE PRECISION SDV
      integer :: STATUS = 0
      INTEGER uniqueSteps(100)
      INTEGER nstepsSP
      INTEGER ngSdv
      CHARACTER*256 stepList, tempStr
C For Frequency
      CHARACTER*256 FREQfile, FREQfilePath
      DOUBLE PRECISION currentStepHour
      INTEGER nstepToken, COUNT_TOKENS
      LOGICAL isStep, isHour, STEP_APPLIES, HOUR_APPLIES
      EXTERNAL COUNT_TOKENS, STEP_APPLIES, HOUR_APPLIES

C
C--------------------------------------------------------------------
C Definition of paths:
C--------------------------------------------------------------------
C

C Write a file with the information of the simulation
      IF (LOP .EQ. 0) THEN                                              !start of the analysis

C Call DATE_AND_TIME
      CALL DATE_AND_TIME(values=date_ini)

C Get the job name
      CALL GETJOBNAME(JOBNAME,LENJOBNAME)
      CALL GETCWD(cwd)
      CALL GETOUTDIR(outdir,lenoutdir)
      lenoutdir = 50
      infoFile_T_D = '/' // TRIM(JOBNAME) // '.info'
      infoFilePath_T_D = TRIM(outdir) // TRIM(infoFile_T_D)
C Get the input file name
      inpFile = '/' // TRIM(JOBNAME) // '.inp'
      inpFilePath = TRIM(outdir) // TRIM(inpFile)

C Get the SDV file or folder name from previous simulation
C NF0: Initial water content
      sdvNF0File = '/' // TRIM(JOBNAME) // '_SDV.NF0'
      sdvNF0FilePath = TRIM(outdir) // TRIM(sdvNF0File)

C SP: Total Stress Pressure
      sdvSPFolder = TRIM(outdir) // '/' // TRIM(JOBNAME) // '_SP/'
      stepList = ''
C Transport folder
      srcFolderTransp = TRIM(outdir) // '/src/IVD/transport/'

C--------------------------------------------------------------------------
C Store the name step in the variable STEPNAMES(NSTEPS), which is allocatable
C Get the NSTEPS from the input file
      CALL GETNSTEPS(NSTEPS, inpFilePath)
      ALLOCATE(STEPNAMES(NSTEPS))
C Get the step names from the input file
      CALL GETSTEPNAMES(STEPNAMES, NSTEPS, inpFilePath)

C--------------------------------------------------------------------------
C Create the info file
      OPEN(15,FILE=infoFilePath_T_D, STATUS='UNKNOWN')
      WRITE(15,*) ' '
      CLOSE(15)

C--------------------------------------------------------------------------
C Write and print the simulation info file
      CALL Initialize_UVARM_NAMES()
      ngSdv = SIZE(UVARM_NAMES)

C--------------------------------------------------------------------------
C Frequency file
      CALL GETENV('FREQfile', FREQfile)
      FREQ_activation = 0
      PN_TIME_file = 1
      IF (FREQfile .NE. '') THEN
            FREQfileSwitch = 1
            FREQfilePath = TRIM(srcFolderTransp) // 'frequency/'
            FREQfilePath = TRIM(FREQfilePath) // TRIM(FREQfile)
C Count the number of activities in the CSV file (counting the commas)
            IF (.NOT. ALLOCATED(STEP_LIST_STR)) THEN
                  CALL NcommasCSV(FREQfilePath, NActivities, 9)
                  CALL NcommasCSV(FREQfilePath, PN_NActivities, 4)
                  NActivities = NActivities + 1
                  PN_NActivities = PN_NActivities + 1
C Allocate the arrays for the frequency file
                  ALLOCATE(STEP_LIST_STR(NActivities))
                  ALLOCATE(HOUR_RANGE_STR(NActivities))
                  ALLOCATE(FREQ_VAL(NActivities))
C for PN step
                  ALLOCATE(PN_HOUR_RANGE_STR(PN_NActivities))
                  ALLOCATE(PN_FREQ_VAL(PN_NActivities))
            END IF
C Read the frequency file to get the activities
            CALL readActivitiesCSV(FREQfilePath, NActivities,
     1                             STEP_LIST_STR, HOUR_RANGE_STR,
     2                             FREQ_VAL, PN_TIME_file,
     3                             PN_NActivities,
     4                             PN_HOUR_RANGE_STR, PN_FREQ_VAL)
      ELSE
            FREQfileSwitch = 0
      END IF
C Define PN_FREQ_file allocated with PN_TIME_file
      ALLOCATE(PN_FREQ_file(PN_TIME_file))
C Fill it with zeros
      PN_FREQ_file = ZERO

      IF (FREQfile .NE. '') THEN
            DO I = 1, PN_NActivities
                  CALL PN_HOURS(PN_FREQ_file, PN_TIME_file,
     1                          PN_HOUR_RANGE_STR(I), PN_FREQ_VAL(I))
            END DO
      END IF

C --------------------------
      END IF

C Write and print the simulation info file
      INCLUDE './src/IVD/transport/functions/write_info_Transport.f'

C--------------------------------------------------------------------------
C Read the initial values of the SDV's
C the file name should be JOBNAME + '_sdv.' + variable

      IF (LOP .EQ. 0) THEN                                              !start of the analysis
      OPEN(15,FILE=infoFilePath_T_D, STATUS='OLD', POSITION='APPEND')
      WRITE(15,*) '****************************************************'
      WRITE(15,*) '                       SDVs'
      WRITE(15,*) '****************************************************'
      WRITE(15,*) ' '
      WRITE(15,*) 'Reading the initial values of the SDVs:'
      WRITE(15,*) ' '
      WRITE(15,*) 'Initial water content (NF0)'
      WRITE(15,*) 'File: ', sdvNF0FilePath
      WRITE(15,*) ' '
      WRITE(15,*) 'Cheking temporal values of the SDVs:'
      WRITE(15,*) ' '
      WRITE(15,*) 'Total Stress Pressure (SP)'
      WRITE(15,*) 'Folder: ', sdvSPFolder
      WRITE(15,*) 'Steps that will include SP from previous simulation:'
      CALL GETUNIQUESTEPS(sdvSPFolder, JOBNAME,
     1                    outdir, uniqueSteps, nstepsSP)
      ! Debug prints: Check nstepsSP and each uniqueSteps value.
      IF (nstepsSP .EQ. 0) THEN
      WRITE(15,*) ' '
      WRITE(15,*) '////////////////////////////////////////////////////'
      WRITE(15,*) 'EXCEPTION: The folder does not exist or the files'
      WRITE(15,*) '           could not be read.'
      WRITE(15,*) 'OBS: As the file could not be read,'
      WRITE(15,*) '     be sure the variable is set as ZERO.'
      WRITE(15,*) '////////////////////////////////////////////////////'
      ELSE
      DO I = 1, nstepsSP
            WRITE(tempStr, '(I3)') uniqueSteps(I)
            tempStr = ADJUSTL(tempStr)
            IF (I .EQ. 1) THEN
                  stepList = TRIM(tempStr)
            ELSE
                  stepList = TRIM(stepList) // ', ' // TRIM(tempStr)
            END IF
      END DO
      WRITE(15, '(A)') stepList
      END IF
      WRITE(15,*) ' '
      CLOSE(15)

      PRINT *, '****************************************************'
      PRINT *, '                       SDVs'
      PRINT *, '****************************************************'
      PRINT *, ' '
      PRINT *, 'Reading the initial values of the SDVs:'
      PRINT *, ' '
      PRINT *, 'Initial water content (NF0)'
      PRINT *, 'File: ', sdvNF0FilePath
      PRINT *, ' '
      PRINT *, 'Cheking temporal values of the SDVs:'
      PRINT *, ' '
      PRINT *, 'Total Stress Pressure (SP)'
      PRINT *, 'Folder: ', sdvSPFolder
      PRINT *, 'Steps that will include SP from previous simulation:'
      ! Debug prints: Check nstepsSP and each uniqueSteps value.
      IF (nstepsSP .EQ. 0) THEN
      PRINT *, ' '
      PRINT *, '////////////////////////////////////////////////////'
      PRINT *, 'EXCEPTION: The folder does not exist or the files'
      PRINT *, '           could not be read.'
      PRINT *, 'OBS: As the file could not be read,'
      PRINT *, '     be sure the variable is set as ZERO.'
      PRINT *, '////////////////////////////////////////////////////'
      ELSE
      DO I = 1, nstepsSP
            WRITE(tempStr, '(I3)') uniqueSteps(I)
            tempStr = ADJUSTL(tempStr)
            IF (I .EQ. 1) THEN
                  stepList = TRIM(tempStr)
            ELSE
                  stepList = TRIM(stepList) // ', ' // TRIM(tempStr)
            END IF
      END DO
      PRINT *, stepList
      END IF
      PRINT *, ' '

C NF0: Initial water content
C --------------------------
C read the file and extract the values for every integration point
      OPEN(15,FILE=sdvNF0FilePath,STATUS='OLD',FORM='FORMATTED'
     1,IOSTAT=error)

      IF (error .EQ. 0) THEN
            iElem = 0
            iIP = 0
            SDV = ZERO
            DO I = 1,numSDVxIdElem
                  DO J = 1,numInt
                        SDV_NF0(I,J) = ZERO
                  END DO
            END DO
                        
            DO I = 1,numElemInt
C Store the values ELEM, IP, SDV
                  READ(15,*,IOSTAT=error) iElem,iIP,SDV
                  SDV_NF0(iElem,iIP) = SDV
            END DO
            CLOSE(15)
      ELSE
      CLOSE(15)
      OPEN(15,FILE=infoFilePath_T_D, STATUS='OLD', POSITION='APPEND')
      WRITE(15,*) ' '
      WRITE(15,*) '////////////////////////////////////////////////////'
      WRITE(15,*) 'EXCEPTION: Could not open file:'
      WRITE(15,*) sdvNF0FilePath
      WRITE(15,*) 'VAR: ', error
      WRITE(15,*) 'STEP: ', KSTEP
      WRITE(15,*) 'INC: ', KINC
      WRITE(15,*) 'OBS: As the file could not be read, be sure to'
      WRITE(15,*) '     implement this variable as input in the PROPS.'
      WRITE(15,*) '////////////////////////////////////////////////////'
      WRITE(15,*) ' '
      WRITE(15,*) '****************************************************'
      WRITE(15,*) ' '
      CLOSE(15)

      PRINT *, ' '
      PRINT *, '////////////////////////////////////////////////////'
      PRINT *, 'EXCEPTION: Could not open file:'
      PRINT *, sdvNF0FilePath
      PRINT *, 'VAR: ', error
      PRINT *, 'STEP: ', KSTEP
      PRINT *, 'INC: ', KINC
      PRINT *, 'OBS: As the file could not be read, be sure to'
      PRINT *, '     implement this variable as input in the PROPS.'
      PRINT *, '////////////////////////////////////////////////////'
      PRINT *, ' '
      PRINT *, '****************************************************'
      PRINT *, ' '
C      CALL XIT
      END IF

C Allocate memory for the SDVs (global and local)
      INCLUDE './src/IVD/transport/functions/UEL_SDVs_Allocate.f'

C --------------------------
      END IF

C--------------------------------------------------------------------------
C Get the temporal SDV's
      IF (LOP .EQ. 1) THEN                                              ! start of the current analysis increment.
C Get the job name
      CALL GETJOBNAME(JOBNAME,LENJOBNAME)
      CALL GETCWD(cwd)
      CALL GETOUTDIR(outdir,lenoutdir)

C SP: Total Stress Pressure
C --------------------------
C Initialize strings for conversion
      kstep_str = '          '
      kinc_str  = '          '
C Convert integers KSTEP and KINC to strings using format I10.
      WRITE(kstep_str, '(I10)') KSTEP
      WRITE(kinc_str,  '(I10)') KINC
      kstep_str = ADJUSTL(kstep_str)
      kinc_str  = ADJUSTL(kinc_str)

C Now concatenate the strings. (TRIM is nonstandard in F77 but often available as an extension)
      sdvSPFolder = '/' // TRIM(JOBNAME) // '_SP'
      sdvSPFile = '/' // TRIM(JOBNAME) // '_SDV.' // TRIM(kstep_str)
      sdvSPFile = TRIM(sdvSPFile) // '.' // TRIM(kinc_str) // '.SP'

      sdvSPFilePath = TRIM(outdir)//TRIM(sdvSPFolder)//TRIM(sdvSPFile)

C read the file and extract the values for every integration point
      OPEN(15,FILE=sdvSPFilePath,STATUS='OLD',FORM='FORMATTED'
     1,IOSTAT=error)

      IF (error .EQ. 0) THEN
            iElem = 0
            iIP = 0
            SDV = ZERO
            DO I = 1,numSDVxIdElem
                  DO J = 1,numInt
                        SDV_SP(I,J) = ZERO
                  END DO
            END DO
                        
            DO I = 1,numElemInt
C Store the values ELEM, IP, SDV
                  READ(15,*,IOSTAT=error) iElem,iIP,SDV
                  SDV_SP(iElem,iIP) = SDV
            END DO
            CLOSE(15)
      ELSE
      CLOSE(15)
            DO I = 1,numSDVxIdElem
                  DO J = 1,numInt
                        SDV_SP(I,J) = ZERO
                  END DO
            END DO
      END IF

C FREQUENCY_CURRENT: Current frequency value
C ------------------------------------------
C If we have a frequency file, interpret it now based on the current step/time
      IF (FREQfileSwitch .EQ. 1) THEN
      FREQ_activation = 0
      FREQUENCY_CURRENT = ZERO
      currentStepHour = TIME(1) / 3600.0D0

      DO i = 1, NActivities
      nstepToken = COUNT_TOKENS(STEP_LIST_STR(i), '-')
      IF (nstepToken .GT. 0) THEN
            isStep = STEP_APPLIES(KSTEP, STEP_LIST_STR(i), nstepToken)
            isHour = HOUR_APPLIES(currentStepHour, HOUR_RANGE_STR(i))
            IF ( isStep .AND. isHour ) THEN
                  FREQUENCY_CURRENT = FREQ_VAL(i)
                  FREQ_activation = 1
                  EXIT
            END IF
      ELSE
            PRINT *, 'ERROR: The frequency file is not well defined.'
            PRINT *, '       Check the number of steps and hours.'
            PRINT *, '       The file will not be read.'
            PRINT *, '       The frequency will be set to zero.'
            PRINT *, ' '
      END IF

      END DO
      END IF

C --------------------------
      END IF


C--------------------------------------------------------------------------
      RETURN
      END SUBROUTINE UEXTERNALDB

C**********************************************************************************************************************************************
C UMAT subroutine input and workspace definitions
C**********************************************************************************************************************************************
C
C UMAT runs for every integration point in the mesh
C
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,
     2 PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,
     3 DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
C ********************************************************************
C
      USE GLOBAL
C
      INCLUDE 'ABA_PARAM.INC'
C
C Define the COMMON block for the diffusion-reaction transport model
      COMMON /comm_share/ run_count_share
      COMMON /comm_info/ infoFilePath_T_D,cwd,outdir
C
      DIMENSION STRESS(NTENS),STATEV(NSTATEV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JAY(NSTATEV),DFGR(3,3)
      CHARACTER*80 MATNAME
      CHARACTER*256 CMNAME, cwd, outdir, infoFile_T_D, infoFilePath_T_D
      INTEGER run_count_share

C
C--------------------------------------------------------------------------
C Variable definitions for the problem
C--------------------------------------------------------------------------
C
      DOUBLE PRECISION DETDFG,CELL_rho, CELL_viab
      DOUBLE PRECISION DFGRD1_inv(NDIM,NDIM),DFGRD1_invT(NDIM,NDIM)
      DOUBLE PRECISION CC,B1,D1,pH_ini,AA,BB,DD
      DOUBLE PRECISION summ1,summ2,exg
      DOUBLE PRECISION ti, deti
      DOUBLE PRECISION glu, NF
      DOUBLE PRECISION TRACE_D, dD_SOLdx

C ALPHA_pH: Death rate due to acidity
C ALPHA_GLUC: Death rate due to glucose.
C pH_val: current pH value
      DOUBLE PRECISION ALPHA_pH, ALPHA_GLUC, pH_val

C Oxygen solubility in water
      DOUBLE PRECISION, PARAMETER :: SO2 = 0.010268D0

C Define constants and parameters
C I, J, M, N: Integer variables
      INTEGER I, J, II, JJ, M, N
      INTEGER stat
      DOUBLE PRECISION D_LAMDA, D_MU, diff, adiff, NF0

C--------------------------------------------------------------------
C    PROPS(1) - DUMMY YOUNG MODULUS [MPa]
C    PROPS(2) - POISSON COEFICIENT (0.3)
C    PROPS(3) - Correction constant for glucose rate
C    PROPS(4) - Glucose rate parameter 1
C    PROPS(5) - Glucose rate parameter 2
C    PROPS(6) - Initial fraction of water content homogeneous
C    PROPS(7) - Initial cell density
C--------------------------------------------------------------------
      DOUBLE PRECISION E, NU, ALPHA_GCF, kk1, kk2, NF0h, CELL_rho_0

C Initialization of variables
C ---------------------------
      DETDFG = ZERO
      CELL_rho = ZERO
      CELL_viab = ZERO
      DFGRD1_inv = ZERO
      DFGRD1_invT = ZERO
      CC = ZERO
      B1 = ZERO
      D1 = ZERO
      pH_ini = ZERO
      AA = ZERO
      BB = ZERO
      DD = ZERO
      summ1 = ZERO
      summ2 = ZERO
      exg = ZERO
      ti = ZERO
      deti = ZERO
      glu = ZERO
      NF = ZERO
      ALPHA_pH = ZERO
      ALPHA_GLUC = ZERO
      pH_val = ZERO
      diff = ZERO
      adiff = ZERO
      NF0 = ZERO
      D_LAMDA = ZERO
      D_MU = ZERO
      TRACE_D = ZERO
      dD_SOLdx = ZERO

C PROPS variables
      E = ZERO
      NU = ZERO
      ALPHA_GCF = ZERO
      kk1 = ZERO
      kk2 = ZERO
      NF0h = ZERO
      CELL_rho_0 = ZERO

C Definitions
      E=PROPS(1)
      NU=PROPS(2)
      ALPHA_GCF = PROPS(3)
      kk1=PROPS(4)
      kk2=PROPS(5)
      NF0h = PROPS(6)
      CELL_rho_0 = PROPS(7)
      MATNAME = CMNAME(4:5)

C ACTIVATE ONLY TO COMPARE O2 BETWEEN UMAT (SHOWN AS TEMPERATURE) AND UEL
      INCLUDE './src/IVD/transport/functions/UMAT_verification.f'

      run_count_share = run_count_share + 1
C--------------------------------------------------------------------------
      RETURN
      END SUBROUTINE UMAT

C**********************************************************************************************************************************************
C UMATHT subroutine for dummy elements (only to compare 02 (temperature here) with UEL)
C**********************************************************************************************************************************************

C
      SUBROUTINE UMATHT(U,DUDT,DUDG,FLUX,DFDT,DFDG,STATEV,TEMP,
     +     DTEMP,DTEMDX,TIME,DTIME,PREDEF,DPRED,CMNAME,NTGRD,NSTATV,
     +     PROPS,NPROPS,COORDS,PNEWDT,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
C********************************************************************
C
      USE GLOBAL
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
C
      DIMENSION DUDG(NTGRD),FLUX(NTGRD),DFDT(NTGRD),
     1     DFDG(NTGRD,NTGRD),STATEV(NSTATV),DTEMDX(NTGRD),TIME(2),
     2     PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3),TEMP(NTGRD)

C Definitions
      DOUBLE PRECISION COND, COND2, SPECHT
C
C--------------------------------------------------------------------------
C Variable definitions for the problem
C--------------------------------------------------------------------------
C

C ACTIVATE ONLY TO COMPARE O2 BETWEEN UMAT (SHOWN AS TEMPERATURE) AND UEL
      INCLUDE './src/IVD/transport/functions/UMATHT_verification.f'

C--------------------------------------------------------------------------
      RETURN
      END SUBROUTINE UMATHT

C**********************************************************************************************************************************************
C UEL subroutine input and workspace definitions
C**********************************************************************************************************************************************

C This code is an ABAQUS User Element (UEL) which describes the 
C diffusion-reaction transport of solutes in an intervertebral disc (IVD).
C The solutes considered: Oxygen (O2), Lactate (lact), and Glucose (gluc).
C
C UEL runs for every element in the mesh
C
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +     PROPS,NPROPS,COORDS,MCRD,NNODE,UALL,dUALL,Vel,Accn,JTYPE,
     +     TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +     PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +     NJPROP,PERIOD)
C
C ********************************************************************
C
      USE GLOBAL
C
      IMPLICIT NONE
C      INCLUDE 'ABA_PARAM.INC'
C
C Define the COMMON block for the diffusion-reaction transport model
      COMMON /comm_info/ infoFilePath_T_D,cwd,outdir
      COMMON /comm_share/ run_count_share
C
C     VARIABLES DEFINED IN UEL, PASSED BACK TO ABAQUS

      DOUBLE PRECISION :: RHS,AMATRX,SVARS,ENERGY

C     VARIABLES PASSED INTO UEL 
      DOUBLE PRECISION :: PROPS,COORDS,UALL,dUALL,Vel,Accn,TIME,
     1  DTIME,PARAMS,ADLMAG,PREDEF,DDLMAG,PNEWDT,PERIOD
      INTEGER :: NDOFEL,NRHS,NSVARS,NPROPS,MCRD,NNODE,JTYPE,KSTEP,KINC,
     1  JELEM,NDLOAD,JDLTYP,NPREDF,LFLAGS,MLVARX,MDLOAD,JPROPS,NJPROP
C
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),UALL(NDOFEL),
     2 dUALL(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

      CHARACTER*80 MATNAME
      CHARACTER*256 infoFilePath_T_D,cwd,outdir

C
C--------------------------------------------------------------------------
C Some necessary variables
C--------------------------------------------------------------------------
C
C NDOFN: Total number of degrees of freedom per node
      INTEGER :: NDOFN

C IdenM: Identity tensor
      DOUBLE PRECISION IdenM(3,3)

C IdenV: Identity vector
      DOUBLE PRECISION IdenV(3,1)

C Cip: Coordinates of the Gauss integration points in a reference or "master" 3D element
C whtG: The corresponding Gauss integration weights
C WHT : Gauss weight multiplied by jacobi of the transformation
      DOUBLE PRECISION Cip(numInt,NDIM), whtG(numInt), WHT

C DFG: Deformation gradient
C DFG_old: Deformation gradient at the previous time step
C DFG_inv: Inverse of the deformation gradient
C DFG_invT: Transpose of the inverse of the deformation gradient
C DETDFG: Determinant of the deformation gradient
C dDFGdx: Derivative of the deformation gradient with respect to the spatial coordinates
      DOUBLE PRECISION DFG(NDIM,NDIM),DFG_old(NDIM,NDIM)
      DOUBLE PRECISION DFG_inv(NDIM,NDIM),DFG_invT(NDIM,NDIM)
      DOUBLE PRECISION DETDFG, dDFGdx(NDIM,NDIM,NDIM)

C sh(nNode)
C N: Shape functions
      DOUBLE PRECISION N(NNODE)

C for ξ-η-ζ (isoparametric) domain:
C dNxi: Derivatives of shape functions
      DOUBLE PRECISION dNxi(NNODE,NDIM)
C d2Nxi: Second derivatives of shape functions
      DOUBLE PRECISION d2Nxi(NNODE,NDIM,NDIM)

C for x-y-z (physical) domain:
C dN: Derivatives of shape functions
C d2N: Second derivatives of shape functions
C DETMapJ: Determinant of the Jacobian matrix
      DOUBLE PRECISION dN(NNODE,NDIM), d2N(NNODE,NDIM,NDIM)
      DOUBLE PRECISION DETMapJ

C for centroid...
C N0: Shape functions
C dN0: Derivatives of shape functions
C d2N0: Second derivatives of shape functions
C Cip0: Coordinates of the centroid
C DETMapJ0: Determinant of the Jacobian matrix at the centroid
      DOUBLE PRECISION N0(NNODE), dN0(NNODE,NDIM), d2N0(NNODE,NDIM,NDIM)
      DOUBLE PRECISION Cip0(NDIM),DETMapJ0

C for current...
C dNC: Derivatives of shape functions
C d2NC: Second derivatives of shape functions
C d2Nmix: Mixed second derivatives of shape functions
C DETMapJC: Determinant of the Jacobian matrix at the current point
      DOUBLE PRECISION dNC(NNODE,NDIM), d2NC(NNODE,NDIM,NDIM)
      DOUBLE PRECISION d2Nmix(NNODE,NDIM,NDIM)
      DOUBLE PRECISION DETMapJC

C for current centroid...
C dN0C: Derivatives of shape functions at the centroid at the current point
C d2N0C: Second derivatives of shape functions at the centroid at the current point
C DETMapJ0C: Determinant of the Jacobian matrix at the centroid at the current point
      DOUBLE PRECISION dN0C(NNODE,NDIM), d2N0C(NNODE,NDIM,NDIM)
      DOUBLE PRECISION DETMapJ0C

C Shape functions vector      
      DOUBLE PRECISION Nvec(1,NNODE) !Shape functions vector

C RHS and AMATRX: Right hand side and stiffness matrix
C O2
      DOUBLE PRECISION R_O2(NNODE,1),K_O2(NNODE,NNODE)
      DOUBLE PRECISION K_O2LACT(NNODE,NNODE),K_O2GLUC(NNODE,NNODE)
C Lactate
      DOUBLE PRECISION R_LACT(NNODE,1),K_LACT(NNODE,NNODE)
      DOUBLE PRECISION K_LACTO2(NNODE,NNODE),K_LACTGLUC(NNODE,NNODE)
C Glucose
      DOUBLE PRECISION R_GLUC(NNODE,1),K_GLUC(NNODE,NNODE)
      DOUBLE PRECISION K_GLUCO2(NNODE,NNODE),K_GLUCLACT(NNODE,NNODE)

C Displacements
      DOUBLE PRECISION R_U(NDIM*NNODE,1) !RHS of displacements (dummy)
      DOUBLE PRECISION K_U(NDIM*NNODE,NDIM*NNODE) !Stiffness matrix of displacements (dummy)

C Residual and Tangent Factors
      DOUBLE PRECISION RF_O2, RF_LACT, RF_GLUC
      DOUBLE PRECISION TF_O2, TF_LACT, TF_GLUC

C
C--------------------------------------------------------------------------
C Variable definitions for the problem
C--------------------------------------------------------------------------
C

C Loop and element number variables
C ---------------------------------
C I, J, K: Integer variables
      INTEGER I, J, K, R, M, Ni, A, JJ

C index: Index of the integration point (1-numElemInt)
C i_intPT: Integration point (1-27)
C lenoutdir: Length of the output directory
      INTEGER index, i_intPT, lenoutdir

C run_count: Static variable to keep track of run count
      INTEGER, SAVE :: i_numElem = 1
      INTEGER run_count_share

C PROPS:
C ------
C 1: E: Young's modulusC
C 2: NU: Poisson's ratio
C 3,4,5: D_WAT: Solute water diffusivities
C 6: ph_val_0: Initial ph value
C 7: CELL_rho_0: Initial cell density
C 8: CELL_viab_0: Initial Cell Viability
C a_K, b_K, c_K: Constants for the K(CELL_rho_0) model
C 10: ALPHA_pH: Death rate due to acidity
C 11: GlucThres: Glucose threshold for cell viability
C 12: pHThres: pH threshold for cell viability
C 13: NF0h: Alternative valuees for the SDVs - Homogeneous Initial water content
      DOUBLE PRECISION E, NU, D_WAT(3)
      DOUBLE PRECISION ph_val_0, CELL_rho_0, CELL_viab_0
      DOUBLE PRECISION a_K, b_K, c_K, ALPHA_pH
      DOUBLE PRECISION GlucThres, pHThres, NF0h

C PN equation Props:
C Previous definition:
C Cell activity (CA): Cell responses (according to research interest / research question),
C                     e.g. in terms of mRNA expression or protein synthesis.
C Cell state (CS): Stimuli (S) coming from an external source that globally regulate
C                  targeted CA. Such stimuli might be mechanical
C                  loads, nutritional factors, radiation etc.
C Variables:
C FREQ: Frequency of load
C PN_TIME: Total desired protein expresion simulation time (hours)
      DOUBLE PRECISION FREQ
      INTEGER PN_TIME

C I PROPS:
C --------
C CNUMBER: Tissue Number
C ngSdv: Number of global SDV's per integration point
C nlSdv: Number of local sdv's per integ pt
C NF0h_d: Use NF0h (Homogeneous NF0)? yes NF0h_d = 1, no NF0h_d = 0
C nD_d: Use ∇D as ZERO (nD_d=0), Isotropic (nD_d=1), Anisotropic (nD_d=2)
      INTEGER CNUMBER, ngSdv, nlSdv, NF0h_d, nD_d

C SDVs from previous sim:
C -----------------------
C NF0: Initial water content
      DOUBLE PRECISION NF0
C SP: Total Stress Pressure
      DOUBLE PRECISION SP

C Variables for the UEL diffusion-reaction transport model
C --------------------------------------------------------

C Initial mechano-transport variables
C U: Current nodal displacements
C dU: Incremental nodal displacements
C U_old: Previous nodal displacements
      DOUBLE PRECISION U(NNODE,NDIM), dU(NNODE,NDOFEL)
      DOUBLE PRECISION U_old(NNODE,NDOFEL)

C COORDSC: Current nodal coordinates
      DOUBLE PRECISION COORDSC(NDIM,NNODE)

C Le: Length of the diagonal of the element
      DOUBLE PRECISION Le

C Concentrations for O2, lactate, and glucose
C Old
      DOUBLE PRECISION O2_CONC_old(NNODE),LACT_CONC_old(NNODE)
      DOUBLE PRECISION GLUC_CONC_old(NNODE)

C Current
      DOUBLE PRECISION O2_CONC(NNODE),LACT_CONC(NNODE)
      DOUBLE PRECISION GLUC_CONC(NNODE)

C Rate of change
      DOUBLE PRECISION dO2_CONC(NNODE),dLACT_CONC(NNODE)
      DOUBLE PRECISION dGLUC_CONC(NNODE)

C Initial concentrations
      DOUBLE PRECISION O2_CONC_0(NNODE),LACT_CONC_0(NNODE)
      DOUBLE PRECISION GLUC_CONC_0(NNODE)

C Time rate of change of the concentrations
      DOUBLE PRECISION dO2_CONCdt(NNODE),dLACT_CONCdt(NNODE)
      DOUBLE PRECISION dGLUC_CONCdt(NNODE)

C Integration point
      DOUBLE PRECISION O2_int,O2_int_old,dO2_dt
      DOUBLE PRECISION dO2dx(NDIM,1), dO2dt
      DOUBLE PRECISION LACT_int,LACT_int_old,dLACT_dt
      DOUBLE PRECISION dLACTdx(NDIM,1), dLACTdt
      DOUBLE PRECISION GLUC_int,GLUC_int_old,dGLUC_dt
      DOUBLE PRECISION dGLUCdx(NDIM,1), dGLUCdt

C Initial Integration point
      DOUBLE PRECISION O2_int_0,LACT_int_0,GLUC_int_0

C D_SOL: Diffusion matrix for solutes within the medium
C R_SOL: Reaction matrix for solutes within the medium
      DOUBLE PRECISION D_SOL(3), R_SOL(3)

C Gradient of the diffusion
      DOUBLE PRECISION dD_O2dx(NDIM,1)
      DOUBLE PRECISION dD_LACTdx(NDIM,1)
      DOUBLE PRECISION dD_GLUCdx(NDIM,1)
      DOUBLE PRECISION dD_SOLdx_MAG

C DERIVATIVE OF REACTION WRT TO O2, LACTATE, AND GLUCOSE
      DOUBLE PRECISION dR_SOL(3,3)

C SDVs Shared variables to be stored in UVARM
C -------------------------------------------
C stat: Status of the element
C NF: water content
C pH_val: pH value
C CELL_rho: Cell Density
C CELL_viab: Cell Viability

      INTEGER stat
      DOUBLE PRECISION NF, pH_val, CELL_rho, CELL_viab

C Variables for the UEL PN equation aproach
C -----------------------------------------
C Previous definition:
C Result of the PN-Equation, i.e. a final CA of a given CS within the system of interest.
C
C Variables:
C Proinflammatory cytokines (PIC): Related to initiations of IVD degeneration,
C                                  as they importantly alter the targeted CA of NP cells
C IL1B: Interleukin 1 beta
C IL1B_PROT: Interleukin 1 beta protein
C TNF: Tumor necrosis factor
C TNF_PROT: Tumor necrosis factor protein
      DOUBLE PRECISION IL1B, IL1B_PROT, TNF, TNF_PROT
C AGG: Aggrecan expression
C COLI: Collagen I expression
C COLII: Collagen II expression
C MMP3: Matrix metalloproteinase expression
C ADM4: A disintegrin and metalloproteinase with thrombospondin motifs expression
C IMN: Immuno negative
C ACC: Accumulation of the protein (depends on PN_TIME in hours)
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

C Old values of the local SDVs
      DOUBLE PRECISION pH_val_old, CELL_rho_old, CELL_viab_old

C Initialization of variables
C ---------------------------
      INCLUDE './src/IVD/transport/functions/UEL_init_var.f'

C
C--------------------------------------------------------------------------
C Definitions of the PROPS array
C--------------------------------------------------------------------------
C

C Mechanical properties from PROPS array
      E = PROPS(1)
      NU = PROPS(2)

C Extracting Water diffusivities from PROPS array. This represents the 
C D_WAT value from the equation. It's given for each solute (O2, lactate, glucose).
C Obtained from the .inp file
      D_WAT(1) = PROPS(3)
      D_WAT(2) = PROPS(4)
      D_WAT(3) = PROPS(5)

C Extracting the initial fraction of water content from PROPS array.
      ph_val_0 = PROPS(6)

C Extracting initial cell density from PROPS array. Represents ρ_cell,0 in the equation.
C Obtained from the .inp file
      CELL_rho_0 = PROPS(7)

C Extracting the initial cell viability from PROPS array.
      CELL_viab_0 = PROPS(8)

C a_K, b_K, c_K: Constants for the K(CELL_rho_0) model
      a_K = PROPS(9)
      b_K = PROPS(10)
      c_K = PROPS(11)

C ALPHA_pH: Death rate due to acidity
      ALPHA_pH = PROPS(12)

C GlucThres: Glucose threshold for cell viability
      GlucThres = PROPS(13)

C pHThres: pH threshold for cell viability
      pHThres = PROPS(14)

C NF0h: Alternative valuees for the SDVs - Homogeneous Initial water content
      NF0h = PROPS(15)

C FREQ: Frequency of load
      FREQ = PROPS(16)

C--------------------------------------------------------------------------
C INTEGER PROPS
C--------------------------------------------------------------------------
C
C set the number of the component of the IVD
C The number is CNUMBER = PROPS(1) and it will be stored in CMNAME
C AF = 1
C CEP = 2
C NP = 3
      CNUMBER = JPROPS(1)

      IF (CNUMBER == 0) THEN
            MATNAME = 'DECOY'

      ELSE IF (CNUMBER == 1) THEN
            MATNAME = 'AF'

      ELSE IF (CNUMBER == 2) THEN
            MATNAME = 'CEP'

      ELSE IF (CNUMBER == 3) THEN
            MATNAME = 'NP'

      END IF

C Number of local sdv's per integ point
      nlSdv  = JPROPS(2)

C Number of global sdv's per integ point
      ngSdv  = JPROPS(3)

C NF0h_d: Use NF0h? yes NF0h_d = 1, no NF0h_d = 0
      NF0h_d = JPROPS(4)

C nD_d: Use ∇D as ZERO (nD_d=0), Isotropic (nD_d=1), Anisotropic (nD_d=2)
      nD_d = JPROPS(5)

C PN_TIME: Total desired protein expresion simulation time (hours)
      PN_TIME = JPROPS(6)

C--------------------------------------------------------------------------
C Allocate memory for the globalSdv's used for viewing
C results on the dummy mesh
C--------------------------------------------------------------------------
C
C Uncomment to print the PROPS and IPROPS before allocate the SDVs
      ! IF(.NOT.allocated(globalSDVs)) THEN
      !       PRINT *, ''
      !       PRINT *, 'PROPS:'
      !       PRINT *, '------'
      !       PRINT *, 'E: ', E
      !       PRINT *, 'NU: ', NU
      !       PRINT *, 'D_WAT: ', D_WAT
      !       PRINT *, 'ph_val_0: ', ph_val_0
      !       PRINT *, 'CELL_rho_0: ', CELL_rho_0
      !       PRINT *, 'CELL_viab_0: ', CELL_viab_0
      !       PRINT *, 'a_K: ', a_K
      !       PRINT *, 'b_K: ', b_K
      !       PRINT *, 'c_K: ', c_K
      !       PRINT *, 'GlucThres: ', GlucThres
      !       PRINT *, 'NF0h: ', NF0h
      !       PRINT *, ' '
      !       PRINT *, 'JPROPS:'
      !       PRINT *, '-------'
      !       PRINT *, 'CNUMBER: ', CNUMBER
      !       PRINT *, 'nlSdv: ', nlSdv
      !       PRINT *, 'ngSdv: ', ngSdv
      !       PRINT *, 'NF0h_d: ', NF0h_d
      !       PRINT *, ''
      ! END IF


C Some internal functions to cosider:
C XIT: can be called from any Abaqus/Standard or Abaqus/Explicit user subroutine, respectively, to terminate an analysis

C
C*******************************************************************************
C Code
C*******************************************************************************
C

C Get the actual iteration number
      IF (LastStep /= KSTEP .OR. LastKINC /= KINC) THEN
          ITERATION = 0
          LastKINC = KINC
          LastStep = KSTEP
          LastElement = JELEM
      ELSE IF (LastElement == JELEM) THEN
          ITERATION = ITERATION + 1
      END IF
      IF (i_numElem .EQ. numElem) THEN  !start of the step
            i_numElem = 0
      END IF

C Identity tensor
      CALL onem(IdenM)

C Identity vector
      CALL onev(IdenV)

C DOFs: Degrees of freedom
      NDOFN = NDOFEL/NNODE

C--------------------------------------------------------------------------
C Nodal displacements, and solutes concentrations
C Obtaining the initial condition from the DOFs

      K = 0
      DO I = 1,NNODE
            DO J = 1,NDIM
C Nodal displacements
                  K = K + 1
                  U(I,J) = UALL(K)
                  dU(I,J) = dUALL(K,1)
                  U_old(I,J) = U(I,J) - dU(I,J)
C Obtain current nodal coodinates
                  COORDSC(J,I) = COORDS(J,I) + U(I,J)
            END DO

C Solutes concentrations (this also consider the BCs)
C
C Temperature (Dummy)
            K = K + 1

C O2
            K = K + 1
            O2_CONC(I) = UALL(K)
            dO2_CONC(I) = dUALL(K,1)
            O2_CONC_old(I) = O2_CONC(I) - dO2_CONC(I)
            dO2_CONCdt(I) = Vel(K)

C Lactate
            K = K + 1
            LACT_CONC(I) = UALL(K)
            dLACT_CONC(I) = dUALL(K,1)
            LACT_CONC_old(I) = LACT_CONC(I) - dLACT_CONC(I)
            dLACT_CONCdt(I) = Vel(K)

C Glucose
            K = K + 1
            GLUC_CONC(I) = UALL(K)
            dGLUC_CONC(I) = dUALL(K,1)
            GLUC_CONC_old(I) = GLUC_CONC(I) - dGLUC_CONC(I)
            dGLUC_CONCdt(I) = Vel(K)
      END DO

C If this is the initial step, read the predifined fields for O2, lactate, and glucose from the .inp file
C The variable is: PREDEF(2,NPREDF,NNODE)

C PREDEF(K1,1,K3)  Temperature.
C PREDEF(K1,2,K3)  First predefined field variable.
C PREDEF(K1,3,K3)  Second predefined field variable.
C Etc.	       Any other predefined field variable.
C PREDEF(K1,K2,K3) Total or incremental value of the K2th predefined field variable at the K3th node of the element.
C PREDEF(1,K2,K3)  Values of the variables at the end of the current increment.
C PREDEF(2,K2,K3)  Incremental values corresponding to the current time increment.

      IF (KSTEP .LT. 2) THEN
            DO I = 1,NNODE
C O2                  
                  O2_CONC_0(I) = PREDEF(1,2,I)

C Lactate
                  LACT_CONC_0(I) = PREDEF(1,3,I)

C Glucose
                  GLUC_CONC_0(I) = PREDEF(1,4,I)

            END DO
      END IF

C Uncomment to print the concentration values
C      IF (JELEM-ElemOffset .EQ. 4703) THEN  !start of the step
C            INCLUDE './src/IVD/transport/functions/infoSim.f'
C      END IF

C--------------------------------------------------------------------------
C Displacement increment, based on element diagonal:
C Impose any time-stepping changes on the increments of
C displacement if you want (based on element diagonal)
      Le = DSQRT(((COORDSC(1,1)-COORDSC(1,7))**TWO) + 
     +     ((COORDSC(2,1)-COORDSC(2,7))**TWO) +
     +     ((COORDSC(3,1)-COORDSC(3,7))**TWO))

C if dU(I,J is superior to 10 times Le, then decrease the time step 0.5 times
      DO I = 1,NNODE
            DO J = 1,NDIM
                  IF (DABS(dU(I,J)) .GT. 10.0D0*Le) THEN
                        PNEWDT = HALF
                  END IF
            END DO
      END DO

C
C--------------------------------------------------------------------------
C Take this opportunity to perform calculations at the element centroid. 
C F-bar method
C--------------------------------------------------------------------------
C
C Obtain shape functions and their local gradients at the element
C centriod, that means ξ=η=ζ=0.0, and nIntPt=1
C shapeFunc_C3D20(numInt,Cip,i_intPT,N,dN,NNODE,NDIM)

C      INCLUDE './src/IVD/transport/functions/f_bar_method.f'

C
C--------------------------------------------------------------------------
C Shape functions and derivatives
C Begin the loop over integration points
C--------------------------------------------------------------------------
C
C Obtain integration point local coordinates and weights
C 27-pt integration, numInt = 27
      CALL CintPtC3D20pt(Cip,whtG,numInt,NDIM)

C Loop over integration points
      JJ = 0 !JJ is used for tracking the state variables
      DO i_intPT = 1, numInt
      index = (i_numElem-1)*numInt + i_intPT ! index of the integration point (1-numElemInt)
C--------------------------------------------------------------------------

C Obtain shape functions and their local gradients
      IF(NNODE .EQ. 20) THEN
      CALL shapeFunc_C3D20(numInt,Cip,i_intPT,N,dNxi,d2Nxi,NNODE,NDIM)

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
      CALL mapShape3D(NNODE,dNxi,COORDS,dN,DETMapJ,stat)
      IF (stat .EQ. 0) THEN
            PNEWDT = HALF
            RETURN
      END IF

C Map shape functions from local to global current coordinate system
      CALL mapShape3D(NNODE,dNxi,COORDSC,dNC,DETMapJC,stat)
      IF (stat .EQ. 0) THEN
            PNEWDT = HALF
            RETURN
      END IF

C Calculate the second derivatives of the shape functions for the reference configuration
      CALL SecondDerN(NNODE,dNxi,d2Nxi,COORDS,d2N,stat)
      IF (stat .EQ. 0) THEN
            PNEWDT = HALF
            RETURN
      END IF

C--------------------------------------------------------------------------
C Calculating the variables

C Obtain the O2 concentration and its derivative's at
C this intPt at the begining and end of the increment
      O2_int = ZERO
      O2_int_old = ZERO
      dO2_dt = ZERO                              !temporal derivative of O2
      dO2dx = ZERO                               !spatial derivative of O2
      dO2dt = ZERO
      DO K = 1,NNODE
            O2_int = O2_int + O2_CONC(K)*N(K)
            O2_int_old = O2_int_old + O2_CONC_old(K)*N(K)
            dO2_dt = dO2_dt + dO2_CONCdt(K)*N(K)
            DO I = 1,NDIM
                  dO2dx(I,1) = dO2dx(I,1) + O2_CONC(K)*dNC(K,I)
            END DO
      END DO
      dO2dt = (O2_int - O2_int_old) / DTIME      !increment of O2

C Obtain the lactate concentration and its derivative's at
C this intPt at the begining and end of the increment
      LACT_int = ZERO
      LACT_int_old = ZERO
      dLACT_dt = ZERO                            !temporal derivative of lactate
      dLACTdx = ZERO                             !spatial derivative of lactate
      dLACTdt = ZERO
      DO K = 1,NNODE
            LACT_int = LACT_int + LACT_CONC(K)*N(K)
            LACT_int_old = LACT_int_old + LACT_CONC_old(K)*N(K)
            dLACT_dt = dLACT_dt + dLACT_CONCdt(K)*N(K)
            DO I = 1,NDIM
                  dLACTdx(I,1) = dLACTdx(I,1) + LACT_CONC(K)*dNC(K,I)
            END DO
      END DO
      dLACTdt = (LACT_int - LACT_int_old) / DTIME !increment of lactate

C Obtain the glucose concentration and its derivative's at
C this intPt at the begining and end of the increment
      GLUC_int = ZERO
      GLUC_int_old = ZERO
      dGLUC_dt = ZERO                            !temporal derivative of glucose
      dGLUCdx = ZERO                             !spatial derivative of glucose
      dGLUCdt = ZERO
      DO K = 1,NNODE
            GLUC_int = GLUC_int + GLUC_CONC(K)*N(K)
            GLUC_int_old = GLUC_int_old + GLUC_CONC_old(K)*N(K)
            dGLUC_dt = dGLUC_dt + dGLUC_CONCdt(K)*N(K)
            DO I = 1,NDIM
                  dGLUCdx(I,1) = dGLUCdx(I,1) + GLUC_CONC(K)*dNC(K,I)
            END DO
      END DO
      dGLUCdt = (GLUC_int - GLUC_int_old) / DTIME !increment of glucose

C--------------------------------------------------------------------------
C Obtain and modify the deformation gradient at this intPt
C   dN = ∂N/∂X  (shape-function gradients in the REFERENCE configuration)
C--------------------------------------------------------------------------
      DFG      = IdenM
      DFG_old  = IdenM
      DO I = 1, NDIM
            DO J = 1, NDIM
                  DO K = 1, NNODE
                  DFG(I,J)     = DFG(I,J)      + dN(K,J) * U(K,I)
                  DFG_old(I,J) = DFG_old(I,J)  + dN(K,J) * U_old(K,I)
                  END DO
            END DO
      END DO

C  --------------------------------------------
C    Determinant of the deformation gradient
C       
C  is determined from the boundary displcement
C  that for each increment correspond to those
C  of the "global" model (poromechanical one)
C
C  --------------------------------------------
C
C DFG: Deformation gradient
C DFG_inv: Inverse of the deformation gradient
C DFG_invT: Transpose of the inverse of the deformation gradient
C DETDFG: Determinant of the deformation gradient
C
      DFG_inv = ZERO
      DFG_invT = ZERO
      DETDFG = ZERO
      CALL matInvDet3D(DFG,DFG_inv,DFG_invT,DETDFG,stat)
      IF (stat .EQ. 0) THEN
      OPEN(15,FILE=infoFilePath_T_D, status='old', position='append')
      WRITE(15,*) ' '
      WRITE(15,*) '//////////////////////////////////////////////'
      WRITE(15,*) 'PROBLEM: DEF GRADIENT DET .LT. 0'
      WRITE(15,*) '//////////////////////////////////////////////'
      WRITE(15,*) ' '
      PRINT *, ' '
      PRINT *, '//////////////////////////////////////////////'
      PRINT *, 'PROBLEM: DEF GRADIENT DET .LT. 0'
      PRINT *, '//////////////////////////////////////////////'
      PRINT *, ' '
      CLOSE(15)
C Terminate the analysis
      CALL XIT
      END IF

C ------------------------------------------------------------
C Push-forward reference Hessian to get the mixed derivative
C d2Nmix(A,K,Ni) = Σ_M  F_inv(M,K) * d2N(A,M,Ni)
C ------------------------------------------------------------
      d2Nmix = ZERO
      DO A = 1, NNODE
            DO K = 1, NDIM                                              ! spatial index  x_k
                  DO Ni = 1, NDIM                                       ! reference index X_n
                        DO M = 1, NDIM                                  ! reference index X_m (summed)
                        d2Nmix(A,K,Ni) = d2Nmix(A,K,Ni) +
     1                                   DFG_inv(M,K) * d2N(A,M,Ni)
                        END DO
                  END DO
            END DO
      END DO

C ------------------------------------------------------------
C    Derivative of the deformation gradient
C    with respect to spatial coordinates
C    ∇ₓF  (3-rd-order tensor stored in dDFGdx(K,M,Ni))
C    dDFGdx(K,M,Ni) = ∂F(M,Ni) / ∂x_K
C    where F(M,Ni) = δ_MNi + Σ_A U(A,M) * ∂N^A/∂X_Ni
C  ⇒ ∂F(M,Ni)/∂x_K = Σ_A U(A,M) * ∂²N^A/(∂x_K ∂X_Ni)
C                   = Σ_A U(A,M) * d2Nmix(A,K,Ni)
C ------------------------------------------------------------
C 4) Then use d2Nmix(A,K,Ni) in your dDFGdx loop.
      dDFGdx = ZERO
      DO K = 1, NDIM                                                    ! gradient direction  (x, y, z)
            DO M = 1, NDIM                                              ! displacement component (row of F)
                  DO Ni = 1, NDIM                                       ! column of F
                        DO A = 1, NNODE
                              dDFGdx(K,M,Ni) = dDFGdx(K,M,Ni) +
     1                                         U(A,M) * d2Nmix(A,K,Ni)
                        END DO
                  END DO
            END DO
      END DO

C--------------------------------------------------------------------------
C SDVs from previous sim

C-------------------------------------------
C NF0: Initial water content
      NF0 = ZERO
C NF0h_d: Use NF0h? yes NF0h_d = 1, no NF0h_d = 0
C Passing the values of NF0 from node to integration point
C Using homogeneous value for NF0
      IF (NF0h_d .EQ. 1) THEN
            DO K = 1,NNODE
                  NF0 = NF0 + NF0h*N(K)
            END DO

C Using the values of NF0 from the previous mechanical simulation
      ELSE IF (NF0h_d .EQ. 0) THEN
            NF0 = SDV_NF0(JELEM-ElemOffset,i_intPT)
      END IF

C-------------------------------------------
C SP: Total Stress Pressure
      SP = ZERO
C Using the values of SP from the previous mechanical simulation
      SP = SDV_SP(JELEM-ElemOffset,i_intPT)

C--------------------------------------------------------------------------
C Obtain state variables from previous sep
      IF (KSTEP .LT. 2) THEN
            ph_val_old = ph_val_0
            CELL_rho_old = CELL_rho_0
            CELL_viab_old = CELL_viab_0

C This is not the first step, read old values
      ELSE
C            CELL_viab_old = localSDVs(JELEM,i_intPT,1)
            CELL_viab_old = SVARS(1 + JJ)

      ENDIF

C @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C Perform the time integration at this integ. point to compute
C all the specific forms and parameters needed for the solution

      CALL integ(STEPNAMES, NSTEPS, DTIME, KSTEP, NDIM,
     +           DFG_inv,DETDFG,dDFGdx,
     +           D_WAT, CELL_rho_0,
     +           ph_val_old, CELL_rho_old, CELL_viab_old,
     +           a_K, b_K, c_K, ALPHA_pH, GlucThres, pHThres, NF0,
     +           O2_int, LACT_int, GLUC_int,
     +           D_SOL, R_SOL, dR_SOL,
     +           dD_O2dx, dD_LACTdx, dD_GLUCdx, dD_SOLdx_MAG, nD_d,
     +           NF, pH_val, CELL_rho, CELL_viab, SP,
     +           PN_TIME, FREQ, PN_TIME_file, PN_FREQ_file,
     +           FREQ_activation, FREQUENCY_CURRENT,
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
     +           MATNAME)

C @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

C--------------------------------------------------------------------------
C Save the state variables at this integ point
C at the end of the increment

C      localSDVs(JELEM,i_intPT,1) = CELL_viab
      SVARS(1 + JJ) = CELL_viab
      JJ = JJ + nlSdv                      ! setup for the next intPt

C--------------------------------------------------------------------------
C Save the state variables at this integ point in the
C global array used for plotting field output
      globalSDVs(JELEM,i_intPT,UVAR_NF0) = NF0                          !1: Initial water content
      globalSDVs(JELEM,i_intPT,UVAR_NF) = NF                            !2: Water content

C Showing the initial concentrations
      IF (KSTEP .LT. 2) THEN
            O2_int_0 = ZERO
            LACT_int_0 = ZERO
            GLUC_int_0 = ZERO
            DO K = 1,NNODE
                  O2_int_0 = O2_int_0 + O2_CONC_0(K)*N(K)
                  LACT_int_0 = LACT_int_0 + LACT_CONC_0(K)*N(K)
                  GLUC_int_0 = GLUC_int_0 + GLUC_CONC_0(K)*N(K)
            END DO
            globalSDVs(JELEM,i_intPT,UVAR_O2) = O2_int_0                !3: O2 concentration
            globalSDVs(JELEM,i_intPT,UVAR_LACT) = LACT_int_0            !4: Lactate concentration
            globalSDVs(JELEM,i_intPT,UVAR_GLUC) = GLUC_int_0            !5: Glucose concentration

      ELSE
            globalSDVs(JELEM,i_intPT,UVAR_O2) = O2_int                  !3: O2 concentration
            globalSDVs(JELEM,i_intPT,UVAR_LACT) = LACT_int              !4: Lactate concentration
            globalSDVs(JELEM,i_intPT,UVAR_GLUC) = GLUC_int              !5: Glucose concentration
      END IF

      globalSDVs(JELEM,i_intPT,UVAR_pH_val) = pH_val                    !6: pH value
      globalSDVs(JELEM,i_intPT,UVAR_CELL_rho) = CELL_rho                !7: Cell density
      globalSDVs(JELEM,i_intPT,UVAR_CELL_viab) = CELL_viab              !8: Cell viability
      globalSDVs(JELEM,i_intPT,UVAR_NabD) = dD_SOLdx_MAG                !9: Gradient of the diffusion term
      globalSDVs(JELEM,i_intPT,UVAR_SP) = SP                            !10: Total Stress Pressure

C PN variables
C TNF & IL1B
      globalSDVs(JELEM,i_intPT,UVAR_IL1B) = IL1B                        !11: Interleukin 1 beta
      globalSDVs(JELEM,i_intPT,UVAR_IL1B_PROT) = IL1B_PROT              !12: Interleukin 1 beta protein
      globalSDVs(JELEM,i_intPT,UVAR_TNF) = TNF                          !13: Tumor necrosis factor
      globalSDVs(JELEM,i_intPT,UVAR_TNF_PROT) = TNF_PROT                !14: Tumor necrosis factor protein
C Aggrecan
      globalSDVs(JELEM,i_intPT,UVAR_AGG_IMN) = AGG_IMN_ACC              !15: Aggrecan expression
      globalSDVs(JELEM,i_intPT,UVAR_AGG_IL1B) = AGG_IL1B_ACC            !16: Aggrecan expression
      globalSDVs(JELEM,i_intPT,UVAR_AGG_TNF) = AGG_TNF_ACC              !17: Aggrecan expression
      globalSDVs(JELEM,i_intPT,UVAR_AGG_IL1B_TNF) = AGG_IL1B_TNF_ACC    !18: Aggrecan expression
C Collagen I
      globalSDVs(JELEM,i_intPT,UVAR_COLI_IMN) = COLI_IMN_ACC            !19: Collagen I expression
      globalSDVs(JELEM,i_intPT,UVAR_COLI_IL1B) = COLI_IL1B_ACC          !20: Collagen I expression
      globalSDVs(JELEM,i_intPT,UVAR_COLI_TNF) = COLI_TNF_ACC            !21: Collagen I expression
      globalSDVs(JELEM,i_intPT,UVAR_COLI_IL1B_TNF) =                    !22: Collagen I expression
     1           COLI_IL1B_TNF_ACC    
C Collagen II    
      globalSDVs(JELEM,i_intPT,UVAR_COLII_IMN) = COLII_IMN_ACC          !23: Collagen II expression
      globalSDVs(JELEM,i_intPT,UVAR_COLII_IL1B) = COLII_IL1B_ACC        !24: Collagen II expression
      globalSDVs(JELEM,i_intPT,UVAR_COLII_TNF) = COLII_TNF_ACC          !25: Collagen II expression
      globalSDVs(JELEM,i_intPT,UVAR_COLII_IL1B_TNF) =                   !26: Collagen II expression
     1           COLII_IL1B_TNF_ACC    
C MMP3    
      globalSDVs(JELEM,i_intPT,UVAR_MMP3_IMN) = MMP3_IMN_ACC            !27: Matrix metalloproteinase expression
      globalSDVs(JELEM,i_intPT,UVAR_MMP3_IL1B) = MMP3_IL1B_ACC          !28: Matrix metalloproteinase expression
      globalSDVs(JELEM,i_intPT,UVAR_MMP3_TNF) = MMP3_TNF_ACC            !29: Matrix metalloproteinase expression
      globalSDVs(JELEM,i_intPT,UVAR_MMP3_IL1B_TNF) =                    !30: Matrix metalloproteinase expression
     1           MMP3_IL1B_TNF_ACC    
C ADM4    
      globalSDVs(JELEM,i_intPT,UVAR_ADM4_IMN) = ADM4_IMN_ACC            !31: A disintegrin and metalloproteinase with thrombospondin motifs expression
      globalSDVs(JELEM,i_intPT,UVAR_ADM4_IL1B) = ADM4_IL1B_ACC          !32: A disintegrin and metalloproteinase with thrombospondin motifs expression
      globalSDVs(JELEM,i_intPT,UVAR_ADM4_TNF) = ADM4_TNF_ACC            !33: A disintegrin and metalloproteinase with thrombospondin motifs expression
      globalSDVs(JELEM,i_intPT,UVAR_ADM4_IL1B_TNF) =                    !34: A disintegrin and metalloproteinase with thrombospondin motifs expression
     1           ADM4_IL1B_TNF_ACC
C Frequency of load
      globalSDVs(JELEM,i_intPT,UVAR_FREQ) = FREQ                        !35: Frequency of load

C
C--------------------------------------------------------------------------
C
C                          RHS and Stiffness matrix
C
C--------------------------------------------------------------------------
C
C Fillin the RHS and AMATRX arrays

C Shape function vector (Nvec = 1xNNODE)
      Nvec = ZERO
      DO I = 1,NNODE
            Nvec(1,I) = N(I)
      END DO

C WHT : Gauss weight multiplied by jacobi of the transformation
      WHT = DETMapJC*whtG(i_intPT)

C D_SOL(3) : Diffusivity matrix for solutes within the medium
C R_SOL(3): Reaction matrix for solutes within the medium
C dR_SOL(3x3): Derivative of the reaction matrix for solutes within the medium
C dD_O2dx(3x1), dD_LACTdx(3x1), dD_GLUCdx(3x1): Gradient of the diffusion term

C dNC(NNODEx3): Gradient of the shape functions
C IdenV(3x1): Identity vector
C dO2dx(3x1), dLACTdx(3x1), dGLUCdx(3x1): Gradient of the solutes concentrations

C--------------------------------------------------------------------
C O2

C RHS concentration of O2
      RF_O2 = dO2_dt - R_SOL(1)

      R_O2 = R_O2 + WHT * 
     1     (
     2     TRANSPOSE(Nvec) * RF_O2                                      ! Time and reaction rate
     3     + MATMUL(dNC, D_SOL(1) * dO2dx)                              ! Isotropic Diffusion term contribution
     4     - TRANSPOSE(Nvec) *                                          ! Gradient of the diffusion term: dD_O2dx(1)
     5     DOT_PRODUCT(dD_O2dx(:,1), dO2dx(:,1))
     6     )

C Stiffness matrix of O2
      TF_O2 = -ONE/DTIME + dR_SOL(1,1)

C Diagonal term of the stiffness matrix      
      K_O2 = K_O2 + WHT *
     1     (
     2     TF_O2 * MATMUL(TRANSPOSE(Nvec), Nvec)                      ! Time derivative and reaction rate
     3     - D_SOL(1) * MATMUL(dNC, TRANSPOSE(dNC))                   ! derivative of the diffusion term contribution
     4     + MATMUL(TRANSPOSE(Nvec),                                  ! derivative of the Gradient of the diffusion term
     5     MATMUL(TRANSPOSE(dD_O2dx), TRANSPOSE(dNC))) 
     6     )

C Off-diagonal terms of the stiffness matrix:
C Stiffness matrix of O2 with Lactate
      K_O2LACT = K_O2LACT + WHT *
     1       (
     2       dR_SOL(1,2) * MATMUL(TRANSPOSE(Nvec), Nvec)                ! Derivative of the reaction term
     3       )

C Stiffness matrix of O2 with Glucose
      K_O2GLUC = K_O2GLUC + WHT *
     1       (
     2       dR_SOL(1,3) * MATMUL(TRANSPOSE(Nvec), Nvec)                ! Derivative of the reaction term
     3       )

C--------------------------------------------------------------------
C Lacate

C RHS concentration of lactate
      RF_LACT = dLACT_dt - R_SOL(2)

      R_LACT = R_LACT + WHT * 
     1     (
     2     RF_LACT * TRANSPOSE(Nvec)                                    ! Time and reaction rate
     3     + MATMUL(dNC, D_SOL(2) * dLACTdx)                            ! Isotropic Diffusion term contribution
     4     - TRANSPOSE(Nvec) *                                          ! Gradient of the diffusion term: dLACTdx(2)
     5     DOT_PRODUCT(dD_LACTdx(:,1), dLACTdx(:,1))
     6     )

C Stiffness matrix of lactate
      TF_LACT = -ONE/DTIME + dR_SOL(2,2)

C Diagonal term of the stiffness matrix
      K_LACT = K_LACT + WHT *
     1     (
     2     TF_LACT * MATMUL(TRANSPOSE(Nvec), Nvec)                      ! Time derivative and reaction rate
     3     - D_SOL(2) * MATMUL(dNC, TRANSPOSE(dNC))                     ! derivative of the diffusion term contribution
     4     + MATMUL(TRANSPOSE(Nvec),                                    ! derivative of the Gradient of the diffusion term
     5     MATMUL(TRANSPOSE(dD_LACTdx), TRANSPOSE(dNC))) 
     6     )


C Off-diagonal terms of the stiffness matrix:
C Stiffness matrix of Lactate with O2
      K_LACTO2 = K_LACTO2 + WHT *
     1       (
     2       dR_SOL(2,1) * MATMUL(TRANSPOSE(Nvec), Nvec)                ! Derivative of the reaction term
     3       )

C Stiffness matrix of Lactate with Glucose
      K_LACTGLUC = K_LACTGLUC + WHT *
     1       (
     2       dR_SOL(2,3) * MATMUL(TRANSPOSE(Nvec), Nvec)                ! Derivative of the reaction term
     3       )

C--------------------------------------------------------------------
C Glucose

C RHS concentration of glucose
      RF_GLUC = dGLUC_dt - R_SOL(3)

      R_GLUC = R_GLUC + WHT * 
     1     (
     2     RF_GLUC * TRANSPOSE(Nvec)                                    ! Time and reaction rate
     3     + MATMUL(dNC, D_SOL(3) * dGLUCdx)                            ! Isotropic Diffusion term contribution
     4     - TRANSPOSE(Nvec) *                                          ! Gradient of the diffusion term: dGLUCdx(3)
     5     DOT_PRODUCT(dD_GLUCdx(:,1), dGLUCdx(:,1))
     6     )

C Stiffness matrix of glucose
      TF_GLUC = -ONE/DTIME + dR_SOL(3,3)

C Diagonal term of the stiffness matrix
      K_GLUC = K_GLUC + WHT *
     1     (
     2     TF_GLUC * MATMUL(TRANSPOSE(Nvec), Nvec)                      ! Time derivative and reaction rate
     3     - D_SOL(3) * MATMUL(dNC, TRANSPOSE(dNC))                     ! derivative of the diffusion term contribution
     4     + MATMUL(TRANSPOSE(Nvec),                                    ! derivative of the Gradient of the diffusion term
     5     MATMUL(TRANSPOSE(dD_GLUCdx), TRANSPOSE(dNC))) 
     6     )

C Off-diagonal terms of the stiffness matrix:
C Stiffness matrix of Glucose with O2
      K_GLUCO2 = K_GLUCO2 + WHT *
     1       (
     2       dR_SOL(3,1) * MATMUL(TRANSPOSE(Nvec), Nvec)                ! Derivative of the reaction term
     3       )

C Stiffness matrix of Glucose with Lactate
      K_GLUCLACT = K_GLUCLACT + WHT *
     1       (
     2       dR_SOL(3,2) * MATMUL(TRANSPOSE(Nvec), Nvec)                ! Derivative of the reaction term
     3       )      

C End the loop over integration points
      END DO
C--------------------------------------------------------------------------

C Return Abaqus the RHS vector and the Stiffness matrix
      CALL AssembleElement(NDIM,NNODE,NDOFEL,NDOFN,
     +           RHS,NRHS,AMATRX,
     +           R_O2,K_O2,
     +           R_LACT,K_LACT,
     +           R_GLUC,K_GLUC,
     +           K_O2LACT,K_O2GLUC,
     +           K_LACTO2,K_LACTGLUC,
     +           K_GLUCO2,K_GLUCLACT)

C End return of RHS and AMATRX
C--------------------------------------------------------------------------

C
C--------------------------------------------------------------------------
C *********************** Exception handling ******************************
C--------------------------------------------------------------------------
C

C If the element number is greater than the total number of elements
C      INCLUDE './src/IVD/transport/functions/UEL_NumElem_Err.f'

C
C--------------------------------------------------------------------------
C ******************** Creating .txt for visualization ********************
C--------------------------------------------------------------------------
C

C Uncomment only to check what subroutine is doing
C it can generate large files
C
      ! IF (i_numElem .EQ. 1) THEN  !start of the step
      !       INCLUDE './src/IVD/transport/functions/infoSim.f'
      ! END IF

C
C--------------------------------------------------------------------------
C ************************* Finish the subroutine *************************
C--------------------------------------------------------------------------
C
      i_numElem = i_numElem + 1
      run_count_share = i_numElem

      RETURN
      END SUBROUTINE UEL
C--------------------------------------------------------------------------




C**********************************************************************************************************************************************
C ADDITIONAL SUBROUTINES
C**********************************************************************************************************************************************

C--------------------------------------------------------------------------
C Shape functions and derivatives
C--------------------------------------------------------------------------

      INCLUDE './src/IVD/transport/functions/UEL_C3D20_Functions.f'

C--------------------------------------------------------------------------
C Calculation subroutines
C--------------------------------------------------------------------------

      INCLUDE './src/IVD/transport/functions/UEL_calcFunc.f'

C --------------------------------------------------------------------------
C PN equation
C --------------------------------------------------------------------------

      INCLUDE './src/IVD/transport/functions/UEL_PNeq.f'
      INCLUDE './src/IVD/transport/functions/UEL_PNeq_Functions.f'