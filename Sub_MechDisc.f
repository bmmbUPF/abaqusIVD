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
C ******    Osmo-Poro-Hyper-Viscoelastic Model for the IVD              ******
C ******    Type: User MATerial (UMAT)                                  ******
C ******                                                                ******
C ******    Auth: Adrea Malandrino, 5/10/2012                           ******
C ******    Updates:                                                    ******
C ******    - Carlos Ruiz, 23/07/2024                                   ******
C ******      Added strain-dependent permeability calculation           ******
C ******    - Estefano Muñoz-Moya, 05/08/2024                           ******
C ******      Adaptation to be compatible with mechano-transport model  ******
C ******                                                                ******
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
C MODULES
C -------------------------

C**********************************************************************************************************************************************
C Module to define the global variables of the simulation
C**********************************************************************************************************************************************

C
      MODULE GLOBAL

C  This module is used to define the global variables of the simulation

C  numElem
C  Total number of elements in the rmesh

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
      INTEGER :: LastInt = -1
      DOUBLE PRECISION :: LastDTIME = -1.D0

C -----------------------------------------------------------
C Number of dimensions
      INTEGER, PARAMETER :: NDIM = 3

C -----------------------------------------------------------
C Set the number of elements used here
C IVD model: 19392 elements
      INTEGER, PARAMETER :: numElem = 19392

C -----------------------------------------------------------
C Set the larget id number for the elements
C IVD model: 52944 largest id element
      INTEGER, PARAMETER :: numIdElem = 52944

C -----------------------------------------------------------
C Set the number of UMAT integration points
      INTEGER, PARAMETER :: numInt = 27

C -----------------------------------------------------------
C Set number of elements * number of integration points (previously called sou)
      INTEGER, PARAMETER :: numElemInt = numElem*numInt

C -----------------------------------------------------------
C Defining the STATEV numbers
c------------------------------------------------------------------------
c    STATEV(1)=PERMEABILITY (mm^4/Ns)
c    STATEV(2)=INITIAL WATER CONTENT (% wet weight)
c    STATEV(3)=INITIAL FIXED CHARGE DENSITY [meq/ml]
c    STATEV(4)=WATER CONTENT (% wet weight)
c    STATEV(5)=FIXED CHARGE DENSITY [meq/ml]
c    STATEV(6)=EXTRAFIBRILLAR FIXED CHARGE DENSITY [meq/ml]
c    STATEV(7)=EXTRAFIBRILLAR WATER CONTENT (% wet weight)
c    STATEV(8)=SWELLING PRESSURE [MPa]
      INTEGER, PARAMETER :: STATEV_K = 1
      INTEGER, PARAMETER :: STATEV_NF0 = 2
      INTEGER, PARAMETER :: STATEV_CF0 = 3
      INTEGER, PARAMETER :: STATEV_NF = 4
      INTEGER, PARAMETER :: STATEV_CF = 5
      INTEGER, PARAMETER :: STATEV_CFext = 6
      INTEGER, PARAMETER :: STATEV_NFext = 7
      INTEGER, PARAMETER :: STATEV_DELTAPI = 8

C -----------------------------------------------------------
C Defining the UVARM numbers
c------------------------------------------------------------------------
c    UVAR(1)=HYDROSTATIC PRESSURE [MPa]
      INTEGER, PARAMETER :: UVAR_SP = 1                                 !Hydrostatic pressure
      INTEGER, PARAMETER :: UVAR_test = 2                               !Test variable

C -----------------------------------------------------------
C Defining vector with STATEV numbers and names
      CHARACTER*256 :: STATEV_NAMES(8)

C -----------------------------------------------------------
C Defining vector with UVARM numbers and names
      CHARACTER*256 :: UVARM_NAMES(2)

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
C STATEV NAMES subroutine
C**********************************************************************************************************************************************

C This subroutine is used to initialize the STATEV_NAMES array
C
      SUBROUTINE Initialize_STATEV_NAMES()
C
C ********************************************************************
C
      USE GLOBAL
C

C Initialize the STATEV_NAMES array
      STATEV_NAMES(STATEV_K) = 'Permeability'
      STATEV_NAMES(STATEV_NF0) = 'Initial Water Content'
      STATEV_NAMES(STATEV_CF0) = 'Initial Fixed Charge Density'
      STATEV_NAMES(STATEV_NF) = 'Water Content'
      STATEV_NAMES(STATEV_CF) = 'Fixed Charge Density'
      STATEV_NAMES(STATEV_CFext) = 'Extrafibrillar Fixed Charge Density'
      STATEV_NAMES(STATEV_NFext) = 'Extrafibrillar Water Content'
      STATEV_NAMES(STATEV_DELTAPI) = 'Swelling Pressure'

C--------------------------------------------------------------------------
      RETURN
      END SUBROUTINE Initialize_STATEV_NAMES

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
      UVARM_NAMES(UVAR_SP) = 'Total Stress Pressure'
      UVARM_NAMES(UVAR_test) = 'User can define more variables'

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
C
      DOUBLE PRECISION :: SP, S11, S22, S33, POR_temp

C Initialize variables
      SP = ZERO
      S11 = ZERO
      S22 = ZERO
      S33 = ZERO
      POR = ZERO

C Get the variable POR from the UEL
      CALL GETVRM('POR',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,
     1      MATLAYO,LACCFLA)
      POR_temp = ARRAY(1)

C Get Total Stress Pressure
      CALL GETVRM('S',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,
     1      MATLAYO,LACCFLA)
      S11 = ARRAY(1)
      S22 = ARRAY(2)
      S33 = ARRAY(3)
      SP = -(S11 + S22 + S33 + 3*POR_temp)/THREE

C Transfer the SDV's from the UEL to the UVARM
      UVAR(UVAR_SP) = SP                                              ! Hydrostatic pressure

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
      USE GLOBAL
C
      INCLUDE 'ABA_PARAM.INC'
C
C Define the COMMON block for the info file
      COMMON /comm_info/ infoFilePath_T_D,cwd,outdir
c
      CHARACTER*256 cwd,outdir,infoFile_T_D,infoFilePath_T_D,JOBNAME
      CHARACTER*256 sdvNF0File,sdvNF0FilePath
      DIMENSION TIME(2)
      INTEGER error, LENJOBNAME
C For SDVs
      INTEGER iElem,iIP
      DOUBLE PRECISION SDV
      integer :: STATUS = 0

C
C--------------------------------------------------------------------
C Definition of paths:
C--------------------------------------------------------------------
C

C write a file with the information of the simulation
      IF (LOP .EQ. 0) THEN            !start of the analysis

C Call DATE_AND_TIME
      CALL DATE_AND_TIME(values=date_ini)

C Get the job name
      CALL GETJOBNAME(JOBNAME,LENJOBNAME)
      CALL GETCWD(cwd)
      CALL GETOUTDIR(outdir,lenoutdir)
      lenoutdir = 50
      infoFile_T_D = '/' // TRIM(JOBNAME) // '.info'
      infoFilePath_T_D = TRIM(outdir) // TRIM(infoFile_T_D)

C--------------------------------------------------------------------------
C Create the info file
      OPEN(15,FILE=infoFilePath_T_D, STATUS='UNKNOWN')
      WRITE(15,*) ' '
      CLOSE(15)

      END IF
C--------------------------------------------------------------------------

C Write and print the simulation info file
      INCLUDE './src/IVD/mechanic/functions/write_info_Mechanic.f'

C--------------------------------------------------------------------------
      RETURN
      END SUBROUTINE UEXTERNALDB

C**********************************************************************************************************************************************
C USDFLD subroutine input and workspace definitions
C**********************************************************************************************************************************************

      SUBROUTINE USDFLD(FIELD,STATEV,PNEWDT,DIRECT,T,CELENT,
     1 TIME,DTIME,CMNAME,ORNAME,NFIELD,NSTATEV,NOEL,NPT,LAYER,
     2 KSPT,KSTEP,KINC,NDI,NSHR,COORD,JMAC,JMATYP,MATLAYO,
     3 LACCFLA)
C
C ********************************************************************
C
      USE GLOBAL
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*256 CMNAME,ORNAME
      CHARACTER*3  FLGRAY(15)
      DIMENSION FIELD(NFIELD),STATEV(NSTATEV),DIRECT(3,3),
     1 T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),
     1 COORD(*)
      DOUBLE PRECISION :: VOID, ALPHA, E0, M, POR, Next

C
C------------------------------------------------------------------------
C                      PREDEFINED FIELD (FV) VARIABLES 
C
C             Fixed paremeters for permeability calculation 
C    (from Schoerder et al. (2008) Journal of Orthopaedic Research)
C
C  - alpha: Initial permeability (Variable=1 in inp file, ARRAY(1))
C  - M: positive constant (Variable=2 in inp file, ARRAY(2))
C
C          Values taken from alpha.inp and emme.inp files
C
C              Last update by Carlos Ruiz 5/20/2013
C
C------------------------------------------------------------------------
C
      CALL GETVRM('FV',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,
     1      MATLAYO,LACCFLA)
c
      IF (KSTEP .EQ. 1 .AND. KINC .EQ. 1) THEN
            ALPHA = ARRAY(1)
      ELSE
            ALPHA = ARRAY(4)
      END IF
c
      M = ARRAY(2)
c
c------------------------------------------------------------------------
c
c                 EXTRAFIBRILLAR WATER CONTENT (Next)
c    Next comes from the UMAT_MATRIX routine as state variavle (SDV) 
c              Next is ARRAY(7) for SDV (STATEV(7) in the UMAT_MATRIX)
c
c------------------------------------------------------------------------
c
      CALL GETVRM('SDV',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,
     1      MATLAYO,LACCFLA)
c
      Next=ARRAY(7)
C
c------------------------------------------------------------------------
c
c         EXTRAFIBRILLAR WATER CONTENT for first step, increment 1
c
c------------------------------------------------------------------------
c
      IF (KSTEP .EQ. 1 .AND. KINC .EQ. 1) THEN
            IF (CMNAME(1:2).EQ.'NP') THEN
                  Next=0.75D0
            ELSEIF (CMNAME(1:2).EQ.'AF') THEN
                  Next=0.70D0
            ELSEIF (CMNAME(1:2).EQ.'CE') THEN
                  Next=0.61D0
            END IF
      END IF
c
c------------------------------------------------------------------------
c
c              Strain-dependent permeability calculation
c From Eq. 4 in Schroeder et al. (2008) Journal of Orthopaedic Research
c
c------------------------------------------------------------------------
c
      FIELD(1) = ALPHA/((1-Next)**(M))
c
c------------------------------------------------------------------------
c
c    Store the permeability (k) as a solution dependent state variable
c
c------------------------------------------------------------------------
c
      STATEV(1) = FIELD(1)

C--------------------------------------------------------------------------
      RETURN
      END SUBROUTINE USDFLD

c
c************************************************************************
c
c************************************************************************
c              M A T R I X    M A T E R I A L   M O D E L 
c 
c         From Wilson et al. (2006) OsteoArthritis and Cartilage
c          Adapted to the IVD by Schoerder et al. (2006-2010)
c
c               Last update by Carlos Ruiz, 23/04/2024 
c************************************************************************
c
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,
     2 PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,
     3 DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
C ********************************************************************
C
      USE GLOBAL
C
      INCLUDE "ABA_PARAM.INC"
C
      DIMENSION STRESS(NTENS),STATEV(NSTATEV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4),JAY(NSTATEV),DFGR(3,3)

      CHARACTER*256 CMNAME, cwd, outdir, infoFile_T_D, infoFilePath_T_D
C

      IF (CMNAME(1:4) .EQ. 'COLL') THEN
            CALL UMAT_COLLAGEN(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,
     2 PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,
     3 DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
      ELSE
            CALL UMAT_MATRIX(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,
     2 PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,
     3 DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
      END IF

C--------------------------------------------------------------------------
      RETURN
      END SUBROUTINE UMAT
c
c************************************************************************
c
c************************************************************************
c            M A T R I X    M A T E R I A L   M O D E L 
c 
c         From Wilson et al. (2006) OsteoArthritis and Cartilage
c          Adapted to the IVD by Schoerder et al. (2006-2010)
c
c               Last update by Carlos Ruiz, 23/04/2024 
c************************************************************************
c
      SUBROUTINE UMAT_MATRIX(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,
     2 PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,
     3 DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
C ********************************************************************
C
      USE GLOBAL
C
      INCLUDE "ABA_PARAM.INC"
C
      DIMENSION STRESS(NTENS),STATEV(NSTATEV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JAY(NSTATEV),DFGR(3,3)

      CHARACTER*256 CMNAME, cwd, outdir, infoFile_T_D, infoFilePath_T_D
c
c    LOCAL ARRAYS
c------------------------------------------------------------------------
c    BBAR   - DEVIATORIC RIGHT CAUCHY-GREEN TENSOR
c    DISTGR - DEVIATORIC DEFORMATION GRADIENT (DISTORTION TENSOR)
c------------------------------------------------------------------------
c
      DOUBLE PRECISION :: fi_int,fi_ext,cF0,cF,NF0,ctot,RR,TT,
     1 Posm,Posm0,dPosm,rhoS,g_ext,dcFexfdPosm,dcFexfdnnf,dcFexfdcF,
     2 Z,dS,dfi_intdcFexf, BBAR(6),DISTGR(3,3),x(4),
     3 dfi_PMdX,dfi_int,dctotdcFexf,dg_intdcFexf,A,B,C,D,E,
     4 dg_MMdcFexf,dg_MMdg_int,dg_MMdctot,dg_PMdX,dg_int,
     5 dctotdg_int,dctot,dXdctot,dXdcFexf,dnnf,
     6 dcF,NF,COLLF_m,ci_phi,NEXTF_m,NEXTF,cFexf,gg0,gg,XX,
     7 g_PM,g_MM,g_int,fi_PM,fi_MM,VOLU1,VOLU2,NS0
c
C PROPS
      DOUBLE PRECISION :: C10,NF_m,ce,cFeq,COLL_rho
C STIFFNESS
      DOUBLE PRECISION :: DFG_inv(NDIM,NDIM),DFG_invT(NDIM,NDIM),DETDFG
      DOUBLE PRECISION :: SCALE,TRBBAR,EG,EG23,PR,EK,G23
      DOUBLE PRECISION :: TERM1,TERM2,TERM3
C LOOP VARIABLES
      INTEGER :: i,j,k,K1,K2,stat

c
c------------------------------------------------------------------------
c    UMAT FOR COMPRESSIBLE NEO-HOOKEAN HYPERELASTICITY
c    CANNOT BE USED FOR PLANE STRESS
c------------------------------------------------------------------------
c
c            Extra-cellular matrix material properties
c   Taken from annulus.inp, nucleus.inp, and endplate.inp files
c
c------------------------------------------------------------------------
c    PROPS(1) - SHEAR MODULUS [MPa]
c    PROPS(2) - INITIAL WATER CONTENT (% wet weight)
c    PROPS(3) - EXTERNAL SALT CONCENTRATION [meq/ml]
c    PROPS(4) - INITIAL FIXED CHARGE DENSITY [meq/ml]
c    PROPS(5) - COLLAGEN CONTENT (% dry weight)
c------------------------------------------------------------------------
c
c   Output state variables (SDV) for visualization in the odb file
c
c------------------------------------------------------------------------
c    STATEV(1)=PERMEABILITY (mm^4/Ns)
c    STATEV(2)=INITIAL WATER CONTENT (% wet weight)
c    STATEV(3)=INITIAL FIXED CHARGE DENSITY [meq/ml]
c    STATEV(4)=WATER CONTENT (% wet weight)
c    STATEV(5)=FIXED CHARGE DENSITY [meq/ml]
c    STATEV(6)=EXTRAFIBRILLAR FIXED CHARGE DENSITY [meq/ml]
c    STATEV(7)=EXTRAFIBRILLAR WATER CONTENT (% wet weight)
c    STATEV(8)=SWELLING PRESSURE [MPa]
c------------------------------------------------------------------------
c

C
C Get the actual iteration number
      IF (LastStep /= KSTEP .OR. LastKINC /= KINC .OR.
     1    LastDTIME /= DTIME) THEN
            ITERATION = 0
            LastKINC = KINC
            LastStep = KSTEP
            LastDTIME = DTIME
            LastElement = NOEL
            LastInt = NPT
      ELSE IF (LastElement == NOEL) THEN
            ITERATION = ITERATION + 1
      END IF

C Initialization of variables
C ---------------------------
      INCLUDE './src/IVD/mechanic/functions/UMATMATRIX_init_var.f'
C
C--------------------------------------------------------------------------
C Definitions of the PROPS array
C--------------------------------------------------------------------------
C
      C10 = PROPS(1)/TWO                                                !Neo hookean constant related to the shear modulus (Gm)
      ce = PROPS(3)                                                     !Initial external concentration of salt 
      COLL_rho = PROPS(5)                                               !Collagen content 
c 
c------------------------------------------------------------------------      
c
c    JACOBIAN AND DISTORTION TENSOR
c
c------------------------------------------------------------------------
c
C Determinant of the deformation gradient
      DETDFG = ZERO
      CALL matInvDet3D(DFGRD1,DFG_inv,DFG_invT,DETDFG,stat)
C 
      SCALE=DETDFG**(-ONE/THREE)

      DO K1=1, 3
            DO K2=1, 3
                  DISTGR(K2, K1)=SCALE*DFGRD1(K2, K1)
            END DO
      END DO
c
c------------------------------------------------------------------------
c
c	          PARAMETER INITIAL CONDITION 
c              Last update by Carlos Ruiz 5/20/2013
c
c------------------------------------------------------------------------
c
      IF (KSTEP.EQ.1) THEN
            rhoS=1.4338D0  ! Mass density of solid matrix [g/m]) (Shapiro et al., Osteoarthritis and Cartilage, 9:533-539 (2001); Basser et al., Arch Biochem Biophys, 351:207-219 (1998))
c
c------------------------------------------------------------------------
c
c     Determine fluid fraction that is required in equilibrium
c                        For Nucleus
c
c------------------------------------------------------------------------
            IF (CMNAME(1:2).EQ.'NP') THEN
                  NF_m = PROPS(2)
                  NF = (NF_m*rhoS)/(-NF_m+rhoS*NF_m+ONE)
                  cFeq = PROPS(4)
c
c------------------------------------------------------------------------
c
c     Determine fluid fraction that is required in equilibrium 
c            For Annulus and Catilage endplates
c
c------------------------------------------------------------------------
c
            ELSE
                  NF_m = PROPS(2)
                  cFeq = PROPS(4)
                  NF = (NF_m*rhoS)/(-NF_m+rhoS*NF_m+ONE)
            END IF
c
c------------------------------------------------------------------------
c
c             Determine initial fluid fraction 
c            based on equilibrium fluid fraction
c   Eq. 16 in Wilson et al. (2005) Journal of Biomechanics
c
c------------------------------------------------------------------------
c
            NF0 = NF*DETDFG-DETDFG+ONE
            STATEV(2) = NF0                                             ! Initial fluid volume fraction storage
            NS0 = ONE-NF0
c
c------------------------------------------------------------------------
c
c           Determine initial fixed charge density
c          based on equilibrium fixed charge density
c   Eq. 17 in Wilson et al. (2005) Journal of Biomechanics
c
c------------------------------------------------------------------------
c
            cF0 = cFeq*(NF0-ONE+DETDFG)/NF0
            STATEV(3) = cF0                                             ! Initial fixed charge density storage
      END IF
c
c------------------------------------------------------------------------
c
c     Definition of the initial water content (NF0) and 
c            initial fixed charged density (cF0) 
c
c------------------------------------------------------------------------
c
      NF0 = STATEV(2)                                                   ! Initial water content (volume fraction)
      cF0 = STATEV(3)                                                   ! Initial fixed charge density (mEq/mL)
c
c------------------------------------------------------------------------
c
c                     OSMOTIC PRESSURE MODEL
c     From Schroeder et al. (2008) Journal of Orthopaedic Research
c                 Adapted by Andrea Malandrino
c
c------------------------------------------------------------------------
c
      RR = 8.3145D0                                                     !Gas constant [N mm /(mmol K)]
      TT = 293.0D0                                                      !Absolute temperature [K]
      rhoS = 1.4338D0                                                   !Mass density of solid matrix [g/m]) (Shapiro et al., Osteoarthritis and Cartilage, 9:533-539 (2001); Basser et al., Arch Biochem Biophys, 351:207-219 (1998))
      x(1) = 0.3507D0                                                   !Material prop for function of activity coeff. [-]
      x(2) = 0.5060D0                                                   !Material prop for function of activity coeff. [-]
      x(3) = 3.1211D0                                                   !Material prop for function of activity coeff. [-]
      x(4) = 0.6447D0                                                   !Material prop for function of activity coeff. [-]
      fi_ext = 0.93D0                                                   !External osmotic coefficient [-]
      g_ext = x(1)*EXP(-cF0**x(2))**x(3)+x(4)                           !External activity coefficient [-]
c
c
c
c------------------------------------------------------------------------
c
c         Determine current fluid fraction (Volume based) (NF)
c  From Eq. 5 in Schroeder et al. (2007) Journal of Orthopaedic Research
c
c------------------------------------------------------------------------
c
      NF=(NF0-ONE+DETDFG)/DETDFG
      STATEV(4) = NF
c
c------------------------------------------------------------------------
c
c     Determine current fluid fraction (Mass based) (NF_m)
c From Eq. 28 in Wilson et al. (2007) Biomechanics and Modelling in Mechanobiology
c
c------------------------------------------------------------------------
c
      NF_m = NF/(NF+(ONE-NF)*rhoS)
c
c------------------------------------------------------------------------
c
c     Define collagen mass fraction t.o.v. total mass (COLLF_m)
c From Eq. 24 in Wilson et al. (2007) Biomechanics and Modelling in Mechanobiology
c
c------------------------------------------------------------------------
c
      COLLF_m = COLL_rho*(ONE-NF_m)
c
c------------------------------------------------------------------------
c
c            Determine current fixed charge density (cF)
c  From Eq. 5 in Schroeder et al. (2007) Journal of Orthopaedic Research
c
c------------------------------------------------------------------------
c
      cF=cF0*NF0/(NF0-ONE+DETDFG)
      STATEV(5) = cF
c
c------------------------------------------------------------------------
c
c       Initial guess for the osmotic pressure gradient (Posm)
c  From Eq. 4 in Schroeder et al. (2007) Journal of Orthopaedic Research
c
c               GG is given by g_ext/g_int
c          g_ext: external activity coefficient
c          g_int: internal activity coefficient 
c
c------------------------------------------------------------------------
c      
      Posm=-RR*TT*(SQRT(cF**TWO+FOUR*ce**TWO)-TWO*ce)*1.e-3
      gg=ONE 
c
c------------------------------------------------------------------------
c
c       Iteratively determine the osmotic pressure gradient (Posm)
c
c------------------------------------------------------------------------
c
      DO i=1,50
            Posm0=Posm
c
c------------------------------------------------------------------------
c
c    Determine the mass of intra-fibrillar water per collagen mass (ci_phi)
c From Eq. 11 and 12 in Schroeder et al. (2007) Journal of Orthopaedic Research
c
c------------------------------------------------------------------------
c	 
            IF (CMNAME(1:2).EQ.'NP') THEN
                  ci_phi = 0.4794D0*EXP(3.218D0*Posm)+0.7410D0 ! Eq. 12
            ELSE
                  ci_phi = 0.5931D0*EXP(3.396D0*Posm)+0.9222D0 ! Eq. 11
            END IF
c
c------------------------------------------------------------------------
c
c    Determine the extrafibrillar water mass fraction (NEXTF_m)
c  From Eq. 7 in Schroeder et al. (2007) Journal of Orthopaedic Research
c
c------------------------------------------------------------------------
c	 
            NEXTF_m = NF_m-ci_phi*COLLF_m
c
c------------------------------------------------------------------------
c
c    Determine the extrafibrillar water fraction (NEXTF)
c From Eq. 22 in Wilson et al. (2007) Biomechanics and Modelling in Mechanobiology
c
c------------------------------------------------------------------------
c
            IF (NEXTF_m .LT. 0.2D0) THEN
                  NEXTF_m = 0.2D0
            END IF
            NEXTF = NEXTF_m*(NF+(ONE-NF)*rhoS) 
c
c------------------------------------------------------------------------
c
c    Determine the extrafibrillar fixed charge density (cexf)
c  From Eq. 6 in Schroeder et al. (2007) Journal of Orthopaedic Research
c
c------------------------------------------------------------------------
c
            cFexf = cF*NF/NEXTF               
            STATEV(6) = cFexf
c
            DO j=1,50
                  gg0 = gg
                  ctot = SQRT(cFexf**TWO+FOUR*gg**TWO*ce**TWO)          ! Total ion concentration [M]
                  XX = cFexf/(cFexf+ctot)                               ! Ratio of fixed charge density to free electrolyte
                  g_PM = EXP(-HALF*0.99D0*XX)                           ! Osmotic coefficient of poly-on/mobile ion interaction in absence of salt
                  g_MM = x(1)*EXP(-(HALF*ctot)**x(2))**x(3)+x(4)        ! Osmotic coefficient of mobile ion/mobile ion interaction
                  g_int = g_MM*g_PM
                  gg = g_ext/g_int
                  IF (ABS(ONE-gg0/gg).LE.1.e-4) THEN
                        GO TO 1071
                  END IF
            END DO
1071        fi_PM = ONE-HALF*0.99D0*XX
            fi_MM = 0.93D0
            fi_int = fi_MM*fi_PM
            Posm = -RR*TT*(fi_int*ctot-TWO*fi_ext*ce)*1.e-3  !Eq.53
            IF (ABS(ONE-Posm0/Posm) .LE. 1.e-4) THEN
                  GO TO 1072
            END IF
      END DO
c
c------------------------------------------------------------------------
c
c     Save extrafibrillar water fraction for strain-dependent permeability
c
c------------------------------------------------------------------------
c
1072  STATEV(7) = NEXTF
      STATEV(8) = Posm

c
c------------------------------------------------------------------------
c      
c     Determine the derivative of cF and NF to JJ
c
c------------------------------------------------------------------------
c
      dcF = -cF0*NF0/(NF0-ONE+DETDFG)**TWO
      dnnf = ONE/DETDFG-(NF0+DETDFG-ONE)/DETDFG**TWO
      dXdcFexf = ctot/(cFexf+ctot)**TWO
      dXdctot = -cFexf/(cFexf+ctot)**TWO

c     dctot=dctotdcFexf*dcFexf+dctotdg_int*dg_int
      dctotdcFexf = cFexf/ctot
      dctotdg_int = -FOUR*g_ext**TWO*ce**TWO/(ctot*g_int**THREE)

c     dg_int=g_MM*dg_PM+g_PM*dg_MM
      dg_PMdX = (-HALF*0.99D0)*EXP(-HALF*0.99D0*XX)
      dg_MMdctot = -x(1)*EXP(-(HALF*ctot)**x(2))**x(3)*x(3)*
     1             (HALF*ctot)**x(2)*x(2)/ctot
      dg_MMdg_int = dg_MMdctot*dctotdg_int
      dg_MMdcFexf = dg_MMdctot*dctotdcFexf

      A = g_PM*dg_MMdg_int
      B = g_PM*dg_MMdcFexf
      C = g_MM*dg_PMdX*dXdcFexf
      D = g_MM*dg_PMdX*dXdctot*dctotdcFexf
      E = g_MM*dg_PMdX*dXdctot*dctotdg_int
      dg_intdcFexf = (B+C+D)/(ONE-(A+E))
      dctotdcFexf = dctotdcFexf+dctotdg_int*dg_intdcFexf

c     dfi_int=fi_MM*dfi_PM
      dfi_PMdX = -HALF*0.99D0
      dfi_intdcFexf = fi_MM*dfi_PMdX*dXdcFexf+
     1                fi_MM*dfi_PMdX*dXdctot*dctotdcFexf

c     dS=Z*dcFexf
      Z = -1.e-3*RR*TT*ctot*dfi_intdcFexf-1.e-3*RR*TT*fi_int*dctotdcFexf
c
c------------------------------------------------------------------------
c
c     Determine the derivative of cFexf
c
c------------------------------------------------------------------------
c
      IF (CMNAME(1:2).EQ.'NP') THEN
            dcFexfdcF = NF/(NF/(NF+(ONE-NF)*rhoS)-ci_phi*COLLF_m)/
     1                  (NF+(ONE-NF)*rhoS)
            dcFexfdnnf = -rhoS*COLL_rho*ci_phi*cF/(NF-ci_phi*
     1                   COLL_rho*rhoS+ci_phi*COLL_rho*rhoS*NF)**TWO
            dcFexfdPosm = 1.5427D0*cF*NF/(NF_m-(0.4794D0*EXP(3.218D0*
     1                    Posm)+0.7410D0)*COLLF_m)**TWO/
     2                    (NF+(ONE-NF)*rhoS)*EXP(3.218D0*Posm)*COLLF_m
            ELSE
            dcFexfdcF = NF/(NF/(NF+(ONE-NF)*rhoS)-ci_phi*COLLF_m)/
     1                  (NF+(ONE-NF)*rhoS)
            dcFexfdnnf = -rhoS*COLL_rho*ci_phi*cF/(NF-ci_phi*
     1                   COLL_rho*rhoS+ ci_phi*COLL_rho*rhoS*NF)**TWO
            dcFexfdPosm = 2.1920D0*cF*NF/(NF_m-(0.5931D0*EXP(3.696D0*
     1                    Posm)+0.9222D0)*COLLF_m)**TWO
     2                    /(NF+(ONE-NF)*rhoS)*EXP(3.696D0*Posm)*COLLF_m
      END IF
c
c------------------------------------------------------------------------
c
c     Determine derivative of osmotic pressure to JJ (Posm)
c
c------------------------------------------------------------------------
c
      dPosm = Z*(dcFexfdcF*dcF+dcFexfdnnf*dnnf)/(ONE-Z*dcFexfdPosm)
c
c------------------------------------------------------------------------
c
c------------------------------------------------------------------------
c
c    CALCULATE LEFT CAUCHY-GREEN TENSOR (Voigt order: 11 22 33 12 13 23)
c
c------------------------------------------------------------------------
c
      BBAR(1)=DISTGR(1, 1)**TWO+DISTGR(1, 2)**2+DISTGR(1, 3)**TWO
      BBAR(2)=DISTGR(2, 1)**TWO+DISTGR(2, 2)**TWO+DISTGR(2, 3)**TWO
      BBAR(3)=DISTGR(3, 1)**TWO+DISTGR(3, 2)**TWO+DISTGR(3, 3)**TWO
      BBAR(4)=DISTGR(1, 1)*DISTGR(2, 1)+DISTGR(1, 2)*DISTGR(2, 2)
     1       +DISTGR(1, 3)*DISTGR(2, 3)
      IF(NSHR.EQ.3) THEN
        BBAR(5)=DISTGR(1, 1)*DISTGR(3, 1)+DISTGR(1, 2)*DISTGR(3, 2)
     1         +DISTGR(1, 3)*DISTGR(3, 3)
        BBAR(6)=DISTGR(2, 1)*DISTGR(3, 1)+DISTGR(2, 2)*DISTGR(3, 2)
     1         +DISTGR(2, 3)*DISTGR(3, 3)
      END IF     
c
c------------------------------------------------------------------------
c
c    CALCULATE THE STRESS AND ALGORITHMIC BULK MODULUS
c
c------------------------------------------------------------------------
c
      TRBBAR = (BBAR(1)+BBAR(2)+BBAR(3))/THREE                          !  I1 / 3
                  
      EG    = TWO*C10/DETDFG                                            !  Gm / J
      EG23  = EG*TWO/THREE                  
      NS0   = ONE - NF0                                                 !  initial solid fraction

      VOLU1 = DETDFG + NS0
      VOLU2 = -DETDFG+NS0

c--- isotropic pressure  p ------------------------------------------------
      G23 = -EG*DETDFG/SIX                                              !  -Gm / 6

      PR = G23*LOG(DETDFG)/DETDFG * ( -ONE
     1      + THREE*VOLU1/VOLU2
     2      + THREE*NS0*DETDFG*LOG(DETDFG)/(VOLU2**TWO) )
     3      + Posm

c--- derivative  dp/dJ  (needed for Keff) ---------------------------------
      TERM1 = (LOG(DETDFG) - ONE) / (DETDFG**TWO)

      TERM2 = ( THREE*VOLU2*(DETDFG*LOG(DETDFG) + VOLU1)
     1       + THREE*VOLU1*LOG(DETDFG)*(DETDFG - VOLU2) )
     2       / (DETDFG**TWO * VOLU2**TWO)

      TERM3 = ( SIX*NS0*LOG(DETDFG)*VOLU2
     1       + SIX*NS0*DETDFG*LOG(DETDFG)**TWO )
     2       / (DETDFG * VOLU2**THREE)

c--- algorithmic bulk modulus --------------------------------------------
      EK = DETDFG * ( G23*(TERM1 + TERM2 + TERM3) + dPosm ) + PR

c--- Cauchy stress  -------------------------------------------------------
      DO K1 = 1, NDI                                                    ! normal components (11-33)
            STRESS(K1) = EG*( BBAR(K1) - DETDFG**(TWO/THREE) ) + PR
      END DO

      DO K1 = NDI+1, NDI+NSHR                                           ! shear components (12-23)
            STRESS(K1) = EG*BBAR(K1)
      END DO

c
c------------------------------------------------------------------------
c
c    CALCULATE THE STIFFNESS
c
c------------------------------------------------------------------------
c
      DDSDDE(1, 1)= EG23*(BBAR(1)+TRBBAR)+EK
      DDSDDE(2, 2)= EG23*(BBAR(2)+TRBBAR)+EK
      DDSDDE(3, 3)= EG23*(BBAR(3)+TRBBAR)+EK
      DDSDDE(1, 2)=-EG23*(BBAR(1)+BBAR(2)-TRBBAR)+EK
      DDSDDE(1, 3)=-EG23*(BBAR(1)+BBAR(3)-TRBBAR)+EK
      DDSDDE(2, 3)=-EG23*(BBAR(2)+BBAR(3)-TRBBAR)+EK
      DDSDDE(1, 4)= EG23*BBAR(4)/TWO
      DDSDDE(2, 4)= EG23*BBAR(4)/TWO
      DDSDDE(3, 4)=-EG23*BBAR(4)
      DDSDDE(4, 4)= EG*(BBAR(1)+BBAR(2))/TWO
      IF(NSHR.EQ.3) THEN
            DDSDDE(1, 5)= EG23*BBAR(5)/TWO
            DDSDDE(2, 5)=-EG23*BBAR(5)
            DDSDDE(3, 5)= EG23*BBAR(5)/TWO
            DDSDDE(1, 6)=-EG23*BBAR(6)
            DDSDDE(2, 6)= EG23*BBAR(6)/TWO
            DDSDDE(3, 6)= EG23*BBAR(6)/TWO
            DDSDDE(5, 5)= EG*(BBAR(1)+BBAR(3))/TWO
            DDSDDE(6, 6)= EG*(BBAR(2)+BBAR(3))/TWO
            DDSDDE(4,5)= EG*BBAR(6)/TWO
            DDSDDE(4,6)= EG*BBAR(5)/TWO
            DDSDDE(5,6)= EG*BBAR(4)/TWO
      END IF
      DO K1=1, NTENS
            DO K2=1, K1-1
                  DDSDDE(K1, K2)=DDSDDE(K2, K1)
            END DO
      END DO

C--------------------------------------------------------------------------
      RETURN
      END SUBROUTINE UMAT_MATRIX
c
c************************************************************************
c     C O L L A G E N    F I B R E S (based on a previous W. Wilson model)
c            Adapted to rebars by Andrea Malandrino
c************************************************************************
c
      SUBROUTINE UMAT_COLLAGEN(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,
     2 PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,
     3 DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
C ********************************************************************
C
      USE GLOBAL
C
      INCLUDE "ABA_PARAM.INC"
C
      DIMENSION STRESS(NTENS),STATEV(NSTATEV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4),JAY(NSTATEV),DFGR(3,3)

      CHARACTER*256 CMNAME, cwd, outdir, infoFile_T_D, infoFilePath_T_D

c      include 'common_cart.inc'

      INTEGER nf,i,j,k,l,m,n
      DOUBLE PRECISION :: de,deps(3,3),bb,cc,sqrtDD,s1,ds,dss(3,3),
     1       dA(3,3),dB(3,3,3,3),dJdF(3,3),invF(3,3),A,B(3,3),
     2       labda,dldF(3,3),TEMP1(3),TEMP2(3,3),TEMP3(3,3),dt,
     3       C(3,3,3,3),AN(3,3,3,3),SIGMA(3,3),dSIGMA(3,3,3,3),
     4       I4s(3,3,3,3),I2s(3,3),I4(3,3,3,3),JJs,s2old,sigmaold,
     5       epsold, epsi, sigmai, s2

C Define unit tensor
      DOUBLE PRECISION :: dil,djk,dik,djl

c
C PROPS
      DOUBLE PRECISION :: NF0,E1,K1,K2,E2

c************************************************************************
c     A         Q in derivation of equations in thesis [MPa]
c     B         P in derivation of equations in thesis [MPa]
c     dA        derivative of A to F [MPa]
c     dB        derivative of B to F [MPa]
c     dens      fibril density [-]
c     deps      derivative of eps to F [-]
c     de        strain rate [-]
c     dldF	derivative of labda to F [-]
c     ds        derivative of sigma to eps [MPa]
c     dSIGMA    derivative of SIGMA to F [MPa]
c     dss       derivative of sigma to F [MPa]
c     E1,k1     material constant (elastic fibrilpart) [MPa][-]
c     E2,k2     material constant (viscoelastic fibrilpart)[MPa][-]
c     elem      element number [-]
c     epsi      current fibril strain [-]
c     epsold    fibril strain from previous increment [-]
c     F 	total deformation tensor [-]
c     I4s       symmetric 4th-order unit tensor [-]
c     ip        integration point number [-]
c     labda     fibril stretch [-]
c     nf        fibril number [-]
c     s1        elastic fibril stress [MPa]
c     s2        viscoelastic fibril stress [MPa]
c     SIGMA     fibril stress tensor [MPa]
c     sigmai    current fibrils stress [MPa]
c     sigmaold  fibril stress from previous increment [MPa]
c     vec       vector representation of the fibril [-]
c************************************************************************
C    PROPS(1) - dashpot stiffness[MPa s]
C    PROPS(2) - E1 material constant (elastic fibrilpart) [MPa][-]
C    PROPS(3) - k1 material constant (elastic fibrilpart) [MPa][-]
C    PROPS(4) - E2 material constant (viscoelastic fibrilpart)[MPa][-]
C    PROPS(5) - k2 material constant (viscoelastic fibrilpart)[MPa][-]
C************************************************************************
c
c------------------------------------------------------------------------
c
c            Collagen I fibres viscoelastic material parameters
c                       From collagen.inp file
c
c------------------------------------------------------------------------
c
C Initialization of variables
C ---------------------------
      INCLUDE './src/IVD/mechanic/functions/UMATCOLLAGEN_init_var.f'
C
C--------------------------------------------------------------------------
C Definitions of the PROPS array
C--------------------------------------------------------------------------
C
      NF0 = PROPS(1)
      E1 = PROPS(2)
      k1 = PROPS(3)
      E2 = PROPS(4)
      k2 = PROPS(5)
c
c------------------------------------------------------------------------
c
c------------------------------------------------------------------------
c
c                         Define unit tensors
c
c------------------------------------------------------------------------
c
      DO i=1,3
            DO j=1,3
                  IF (i.EQ.j) THEN
                        I2s(i,j)=ONE
                  ELSE
                        I2s(i,j)=ZERO
                  END IF

                  DO k=1,3
                        DO l=1,3
                              dil=ZERO
                              djk=ZERO
                              dik=ZERO
                              djl=ZERO
                              IF (i.EQ.l) dil=ONE
                              IF (j.EQ.k) djk=ONE
                              IF (i.EQ.k) dik=ONE
                              IF (j.EQ.l) djl=ONE
                              I4s(i,j,k,l)=HALF*(dil*djk+dik*djl)
                              I4(i,j,k,l)=dil*djk
                        END DO
                  END DO
            END DO
      END DO
c
c------------------------------------------------------------------------
c
c                   Jacobian for rebar element
c
c------------------------------------------------------------------------
c
      JJs = DFGRD1(1,1)
c
c
      dJdF(1,1) = ONE
c
c------------------------------------------------------------------------
c
c 	               Current fibril strain      
c
c------------------------------------------------------------------------
c
      s2old = STATEV(1)
      sigmaold = STRESS(1)
      epsold = STRAN(1)
      epsi = STRAN(1) + DSTRAN(1)
c
c------------------------------------------------------------------------
c
c    Determine the stress matrix and derivative of the stress matrix 
c                        (for single fibril)
c
c------------------------------------------------------------------------
c
      IF (epsi.LE.1.e-6) THEN   !fibrils can only resist postive strains
            epsi=ZERO
            sigmai=ZERO
            s2=ZERO
      ELSE
c
c------------------------------------------------------------------------
c
c              UPDATE STRESS FOLLOWING ZENER MODEL  
c
c------------------------------------------------------------------------
c
c------------------------------------------------------------------------
c
c                Determine fibril stretch (labda)
c
c------------------------------------------------------------------------
c
            labda=EXP(epsi)  !labda  (verlenging)
c
c------------------------------------------------------------------------
c
c      Determine fibril stress (sigmai) and fibril stress tensor (A)
c         From Wilson et al. (2006) OsteoArthritis and Cartilage
c
c------------------------------------------------------------------------
c
c      IF (S2choice.EQ.0) THEN
c
            dt = DTIME
            de = (epsi-epsold)/dt
            bb = -de*NF0+NF0/(k2*dt)+E2
            cc = -de*E2*NF0-NF0*s2old/(k2*dt)
            sqrtDD = SQRT(bb**TWO-FOUR*cc)
            s1 = E1*(EXP(k1*epsi)-ONE)                                    !Stress of elastic part (Eq. 3)
            s2 = -HALF*bb+HALF*sqrtDD
            IF (s2 .LT. ZERO) s2=ZERO                                   !Stress of viscoelastic part (Eq. 4)
            sigmai = s1+s2                                                !Total fibril stress (Eq. 5)
            A = sigmai/(labda*JJs)	                                    !Fibril stress tensor (Eq. 2)
            STATEV(1) = s2
c
c------------------------------------------------------------------------
c   
c         Determine dlabda/dF (dldF) and deps/dF (deps)
c
c------------------------------------------------------------------------
c
            TEMP1(1) = ZERO
c
            TEMP1(1) = TEMP1(1)+DFGRD1(1,1)
c
            dldF(1,1) = TEMP1(1)/labda
            deps(1,1) = (ONE/labda)*dldF(1,1)

c
c------------------------------------------------------------------------
c
c                  Determime dsigmaf/deps (ds)
c
c------------------------------------------------------------------------
c
            IF (s2.LE.ZERO) THEN
                  ds = E1*k1*EXP(k1*epsi)
            ELSE
                  ds = E1*k1*EXP(k1*epsi)+
     1                 HALF*NF0/dt+(-HALF*bb*NF0/dt+E2*NF0/dt)/sqrtDD
            END IF
c
c------------------------------------------------------------------------
c
c           Determime dsigmaf/dF (dss) and dA/dF (dA)
c
c------------------------------------------------------------------------
c
            dss(1,1) = ds*deps(1,1)
C
            dA(1,1) = ONE/(labda*JJs)*dss(1,1)-
     1                sigmai/(labda**TWO*JJs)*dldF(1,1)-
     2                sigmai/(labda*JJs**TWO)*dJdF(1,1)
            TEMP2(1,1) = I2s(1,1)
c
c------------------------------------------------------------------------
c
c               Determine B (P in the derivations)
c
c------------------------------------------------------------------------
c
            B(1,1) = ZERO
            TEMP3(1,1) = ZERO
c
            TEMP3(1,1) = TEMP3(1,1)+DFGRD1(1,1)*TEMP2(1,1)
            B(1,1) = B(1,1)+DFGRD1(1,1)*TEMP2(1,1)*DFGRD1(1,1)
c
c------------------------------------------------------------------------
c
c      Determine dB/dF (dB) and fibril stress tensor (SIGMA)
c
c------------------------------------------------------------------------
c
            SIGMA(1,1) = A*B(1,1)
            dB(1,1,1,1) = ZERO
c
            dB(1,1,1,1) = dB(1,1,1,1)+TWO*TEMP3(1,1)
c
c------------------------------------------------------------------------
c
c                   Determine dSIGMA/dF (dSIGMA)
c
c------------------------------------------------------------------------
c
            dSIGMA(1,1,1,1) = (B(1,1)*dA(1,1)+dB(1,1,1,1)*A)
c
c------------------------------------------------------------------------
c
c Assign components from the tensors to the UMAT vectors (STRESS and DDSDDE)
c
c------------------------------------------------------------------------
c
c------------------------------------------------------------------------
c
c                           Determine A
c
c------------------------------------------------------------------------
c
            AN(1,1,1,1) = SIGMA(1,1)*dJdF(1,1)+JJs*dSIGMA(1,1,1,1)
c
c------------------------------------------------------------------------
c
c             Determine the 4-order fibril stiffness tensor C
c
c------------------------------------------------------------------------
c
            C(1,1,1,1) = ZERO
c
            C(1,1,1,1) = C(1,1,1,1)+AN(1,1,1,1)*DFGRD1(1,1)/JJs

c
c------------------------------------------------------------------------
c
c          Compute the fibril stress vector & stiffness tensor
c
c------------------------------------------------------------------------
c
            STRESS(1) = SIGMA(1,1)
            DDSDDE(1,1) = C(1,1,1,1)
c
      END IF
C--------------------------------------------------------------------------
      RETURN
      END SUBROUTINE UMAT_COLLAGEN
C--------------------------------------------------------------------------




C**********************************************************************************************************************************************
C ADDITIONAL SUBROUTINES
C**********************************************************************************************************************************************

C--------------------------------------------------------------------------
C Calculation subroutines
C--------------------------------------------------------------------------

      INCLUDE './src/IVD/mechanic/functions/UMAT_calcFunc.f'