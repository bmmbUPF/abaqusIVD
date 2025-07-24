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

C
C ***************************************************************************
C ******    ABAQUS UEL: Useful functions                               ******
C ******    auth: Estefano Muñoz-Moya                                  ******
C ******    LinkTree: https://linktr.ee/estefano23                     ******
C ******    webPage: https://estefano23.github.io/                     ******
C ******    github: estefano23                                         ******
C ******    email: estefano.munoz.moya@gmail.com                       ******
C ***************************************************************************
C

C
C-------------------------------------------------------------------------------------------------------------------------------------------------------------------
C Code
C ------------

      SUBROUTINE GETNSTEPS(NSTEPS, inpFilePath)
C-----------------------------------------------------------------------
C  GETNSTEPS:
C    Opens the input file (inpFilePath) and counts the number of
C    step definitions. It assumes that each step line begins with
C    "*STEP" (e.g., "*STEP, NLGEOM=YES, NAME=Lo_D_02_01").
C  Arguments:
C    NSTEPS      (OUT) - Number of steps found.
C    inpFilePath (IN)  - Full path to the inp file.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: NSTEPS
      CHARACTER*(*) , INTENT(IN) :: inpFilePath
      CHARACTER*256 line
      INTEGER unit, ierr, count

      unit = 20
      count = 0

      OPEN(UNIT=unit, FILE=inpFilePath, STATUS='OLD', IOSTAT=ierr)
      IF (ierr .NE. 0) THEN
         PRINT *, 'ERROR in GETNSTEPS: Unable to open file:', TRIM(inpFilePath)
         NSTEPS = 0
         RETURN
      END IF

      DO
         READ(unit, '(A)', IOSTAT=ierr) line
         IF (ierr .NE. 0) EXIT
         IF (LEN_TRIM(line) .GE. 5) THEN
            IF (line(1:5) .EQ. '*STEP') THEN
               count = count + 1
            END IF
         END IF
      END DO

      NSTEPS = count
      CLOSE(unit)

C--------------------------------------------------------------------------
      RETURN
      END SUBROUTINE GETNSTEPS

C----------------------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE GETSTEPNAMES(STEPNAMES, NSTEPS, inpFilePath)
C-----------------------------------------------------------------------
C  GETSTEPNAMES:
C    Opens the input file (inpFilePath) and for each line that begins 
C    with "*STEP" extracts the step name. It looks for the substring 
C    "NAME=" and extracts everything after "NAME=" up to the first comma 
C    (if present) or to the end of the line.
C  Arguments:
C    STEPNAMES   (OUT) - Allocated array of character*256 that will hold 
C                          the step names.
C    NSTEPS      (IN)  - Number of steps (rows) expected.
C    inpFilePath (IN)  - Full path to the inp file.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NSTEPS
      CHARACTER*256, DIMENSION(NSTEPS), INTENT(OUT) :: STEPNAMES
      CHARACTER*(*) , INTENT(IN) :: inpFilePath
      CHARACTER*256 line
      INTEGER unit, ierr, count, pos, posComma
      CHARACTER*256 stepName

      unit = 20
      count = 0

      OPEN(UNIT=unit, FILE=inpFilePath, STATUS='OLD', IOSTAT=ierr)
      IF (ierr .NE. 0) THEN
         PRINT *, 'ERROR in GETSTEPNAMES: Unable to open file:', TRIM(inpFilePath)
         RETURN
      END IF

      DO
         READ(unit, '(A)', IOSTAT=ierr) line
         IF (ierr .NE. 0) EXIT
         IF (LEN_TRIM(line) .GE. 5) THEN
            IF (line(1:5) .EQ. '*STEP') THEN
               pos = INDEX(line, 'NAME=')
               IF (pos .GT. 0) THEN
                  ! Extract substring after "NAME="
                  stepName = ADJUSTL(line(pos+5:))
                  ! If a comma is present, take substring before comma
                  posComma = INDEX(stepName, ',')
                  IF (posComma .GT. 0) THEN
                     stepName = stepName(1:posComma-1)
                  END IF
                  stepName = TRIM(stepName)
                  count = count + 1
                  IF (count .LE. NSTEPS) THEN
                     STEPNAMES(count) = stepName
                  END IF
               END IF
            END IF
         END IF
         IF (count .GE. NSTEPS) EXIT
      END DO

      CLOSE(unit)

C ------------------------------------------------------------------------
      RETURN
      END SUBROUTINE GETSTEPNAMES

C----------------------------------------------------------------------------------------------------------------------------------------------

C Calculate the identity matrix
C -----------------------------

      SUBROUTINE onem(A)
C This subroutine stores the identity matrix in the
C 3 by 3 matrix [A]
      IMPLICIT NONE
C
      INTEGER I,J

      DOUBLE PRECISION A(3,3)

C some constant parameters
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.D0, ONE=1.0D0

      DO I=1,3
            DO J=1,3
                  IF (I .EQ. J) THEN
                              A(I,J) = ONE
                  else
                              A(I,J) = ZERO
                  END IF
            END DO
      END DO

C--------------------------------------------------------------------------
      RETURN
      END SUBROUTINE onem

C----------------------------------------------------------------------------------------------------------------------------------------------

C Calculate the identity vector
C -----------------------------

      SUBROUTINE onev(A)
C This subroutine stores the identity vector in the
C 3 by 1 vector [A]
      IMPLICIT NONE
C
      INTEGER I

      DOUBLE PRECISION A(3,1)

C some constant parameters
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.D0, ONE=1.0D0
      DO I=1,3
            A(I,1) = ONE
      END DO

C--------------------------------------------------------------------------
      RETURN
      END SUBROUTINE onev

C----------------------------------------------------------------------------------------------------------------------------------------------

C Calculate the determinant of the deformation gradient
C -----------------------------------------------------------

      SUBROUTINE mDETDFG(DFGRD1,DETDFG,NSHR)
C This subroutine calculates the determinant of a 3 by 3 matrix [A]
      IMPLICIT NONE
C
      INTEGER :: NSHR
      DOUBLE PRECISION  DFGRD1(3,3), DETDFG

      DETDFG=DFGRD1(1, 1)*DFGRD1(2, 2)*DFGRD1(3, 3)
     1   -DFGRD1(1, 2)*DFGRD1(2, 1)*DFGRD1(3, 3)
      IF(NSHR.EQ.3) THEN
            DETDFG=DETDFG+DFGRD1(1, 2)*DFGRD1(2, 3)*DFGRD1(3, 1)
     1         +DFGRD1(1, 3)*DFGRD1(3, 2)*DFGRD1(2, 1)
     2         -DFGRD1(1, 3)*DFGRD1(3,1)*DFGRD1(2, 2)
     3         -DFGRD1(2, 3)*DFGRD1(3, 2)*DFGRD1(1, 1)
      END IF

C--------------------------------------------------------------------------
      RETURN
      END SUBROUTINE mDETDFG

C----------------------------------------------------------------------------------------------------------------------------------------------
C Calculate the determinant of a given 3x3 matrix
C -----------------------------------------------------------

      SUBROUTINE mDet(A,det_A)
C This subroutine calculates the determinant of a 3 by 3 matrix [A]
      IMPLICIT NONE
C
      DOUBLE PRECISION  A(3,3),det_A
C
      det_A = A(1,1)*A(2,2)*A(3,3)
     +  + A(1,2)*A(2,3)*A(3,1)
     +  + A(1,3)*A(2,1)*A(3,2)
     +  - A(3,1)*A(2,2)*A(1,3)
     +  - A(3,2)*A(2,3)*A(1,1)
     +  - A(3,3)*A(2,1)*A(1,2)

C--------------------------------------------------------------------------
      RETURN
      END SUBROUTINE mDet

C----------------------------------------------------------------------------------------------------------------------------------------------
C Calculate the inverse and determinant of a given 3x3 matrix
C -----------------------------------------------------------

      SUBROUTINE matInvDet3D(A,A_inv,A_invT,det_A,istat)
C Returns A_inv, the inverse and det_A, the determinant
C Note that the det is of the original matrix, not the
C inverse
      IMPLICIT NONE
C
      INTEGER istat
C
      DOUBLE PRECISION A(3,3), A_inv(3,3), A_invT(3,3)
      DOUBLE PRECISION det_A, det_A_inv

C initialize
      A_inv = 0.D0
      A_invT = 0.D0
      det_A = 0.D0
      istat = 1

C Calculate the determinant of the matrix
      det_A = ABS(A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3)) -
     +        A(2,1)*(A(1,2)*A(3,3) - A(3,2)*A(1,3)) +
     +        A(3,1)*(A(1,2)*A(2,3) - A(2,2)*A(1,3)))

C Check if the determinant is zero
      if (det_A .le. 0.D0) then
            write(*,*) 'WARNING: subroutine matInvDet3D:'
            write(*,*) 'WARNING: det of mat=',det_A
            istat = 0
            return
      end if

C Calculate the inverse of the matrix
      det_A_inv = ABS(1.D0/det_A)

C Calculate the inverse of the matrix
      A_inv(1,1) = det_A_inv*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
      A_inv(1,2) = det_A_inv*(A(3,2)*A(1,3)-A(1,2)*A(3,3))
      A_inv(1,3) = det_A_inv*(A(1,2)*A(2,3)-A(2,2)*A(1,3))
      A_inv(2,1) = det_A_inv*(A(3,1)*A(2,3)-A(2,1)*A(3,3))
      A_inv(2,2) = det_A_inv*(A(1,1)*A(3,3)-A(3,1)*A(1,3))
      A_inv(2,3) = det_A_inv*(A(2,1)*A(1,3)-A(1,1)*A(2,3))
      A_inv(3,1) = det_A_inv*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
      A_inv(3,2) = det_A_inv*(A(3,1)*A(1,2)-A(1,1)*A(3,2))
      A_inv(3,3) = det_A_inv*(A(1,1)*A(2,2)-A(2,1)*A(1,2))

C Calculate the transpose of the inverse of the matrix
      A_invT = TRANSPOSE(A_inv)

C--------------------------------------------------------------------------
      RETURN
      END SUBROUTINE matInvDet3D

C----------------------------------------------------------------------------------------------------------------------------------------------
C The D/H/M/S between two dates
C -----------------------------------------------------------

      SUBROUTINE diff_DATE_AND_TIME(date_ini, date_fin, date_diff)
C This subroutine calculates the difference between two dates and times
      IMPLICIT NONE
C
C Declare variables to store date and time components:
C date_ini: Initial date and time values
C date_fin: Current date and time values
C date_diff: Difference between the initial and current date and time values
      INTEGER :: date_ini(8), date_fin(8), date_diff(5)

C date_diff(1): Days
C date_diff(2): Hours
C date_diff(3): Minutes
C date_diff(4): Seconds
C date_diff(5): Total seconds

C Initialize date_diff
      date_diff = 0

C Compute elapsed time in total seconds
      date_diff(5) = (date_fin(1) - date_ini(1))
     1 * 31536000 + (date_fin(2) - date_ini(2))
     2 * 2628000 + (date_fin(3) - date_ini(3))
     3 * 86400 + (date_fin(5) - date_ini(5))
     4 * 3600 + (date_fin(6) - date_ini(6))
     5 * 60 + (date_fin(7) - date_ini(7))
      
      ! Convert total seconds to D/H/M/S
      date_diff(1) = date_diff(5) / (24 * 3600)                     ! Days
      date_diff(2) = (date_diff(5) - date_diff(1) * 24 * 3600)
     1 / 3600  ! Hours
      date_diff(3) = (date_diff(5) - date_diff(1) * 24 * 3600
     1 - date_diff(2) * 3600) / 60  ! Minutes
      date_diff(4) = date_diff(5) - date_diff(1) * 24 * 3600
     1 - date_diff(2) * 3600 - date_diff(3) * 60  ! Seconds

C--------------------------------------------------------------------------
      RETURN
      END SUBROUTINE diff_DATE_AND_TIME

C----------------------------------------------------------------------------------------------------------------------------------------------
      SUBROUTINE GETUNIQUESTEPS(sdvFolder, JOBNAME,
     1                          outdir, uniqueSteps, nstepsSP)
C***********************************************************************
C  SUBROUTINE GETUNIQUESTEPS
C
C  This subroutine examines all files in the folder sdvFolder and
C  collects the unique step numbers from filenames that begin with
C  JOBNAME//'_SDV.' (e.g. GENERIC_Transport_DiffReac_SDV.2.3.SP).
C  The extracted step numbers (assumed to be integers) are returned in
C  the array uniqueSteps and the count in nstepsSP.
C
C  Inputs:
C     sdvFolder  - Character string containing the folder path.
C     JOBNAME    - Character string with the job name.
C     outdir     - Absolute path where the temporary file list will be written.
C
C  Outputs:
C     uniqueSteps  - Integer array (dimensioned to at least 100) of unique step numbers.
C     nstepsSP       - Number of unique step numbers found.
C***********************************************************************
      IMPLICIT NONE
      CHARACTER*(*) sdvFolder, JOBNAME, outdir
      CHARACTER*256 command, fileName, filelistName
      INTEGER uniqueSteps(100)
      INTEGER nstepsSP
      INTEGER ios, pos, posPeriod, stepValue, i, j

      nstepsSP = 0

C Ensure outdir ends with a slash.
      IF (outdir(LEN(outdir):LEN(outdir)) .NE. '/') THEN
         outdir = TRIM(outdir) // '/'
      END IF

C Build the file list file name using JOBNAME to avoid collisions.
      filelistName = TRIM(outdir) // TRIM(JOBNAME) // '_filelist.txt'

C Use ls -1 to force one file per line.
      command = 'ls -1 ' // TRIM(sdvFolder) // ' > ' // filelistName
      CALL SYSTEM(command)

      OPEN(20, FILE=filelistName, STATUS='OLD', IOSTAT=ios)
      IF (ios .NE. 0) THEN
         WRITE(*,*) 'Error: Could not open ', filelistName
         RETURN
      END IF

 10   CONTINUE
      READ(20, '(A)', END=100, IOSTAT=ios) fileName
      IF (ios .NE. 0) GOTO 10

      IF (INDEX(fileName, TRIM(JOBNAME)//'_SDV.') .EQ. 1) THEN
            ! Determine the position where the step number starts.
            pos = LEN(TRIM(JOBNAME)) + LEN('_SDV.') + 1
            ! Find the first period in the substring starting at pos.
            posPeriod = INDEX(fileName(pos:), '.')
            IF (posPeriod .GT. 0) THEN
                  READ(fileName(pos: pos+posPeriod-2), '(I10)',
     1                 IOSTAT=ios) stepValue
            ELSE
            stepValue = -1
            END IF
            IF (ios .EQ. 0) THEN
                  j = 0
                  DO i = 1, nstepsSP
                        IF (uniqueSteps(i) .EQ. stepValue) THEN
                              j = 1
                              EXIT
                        END IF
                  END DO
                  IF (j .EQ. 0) THEN
                        nstepsSP = nstepsSP + 1
                        uniqueSteps(nstepsSP) = stepValue
                  END IF
            END IF
      END IF
      GOTO 10
100   CONTINUE
      CLOSE(20)

C ----- Now sort the uniqueSteps array in ascending order (simple bubble sort) -----
      DO i = 1, nstepsSP-1
         DO j = i+1, nstepsSP
            IF (uniqueSteps(i) > uniqueSteps(j)) THEN
               stepValue = uniqueSteps(i)
               uniqueSteps(i) = uniqueSteps(j)
               uniqueSteps(j) = stepValue
            END IF
         END DO
      END DO

C--------------------------------------------------------------------------
      RETURN
      END SUBROUTINE GETUNIQUESTEPS

C----------------------------------------------------------------------------------------------------------------------------------------------
      SUBROUTINE NcommasCSV(fileName, nCommas, rowWanted)
C=======================================================================
C  NcommasCSV:
C    Reads the specified row ("rowWanted") from the given file ("fileName")
C    and counts how many commas (',') are in that line.
C    The result is returned in "nCommas".
C
C    Example usage in UEXTERNALDB:
C      CALL NcommasCSV(FREQfilePath, NActivities, 1)
C      NActivities = NActivities + 1
C    (since columns = commas + 1)
C=======================================================================

      CHARACTER*(*) fileName
      INTEGER nCommas
      INTEGER rowWanted

C     Local variables
      INTEGER i, iErr, iLine
      CHARACTER*256 line

C     Initialize
      nCommas = 0

C     Open the file
      OPEN(UNIT=15, FILE=fileName, STATUS='OLD', IOSTAT=iErr)
      IF (iErr .NE. 0) THEN
         PRINT *, '*** NcommasCSV Could not open file: ', TRIM(fileName)
         RETURN
      END IF

C     Skip to the desired row
      DO iLine = 1, rowWanted
         READ(15, '(A)', IOSTAT=iErr, END=999) line
         IF (iErr .NE. 0) THEN
            PRINT *, '*** NcommasCSV Error reading line ', iLine
            GOTO 999
         END IF
      END DO

C     Now "line" holds the desired row.
C     Count commas:
      DO i = 1, LEN_TRIM(line)
         IF (line(i:i) .EQ. ',') THEN
            nCommas = nCommas + 1
         END IF
      END DO

 999  CONTINUE
      CLOSE(UNIT=15)

C--------------------------------------------------------------------------
      RETURN
      END SUBROUTINE NcommasCSV

C----------------------------------------------------------------------------------------------------------------------------------------------
      SUBROUTINE readActivitiesCSV(filePath, NActivities,
     1                             STEP_LIST_STR, HOUR_RANGE_STR,
     2                             FREQ_VAL, PN_TIME_file,
     3                             PN_NActivities,
     4                             PN_HOUR_RANGE_STR, PN_FREQ_VAL)

C---------------------------------------------------------------------
C  readActivitiesCSV:
C    Reads four lines from the CSV at "filePath". The file format is
C      1) Activity names (e.g., "Walking  Standing  Sitting")
C      2) Steps (e.g., "5-9-13  5-9-13  5-9-13")
C      3) Hour ranges (e.g., "4-5  5-9  9-12")
C      4) Frequencies (e.g., "1.8  0.9  0.0001")
C
C    We assume "NActivities" columns in each line. This subroutine:
C      - opens the file,
C      - reads these 4 lines,
C      - discards line 1 (or you can store it if needed),
C      - parses line 2 into STEP_LIST_STR(*),
C      - parses line 3 into HOUR_RANGE_STR(*),
C      - parses line 4 into FREQ_VAL(*),
C
C  Arguments:
C    filePath        (IN)   full path to CSV
C    NActivities     (IN)   number of columns (fields) in each line
C    STEP_LIST_STR   (OUT)  e.g., "5-9-13"
C    HOUR_RANGE_STR  (OUT)  e.g., "4-5"
C    FREQ_VAL        (OUT)  numeric frequencies
C
C---------------------------------------------------------------------
      CHARACTER*(*) filePath
      INTEGER NActivities, PN_NActivities
      CHARACTER*256 STEP_LIST_STR(*)
      CHARACTER*256 HOUR_RANGE_STR(*)
      DOUBLE PRECISION FREQ_VAL(*)
      CHARACTER*256 PN_HOUR_RANGE_STR(*)
      DOUBLE PRECISION PN_FREQ_VAL(*)
      INTEGER PN_TIME_file

      INTEGER ioErr
      CHARACTER*256 line1, line2, line3, line4, line5, line6
      CHARACTER*256 line7, line8, line9, line10, line11, line12

      ! Open the file
      OPEN(UNIT=30, FILE=filePath, STATUS='OLD', IOSTAT=ioErr)
      IF (ioErr .NE. 0) THEN
         PRINT *, '*** Could not open file:', TRIM(filePath)
         PRINT *, ' '
         RETURN
      END IF

      ! Read 6 lines
      READ(30, '(A)', IOSTAT=ioErr) line1
      READ(30, '(A)', IOSTAT=ioErr) line2
      READ(30, '(A)', IOSTAT=ioErr) line3
      READ(30, '(A)', IOSTAT=ioErr) line4
      READ(30, '(A)', IOSTAT=ioErr) line5
      READ(30, '(A)', IOSTAT=ioErr) line6
      READ(30, '(A)', IOSTAT=ioErr) line7
      READ(30, '(A)', IOSTAT=ioErr) line8
      READ(30, '(A)', IOSTAT=ioErr) line9
      READ(30, '(A)', IOSTAT=ioErr) line10
      READ(30, '(A)', IOSTAT=ioErr) line11
      READ(30, '(A)', IOSTAT=ioErr) line12
      CLOSE(30)

      IF (ioErr .NE. 0) THEN
         PRINT *, '*** Error reading lines from', TRIM(filePath)
         PRINT *, ' '
         RETURN
      END IF

C     -----------------------------------------------------------
C     We are ignoring the third line (activity names). If you want
C     to store or process them, you can parse line3 similarly.
C     e.g., READ(line3, *) (ACTIVITY_NAMES(i), i=1,NActivities)

C     Parse the line 3 into PN_TIME_file (contains only one value)
      READ(line3, *, IOSTAT=ioErr) PN_TIME_file
      IF (ioErr .NE. 0) THEN
         PRINT *, '*** Error parsing time values, line3'
         PRINT *, ' '
         RETURN
      END IF

C     Parse the line 5 into PN_HOUR_RANGE_STR
      READ(line5, *, IOSTAT=ioErr)
     1     (PN_HOUR_RANGE_STR(i), i=1,PN_NActivities)
      IF (ioErr .NE. 0) THEN
         PRINT *, '*** Error parsing hour ranges, line5'
         PRINT *, ' '
         RETURN
      END IF

C    Parse the line 6 into PN_FREQ_VAL
      READ(line6, *, IOSTAT=ioErr)
     1     (PN_FREQ_VAL(i), i=1,PN_NActivities)
      IF (ioErr .NE. 0) THEN
         PRINT *, '*** Error parsing frequency values, line6'
         PRINT *, ' '
         RETURN
      END IF

C     Parse the line 10 into STEP_LIST_STR
      READ(line10, *, IOSTAT=ioErr) (STEP_LIST_STR(i), i=1,NActivities)
      IF (ioErr .NE. 0) THEN
         PRINT *, '*** Error parsing step strings, line10'
         PRINT *, ' '
         RETURN
      END IF

C     Parse the line 11 into HOUR_RANGE_STR
      READ(line11, *, IOSTAT=ioErr) (HOUR_RANGE_STR(i), i=1,NActivities)
      IF (ioErr .NE. 0) THEN
         PRINT *, '*** Error parsing hour ranges, line11'
         PRINT *, ' '
         RETURN
      END IF

C     Parse the line 12 into FREQ_VAL
      READ(line12, *, IOSTAT=ioErr) (FREQ_VAL(i), i=1,NActivities)
      IF (ioErr .NE. 0) THEN
         PRINT *, '*** Error parsing frequency values, line12'
         PRINT *, ' '
         RETURN
      END IF

C ------------------------------------------------------------
      RETURN
      END SUBROUTINE readActivitiesCSV

C -----------------------------------------------------------------------------------------------------------------------------------------------
      LOGICAL FUNCTION STEP_APPLIES(currentStep, stepStr, nstepToken)

C---------------------------------------------------------------------
C  STEP_APPLIES:
C    Checks if the integer "currentStep" is present
C    in the dash-separated string "stepStr".
C    Example: stepStr="5-9-13". Replaces '-' with ' '
C    then reads tokens "5", "9", "13".
C---------------------------------------------------------------------

      INTEGER currentStep
      CHARACTER*(*) stepStr
      INTEGER nstepToken

      INTEGER i, ioErr
      INTEGER steps(nstepToken)  ! read up to nstepToken integer tokens
      CHARACTER*256 localStr
      localStr = stepStr

C Replace dash with space, so list-directed read sees separate tokens
      CALL REPLACE_CHAR(localStr, '-', ' ')

C Initialize steps array to a sentinel value
      DO i=1,nstepToken
         steps(i) = -1
      END DO

C Parse up to nstepToken integers from localStr
      READ(localStr, *, IOSTAT=ioErr) (steps(i), i=1,nstepToken)

C Default is .FALSE.
      STEP_APPLIES = .FALSE.

C If read succeeded, check if currentStep matches any token
      IF (ioErr .EQ. 0) THEN
         DO i=1,nstepToken
            IF (steps(i) .EQ. currentStep) THEN
               STEP_APPLIES = .TRUE.
               RETURN
            END IF
         END DO
      END IF

C ---------------------------------------------------------------------------
      RETURN
      END FUNCTION STEP_APPLIES

C -----------------------------------------------------------------------------------------------------------------------------------------------
      LOGICAL FUNCTION HOUR_APPLIES(currentHour, hourStr)

C---------------------------------------------------------------------
C  HOUR_APPLIES:
C    Checks if the double precision "currentHour" is within the
C    dash-separated hour range in the string "hourStr".
C    Example: hourStr="4-5" -> range [4,5).
C---------------------------------------------------------------------

      DOUBLE PRECISION currentHour
      CHARACTER*(*) hourStr

      CHARACTER*256 localStr
      DOUBLE PRECISION startHr, endHr
      INTEGER ioErr

      localStr = hourStr
      CALL REPLACE_CHAR(localStr, '-', ' ')

C Default range = 0..0 if we fail
      startHr = ZERO
      endHr   = ZERO

C Parse 2 doubles from localStr
      READ(localStr, *, IOSTAT=ioErr) startHr, endHr

      HOUR_APPLIES = .FALSE.
      IF (ioErr .EQ. 0) THEN
         IF ( (currentHour .GE. startHr) .AND.
     1        (currentHour .LE. endHr) ) THEN
            HOUR_APPLIES = .TRUE.
         END IF
      END IF

C ---------------------------------------------------------------------------
      RETURN
      END FUNCTION HOUR_APPLIES

C -----------------------------------------------------------------------------------------------------------------------------------------------
      SUBROUTINE REPLACE_CHAR(str, oldC, newC)
      CHARACTER*(*) str
      CHARACTER*1 oldC, newC
      INTEGER i, L

      L = LEN_TRIM(str)
      DO i=1,L
         IF (str(i:i) .EQ. oldC) THEN
            str(i:i) = newC
         END IF
      END DO

C ---------------------------------------------------------------------------
      RETURN
      END SUBROUTINE REPLACE_CHAR

C -----------------------------------------------------------------------------------------------------------------------------------------------
      INTEGER FUNCTION COUNT_TOKENS(stepStr, delimiter)
      IMPLICIT NONE
      CHARACTER*(*) stepStr
      CHARACTER*1   delimiter

      INTEGER i, L, nDash
      L = LEN_TRIM(stepStr)

      nDash = 0
      DO i = 1, L
            IF (stepStr(i:i) .EQ. delimiter) THEN
                  nDash = nDash + 1
            END IF
      END DO

      IF (L .GT. 0) THEN
      COUNT_TOKENS = nDash + 1
      ELSE
      COUNT_TOKENS = 0
      END IF

C ---------------------------------------------------------------------------
      RETURN
      END FUNCTION COUNT_TOKENS

C -----------------------------------------------------------------------------------------------------------------------------------------------
      SUBROUTINE PN_HOURS(PN_FREQ_file, PN_TIME_file, hourStr, FREQ)

C---------------------------------------------------------------------
C  PN_HOURS:
C    Set the PN_FREQ_file array to the value FREQ for the corresponding hours
C---------------------------------------------------------------------

      INTEGER PN_TIME_file
      DOUBLE PRECISION PN_FREQ_file(PN_TIME_file), FREQ
      CHARACTER*(*) hourStr

      CHARACTER*256 localStr
      INTEGER startHr, endHr
      INTEGER ioErr, I

      localStr = hourStr
      CALL REPLACE_CHAR(localStr, '-', ' ')

C Default range = 0..0 if we fail
      startHr = 0
      endHr   = 0

C Parse 2 doubles from localStr
      READ(localStr, *, IOSTAT=ioErr) startHr, endHr

      IF (ioErr .EQ. 0) THEN
         DO I = startHr, endHr - 1
            PN_FREQ_file(I) = FREQ
         END DO
      END IF

C ---------------------------------------------------------------------------
      RETURN
      END SUBROUTINE PN_HOURS