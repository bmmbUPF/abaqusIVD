C
C Write simulation info to file at the stat of the analisys
      IF (LOP .EQ. 0) THEN  
      OPEN(15,FILE=infoFilePath_T_D, STATUS='OLD', POSITION='APPEND')

C Write the design to the file
      WRITE(15,*) ' '
      WRITE(15,*) '  ___  ___  ________     ___    ___'
      WRITE(15,*) ' |\  \|\  \|\   __  \   |\  \  /  /|'
      WRITE(15,*) ' \ \  \\\  \ \  \|\  \  \ \  \/  / /'
      WRITE(15,*) '  \ \   __  \ \   __  \  \ \    / /'
      WRITE(15,*) '   \ \  \ \  \ \  \ \  \  /     \/'
      WRITE(15,*) '    \ \__\ \__\ \__\ \__\/  /\   \'
      WRITE(15,*) '     \|__|\|__|\|__|\|__/__/ /\ __\'
      WRITE(15,*) '                        |__|/ \|__|'
      WRITE(15,*) '      01001000 01100001 01111000'
      WRITE(15,*) '  ___  ___  ________     ___    ___ ___    ___'
      WRITE(15,*) ' |\  \|\  \|\   __  \   |\  \  /  /|\  \  /  /|'
      WRITE(15,*) ' \ \  \\\  \ \  \|\  \  \ \  \/  / | \  \/  / /'
      WRITE(15,*) '  \ \   __  \ \   __  \  \ \    / / \ \    / /'
      WRITE(15,*) '   \ \  \ \  \ \  \ \  \  /     \/   /     \/'
      WRITE(15,*) '    \ \__\ \__\ \__\ \__\/  /\   \  /  /\   \'
      WRITE(15,*) '     \|__|\|__|\|__|\|__/__/ /\ __\/__/ /\ __\'
      WRITE(15,*) '                       |__|/ \|__||__|/ \|__|'
      WRITE(15,*) '        01001000 01100001 01111000 01111000'
      WRITE(15,*) '   _______  ________'
      WRITE(15,*) '  /  ___  \|\_____  \'
      WRITE(15,*) ' /__/|_/  /\|____|\ /_'
      WRITE(15,*) ' |__|//  / /     \|\  \'
      WRITE(15,*) '     /  /_/__   __\_\  \'
      WRITE(15,*) '    |\________\|\_______\'
      WRITE(15,*) '     \|_______|\|_______|'
      WRITE(15,*) '      00110010 00110011'      
      WRITE(15,*) ' '
      WRITE(15,*) '****************************************************'
      WRITE(15,*) ' ABAQUS SUBROUTINE'
      WRITE(15,*) ' '
      WRITE(15,*) ' Diffusion-Reaction Transport IVD Simulation'
      WRITE(15,*) ' Obs: Coupled with a Mechanical IVD Simulation'
      WRITE(15,*) ' Type: User Element (UEL)'
      WRITE(15,*) ' '
      WRITE(15,*) ' Auth: Estefano Muñoz-Moya'
      WRITE(15,*) ' LinkTree: https://linktr.ee/estefano23'
      WRITE(15,*) ' Web Page: https://estefano23.github.io/'
      WRITE(15,*) ' GitHub: estefano23'
      WRITE(15,*) ' Email: estefano.munoz.moya@gmail.com'
      WRITE(15,*) ' HaxHaxx23'
      WRITE(15,*) ' 01001000 01100001 01111000'
      WRITE(15,*) ' 01001000 01100001 01111000 01111000'
      WRITE(15,*) ' 00110010 00110011'
      WRITE(15,*) '****************************************************'
      WRITE(15,*) ' '
      WRITE(15,*) '---------------'
      WRITE(15,*) 'SIMULATION INFO'
      WRITE(15,*) '---------------'
      WRITE(15,*) ' '
      WRITE(15,'(A, I4.4, A, I2.2, A, I2.2, A)') ' Date: '
     1, date_ini(1), '-', date_ini(2), '-'
     2, date_ini(3)
      WRITE(15,'(A, I2.2, A, I2.2, A, I2.2, A, I3.3)') ' Time: '
     1, date_ini(5), ':', date_ini(6), ':'
     2, date_ini(7)
      WRITE(15,*) 'The time difference from UTC in minutes: '
     1, date_ini(4)
      WRITE(15,*) 'Current job name: ', JOBNAME
      WRITE(15,*) 'Current job directory: '
      WRITE(15,*) outdir
      WRITE(15,*) 'Temporary directory: '
      WRITE(15,*) cwd
      WRITE(15,*) ' '
      WRITE(15,*) '----------------------------------------------------'
      WRITE(15,*) ' '
      WRITE(15,*) 'TYPE OF ELEMENT: C3D20'
      WRITE(15,*) '----------------------'
      WRITE(15,*) ' '
      WRITE(15,*) ' Normalized Nodal Coordinates'
      WRITE(15,*) ' for the shape functions (N)'
      WRITE(15,*) '       x   y    z'
      WRITE(15,*) 'N(1)  -1  -1   -1'
      WRITE(15,*) 'N(2)   1  -1   -1                          \y'
      WRITE(15,*) 'N(3)   1   1   -1                \z       /'
      WRITE(15,*) 'N(4)  -1   1   -1                |       /'
      WRITE(15,*) 'N(5)  -1  -1    1          8-----|15--------7'
      WRITE(15,*) 'N(6)   1  -1    1         /|     |     /   /|'
      WRITE(15,*) 'N(7)   1   1    1        / |     |    /   / |'
      WRITE(15,*) 'N(8)  -1   1    1      16  |     |   /  14  |'
      WRITE(15,*) 'N(9)   0  -1   -1      /  20     |  /   /   19'
      WRITE(15,*) 'N(10)  1   0   -1     /    |     | /   /    |'
      WRITE(15,*) 'N(11)  0   1   -1    5-------13-------6     |'
      WRITE(15,*) 'N(12) -1   0   -1    |     |     +----|--------\x'
      WRITE(15,*) 'N(13)  0  -1    1    |     4-------11-|-----3'
      WRITE(15,*) 'N(14)  1   0    1    |    /           |    /'
      WRITE(15,*) 'N(15)  0   1    1   17   /           18   /'
      WRITE(15,*) 'N(16) -1   0    1    |  12            |  10'
      WRITE(15,*) 'N(17) -1  -1    0    | /              | /'
      WRITE(15,*) 'N(18)  1  -1    0    |/               |/'
      WRITE(15,*) 'N(19)  1   1    0    1--------9-------2'
      WRITE(15,*) 'N(20) -1   1    0'
      WRITE(15,*) ' '
      WRITE(15,*) ' Normalized Integration Point coordinates (Cip)'
      WRITE(15,*) '         ξ   η    ζ'
      WRITE(15,*) 'Cip(1)  -1  -1   -1'
      WRITE(15,*) 'Cip(2)   0  -1   -1'
      WRITE(15,*) 'Cip(3)   1  -1   -1'
      WRITE(15,*) 'Cip(4)  -1   0   -1                          \η'
      WRITE(15,*) 'Cip(5)   0   0   -1                \ζ       /'
      WRITE(15,*) 'Cip(6)   1   0   -1                |       /'
      WRITE(15,*) 'Cip(7)  -1   1   -1         25-----|26-------27'
      WRITE(15,*) 'Cip(8)   0   1   -1         /|     |     /   /|'
      WRITE(15,*) 'Cip(9)   1   1   -1        / |     |    /   / |'
      WRITE(15,*) 'Cip(10) -1  -1    0      22--|----23------24  |'
      WRITE(15,*) 'Cip(11)  0  -1    0      /  16-----|17/---/---18'
      WRITE(15,*) 'Cip(12)  1  -1    0     /   /|     | /   /   /|'
      WRITE(15,*) 'Cip(13) -1   0    0   19-------20------21   / |'
      WRITE(15,*) 'Cip(14)  0   0    0    |  13------14----|-15-----\ξ'
      WRITE(15,*) 'Cip(15)  1   0    0    | /   7-------8--|-/---9'
      WRITE(15,*) 'Cip(16) -1   1    0    |/   /           |/   /'
      WRITE(15,*) 'Cip(17)  0   1    0   10--------11-----12   /'
      WRITE(15,*) 'Cip(18)  1   1    0    |  4 -------5----|--6'
      WRITE(15,*) 'Cip(19) -1  -1    1    | /              | /'
      WRITE(15,*) 'Cip(20)  0  -1    1    |/               |/'
      WRITE(15,*) 'Cip(21)  1  -1    1    1--------2-------3'
      WRITE(15,*) 'Cip(22) -1   0    1   Multiplicative Factor:'
      WRITE(15,*) 'Cip(23)  0   0    1      α = SQR(3/5)'
      WRITE(15,*) 'Cip(24)  1   0    1   ALPHA=0.774596669241483'
      WRITE(15,*) 'Cip(25) -1   1    1'
      WRITE(15,*) 'Cip(26)  0   1    1'
      WRITE(15,*) 'Cip(27)  1   1    1'
      WRITE(15,*) ' '
      WRITE(15,*) '----------------------------------------------------'
      WRITE(15,*) ' '
      WRITE(15,*) 'PROBLEM INFO'
      WRITE(15,*) '------------'
      WRITE(15,*) ' '
      WRITE(15,*) 'Number of dimensions: ', NDIM
      WRITE(15,*) 'Number of elements: ', numElem
      WRITE(15,*) 'ID Offset original elements vs UEL: ', ElemOffset
      WRITE(15,*) 'Number of integration points per element: ', numInt
      WRITE(15,*) 'Total number of integration points: ', numElemInt
      WRITE(15,*) 'Larget id number for the UEL elements: ', numIdElem
      IF (FREQfileSwitch .EQ. 1) THEN
            WRITE(15,*) 'Frequency file: ', FREQfile
            WRITE(15,*) 'Frequency file path: ', FREQfilePath
            WRITE(15,*) ' '
            WRITE(15,*) 'For PN eq. during the simulation:'
            WRITE(15,*) 'Number of activities: ', NActivities
            WRITE(15,*) 'STEP_LIST_STR:  ', STEP_LIST_STR
            WRITE(15,*) 'HOUR_RANGE_STR: ', HOUR_RANGE_STR
            WRITE(15,*) 'FREQ_VAL:       ', FREQ_VAL
            WRITE(15,*) ' '
            WRITE(15,*) 'For PN eq. STEP:'
            WRITE(15,*) 'Number of activities: ', PN_NActivities
            WRITE(15,*) 'Hours of simulation: ', PN_TIME_file
            WRITE(15,*) 'Frequency per hour:'
            DO I = 1, PN_TIME_file
                  WRITE(15,*) ' ', I, ': ', PN_FREQ_file(I)
            ENDDO
      ELSE
            WRITE(15,*) 'Frequency file: ', 'None'
      END IF
      WRITE(15,*) ' '
      WRITE(15,*) 'Number of Step: ', NSTEPS
C NOW THE STEP NAMES STORED IN STEPNAMES(NSTEPS)
      DO I=1,NSTEPS
            WRITE(15,'(A,I3,A)') ' Step(',I,'): ', STEPNAMES(I)
      ENDDO
      WRITE(15,*) ' '
      WRITE(15,*) '****************************************************'
      WRITE(15,*) '   THE SDVs ARE STORED IN THE UVARM SUBROUTINE'
      WRITE(15,*) '****************************************************'
      WRITE(15,*) ' '
      CALL Initialize_UVARM_NAMES()
      DO I=1,ngSdv
            WRITE(15,'(A,I3,A)') ' UVARM(',I,'):'
            WRITE(15,*) UVARM_NAMES(I)
      ENDDO
      WRITE(15,*) ' '
      CLOSE(15)

      PRINT *, ' '
      PRINT *, '  ___  ___  ________     ___    ___'
      PRINT *, ' |\  \|\  \|\   __  \   |\  \  /  /|'
      PRINT *, ' \ \  \\\  \ \  \|\  \  \ \  \/  / /'
      PRINT *, '  \ \   __  \ \   __  \  \ \    / /'
      PRINT *, '   \ \  \ \  \ \  \ \  \  /     \/'
      PRINT *, '    \ \__\ \__\ \__\ \__\/  /\   \'
      PRINT *, '     \|__|\|__|\|__|\|__/__/ /\ __\'
      PRINT *, '                        |__|/ \|__|'
      PRINT *, '      01001000 01100001 01111000'
      PRINT *, '  ___  ___  ________     ___    ___ ___    ___'
      PRINT *, ' |\  \|\  \|\   __  \   |\  \  /  /|\  \  /  /|'
      PRINT *, ' \ \  \\\  \ \  \|\  \  \ \  \/  / | \  \/  / /'
      PRINT *, '  \ \   __  \ \   __  \  \ \    / / \ \    / /'
      PRINT *, '   \ \  \ \  \ \  \ \  \  /     \/   /     \/'
      PRINT *, '    \ \__\ \__\ \__\ \__\/  /\   \  /  /\   \'
      PRINT *, '     \|__|\|__|\|__|\|__/__/ /\ __\/__/ /\ __\'
      PRINT *, '                       |__|/ \|__||__|/ \|__|'
      PRINT *, '        01001000 01100001 01111000 01111000'
      PRINT *, '   _______  ________'
      PRINT *, '  /  ___  \|\_____  \'
      PRINT *, ' /__/|_/  /\|____|\ /_'
      PRINT *, ' |__|//  / /     \|\  \'
      PRINT *, '     /  /_/__   __\_\  \'
      PRINT *, '    |\________\|\_______\'
      PRINT *, '     \|_______|\|_______|'
      PRINT *, '      00110010 00110011'      
      PRINT *, ' '
      PRINT *, '****************************************************'
      PRINT *, ' ABAQUS SUBROUTINE'
      PRINT *, ' '
      PRINT *, ' Diffusion-Reaction Transport IVD Simulation'
      PRINT *, ' Obs: Coupled with a Mechanical IVD Simulation'
      PRINT *, ' Type: User Element (UEL)'
      PRINT *, ' '
      PRINT *, ' Auth: Estefano Muñoz-Moya'
      PRINT *, ' LinkTree: https://linktr.ee/estefano23'
      PRINT *, ' Web Page: https://estefano23.github.io/'
      PRINT *, ' GitHub: estefano23'
      PRINT *, ' Email: estefano.munoz.moya@gmail.com'
      PRINT *, ' HaxHaxx23'
      PRINT *, ' 01001000 01100001 01111000'
      PRINT *, ' 01001000 01100001 01111000 01111000'
      PRINT *, ' 00110010 00110011'
      PRINT *, '****************************************************'
      PRINT *, ' '
      PRINT *, '---------------'
      PRINT *, 'SIMULATION INFO'
      PRINT *, '---------------'
      PRINT *, ' '
      PRINT '(A, I4.4, A, I2.2, A, I2.2, A)', ' Date: '
     1, date_ini(1), '-', date_ini(2), '-'
     2, date_ini(3)
      PRINT '(A, I2.2, A, I2.2, A, I2.2, A, I3.3)', ' Time: '
     1, date_ini(5), ':', date_ini(6), ':'
     2, date_ini(7)
      PRINT *, 'The time difference from UTC in minutes: '
     1, date_ini(4)
      PRINT *, 'Current job name: ', JOBNAME
      PRINT *, 'Current job directory: '
      PRINT *, outdir
      PRINT *, 'Temporary directory: '
      PRINT *, cwd
      PRINT *, ' '
      PRINT *, '----------------------------------------------------'
      PRINT *, ' '
      PRINT *, 'TYPE OF ELEMENT: C3D20'
      PRINT *, '----------------------'
      PRINT *, ' '
      PRINT *, ' Normalized Nodal Coordinates'
      PRINT *, ' for the shape functions (N)'
      PRINT *, '       x   y    z'
      PRINT *, 'N(1)  -1  -1   -1'
      PRINT *, 'N(2)   1  -1   -1                          \y'
      PRINT *, 'N(3)   1   1   -1                \z       /'
      PRINT *, 'N(4)  -1   1   -1                |       /'
      PRINT *, 'N(5)  -1  -1    1          8-----|15--------7'
      PRINT *, 'N(6)   1  -1    1         /|     |     /   /|'
      PRINT *, 'N(7)   1   1    1        / |     |    /   / |'
      PRINT *, 'N(8)  -1   1    1      16  |     |   /  14  |'
      PRINT *, 'N(9)   0  -1   -1      /  20     |  /   /   19'
      PRINT *, 'N(10)  1   0   -1     /    |     | /   /    |'
      PRINT *, 'N(11)  0   1   -1    5-------13-------6     |'
      PRINT *, 'N(12) -1   0   -1    |     |     +----|--------\x'
      PRINT *, 'N(13)  0  -1    1    |     4-------11-|-----3'
      PRINT *, 'N(14)  1   0    1    |    /           |    /'
      PRINT *, 'N(15)  0   1    1   17   /           18   /'
      PRINT *, 'N(16) -1   0    1    |  12            |  10'
      PRINT *, 'N(17) -1  -1    0    | /              | /'
      PRINT *, 'N(18)  1  -1    0    |/               |/'
      PRINT *, 'N(19)  1   1    0    1--------9-------2'
      PRINT *, 'N(20) -1   1    0'
      PRINT *, ' '
      PRINT *, ' Normalized Integration Point coordinates (Cip)'
      PRINT *, '         ξ   η    ζ'
      PRINT *, 'Cip(1)  -1  -1   -1'
      PRINT *, 'Cip(2)   0  -1   -1'
      PRINT *, 'Cip(3)   1  -1   -1'
      PRINT *, 'Cip(4)  -1   0   -1                          \η'
      PRINT *, 'Cip(5)   0   0   -1                \ζ       /'
      PRINT *, 'Cip(6)   1   0   -1                |       /'
      PRINT *, 'Cip(7)  -1   1   -1         25-----|26-------27'
      PRINT *, 'Cip(8)   0   1   -1         /|     |     /   /|'
      PRINT *, 'Cip(9)   1   1   -1        / |     |    /   / |'
      PRINT *, 'Cip(10) -1  -1    0      22--|----23------24  |'
      PRINT *, 'Cip(11)  0  -1    0      /  16-----|17/---/---18'
      PRINT *, 'Cip(12)  1  -1    0     /   /|     | /   /   /|'
      PRINT *, 'Cip(13) -1   0    0   19-------20------21   / |'
      PRINT *, 'Cip(14)  0   0    0    |  13------14----|-15-----\ξ'
      PRINT *, 'Cip(15)  1   0    0    | /   7-------8--|-/---9'
      PRINT *, 'Cip(16) -1   1    0    |/   /           |/   /'
      PRINT *, 'Cip(17)  0   1    0   10--------11-----12   /'
      PRINT *, 'Cip(18)  1   1    0    |  4 -------5----|--6'
      PRINT *, 'Cip(19) -1  -1    1    | /              | /'
      PRINT *, 'Cip(20)  0  -1    1    |/               |/'
      PRINT *, 'Cip(21)  1  -1    1    1--------2-------3'
      PRINT *, 'Cip(22) -1   0    1   Multiplicative Factor:'
      PRINT *, 'Cip(23)  0   0    1      α = SQR(3/5)'
      PRINT *, 'Cip(24)  1   0    1   ALPHA=0.774596669241483'
      PRINT *, 'Cip(25) -1   1    1'
      PRINT *, 'Cip(26)  0   1    1'
      PRINT *, 'Cip(27)  1   1    1'
      PRINT *, ' '
      PRINT *, '----------------------------------------------------'
      PRINT *, ' '
      PRINT *, 'PROBLEM INFO'
      PRINT *, '------------'
      PRINT *, ' '
      PRINT *, 'Number of dimensions: ', NDIM
      PRINT *, 'Number of elements: ', numElem
      PRINT *, 'ID Offset original elements vs UEL: ', ElemOffset
      PRINT *, 'Number of integration points per element: ', numInt
      PRINT *, 'Total number of integration points: ', numElemInt
      PRINT *, 'Larget id number for the UEL elements: ', numIdElem
      IF (FREQfileSwitch .EQ. 1) THEN
            PRINT *, 'Frequency file: ', FREQfile
            PRINT *, 'Frequency file path: ', FREQfilePath
            PRINT *, ' '
            PRINT *, 'For PN eq. during the simulation:'
            PRINT *, 'Number of activities: ', NActivities
            PRINT *, 'STEP_LIST_STR:  ', STEP_LIST_STR
            PRINT *, 'HOUR_RANGE_STR: ', HOUR_RANGE_STR
            PRINT *, 'FREQ_VAL:       ', FREQ_VAL
            PRINT *, ' '
            PRINT *, 'For PN eq. STEP:'
            PRINT *, 'Number of activities: ', PN_NActivities
            PRINT *, 'Hours of simulation: ', PN_TIME_file
            PRINT *, 'Frequency per hour:'
            DO I = 1, PN_TIME_file
                  PRINT *, ' ', I, ': ', PN_FREQ_file(I)
            ENDDO
      ELSE
            PRINT *, 'Frequency file: ', 'None'
      END IF
      PRINT *, ' '
      PRINT *, 'Number of Step: ', NSTEPS
C NOW THE STEP NAMES STORED IN STEPNAMES(NSTEPS)
      DO I=1,NSTEPS
            PRINT '(A,I3,A)', ' Step(',I,'): ', STEPNAMES(I)
      ENDDO
      PRINT *, ' '
      PRINT *, '****************************************************'
      PRINT *, '   THE SDVs ARE STORED IN THE UVARM SUBROUTINE'
      PRINT *, '****************************************************'
      PRINT *, ' '
      DO I=1,ngSdv
            PRINT '(A,I3,A)', ' UVARM(',I,'):'
            PRINT *, UVARM_NAMES(I)
      ENDDO
      PRINT *, ' '
      PRINT *, '****************************************************'
      PRINT *, ' '

C--------------------------------------------------------------------------
      ELSEIF (LOP .EQ. 2) THEN  !end of the increment

C Call DATE_AND_TIME
C date_ini: Initial date and time values
C date_fin: Current date and time values
      CALL DATE_AND_TIME(values=date_fin)

C Calculate the difference between the initial and current date and time
      CALL diff_DATE_AND_TIME(date_ini, date_fin, date_diff)

C Opening file and write summary
      OPEN(15,FILE=infoFilePath_T_D, STATUS='OLD', POSITION='APPEND')

      IF (KINC .EQ. 1 .AND. KSTEP .EQ. 1) THEN
      PRINT *, ' '
      PRINT *, 'SUMMARY OF JOB INFORMATION:'
      PRINT *, TRIM('STEP  INC  TOTAL    TOTAL        STEP')
     1//TRIM('       INC OF     REAL TIME')
      PRINT *, TRIM('           ITERS    TIME/      TIME/LPF')
     1//TRIM('    TIME/LPF     D/H/M/S ')

      WRITE(15,*) ' '
      WRITE(15,*) 'SUMMARY OF JOB INFORMATION:'
      WRITE(15,*) TRIM('STEP  INC  TOTAL    TOTAL        STEP')
     1//TRIM('       INC OF     REAL TIME')
      WRITE(15,*) TRIM('           ITERS    TIME/      TIME/LPF')
     1//TRIM('    TIME/LPF     D/H/M/S ')

      END IF

      PRINT '(2I5, I6, 3E12.3, A, I3, A, I2, A, I2, A, I2, A)', 
     1KSTEP, KINC, ITERATION, TIME(2), TIME(1), DTIME, ' '
     2, date_diff(1), 'd', date_diff(2), 'h', date_diff(3), 'm'
     3, date_diff(4), 's'

      WRITE(15, '(2I5, I6, 3E12.3, A, I3, A, I2, A, I2, A, I2, A)') 
     1KSTEP, KINC, ITERATION, TIME(2), TIME(1), DTIME, ' '
     2, date_diff(1), 'd', date_diff(2), 'h', date_diff(3), 'm'
     3, date_diff(4), 's'
      CLOSE(15)

C--------------------------------------------------------------------------
      ELSEIF (LOP .EQ. 3) THEN  !end of the analisys

C Call DATE_AND_TIME
      CALL DATE_AND_TIME(values=date_fin)
C Calculate the difference between the initial and current date and time
      CALL diff_DATE_AND_TIME(date_ini, date_fin, date_diff)

      OPEN(15,FILE=infoFilePath_T_D, STATUS='OLD', POSITION='APPEND')
      WRITE(15,*) ' '
      WRITE(15,*) '****************************************************'
      WRITE(15,*) ' '
      WRITE(15,*) ' End of the püdü analisys'
      WRITE(15,*) ' '
      WRITE(15,'(A, I4.4, A, I2.2, A, I2.2, A)') ' Date: '
     1, date_fin(1), '-', date_fin(2), '-'
     2, date_fin(3)
      WRITE(15,'(A, I2.2, A, I2.2, A, I2.2, A, I3.3)') ' Time: '
     1, date_fin(5), ':', date_fin(6), ':'
     2, date_fin(7)
      WRITE(15,*) ' '
      WRITE(15,*) 'Total simulation time: ', TIME(2), 's'
      WRITE(15,*) 'Total real time: ', date_diff(1), 'd', date_diff(2)
     1, 'h', date_diff(3), 'm', date_diff(4), 's'
      WRITE(15,*) 'Total real seconds: ', date_diff(5)
      WRITE(15,*) ' '
      WRITE(15,*) ' 01001000 01100001 01111000'
      WRITE(15,*) ' 01001000 01100001 01111000 01111000'
      WRITE(15,*) ' 00110010 00110011'
      WRITE(15,*) ' '
      WRITE(15,*) ' Auth: Estefano Muñoz-Moya'
      WRITE(15,*) ' LinkTree: https://linktr.ee/estefano23'
      WRITE(15,*) ' Web Page: https://estefano23.github.io/'
      WRITE(15,*) ' GitHub: estefano23'
      WRITE(15,*) ' Email: estefano.munoz.moya@gmail.com'
      WRITE(15,*) ' HaxHaxx23'
      WRITE(15,*) ' '
      WRITE(15,*) '          :   .:.                                   '
      WRITE(15,*) '         ^~  :^:                                    '
      WRITE(15,*) '   :::. :~^::~..  ...                               '
      WRITE(15,*) '   ^::^^!^~^^~~~:~^.^.                              '
      WRITE(15,*) '   .~.^:..  .^~ ~~.:^                               '
      WRITE(15,*) '    :~.  .^: .:.^^~^                                '
      WRITE(15,*) '    :    !5?    .~                                  '
      WRITE(15,*) '  ~!:            .^            ......:....          '
      WRITE(15,*) ' .5P~             .:.    ..::::....  .  .....       '
      WRITE(15,*) '  .~^:::.....       ::::...                 .:.     '
      WRITE(15,*) '       ..~:..                                 .:.   '
      WRITE(15,*) '          ^.        ..                 ..       :^  '
      WRITE(15,*) '           ^       ..        .       .:          !  '
      WRITE(15,*) '           .^.                      .:          ^~  '
      WRITE(15,*) '             :.                     ..        .~~   '
      WRITE(15,*) '              ::.. .       .:      ..        :~:    '
      WRITE(15,*) '                :~::^.     .~...^.^^.        :~     '
      WRITE(15,*) '                 ^. .^^    ~:...:^..:::..     ::    '
      WRITE(15,*) '                 :  :^^   ^.      .:.  .~::.   ~:   '
      WRITE(15,*) '                :: .:::  :.         .:. .^..:: .^.  '
      WRITE(15,*) '                : ::  ^ ::            ~  ^   ~ .:   '
      WRITE(15,*) '               : :^   ^ ::           :: ^.  ~^.^    '
      WRITE(15,*) '              ~^:~   :  !.          :^.~:  ^^ !.    '
      WRITE(15,*) '             .?J?.  ^?~^:          .7J7.  :?J?:     '
      WRITE(15,*) '                    ^7!^                   .:.      '
      WRITE(15,*) ' '
      WRITE(15,*) '****************************************************'
      WRITE(15,*) ' '
      CLOSE(15)

      PRINT *, ' '
      PRINT *, '****************************************************'
      PRINT *, ' '
      PRINT *, ' End of the püdü analisys'
      PRINT *, ' '
      PRINT '(A, I4.4, A, I2.2, A, I2.2, A)', ' Date: '
     1, date_fin(1), '-', date_fin(2), '-'
     2, date_fin(3)
      PRINT '(A, I2.2, A, I2.2, A, I2.2, A, I3.3)', ' Time: '
     1, date_fin(5), ':', date_fin(6), ':'
     2, date_fin(7)
      PRINT *, ' '
      PRINT *, 'Total simulation time: ', TIME(2), 's'
      PRINT *, 'Total real time: ', date_diff(1), 'd', date_diff(2)
     1, 'h', date_diff(3), 'm', date_diff(4), 's'
      PRINT *, 'Total real seconds: ', date_diff(5)
      PRINT *, ' '
      PRINT *, ' 01001000 01100001 01111000'
      PRINT *, ' 01001000 01100001 01111000 01111000'
      PRINT *, ' 00110010 00110011'
      PRINT *, ' '
      PRINT *, ' Auth: Estefano Muñoz-Moya'
      PRINT *, ' LinkTree: https://linktr.ee/estefano23'
      PRINT *, ' Web Page: https://estefano23.github.io/'
      PRINT *, ' GitHub: estefano23'
      PRINT *, ' Email: estefano.munoz.moya@gmail.com'
      PRINT *, ' HaxHaxx23'
      PRINT *, ' '
      PRINT *, '          :   .:.                                   '
      PRINT *, '         ^~  :^:                                    '
      PRINT *, '   :::. :~^::~..  ...                               '
      PRINT *, '   ^::^^!^~^^~~~:~^.^.                              '
      PRINT *, '   .~.^:..  .^~ ~~.:^                               '
      PRINT *, '    :~.  .^: .:.^^~^                                '
      PRINT *, '    :    !5?    .~                                  '
      PRINT *, '  ~!:            .^            ......:....          '
      PRINT *, ' .5P~             .:.    ..::::....  .  .....       '
      PRINT *, '  .~^:::.....       ::::...                 .:.     '
      PRINT *, '       ..~:..                                 .:.   '
      PRINT *, '          ^.        ..                 ..       :^  '
      PRINT *, '           ^       ..        .       .:          !  '
      PRINT *, '           .^.                      .:          ^~  '
      PRINT *, '             :.                     ..        .~~   '
      PRINT *, '              ::.. .       .:      ..        :~:    '
      PRINT *, '                :~::^.     .~...^.^^.        :~     '
      PRINT *, '                 ^. .^^    ~:...:^..:::..     ::    '
      PRINT *, '                 :  :^^   ^.      .:.  .~::.   ~:   '
      PRINT *, '                :: .:::  :.         .:. .^..:: .^.  '
      PRINT *, '                : ::  ^ ::            ~  ^   ~ .:   '
      PRINT *, '               : :^   ^ ::           :: ^.  ~^.^    '
      PRINT *, '              ~^:~   :  !.          :^.~:  ^^ !.    '
      PRINT *, '             .?J?.  ^?~^:          .7J7.  :?J?:     '
      PRINT *, '                    ^7!^                   .:.      '
      PRINT *, ' '
      PRINT *, '****************************************************'
      PRINT *, ' '
      END IF
