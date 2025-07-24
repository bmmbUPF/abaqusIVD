      IF (i_numElem .GT. numElem) THEN  !If the element number is greater than the total number of elements
      
C Write the error message to the output file
      OPEN(15,FILE=infoFilePath_T_D, STATUS='OLD', POSITION='APPEND')
      WRITE(15,*) ' '
      WRITE(15,*) '////////////////////////////////////////////////////'
      WRITE(15,*) '     ERROR: ELEMENT NUMBER IS GREATER THAN THE      '
      WRITE(15,*) '     TOTAL NUMBER OF ELEMENTS                       '
      WRITE(15,*) 'ACTUAL ELEMENT NUMBER:    ', i_numElem
      WRITE(15,*) 'TOTAL NUMBER OF ELEMENTS: ', numElem
      WRITE(15,*) '////////////////////////////////////////////////////'
      WRITE(15,*) ' '

C Print the error message to the console
      PRINT *, ' '
      PRINT *, '////////////////////////////////////////////////////'
      PRINT *, '     ERROR: ELEMENT NUMBER IS GREATER THAN THE      '
      PRINT *, '     TOTAL NUMBER OF ELEMENTS                       '
      PRINT *, 'ACTUAL ELEMENT NUMBER:    ', i_numElem
      PRINT *, 'TOTAL NUMBER OF ELEMENTS: ', numElem
      PRINT *, '////////////////////////////////////////////////////'
      PRINT *, ' '

C Terminate the analysis
      CALL XIT
      END IF