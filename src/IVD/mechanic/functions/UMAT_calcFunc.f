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
C ******    ABAQUS UMAT: Useful functions                              ******
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

C Calculate the identity matrix
C -----------------------------

      SUBROUTINE onem(A)
C This subroutine stores the identity matrix in the
C 3 by 3 matrix [A]
      IMPLICIT NONE

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

      RETURN
      END SUBROUTINE onem

C----------------------------------------------------------------------------------------------------------------------------------------------

C Calculate the identity vector
C -----------------------------

      SUBROUTINE onev(A)
C This subroutine stores the identity vector in the
C 3 by 1 vector [A]
      IMPLICIT NONE

      INTEGER I

      DOUBLE PRECISION A(3,1)

C some constant parameters
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.D0, ONE=1.0D0
      DO I=1,3
            A(I,1) = ONE
      END DO

      RETURN
      END SUBROUTINE onev

C----------------------------------------------------------------------------------------------------------------------------------------------

C Calculate the determinant of the deformation gradient
C -----------------------------------------------------------

      SUBROUTINE mDETDFG(DFGRD1,DETDFG,NSHR)
C This subroutine calculates the determinant of a 3 by 3 matrix [A]
      IMPLICIT NONE

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

      RETURN
      END SUBROUTINE mDETDFG

C----------------------------------------------------------------------------------------------------------------------------------------------

C Calculate the determinant of a given 3x3 matrix
C -----------------------------------------------------------

      SUBROUTINE mDet(A,det_A)
C This subroutine calculates the determinant of a 3 by 3 matrix [A]
      IMPLICIT NONE

      DOUBLE PRECISION  A(3,3),det_A

      det_A = A(1,1)*A(2,2)*A(3,3)
     +  + A(1,2)*A(2,3)*A(3,1)
     +  + A(1,3)*A(2,1)*A(3,2)
     +  - A(3,1)*A(2,2)*A(1,3)
     +  - A(3,2)*A(2,3)*A(1,1)
     +  - A(3,3)*A(2,1)*A(1,2)

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

      INTEGER istat

      DOUBLE PRECISION A(3,3), A_inv(3,3), A_invT(3,3)
      DOUBLE PRECISION det_A, det_A_inv


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

      RETURN
      END SUBROUTINE matInvDet3D

C----------------------------------------------------------------------------------------------------------------------------------------------

C The D/H/M/S between two dates
C -----------------------------------------------------------

      SUBROUTINE diff_DATE_AND_TIME(date_ini, date_fin, date_diff)
C This subroutine calculates the difference between two dates and times
      IMPLICIT NONE

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


      RETURN
      END SUBROUTINE diff_DATE_AND_TIME
