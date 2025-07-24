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
C ******    - Integration point locations for a 3D element              ******
C ******    - Shape functions and derivatives for a 3D element          ******
C ******    - Mapping of shape functions from the isoparametric domain  ******
C ******      to the physical domain                                    ******
C ******    - Time integration of the problem                           ******
C ******    - Assemble the global stiffness and force matrices          ******
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

C This subroutine will get the integration point locations for a 3D element
C Hexahedral element with 20 nodes (C3D20) using 27 gauss points for integration

C Normalized Nodal Coordinates for the shape functions (N)
C        x   y    z
C N(1)  -1  -1   -1
C N(2)   1  -1   -1                                  \y
C N(3)   1   1   -1                        \z       /
C N(4)  -1   1   -1                        |       /
C N(5)  -1  -1    1                  8-----|15--------7
C N(6)   1  -1    1                 /|     |     /   /|
C N(7)   1   1    1                / |     |    /   / |
C N(8)  -1   1    1              16  |     |   /  14  |
C N(9)   0  -1   -1              /  20     |  /   /   19
C N(10)  1   0   -1             /    |     | /   /    |
C N(11)  0   1   -1            5-------13-------6     |
C N(12) -1   0   -1            |     |     +----|------------\x
C N(13)  0  -1    1            |     4-------11-|-----3
C N(14)  1   0    1            |    /           |    /
C N(15)  0   1    1           17   /           18   /
C N(16) -1   0    1            |  12            |  10
C N(17) -1  -1    0            | /              | /
C N(18)  1  -1    0            |/               |/
C N(19)  1   1    0            1--------9-------2
C N(20) -1   1    0
C

C Normalized Integration Point Coordinates (Cip)
C          ξ   η    ζ
C Cip(1)  -1  -1   -1
C Cip(2)   0  -1   -1
C Cip(3)   1  -1   -1
C Cip(4)  -1   0   -1                                  \η
C Cip(5)   0   0   -1                        \ζ       /
C Cip(6)   1   0   -1                        |       /
C Cip(7)  -1   1   -1                 25-----|26-------27
C Cip(8)   0   1   -1                 /|     |     /   /|
C Cip(9)   1   1   -1                / |     |    /   / |
C Cip(10) -1  -1    0              22--|----23------24  |
C Cip(11)  0  -1    0              /  16-----|17/---/---18
C Cip(12)  1  -1    0             /   /|     | /   /   /|
C Cip(13) -1   0    0           19-------20------21   / |
C Cip(14)  0   0    0            |  13------14----|-15---------\ξ
C Cip(15)  1   0    0            | /   7-------8--|-/---9
C Cip(16) -1   1    0            |/   /           |/   /
C Cip(17)  0   1    0           10--------11-----12   /
C Cip(18)  1   1    0            |  4 -------5----|--6
C Cip(19) -1  -1    1            | /              | /
C Cip(20)  0  -1    1            |/               |/
C Cip(21)  1  -1    1            1--------2-------3
C Cip(22) -1   0    1           Multiplicative Factor:
C Cip(23)  0   0    1              α = SQR(3/5)
C Cip(24)  1   0    1          ALPHA=0.774596669241483
C Cip(25) -1   1    1
C Cip(26)  0   1    1
C Cip(27)  1   1    1
C


! 7------8------9
! |             |
! |             |
! 4------5------6
! |             |
! |             |
! 1------2------3

! 16-----17-----18
! |             |
! |             |
! 13-----14-----15
! |             |
! |             |
! 10-----11-----12

! 25-----26-----27
! |             |
! |             |
! 22-----23-----24
! |             |
! |             |
! 19-----20-----21

C
C-------------------------------------------------------------------------------------------------------------------------------------------------------------------
C Code
C ------------

C Integration points and weights
C ------------------------------

      SUBROUTINE CintPtC3D20pt(Cip,whtG,i_intPT,NDIM)
      
C This subroutine will get the integration point locations
C and corresponding gauss quadrature weights for 3D elements
C using 27 gauss points for integration
C
C Cip(i_intPT,3): Two-dimensional array representing the coordinates of the Gauss integration points
C             in a reference or "master" 3D element. Each row of Cip corresponds to an integration point,
C             and each column corresponds to a coordinate direction: ξ (xi), η (eta), and ζ (zeta).
C             These are the natural (or isoparametric) coordinates,
C             which are used to define positions within the reference element.

C whtG(i_intPT): The corresponding Gauss integration weights

      IMPLICIT NONE

      INTEGER i_intPT,NDIM

      DOUBLE PRECISION Cip(i_intPT,NDIM),whtG(i_intPT)

C CONSTANT PARAMETERS
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.D0, ONE=1.0D0, TWO=2.0D0
      DOUBLE PRECISION, PARAMETER :: HALF=0.5D0, FOURTH=0.25D0

C ALPHA (SQR(3/5)), C1 (5/9), AND C2 (8/9) ARE THE GAUSS INTEGRATION POINTS
C      DOUBLE PRECISION, PARAMETER :: ALPHA=0.774596669241483D0
      DOUBLE PRECISION, PARAMETER :: ALPHA=DSQRT(3.D0/5.D0)
      DOUBLE PRECISION, PARAMETER :: C1=125.D0
      DOUBLE PRECISION, PARAMETER :: C2=200.D0
      DOUBLE PRECISION, PARAMETER :: C3=320.D0
      DOUBLE PRECISION, PARAMETER :: C4=512.D0
      DOUBLE PRECISION, PARAMETER :: D1=729.D0

C Initialize
      whtG = ZERO
      Cip = ZERO

C--------------------------------------------------------------------------
C Gauss weights
      whtG(1) = C1/D1 !0.171467764060356
      whtG(2) = C2/D1
      whtG(3) = C1/D1
      whtG(4) = C2/D1
      whtG(5) = C3/D1
      whtG(6) = C2/D1
      whtG(7) = C1/D1
      whtG(8) = C2/D1
      whtG(9) = C1/D1
      whtG(10) = C2/D1
      whtG(11) = C3/D1
      whtG(12) = C2/D1
      whtG(13) = C3/D1
      whtG(14) = C4/D1
      whtG(15) = C3/D1
      whtG(16) = C2/D1
      whtG(17) = C3/D1
      whtG(18) = C2/D1
      whtG(19) = C1/D1
      whtG(20) = C2/D1
      whtG(21) = C1/D1
      whtG(22) = C2/D1
      whtG(23) = C3/D1
      whtG(24) = C2/D1
      whtG(25) = C1/D1
      whtG(26) = C2/D1
      whtG(27) = C1/D1

C--------------------------------------------------------------------------
C Gauss pt locations in master element
C Cip(integrationPoint,Coord(ξ,η,ζ))
      Cip(1,1) = -ALPHA
      Cip(1,2) = -ALPHA
      Cip(1,3) = -ALPHA

      Cip(2,1) = ZERO
      Cip(2,2) = -ALPHA
      Cip(2,3) = -ALPHA

      Cip(3,1) = ALPHA
      Cip(3,2) = -ALPHA
      Cip(3,3) = -ALPHA

      Cip(4,1) = -ALPHA
      Cip(4,2) = ZERO
      Cip(4,3) = -ALPHA

      Cip(5,1) = ZERO
      Cip(5,2) = ZERO
      Cip(5,3) = -ALPHA

      Cip(6,1) = ALPHA
      Cip(6,2) = ZERO
      Cip(6,3) = -ALPHA

      Cip(7,1) = -ALPHA
      Cip(7,2) = ALPHA
      Cip(7,3) = -ALPHA

      Cip(8,1) = ZERO
      Cip(8,2) = ALPHA
      Cip(8,3) = -ALPHA

      Cip(9,1) = ALPHA
      Cip(9,2) = ALPHA
      Cip(9,3) = -ALPHA

      Cip(10,1) = -ALPHA
      Cip(10,2) = -ALPHA
      Cip(10,3) = ZERO

      Cip(11,1) = ZERO
      Cip(11,2) = -ALPHA
      Cip(11,3) = ZERO

      Cip(12,1) = ALPHA
      Cip(12,2) = -ALPHA
      Cip(12,3) = ZERO

      Cip(13,1) = -ALPHA
      Cip(13,2) = ZERO
      Cip(13,3) = ZERO

      Cip(14,1) = ZERO
      Cip(14,2) = ZERO
      Cip(14,3) = ZERO

      Cip(15,1) = ALPHA
      Cip(15,2) = ZERO
      Cip(15,3) = ZERO

      Cip(16,1) = -ALPHA
      Cip(16,2) = ALPHA
      Cip(16,3) = ZERO

      Cip(17,1) = ZERO
      Cip(17,2) = ALPHA
      Cip(17,3) = ZERO

      Cip(18,1) = ALPHA
      Cip(18,2) = ALPHA
      Cip(18,3) = ZERO

      Cip(19,1) = -ALPHA
      Cip(19,2) = -ALPHA
      Cip(19,3) = ALPHA

      Cip(20,1) = ZERO
      Cip(20,2) = -ALPHA
      Cip(20,3) = ALPHA

      Cip(21,1) = ALPHA
      Cip(21,2) = -ALPHA
      Cip(21,3) = ALPHA

      Cip(22,1) = -ALPHA
      Cip(22,2) = ZERO
      Cip(22,3) = ALPHA

      Cip(23,1) = ZERO
      Cip(23,2) = ZERO
      Cip(23,3) = ALPHA

      Cip(24,1) = ALPHA
      Cip(24,2) = ZERO
      Cip(24,3) = ALPHA

      Cip(25,1) = -ALPHA
      Cip(25,2) = ALPHA
      Cip(25,3) = ALPHA

      Cip(26,1) = ZERO
      Cip(26,2) = ALPHA
      Cip(26,3) = ALPHA

      Cip(27,1) = ALPHA
      Cip(27,2) = ALPHA
      Cip(27,3) = ALPHA

C--------------------------------------------------------------------------
      RETURN
      END SUBROUTINE CintPtC3D20pt

C----------------------------------------------------------------------------------------------------------------------------------------------
C Shape functions and derivatives
C -------------------------------

      SUBROUTINE shapeFunc_C3D20(i_intPT,Cip,intpt,N,dN,d2N,NNODE,NDIM)

C Define the shape functions for C3D20 Element
      
      IMPLICIT NONE

      INTEGER i_intPT,intpt,NDIM,NNODE

      DOUBLE PRECISION Cip(i_intPT,NDIM),N(NNODE)
      DOUBLE PRECISION dN(NNODE,NDIM), d2N(NNODE,NDIM,NDIM)
      DOUBLE PRECISION xi,eta,zeta
      DOUBLE PRECISION omg,omh,omr,opg,oph,opr, tpgphpr,tmgphpr,tmgmhpr
      DOUBLE PRECISION tpgmhpr,tpgphmr,tmgphmr,tmgmhmr,tpgmhmr, omgopg
      DOUBLE PRECISION omhoph,omropr,omgmopg,omhmoph,omrmopr


C some constant parameters
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.D0, ONE=1.0D0, TWO=2.0D0
      DOUBLE PRECISION, PARAMETER :: HALF=0.5D0, FOURTH=0.25D0
      DOUBLE PRECISION, PARAMETER :: EIGHTH=0.125D0

C Location in the master element
      xi = Cip(intpt,1)
      eta = Cip(intpt,2)
      zeta = Cip(intpt,3)

C Defining initial values to N and dN
      N = ZERO
      dN = ZERO
      d2N = ZERO

C Helping variables
      omg = ONE-xi
      omh = ONE-eta
      omr = ONE-zeta
      opg = ONE+xi
      oph = ONE+eta
      opr = ONE+zeta

      tpgphpr = opg+oph+zeta
      tmgphpr = omg+oph+zeta
      tmgmhpr = omg+omh+zeta
      tpgmhpr = opg+omh+zeta

      tpgphmr = opg+oph-zeta
      tmgphmr = omg+oph-zeta
      tmgmhmr = omg+omh-zeta
      tpgmhmr = opg+omh-zeta

      omgopg = omg*opg*FOURTH
      omhoph = omh*oph*FOURTH
      omropr = omr*opr*FOURTH

      omgmopg = (omg-opg)*FOURTH
      omhmoph = (omh-oph)*FOURTH
      omrmopr = (omr-opr)*FOURTH

C--------------------------------------------------------------------------
C Fill in N(1:20) for all 20 nodes based on their location

C Nodes 1 - 8: These are the corner nodes of the hexahedral element.
C Their shape functions are derived from the trilinear interpolation within the element.
C They are associated with the vertices of the brick.
      N(1) = -omg*omh*omr*tpgphpr*EIGHTH
      N(2) = -opg*omh*omr*tmgphpr*EIGHTH
      N(3) = -opg*oph*omr*tmgmhpr*EIGHTH
      N(4) = -omg*oph*omr*tpgmhpr*EIGHTH
      N(5) = -omg*omh*opr*tpgphmr*EIGHTH
      N(6) = -opg*omh*opr*tmgphmr*EIGHTH
      N(7) = -opg*oph*opr*tmgmhmr*EIGHTH
      N(8) = -omg*oph*opr*tpgmhmr*EIGHTH

C Nodes 9 - 16: These are the midside nodes located on the edges of the hexahedron
C in the xi, eta, and zeta directions, but not at the corners. 
      N(9) = omgopg*omh*omr
      N(10) = omhoph*opg*omr
      N(11) = omgopg*oph*omr
      N(12) = omhoph*omg*omr
      N(13) = omgopg*omh*opr
      N(14) = omhoph*opg*opr
      N(15) = omgopg*oph*opr
      N(16) = omhoph*omg*opr


C Nodes 17 - 20: These nodes are located on the center of the faces of the vertical edges.
      N(17) = omropr*omg*omh
      N(18) = omropr*opg*omh
      N(19) = omropr*opg*oph
      N(20) = omropr*omg*oph

C--------------------------------------------------------------------------
C calculating the derivatives of the shape functions with respect to xi, eta, and zeta

C xi derivatives
      dN(1,1) = omh*omr*(tpgphpr-omg)*EIGHTH
      dN(2,1) = (opg-tmgphpr)*omh*omr*EIGHTH
      dN(3,1) = (opg-tmgmhpr)*oph*omr*EIGHTH
      dN(4,1) = oph*omr*(tpgmhpr-omg)*EIGHTH
      dN(5,1) = omh*opr*(tpgphmr-omg)*EIGHTH
      dN(6,1) = (opg-tmgphmr)*omh*opr*EIGHTH
      dN(7,1) = (opg-tmgmhmr)*oph*opr*EIGHTH
      dN(8,1) = oph*opr*(tpgmhmr-omg)*EIGHTH
      dN(9,1) = omgmopg*omh*omr
      dN(10,1) = omhoph*omr
      dN(11,1) = omgmopg*oph*omr
      dN(12,1) = -omhoph*omr
      dN(13,1) = omgmopg*omh*opr
      dN(14,1) = omhoph*opr
      dN(15,1) = omgmopg*oph*opr
      dN(16,1) = -omhoph*opr
      dN(17,1) = -omropr*omh
      dN(18,1) = omropr*omh
      dN(19,1) = omropr*oph
      dN(20,1) = -omropr*oph

C eta derivatives
      dN(1,2) = omg*omr*(tpgphpr-omh)*EIGHTH
      dN(2,2) = opg*omr*(tmgphpr-omh)*EIGHTH
      dN(3,2) = opg*(oph-tmgmhpr)*omr*EIGHTH
      dN(4,2) = omg*(oph-tpgmhpr)*omr*EIGHTH
      dN(5,2) = omg*opr*(tpgphmr-omh)*EIGHTH
      dN(6,2) = opg*opr*(tmgphmr-omh)*EIGHTH
      dN(7,2) = opg*(oph-tmgmhmr)*opr*EIGHTH
      dN(8,2) = omg*(oph-tpgmhmr)*opr*EIGHTH
      dN(9,2) = -omgopg*omr
      dN(10,2) = omhmoph*opg*omr
      dN(11,2) = omgopg*omr
      dN(12,2) = omhmoph*omg*omr
      dN(13,2) = -omgopg*opr
      dN(14,2) = omhmoph*opg*opr
      dN(15,2) = omgopg*opr
      dN(16,2) = omhmoph*omg*opr
      dN(17,2) = -omropr*omg
      dN(18,2) = -omropr*opg
      dN(19,2) = omropr*opg
      dN(20,2) = omropr*omg

C zeta derivatives
      dN(1,3) = omg*omh*(tpgphpr-omr)*EIGHTH
      dN(2,3) = opg*omh*(tmgphpr-omr)*EIGHTH
      dN(3,3) = opg*oph*(tmgmhpr-omr)*EIGHTH
      dN(4,3) = omg*oph*(tpgmhpr-omr)*EIGHTH
      dN(5,3) = omg*omh*(opr-tpgphmr)*EIGHTH
      dN(6,3) = opg*omh*(opr-tmgphmr)*EIGHTH
      dN(7,3) = opg*oph*(opr-tmgmhmr)*EIGHTH
      dN(8,3) = omg*oph*(opr-tpgmhmr)*EIGHTH
      dN(9,3) = -omgopg*omh
      dN(10,3) = -omhoph*opg
      dN(11,3) = -omgopg*oph
      dN(12,3) = -omhoph*omg
      dN(13,3) = omgopg*omh
      dN(14,3) = omhoph*opg
      dN(15,3) = omgopg*oph
      dN(16,3) = omhoph*omg
      dN(17,3) = omrmopr*omg*omh
      dN(18,3) = omrmopr*opg*omh
      dN(19,3) = omrmopr*opg*oph
      dN(20,3) = omrmopr*omg*oph

C--------------------------------------------------------------------------
C calculating the second derivatives of the shape functions with respect to xi, eta, and zeta

C second derivatives wrt xi^2
      d2N(1,1,1) = (eta - ONE)*(zeta - ONE)*FOURTH
      d2N(2,1,1) = d2N(1,1,1)
      d2N(3,1,1) = -(eta + ONE)*(zeta - ONE)*FOURTH
      d2N(4,1,1) = d2N(3,1,1)
      d2N(5,1,1) = -(eta - ONE)*(zeta + ONE)*FOURTH
      d2N(6,1,1) = d2N(5,1,1)
      d2N(7,1,1) = (eta + ONE)*(zeta + ONE)*FOURTH
      d2N(8,1,1) = d2N(7,1,1)
      d2N(9,1,1) = -(eta - ONE)*(zeta - ONE)/TWO
!       d2N(10,1,1) = ZERO
      d2N(11,1,1) = (eta + ONE)*(zeta - ONE)/TWO
!       d2N(12,1,1) = ZERO
      d2N(13,1,1) = (eta - ONE)*(zeta + ONE)/TWO
!       d2N(14,1,1) = ZERO
      d2N(15,1,1) = -(eta + ONE)*(zeta + ONE)/TWO
!       d2N(16,1,1) = ZERO
!       d2N(17,1,1) = ZERO
!       d2N(18,1,1) = ZERO
!       d2N(19,1,1) = ZERO
!       d2N(20,1,1) = ZERO

C second derivatives wrt eta^2
      d2N(1,2,2) = (xi - ONE)*(zeta - ONE)*FOURTH
      d2N(2,2,2) = -(xi + ONE)*(zeta - ONE)*FOURTH
      d2N(3,2,2) = d2N(2,2,2)
      d2N(4,2,2) = d2N(1,2,2)
      d2N(5,2,2) = -(xi - ONE)*(zeta + ONE)*FOURTH
      d2N(6,2,2) = (xi + ONE)*(zeta + ONE)*FOURTH
      d2N(7,2,2) = d2N(6,2,2)
      d2N(8,2,2) = d2N(5,2,2)
!       d2N(9,2,2) = ZERO
      d2N(10,2,2) = (xi + ONE)*(zeta - ONE)/TWO
!       d2N(11,2,2) = ZERO
      d2N(12,2,2) = -(xi - ONE)*(zeta - ONE)/TWO
!       d2N(13,2,2) = ZERO
      d2N(14,2,2) = -(xi + ONE)*(zeta + ONE)/TWO
!       d2N(15,2,2) = ZERO
      d2N(16,2,2) = (xi - ONE)*(zeta + ONE)/TWO
!       d2N(17,2,2) = ZERO
!       d2N(18,2,2) = ZERO
!       d2N(19,2,2) = ZERO
!       d2N(20,2,2) = ZERO

C second derivatives wrt zeta^2
      d2N(1,3,3) = (eta - ONE)*(xi - ONE)*FOURTH
      d2N(2,3,3) = -(eta - ONE)*(xi + ONE)*FOURTH
      d2N(3,3,3) = (eta + ONE)*(xi + ONE)*FOURTH
      d2N(4,3,3) = -(eta + ONE)*(xi - ONE)*FOURTH
      d2N(5,3,3) = d2N(1,3,3)
      d2N(6,3,3) = d2N(2,3,3)
      d2N(7,3,3) = d2N(3,3,3)
      d2N(8,3,3) = d2N(4,3,3)
!       d2N(9,3,3) = ZERO
!       d2N(10,3,3) = ZERO
!       d2N(11,3,3) = ZERO
!       d2N(12,3,3) = ZERO
!       d2N(13,3,3) = ZERO
!       d2N(14,3,3) = ZERO
!       d2N(15,3,3) = ZERO
!       d2N(16,3,3) = ZERO
      d2N(17,3,3) = -(eta - ONE)*(xi - ONE)/TWO
      d2N(18,3,3) = (eta - ONE)*(xi + ONE)/TWO
      d2N(19,3,3) = -(eta + ONE)*(xi + ONE)/TWO
      d2N(20,3,3) = (eta + ONE)*(xi - ONE)/TWO

C second derivatives wrt xi,eta
      d2N(1,1,2) = (zeta - ONE)*(TWO*eta + TWO*xi + zeta)*EIGHTH
      d2N(2,1,2) = -(zeta - ONE)*(TWO*eta - TWO*xi + zeta)*EIGHTH
      d2N(3,1,2) = -(zeta - ONE)*(TWO*eta + TWO*xi - zeta)*EIGHTH
      d2N(4,1,2) = (zeta - ONE)*(TWO*eta - TWO*xi - zeta)*EIGHTH
      d2N(5,1,2) = -(zeta + ONE)*(TWO*eta + TWO*xi - zeta)*EIGHTH
      d2N(6,1,2) = (zeta + ONE)*(TWO*eta - TWO*xi - zeta)*EIGHTH
      d2N(7,1,2) = (zeta + ONE)*(TWO*eta + TWO*xi + zeta)*EIGHTH
      d2N(8,1,2) = -(zeta + ONE)*(TWO*eta - TWO*xi + zeta)*EIGHTH
      d2N(9,1,2) = -xi*(zeta - ONE)/TWO
      d2N(10,1,2) = eta*(zeta - ONE)/TWO
      d2N(11,1,2) = -d2N(9,1,2)
      d2N(12,1,2) = -d2N(10,1,2)
      d2N(13,1,2) = xi*(zeta + ONE)/TWO
      d2N(14,1,2) = -eta*(zeta + ONE)/TWO
      d2N(15,1,2) = -d2N(13,1,2)
      d2N(16,1,2) = -d2N(14,1,2)
      d2N(17,1,2) = -(zeta - ONE)*(zeta + ONE)*FOURTH
      d2N(18,1,2) = -d2N(17,1,2)
      d2N(19,1,2) = d2N(17,1,2)
      d2N(20,1,2) = -d2N(17,1,2)

C second derivatives wrt eta,xi
      d2N(1,2,1) = d2N(1,1,2)
      d2N(2,2,1) = d2N(2,1,2)
      d2N(3,2,1) = d2N(3,1,2)
      d2N(4,2,1) = d2N(4,1,2)
      d2N(5,2,1) = d2N(5,1,2)
      d2N(6,2,1) = d2N(6,1,2)
      d2N(7,2,1) = d2N(7,1,2)
      d2N(8,2,1) = d2N(8,1,2)
      d2N(9,2,1) = d2N(9,1,2)
      d2N(10,2,1) = d2N(10,1,2)
      d2N(11,2,1) = d2N(11,1,2)
      d2N(12,2,1) = d2N(12,1,2)
      d2N(13,2,1) = d2N(13,1,2)
      d2N(14,2,1) = d2N(14,1,2)
      d2N(15,2,1) = d2N(15,1,2)
      d2N(16,2,1) = d2N(16,1,2)
      d2N(17,2,1) = d2N(17,1,2)
      d2N(18,2,1) = d2N(18,1,2)
      d2N(19,2,1) = d2N(19,1,2)
      d2N(20,2,1) = d2N(20,1,2)
      
C second derivatives wrt xi,zeta
      d2N(1,1,3) = (eta - ONE)*(eta + TWO*xi + TWO*zeta)*EIGHTH
      d2N(2,1,3) = -(eta - ONE)*(eta - TWO*xi + TWO*zeta)*EIGHTH
      d2N(3,1,3) = -(eta + ONE)*(eta + TWO*xi - TWO*zeta)*EIGHTH
      d2N(4,1,3) = (eta + ONE)*(eta - TWO*xi - TWO*zeta)*EIGHTH
      d2N(5,1,3) = -(eta - ONE)*(eta + TWO*xi - TWO*zeta)*EIGHTH
      d2N(6,1,3) = (eta - ONE)*(eta - TWO*xi - TWO*zeta)*EIGHTH
      d2N(7,1,3) = (eta + ONE)*(eta + TWO*xi + TWO*zeta)*EIGHTH
      d2N(8,1,3) = -(eta + ONE)*(eta - TWO*xi + TWO*zeta)*EIGHTH
      d2N(9,1,3) = -xi*(eta - ONE)/TWO
      d2N(10,1,3) = (eta - ONE)*(eta + ONE)*FOURTH
      d2N(11,1,3) = xi*(eta + ONE)/TWO
      d2N(12,1,3) = -d2N(10,1,3)
      d2N(13,1,3) = -d2N(9,1,3)
      d2N(14,1,3) = -d2N(10,1,3)
      d2N(15,1,3) = -d2N(11,1,3)
      d2N(16,1,3) = d2N(10,1,3)
      d2N(17,1,3) = -zeta*(eta - ONE)/TWO
      d2N(18,1,3) = -d2N(17,1,3)
      d2N(19,1,3) = -zeta*(eta + ONE)/TWO
      d2N(20,1,3) = -d2N(19,1,3)

C second derivatives wrt zeta,xi
      d2N(1,3,1) = d2N(1,1,3)
      d2N(2,3,1) = d2N(2,1,3)
      d2N(3,3,1) = d2N(3,1,3)
      d2N(4,3,1) = d2N(4,1,3)
      d2N(5,3,1) = d2N(5,1,3)
      d2N(6,3,1) = d2N(6,1,3)
      d2N(7,3,1) = d2N(7,1,3)
      d2N(8,3,1) = d2N(8,1,3)
      d2N(9,3,1) = d2N(9,1,3)
      d2N(10,3,1) = d2N(10,1,3)
      d2N(11,3,1) = d2N(11,1,3)
      d2N(12,3,1) = d2N(12,1,3)
      d2N(13,3,1) = d2N(13,1,3)
      d2N(14,3,1) = d2N(14,1,3)
      d2N(15,3,1) = d2N(15,1,3)
      d2N(16,3,1) = d2N(16,1,3)
      d2N(17,3,1) = d2N(17,1,3)
      d2N(18,3,1) = d2N(18,1,3)
      d2N(19,3,1) = d2N(19,1,3)
      d2N(20,3,1) = d2N(20,1,3)

C second derivatives wrt eta,zeta
      d2N(1,2,3) = (xi - ONE)*(TWO*eta + xi + TWO*zeta)*EIGHTH
      d2N(2,2,3) = -(xi + ONE)*(TWO*eta - xi + TWO*zeta)*EIGHTH
      d2N(3,2,3) = -(xi + ONE)*(TWO*eta + xi - TWO*zeta)*EIGHTH
      d2N(4,2,3) = (xi - ONE)*(TWO*eta - xi - TWO*zeta)*EIGHTH
      d2N(5,2,3) = -(xi - ONE)*(TWO*eta + xi - TWO*zeta)*EIGHTH
      d2N(6,2,3) = (xi + ONE)*(TWO*eta - xi - TWO*zeta)*EIGHTH
      d2N(7,2,3) = (xi + ONE)*(TWO*eta + xi + TWO*zeta)*EIGHTH
      d2N(8,2,3) = -(xi - ONE)*(TWO*eta - xi + TWO*zeta)*EIGHTH
      d2N(9,2,3) = -(xi - ONE)*(xi + ONE)*FOURTH
      d2N(10,2,3) = eta*(xi + ONE)/TWO
      d2N(11,2,3) = -d2N(9,2,3)
      d2N(12,2,3) = -eta*(xi - ONE)/TWO
      d2N(13,2,3) = -d2N(9,2,3)
      d2N(14,2,3) = -d2N(10,2,3)
      d2N(15,2,3) = d2N(9,2,3)
      d2N(16,2,3) = -d2N(12,2,3)
      d2N(17,2,3) = -zeta*(xi - ONE)/TWO
      d2N(18,2,3) = zeta*(xi + ONE)/TWO
      d2N(19,2,3) = -d2N(18,2,3)
      d2N(20,2,3) = -d2N(17,2,3)

C second derivatives wrt zeta,eta
      d2N(1,3,2) = d2N(1,2,3)
      d2N(2,3,2) = d2N(2,2,3)
      d2N(3,3,2) = d2N(3,2,3)
      d2N(4,3,2) = d2N(4,2,3)
      d2N(5,3,2) = d2N(5,2,3)
      d2N(6,3,2) = d2N(6,2,3)
      d2N(7,3,2) = d2N(7,2,3)
      d2N(8,3,2) = d2N(8,2,3)
      d2N(9,3,2) = d2N(9,2,3)
      d2N(10,3,2) = d2N(10,2,3)
      d2N(11,3,2) = d2N(11,2,3)
      d2N(12,3,2) = d2N(12,2,3)
      d2N(13,3,2) = d2N(13,2,3)
      d2N(14,3,2) = d2N(14,2,3)
      d2N(15,3,2) = d2N(15,2,3)
      d2N(16,3,2) = d2N(16,2,3)
      d2N(17,3,2) = d2N(17,2,3)
      d2N(18,3,2) = d2N(18,2,3)
      d2N(19,3,2) = d2N(19,2,3)
      d2N(20,3,2) = d2N(20,2,3)
      
C--------------------------------------------------------------------------
      RETURN
      END SUBROUTINE shapeFunc_C3D20

C----------------------------------------------------------------------------------------------------------------------------------------------
C Calculates the derivatives of shape functions
C as they are mapped from the ξ-η-ζ (isoparametric) domain
C to the x-y-z (physical) domain
C --------------------------------------------------------

      SUBROUTINE mapShape3D(NNODE,dNxi,COORDS,dN,mapJ_det,stat)
C----------------------------------------------------------------------
C  Map derivatives of shape fns from xi-eta-zeta domain (dNxi)
C  to x-y-z domain (dN, d2N).
C  * dNxi = partials in param. space (xi,eta,zeta)
C  * dN    = partials in physical domain (x,y,z)
C----------------------------------------------------------------------
C
      IMPLICIT NONE

      INTEGER NNODE, stat
      INTEGER K, P, A
C dNxi(NNODE,3): ∂N_K / ∂ξ_i
      DOUBLE PRECISION dNxi(NNODE,3)
      DOUBLE PRECISION dN(NNODE,3)

C COORDS(3,NNODE): coordinates of each node in physical space
C COORDS(V,I) = x_V^(I) (V = x,y,z)
      DOUBLE PRECISION COORDS(3,NNODE)

C → mapJ = Jᵀ (transpose of the Jacobian)
      DOUBLE PRECISION mapJ(3,3)
C mapJ_inv = inverse of mapJ = (Jᵀ)⁻¹ = J⁻ᵀ
C mapJ_invt = transpose of mapJ_inv = (J⁻ᵀ)ᵀ = J⁻¹
C mapJ_det = determinant of mapJ = |J|
      DOUBLE PRECISION mapJ_inv(3,3), mapJ_invt(3,3)
      DOUBLE PRECISION mapJ_det
C some constant parameters
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.D0, ONE=1.0D0, TWO=2.0D0
      DOUBLE PRECISION, PARAMETER :: HALF=0.5D0, FOURTH=0.25D0
      DOUBLE PRECISION, PARAMETER :: EIGHTH=0.125D0

C  Initialize
      mapJ = ZERO
      mapJ_inv = ZERO
      mapJ_invt = ZERO
      mapJ_det = ZERO
      dN = ZERO
C Integer variables
      K = 0
      A = 0
      P = 0

C--------------------------------------------------------------------
C 1) Build the Jacobian: mapJ(m,p) = ∂x_m / ∂ξ_p = Jᵀ
C--------------------------------------------------------------------
      DO K=1,3
            DO P=1,3
                  DO A=1,NNODE
                        mapJ(K,P) = mapJ(K,P) + dNxi(A,K)*COORDS(P,A)   ! Jᵀ
                  END DO
            END DO
      END DO

C--------------------------------------------------------------------
C 2) Invert Jacobian: mapJ_inv = J⁻ᵀ, mapJ_invt = J⁻¹
C    Invert Jacobian: jacob_inv = J⁻¹, jacob_invt = (J⁻¹)ᵀ
C--------------------------------------------------------------------
      CALL matInvDet3D(mapJ,mapJ_inv,mapJ_invt,mapJ_det,stat)

C--------------------------------------------------------------------
C 3) First derivatives transform:
C    dN = (J⁻ᵀ) ⋅ dN/dξ
C    → using matrix notation: dN = J⁻ᵀ ⋅ dNxi
C    where:
C      mapJ_inv = J⁻ᵀ
C      mapJ_invt = J⁻¹
C--------------------------------------------------------------------
      dN = TRANSPOSE(MATMUL(mapJ_inv,TRANSPOSE(dNxi)))

C--------------------------------------------------------------------------
      RETURN
      END SUBROUTINE mapShape3D

C----------------------------------------------------------------------------------------------------------------------------------------------
C Calculates the derivatives of shape functions
C as they are mapped from the ξ-η-ζ (isoparametric) domain
C to the x-y-z (physical) domain
C --------------------------------------------------------

      SUBROUTINE SecondDerN(NNODE,dNxi,d2Nxi,COORDS,d2N,stat)
C----------------------------------------------------------------------
C  Map derivatives of shape fns from xi-eta-zeta domain (d2Nxi)
C  to x-y-z domain (dN, d2N).
C  * dNxi, d2Nxi = partials in param. space (xi,eta,zeta)
C  * d2N    = partials in physical domain (x,y,z)
C----------------------------------------------------------------------
C
      IMPLICIT NONE

      INTEGER NNODE, stat
      INTEGER K, A
      INTEGER S, R, P, Q, S2, T
      INTEGER MPHYS, NPHYS, PREF, QREF
C dNxi(NNODE,3): ∂N_K / ∂ξ_i
C d2Nxi(NNODE,3,3): ∂²N_K / ∂ξ_i∂ξ_j
      DOUBLE PRECISION dNxi(NNODE,3), d2Nxi(NNODE,3,3)
      DOUBLE PRECISION d2N(NNODE,3,3)

C COORDS(3,NNODE): coordinates of each node in physical space
C COORDS(V,I) = x_V^(I) (V = x,y,z)
      DOUBLE PRECISION COORDS(3,NNODE)
C jacob(i,j) = ∂x_i / ∂ξ_j
C → jacob = J
      DOUBLE PRECISION jacob(3,3)
C jacob_inv = inverse of jacob = J⁻¹
C jacob_invt = transpose of jacob_inv = (J⁻¹)ᵀ
C jacob_det = determinant of jacob = |J|
      DOUBLE PRECISION jacob_inv(3,3), jacob_invt(3,3)
      DOUBLE PRECISION jacob_det
C tmpJ(U,V,R) ≈ ∂J_{VU}/∂x_R = ∂²x_V / ∂ξ_U∂ξ_S * ∂ξ_S / ∂x_R
C dMapJdx(P,Q,R) = ∂(J⁻¹)_{P Q} / ∂x_R
      DOUBLE PRECISION tmpJ(3,3,3), dMapJdx(3,3,3)

C some constant parameters
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.D0, ONE=1.0D0, TWO=2.0D0
      DOUBLE PRECISION, PARAMETER :: HALF=0.5D0, FOURTH=0.25D0
      DOUBLE PRECISION, PARAMETER :: EIGHTH=0.125D0

C  Initialize
      jacob = ZERO
      jacob_inv = ZERO
      jacob_invt = ZERO
      jacob_det = ZERO
      tmpJ = ZERO
      dMapJdx = ZERO
      d2N = ZERO
C Integer variables
      K = 0
      A = 0
      S = 0
      R = 0
      P = 0
      Q = 0
      S2 = 0
      T = 0
      MPHYS = 0
      NPHYS = 0
      PREF = 0
      QREF = 0

C--------------------------------------------------------------------
C 1) Build the Jacobian: mapJ(m,p) = ∂x_m / ∂ξ_p = Jᵀ
C--------------------------------------------------------------------
      DO K=1,3
            DO P=1,3
                  DO A=1,NNODE
                        jacob(K,P) = jacob(K,P) + COORDS(K,A)*dNxi(A,P) ! J
                  END DO
            END DO
      END DO

C--------------------------------------------------------------------
C 2) Invert Jacobian: jacob_inv = J⁻¹, jacob_invt = (J⁻¹)ᵀ
C--------------------------------------------------------------------
      CALL matInvDet3D(jacob,jacob_inv,jacob_invt,jacob_det,stat)
C--------------------------------------------------------------------
C 3) Second derivative transform:
C    We apply the chain rule:
C    ∂²N/∂x_i∂x_j =
C        Σ_{m,l} ∂²N/∂ξ_m∂ξ_l ⋅ J⁻¹_{m i} ⋅ J⁻¹_{l j}
C      + Σ_m ∂N/∂ξ_m ⋅ [ ∂(J⁻¹_{m i})/∂x_j + ∂(J⁻¹_{m j})/∂x_i ]
C--------------------------------------------------------------------

C── 3.1) Build tmpJ(R,S,T): ∂J_{R,S} / ∂x_T
      DO A = 1, NNODE
            DO S = 1, 3                                                 ! reference direction
                  DO R = 1, 3                                           ! physical direction
                        DO T = 1, 3                                     ! physical derivative direction
                              DO S2 = 1, 3                              ! second reference direction
                                    tmpJ(R,S,T) = tmpJ(R,S,T)
     1                                          + COORDS(R,A)
     2                                          * d2Nxi(A,S,S2)
     3                                          * jacob_inv(S2,T)
                              END DO
                        END DO
                  END DO
            END DO
      END DO

C── 3.2) Build dMapJdx(P,Q,T) = ∂(J⁻¹)_{P,Q} / ∂x_T
      DO P = 1,3
            DO Q = 1,3
                  DO T = 1,3
                        DO R = 1,3
                              DO S = 1,3
                                    ! Inverse matrix derivative identity
                                    dMapJdx(P,Q,T) = dMapJdx(P,Q,T)
     1                                             - jacob_inv(P,R)     ! = J⁻¹(P,R)
     2                                             * tmpJ(R,S,T)        ! = ∂J_{R,S}/∂x_T
     3                                             * jacob_inv(S,Q)     ! = J⁻¹(S,Q)
                              END DO
                        END DO
                  END DO
            END DO
      END DO

C── 3.3) Compute full second derivatives in physical space: d2N(K,I,J)
      DO A = 1, NNODE                                                   ! Shape function index: A = N^A
            DO MPHYS = 1, 3                                             ! Physical direction x_m
                  DO NPHYS = MPHYS, 3                                   ! Physical direction x_n (symmetry: n ≥ m)
                        d2N(A,MPHYS,NPHYS) = ZERO                       ! clear accumulator
                        DO PREF = 1, 3                                  ! Reference ξ_p
                              DO QREF = 1, 3                            ! Reference ξ_q
                                    ! (a) pure chain‐rule piece
                                    d2N(A,MPHYS,NPHYS) = 
     1                                           d2N(A,MPHYS,NPHYS)
     2                                         + d2Nxi(A,PREF,QREF)
     3                                         * jacob_inv(PREF,MPHYS)  ! = J⁻¹_{p m}
     4                                         * jacob_inv(QREF,NPHYS)  ! = J⁻¹_{q n}
                              END DO                                    ! QREF
                              ! (b+c) correction terms:
                              ! = ∂(J⁻¹_{p m})/∂x_n + ∂(J⁻¹_{p n})/∂x_m
                              d2N(A,MPHYS,NPHYS) = d2N(A,MPHYS,NPHYS)
     1                                   + dNxi(A,PREF) * (
                                           ! (b) ∂(J⁻¹_{p m})/∂x_n
     2                                     dMapJdx(PREF,MPHYS,NPHYS) +
                                           ! (c) ∂(J⁻¹_{p n})/∂x_m
     3                                     dMapJdx(PREF,NPHYS,MPHYS)
     4                                     )
                        END DO                                          ! PREF
                        ! Symmetry by Clairaut’s theorem (also called Schwarz’s theorem)
                        d2N(A,NPHYS,MPHYS) = d2N(A,MPHYS,NPHYS)         ! ∂²N/∂x_n∂x_m = ∂²N/∂x_m∂x_n
                  END DO                                                ! NPHYS
            END DO                                                      ! MPHYS
      END DO                                                            ! A

C--------------------------------------------------------------------------
      RETURN
      END SUBROUTINE SecondDerN

C----------------------------------------------------------------------------------------------------------------------------------------------
C Computes everything required for the time integration of the problem
C --------------------------------------------------------------------
C This subroutine combines different behaviors to calculate the new state
C of these variables based on increments over a given time step
C

      SUBROUTINE integ(STEPNAMES, NSTEPS, DTIME, KSTEP, NDIM,
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

C VARIABLES:
C I,J,K: Loop integers
C STEPNAMES: Names of the steps
C NSTEPS: Number of steps
C DTIME: Time increment
C KSTEP: Current step
C NDIM: Number of dimensions

C NECESSARY VARIABLES OF THE PROBLEM:
C K_CELL_rho_0: value for the cell viability model

C INPUTS:
C PROPS (E, NU, D_WAT, CELL_rho_0, GlucThres, pHThres,C10)
C
C E: YOUNG'S MODULUS
C NU: POISSON'S RATIO
C CELL_rho_0: INITIAL CELL DENSITY
C a_K, b_K, c_K: Constants for the K(CELL_rho_0) model
C ALPHA_pHLACT: ALPHA_pHLACT: nmol/mL is a constant that quantifies change of pH per unit of lactate concentration
C D_WAT: WATER DIFFUSION COEFFICIENT OF THE SOLUTES
C ALPHA_pH: Death rate due to acidity
C GlucThres: Glucose threshold for cell viability
C pHThres: pH threshold for cell viability
C DFG_inv: Inverse of the deformation gradient
C DETDFG: Determinant of the deformation gradient
C dDFGdx: Gradient of the deformation gradient
C

C SDVs FROM PREVIOUS SIMULATION OF THE GLOBAL MODEL
C NF0: INITIAL WATER CONTENT (POROSITY BECAUSE THE PORUS MATRIX IS FULLY SATURATED)
C
C Old values of the local SDVs
C CELL_rho_old: OLD CELL DENSITY
C pH_val_old: OLD pH
C CELL_viab_old: OLD CELL VIABILITY
C
C SOLUTE CONCENTRATIONS (O2_int, LACT_int, GLUC_int)
C THE DERIVATIVES WERE CALCULATED IN THE UEL AS:
C DO I=1,NDIM
C       dSOLUTEdX(I,1) = dSOLUTEdX(I,1) + SOLUTE_CONC(K)*dNC(K,I)
C END DO
C REACTION RATES (ALPHA_pH, ALPHA_GLUC)

C OUTPUTS:
C DIFFUSION IN THE MEDIUM AND REACTION (CONSUMPTION) OF THE SOLUTES (D_SOL, R_SOL)
C GRADIENT OF THE DIFFUSION (dD_O2dx, dD_LACTdx, dD_GLUCdx)
C MAGNITUDE OF THE GRADIENT (dD_SOLdx_MAG)
C DERIVATIVE OF REACTION WRT TO O2, LACTATE, AND GLUCOSE (dR_SOL)
C ACTUAL WATER CONTENT OR POROSITY (NF)
C CELL_rho: ACTUAL CELL DENSITY
C pH_val: ACTUAL pH
C CELL_viab: ACTUAL CELL VIABILITY
C PN equation variables:
C SP: Total Stress Pressure
C IL1B: Interleukin 1 beta
C IL1B_PROT: Interleukin 1 beta protein
C TNF: Tumor necrosis factor
C TNF_PROT: Tumor necrosis factor protein
C AGG: Aggrecan expression
C COLI: Collagen I expression
C COLII: Collagen II expression
C MMP3: Matrix metalloproteinase expression
C ADM4: A disintegrin and metalloproteinase with thrombospondin motifs expression
C ACC: Accumulation of the protein (depends on PN_TIME in hours)

C
C ********************************************************************
C
      IMPLICIT NONE
C
C FUNCTION
      DOUBLE PRECISION ROUND
C
C VARIABLES:
C LOOP INTEGERS
      INTEGER I, J, K, R, A, B
C SIMULATION VALUES
      INTEGER NSTEPS
      CHARACTER*256 STEPNAMES(NSTEPS)
      DOUBLE PRECISION DTIME
      INTEGER KSTEP, NDIM

C NECESSARY VARIABLES OF THE PROBLEM:
C DEFORMATION GRADIENT
      DOUBLE PRECISION DFG_inv(NDIM,NDIM),DETDFG,dDFGdx(NDIM,NDIM,NDIM)

C INPUTS:
C MATERIAL PROPERTIES AS PROPS
      DOUBLE PRECISION D_WAT(3), CELL_rho_0, ALPHA_pH
      DOUBLE PRECISION a_K, b_K, c_K
      DOUBLE PRECISION GlucThres, pHThres
C SDVs FROM PREVIOUS SIMULATION OF THE GLOBAL MODEL
C NF0: INITIAL WATER CONTENT (POROSITY BECAUSE THE PORUS MATRIX IS FULLY SATURATED)
      DOUBLE PRECISION NF0
C SOLUTE CONCENTRATIONS
      DOUBLE PRECISION O2_int, LACT_int, GLUC_int
C REACTION RATES
      DOUBLE PRECISION ALPHA_GLUC
C PN equation variables
C FREQ: Frequency of load
C PN_TIME: Total desired protein expresion simulation time (hours)
      DOUBLE PRECISION FREQ
      INTEGER PN_TIME, PN_TIME_file
      INTEGER FREQ_activation
      DOUBLE PRECISION FREQUENCY_CURRENT, PN_FREQ_file(PN_TIME_file)

C OUTPUTS:
C DIFFUSION IN THE MEDIUM AND REACTION (CONSUMPTION) OF THE SOLUTES:
      DOUBLE PRECISION D_SOL(3), R_SOL(3)
C GRADIENT OF THE DIFFUSION
      DOUBLE PRECISION prefdD_SOLdx
      DOUBLE PRECISION dD_O2dx(NDIM,1)
      DOUBLE PRECISION dD_LACTdx(NDIM,1)
      DOUBLE PRECISION dD_GLUCdx(NDIM,1)
      DOUBLE PRECISION dD_SOLdx_MAG
C nD_d: Use ∇D as ZERO (nD_d=0), Isotropic (nD_d=1), Anisotropic (nD_d=2)
      INTEGER nD_d
C DERIVATIVE OF REACTION WRT TO O2, LACTATE, AND GLUCOSE
      DOUBLE PRECISION dR_SOL(3,3)
c ACTUAL WATER CONTENT OF POROSITY, CELL DENSITY, CELL VIABILITY pH
      DOUBLE PRECISION NF, CELL_rho, CELL_viab, pH_val
C OLD CELL VIABILITY
      DOUBLE PRECISION CELL_rho_old, pH_val_old, CELL_viab_old
C PN equation variables
      DOUBLE PRECISION SP
      DOUBLE PRECISION IL1B, IL1B_PROT, TNF, TNF_PROT
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
      INTEGER PNeq_STEP

C OTHER VARIABLES:
C MATERIAL NAME
      CHARACTER*80 MATNAME

C TEMPORAL VARIABLES
C DIFFUSION TRACE_D: tr(F^{-1} ∂F/∂x_k)
      DOUBLE PRECISION TRACE_D(NDIM)
C CELL VIABILITY
      DOUBLE PRECISION DIFF, aDIFF, SUM_var, K_CELL_rho_0
C O2 
      DOUBLE PRECISION AA,BB,CC,DD,B1,D1
C LACTATE
      DOUBLE PRECISION AL, BL, CL, DL

C CONSTANT PARAMETERS:
C SIMPLE DOUBLE PRECISION NUMBERS
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0, ONE=1.0D0, TWO=2.0D0
      DOUBLE PRECISION, PARAMETER :: THREE=3.0D0, FOUR=4.0D0, SIX=6.0D0
      DOUBLE PRECISION, PARAMETER :: HALF=0.5D0, FOURTH=0.25D0
      DOUBLE PRECISION, PARAMETER :: EIGHTH=0.125D0

C Pressure necessary variables
C ----------------------------
      DOUBLE PRECISION, PARAMETER :: RR=8.3145D0                        !Gas constant [N mm /(mmol K)]
      DOUBLE PRECISION, PARAMETER :: TT=293.0D0                         !Absolute temperature [K]
      DOUBLE PRECISION, PARAMETER :: rhoS=1.4338D0                      !Mass density of solid matrix [g/m]) (Shapiro et al., Osteoarthritis and Cartilage, 9:533-539 (2001); Basser et al., Arch Biochem Biophys, 351:207-219 (1998))

C Nutrient transport necessary variables
C --------------------------------------
C ALPHA_pHLACT: nmol/mL is a constant that quantifies change of pH per unit of lactate concentration
C 1/A = 11.11 nmol/mL is a constant that quantifies change of pH per unit of lactate concentration
      DOUBLE PRECISION, PARAMETER :: ALPHA_pHLACT = 1.0D0/11.11D0
C OXYGEN SOLUBILITY IN WATER
      DOUBLE PRECISION, PARAMETER :: SO2 = 0.010268D0
C ALPHA_GTF: ALPHA GLUCOSE TIME FACTOR
      DOUBLE PRECISION, PARAMETER :: ALPHA_GTF = 1.0D0/86400.0D0 ! 1/day = 0.0000115741D0
C dCELL_rhodGluc: Derivative of cell density wrt glucose
      DOUBLE PRECISION dCELL_rhodGluc

C--------------------------------------------------------------------------
C Initializing variables
      I = 0
      J = 0
      K = 0
      R = 0
      A = 0
      B = 0
      NF = ZERO
      D_SOL = ZERO
      prefdD_SOLdx = ZERO
      dD_O2dx = ZERO
      dD_LACTdx = ZERO
      dD_GLUCdx = ZERO
      dD_SOLdx_MAG = ZERO
      DIFF = ZERO
      aDIFF = ZERO
      SUM_var = ZERO
      K_CELL_rho_0 = ZERO
      R_SOL = ZERO
      dR_SOL = ZERO
      ALPHA_GLUC = ZERO
      PNeq_STEP = 0
C DIFFUSION
      TRACE_D = ZERO
C O2
      AA = ZERO
      BB = ZERO
      CC = ZERO
      DD = ZERO
      B1 = ZERO
      D1 = ZERO
C LACTATE
      AL = ZERO
      BL = ZERO
      CL = ZERO
      DL = ZERO
C GLUCOSE
      dCELL_rhodGluc = ZERO

C--------------------------------------------------------------------------
C Uncomment to use the determinant function developed by Malandrino
C      DETDFG = ZERO
C      CALL mDETDFG(DFG,DETDFG,3)

C--------------------------------------------------------------------------
C Check if O2 is lower than ZERO, if so, set it to ZERO
      IF (O2_int .LT. ZERO) THEN
            O2_int = ZERO
      END IF

C Check if Lactate is lower than ZERO, if so, set it to ZERO
      IF (LACT_int .LT. ZERO) THEN
            LACT_int = ZERO
      END IF

C Check if Glucose is lower than ZERO, if so, set it to ZERO
      IF (GLUC_int .LT. ZERO) THEN
            GLUC_int = ZERO
      END IF

C--------------------------------------------------------------------------
C     Calculating the Water Content (NF) to be compared between the
C     poro-mechanical model, UMAT, and UEL
C     here we use the Initial Water Content (NF0)
      IF (DETDFG .EQ. ZERO) THEN
            NF = NF0
            DETDFG = ONE
      ELSE
            NF = (NF0-ONE+DETDFG)/DETDFG
      END IF

C
C********************************************************************
C
C               Diffusion and flux calculation
C              Based on Mackie and Meares (1955)
C  
C********************************************************************
C

C O2, Lactate, and Glucose Diffusivity in the medium D_SOL(1,2,3) considering the porosity
      DO I = 1,3
            D_SOL(I) = D_WAT(I)*(NF/(TWO-NF))**TWO
      ENDDO

C
C********************************************************************
C
C                   Gradient of the diffusion
C           Based on Workineh & Muñoz-Moya et al. (2025)
C  
C********************************************************************
C
C--------------------------------------------------------------
C Prefactor common to all solutes (depends only on n_f)
C--------------------------------------------------------------
      prefdD_SOLdx = FOUR * NF * (ONE - NF) / ( (TWO - NF)**THREE )

C--------------------------------------------------------------
C Compute slice-wise trace  tr( F^{-1} * ∂F/∂x_k )
C TRACE_D(k) holds the value for gradient direction k = 1..NDIM
C--------------------------------------------------------------
      DO K = 1, NDIM                                                    ! gradient direction  x, y, z
            DO A = 1, NDIM
                  DO B = 1, NDIM
                        TRACE_D(K) = TRACE_D(K) +
     1                               DFG_inv(A,B) * dDFGdx(K,B,A)
                  END DO
            END DO                                                      !  (notice J,I to pick up the transpose inside the trace)
C --------------------------------------------------------------------------
C Now fill dD_SOLdx(K, I) for each solute using your prefactor
C Fill the diffusivity gradients  dD_SOLdx(K, I)
C I = 1:O2 , 2:Lact , 3:Gluc
C--------------------------------------------------------------
            dD_O2dx(K,1) = D_WAT(1) * prefdD_SOLdx * TRACE_D(K)
            dD_LACTdx(K,1) = D_WAT(2) * prefdD_SOLdx * TRACE_D(K)
            dD_GLUCdx(K,1) = D_WAT(3) * prefdD_SOLdx * TRACE_D(K)
C Magnitude of the gradient
            dD_SOLdx_MAG = dD_SOLdx_MAG +
     1                     (prefdD_SOLdx * TRACE_D(K))**TWO
      END DO

      dD_SOLdx_MAG = DSQRT(dD_SOLdx_MAG)

C nD_d = 0: Zero gradient
      IF (nD_d .EQ. 0) THEN
            dD_O2dx = ZERO
            dD_LACTdx = ZERO
            dD_GLUCdx = ZERO
            dD_SOLdx_MAG = ZERO
C nD_d = 1: Isotropic gradient
      ELSE IF (nD_d .EQ. 1) THEN
            DO K = 1,3
                  dD_O2dx(K,1) = D_WAT(1) * dD_SOLdx_MAG
     1                           / DSQRT(THREE)
                  dD_LACTdx(K,1) = D_WAT(2) * dD_SOLdx_MAG
     1                           / DSQRT(THREE)
                  dD_GLUCdx(K,1) = D_WAT(3) * dD_SOLdx_MAG
     1                           / DSQRT(THREE)
            END DO
      ELSE
      ! nD_d = 2 → keep anisotropic result (no change)
      END IF
C
C--------------------------------------------------------------------
C
C             Parameters for cell viability calculation
C
C--------------------------------------------------------------------
C

C FIRST STEP: INITIAL CONDITIONS
      IF (KSTEP .LT. 2) THEN
            pH_val = ph_val_old
            CELL_rho = CELL_rho_old
            CELL_viab = CELL_viab_old

C NEXT STEPS
      ELSE

C  Bibby et al. (2005), Soukane et al. (2005, 2007), and Magnier et al. (2009)
C Calculation of pH values using lactate concentration
      pH_val = 7.4D0 - ALPHA_pHLACT*LACT_int

C Muñoz-Moya et al. (2025)
C Defining the value of K(CELL_rho_0)
      K_CELL_rho_0 = -(CELL_rho_0 * 10.0D0**THREE - a_K) /
     1          DSQRT((CELL_rho_0 * 10.0D0**THREE
     2          - a_K)**TWO+ 10.0D0**c_K)
     3          * ((GlucThres-0.001D0)-b_K)/TWO  +
     4          ((GlucThres-0.001D0)-b_K)/TWO  + b_K

C Zhu et al. (2012)
C ALPHA_GLUC for glucose death rate
      SUM_var = GLUC_int + K_CELL_rho_0
      DIFF = GLUC_int - GlucThres
      aDIFF = ABS(DIFF)
      ALPHA_GLUC = ALPHA_GTF*((DIFF/SUM_var)-(aDIFF/SUM_var))           ! glucose death rate

C--------------------------------------------------------------------------
C FIRST CONDITION: NO ACTIVATION      
      IF ((GLUC_int .GE. GlucThres).AND.(pH_val .GE. pHThres)) THEN
      !CELL VIABILITY
      CELL_viab = CELL_viab_old
      !DERIVATIVE OF CELL DENSITY WRT GLUCOSE
      dCELL_rhodGluc = ZERO

C--------------------------------------------------------------------------
C Zhu et al. (2012)
C SECOND CONDITION: GLUCOSE ACTIVATE THE CELL
      ELSE IF ((GLUC_int .LT. GlucThres).AND.(pH_val .GE. pHThres)) THEN
      !CELL VIABILITY
      CELL_viab = CELL_viab_old*EXP(ALPHA_GLUC*DTIME)
      !DERIVATIVE OF CELL DENSITY WRT GLUCOSE
      dCELL_rhodGluc = ALPHA_GTF *
     1                 ((aDIFF - DIFF) * (aDIFF + SUM_var))/
     2                 (aDIFF * SUM_var**TWO)
      dCELL_rhodGluc = CELL_viab*DTIME*dCELL_rhodGluc

C--------------------------------------------------------------------------
C Malandrino et al. (2015b) and Muñoz-Moya et al. (2025)
C THIRD CONDITION: pH ACTIVATE THE CELL
      ELSE IF ((GLUC_int .GE. GlucThres).AND.(pH_val .LT. pHThres)) THEN
      CELL_viab = CELL_viab_old*EXP(ALPHA_pH*DTIME)
      dCELL_rhodGluc = ZERO

C--------------------------------------------------------------------------
C Zhu et al. (2012), Malandrino et al. (2015b) and Muñoz-Moya et al. (2025)
C FOURTH CONDITION: BOTH ACTIVATE THE CELL
      ELSE IF ((GLUC_int .LT. GlucThres).AND.(pH_val .LT. pHThres)) THEN
      !CELL VIABILITY
      CELL_viab = CELL_viab_old*EXP((ALPHA_GLUC+ALPHA_pH)*DTIME)
      !DERIVATIVE OF CELL DENSITY WRT GLUCOSE
      dCELL_rhodGluc = ALPHA_GTF *
     1                 ((aDIFF - DIFF) * (aDIFF + SUM_var))/
     2                 (aDIFF * SUM_var**TWO)
      dCELL_rhodGluc = CELL_viab*DTIME*dCELL_rhodGluc

      END IF

C--------------------------------------------------------------------------
C Uncomment the following lines to use the Initial Contidions of the variables for testing
C      pH_val = 7.4D0
C      CELL_viab = ONE
C      dCELL_rhodGluc = ZERO

C CELL DENSITY
      CELL_rho = CELL_rho_0*CELL_viab/DETDFG
      dCELL_rhodGluc = CELL_rho_0*dCELL_rhodGluc/DETDFG

      ENDIF

C
C********************************************************************
C
C                   Metabolic reactions calculation
C                    Based on Bibby et al. (2005)
C  
C********************************************************************
C

C
C--------------------------------------------------------------------
C
C                    Reaction for Oxygen by tissue
C
C--------------------------------------------------------------------
C

C Oxygen rate
      CC = 1.46D0
      B1 = 4.95D0
      D1 = 4.03D0

      BB = pH_val-B1
      DD = D1*BB

      AA = (7.28D0*CELL_rho*NF)/(3600.D0*SO2)
      R_SOL(1) = -(AA*BB*O2_int)/(CC+DD+O2_int)

C The derivatives were calculated considering this reaction:
C R_SOL(1) = -(AA*(7.4D0 - ALPHA_pHLACT*LACT_int-B1)*O2_int)/(CC+D1*(7.4D0 - ALPHA_pHLACT*LACT_int-B1)+O2_int)

C Derivative of O2 Reaction
      dR_SOL(1,1) = -(AA*BB*(CC+DD))/((CC+DD+O2_int)**TWO)

C Derivative of R_SOL(1) wrt Lactate Reaction
      dR_SOL(1,2) = (AA*O2_int*ALPHA_pHLACT*(CC+O2_int))/
     1              (CC+DD+O2_int)**TWO

C Derivative of R_SOL(1) wrt Glucose Reaction
      dR_SOL(1,3) = R_SOL(1)/CELL_rho*dCELL_rhodGluc

C
C--------------------------------------------------------------------
C
C                   Reaction for Lactate by tissue
C
C--------------------------------------------------------------------
C

C Lactate rate
      AL = -2.47D0
      BL = 0.93D0
      CL = 0.16D0
      DL = -0.0058D0

      R_SOL(2) = (CELL_rho/3600.D0) *  
     1           EXP(AL + BL*pH_val + CL*O2_int + DL*O2_int*O2_int)

C Derivative R_SOL(2) wrt O2 Reaction
      dR_SOL(2,1) = R_SOL(2)*(CL + 2*DL*O2_int)

C Derivative of Lactate Reaction
      dR_SOL(2,2) = -R_SOL(2)*BL*ALPHA_pHLACT

C Derivative R_SOL(2) wrt Glucose Reaction
      dR_SOL(2,3) = R_SOL(2)/CELL_rho*dCELL_rhodGluc

C
C--------------------------------------------------------------------
C
C                   Reaction for Glucose by tissue
C
C--------------------------------------------------------------------
C

C Glucose rate
      R_SOL(3) = -HALF*R_SOL(2)

C Derivative R_SOL(2) wrt O2 Reaction
      dR_SOL(3,1) = R_SOL(3)*(CL + 2*DL*O2_int)

C Derivative of Lactate Reaction
      dR_SOL(3,2) = -R_SOL(3)*BL*ALPHA_pHLACT

C Derivative of Glucose Reaction
      dR_SOL(3,3) = R_SOL(3)/CELL_rho*dCELL_rhodGluc

C
C********************************************************************
C
C                      PN equation calculation
C                     Baumgartner et al. (2025)
C  
C********************************************************************
C

C
C--------------------------------------------------------------------
C
C           Protein expression calculation in nucleus cells
C
C--------------------------------------------------------------------
C
      IF (MATNAME .EQ. 'NP') THEN

C Check if the frequency file is activated for this step and we are not in the PN equation step
      IF ((FREQ_activation .EQ. 1) .AND.
     1    (STEPNAMES(KSTEP) .NE. 'PN_eq')) THEN
            FREQ = FREQUENCY_CURRENT
C If the frequency file is activated for this step and we are in the PN equation step
      ELSE IF ((FREQ_activation .EQ. 1) .AND.
     1    (STEPNAMES(KSTEP) .EQ. 'PN_eq')) THEN
            PNeq_STEP = 1
            PN_TIME = PN_TIME_file
      END IF

      CALL PNeq(
     +           ROUND(ABS(GLUC_int), 4), ROUND(ABS(pH_val), 4),
     +           ROUND(ABS(SP), 4), ROUND(ABS(FREQ), 4),
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

      END IF

C--------------------------------------------------------------------------
      RETURN
      END SUBROUTINE integ

C----------------------------------------------------------------------------------------------------------------------------------------------
C Subroutine to assemble the local elements residual and tangent
C --------------------------------------------------------------------
C
      SUBROUTINE AssembleElement(NDIM,NNODE,NDOFEL,NDOFN,
     +           RHS,NRHS,AMATRX,
     +           R_O2,K_O2,
     +           R_LACT,K_LACT,
     +           R_GLUC,K_GLUC,
     +           K_O2LACT,K_O2GLUC,
     +           K_LACTO2,K_LACTGLUC,
     +           K_GLUCO2,K_GLUCLACT)
C
      IMPLICIT NONE
C
      INTEGER I, J
      INTEGER A11, A12, B11, B12
      INTEGER NDIM, NNODE, NDOFEL, NDOFN, NRHS

C RHS AND AMATRX ARE THE LOCAL ELEMENT RESIDUAL AND TANGENT
      DOUBLE PRECISION RHS(NDOFEL,NRHS),AMATRX(NDOFEL,NDOFEL)

C RHS AND STIFFNESS MATRX FOR EACH SOLUTE
C O2
      DOUBLE PRECISION R_O2(NNODE,1),K_O2(NNODE,NNODE)
      DOUBLE PRECISION K_O2LACT(NNODE,NNODE),K_O2GLUC(NNODE,NNODE)
C Lactate
      DOUBLE PRECISION R_LACT(NNODE,1),K_LACT(NNODE,NNODE)
      DOUBLE PRECISION K_LACTO2(NNODE,NNODE),K_LACTGLUC(NNODE,NNODE)
C Glucose
      DOUBLE PRECISION R_GLUC(NNODE,1),K_GLUC(NNODE,NNODE)
      DOUBLE PRECISION K_GLUCO2(NNODE,NNODE),K_GLUCLACT(NNODE,NNODE)

C CONSTANT PARAMETERS
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.D0, ONE=1.0D0, TWO=2.0D0
      DOUBLE PRECISION, PARAMETER :: HALF=0.5D0, FOURTH=0.25D0
      DOUBLE PRECISION, PARAMETER :: EIGHTH=0.125D0

C--------------------------------------------------------------------------
C Initializing RHS and AMATRX
      RHS(:,NRHS) = ZERO
      AMATRX = ZERO

C--------------------------------------------------------------------------
C Loop over nodes to assemble the element level residual
C Assemble the element level residual
      DO I=1,NNODE
            A11 = NDOFN*(I-1)+1
            A12 = NDIM*(I-1)+1

C Displacement: Uncomment if you want to edit the displacement
C            RHS(A11,NRHS) = ZERO
C            RHS(A11+1,NRHS) = ZERO
C            RHS(A11+2,NRHS) = ZERO

C Dummy Temperature: It is here to perform a Coupled temperature-displacement analysis, but it is not used in this model
C            RHS(A11+3,NRHS) = ZERO

C--------------------------------------------------------------------------
C Extra DoF
C O2 Concentration
            RHS(A11+4,NRHS) = -R_O2(I,1)

C Lactate Concentration
            RHS(A11+5,NRHS) = -R_LACT(I,1)

C Glucose Concentration
            RHS(A11+6,NRHS) = -R_GLUC(I,1)
      ENDDO

C--------------------------------------------------------------------------
C Loop over nodes to assemble the element level tangent matrix
C Assembly the element level tangent matrix
      DO I=1,NNODE
            DO J=1,NNODE
                  A11 = NDOFN*(I-1)+1
                  A12 = NDIM*(I-1)+1                  
                  B11 = NDOFN*(J-1)+1
                  B12 = NDIM*(J-1)+1

C Displacement: Uncomment if you want to edit the displacement
                  ! AMATRX(A11,B11) = ZERO
                  ! AMATRX(A11,B11+1) = ZERO
                  ! AMATRX(A11,B11+2) = ZERO
                  ! AMATRX(A11+1,B11) = ZERO
                  ! AMATRX(A11+1,B11+1) = ZERO
                  ! AMATRX(A11+1,B11+2) = ZERO
                  ! AMATRX(A11+2,B11) = ZERO
                  ! AMATRX(A11+2,B11+1) = ZERO
                  ! AMATRX(A11+2,B11+2) = ZERO

C Dummy Temperature: It is here to perform a Coupled temperature-displacement analysis, but it is not used in this model
                  ! AMATRX(A11+3,B11+3) = ZERO

C--------------------------------------------------------------------------
C Extra DoF: Diagonal terms

C O2 Concentration
                  AMATRX(A11+4,B11+4) = -K_O2(I,J)

C Lactate Concentration
                  AMATRX(A11+5,B11+5) = -K_LACT(I,J)

C Glucose Concentration
                  AMATRX(A11+6,B11+6) = -K_GLUC(I,J)

C--------------------------------------------------------------------------
C Extra DoF: Off-diagonal terms:

C O2 - Lactate
                  AMATRX(A11+4,B11+5) = -K_O2LACT(I,J)

C O2 - Glucose
                  AMATRX(A11+4,B11+6) = -K_O2GLUC(I,J)
      
C Lactate - O2
                  AMATRX(A11+5,B11+4) = -K_LACTO2(I,J)

C Lactate - Glucose
                  AMATRX(A11+5,B11+6) = -K_LACTGLUC(I,J)

C Glucose - O2
                  AMATRX(A11+6,B11+4) = -K_GLUCO2(I,J)

C Glucose - Lactate
                  AMATRX(A11+6,B11+5) = -K_GLUCLACT(I,J)

            END DO

      END DO


      END SUBROUTINE AssembleElement

C----------------------------------------------------------------------------------------------------------------------------------------------