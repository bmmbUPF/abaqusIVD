C***********************************************************************

C GLOBAL SDVs

      IF(.NOT.allocated(globalSDVs)) THEN
C
C allocate memory for the globalSdv's
C
C numElem needs to be set in the MODULE
C numInt needs to be set in the UEL
C

      stat=0
      ALLOCATE(globalSDVs(numIdElem, numInt, ngSdv),stat=err)
      globalSDVs = ZERO
      OPEN(15,FILE=infoFilePath_T_D, status='old', position='append')
      IF(stat.ne.0) THEN
      WRITE(15,*) '//////////////////////////////////////////////'
      WRITE(15,*) 'error when allocating globalSdv'
      WRITE(15,*) '//////////////////////////////////////////////'
      WRITE(15,*) 'stat =',stat
      WRITE(15,*) 'ngSdv =',ngSdv
      WRITE(15,*) 'numInt =',numInt
      WRITE(15,*) 'numElem =',numElem
      WRITE(15,*) 'NNODE =',NNODE
      WRITE(15,*) 'lbound(globalSDVs)',lbound(globalSDVs)
      WRITE(15,*) 'ubound(globalSDVs)',ubound(globalSDVs)
      WRITE(15,*) '//////////////////////////////////////////////'
      WRITE(15,*) ' '

      PRINT *, '//////////////////////////////////////////////'
      PRINT *, 'error when allocating globalSdv'
      PRINT *, '//////////////////////////////////////////////'
      PRINT *, 'stat =',stat
      PRINT *, 'ngSdv =',ngSdv
      PRINT *, 'numInt =',numInt
      PRINT *, 'numElem =',numElem
      PRINT *, 'NNODE =',NNODE
      PRINT *, 'lbound(globalSDVs)',lbound(globalSDVs)
      PRINT *, 'ubound(globalSDVs)',ubound(globalSDVs)
      PRINT *, '//////////////////////////////////////////////'
      PRINT *, ' '
      CALL xit
      ENDIF
C
      WRITE(15,*) '-------------------------------------------------'
      WRITE(15,*) '-------------- SDVs ALLOCATED -------------------'
      WRITE(15,*) '-------------------------------------------------'
      WRITE(15,*) '---------- YOU PUT NUMBER OF ELEMENTS -----------'
      WRITE(15,*) '---------- numElem =',numElem
      WRITE(15,*) '---------- U3D20 ELEMENTS -----------------------'
      WRITE(15,*) '-------------------------------------------------'
      WRITE(15,*) '---------- YOU PUT NUMBER OF POINTS -------------'
      WRITE(15,*) '---------- numInt =',numInt
      WRITE(15,*) '-------------------------------------------------'
      WRITE(15,*) '---------- YOU PUT NUMBER OF SDVs ---------------'
      WRITE(15,*) '---------- Global SDVs (ngSdv) =',ngSdv
      WRITE(15,*) '-------------------------------------------------'
      WRITE(15,*) 'lbound(globalSDVs)',lbound(globalSDVs)
      WRITE(15,*) 'ubound(globalSDVs)',ubound(globalSDVs)
      WRITE(15,*) '-------------------------------------------------'
      WRITE(15,*) ' '

      PRINT *, '-------------------------------------------------'
      PRINT *, '-------------- SDVs ALLOCATED -------------------'
      PRINT *, '-------------------------------------------------'
      PRINT *, '---------- YOU PUT NUMBER OF ELEMENTS -----------'
      PRINT *, '---------- numElem =',numElem
      PRINT *, '---------- U3D20 ELEMENTS -----------------------'
      PRINT *, '-------------------------------------------------'
      PRINT *, '---------- YOU PUT NUMBER OF POINTS -------------'
      PRINT *, '---------- numInt =',numInt
      PRINT *, '-------------------------------------------------'
      PRINT *, '---------- YOU PUT NUMBER OF SDVs ---------------'
      PRINT *, '---------- Global SDVs (ngSdv) =',ngSdv
      PRINT *, '-------------------------------------------------'
      PRINT *, 'lbound(globalSDVs)',lbound(globalSDVs)
      PRINT *, 'ubound(globalSDVs)',ubound(globalSDVs)
      PRINT *, '-------------------------------------------------'
      PRINT *, ' '

      ENDIF
      CLOSE(15)

C LOCAL SDVs

!       IF(.NOT.allocated(localSDVs)) THEN
! C
! C allocate memory for the localSdv's
! C
! C numElem needs to be set in the MODULE
! C numInt needs to be set in the UEL
! C
!       stat=0
!       ALLOCATE(localSDVs(numIdElem, numInt, nlSdv),stat=err)
!       localSDVs = ZERO
!       OPEN(15,FILE=infoFilePath_T_D, status='old', position='append')
!       IF(stat.ne.0) THEN
!       WRITE(15,*) '//////////////////////////////////////////////'
!       WRITE(15,*) 'error when allocating localSdv'
!       WRITE(15,*) '//////////////////////////////////////////////'
!       WRITE(15,*) 'stat=',stat
!       WRITE(15,*) 'nlSdv=',nlSdv
!       WRITE(15,*) 'numInt=',numInt
!       WRITE(15,*) 'numElem=',numElem
!       WRITE(15,*) 'NNODE=',NNODE
!       WRITE(15,*) 'lbound(localSDVs)',lbound(localSDVs)
!       WRITE(15,*) 'ubound(localSDVs)',ubound(localSDVs)
!       WRITE(15,*) '//////////////////////////////////////////////'
!       WRITE(15,*) ' '
!       CALL xit
!       ENDIF
! C
!       WRITE(15,*) '-------------------------------------------------'
!       WRITE(15,*) '----------- localSDV ALLOCATED ------------------'
!       WRITE(15,*) '-------------------------------------------------'
!       WRITE(15,*) '---------- YOU PUT NUMBER OF ELEMENTS -----------'
!       WRITE(15,*) '---------- numElem=',numElem
!       WRITE(15,*) '---------- U3D20 ELEMENTS -----------------------'
!       WRITE(15,*) '-------------------------------------------------'
!       WRITE(15,*) '---------- YOU PUT NUMBER OF POINTS -------------'
!       WRITE(15,*) '---------- numInt =',numInt
!       WRITE(15,*) '-------------------------------------------------'
!       WRITE(15,*) '---------- YOU PUT NUMBER OF SDVs ---------------'
!       WRITE(15,*) '---------- nlSdv=',nlSdv
!       WRITE(15,*) '-------------------------------------------------'
!       WRITE(15,*) 'lbound(localSDVs)',lbound(localSDVs)
!       WRITE(15,*) 'ubound(localSDVs)',ubound(localSDVs)
!       WRITE(15,*) '-------------------------------------------------'
!       WRITE(15,*) ' '

!       PRINT *, '-------------------------------------------------'
!       PRINT *, '----------- localSDV ALLOCATED ------------------'
!       PRINT *, '-------------------------------------------------'
!       PRINT *, '---------- YOU PUT NUMBER OF ELEMENTS -----------'
!       PRINT *, '---------- numElem=',numElem
!       PRINT *, '---------- U3D20 ELEMENTS -----------------------'
!       PRINT *, '-------------------------------------------------'
!       PRINT *, '---------- YOU PUT NUMBER OF POINTS -------------'
!       PRINT *, '---------- numInt =',numInt
!       PRINT *, '-------------------------------------------------'
!       PRINT *, '---------- YOU PUT NUMBER OF SDVs ---------------'
!       PRINT *, '---------- nlSdv=',nlSdv
!       PRINT *, '-------------------------------------------------'
!       PRINT *, 'lbound(localSDVs)',lbound(localSDVs)
!       PRINT *, 'ubound(localSDVs)',ubound(localSDVs)
!       PRINT *, '-------------------------------------------------'
!       PRINT *, ' '

!       ENDIF
!       CLOSE(15)



C***********************************************************************