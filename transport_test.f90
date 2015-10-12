PROGRAM transport_test 
  USE chaco_data
  USE parameters

IMPLICIT NONE

!******************************************************************************
! Data dictionary
!******************************************************************************
INTEGER(KIND=SHORT) , ALLOCATABLE, DIMENSION(:,:,:), TARGET  :: matrix         ! Ice-mantle matrix
INTEGER(KIND=SHORT)              , DIMENSION(:,:,:), POINTER :: matrix_ptr     ! Pointer to the matrix
INTEGER                                                      :: n,i,j,k          ! Counters
INTEGER, DIMENSION(3) :: prev,curr,next

i = 1000
j = 10000
ALLOCATE ( matrix(i,i,i) )
matrix_ptr => matrix
PRINT *, 'Allocated matrix'

CALL RANDOM_SEED()

!******************************************************************************
! Begin the simulation 
!******************************************************************************

OPEN(UNIT=1012,FILE="transport_track.csv",POSITION='APPEND', STATUS='REPLACE')
PRINT *, 'Opened file...'

prev = i/2 
curr = (i/2) + 1 
PRINT *, 'prev=',prev
PRINT *, 'curr=',curr
PRINT *, 'Starting simulation...'
DO n=1,j
!  PRINT *, n
  CALL transport(prev,curr,next,matrix_ptr)
  WRITE(1012,*) next(1),',',next(2),',',next(3) 
!  PRINT *, next(1),',',next(2),',',next(3)
  prev = curr
  curr = next
END DO
PRINT *, 'Ending simulation...'

CLOSE(1012)
END PROGRAM transport_test 
