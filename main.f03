PROGRAM main

USE subroutines
USE parameters
USE typedefs
USE functiondefs
IMPLICIT NONE

!******************************************************************************
! Data dictionary
!******************************************************************************
INTEGER                                                      :: n
INTEGER             , ALLOCATABLE, DIMENSION(:,:,:), TARGET  :: matrix         ! Ice-mantle matrix
INTEGER(KIND=SHORT) , ALLOCATABLE, DIMENSION(:,:,:), TARGET  :: qube
INTEGER             , ALLOCATABLE, DIMENSION(:)    , TARGET  :: anion_target
INTEGER                          , DIMENSION(:)    , POINTER :: anion_list     ! List of anionic species
INTEGER                          , DIMENSION(:,:,:), POINTER :: matrix_ptr     ! Pointer to the matrix
INTEGER(KIND=SHORT)              , DIMENSION(:,:,:), POINTER :: qube_ptr       ! Pointer to the reaction cube
INTEGER                                            , TARGET  :: wlen_target
INTEGER                                            , POINTER :: wait_len       ! Length of nonzero waitlist elements
INTEGER(KIND=SHORT)              , DIMENSION(3)              :: ev_nums
INTEGER                                                      :: mindex
INTEGER                          , DIMENSION(3)              :: dimens         ! Dimensions of the matrix
INTEGER                                                      :: i,j,k          ! Counters
INTEGER                                                      :: err1, err2     ! Error numbers for the files
INTEGER                                                      :: cr_num         ! species number for cosmic rays
INTEGER(KIND=SHORT)                                          :: lines_spec     ! Number of lines in species file
INTEGER(KIND=SHORT)                                          :: lines_react    ! Number of lines in reactions file
INTEGER                                                      :: count_num      ! Number of abundance file
INTEGER                                                      :: result
INTEGER                                                      :: time_check     ! DEBUGGING VAR
REAL                , ALLOCATABLE, DIMENSION(:,:)  , TARGET  :: en_list
REAL                             , DIMENSION(:,:)  , POINTER :: en_ptr         ! Pointer to en_list
REAL(KIND=DBL)                                               :: cr_time        ! Time till next proton collision
REAL(KIND=DBL)                                               :: rndnum         ! Random number
REAL(KIND=DBL)                                     , TARGET  :: time_target    ! Target for time pointer
REAL(KIND=DBL)                                     , POINTER :: time           ! Current simulation time
REAL(KIND=DBL)                                               :: time_step      ! Time between abundance checks
!REAL(KIND=DBL)                                               :: time_check     ! Time used to determine ab. checks
REAL(KIND=DBL)                                               :: time_diff !DEBUGGING VAR
REAL(KIND=DBL)                                               :: t1,t2
REAL(KIND=DBL)                                               :: cpu_total
REAL(KIND=DBL)                                               :: cpu_max_time
CHARACTER(len=10)   , ALLOCATABLE, DIMENSION(:)    , TARGET  :: sp_list        !  List of species
CHARACTER(len=10)                , DIMENSION(:)    , POINTER :: sp_ptr         ! Pointer to species list
CHARACTER(len=80)                                            :: species_file   ! Name of species file
CHARACTER(len=80)                                            :: reactions_file ! Name of reactions file
CHARACTER(len=80)                                            :: hopping_file   ! File containing hopping data
TYPE (wait_info)    , ALLOCATABLE, DIMENSION(:)    , TARGET  :: wait_target
TYPE (wait_info)                 , DIMENSION(:)    , POINTER :: wait_list      !Derived data type described in chaco_data.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!! DEBUGGING/ANALYTICS VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER(KIND=LONG)                                           :: numprotons

!To enable debugging outputs, set debug to true


PRINT *, "*************************"
PRINT *, "***STARTING SIMULATION***"
PRINT *, "*************************"

! Initialize analytics and debugging vals
numprotons = 0
!numo2      = 0
!numo3      = 0
!p_e_loss   = 0D0
!disc_fluence = 0D0

time_check = 0
time_diff  = 0
cpu_total  = 0
t2         = 0
t1         = 0

! Open files
hopping_file = "hopping_data.txt"
OPEN(UNIT=1009,FILE="abundance.csv",POSITION='APPEND', STATUS='REPLACE')
!OPEN(UNIT=1011,FILE=hopping_file)
!OPEN(UNIT=1013,FILE="time_data.csv")
!Below for debugging and analytics
IF ( DEBUG .EQV. .TRUE. ) OPEN(UNIT=777,FILE='reaction_analytics.csv',STATUS='REPLACE',POSITION='APPEND')

! Nullify pointers
NULLIFY ( wait_list,anion_list,matrix_ptr,qube_ptr,sp_ptr,time,wait_len )


! Associate time value
time => time_target

! Initialize time step between abundance checks
time_step = time_total/time_counts
t1 = 0

! Initialize count_num
count_num = 1
cpu_max_time = 100.0
cpu_total = 0

species_file = 'species.dat'
reactions_file = 'reactions.dat'
OPEN (UNIT=1, FILE=species_file, STATUS='OLD', ACTION='READ', IOSTAT=err1)
OPEN (UNIT=2, FILE=reactions_file, STATUS='OLD', ACTION='READ', IOSTAT=err2)
CALL linecount(1,err1,lines_spec)
CALL linecount(2,err2,lines_react)
CLOSE(1)
CLOSE(2)

!******************************************************************************
! Create the Reaction Array
!******************************************************************************
! Note, currently ions is somewhat of a magic number
ALLOCATE( qube(lines_spec,lines_spec,3), sp_list(lines_spec), en_list(lines_spec,3), anion_target(ions) )
sp_list = "0"
anion_list => anion_target
CALL qbert(qube,lines_spec,lines_react,en_list,species_file,reactions_file,sp_list ,ions, anion_list)

PRINT *, 'The anion list is:',anion_list

! Associate the pointer to the qube
qube_ptr => qube

! Associate the pointer to en_list
en_ptr => en_list

! Associate the species list pointer
sp_ptr => sp_list


! Lookup to numbers of CRP and electron in the listj
CALL lookup( "CRP", lines_spec, sp_list, ev_nums(2))
CALL lookup( '*'  , lines_spec, sp_list, ev_nums(1) )
CALL lookup( 'e'  , lines_spec, sp_list, ev_nums(3) )

!******************************************************************************
! Calculate the dimensions of the matrix
!******************************************************************************
dimens(1) = NTHICK
dimens(2) = NEDGE
dimens(3) = dimens(2)
PRINT *, 'In main, the dimens are:',dimens

!******************************************************************************
! Create the matrix
!******************************************************************************
ALLOCATE ( matrix( dimens(1),dimens(2),dimens(3) ) )
matrix = 0

FORALL ( i=1:dimens(1),j=1:dimens(2),k=1:dimens(3), MOD(k,2) .EQ. 1 .AND. MOD(j,2) .EQ. 1)
  matrix(i,j,k) = -1
END FORALL

! Associate the pointer to the matrix
matrix_ptr => matrix

! Initialize wait list to have nothing in it
ALLOCATE( wait_target(SIZE(matrix)/3) )
wlen_target = 0
wait_len  => wlen_target
wait_list => wait_target
wait_list%wait_time = 0.0
wait_list%i         = 0
wait_list%j         = 0
wait_list%k         = 0
wait_list%sp_num    = 0
wait_list%act_type  = 0

!******************************************************************************
! Calculate initial waiting times
!******************************************************************************
time = 0.0

! Call random seed
CALL RANDOM_SEED()

! Find species number for cosmic ray
cr_num  = ev_nums(2)
rndnum  = RAND()
cr_time = -1*( DLOG(rndnum)/CR_RATE )

! Populate wait_list with cr arrival time
wait_len                      = wait_len + 1
wait_list(wait_len)%wait_time = cr_time + time
wait_list(wait_len)%sp_num    = cr_num

! Write first line in abundance.out file
!WRITE(1009,*) '  [TIME]                        ','[FLUENCE]                                    ','[O]            ','[O3]'
CALL COUNTER(numprotons,time,AB_UNIT_NUM,matrix_ptr,wait_list,4,7,wait_len)

!******************************************************************************
! Begin the simulation
!******************************************************************************
!counter = 0
DO WHILE ( time .LE. time_total )
  ! Read the top of the waiting list
  CALL roll_call( wait_list, time, wait_len, mindex )
  IF ( wait_list(mindex)%sp_num .EQ. cr_num ) THEN
    ! If the event is a proton collision, call Fallout
    numprotons = numprotons + 1
    CALL reactant_remove(wait_list,mindex,matrix_ptr,wait_len)
    CALL fallout( qube_ptr,matrix_ptr,en_ptr,anion_list,wait_list, &
                  wait_len,time,ev_nums)
    ! Calculate time to next cosmic-ray event
    rndnum  = RAND()
    cr_time = -1*( DLOG(rndnum)/cr_rate )
    ! Populate wait_list with new time
    wait_len = wait_len + 1
    wait_list(wait_len)%wait_time   = cr_time + time
    wait_list(wait_len)%sp_num   = cr_num
  ELSE
    ! If it is a regular species, decide it hopping or desorption
    SELECT CASE (wait_list(mindex)%act_type)
    CASE(1) ! The species hops
      CALL meta_hop( mindex,qube_ptr,matrix_ptr,en_ptr,wait_list,result,wait_len,time )
    CASE(2) ! The species desorbs
      matrix_ptr( wait_list(mindex)%i,wait_list(mindex)%j,wait_list(mindex)%k ) = 0
      CALL reactant_remove(wait_list,mindex,matrix_ptr,wait_len)
    CASE(3) ! The species reacts quickly
      CALL fast_reaction(mindex,wait_len,matrix_ptr,qube_ptr,en_ptr,time,wait_list)
    END SELECT
  END IF


  time_check = time_check + 1
  IF ( MOD(time_check,100) .EQ. 0 ) THEN
    CALL counter( numprotons,time, AB_UNIT_NUM, matrix_ptr,wait_list,4,7,wait_len)
    CALL CPU_TIME(t2)
    cpu_total = cpu_total + (t2-t1)
    !PRINT *, "Time =", time, "|*| Fluence =", numprotons/AREA !, "|*| Clock_diff =",t2-t1
!    WRITE(1013, *) time,",", time-time_diff,",", t2-t1,",",cpu_total
    time_diff = time
    t1 = t2
  END IF
!  IF ( MOD(time_check,1000) .EQ. 0 ) CALL counter( count_num, matrix_ptr )
END DO

!OPEN(UNIT=1013,FILE="wait_list_flaw.txt")
!DO n=1,wait_len
!  IF ( wait_list(n)%sp_num .NE. 20 ) THEN
!    IF ( matrix_ptr(wait_list(n)%i,wait_list(n)%j,wait_list(n)%k) .NE. n ) THEN
!      WRITE(1013,*) matrix_ptr(wait_list(n)%i,wait_list(n)%j,wait_list(n)%k),", ",n,",",wait_len
!    END IF
!  END IF
!END DO
!CLOSE(1013)

OPEN(UNIT=1012,FILE="final_wait_list.txt")
DO n=1,wait_len+1
  WRITE(1012,*) wait_list(n)
END DO
CLOSE(1012)

PRINT *, "****************"
PRINT *, "ENDING LOSALAMOS"
PRINT *, "****************"

CALL counter( numprotons,time, AB_UNIT_NUM,matrix_ptr,wait_list,4,7,wait_len)
CLOSE(1009)
!CLOSE(1011)
!CLOSE(1013)
IF ( DEBUG .EQV. .TRUE. ) CLOSE(777)

!PRINT *, "wait_len is: ",wait_len
!PRINT *, "Number of normal sites is:",SIZE(matrix)/3
!PRINT *, "wait_list is size ",SIZE(wait_list)
!PRINT *, 'Area is:',AREA
NULLIFY ( wait_list,anion_list,matrix_ptr,qube_ptr,sp_ptr,time,wait_len )
DEALLOCATE( qube, sp_list, en_list, matrix,wait_target )
END PROGRAM main
