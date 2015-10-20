PROGRAM moonbase 
  USE chaco_data
  USE parameters
  USE typedefs
  USE functiondefs

IMPLICIT NONE

!******************************************************************************
! Data dictionary
!******************************************************************************
INTEGER             , ALLOCATABLE, DIMENSION(:,:,:), TARGET  :: matrix         ! Ice-mantle matrix
INTEGER             , ALLOCATABLE, DIMENSION(:)    , TARGET  :: mobile_list    ! List of mobile species
INTEGER(KIND=SHORT) , ALLOCATABLE, DIMENSION(:,:,:), TARGET  :: qube
INTEGER                          , DIMENSION(:)    , POINTER :: mobile_ptr     ! Pointer of mobile species 
INTEGER             , ALLOCATABLE, DIMENSION(:)    , TARGET  :: anion_target
INTEGER                          , DIMENSION(:)    , POINTER :: anion_list     ! List of anionic species
INTEGER                          , DIMENSION(:,:,:), POINTER :: matrix_ptr     ! Pointer to the matrix
INTEGER(KIND=SHORT)              , DIMENSION(:,:,:), POINTER :: qube_ptr       ! Pointer to the reaction cube
INTEGER                                            , TARGET  :: wlen_target
INTEGER                                            , POINTER :: wait_len       ! Length of nonzero waitlist elements
INTEGER(KIND=SHORT)              , DIMENSION(3)              :: ev_nums
INTEGER                                                      :: mindex
INTEGER                          , DIMENSION(3)              :: dimens         ! Dimensions of the matrix
INTEGER                                                      :: n,i,j,k          ! Counters
INTEGER                                                      :: err1, err2     ! Error numbers for the files
INTEGER                                                      :: cr_num         ! species number for cosmic rays
INTEGER(KIND=SHORT)                                          :: lines_spec     ! Number of lines in species file
INTEGER(KIND=SHORT)                                          :: lines_react    ! Number of lines in reactions file
INTEGER                                                      :: count_num      ! Number of abundance file
INTEGER                                                      :: result
INTEGER                                                      :: time_check     ! DEBUGGING VAR
REAL                , ALLOCATABLE, DIMENSION(:,:)  , TARGET  :: en_list
REAL                             , DIMENSION(:,:)  , POINTER :: en_ptr         ! Pointer to en_list
REAL(KIND=DBL)                                               :: cr_rate        ! Rate of proton arrival
REAL(KIND=DBL)                                               :: cr_time        ! Time till next proton collision
REAL(KIND=DBL)                                               :: rand           ! Random number
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
!*****************************Cross-Sections************************************
!TYPE (sigma_box)    , ALLOCATABLE, DIMENSION(:)    , TARGET  :: psigmas_target ! Proton cross-sections
!TYPE (sigma_box)                 , DIMENSION(:)    , POINTER :: psigmas
!DOUBLE PRECISION    , ALLOCATABLE, DIMENSION(:)    , TARGET  :: psigij_target
!DOUBLE PRECISION                 , DIMENSION(:)    , POINTER :: psigij
!DOUBLE PRECISION    , ALLOCATABLE, DIMENSION(:)    , TARGET  :: psigexj_target !individual state cross-sections
!DOUBLE PRECISION                 , DIMENSION(:)    , POINTER :: psigexj
!DOUBLE PRECISION                                   , POINTER :: ionenergy
!DOUBLE PRECISION                                   , TARGET  :: energy_target

PRINT *, "here we go???"




time_check = 0
time_diff = 0
cpu_total = 0
t2 = 0
t1 = 0

! Open files 
hopping_file = "hopping_data.txt"
OPEN(UNIT=1009,FILE="abundance.csv",POSITION='APPEND', STATUS='REPLACE')
OPEN(UNIT=1011,FILE=hopping_file)
OPEN(UNIT=1013,FILE="time_data.csv")


! Nullify pointers
NULLIFY ( wait_list,mobile_ptr,anion_list,matrix_ptr,qube_ptr,sp_ptr,time,wait_len ) 


! Allocate and associate mobile list and ptr
ALLOCATE ( mobile_list(1) ) 
mobile_list = (/ 4 /)
!ALLOCATE( mobile_ptr(SIZE(mobile_list)) )
mobile_ptr => mobile_list

! Associate time value
time => time_target

! Initialize time step between abundance checks
time_step = time_total/time_counts
t1 = 0

! Initialize count_num
count_num = 1
cpu_max_time = 100.0 
cpu_total = 0

species_file = 'species.d'                                                   
reactions_file = 'reactions.d'                                             
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

! Associate the pointer to the qube 
!ALLOCATE( qube_ptr(SIZE(qube,1),SIZE(qube,2),SIZE(qube,3) ) )
qube_ptr => qube

! Associate the pointer to en_list
!ALLOCATE( en_ptr(SIZE(en_list,1),SIZE(en_list,2)) )
en_ptr => en_list

! Associate the species list pointer
!ALLOCATE( sp_ptr(lines_spec) )
sp_ptr => sp_list


! Lookup to numbers of CRP and electron in the listj
CALL lookup( "CRP", lines_spec, sp_list, ev_nums(2))  
CALL lookup( '*'  , lines_spec, sp_list, ev_nums(1) ) 
CALL lookup( 'e'  , lines_spec, sp_list, ev_nums(3) )
 
!******************************************************************************
! Calculate the dimensions of the matrix 
!******************************************************************************

dimens(1) = FLOOR(THICK/C_PR)
dimens(2) = FLOOR(EDGE/BDIM)
dimens(3) = dimens(2) !FLOOR(edge/a)
!******************************************************************************
! Create the matrix 
!******************************************************************************

ALLOCATE ( matrix( dimens(1),dimens(2),dimens(3) ) )
matrix = 0

FORALL ( i=1:dimens(1),j=1:dimens(2),k=1:dimens(3), MOD(k,2) .EQ. 1 .AND. MOD(j,2) .EQ. 1) 
  matrix(i,j,k) = -1
END FORALL

! Associate the pointer to the matrix
!ALLOCATE( matrix_ptr(dimens(1),dimens(2),dimens(3)) )
matrix_ptr => matrix

! Initialize wait list to have nothing in it
ALLOCATE( wait_target(SIZE(matrix)/3) )
!ALLOCATE( wait_list(SIZE(matrix)/3) ) 
wlen_target = 0
wait_len => wlen_target
wait_list => wait_target
wait_list%wait_time = 0.0

wait_list%i = 0
wait_list%j = 0
wait_list%k = 0
wait_list%sp_num = 0
wait_list%act_type = 0
!******************************************************************************
! Calculate initial waiting times 
!******************************************************************************

! Initialize time
time = 0.0

! Call random seed
CALL RANDOM_SEED()

! Find species number for cosmic ray
cr_num = ev_nums(2)

cr_rate = cr_flux*area
CALL RANDOM_NUMBER(rand)
cr_time = -1*( LOG(rand)/cr_rate )

! Populate wait_list with cr arrival time
wait_len = wait_len + 1
wait_list(wait_len)%wait_time   = cr_time + time
wait_list(wait_len)%sp_num   = cr_num

! Write first line in abundance.out file
!WRITE(1009,*) '  [TIME]                        ','[FLUENCE]                                    ','[O]            ','[O3]'
CALL COUNTER(time,AB_UNIT_NUM,matrix_ptr,wait_list,4,7)


!******************************************************************************
! Begin the simulation 
!******************************************************************************

! Initialize time to 0

! To test, call Fallout once
!CALL fallout( qube_ptr,matrix_ptr,en_ptr,anion_list,mobile_ptr,wait_list,wait_len,time, ev_nums )
! Write wait_list data to file
!OPEN(UNIT=1012,FILE="initial_wait_list.txt")
!DO n=1,wait_len+1
!  WRITE(1012,*) wait_list(n)
!END DO
!CLOSE(1012)

 
! Begin simulation
!counter = 0
DO WHILE ( time .LE. time_total )
  ! Read the top of the waiting list
  CALL roll_call( wait_list, time, wait_len, mindex )
  IF ( wait_list(mindex)%sp_num .EQ. cr_num ) THEN
!    ! If the event is a proton collision, call Fallout
    CALL reactant_remove(wait_list,mindex,matrix_ptr,wait_len)
    CALL fallout( qube_ptr,matrix_ptr,en_ptr,anion_list,mobile_ptr,wait_list, &
                  wait_len,time,ev_nums)
    ! Calculate time to next cosmic-ray event 
    CALL RANDOM_NUMBER(rand)
    cr_time = -1*( LOG(rand)/cr_rate )
    ! Populate wait_list with new time
    wait_len = wait_len + 1
    wait_list(wait_len)%wait_time   = cr_time + time
    wait_list(wait_len)%sp_num   = cr_num
  ELSE 
    ! If it is a regular species, decide it hopping or desorption
    IF ( wait_list(mindex)%act_type .EQ. 1 ) THEN
      ! The species hops
      CALL meta_hop( mindex,qube_ptr,matrix_ptr,en_ptr,wait_list,mobile_ptr,result,wait_len,time )
    ELSE
      ! The species desorbs
      matrix_ptr( wait_list(mindex)%i,wait_list(mindex)%j,wait_list(mindex)%k ) = 0
      CALL reactant_remove(wait_list,mindex,matrix_ptr,wait_len)  
    END IF
  END IF 


  time_check = time_check + 1  
  IF ( MOD(time_check,100) .EQ. 0 ) THEN
    CALL counter( time, AB_UNIT_NUM, matrix_ptr,wait_list,4,7)
    CALL CPU_TIME(t2)
    cpu_total = cpu_total + (t2-t1)
    PRINT *, "Time =", time, "|*| cpu_total =", cpu_total !, "|*| Clock_diff =",t2-t1 
!    WRITE(1013, *) time,",", time-time_diff,",", t2-t1,",",cpu_total 
    time_diff = time
    t1 = t2
  END IF
!  IF ( MOD(time_check,1000000) .EQ. 0 ) CALL counter( count_num, matrix_ptr )
END DO

OPEN(UNIT=1013,FILE="wait_list_flaw.txt")
DO n=1,wait_len
  IF ( wait_list(n)%sp_num .NE. 20 ) THEN
    IF ( matrix_ptr(wait_list(n)%i,wait_list(n)%j,wait_list(n)%k) .NE. n ) THEN
      WRITE(1013,*) matrix_ptr(wait_list(n)%i,wait_list(n)%j,wait_list(n)%k),", ",n,",",wait_len
    END IF
  END IF
END DO
CLOSE(1013)


!OPEN(UNIT=1012,FILE="final_wait_list.txt")
!DO n=1,wait_len+1
!  WRITE(1012,*) wait_list(n)
!END DO
!CLOSE(1012)



PRINT *, "LEAVING THE MOONBASE!!!!!!"
PRINT *, "3"
PRINT *, "2"
PRINT *, "1"
PRINT *, "BLASTOFF!!!!!!!"
CALL counter( time, AB_UNIT_NUM,matrix_ptr,wait_list,4,7)
CLOSE(1009)
CLOSE(1011)
CLOSE(1013)

!CALL chess(10,matrix_ptr)
PRINT *, "wait_len is: ",wait_len
PRINT *, "matrix is size ",SIZEOF(matrix)
PRINT *, "wait_list is size ",SIZEOF(wait_list)
!NULLIFY(wait_list,matrix_ptr)
!DEALLOCATE(wait_target,matrix_ptr)
NULLIFY ( wait_list,mobile_ptr,anion_list,matrix_ptr,qube_ptr,sp_ptr,time,wait_len ) 
DEALLOCATE( qube, sp_list, en_list, matrix,mobile_list,wait_target )
END PROGRAM moonbase
