PROGRAM main
  USE subroutines
  USE parameters
  USE typedefs
  USE functiondefs
  USE gp
  IMPLICIT NONE

  !******************************************************************************
  ! Data dictionary
  !******************************************************************************
  INTEGER                                                      :: n
  INTEGER(KIND=SHORT) , ALLOCATABLE, DIMENSION(:,:,:), TARGET  :: qube
  INTEGER             , ALLOCATABLE, DIMENSION(:)    , TARGET  :: anion_target
  INTEGER                          , DIMENSION(:)    , POINTER :: anion_list     ! List of anionic species
  INTEGER(KIND=SHORT)              , DIMENSION(:,:,:), POINTER :: qube_ptr       ! Pointer to the reaction cube
  INTEGER                                                      :: i,j,k          ! Counters
  INTEGER                                                      :: err1, err2     ! Error numbers for the files
  INTEGER(KIND=SHORT)                                          :: lines_spec     ! Number of lines in species file
  INTEGER(KIND=SHORT)                                          :: lines_react    ! Number of lines in reactions file
  INTEGER                                                      :: count_num      ! Number of abundance file
  INTEGER                                                      :: result
  INTEGER                                                      :: time_check     ! DEBUGGING VAR
  REAL                , ALLOCATABLE, DIMENSION(:)    , TARGET  :: en_list
  REAL                             , DIMENSION(:)    , POINTER :: en_ptr         ! Pointer to en_list
  REAL(KIND=DBL)                                               :: cr_time        ! Time till next proton collision
  REAL(KIND=DBL)                                               :: rndnum         ! Random number
  REAL(KIND=DBL)                                               :: time_step      ! Time between abundance checks
  !REAL(KIND=DBL)                                               :: time_check     ! Time used to determine ab. checks
  REAL(KIND=DBL)                                               :: fluence        ! Run simulation until some max fluence
  REAL(KIND=DBL)                                               :: time_diff !DEBUGGING VAR
  REAL(KIND=DBL)                                               :: t1,t2
  REAL(KIND=DBL)                                               :: cpu_total
  REAL(KIND=DBL)                                               :: cpu_max_time
  CHARACTER(len=10)   , ALLOCATABLE, DIMENSION(:)    , TARGET  :: sp_list        !  List of species
  CHARACTER(len=10)                , DIMENSION(:)    , POINTER :: sp_ptr         ! Pointer to species list
  CHARACTER(len=80)                                            :: hopping_file   ! File containing hopping data
  INTEGER                                            , TARGET  :: o3_prod_target,o3_dest_target
  INTEGER                                            , POINTER :: o3_prod,o3_dest
  INTEGER                                                      :: spec_header,reac_header
  INTEGER                                                      :: num_species,num_reacts
  LOGICAL                                                      :: not_infty
  REAL(KIND=DBL)                                               :: total_fitness   ! Total fitness
  LOGICAL                                                      :: unfit           ! TRUE if solution is too unfit -> stop simulation
  INTEGER(KIND=LONG)                                           :: numprotons
  INTEGER                                                      :: o_temp,o2_temp,o3_temp
  REAL(KIND=DBL)                                               :: temp_time
  REAL(KIND=DBL)                                               :: j2, j3
  REAL(KIND=DBL)                                               :: k12,k13
  REAL(KIND=DBL)                                               :: geminacy


  !To enable debugging outputs, set debug to true


  PRINT *, "*************************"
  PRINT *, "***STARTING SIMULATION***"
  PRINT *, "*************************"

  ! Read in constants and save seed
  CALL SYSTEM("/bin/bash pre.sh")
  CALL initconstants()
  CALL store_rand()

  ! Initialize total_fitness
  total_fitness = 0

  ! Initialize analytics and debugging vals
  numprotons = 0
  o3_prod_target = 0
  o3_dest_target = 0
  !numo2      = 0
  !numo3      = 0
  !p_e_loss   = 0D0
  !disc_fluence = 0D0

  temp_time  = 0
  time_check = 0
  time_diff  = 0
  t2         = 0
  t1         = 0


  ! Open files
  hopping_file = "hopping_data.txt"
  OPEN(UNIT=AB_UNIT_NUM,FILE="abundance.csv",POSITION='APPEND', STATUS='REPLACE')
  OPEN(UNIT=RATE_UNIT_NUM,FILE="rates.csv",POSITION='APPEND', STATUS='REPLACE')

  !OPEN(UNIT=1011,FILE=hopping_file)
  !OPEN(UNIT=1013,FILE="time_data.csv")
  !Below for debugging and analytics
  IF ( DEBUG .EQV. .TRUE. ) OPEN(UNIT=777,FILE='reaction_analytics.csv',STATUS='REPLACE',POSITION='APPEND')
  IF ( O3_ANALYTICS .EQV. .TRUE. ) OPEN(UNIT=O3_NUM,FILE='ozone_reactions.csv', STATUS='REPLACE',POSITION='APPEND')
  ! Nullify pointers
  NULLIFY ( anion_list,qube_ptr,time,o3_prod,o3_dest )
  NULLIFY( root, temp, prevNode, nextNode )

  ! Associate analytics variables
  o3_prod => o3_prod_target
  o3_dest => o3_dest_target

  ! Initialize time step between abundance checks
  time_step = time_total/time_counts
  t1 = 0

  ! Initialize count_num
  count_num = 1
  cpu_max_time = 100.0
  cpu_total = 0

  OPEN (UNIT=1, FILE=SPECIES_FILE, STATUS='OLD', ACTION='READ', IOSTAT=err1)
  OPEN (UNIT=2, FILE=REACTIONS_FILE, STATUS='OLD', ACTION='READ', IOSTAT=err2)
  CALL linecount(1,err1,lines_spec,spec_header)
  CALL linecount(2,err2,lines_react,reac_header)
  CLOSE(1)
  CLOSE(2)

  PRINT *, 'files now closed'
  num_species = lines_spec-spec_header
  num_reacts = lines_react-reac_header

  !******************************************************************************
  ! Create the Reaction Array
  !******************************************************************************
  ! Note, currently ions is somewhat of a magic number
  ALLOCATE( qube(num_species,num_species,3), sp_list(num_species), en_list(num_species), anion_target(ions) )
  sp_list = "0"
  PRINT *, 'size of species list=',SIZE(sp_list)
  anion_list => anion_target
  CALL qbert(num_species,num_reacts,qube,en_list,sp_list ,anion_list)

  PRINT *, 'The binding energy of O=',en_list(4),' and O3=',en_list(7)

  ! Associate the pointer to the qube
  qube_ptr => qube

  ! Associate the pointer to en_list
  en_ptr => en_list

  ! Associate the species list pointer
  sp_ptr => sp_list


  ! Lookup to numbers of CRP and electron in the listj
  CALL lookup( "CRP", num_species, sp_list,  CRPNUM)
  CALL lookup( '*'  , num_species, sp_list,  EXCNUM )
  CALL lookup( 'e'  , num_species, sp_list,  ELECNUM )

  SPECIAL_LIST = (/ CRPNUM, EXCNUM /)

! Initialize rate info
  ! Set all counts initially to 0
  RATEINFO%count = 0
  ! CRP + O2
  RATEINFO(1)%r1 = CRPNUM
  RATEINFO(1)%r2 = 1
  ! CRP + O3
  RATEINFO(1)%r1 = CRPNUM
  RATEINFO(1)%r2 = 7
  ! * + O2
  RATEINFO(1)%r1 = EXCNUM
  RATEINFO(1)%r2 = 1
  ! * + O3
  RATEINFO(1)%r1 = EXCNUM
  RATEINFO(1)%r2 = 7
  ! e + O2
  RATEINFO(1)%r1 = ELECNUM
  RATEINFO(1)%r2 = 1
  ! e + O3
  RATEINFO(1)%r1 = ELECNUM
  RATEINFO(1)%r2 = 7
  ! O + O2
  RATEINFO(1)%r1 = 4
  RATEINFO(1)%r2 = 1
  ! O + O3
  RATEINFO(1)%r1 = 4
  RATEINFO(1)%r2 = 7

  !******************************************************************************
  ! Calculate the dimensions of the matrix
  !******************************************************************************
  IF ( FIXED_SIZE .EQV. .TRUE. ) THEN
    dimens(1) = FIX1
    dimens(2) = FIX2
    dimens(3) = FIX3
  ELSE
    DIMENS(1) = NTHICK
    DIMENS(2) = NEDGE
    DIMENS(3) = DIMENS(2)
  END IF

  !******************************************************************************
  ! Create the matrix
  !******************************************************************************
  PRINT *, 'In main, the dimens are:',dimens
  ALLOCATE ( MATRIX( DIMENS(1),DIMENS(2),DIMENS(3) ) )

  PRINT *, "In main, matrix has been initialized to 0, now assigning -1 to O2"
  DO k=1,DIMENS(3)
    DO j=1,DIMENS(2)
      DO i=1,DIMENS(1)
        temp => MATRIX(i,j,k)
        CALL init_node(temp,i,j,k)
      END DO
    END DO
  END DO
  PRINT *, "Matrix has been fully initialized"

  !******************************************************************************
  ! Calculate initial waiting times
  !******************************************************************************
  time = 0.0

  ! Call random seed
  CALL RANDOM_SEED()

  ! Find species number for cosmic ray
  cr_num  = ev_nums(2)
  !rndnum  = RAND()
  CALL RANDOM_NUMBER(rndnum)
  cr_time = -1*( DLOG(rndnum)/CR_RATE )

  CALL COUNTER(o3_prod,o3_dest,numprotons,time,matrix_ptr,wait_list,4,7,wait_len)

  !******************************************************************************
  ! Begin the simulation
  !******************************************************************************
  !counter = 0
  fluence = time * CR_FLUX
  unfit = .FALSE.
  PRINT *, "Now  beginning loop"
  DO WHILE ( fluence .LE. FLUENCE_TOTAL .AND. .NOT. unfit)
    IF ((.NOT. ASSOCIATED(root) .OR. (TIME .LE. cr_time))) THEN
      ! Calculate time to next cosmic-ray event
      not_infty = .FALSE.
      DO WHILE ( not_infty .EQV. .FALSE. )
        !      rndnum  = RAND()
        CALL RANDOM_NUMBER(rndnum)
        cr_time = -1*( DLOG(rndnum)/cr_rate )
        IF ( cr_time + TIME .LT. 9E6 ) not_infty = .TRUE.
      END DO
      cr_time = cr_time + TIME
      ! Initialize variables for model analytics
      o_temp       = O_ABUNDANCE
      o2_temp      = O2_ABUNDANCE
      o3_temp      = O3_ABUNDANCE
      PROTON_ELOSS = 0.d0 !Reset protpn energy loss to 0
      CALL fallout( o3_prod,o3_dest,qube_ptr,en_ptr,anion_list )
      IF ( PROTON_ELOSS .GT. 0 ) numprotons = numprotons + 1
      ALTFLUENCE = numprotons/AREA
    ELSE
      CALL find_min(root, temp)
      ! If it is a regular species, decide it hopping or desorption
      SELECT CASE (temp%act_type)
      CASE(1) ! The species hops
        CALL meta_hop( o3_prod,o3_dest,qube_ptr,en_ptr,result,root,temp,prevNode,nextNode )
      CASE(2) ! The species desorbs
        x = temp%coord1
        y = temp%coord2
        z = temp%coord3
        CALL delete_node(root,temp,prevNode,nextNode,error)
        temp => matrix(x,y,z)
        CALL init_node(temp,x,y,z)
      CASE(3) ! The species reacts quickly
        CALL fast_reaction(o3_prod,o3_dest,qube_ptr,en_ptr,root,temp,prevNode,nextNode)
      END SELECT
    END IF

!    IF ( fluence .LE. 5.0E12 ) TIME_FREQ = 100000
!    IF ( fluence .GT. 5.0E12 .AND. fluence .LE. 5.0e14 ) TIME_FREQ = 10000
!    IF ( fluence .GT. 5.0E14  )  TIME_FREQ = 100000
    time_check = time_check + 1
    IF ( (MOD(time_check,TIME_FREQ) .EQ. 0) .AND. (wait_list(mindex)%sp_num .EQ. cr_num) ) THEN
      CALL counter( o3_prod,o3_dest,numprotons,time, matrix_ptr,wait_list,4,7,wait_len)

      ! Testing out the new fitness function
      CALL fitness(unfit,o3_prod,o3_dest,ALTFLUENCE,total_fitness,dimens,matrix_ptr,wait_list)

      ! Perform reaction analytics
      IF ( (FLOAT(O_ABUNDANCE-o_temp) .GT. 0) .AND. (PROTON_ELOSS .GT. 0.0) ) THEN
        geminacy = FLOAT(O_ABUNDANCE-o_temp)/PROTON_ELOSS
      ELSE
        geminacy = 0.d0
      END IF
      DELTA_TIME = time - temp_time
      ! Calculate rate-coefficient for the following reactions:
      ! (1) O + O2 -> O3
      ! This value should have units of cm^6 s^-1
      k12 = (FLOAT(RATEINFO(7)%count)/VOLUME)/DELTA_TIME
      k12 = k12/((FLOAT(O_ABUNDANCE)/VOLUME)*((FLOAT(O2_ABUNDANCE)/VOLUME)**2))
      IF ( ISNAN(k12) .OR. (k12 .GT. 1e30)) k12 = 0
      ! (2) O + O3 -> O2 + O2
      ! This value should be in units of cm^3 s^-1
      k13 = (FLOAT(RATEINFO(8)%count)/VOLUME)/DELTA_TIME
      k13 = k13/((FLOAT(O_ABUNDANCE)/VOLUME)*((FLOAT(O3_ABUNDANCE)/VOLUME)))
      IF ( ISNAN(k13) .OR. (k13 .GT. 1e30)) k13 = 0
      ! (3) X + O2 -> O + O
      ! This value is in units of s^-1
!      j2 = ABS(FLOAT(O2_ABUNDANCE-o2_temp)/(DELTA_TIME*FLOAT(O2_ABUNDANCE)))
      j2 = ABS(FLOAT(RATEINFO(1)%count + RATEINFO(3)%count + RATEINFO(5)%count)/(DELTA_TIME*FLOAT(O2_ABUNDANCE)))
      IF ( ISNAN(j2) ) j2 = 0
      ! (4) X + O3 -> O2 + O
      ! This value is in units of s^-1
!      j3 = ABS(FLOAT(O3_ABUNDANCE-o3_temp)/(DELTA_TIME*FLOAT(O3_ABUNDANCE)))
      j3 = ABS(FLOAT(RATEINFO(2)%count + RATEINFO(4)%count + RATEINFO(6)%count)/(DELTA_TIME*FLOAT(O2_ABUNDANCE)))
      IF ( ISNAN(j3) ) j3 = 0
      IF ( NO_OUTPUT .EQV. .FALSE.) WRITE(RATE_UNIT_NUM,*) geminacy,',',j2,',',j3,',',k12,',',k13
      RATEINFO%count = 0

      CALL CPU_TIME(t2)
      cpu_total = cpu_total + (t2-t1)
      time_diff = time
      t1 = t2
    END IF
    ! update fluence
    fluence = time * CR_FLUX
  END DO

  PRINT *, "****************"
  PRINT *, "ENDING LOSALAMOS"
  PRINT *, "****************"

  CALL counter( o3_prod,o3_dest,numprotons,time,matrix_ptr,wait_list,4,7,wait_len)
  CLOSE(AB_UNIT_NUM)
  !CLOSE(1011)
  !CLOSE(1013)
  IF ( DEBUG .EQV. .TRUE. ) CLOSE(777)

  CALL SYSTEM("/bin/bash post.sh")
END PROGRAM main
