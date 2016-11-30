PROGRAM main
  USE bsimple
  USE subroutines
  USE parameters
  USE typedefs
  USE functiondefs
  USE gp
  IMPLICIT NONE

  !******************************************************************************
  ! Data dictionary
  !******************************************************************************
  INTEGER             :: tmp_sp_num
  INTEGER             :: loop_count
  INTEGER             :: xx,yy,zz
  INTEGER             :: i(3),j(3),k(3)
  INTEGER             :: count_num      ! Number of abundance file
  INTEGER             :: time_check     ! DEBUGGING VAR
  INTEGER             :: error
  DOUBLE PRECISION      :: cr_time        ! Time till next proton collision
  DOUBLE PRECISION      :: rndnum         ! Random number
  DOUBLE PRECISION      :: time_step      ! Time between abundance checks
  DOUBLE PRECISION      :: time_diff !DEBUGGING VAR
  DOUBLE PRECISION      :: t1,t2
  DOUBLE PRECISION      :: cpu_total
  DOUBLE PRECISION      :: cpu_max_time
  CHARACTER(len=10)   :: chartime
  CHARACTER(len=80)   :: hopping_file   ! File containing hopping data
  INTEGER, TARGET     :: o3_prod_target,o3_dest_target
  INTEGER, POINTER    :: o3_prod,o3_dest
  LOGICAL             :: not_infty
  LOGICAL             :: cr_arrival
  DOUBLE PRECISION      :: total_fitness   ! Total fitness
  LOGICAL             :: unfit           ! TRUE if solution is too unfit -> stop simulation
  INTEGER             :: o_temp,o2_temp,o3_temp
  DOUBLE PRECISION      :: temp_time
  DOUBLE PRECISION      :: j2, j3
  DOUBLE PRECISION      :: k12,k13
  DOUBLE PRECISION      :: geminacy
  TYPE(node), POINTER :: root,temp,prevNode,nextNode

  PRINT *, "********************"
  PRINT *, "***STARTING CIRIS***"
  PRINT *, "********************"


  ! Read in constants and save seed
  CALL SYSTEM("/bin/bash pre.sh")
  CALL initconstants()
  CALL store_rand()
  SPECIAL_LIST = (/ CRPNUM, EXCNUM, ELECNUM /)

  ! Initialize total_fitness
  total_fitness = 0

  ! Initialize analytics and debugging vals
  numprotons = 0
  o3_prod_target = 0
  o3_dest_target = 0
  temp_time  = 0
  time_check = 0
  time_diff  = 0
  t2         = 0
  t1         = 0

  ! Open files
  hopping_file = "hopping_data.txt"

  OPEN(UNIT=AB_UNIT_NUM,&
       FILE="abundance.wsv",&
       POSITION='APPEND', &
       STATUS='REPLACE')

  OPEN(UNIT=RATE_UNIT_NUM,&
       FILE="rates.csv",&
       POSITION='APPEND', &
       STATUS='REPLACE')

  IF ( O3_ANALYTICS .EQV. .TRUE. ) THEN
     OPEN(UNIT=REACTIONS_UNIT_NUM,&
          FILE='ozone_reactions.csv', &
          STATUS='REPLACE',&
          POSITION='APPEND')
     CLOSE(REACTIONS_UNIT_NUM)
  END IF

  ! Nullify pointers
  NULLIFY ( o3_prod,o3_dest )
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

  !******************************************************************************
  ! Create the Reaction Array
  !******************************************************************************
  CALL qbert()

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

  DO zz=1,DIMENS(3)
     DO yy=1,DIMENS(2)
        DO xx=1,DIMENS(1)
           temp => MATRIX(xx,yy,zz)
           CALL init_node(root,temp,xx,yy,zz)
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
  CALL RANDOM_NUMBER(rndnum)
  cr_time = -1*( DLOG(rndnum)/CR_RATE )

  ! Get initial abundances
  CALL counter()

  !******************************************************************************
  ! Print fitting parameters 
  !******************************************************************************
  PRINT *, "Ed(O)=",EN_LIST(ONUM)
  PRINT *, "Ed(O3)=",EN_LIST(O3NUM)
  PRINT *, "Pdis(O2)=",O2_DISPROB
  PRINT *, "Pdis(O3)=",O3_DISPROB
  PRINT *, "Aval=",AVAL

  !******************************************************************************
  ! Begin the simulation
  !******************************************************************************
  loop_count = 0
  FLUENCE    = time * CR_FLUX
  unfit      = .FALSE.
  cr_arrival = .FALSE.
  PRINT *, "Now  beginning loop"
  DO WHILE ( FLUENCE .LE. FLUENCE_TOTAL .AND. .NOT. unfit)
     loop_count = loop_count + 1
     ! At the start of the simulation, or whenever it's time for a particle
     CALL find_min(root, temp)

     IF ( temp%wait_time .LT. TIME ) THEN
        PRINT *, "ERROR!!!!! temptime < TIME!!!!!"
        CALL EXIT()
     END IF

     IF ( temp%sp_num .EQ. ONUM) THEN
        N_OHOP = N_OHOP + 1
     ElSE IF ( temp%sp_num .EQ. O3NUM ) THEN
        N_O3HOP = N_O3HOP + 1
     END IF

     IF (temp%wait_time .GT. cr_time) THEN
        NUMPROTONS = NUMPROTONS + 1
        time_check = time_check + 1
        cr_arrival = .TRUE.
        TIME = cr_time
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

        ! Calculate track/damage
        CALL fallout( root,temp,prevNode,nextNode )
        IF ( PROTON_ELOSS .GT. 0 ) numprotons = numprotons + 1
     ELSE
        IF ( FIX_FREQ .EQV. .TRUE. ) time_check = time_check + 1
        TIME = temp%wait_time
        !        PRINT *, "Min time is:",temp%wait_time,"Min species is:",temp%sp_num
        i(1) = temp%coord1
        i(2) = temp%coord2
        i(3) = temp%coord3
        ! If it is a regular species, decide it hopping or desorption
        SELECT CASE (temp%act_type)
        CASE(1) ! The species hops
           ! Call hopping to get new coords
           CALL hopping(temp%coord1,temp%coord2,temp%coord3,&
                j(1),j(2),j(3),temp%hop_dir)
           nextNode => MATRIX(j(1),j(2),j(3))
           ! If the site is empty, move product there
           IF ( (nextNode%sp_num .EQ. 0) .OR. (ALL(ABS(j-i) .EQ. 0)) ) THEN
              ! Save sp_num
              tmp_sp_num = temp%sp_num
              CALL delete_node(root,temp,prevNode,nextNode,error)
              CALL wipe_node(i(1),i(2),i(3))
              ! Point temp to new location
              temp => MATRIX(j(1),j(2),j(3))
              ! Update species number
              temp%sp_num = tmp_sp_num
              ! Get new hopping time
              CALL wait_calc(temp)
              ! Add back to tree
              CALL add_node(root,temp)
           ELSE
              ! If the site isn't empty, call new_reaction
              ! First, initialize k to 1
              ! i => the original location of hopping species
              ! j => the site to which the species is hopping
              k = 1
              CALL new_reaction(i,j,k,root,temp,prevNode,nextNode,error)
           END IF
        CASE(2) ! The species desorbs
           CALL delete_node(root,temp,prevNode,nextNode,error)
           CALL wipe_node(i(1),i(2),i(3))
           IF ( DEBUG .EQV. .TRUE. ) THEN
              ! Test to make sure temp still points to the
              ! x,y,z coords of the matrix
              IF ( (temp%coord1 .EQ. i(1)) .AND. &
                   (temp%coord2 .EQ. i(2)) .AND. &
                   (temp%coord3 .EQ. i(3))) THEN
                 CONTINUE
              ELSE
                 PRINT *, "Temp no longer points to i!"
                 CALL EXIT()
              END IF
           END IF
        END SELECT
     END IF

     ! update FLUENCE
     FLUENCE = TIME*CR_FLUX
!     FLUENCE = DBLE(NUMPROTONS)/AREA

     ! If event this loop is a collision...
     crarrive: IF ( (cr_arrival .EQV. .TRUE.) .OR. (FIX_FREQ .EQV. .TRUE. ) )THEN
        CALL date_and_time(TIME=chartime)
        READ(chartime,*) t2
        cpu_total = cpu_total + (t2-t1)
        time_diff = t2-t1
        cr_arrival = .FALSE.

        IF ( FIX_FREQ .EQV. .FALSE. ) THEN
           ! Set time_freq
           IF ( FLUENCE .GT. 1.0d13 .AND. FLUENCE .LT. 1.0d14 ) THEN 
              TIME_FREQ = 10 
           ELSE IF ( FLUENCE .GT. 1.0d14 .AND. FLUENCE .LT. 1.0d15) THEN 
              TIME_FREQ = 100
           ELSE IF ( FLUENCE .GT. 1.0d15 .AND. FLUENCE .LT. 1.0d16) THEN 
              TIME_FREQ = 1000
           ELSE IF ( FLUENCE .GT. 1.0d16 ) THEN 
              TIME_FREQ = 10000
           END IF
        END IF

        checktime: IF ( MOD(time_check,TIME_FREQ) .EQ. 0 ) THEN
           CALL counter()
           CALL fitness(unfit,FLUENCE,total_fitness)
           t1 = t2
           N_OHOP = 0
           N_O3HOP = 0
        END IF checktime

        ratecalc: IF ( CALC_RATES .EQV. .TRUE. ) THEN
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
           j2 = ABS(FLOAT(RATEINFO(1)%count + RATEINFO(3)%count + RATEINFO(5)%count)/&
                (DELTA_TIME*FLOAT(O2_ABUNDANCE)))
           IF ( ISNAN(j2) ) j2 = 0
           ! (4) X + O3 -> O2 + O
           ! This value is in units of s^-1
           !      j3 = ABS(FLOAT(O3_ABUNDANCE-o3_temp)/(DELTA_TIME*FLOAT(O3_ABUNDANCE)))
           j3 = ABS(FLOAT(RATEINFO(2)%count + RATEINFO(4)%count + RATEINFO(6)%count)/&
                (DELTA_TIME*FLOAT(O2_ABUNDANCE)))
           IF ( ISNAN(j3) ) j3 = 0
           IF ( NO_OUTPUT .EQV. .FALSE.) WRITE(RATE_UNIT_NUM,*) geminacy,',',j2,',',j3,',',k12,',',k13
           RATEINFO%count = 0
        END IF ratecalc

     END IF crarrive
  END DO

  PRINT *, "************"
  PRINT *, "ENDING CIRIS"
  PRINT *, "************"

  CALL counter()
  CLOSE(AB_UNIT_NUM)
  !CLOSE(1011)
  !CLOSE(1013)
  IF ( DEBUG .EQV. .TRUE. ) CLOSE(777)

  CALL SYSTEM("/bin/bash post.sh")
END PROGRAM main
