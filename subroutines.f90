MODULE subroutines
  USE bsimple
  USE branchmod
  USE parameters
  USE typedefs
  USE functiondefs
  USE mc_toolbox
  USE specdata
CONTAINS
  ! *********************************************************
  ! ******* SUBROUTINES *************************************
  ! *********************************************************
  RECURSIVE SUBROUTINE lookup(name,node,id)
    !
    ! Purpose:
    !   This is a subroutine that compares a string value to values
    !  in a list and gives the index of a matching result and an
    !  error if there is no match.
    !
    !! LOOKUP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE
    INTEGER            :: id
    CHARACTER(len=10)  :: name
    TYPE(species) :: node

    IF ( TRIM(name) .EQ. TRIM(node%name) ) THEN
       id = node%id
       RETURN
    ELSE
       IF ( ASSOCIATED(node%next)) THEN
          CALL lookup(name,node%next,id)
       ELSE
          id = -1
          PRINT *, "ERROR! No match!"
          CALL EXIT()
       END IF
    END IF
    RETURN
  END SUBROUTINE lookup

  RECURSIVE SUBROUTINE writereactions(node)
    !
    ! Purpose:
    !   This is a subroutine that writes out the number of times each reaction
    !  has been called so far in the model
    !! WRITEREACTIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE
    type(reaction) :: node
    character(80) :: varfmt
    varfmt = "(I20,2ES20.4,I20)"

    OPEN(FILE="out_reactions.wsv", &
         UNIT=REACTIONS_UNIT_NUM, &
         STATUS="UNKNOWN", &
         POSITION="APPEND")
    WRITE(REACTIONS_UNIT_NUM,varfmt) node%count,TIME,FLUENCE, node%id
    node%count = 0
    IF ( ASSOCIATED(node%next)) THEN
       CALL writereactions(node%next)
    ELSE
       return
    END IF
    RETURN
  END SUBROUTINE writereactions

  SUBROUTINE hopping ( i_in, j_in, k_in, i_out, j_out, k_out, prob )
    !
    ! Purpose:
    !   The purpose of this  is to move a species from one site to another.
    !  It can move fr/ba/le/ri and up/down. Periodic boundary conditions are in
    !  place such that lateral motion moves to the other side of the lattice if
    !  it goes "overboard"
    !
    !! HOPPING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE
    !*****************
    ! Input and output
    !*****************
    INTEGER            , INTENT(IN)                          :: i_in, j_in, k_in
    INTEGER            , INTENT(OUT)                         :: i_out, j_out, k_out
    INTEGER            , INTENT(IN)                          :: prob

    ! Determine the direction of travel based on input number
    ! note that the second and third indices, j and k, are
    ! incremented by +- 2
    SELECT CASE (prob)
    CASE (1)
       ! hop back => k-2
       IF ( k_in .EQ. 1 .OR. k_in .EQ. 2 ) THEN
          IF ( MOD(DIMENS(3),2) .EQ. 1 ) THEN
             IF ( k_in .EQ. 1 ) k_out = DIMENS(3)
             IF ( k_in .EQ. 2 ) k_out = DIMENS(3)-1
          ELSE
             IF ( k_in .EQ. 1 ) k_out = DIMENS(3)-1
             IF ( k_in .EQ. 2 ) k_out = DIMENS(3)
          END IF
       ELSE
          k_out = k_in-2
       END IF
       i_out = i_in
       j_out = j_in
    CASE (2)
       ! hop forward => k+2
       IF ( k_in .EQ. DIMENS(3) .OR. k_in .EQ. DIMENS(3)-1 ) THEN
          IF ( MOD(DIMENS(3),2) .EQ. 1 ) THEN
             IF ( k_in .EQ. DIMENS(3) ) k_out = 1
             IF ( k_in .EQ. DIMENS(3)-1 ) k_out = 2
          ELSE
             IF ( k_in .EQ. DIMENS(3) ) k_out = 2
             IF ( k_in .EQ. DIMENS(3)-1) k_out = 1
          END IF
       ELSE
          k_out = k_in + 2
       END IF
       i_out = i_in
       j_out = j_in
    CASE (3)
       ! hop left => j-2
       IF ( j_in .EQ. 1 .OR. j_in .EQ. 2 ) THEN
          IF ( MOD(DIMENS(2),2) .EQ. 1 ) THEN
             IF ( j_in .EQ. 1 ) j_out = DIMENS(2)
             IF ( j_in .EQ. 2 ) j_out = DIMENS(2)-1
          ELSE
             IF ( j_in .EQ. 1 ) j_out = DIMENS(2)-1
             IF ( j_in .EQ. 2 ) j_out = DIMENS(2)
          END IF
       ELSE
          j_out = j_in-2
       END IF
       i_out = i_in
       k_out = k_in
    CASE (4)
       ! hop right => j+2
       IF ( j_in .EQ. DIMENS(2) .OR. j_in .EQ. DIMENS(2)-1 ) THEN
          IF ( MOD(DIMENS(2),2) .EQ. 1 ) THEN
             IF ( j_in .EQ. DIMENS(2) ) j_out = 1
             IF ( j_in .EQ. DIMENS(2)-1) j_out = 2
          ELSE
             IF ( j_in .EQ. DIMENS(2)) j_out = 2
             IF ( j_in .EQ. DIMENS(2)-1) j_out = 1
          END IF
       ELSE
          j_out = j_in+2
       END IF
       i_out = i_in
       k_out = k_in
    CASE (5)
       ! hop down => i+1
       ! Hopping to the monolayer above or below the current
       ! one involves +- 1 to the first DIMENSion (i)
       IF ( i_in .EQ. DIMENS(1) ) THEN
          i_out = i_in
       ELSE
          i_out = i_in+1
       END IF
       j_out = j_in
       k_out = k_in
    CASE (6)
       ! hop up => i-1
       IF ( i_in .EQ. 1 ) THEN
          i_out = i_in
       ELSE
          i_out = i_in-1
       END IF
       j_out = j_in
       k_out = k_in
    END SELECT
  END SUBROUTINE hopping

  SUBROUTINE fallout ( root,temp,prevNode,nextNode)
    ! Purpose:
    !   To calculate the track of a particle of ionizing radiation through a
    !  crystaline solid.
    !
    ! Note:
    !   The input array, sigmas, contains the proton collision cross-sections.
    !
    ! Note:
    !   sgse is short for second-generation secondary eletron.
    !
    !! FALLOUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE
    INTEGER          :: null
    INTEGER                      :: num_elecs ! number of secondary electrons pruduced
    INTEGER                      :: num_izns
    INTEGER                      :: num_exs, num_els
    INTEGER                      :: x,y,z !coordinates of cosmic-ray along track
    INTEGER                      :: step !distance the track is incremented
    INTEGER                      :: switch
    INTEGER                      :: ev_coords(3)
    INTEGER                      :: prcoords(3)
    DOUBLE PRECISION               :: p,u,rand1 ! rand num
    DOUBLE PRECISION               :: sigma_tot !total cross-section
    DOUBLE PRECISION               :: mfp ! mean free path
    DOUBLE PRECISION               :: dz ! move dist
    DOUBLE PRECISION               :: dist_trav !distance travelled since last collision
    DOUBLE PRECISION               :: e_loss,e_ion,e_exc
    DOUBLE PRECISION               :: labtheta
    DOUBLE PRECISION     , TARGET  :: e_se
    DOUBLE PRECISION   , TARGET  :: energy_target
    DOUBLE PRECISION   , POINTER :: ione
    DOUBLE PRECISION   , POINTER :: psigij(:),psigexj(:)
    CHARACTER(len=15)            :: nature
    TYPE(SIGMA_BOX)    , POINTER :: psigmas(:)
    TYPE(SE_INFO)                :: se_box
    !*************************************************************************
    !Proton cross-section data, to be phased out and replaced with a struct as
    !with se_box
    !*************************************************************************
    DOUBLE PRECISION, ALLOCATABLE, TARGET :: psigij_target(:),psigexj_target(:)
    TYPE(SIGMA_BOX),  ALLOCATABLE, TARGET :: psigmas_target(:)
    INTEGER :: thinghit
    LOGICAL :: proceed
    TYPE(node), POINTER :: root,temp,prevNode,nextNode

    IF ( DEBUG .EQV. .TRUE. ) PRINT *, '*****Starting Fallout*****'
    IF ( TRACKPLOT .EQV. .TRUE. ) THEN
       CLOSE(TRACKPLOT_UNIT_NUM)
       OPEN(UNIT=TRACKPLOT_UNIT_NUM,FILE="trackplot.csv", STATUS='REPLACE')
    END IF

    !****************************************************************************!
    ! Preliminary  calculations                                                  !
    !****************************************************************************!
    count_count = 0
    BI_CALLS = 0
    PROTON_ELOSS = 0

    ! Calculate the initial cross-sections based on the initial ion energy
    ALLOCATE(psigmas_target(3))
    psigmas => psigmas_target
    ALLOCATE(psigij_target(SIZE(o2_p_ion)))
    psigij => psigij_target
    ALLOCATE(psigexj_target(SIZE(o2_p_ex)))
    psigexj => psigexj_target
    ione => energy_target
    ione = EINIT
    psigmas%cross_section = 0D0
    psigij  = 0D0
    psigexj = 0D0
    CALL psigma_suite(ione,psigmas,psigij,psigexj)

    ! Calculate the total cross-section and the mean free path
    ! Note, sigma_i is the inelastic ionization cross-section and
    ! sigma_e is the inelastic excitation cross section
    sigma_tot = SUM(psigmas%cross_section)
    mfp       = 1./(rho*sigma_tot)

    !****************************************************************************!
    ! Determine random entry site                                                !
    !****************************************************************************!
    proceed = .FALSE.
    IF ( DEBUG .EQV. .TRUE. ) PRINT *, 'Now selecting initial site'
    DO WHILE ( proceed .EQV. .FALSE. )
       CALL RANDOM_NUMBER(p)
       CALL RANDOM_NUMBER(u)
       ! the value will be in range [1,bound]
       x = 1 + FLOOR( DIMENS(3)*p )
       y = 1 + FLOOR( DIMENS(2)*u )
       z = 1
       IF ( TRACKPLOT .EQV. .TRUE. ) THEN
          PRINT *, "Trackplot on"
          IF ( x .GT. ((DIMENS(3)/2)-(DIMENS(3)*0.1)) .AND. &
               x .LT. ((DIMENS(3)/2)+(DIMENS(3)*0.1)) .AND. &
               y .GT. ((DIMENS(2)/2)-(DIMENS(2)*0.1)) .AND. &
               y .LT. ((DIMENS(2)/2)+(DIMENS(2)*0.1)) ) proceed = .TRUE.
       ELSE
          proceed = .TRUE.
       END IF
    END DO

    IF ( DEBUG .EQV. .TRUE. ) PRINT *, "The entry site is:",x,y,z
    !****************************************************************************!
    ! Beginning of track event calculation                                       !
    !****************************************************************************!
    dist_trav = 0
    dz        = 0
    step      = 0

    ! Repeat the section below until z + step is greater (in ml) than the
    ! thickness of the ice
    count_count = 0
    num_elecs   = 0
    num_izns    = 0
    num_exs     = 0
    num_els     = 0

    main_loop: DO WHILE ((dist_trav .LE. THICK) .AND. (ione .GE. 5.0 ))
       IF ( DEBUG .EQV. .TRUE. ) PRINT *, 'Now entering loop: z=',z,&
            ' and DIMENS(1)=',DIMENS(1), &
            ' and step=',step
       count_count = count_count + 1

       ! Define event coords
       ev_coords(1) = z+step
       ev_coords(2) = y
       ev_coords(3) = x
       IF ( DEBUG .EQV. .TRUE. ) PRINT *, "The event coords in Fallout are:",ev_coords

       thinghit = MATRIX(ev_coords(1),ev_coords(2),ev_coords(3))%sp_num

       IF ( thinghit .NE. 0 ) THEN
          IF ( DEBUG .EQV. .TRUE. ) PRINT *, "The value of the matrix is:",thinghit
          ! If the site is occupied, then determine the type of event to occur
          CALL RANDOM_NUMBER(u)
          CALL RANDOM_NUMBER(rand1)
          !****************************************************************************!
          ! Determine the nature of the collision and the energy lost                  !
          !****************************************************************************!
          switch = 0
          e_loss = 0
          e_ion  = 0
          e_exc  = 0
          e_se   = 0
          ASSOCIATE ( sigma_i => psigmas(2)%cross_section, &
               sigma_e => psigmas(3)%cross_section )
            IF ( (u .GT. 0.0) .AND. (u .LE. (sigma_i + sigma_e)/sigma_tot) ) THEN
               IF ( (u .GT. 0) .AND. (u .LE. sigma_i/sigma_tot )) THEN
                  ! Ionization will occur
                  num_izns = num_izns + 1
                  switch = 2
                  CALL p_ion_select(psigij,e_ion,e_se)
                  e_loss = e_ion + e_se
                  nature = "Ionization"
                  PROTON_ELOSS = PROTON_ELOSS + e_loss
               ELSE
                  ! Excitation will occur
                  num_exs = num_exs + 1
                  switch = 1
                  CALL p_ex_select(psigexj,e_exc)
                  e_loss = e_exc
                  nature = "Excitation"
               END IF
            ELSE
               ! Elastic Collision will occur
               num_els = num_els + 1
               switch = 0
               CALL elastic_event(ione,e_loss,labtheta)
               nature = "Elastic"
            END IF
          END ASSOCIATE
          ione = ione - e_loss
          CALL psigma_suite(ione,psigmas,psigij,psigexj)

          IF ( TRACKPLOT .EQV. .TRUE. ) THEN
             count_count = count_count + 1
             WRITE(TRACKPLOT_UNIT_NUM,*) ev_coords(1),',',ev_coords(2),',',ev_coords(3),&
                  ', proton,',nature
          END IF

          IF ( switch .EQ. 0 ) THEN
             !****************************************************************************!
             ! Elastic collision                                                          !
             !****************************************************************************!
             !NB: FUTURE WORK TO ADD LATTICE DAMAGE
             CONTINUE
          ELSE IF ( switch .EQ. 1 ) THEN ! Dissociate target species on track and pla       !
             !****************************************************************************!
             ! Place excitation on site
             null = 0
             MATRIX(ev_coords(1),ev_coords(2),ev_coords(3))%sec_sp_num = EXCNUM
             prcoords = 1 ! Initialize product coordinates to 1 for recursive subroutine
             CALL new_reaction(ev_coords,ev_coords,prcoords,root,temp,prevNode,nextNode,null)
             IF ( null .EQ. 1 ) RETURN
          ELSE IF ( switch .EQ. 2 .AND. z+step .NE. 1 .AND. z+step .NE. 2 ) THEN
             !****************************************************************************!
             ! Ionization                                                                 !
             !****************************************************************************!
             temp => MATRIX(ev_coords(1),ev_coords(2),ev_coords(3))
             ! Place CRP at secondary site for reaction
             temp%sec_sp_num = CRPNUM
             ! Initialize product coords
             prcoords = 1
             ! Call new_reaction to form electron
             CALL new_reaction(ev_coords,ev_coords,prcoords,root,temp,prevNode,nextNode,null)
             IF ( DEBUG .EQV. .TRUE. ) PRINT *, temp%sp_num, temp%sec_sp_num
             ! Test to make sure that the new reaction worked properly
             IF ( DEBUG .EQV. .TRUE. ) THEN
                IF ( (.NOT. ANY(IONLIST .EQ. temp%sp_num)) .OR. &
                     (ELECNUM .NE. temp%sec_sp_num) ) THEN
                   PRINT *, "Ionization not successful!! Products not as expected!"
                   PRINT *, temp%sp_num
                   PRINT *, temp%sec_sp_num
                   CALL EXIT()
                END IF
             END IF
             ! Initialize se_box with initial energy and parent coords
             se_box%se_energy = e_se
             IF ( SECELEC .EQV. .FALSE. ) se_box%se_energy = 0D0
             se_box%parent_coords = ev_coords
             CALL se_info_init(se_box)
             ! Call new_electron
             CALL new_electron(se_box,root,temp,prevNode,nextNode)
             !Manual garbage collection
             CALL se_info_garbage(se_box)
          END IF
       END IF

       !********************!
       ! Track Plotting bit !
       !********************!
       dist_trav = dist_trav + dz
       IF ( DEBUG .EQV. .TRUE. ) PRINT *, "z=",z,"of",DIMENS(1),&
            " which is",(REAL(z)/REAL(DIMENS(1)))*100,"% of thickness"
       ! Increment z for next cycle
       z = z + step

       !********************!
       ! GOTO jumps to here !
       !********************!
       ! Call a random number between [0,1)
100    CALL RANDOM_NUMBER(p)

       ! Make sure the random number does not equal 1
       IF ( p .EQ. 1.0 ) THEN
          DO
             IF ( p .NE. 1.0 ) EXIT
             CALL RANDOM_NUMBER(p)
          END DO
       END IF

       ! Determine the travel distance
       dz = -mfp*LOG(1-p)
       dz = (dz*STEPFAC)

       ! Determine whether or not the site is occupied by dividing the
       ! Delta z by the height of the crystal cube, i.e. \Delta ml =
       ! \Delta z(m) * (1ml/c(m))
       step = INT((dz/(THICK/REAL(DIMENS(1)))))

       ! Make sure the next site is different than the previous one
       IF ( z+step .EQ. z ) THEN
          z = z + 1
       END IF

       ! Make sure the site is greater than the previous one
       IF (step .LE. 0. ) GOTO 100
       IF (z+step .GE. DIMENS(1) ) THEN
          IF ( TRACKPLOT .EQV. .TRUE. ) THEN
             IF ( (count_count .GT. TRACKMIN) .AND. (count_count .LT. TRACKMAX) ) THEN
                CALL EXIT()
             ELSE
                RETURN
             END IF
          ELSE
             RETURN
          END IF
       END IF

       IF ( DEBUG .EQV. .TRUE. ) PRINT *, 'Now at end of main loop in Fallout'
       IF ( TRACKPLOT .EQV. .TRUE. ) PRINT *, "count_count=",count_count
    END DO main_loop

    IF ( ( TRACKPLOT .EQV. .TRUE. ) .AND. (count_count .GE. TRACKMIN ) ) CALL EXIT()
    IF ( DEBUG .EQV. .TRUE. ) PRINT *, '*****Ending Fallout*****'
  END SUBROUTINE fallout

  RECURSIVE SUBROUTINE find_empty_site( in_coords,out_coords,null )
    ! Purose:
    !   This subtroutine takes some ion/bulk interaction site and finds a nearby
    !  empty site to put a second product. The return of the function is a set of
    !  coordinates.
    !
    ! Note:
    !   null = 1 => no empty sites
    !   null = 0 => empty site available
    !
    !! FIND_EMPTY_SITE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE
    !*****************
    ! Input and output
    !*****************
    INTEGER            , INTENT(IN) , DIMENSION(3)               :: in_coords !coords of reaction site
    INTEGER            , INTENT(OUT), DIMENSION(3)               :: out_coords !coords for product
    INTEGER                                          :: null !null error flag
    !****************
    ! Local variables
    !****************
    INTEGER                                                      :: i, j, k
    INTEGER                                                      :: prev(3),curr(3),next(3)
    INTEGER                                                      :: large_count
    INTEGER                                                      :: small_count
    INTEGER                                                      :: lucky !index of selected site, from rand
    INTEGER                                                      :: n !counters
    INTEGER                                                      :: i_re, j_re, k_re !in coords
    INTEGER                                                      :: i_re2,j_re2,k_re2 !out coords
    INTEGER                         , DIMENSION(6,3)             :: large_temp
    INTEGER                         , DIMENSION(4,3)             :: small_temp
    INTEGER            , ALLOCATABLE, DIMENSION(:,:)             :: temp_arr !temporary empty site array
    REAL                                                         :: rand !random number
    LOGICAL                                                      :: vacant

    ! Initialize counters and arrays
    large_count   = 0
    small_count  = 0

    ! Assign i_re, j_re, k_re to in_coords
    i_re = in_coords(1)
    j_re = in_coords(2)
    k_re = in_coords(3)

    large_temp = 0
    small_temp = 0
    scc: DO n=1,6
       layercond: IF ( i_re .EQ. 1 .AND. ( n .EQ. 5 .OR. n .EQ. 6 ) ) THEN
          ! If on top layer, stay on top layer
          CONTINUE
       ELSE IF ( i_re .EQ. DIMENS(1) .AND. n .EQ. 6 ) THEN
          ! Don't hop down if on bottom layer
          CONTINUE
       ELSE
          CALL hopping(i_re,j_re,k_re,i_re2,j_re2,k_re2,n )
          IF ( (MATRIX(i_re2,j_re2,k_re2)%sp_num .EQ. 0) .AND. &
               (MATRIX(i_re2,j_re2,k_re2)%sec_sp_num .EQ. 0) ) THEN
             large_temp(n,1)=i_re2
             large_temp(n,2)=j_re2
             large_temp(n,3)=k_re2
             large_count = large_count + 1
          END IF
       END IF layercond
    END DO scc

    largecount: IF ( large_count .GT. 0 ) THEN
       CONTINUE
    ELSE
       phantom: DO n=1,4
          ! Go to a phantom position to hop to nearest neighbors
          SELECT CASE (n)
          CASE (1)
             IF ( (j_re-1 .GT. 0) .AND. (k_re-1 .GT. 0) .AND. (k_re+1 .LE. DIMENS(3))) THEN
                CALL hopping(i_re,j_re-1,k_re+1,i_re2,j_re2,k_re2,1)
             ELSE
                GOTO 1944
             END IF

          CASE (2)
             IF ( (j_re-1 .GT. 0) .AND. (k_re-1 .GT. 0) .AND. (k_re+1 .LE. DIMENS(3))) THEN
                CALL hopping(i_re,j_re-1,k_re-1,i_re2,j_re2,k_re2,2)
             ELSE
                GOTO 1944
             END IF

          CASE (3)
             IF ( (j_re+1 .LE. DIMENS(2)) .AND. (k_re-1 .GT. 0) .AND. (k_re+1 .LE. DIMENS(3))) THEN
                CALL hopping(i_re,j_re+1,k_re-1,i_re2,j_re2,k_re2,2)
             ELSE
                GOTO 1944
             END IF

          CASE (4)
             IF ( (j_re+1 .LE. DIMENS(2)) .AND. (k_re-1 .GT. 0) .AND. (k_re+1 .LE. DIMENS(3))) THEN
                CALL hopping(i_re,j_re+1,k_re+1,i_re2,j_re2,k_re2,1)
             ELSE
                GOTO 1944
             END IF
          END SELECT

          IF ( (MATRIX(i_re2,j_re2,k_re2)%sp_num .EQ. 0) .AND. &
               (MATRIX(i_re2,j_re2,k_re2)%sec_sp_num .EQ. 0) ) THEN
             small_temp(n,1)=i_re2
             small_temp(n,2)=j_re2
             small_temp(n,3)=k_re2
             small_count = small_count + 1
          END IF
1944      CONTINUE
       END DO phantom
    END IF largecount

    ! Determine if there has been a null event
    nullevent: IF ( large_count .EQ. 0 .AND. small_count .EQ. 0 ) THEN
       checkall: IF ( FIND_EMPTY_COUNT .GT. FINDMAX ) THEN
          DO i = 1,DIMENS(3)
             DO j = 1,DIMENS(2)
                DO k = 1,DIMENS(1)
                   isempty: IF ( MATRIX(i,j,k)%sp_num .EQ. 0 ) THEN
                      out_coords(1) = i
                      out_coords(2) = j
                      out_coords(3) = k
                      null = 0
                      RETURN
                   END IF isempty
                END DO
             END DO
          END DO
          PRINT *, "Could not find site in all matrix!!!!"
          CALL EXIT()
       ELSE
          IF ( DEBUG .EQV. .TRUE.) PRINT *, "314159!"
          vacant = .TRUE.
          curr = in_coords
          next = in_coords
          movefind: DO WHILE ( vacant .EQV. .TRUE. )
             FIND_EMPTY_COUNT = FIND_EMPTY_COUNT + 1
             prev = curr
             curr = next
             CALL transport(prev,curr,next)
             CALL find_empty_site(next,out_coords,null)
             IF ( null .EQ. 0 ) vacant = .FALSE.
          END DO movefind
          RETURN
       END IF checkall
    ELSE
       FIND_EMPTY_COUNT = 0
       null = 0
    END IF nullevent

    ! populate the temp_arr such that it consists of only
    ! coordinates where there are empty spaces
    i = 1
    IF ( large_count .EQ. 0 ) THEN
       ALLOCATE( temp_arr(small_count,3) )
       temp_arr = 0
       n = 0
       DO n=1,SIZE(small_temp,1)
          IF ( ANY( small_temp(n,:) .NE. 0 ) ) THEN
             temp_arr(i,1) = small_temp(n,1)
             temp_arr(i,2) = small_temp(n,2)
             temp_arr(i,3) = small_temp(n,3)
             i = i + 1
          ELSE
             CONTINUE
          END IF
       END DO
    ELSE
       ALLOCATE( temp_arr(large_count,3) )
       temp_arr = 0
       n = 0
       DO n=1,SIZE(large_temp,1)
          IF ( ANY( large_temp(n,:) .NE. 0 ) ) THEN
             temp_arr(i,1) = large_temp(n,1)
             temp_arr(i,2) = large_temp(n,2)
             temp_arr(i,3) = large_temp(n,3)
             i = i + 1
          ELSE
             CONTINUE
          END IF
       END DO
    END IF

    ! choose one at random
    CALL RANDOM_NUMBER(rand)
    out_coords=0
    IF (SIZE(temp_arr,1) .EQ. 1) THEN
       n=0
       DO n=1,3
          out_coords(n) = temp_arr(1,n)
       END DO
    ELSE
       IF ( large_count .NE. 0 ) THEN
          lucky = INT(rand*large_count) + 1
          n=0
          DO n=1,3
             out_coords(n) = temp_arr(lucky,n)
          END DO
       ELSE
          lucky = INT(rand*small_count) + 1
          n=0
          DO n=1,3
             out_coords(n) = temp_arr(lucky,n)
          END DO
       END IF
    END IF
  END SUBROUTINE find_empty_site

  SUBROUTINE wait_calc ( temp )
    !
    ! Purpose:
    !  The purpose of this subroutine is to calculate
    ! waiting times for some species and save these
    ! values to the waiting list
    !
    ! Note: for information on rates, see the description
    !  in Chang and Herbst 2014, hereafter CH14
    !
    ! Note: As described in CH14, the rate, b, is
    !  different for surface and bulk species.
    !  For surface species, which can either
    !  desorb or diffuse, b = b_desorb + b_diff. On the
    !  other hand, for bulk species, b = b_bulkdiff.
    !
    ! Documentation:
    !  DATE          PROGRAMMER           DESCRIPTION
    !  ========      ==========           ===========
    !  20150414      C. Shingledecker     Original code
    !
    !
    !  This subroutine changes the following properties in the node
    !  1) wait_time
    !  2) hop_dir
    !  3) act_type
    !
    !! WAIT_CALC !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE
    TYPE(node), POINTER :: temp
    INTEGER             :: n
    INTEGER             :: ix,iy,iz
    DOUBLE PRECISION                                       :: rand_num  !pseudorandom number
    DOUBLE PRECISION                                       :: b_1       !surface thermal hopping rate
    DOUBLE PRECISION                                       :: b_2       !surface desorption rate
    DOUBLE PRECISION                                       :: b_3       !bulk diffusion rate
    DOUBLE PRECISION                                       :: b         !total rate, from CH14
    DOUBLE PRECISION                                       :: comp_val
    DOUBLE PRECISION                                       :: el_tmp

    ! Initialize variables
    el_tmp = 0

    ! Decide whether or not the species is on the surface
    IF ( temp%coord1 .EQ. 1 ) THEN
       DO n=1,5
          CALL hopping(temp%coord1,temp%coord2,temp%coord3,ix,iy,iz,n)
          IF ( MATRIX(ix,iy,iz)%sp_num .NE. 0 ) el_tmp = el_tmp + 0.1*EN_LIST(MATRIX(ix,iy,iz)%sp_num)
       END DO
       ! Surface species, separate rates for
       ! desorption and diffusion
       b_1 = trl_nu*EXP( -1*((EN_LIST(temp%sp_num)*E_SURF + el_tmp) / kin_temp ) )
       b_2 = trl_nu*EXP( -1*((EN_LIST(temp%sp_num) + el_tmp) / kin_temp ) )
       b = b_1 + b_2
       comp_val = b_1 / (b_1 + b_2)
       CALL RANDOM_NUMBER(rand_num)
       ! Decide whether desorption or hopping occurs
       IF ( rand_num .LT. comp_val ) THEN
          ! Diffusion occurs
          temp%act_type = 1
       ELSE
          ! Desorption occurs
          temp%act_type = 2
       END IF
    ELSE
       ! Bulk species, only bulk diffusion
       IF ( (matrix(temp%coord1,temp%coord2,temp%coord3)%normal .eqv. .true.) .and. &
            (.not. any(fast_reacts .eq. temp%sp_num)) ) THEN
          DO n=1,5
             CALL hopping(temp%coord1,temp%coord2,temp%coord3,ix,iy,iz,n)
             IF ( MATRIX(ix,iy,iz)%sp_num .NE. 0 ) el_tmp = el_tmp + 0.1*EN_LIST(MATRIX(ix,iy,iz)%sp_num)
          END DO
          b_3 = trl_nu*EXP( -1*((EN_LIST(temp%sp_num)*E_BULK + el_tmp)/ kin_temp ) )
       ELSE
          b_3 = trl_nu*EXP( -1*( EN_LIST(temp%sp_num)*E_BULK     / kin_temp ) )
       END IF
       b = b_3
       ! Only hopping (diffusion) can occur
       temp%act_type = 1
    END IF

    ! Calculate waiting time
    CALL RANDOM_NUMBER(rand_num)
    temp%wait_time = (-1*LOG(rand_num) / b) + TIME

    ! Get hopping direction
    temp%hop_dir = 1 + FLOOR((6+1-1)*rand_num)
    IF ( DEBUG .EQV. .TRUE. ) THEN
       ! Test for hopping in ranges
       IF ( (temp%hop_dir .GT. 6) .OR. (temp%hop_dir .LT. 1) ) THEN
          PRINT *, "Hopping out of ranges!"
          CALL EXIT()
       END IF
    END IF
  END SUBROUTINE wait_calc

  SUBROUTINE counter()
    !
    ! Purpose:
    !   The purpose of this subroutine is to count
    !  the number of each species in the matrix and
    !  print this information out to a file. These
    !  data shows the abundance as a function of time.
    !
    ! Documentation:
    !  DATE          PROGRAMMER           DESCRIPTION
    !  ========      ==========           ===========
    !  20150427      C. Shingledecker     Original code
    !
    !! COUNTER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE
    INTEGER                                                    :: i,j,k
    INTEGER                                                    :: o_count,o2_count,o3_count
    DOUBLE PRECISION                                             :: denom
    DOUBLE PRECISION                                             :: area
    DOUBLE PRECISION                                             :: fluence
    CHARACTER(len=80)                                          :: varfmt

    area   = EDGE*EDGE
    denom = THICK*EDGE*EDGE*1.0E20
    o_count = 0
    o2_count = 0
    o3_count = 0
    ! Method 1 of fluence calculation
    fluence  = CR_FLUX*TIME ! Note: This is the x-value for the objective function
    ! Method 2 of fluence calculation (only use 1 at a time )
    ! fluence = numprotons/area

    DO k = 1,DIMENS(3)
       DO j = 1,DIMENS(2)
          DO i = 1,DIMENS(1)
             IF ( MATRIX(i,j,k)%sp_num .EQ. O3NUM ) THEN
                o3_count = o3_count + 1
             ELSE IF ( MATRIX(i,j,k)%sp_num .EQ. O2NUM ) THEN
                o2_count = o2_count + 1
             ELSE IF ( MATRIX(i,j,k)%sp_num .EQ. ONUM ) THEN
                o_count = o_count + 1
             END IF
          END DO
       END DO
    END DO

    O_ABUNDANCE  = o_count
    O2_ABUNDANCE = o2_count
    O3_ABUNDANCE = o3_count

    IF ( NO_OUTPUT .EQV. .FALSE. ) THEN
       varfmt = "(2ES15.4,3I10)"
       WRITE(AB_UNIT_NUM,varfmt) &
            FLUENCE, & ! 1. Float64
            TIME,    & ! 3. Float64
            o2_count, & ! 4. Int64
            o_count, & ! 5. Int64
            o3_count          ! 6. Int64
    END IF

    if ( o3_analytics .eqv. .true. ) call writereactions(RE_HEAD)

    IF ( QUIET .EQV. .FALSE. ) THEN
       varfmt = "(A6,ES10.4,A9,ES10.4)"
       PRINT varfmt, " TIME=",TIME,"FLUENCE=",fluence
       varfmt = "(A5,ES10.4,A6,ES10.4)"
       PRINT varfmt, " [O]=",o_count/denom," [O3]=",o3_count/denom
       PRINT *, '***********************************************************************'
       PRINT *, "O=",O_ABUNDANCE,"O3=",O3_ABUNDANCE
       PRINT *, '***********************************************************************'
    END IF
  END SUBROUTINE counter

  SUBROUTINE transport(prev,curr,next)
    !
    ! Purpose:
    !    The purpose of this subroutine is to calculate the next step in the bulk
    !   scattering and diffusion of secondary electrons in a solid.
    !
    ! Documentation:
    !  DATE          PROGRAMMER           DESCRIPTION
    !  ========      ==========           ===========
    !  20150811      C. Shingledecker     Original code
    !
    !! TRANSPORT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE
    INTEGER            , INTENT(IN)          , DIMENSION(3)     :: prev !coordinates of previous location
    INTEGER            , INTENT(IN)          , DIMENSION(3)     :: curr !current coordinates
    INTEGER            , INTENT(OUT)         , DIMENSION(3)     :: next !coordinates of next position

    !****************
    ! Local variables
    !****************
    INTEGER                                                      :: n !counters
    REAL                                                         :: insides
    REAL                                                         :: hopdist,nextdist
    REAL                                                         :: rand,rand2 !random number
    REAL                                                         :: sigma
    LOGICAL                                                      :: carnap

    insides = (curr(2)-prev(2))**2 + (curr(3)-prev(3))**2 + (curr(1)-prev(1))**2
    hopdist = SQRT( insides  )

    sigma = 0.0
    carnap = .FALSE.
    rand = 0.0
    rand2 = 0.0
    n = 0
    next = curr
    DO WHILE ( carnap .EQV. .FALSE. )
       !    PRINT *, 'prev=',prev
       !    PRINT *, 'curr=',curr
       CALL RANDOM_NUMBER( rand )
       CALL RANDOM_NUMBER( rand2 )
       IF ( rand .LE. 0.6 ) THEN
          n = 1 + FLOOR(6*rand2)
          IF ( curr(1) .EQ. 1 .AND. ( n .EQ. 5 .OR. n .EQ. 6 ) ) THEN
             ! If on top layer, stay on top layer
             CONTINUE
          ELSE IF ( curr(1) .EQ. DIMENS(1) .AND. n .EQ. 6 ) THEN
             ! Don't hop down if on bottom layer
             CONTINUE
          ELSE
             CALL hopping(curr(1),curr(2),curr(3),next(1),next(2),next(3),n )
          END IF
       ELSE
          n = 1 + FLOOR(4*rand2)
          ! Go to a phantom position to hop to nearest neighbors
          SELECT CASE (n)
          CASE (1)
             IF ( curr(2)-1 .GT. 0          .AND. &
                  curr(3)-1 .GT. 0          .AND. &
                  curr(3)+1 .LE. DIMENS(3) ) THEN
                CALL hopping(curr(1),curr(2)-1,curr(3)+1,next(1),next(2),next(3),1)
             ELSE
                CONTINUE
             END IF
          CASE (2)
             IF ( curr(2)-1 .GT. 0          .AND. &
                  curr(3)-1 .GT. 0          .AND. &
                  curr(3)+1 .LE. DIMENS(3) ) THEN
                CALL hopping(curr(1),curr(2)-1,curr(3)-1,next(1),next(2),next(3),2)
             ELSE
                CONTINUE
             END IF
          CASE (3)
             IF ( curr(2)+1 .LE. DIMENS(2)  .AND. &
                  curr(3)-1 .GT. 0          .AND. &
                  curr(3)+1 .LE. DIMENS(3) ) THEN
                CALL hopping(curr(1),curr(2)+1,curr(3)-1,next(1),next(2),next(3),2)
             ELSE
                CONTINUE
             END IF
          CASE (4)
             IF ( curr(2)+1 .LE. DIMENS(2)  .AND. &
                  curr(3)-1 .GT. 0          .AND. &
                  curr(3)+1 .LE. DIMENS(3) ) THEN
                CALL hopping(curr(1),curr(2)+1,curr(3)+1,next(1),next(2),next(3),1)
             ELSE
                CONTINUE
             END IF
          END SELECT
       END IF
       IF ( n .NE. 5 .AND. n .NE. 6 ) THEN
          insides = (next(2)-prev(2))**2 + (next(3)-prev(3))**2
          sigma = SQRT( insides )
          insides = (next(2)-curr(2))**2 + (next(3)-curr(3))**2 + (next(1)-curr(1))**2
          nextdist = SQRT( insides )
          IF ( nextdist .NE. 0.0 ) THEN
             !        PRINT *, 'hopdist=',hopdist,'nextdist=',nextdist,'sigma=',sigma
             IF ( hopdist .LT. 2.0 .AND. sigma .GE. 2.0 ) THEN
                carnap = .TRUE.
             ELSE IF ( hopdist .GE. 2.0 .AND. sigma .GT. 3.0 ) THEN
                carnap = .TRUE.
             ELSE
                CONTINUE
             END IF
          ELSE
             CONTINUE
          END IF
       ELSE IF ( n .EQ.5 .OR. n .EQ. 6 ) THEN
          carnap = .TRUE.
       END IF
       !    PRINT *, 'next=',next
       !    PRINT *, 'hopdist=',hopdist,'sigma=',sigma,'carnap=',carnap
    END DO
    !  PRINT *, n
    RETURN
  END SUBROUTINE transport

  SUBROUTINE psigma_suite(energy,psigmas,psigij,psigexj)
    !
    ! Purpose:
    !   This subroutine is to calculate a new batch of PROTON cross-sections
    ! either at the beginning of the simulation, or after an energy-loss event,
    ! i.e. a collision.
    !
    ! Note:
    !   This subroutine calculates three cross-sections: elastic, ionization, and
    ! excitation.
    !
    ! Note:
    !   The Green-McNeal formalism is used for both proton ionization and
    ! excitation.
    !
    !! PSIGMA_SUITE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    !Data dictionary: Input parameters
    DOUBLE PRECISION, POINTER               :: energy
    TYPE(SIGMA_BOX) , POINTER, DIMENSION(:) :: psigmas
    DOUBLE PRECISION, POINTER, DIMENSION(:) :: psigij,psigexj

    !Data dictionary: Local variables
    INTEGER                                 :: n
    DOUBLE PRECISION                        :: scr_len
    DOUBLE PRECISION                        :: lss_en
    DOUBLE PRECISION                        :: massfac
    DOUBLE PRECISION                        :: sn_e
    DOUBLE PRECISION                        :: sn_eps

    !(1) Calculate elastic cross-section
    !i. calculate screening length
    scr_len = au(ZP,ZO2)
    !ii. calculate mass-factor
    massfac = mass_fac(MP,MO2)
    !iii. calculate reduced energy
    lss_en = eps(energy,ZP,ZO2,MP,MO2,scr_len)
    !iv. calculate reduced energy stopping cross-section
    sn_eps = sneps(lss_en)
    !v. calculate stopping cross-section
    sn_e = sne(sn_eps,ZP,ZO2,MP,MO2)
    !vi. calculate elastic collision cross-section
    psigmas(1)%cross_section = pelsig(energy,sn_e,massfac)
    psigmas(1)%description   = 'Elastic'

    !(2) Calculate ionization cross-section
    DO n=1,SIZE(o2_p_ion)
       ASSOCIATE( a => o2_p_ion(n)%a_epgion, &
            j => o2_p_ion(n)%j_epgion, &
            v => o2_p_ion(n)%nu_epgion, &
            o => o2_p_ion(n)%omega_epgion, &
            i => o2_p_ion(n)%i_epgion )
         psigij(n) = green_mcneal(energy,a,j,v,o,ZO2,i)
       END ASSOCIATE
    END DO
    psigmas(2)%cross_section = SUM(psigij)
    psigmas(2)%description = 'Ionization'

    !(3) Calculate excitation cross-section
    DO n=1,SIZE(o2_p_ex)
       ASSOCIATE( a => o2_p_ex(n)%a_epgex, &
            j => o2_p_ex(n)%j_epgex, &
            v => o2_p_ex(n)%nu_epgex, &
            o => o2_p_ex(n)%omega_epgex )
         psigexj(n) = green_mcneal(energy,a,j,v,o,ZO2,0D0)
       END ASSOCIATE
    END DO
    psigmas(3)%cross_section = SUM(psigexj)
    psigmas(3)%description    = 'Excitation'
  END SUBROUTINE psigma_suite

  SUBROUTINE esigma_suite(se_box)
    !
    ! Purpose:
    !   This subroutine is to calculate a set of ELECTRON cross-sections
    ! either at the beginning of the simulation, or after an energy-loss event,
    ! i.e. a collision.
    !
    ! Note:
    !   This subroutine calculates three cross-sections: elastic, ionization, and
    ! excitation.
    !
    ! Note:
    !   The ionization cross-section makes use of the formalism described in
    ! Green and Sawada (1972)
    !
    ! Note:
    !   The Porter, Jackman, and Green formalism is used for allowed transitions
    ! and the Green-Dutta (1967) formalism is used for forbidden transitions.
    !
    !! ESIGMA_SUITE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    !Data dictionary: Input parameters
    TYPE(se_info)                      :: se_box

    !Data dictionary: Local variables
    INTEGER                                 :: n
    DOUBLE PRECISION                        :: ae,ge,tnaught,tmax


    !(0) Initialize se_box
    !(1) Initialize variables
    n = 0
    ae = 0
    ge = 0
    tnaught = 0
    tmax = 0

    !(2) Calculate ionization cross-section
    DO n=1,SIZE(se_box%se_ionst)
       ASSOCIATE ( en  => se_box%se_energy           , &
            i   => se_box%se_ionst(n)%i_energy, &
            k   => se_box%se_ionst(n)%k_ion   , &
            kb  => se_box%se_ionst(n)%kb_ion  , &
            j   => se_box%se_ionst(n)%j_ion   , &
            jb  => se_box%se_ionst(n)%jb_ion  , &
            jc  => se_box%se_ionst(n)%jc_ion  , &
            gs  => se_box%se_ionst(n)%gams_ion, &
            gb  => se_box%se_ionst(n)%gamb_ion, &
            ts  => se_box%se_ionst(n)%ts_ion  , &
            ta  => se_box%se_ionst(n)%ta_ion  , &
            tb  => se_box%se_ionst(n)%tb_ion  , &
            sig => se_box%se_ionsigs(n)          )
         ! PRINT *, '#',n,'For energy:',en
         ! PRINT *, 'i=',i
         ! PRINT *, 'k=',k
         ! PRINT *, 'kb=',kb
         ! PRINT *, 'j=',j
         ! PRINT *, 'jb=',jb
         ! PRINT *, 'jc=',jc
         ! PRINT *, 'gs=',gs
         ! PRINT *, 'gb=',gb
         ! PRINT *, 'ts=',ts
         ! PRINT *, 'ta=',ta
         ! PRINT *, 'tb=',tb
         IF ( en .LT. i ) THEN
            sig = 0D0
         ELSE
            !i. Calculate the A(E) value from Green & Sawada
            ae = a_gs(en,k,kb,j,jb,jc)
            ! PRINT *, '#',n,' ae=',ae
            !ii. Calculate the \Gamma(E) factor
            ge = gamma_gs(en,gs,gb)
            ! PRINT *, '#',n,' ge=',ge
            !iii. Calculate the T_0 value
            tnaught = t_0_gs(en,ta,tb,ts)
            ! PRINT *, '#',n,' tnaught=',tnaught
            !iv. Calculate the Tmas value
            tmax = t_max_gs(en,i)
            ! PRINT *, '#',n,' tmax=',tmax
            !v. Calculate the cross-section for the state
            sig = green_sawada(ae,ge,tmax,tnaught)
         END IF
         ! PRINT *, 'the',n,' value of sig is:',sig
       END ASSOCIATE
    END DO
    !The total electron impact cross-section is the sum over the
    !cross-sections for the individual states.
    se_box%se_iontot = SUM(se_box%se_ionsigs)

    !(3) Calculate allowed excitation cross-sections
    !NB: the subroutine returns an array of values, so no
    !    loop is required
    DO n=1,SIZE(o2_e_ex_alwd)
       ASSOCIATE( e => se_box%se_energy            , &
            w => se_box%se_alwd(n)%wj_alwd   , &
            f => se_box%se_alwd(n)%fj_alwd   , &
            c => se_box%se_alwd(n)%cj_alwd   , &
            a => se_box%se_alwd(n)%alpha_alwd, &
            b => se_box%se_alwd(n)%beta_alwd    )
         se_box%se_alwdsigs(n) = pjgsigma(e,f,w,c,a,b)
         ! PRINT *, "Allowed sig #",n," =", se_box%se_alwdsigs(n)
       END ASSOCIATE
    END DO
    se_box%se_alwd_extot = SUM(se_box%se_alwdsigs)

    !(4) The total electron impact excitation is the sum of the
    !    allowed and forbidden transition cross-sections
    se_box%se_extot = se_box%se_alwd_extot !+ se_box%se_fbdn_extot

    !(5) Calculate the total cross-section as the sum of the
    ! constituent cross-sections
    se_box%se_ineltot = se_box%se_iontot + se_box%se_extot
    RETURN
  END SUBROUTINE esigma_suite

  SUBROUTINE p_ion_select(psigij,e_ion,e_se)
    !
    ! Purpose:
    !   This subroutine is to determine the specific ionization state that an
    !  inelastic collision ionizes from.
    !
    ! INPUT:
    !   An array containing the cross-sections for the distinct continuum states.
    !
    ! OUTPUT:
    !   Two energies, in eV: the ionization energy from the selected continuum
    !  state and the kinetic energy of the secondary electron.
    !
    !! P_ION_SELECT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    !Data dictionary: Input and output parameters
    DOUBLE PRECISION             , POINTER, DIMENSION(:) :: psigij
    DOUBLE PRECISION, INTENT(OUT)                        :: e_ion,e_se

    !Data dictionary: Local variables
    DOUBLE PRECISION                                     :: prob,prevprob
    DOUBLE PRECISION                                     :: sigtot
    DOUBLE PRECISION                                     :: rn
    INTEGER                                              :: n

    !(0) Initialize variables
    prob     = 0D0
    prevprob = 0d0
    sigtot   = 0D0
    rn       = 0D0
    n        = 0
    e_ion    = 0D0
    e_se     = 0D0


    !(1) Calculate the probabilities of each state based on the relative size
    !    of the cross-sections
    sigtot   = SUM(psigij)
    CALL RANDOM_NUMBER(rn)

    e_ion = 1234567d0 !Just to know what's happening for debugging
    DO n=1,SIZE(psigij)
       prob = (psigij(n)/sigtot) + prevprob
       IF ( rn .GT. prevprob .AND. rn .LT. prob ) THEN
          e_ion = o2_p_ion(n)%i_epgion
       END IF
       prevprob = prob
    END DO

    !(4) Draw another pseudo-random number, this time from a Gamma distribution
    !    to determine the kinetic energy of the low-energy electron.
    !
    !NB: The input to rgamma, aval, is a global parameter
    e_se = rgamma(AVAL)
    RETURN
  END SUBROUTINE p_ion_select

  SUBROUTINE p_ex_select(psigexj,e_exc)
    !
    ! Purpose:
    !   This subroutine is to determine the specific excited state that an
    !  inelastic collision results in the target species being promoted to.
    !
    ! Note:
    !   The formalism in Edgar, Porter, & Green (1974) is used
    !
    ! INPUT:
    !   An array containing the cross-sections for the discrete states.
    !
    ! OUTPUT:
    !   In eV: the excitation energy from the selected discrete state.
    !
    !! P_EX_SELECT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    !Data dictionary: Calling parameters
    DOUBLE PRECISION             , POINTER, DIMENSION(:) :: psigexj
    DOUBLE PRECISION, INTENT(OUT)                        :: e_exc
    !Data dictionary: Local variables
    DOUBLE PRECISION                                     :: sigtot,rn
    DOUBLE PRECISION                                     :: prevprob,prob
    INTEGER                                              :: n

    !(0) Initialize variables
    e_exc    = 0D0
    sigtot   = 0D0
    rn       = 0D0
    prevprob = 0D0
    prob     = 0D0
    n        = 0

    !(1) Calculate the probabilities of each state based on the relative size
    !    of the cross-sections
    sigtot   = SUM(psigexj)
    CALL RANDOM_NUMBER(rn)
    prevprob = 0d0
    e_exc = 1234567d0 !Just to know what's happening for debugging
    DO n=1,SIZE(psigexj)
       prob = (psigexj(n)/sigtot) + prevprob
       IF ( rn .GT. prevprob .AND. rn .LT. prob ) THEN
          e_exc = o2_p_ex(n)%w_epgex
          RETURN
       END IF
       prevprob = prob
    END DO
    RETURN
  END SUBROUTINE p_ex_select

  SUBROUTINE e_ion_select(se_box,e_ion,e_se,null)
    !
    ! Purpose:
    !   This subroutine is to determine the specific ionization state that an
    !  inelastic collision ionizes from.
    !
    !! E_ION_SELECT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    !Data dictionary: Calling parameters
    TYPE(se_info)                                  :: se_box
    DOUBLE PRECISION   , INTENT(OUT)                     :: e_ion, e_se
    INTEGER            , INTENT(OUT)                     :: null

    !Data dictionary: Local variables
    DOUBLE PRECISION                                     :: sigtot
    DOUBLE PRECISION                                     :: prob,prevprob
    DOUBLE PRECISION                                     :: rn
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)        :: arr,temparr
    INTEGER                                              :: n, arrcount,i
    INTEGER                                              :: incount

    !Ensure that there is not a null event
    null = 0
    IF ( se_box%se_ineltot .EQ. 0.0 ) THEN
       null = 1
       RETURN
    ELSE IF ( se_box%se_iontot .EQ. 0.0 ) THEN
       null = 1
       RETURN
    END IF

    !Initialize values
    e_ion = 0D0
    e_se  = 0D0

    !Initialize integers
    n        = 0
    arrcount = 0
    i        = 0
    incount  = 0

    !DETERMINE WHICH TYPE OF TRANSITION WILL OCCUR
    ALLOCATE(arr(SIZE(se_box%se_ionsigs),2))
    arr = 0
    arr(:,1) = se_box%se_ionsigs
    arr(:,2) = se_box%se_ionst%i_energy
    sigtot   = se_box%se_iontot

    IF ( SIZE(arr,1) .EQ. 1 ) THEN
       e_ion = arr(1,2)
    ELSE
       !Populate a new array with possible transitions
       DO n=1,SIZE(arr,1)
          IF ( arr(n,1) .NE. 0.0 ) arrcount = arrcount + 1
          IF ( n .EQ. SIZE(arr,1) ) THEN
             ALLOCATE(temparr(arrcount,2))
             temparr = 0
             incount = 1
             DO i=1,SIZE(arr,1)
                IF ( arr(i,1) .NE. 0.0 ) THEN
                   temparr(incount,1) = arr(i,1)
                   temparr(incount,2) = arr(i,2)
                   incount = incount + 1
                END IF
             END DO
          END IF
       END DO

       !Draw a random number and determine the precise amount of energy lost.
       CALL RANDOM_NUMBER(rn)
       prevprob = 0d0
       prob     = 0D0
       DO n=1,SIZE(temparr,1)
          prob = (temparr(n,1)/sigtot) + prevprob
          IF ( rn .GT. prevprob .AND. rn .LE. prob ) THEN
             e_ion = temparr(n,2)
             RETURN
          END IF
          prevprob = prob
       END DO
    END IF

    !Draw another pseudo-random number, this time from a Gamma distribution
    !to determine the kinetic energy of the low-energy electron.
    CALL RANDOM_NUMBER(rn)
    e_se = rn*(se_box%se_energy - e_ion)
    RETURN
  END SUBROUTINE e_ion_select

  SUBROUTINE e_ex_select(se_box,e_exc,null)
    !
    ! Purpose:
    !   This subroutine is to determine the specific excited state that an
    !  inelastic collision results in the target species being promoted to.
    !
    !! E_EX_SELECT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE
    !Data dictionary: Calling parameters
    TYPE(se_info)                                  :: se_box
    DOUBLE PRECISION, INTENT(OUT)                        :: e_exc
    INTEGER, INTENT(OUT)                                 :: null

    !Data dicitonary: Local variables
    DOUBLE PRECISION                                     :: prob,prevprob,rn
    DOUBLE PRECISION                                     :: sigtot
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)        :: arr,temparr
    INTEGER                                              :: n, arrcount,i
    INTEGER                                              :: incount

    !Initialize values
    e_exc = 0D0
    prob  = 0D0
    prevprob = 0D0
    rn       = 0D0
    n        = 0
    arrcount = 0
    i        = 0
    incount  = 0

    !Ensure that there is not a null event
    IF ( se_box%se_ineltot .EQ. 0.0 ) THEN
       null = 1
       RETURN
    ELSE IF ( se_box%se_extot .EQ. 0.0 ) THEN
       null = 1
       RETURN
    END IF

    ALLOCATE(arr(SIZE(se_box%se_alwdsigs),2))
    arr = 0
    arr(:,1) = se_box%se_alwdsigs
    arr(:,2) = se_box%se_alwd%wj_alwd
    sigtot   = se_box%se_alwd_extot

    !Populate a new array with possible transitions
    DO n=1,SIZE(arr,1)
       IF ( arr(n,1) .NE. 0.0 ) arrcount = arrcount + 1
       IF ( n .EQ. SIZE(arr,1) ) THEN
          ALLOCATE(temparr(arrcount,2))
          temparr = 0
          incount = 1
          DO i=1,SIZE(arr,1)
             IF ( arr(i,1) .NE. 0.0 ) THEN
                temparr(incount,1) = arr(i,1)
                temparr(incount,2) = arr(i,2)
                incount = incount + 1
             END IF
          END DO
       END IF
    END DO

    !Draw a random number and determine the precise amount of energy lost.
    CALL RANDOM_NUMBER(rn)
    prevprob = 0d0
    e_exc = 0d0 !Just to know what's happening for debugging
    DO n=1,SIZE(temparr,1)
       prob = (temparr(n,1)/sigtot) + prevprob
       IF ( rn .GT. prevprob .AND. rn .LE. prob ) THEN
          e_exc = temparr(n,2)
          RETURN
       END IF
       prevprob = prob
    END DO
  END SUBROUTINE e_ex_select

  SUBROUTINE elastic_event(energy,e_loss,labtheta)
    !
    ! Purpose:
    !   This subroutine is to determine the specific amount of energy lost by an
    !  ion in an elastic collisional event.
    !
    !! ELASTIC_EVENT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    ! Data dictionary: Define calling parameters
    DOUBLE PRECISION             , POINTER :: energy !ion energy in eV
    DOUBLE PRECISION, INTENT(OUT)          :: e_loss !energy lost in the collision
    DOUBLE PRECISION, INTENT(OUT)          :: labtheta !lab scattering angle

    ! Data dictionary: Define local variables
    DOUBLE PRECISION                       :: afac !screening length
    DOUBLE PRECISION                       :: gamfac !mass factor
    DOUBLE PRECISION                       :: e_lss !Lindhard-Scharff-Sigmund red. energy
    DOUBLE PRECISION                       :: bfac !reduced impact parameter
    DOUBLE PRECISION                       :: c2,s2 !cos2(cmtheta/2) and sin2(cmtheta/2)
    DOUBLE PRECISION                       :: cmtheta !center-of-mass theta
    DOUBLE PRECISION                       :: rn !random number

    afac = au(ZP,ZO2)
    gamfac = mass_fac(MP,MO2)
    e_lss = eps(energy,ZP,ZO2,MP,MO2,afac)
    CALL RANDOM_NUMBER(rn)
    bfac = b_magic(rn,afac,RHO2)
    CALL magic(e_lss,bfac,c2,s2,cmtheta)
    e_loss = t_coll(energy,gamfac,s2)
    labtheta = lab_theta(cmtheta,MP,MO2)
    RETURN
  END SUBROUTINE elastic_event

  SUBROUTINE se_info_init(se_box)
    !
    ! Purpose:
    !   This subroutine is to set up the se_info struct, which should only have
    !  an energy and parent coords assigned at the time of calling.
    !
    !! SE_INFO_INIT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    TYPE(se_info) :: se_box

    !(1) Initialize scalar values
    se_box%se_iontot     = 0
    se_box%se_extot      = 0
    se_box%se_ineltot    = 0
    se_box%se_alwd_extot = 0

    !(2) Initialize vector values
    IF ( ALLOCATED(se_box%se_ionst) .EQV. .FALSE. ) THEN
       !Initialize the ionization arrays
       ALLOCATE(se_box%se_ionst(SIZE(o2_e_ion)))
       se_box%se_ionst = o2_e_ion
       ALLOCATE(se_box%se_ionsigs(SIZE(se_box%se_ionst)))
       se_box%se_ionsigs = 0

       !Initialize allowed excitation arrays
       ALLOCATE(se_box%se_alwd(SIZE(o2_e_ex_alwd)))
       se_box%se_alwd = o2_e_ex_alwd
       ALLOCATE(se_box%se_alwdsigs(SIZE(o2_e_ex_alwd)))
       se_box%se_alwdsigs = 0
       RETURN
    ELSE
       RETURN
    END IF
  END SUBROUTINE se_info_init

  SUBROUTINE se_info_garbage(se_box)
    !
    ! Purpose:
    !   This subroutine is to garbage collect the memory used in the se_info struct
    !
    !! SE_INFO_GARBAGE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    TYPE(se_info) :: se_box

    IF ( ALLOCATED(se_box%se_ionst) .EQV. .TRUE. ) THEN
       DEALLOCATE(se_box%se_ionst)
       DEALLOCATE(se_box%se_ionsigs)
       DEALLOCATE(se_box%se_alwd)
       DEALLOCATE(se_box%se_alwdsigs)
       RETURN
    ELSE
       RETURN
    END IF
  END SUBROUTINE se_info_garbage

  SUBROUTINE init_node(root,temp_node,x,y,z)
    !
    ! Purpose
    !   This is a subroutine that compares a string value to values
    !  in a list and gives the index of a matching result and an
    !  error if there is no match.
    !
    !! INIT_NODE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE
    TYPE(node), POINTER :: root, temp_node
    INTEGER :: x,y,z

    temp_node%wait_time = 0.0
    temp_node%coord1 = x
    temp_node%coord2 = y
    temp_node%coord3 = z
    temp_node%sec_sp_num = 0
    temp_node%act_type = 0
    temp_node%hop_dir = 0
    temp_node%leftRight = 0
    temp_node%normal = .FALSE.
    temp_node%interstitial = .FALSE.

    NULLIFY(temp_node%before,temp_node%after,temp_node%parent)

    IF ( (MOD(y,2) .EQ. 1) .AND. (MOD(z,2) .EQ. 1) ) THEN
       temp_node%normal = .TRUE.
       temp_node%sp_num = 1
       CALL wait_calc(temp_node)
       CALL add_node(root,temp_node)
       ! Test for parent association
       IF ( DEBUG .EQV. .TRUE. ) THEN
          IF ( (root%coord1 .EQ. temp_node%coord1) .AND. &
               (root%coord2 .EQ. temp_node%coord2) .AND. &
               (root%coord3 .EQ. temp_node%coord3) ) THEN
             CONTINUE
          ELSE
             IF ( .NOT. ASSOCIATED(temp_node%parent)) THEN
                PRINT *, "temp_node parent not associated!"
                CALL EXIT()
             END IF
          END IF
       END IF
    ELSE IF ( (MOD(y,2) .EQ. 0) .AND. (MOD(z,2) .EQ. 0) ) THEN
       temp_node%interstitial = .true.
    END IF
  END SUBROUTINE init_node

  SUBROUTINE wipe_node(x,y,z)
    !
    ! Purpose
    !   This is a subroutine that wipes a lattice site of species specific
    !  information, while leaving species independent information intact
    !
    !! WIPE_NODE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE
    TYPE(node), POINTER :: temp
    INTEGER :: x,y,z

    temp => MATRIX(x,y,z)
    temp%wait_time = 0.0
    temp%sec_sp_num = 0
    temp%sp_num = 0
    temp%act_type = 0
    temp%hop_dir = 0
    temp%leftRight = -1
    IF ( ASSOCIATED(temp%parent) ) NULLIFY(temp%parent)
    IF ( ASSOCIATED(temp%before) ) NULLIFY(temp%before)
    IF ( ASSOCIATED(temp%after)  ) NULLIFY(temp%after)
  END SUBROUTINE wipe_node

  RECURSIVE SUBROUTINE new_reaction(i,j,k,root,temp,prevNode,nextNode,error)
    !
    ! Purpose
    !   This is a subroutine that handles reactions between two species and yields
    !  up to three products.
    !
    ! INPUT:
    !  i -> the coords of reactant 1 in matrix
    !  j -> the coords of reactant 2 in matrix (can be same as i in some cases)
    !  k -> the coords of the product to be placed in react cube
    !
    !  OUTPUT:
    !   i -> the coords of product 2 in matrix
    !   j -> the coords of product 1 in matrix
    !   k -> the coords of product 3 in matrix
    !
    !  NOTE:
    !   When this subroutine gets first called, k should be 1
    !
    !! NEW_REACTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE
    INTEGER, INTENT(INOUT)          :: i(3),j(3),k(3)
    INTEGER                         :: pr_coords(3)
    INTEGER                         :: r1, r2, pr
    INTEGER                         :: error
    TYPE(node)            , POINTER :: root, temp, prevNode, nextNode

    error = 0
    SELECT CASE (k(3))
    CASE(1)
       r1 = MATRIX(i(1),i(2),i(3))%sp_num
       r2 = MERGE(MATRIX(i(1),i(2),i(3))%sec_sp_num,&
            MATRIX(j(1),j(2),j(3))%sp_num,&
            ALL(ABS(i-j) .EQ. 0))
       k(1) = r1
       k(2) = r2
       pr = branching(k(1),k(2),k(3))
       pr_coords = j
       IF ( pr .EQ. 0 ) THEN
          ! If there can be no reaction, delete species from tree
          temp => MATRIX(i(1),i(2),i(3))
          ! Delete node but DO NOT wipe it
          CALL delete_node(root,temp,prevNode,nextNode,error)
          ! Calculate new waiting time
          CALL wait_calc(temp)
          ! Re-add node to tree
          CALL add_node(root,temp)
          error = 1
          !          IF ( DEBUG .EQV. .TRUE. ) &
          PRINT *, "Products = 0:",r1,"+",r2,"=",pr
          CALL EXIT()
          RETURN
       END IF
       k(3) = 2
    CASE(2)
       pr = branching(k(1),k(2),k(3))
       ! If product is 0, then re-add hopping species to list
       ! DO NOT do this if i=j, i.e. an excitation with 1 product
       ! has occured
       IF ( (pr .EQ. 0) .AND. (.NOT. ALL(ABS(i-j) .EQ. 0)) ) THEN
          ! If there is only one product, clear the site of R1
          ! which now becomes a lattice vacancy
          temp => MATRIX(i(1),i(2),i(3))
          CALL delete_node(root,temp,prevNode,nextNode,error)
          CALL wipe_node(i(1),i(2),i(3))
          error = 0
          RETURN
       ELSE IF (pr .EQ. 0 ) THEN
          error = 0
          RETURN
       END IF
       ! if i=j, then non-special P2 overwrites P1
       IF ( (ALL(ABS(i-j) .EQ. 0)) .AND. (.NOT. ANY(SPECIAL_LIST .EQ. pr)) ) THEN
          ! Test to see if k(2) is EXCNUM, if so, and pr!=0, get new site
          CALL find_empty_site(j,i,error) ! save out-coords to i, NOT j
       END IF
       ! For all other cases (the second product can be at site i
       ! even when i=j
       pr_coords = i ! For species
       k(3) = 3
    CASE(3)
       pr = branching(k(1),k(2),k(3))
       IF ( pr .EQ. 0 ) THEN
          error = 0
          RETURN
       END IF
       IF ( DEBUG .EQV. .TRUE. ) PRINT *, "Third product=:",pr
       CALL find_empty_site(i,pr_coords,error)
       IF ( (error .EQ. 1) .AND. (.NOT. ALL(ABS(i-j) .EQ. 0)) ) THEN
          ! If there are no empty sites around i, and i!=j, try j
          error = 0
          CALL find_empty_site(j,pr_coords,error)
          IF ( error .EQ. 1 ) THEN
             ! If there are no empty sites around either i or j, quit...
             PRINT *, "Unable to find a site for the third product!"
             CALL EXIT()
          END IF
       END IF
       ! At this points, when the subroutine quits, i=P1,j=p2,k=P3...
       IF ( DEBUG .EQV. .TRUE. ) PRINT *, "Third product coords are:",pr_coords
    END SELECT

    ! Once the case is selected, point to coords to place
    ! product
    temp => MATRIX(pr_coords(1),pr_coords(2),pr_coords(3))
    ! Place product
    IF ( pr .EQ. ELECNUM ) THEN
       ! Ensure that the electron is placed on the same site
       ! as the cation
       IF ( DEBUG .EQV. .TRUE. ) THEN
          IF ( .NOT. ALL(ABS(i-j) .EQ. 0) ) THEN
             PRINT *, "Electron not placed on same site as cation!"
             CALL EXIT()
          END IF
       END IF
       ! The rest of the node should have a species on it already
       ! just add the electron as the secondary node species
       temp%sec_sp_num = pr
    ELSE
       ! Delete whatever is there, if the site isn't empty
       IF ( temp%sp_num .NE. 0) THEN
          CALL delete_node(root,temp,prevNode,nextNode,error)
          CALL wipe_node(pr_coords(1),pr_coords(2),pr_coords(3))
          temp => MATRIX(pr_coords(1),pr_coords(2),pr_coords(3))
       END IF
       temp%sp_num = pr
       CALL wait_calc(temp)
       CALL add_node(root,temp)
    END IF

    ! Call subroutine again, and check for next product:
    ! if the next product is zero, return
    IF ( k(3) .LT. 3 ) THEN
       CALL new_reaction(i,j,k,root,temp,prevNode,nextNode,error)
    ELSE
       k = pr_coords
       RETURN
    END IF
  END SUBROUTINE new_reaction

  RECURSIVE SUBROUTINE new_electron(se_box,root,temp,prevNode,nextNode)
    ! Takes as input a "se-box" struct containing the initial electron energy
    ! and the coordinates at which it formed and does the following:
    !   1) Calculates cross-sections
    !   2) Hops the electron from inelastic collision at which the following can occur
    !      i) Electron impact ionization -> form a new electron and call self
    !     ii) Electron impact excitation -> chance of dissociating molecule
    !   3) When the electron falls below ECUTOFF, it reacts to form a negative ion
    !      which then reacts dissociatively with it's parent cation
    !
    ! NEW_ELECTRON !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE
    TYPE(se_info) :: new_se_box
    TYPE(se_info) :: se_box
    INTEGER :: n,m
    INTEGER :: curr(3),next(3),prev(3),prcoords(3)
    INTEGER :: estep
    INTEGER :: eswitch
    INTEGER :: error
    INTEGER :: cation, anion
    INTEGER :: nhops
    INTEGER :: hopCoords(6,3)
    REAL               :: erand
    DOUBLE PRECISION :: emfp
    DOUBLE PRECISION :: new_e_energy
    DOUBLE PRECISION :: ee_loss
    DOUBLE PRECISION :: de
    DOUBLE PRECISION :: e_ion, e_exc
    TYPE(node), POINTER :: root,temp,prevNode,nextNode
    LOGICAL :: vacant
    LOGICAL :: react_with_parent

    ! Initialize energy losses
    e_ion = 0.0
    e_exc = 0.0
    ee_loss = 0

    ! Initialize coordinates
    prev = 0
    curr = se_box%parent_coords
    next = curr
    estep = 0

    ! Initialize switches
    eswitch = 0

    ! Initialize real variables
    emfp = 0
    erand = 0
    de = 0
    prevNode => MATRIX(se_box%parent_coords(1),se_box%parent_coords(2),se_box%parent_coords(3))

    ! Calculate track until the electron's energy is depleted
    DO WHILE ( se_box%se_energy .GE. ECUTOFF )
       ! Calculate cross-sections for transport calculations
       CALL esigma_suite(se_box)
       ! The electron's mean-free-path is a function of the total cross sections.
       ! Note: here, we have explicitly calculated the inelastic cross section and
       ! have approximated the elastic cross section to be 1.0E-17 cm^2
       CALL RANDOM_NUMBER(erand)
       emfp = 1./(RHO*(se_box%se_ineltot+1.0E-17))
       ! Determine the actual distance travelled
       de = -1.*emfp*LOG(1.-erand)

       ! Convert to integer value
       estep = INT(de/C_PR)

       ! If the estep = 0, force it to 1
       IF ( estep .EQ. 0 ) estep = 1

       ! Hop estep number of times until the next inelastic collision
       DO n = 1,estep
          prev = curr
          curr = next
          CALL transport(prev,curr,next)
          ! Calculate approximate elastic energy loss
          ee_loss = se_box%se_energy*ELASTIC_LOSS
          se_box%se_energy = se_box%se_energy - ee_loss
          IF ( TRACKPLOT .EQV. .TRUE. ) THEN
             count_count = count_count + 1
             WRITE(TRACKPLOT_UNIT_NUM,*) curr(1),',',curr(2),',',curr(3),', electron, movement'
          END IF
       END DO

       ! Determine the nature of the inelastic event
       CALL RANDOM_NUMBER(erand)
       IF ( (erand .GT. 0) .AND. (erand .LE. se_box%se_iontot/se_box%se_ineltot) ) THEN
          eswitch = 1
       ELSE
          eswitch = 0
       END IF

       ! Based on collision type, perform action
       SELECT CASE (eswitch)
       CASE(1)
          ! Electron impact ionization
          temp => MATRIX(next(1),next(2),next(3))
          IF ((temp%sp_num .NE. 0) .AND. (.NOT. ANY(IONLIST .EQ. temp%sp_num))) THEN
             ! Calculate energy loss
             CALL e_ion_select(se_box,e_ion,new_e_energy,error)
             ! Call random number
             CALL RANDOM_NUMBER(erand)
             ! Calculate new electron energy based on \DeltaE
             !             new_e_energy = erand*(se_box%se_energy - e_ion)
             ! Energy lost is sum of ionization energy + new electron energy
             ee_loss = e_ion + new_e_energy
             ! Ionize the species and call new_electron again
             ! Place CRP pseudo reactant at site to indicate ionization
             temp%sec_sp_num = CRPNUM
             ! Initialize 3rd set of coordinates to 1
             prcoords = 1
             CALL new_reaction(next,next,prcoords,root,temp,prevNode,nextNode,error)
             temp => MATRIX(next(1),next(2),next(3))
             ! Now this site should have a cation and an electron as the secondary species
             ! at the same site
             IF ( DEBUG .EQV. .TRUE. ) THEN
                IF ( temp%sec_sp_num .NE. ELECNUM ) THEN
                   PRINT *, temp%sec_sp_num, &
                        temp%sp_num
                END IF
             END IF
             ! 1) Populate the se_box with initial energy and parent coords
             new_se_box%se_energy = new_e_energy
             new_se_box%parent_coords = next
             ! 2) Call se_info_init to initialize the arrays in the se_box
             CALL se_info_init(new_se_box)
             ! 3) Call new_electron and send it on its merry way
             CALL new_electron(new_se_box,root,temp,prevNode,nextNode)
             ! 4) When it has lost energy and reacted back with its parent cation (above)
             !    we can free-up the arrays in the se_box
             CALL se_info_garbage(new_se_box)
             !    the species at the site above should now be neutral
             IF ( DEBUG .EQV. .TRUE. ) THEN
                PRINT *, "At the end of the EII cycle, the site now has:",&
                     MATRIX(next(1),next(2),next(3))%sp_num
             END IF
          END IF
       CASE(0)
          ! Electron impact excitation
          temp => MATRIX(next(1),next(2),next(3))
          IF ((temp%sp_num .NE. 0) .AND. (.NOT. ANY(IONLIST .EQ. temp%sp_num))) THEN
             ! Place special excitation reactant at site
             temp%sec_sp_num = EXCNUM
             ! Initialize product coords to 1
             prcoords = 1
             CALL new_reaction(next,next,prcoords,root,temp,prevNode,nextNode,error)
          END IF
          CALL e_ex_select(se_box,e_exc,error)
          ! Energy lost is transition energy
          ee_loss = e_exc
       END SELECT
       ! Deduct energy lost by electron
       se_box%se_energy = se_box%se_energy - ee_loss
    END DO

    vacant = .TRUE.
    prevNode => MATRIX(se_box%parent_coords(1),se_box%parent_coords(2),se_box%parent_coords(3))
    nhops = 0
    react_with_parent = .FALSE.
    hopCoords = 0 
    curr = next
    ! Once the electron has fallen below the energy threshold, make
    ! Call transport until a non-vacant site is found
    DO n=1,6
       CALL hopping(curr(1),curr(2),curr(3),next(1),next(2),next(3),n)
       hopCoords(n,:) = next
       temp => MATRIX(next(1),next(2),next(3))
       ! Test to make sure that the species isn't some other electron's cation
       IF ( (temp%sp_num .NE. 0) .AND. &
            (.NOT. ANY(IONLIST .EQ. temp%sp_num))) THEN
          vacant = .FALSE.
          EXIT
       END IF
    END DO

    ! If no site has been found, expand the search region
    IF ( vacant .EQV. .TRUE. ) THEN
       DO n=1,6
          prev = hopCoords(n,:)
          DO m=1,6
             CALL hopping(prev(1),prev(2),prev(3),next(1),next(2),next(3),m)
             temp => MATRIX(next(1),next(2),next(3))
             ! Test to make sure that the species isn't some other electron's cation
             IF ( (temp%sp_num .NE. 0) .AND. &
                  (.NOT. ANY(IONLIST .EQ. temp%sp_num))) THEN
                vacant = .FALSE.
                EXIT
             END IF
          END DO
          IF ( vacant .EQV. .FALSE. ) EXIT
       END DO
    END IF

    ! IF still no occupied site has been found, have electron react with parent
    IF ( vacant .EQV. .TRUE. ) THEN
       react_with_parent = .TRUE.
       next = se_box%parent_coords
       temp => MATRIX(next(1),next(2),next(3))
       temp%sec_sp_num = ELECNUM
    END IF

    ! Form anion, if a suitable target was found
    IF ( react_with_parent .EQV. .FALSE. ) THEN
       ! Electron reacts to form an anion with a surrounding species
       ! i=next,j=next,form negative anion
       prcoords = 1 ! Initialize k(3)=1 for the recursive subroutine
       ! Place electron at same site to form anion
       temp%sec_sp_num = ELECNUM
       CALL new_reaction(next,next,prcoords,root,temp,prevNode,nextNode,error)
       anion = temp%sp_num
       IF ( (error .EQ. 1) .OR. (.NOT. ANY(IONLIST .EQ. anion)) ) THEN
          PRINT *, "Electron couldn't form anion!"
          PRINT *, temp%sp_num
          PRINT *, temp%sec_sp_num
          PRINT *, next
          next = se_box%parent_coords
          temp => MATRIX(next(1),next(2),next(3))
          temp%sec_sp_num = ELECNUM
       END IF
    END IF

    ! Make newly formed anion react with parent cation at %parent_coords
    ! i=next,j=parent_coords
    prcoords = 1
    prevNode => MATRIX(se_box%parent_coords(1),se_box%parent_coords(2),se_box%parent_coords(3))
    cation = prevNode%sp_num
    IF ( .NOT. ANY(IONLIST .EQ. cation)) THEN
       PRINT *, "Parent coords not a cation!!"
       CALL EXIT()
    END IF

    CALL new_reaction(next,se_box%parent_coords,prcoords,root,temp,nextNode,prevNode,error)

    IF ( error .EQ. 1) THEN
       PRINT *, "Ions couldn't recombine!"
       temp => MATRIX(next(1),next(2),next(3))
       PRINT *, temp%sp_num
       PRINT *, temp%sec_sp_num
       PRINT *, "ELECNUM is:",ELECNUM," and EXCNUM is:",EXCNUM
       PRINT *, next
       PRINT *, se_box%parent_coords
       CALL EXIT()
    END IF
  END SUBROUTINE new_electron
END MODULE subroutines
