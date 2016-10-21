MODULE subroutines
  USE parameters
  USE typedefs
  USE functiondefs
  USE mc_toolbox
  USE specdata


CONTAINS
  ! *********************************************************
  ! ******* SUBROUTINES *************************************
  ! *********************************************************
  SUBROUTINE lookup(string, nlines, array, n)
    !
    ! Purpose:
    !   This is a subroutine that compares a string value to values
    !  in a list and gives the index of a matching result and an
    !  error if there is no match.
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! LOOKUP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    !*****************
    ! Input and output
    !*****************

    INTEGER            , INTENT(IN)                    :: nlines
    INTEGER            , INTENT(OUT)                   :: n
    CHARACTER(*)       , INTENT(IN)                    :: string
    CHARACTER(len=10)  , INTENT(IN), DIMENSION(nlines) :: array

    !****************
    ! Local variables
    !****************

    INTEGER(KIND=SHORT)                                :: i
    CHARACTER(len=10)              , DIMENSION(1)      :: string_arr

    ! Go through the array and compare the supplied string with the
    ! strings in the array
    n = 0
    string_arr = (/ string /)
    DO i=1,nlines
       IF ( TRIM(string_arr(1)) .EQ. TRIM(array(i)) )THEN
          n=i
       END IF
    END DO
  END SUBROUTINE lookup

  SUBROUTINE linecount(unitnum, errcode, lines, header_num)
    !
    ! Purpose:
    !   This subroutine simply counts the number of lines in a file.
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! LINECOUNT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    !*****************
    ! Input and output
    !*****************
    INTEGER            , INTENT(IN)  :: unitnum
    INTEGER            , INTENT(OUT) :: errcode
    INTEGER(KIND=SHORT), INTENT(OUT) :: lines
    INTEGER(KIND=SHORT)              :: adv
    CHARACTER(LEN=100)               :: line
    INTEGER                          :: header_num

    adv = 1
    lines = 0
    header_num = 0
    DO
       READ(unitnum,*,IOSTAT=errcode) line
       !        PRINT *, line
       IF ( errcode .NE. 0 ) EXIT
       IF ( line(1:1) .EQ. '!' ) header_num = header_num + 1
       lines = lines + adv
    END DO
    REWIND(unitnum)
  END SUBROUTINE linecount

  SUBROUTINE hopping ( i_in, j_in, k_in, i_out, j_out, k_out, prob, dimens )
    !
    ! Purpose:
    !   The purpose of this  is to move a species from one site to another.
    !  It can move fr/ba/le/ri and up/down. Periodic boundary conditions are in
    !  place such that lateral motion moves to the other side of the lattice if
    !  it goes "overboard"
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! HOPPING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    !*****************
    ! Input and output
    !*****************

    INTEGER            , INTENT(IN)                          :: i_in, j_in, k_in
    INTEGER            , INTENT(OUT)                         :: i_out, j_out, k_out
    INTEGER            , INTENT(IN)                          :: prob
    INTEGER            , INTENT(IN), DIMENSION(3)            :: dimens !dimension of matrix


    ! Determine the direction of travel based on input number
    ! note that the second and third indices, j and k, are
    ! incremented by +- 2
    SELECT CASE (prob)


    CASE (1)
       ! hop back => k-2
       IF ( k_in .EQ. 1 .OR. k_in .EQ. 2 ) THEN
          IF ( MOD(dimens(3),2) .EQ. 1 ) THEN
             IF ( k_in .EQ. 1 ) k_out = dimens(3)
             IF ( k_in .EQ. 2 ) k_out = dimens(3)-1
          ELSE
             IF ( k_in .EQ. 1 ) k_out = dimens(3)-1
             IF ( k_in .EQ. 2 ) k_out = dimens(3)
          END IF
       ELSE
          k_out = k_in-2
       END IF
       i_out = i_in
       j_out = j_in

    CASE (2)
       ! hop forward => k+2
       IF ( k_in .EQ. dimens(3) .OR. k_in .EQ. dimens(3)-1 ) THEN
          IF ( MOD(dimens(3),2) .EQ. 1 ) THEN
             IF ( k_in .EQ. dimens(3) ) k_out = 1
             IF ( k_in .EQ. dimens(3)-1 ) k_out = 2
          ELSE
             IF ( k_in .EQ. dimens(3) ) k_out = 2
             IF ( k_in .EQ. dimens(3)-1) k_out = 1
          END IF
       ELSE
          k_out = k_in + 2
       END IF
       i_out = i_in
       j_out = j_in

    CASE (3)
       ! hop left => j-2
       IF ( j_in .EQ. 1 .OR. j_in .EQ. 2 ) THEN
          IF ( MOD(dimens(2),2) .EQ. 1 ) THEN
             IF ( j_in .EQ. 1 ) j_out = dimens(2)
             IF ( j_in .EQ. 2 ) j_out = dimens(2)-1
          ELSE
             IF ( j_in .EQ. 1 ) j_out = dimens(2)-1
             IF ( j_in .EQ. 2 ) j_out = dimens(2)
          END IF
       ELSE
          j_out = j_in-2
       END IF
       i_out = i_in
       k_out = k_in

    CASE (4)
       ! hop right => j+2
       IF ( j_in .EQ. dimens(2) .OR. j_in .EQ. dimens(2)-1 ) THEN
          IF ( MOD(dimens(2),2) .EQ. 1 ) THEN
             IF ( j_in .EQ. dimens(2) ) j_out = 1
             IF ( j_in .EQ. dimens(2)-1) j_out = 2
          ELSE
             IF ( j_in .EQ. dimens(2)) j_out = 2
             IF ( j_in .EQ. dimens(2)-1) j_out = 1
          END IF
       ELSE
          j_out = j_in+2
       END IF
       i_out = i_in
       k_out = k_in


    CASE (5)
       ! hop down => i+1
       ! Hopping to the monolayer above or below the current
       ! one involves +- 1 to the first dimension (i)
       IF ( i_in .EQ. dimens(1) ) THEN
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

  SUBROUTINE lookaroundyou ( react_cube, coords, null, small_count, &
       large_count, small_arr, large_arr, wait_list )
    !
    ! Purpose:
    !   This subroutine is designed to look at the surrounding spaces in a
    !  matrix and determine if any of them are possible co-reactants for
    !  any reactions in the network used. If there is a match, null=0,
    !  if there is no match found, i.e. no reacting partners, null=1.
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! LOOKAROUNDYOU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IMPLICIT NONE

    !*****************
    ! Input and output
    !*****************

    INTEGER(KIND=SHORT), INTENT(OUT)                             :: null
    INTEGER            , INTENT(OUT)                             :: small_count
    INTEGER            , INTENT(OUT)                             :: large_count
    INTEGER(KIND=SHORT)             , DIMENSION(:,:,:), POINTER  :: react_cube
    INTEGER            , INTENT(OUT), DIMENSION(6,4)             :: large_arr
    INTEGER            , INTENT(OUT), DIMENSION(4,4)             :: small_arr
    INTEGER            , INTENT(IN) , DIMENSION(3)               :: coords
    INTEGER                         , DIMENSION(3)               :: new_coords
    TYPE(wait_info)                 , DIMENSION(:)    , POINTER  :: wait_list


    !****************
    ! Local variables
    !****************
    INTEGER                                                      :: r1, r2
    INTEGER                                                      :: i_re,j_re,k_re
    INTEGER                                                      :: i_re2,j_re2,k_re2
    INTEGER                                                      :: n


    ! Initialize coordinates
    i_re = coords(1)
    j_re = coords(2)
    k_re = coords(3)

    ! Initialize counters
    large_count  = 0
    small_count = 0

    ! Initialize arrays
    large_arr = 0
    small_arr = 0

    ! If the site is on the top or bottom layers, limit the hopping
    n=0
    initial_loop: DO n=1,6
       top_bottom: IF ( i_re .EQ. 1 .AND. ( n .EQ. 5 .OR. n .EQ. 6 ) ) THEN
          ! If on top layer, stay on top layer
          CONTINUE
       ELSE IF ( i_re .EQ. dimens(1) .AND. n .EQ. 6 ) THEN
          ! Don't hop down if on bottom layer
          CONTINUE
       ELSE
          CALL hopping(i_re,j_re,k_re,i_re2,j_re2,k_re2,n,dimens)
          IF ( DEBUG .EQV. .TRUE. ) THEN
             PRINT *, 'i_re,j_re,k_re=',i_re,j_re,k_re
             PRINT *, 'i_re2,j_re2,k_re2=',i_re2,j_re2,k_re2
          END IF
          new_coords(1) = i_re2
          new_coords(2) = j_re2
          new_coords(3) = k_re2
          IF ( matrix(i_re2,j_re2,k_re2)%sp_num .NE. 0 ) THEN
             IF ( ALL(coords .EQ. new_coords) .EQV. .FALSE.) THEN
                ! Determine if the hopped to species can react with the hopping species
                r1 = matrix(i_re,j_re,k_re)%sp_num
                r2 = matrix(i_re2,j_re2,k_re2)%sp_num
                CALL canreact(r1,r2,react_cube,null)
                ! Determine if site is occupied and if species can react
                IF ( null .EQ. 0 ) THEN
                   large_arr(n,4) = 1
                   large_count = large_count + 1
                ELSE
                   large_arr(n,4) = 0
                END IF

                large_arr(n,1)=i_re2
                large_arr(n,2)=j_re2
                large_arr(n,3)=k_re2
             END IF
          END IF
       END IF top_bottom
    END DO initial_loop

    ! Only look at normal sites if there are no interstitial reactants
    IF ( large_count .GT. 0 ) THEN
       CONTINUE
    ELSE
       DO n=1,4 ! Go to a phantom position
          IF ( n .EQ. 1 ) THEN
             IF ( (j_re-1 .GT. 0) .AND. (k_re-1 .GT. 0) .AND. (k_re+1 .LE. dimens(3))) THEN
                CALL hopping(i_re,j_re-1,k_re+1,i_re2,j_re2,k_re2,1,dimens)
             ELSE
                i_re2 = i_re
                j_re2 = 1
                k_re2 = 2
             END IF
          ELSE IF ( n .EQ. 2 ) THEN
             IF ( (j_re-1 .GT. 0) .AND. (k_re-1 .GT. 0) .AND. (k_re+1 .LE. dimens(3))) THEN
                CALL hopping(i_re,j_re-1,k_re-1,i_re2,j_re2,k_re2,2,dimens)
             ELSE
                i_re2 = i_re
                j_re2 = 1
                k_re2 = 2
             END IF
          ELSE IF ( n .EQ. 3 ) THEN
             IF ( (j_re+1 .LE. dimens(2)) .AND. (k_re-1 .GT. 0) .AND. (k_re+1 .LE. dimens(3))) THEN
                CALL hopping(i_re,j_re+1,k_re+1,i_re2,j_re2,k_re2,1,dimens)
             ELSE
                i_re2 = i_re
                j_re2 = 1
                k_re2 = 2
             END IF
          ELSE
             IF ( (j_re+1 .LE. dimens(2)) .AND. (k_re-1 .GT. 0) .AND. (k_re+1 .LE. dimens(3))) THEN
                CALL hopping(i_re,j_re+1,k_re-1,i_re2,j_re2,k_re2,2,dimens)
             ELSE
                i_re2 = i_re
                j_re2 = 1
                k_re2 = 2
             END IF
          END IF

          ! Determine if the matrix site is occupied and can react
          new_coords(1) = i_re2
          new_coords(2) = j_re2
          new_coords(3) = k_re2
          IF ( matrix(i_re2,j_re2,k_re2)%sp_num .NE. 0 ) THEN
             IF ( ALL(coords .EQ. new_coords) .EQV. .FALSE. ) THEN
                r1 = matrix(i_re,j_re,k_re)%sp_num
                r2 = matrix(i_re2,j_re2,k_re2)%sp_num
                !             PRINT *, 'r1=',r1,'and r2=',r2
                CALL canreact(r1,r2,react_cube,wait_list,null)
                IF ( null .EQ. 0 ) THEN
                   small_arr(n,4) = 1
                   small_count = small_count + 1
                ELSE
                   small_arr(n,4) = 0
                END IF
                small_arr(n,1)=i_re2
                small_arr(n,2)=j_re2
                small_arr(n,3)=k_re2
             ELSE
                small_arr(n,1) = i_re
                small_arr(n,2) = j_re
                small_arr(n,3) = k_re
                small_arr(n,4) = 0
             END IF
          END IF
       END DO
    END IF

    ! Determine if there has been a null event
    IF ( large_count .EQ. 0 .AND. small_count .EQ. 0 ) THEN
       null = 1
    ELSE
       null = 0
    END IF
  END SUBROUTINE lookaroundyou

  SUBROUTINE solarlottery( small_count,large_count, small_temp,large_temp, coords )
    !
    ! Purpose:
    !   This subroutine looks at either the interstitial or normal neighbors
    !  that contain a potential reacion partner and chooses one at random.
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! SOLARLOTTERY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IMPLICIT NONE

    !*****************
    ! Input and output
    !*****************
    INTEGER            , INTENT(IN)                                :: large_count   !number of interstitial reactants
    INTEGER            , INTENT(IN)                                :: small_count  !number of normal reactants
    INTEGER            , INTENT(OUT), DIMENSION(3)                 :: coords   !coordinates of selected site
    INTEGER                         , DIMENSION(6,4)               :: large_temp  !array of interstitial neighbor coords
    INTEGER                         , DIMENSION(4,4)               :: small_temp !array of normal neighbor coords

    !****************
    ! Local variables
    !****************
    INTEGER            , ALLOCATABLE, DIMENSION(:,:)               :: temp_arr !temporary array of coordinates
    INTEGER                                                        :: i,n   !counters
    INTEGER                                                        :: lucky    !index of selected coords
    REAL                                                           :: rand     !random number

    !    PRINT *, "Started solarlottery"
    !    PRINT *, "large_count is: ", large_count
    !    PRINT *, "small_count is: ", small_count

    i = 1
    IF ( large_count .NE. 0 ) THEN ! Check if any same type reactants exist
       ALLOCATE ( temp_arr(large_count,3) )
       DO n=1,SIZE(large_temp,1)
          IF ( large_temp(n,4) .NE. 0 ) THEN
             temp_arr(i,1) = large_temp(n,1)
             temp_arr(i,2) = large_temp(n,2)
             temp_arr(i,3) = large_temp(n,3)
             i = i + 1
          ELSE
             CONTINUE
          END IF
       END DO
    ELSE ! See if any opposite type reactants exist
       ALLOCATE ( temp_arr(small_count,3) )
       DO n=1,SIZE(small_temp,1)
          IF ( small_temp(n,4) .NE. 0 ) THEN
             temp_arr(i,1) = small_temp(n,1)
             temp_arr(i,2) = small_temp(n,2)
             temp_arr(i,3) = small_temp(n,3)
             i = i + 1
          ELSE
             CONTINUE
          END IF
       END DO
    END IF

    !    PRINT *, "temp_arr size is: ", SIZE(temp_arr)
    !    PRINT *, "The contents of temp_arr are: "
    !    DO n=1,SIZE(temp_arr,1)
    !      PRINT *, temp_arr(n,:)
    !    END DO

    CALL RANDOM_NUMBER(rand) ! Choose a random temp_arr element
    coords = 0
    IF ( SIZE(temp_arr,1) .EQ. 1 ) THEN
       coords = temp_arr(1,:)
       !      DO n=1,3
       !        coords(n) = temp_arr(1,n)
       !      END DO
    ELSE
       lucky = INT(rand*SIZE(temp_arr,1)) + 1
       DO n=1,3
          coords(n) = temp_arr(lucky,n)
       END DO
    END IF
    !    PRINT *, "The coords are: ",coords, "ending Solarlottery"

  END SUBROUTINE solarlottery

  SUBROUTINE thirdman( prod,i_re,j_re,k_re, null, prod_coords)
    !
    ! Purpose:
    !   This subroutine attempts to find an empty site to place a third product
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! THIRDMAN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    !******************!
    ! Input and output !
    !******************!

    INTEGER            , INTENT(IN)                                       :: prod !product to be placed
    INTEGER(KIND=SHORT), INTENT(OUT)                                      :: null !error flag
    INTEGER            , INTENT(IN)                                       :: i_re,j_re,k_re !coords of original site
    INTEGER            , INTENT(OUT), DIMENSION(3)             , OPTIONAL :: prod_coords !product placement coords
    REAL                            , DIMENSION(:)    , POINTER           :: en_list

    !*****************!
    ! Local variables !
    !*****************!

    INTEGER                                                               :: n,m !counters
    INTEGER                                                               :: i_pr,j_pr,k_pr !product coords

    DO n=1,6
       IF ( i_re .EQ. 1 .AND. ( n .EQ. 5 .OR. n .EQ. 6 ) ) THEN
          CONTINUE
       ELSE IF ( i_re .EQ. dimens(1) .AND. n .EQ. 6 ) THEN
          ! Don't hop down if on bottom layer
          CONTINUE
       ELSE
          CALL hopping(i_re,j_re,k_re,i_pr,j_pr,k_pr,n,dimens)
          IF ( matrix(i_pr,j_pr,k_pr)%sp_num .EQ. 0 ) THEN
             GOTO 1985
          END IF
       END IF
    END DO

    DO m=1,4 ! Go to a phantom position
       SELECT CASE (m)
       CASE (1)
          IF ((j_re-1 .GT. 0) .AND. (k_re-1 .GT. 0) .AND. (k_re+1 .LE. dimens(3))) THEN
             CALL hopping(i_re,j_re-1,k_re+1,i_pr,j_pr,k_pr,1,dimens)
             IF (matrix(i_pr,j_pr,k_pr)%sp_num .EQ. 0 ) THEN
                GOTO 1985
             ELSE
                CONTINUE
             END IF
          ELSE
             CONTINUE
          END IF
       CASE (2)
          IF (((j_re-1 .GT. 0) .AND. (k_re-1 .GT. 0) .AND. (k_re+1 .LE. dimens(3)))) THEN
             CALL hopping(i_re,j_re-1,k_re-1,i_pr,j_pr,k_pr,2,dimens)
             IF (matrix(i_pr,j_pr,k_pr)%sp_num .EQ. 0 ) THEN
                GOTO 1985
             ELSE
                CONTINUE
             END IF
          ELSE
             CONTINUE
          END IF
       CASE (3)
          IF ((j_re+1 .LE. dimens(2)) .AND. (k_re-1 .GT. 0) .AND. (k_re+1 .LE. dimens(3))) THEN
             CALL hopping(i_re,j_re+1,k_re+1,i_pr,j_pr,k_pr,1,dimens)
             IF (matrix(i_pr,j_pr,k_pr)%sp_num .EQ. 0 ) THEN
                GOTO 1985
             ELSE
                CONTINUE
             END IF
          ELSE
             CONTINUE
          END IF
       CASE (4)
          IF ((j_re+1 .LE. dimens(2)) .AND. (k_re-1 .GT. 0) .AND. (k_re+1 .LE. dimens(3))) THEN
             CALL hopping(i_re,j_re+1,k_re-1,i_pr,j_pr,k_pr,2,dimens)
             IF (matrix(i_pr,j_pr,k_pr)%sp_num .EQ. 0 ) THEN
                GOTO 1985
             ELSE
                GOTO 2001
             END IF
          ELSE
             GOTO 2001
          END IF
       END SELECT

2001   IF ( m .EQ. 4 ) THEN
          null = 1
          !          PRINT *, 'No reaction possible! ERROR!!!'
          RETURN
       END IF
    END DO

    ! Place reactant at chosen site
1985 matrix(i_pr,j_pr,k_pr)%sp_num = prod
    CALL wait_calc(i_pr,j_pr,k_pr)
    temp => matrix(i_pr,j_pr,k_pr)
    CALL add_node(root,temp)

    ! Assign output coordinates array
    IF ( PRESENT(prod_coords) ) THEN
       prod_coords(1) = i_pr
       prod_coords(2) = j_pr
       prod_coords(3) = k_pr
    END IF
  END SUBROUTINE thirdman

  SUBROUTINE fallout ( o3_prod,o3_dest)
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! FALLOUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE
    !******************!
    ! Input and output !
    !******************!
    INTEGER                              , POINTER :: o3_prod,o3_dest

    !*****************!
    ! Local variables !
    !*****************!
    INTEGER(KIND=SHORT)                            :: null
    INTEGER                                        :: n,nn,jj
    INTEGER                                        :: num_elecs ! number of secondary electrons pruduced
    INTEGER                                        :: num_izns
    INTEGER                                        :: num_exs, num_els
    INTEGER                                        :: x,y,z !coordinates of cosmic-ray along track
    INTEGER                                        :: exitcount
    INTEGER                                        :: count_count
    INTEGER                                        :: step,estep !distance the track is incremented
    INTEGER                                        :: switch,eswitch
    INTEGER            , DIMENSION(3)              :: ev_coords !coordiantes of collision
    INTEGER            , DIMENSION(3)              :: elec_coords
    INTEGER            , DIMENSION(3)              :: prev, curr, next
    INTEGER                                        :: enull1,enull2
    REAL(KIND=DBL)                                 :: erand, emfp,de
    REAL(KIND=DBL)                                 :: p,u,rand1 ! rand num
    REAL(KIND=DBL)                                 :: ion_dist ! distance from last ionization
    REAL(KIND=DBL)                                 :: sigma_tot !total cross-section
    REAL(KIND=DBL)                                 :: mfp ! mean free path
    REAL(KIND=DBL)                                 :: dz ! move dist
    REAL(KIND=DBL)                                 :: dist_trav !distance travelled since last collision
    REAL(KIND=DBL)                                 :: e_loss,e_ion,e_exc,ee_loss
    REAL(KIND=DBL)                                 :: labtheta
    REAL(KIND=DBL)                       , TARGET  :: e_se
    REAL(KIND=DBL)                                 :: subexrand
    DOUBLE PRECISION                     , TARGET  :: energy_target
    DOUBLE PRECISION                     , POINTER :: ione
    DOUBLE PRECISION   , DIMENSION(:)    , POINTER :: psigij,psigexj
    CHARACTER(len=15)                              :: nature
    TYPE(SIGMA_BOX)    , DIMENSION(:)    , POINTER :: psigmas
    TYPE(SE_INFO)                                  :: se_box
    !*************************************************************************
    !Proton cross-section data, to be phased out and replaced with a struct as
    !with se_box
    !*************************************************************************
    DOUBLE PRECISION   , DIMENSION(:), ALLOCATABLE, TARGET :: psigij_target,psigexj_target
    TYPE(SIGMA_BOX)    , DIMENSION(:), ALLOCATABLE, TARGET :: psigmas_target
    INTEGER :: thinghit
    LOGICAL :: proceed


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
       x = 1 + FLOOR( dimens(3)*p )
       y = 1 + FLOOR( dimens(2)*u )
       z = 1
       IF ( TRACKPLOT .EQV. .TRUE. ) THEN
          PRINT *, "Trackplot on"
          IF ( x .GT. ((dimens(3)/2)-(dimens(3)*0.1)) .AND. x .LT. ((dimens(3)/2)+(dimens(3)*0.1)) &
               .AND. y .GT. ((dimens(2)/2)-(dimens(2)*0.1)) .AND. y .LT. ((dimens(2)/2)+(dimens(2)*0.1)) ) proceed = .TRUE.
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
       IF ( DEBUG .EQV. .TRUE. ) PRINT *, 'Now entering loop: z=',z,' and dimens(1)=',dimens(1),' and step=',step
       count_count = count_count + 1
       !      IF ( MOD(count_count,1000) .EQ. 0 ) CALL counter(time, AB_UNIT_NUM, matrix, wait_list, 4,7)
       ! Define event coords
       ev_coords(1) = z+step
       ev_coords(2) = y
       ev_coords(3) = x
       !      PRINT *, 'The event coords are:',ev_coords

       thinghit = matrix(ev_coords(1),ev_coords(2),ev_coords(3))%sp_num


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
                  se_box%se_energy = e_se
                  !**********************
                  ! For debugging \/ \/
                  !**********************
                  IF ( SECELEC .EQV. .FALSE. ) se_box%se_energy = 0D0
                  !            ELSE IF ( rand1 .LE. DISPROB .OR. ANY(FRAGILE .EQ. -1*thinghit) ) THEN
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
            PROTON_ELOSS = e_loss
          END ASSOCIATE
          ione = ione - e_loss
          CALL psigma_suite(ione,psigmas,psigij,psigexj)
          IF ( TRACKPLOT .EQV. .TRUE. ) THEN
             count_count = count_count + 1
             WRITE(TRACKPLOT_UNIT_NUM,*) ev_coords(1),',',ev_coords(2),',',ev_coords(3),', proton,',nature
          END IF

          !****************************************************************************!
          ! Elastic collision                                                          !
          !****************************************************************************!
          IF ( switch .EQ. 0 ) THEN
             !          PRINT *, "Elastic collision"
             !NB: FUTURE WORK TO ADD LATTICE DAMAGE
             CONTINUE
             !****************************************************************************!
             ! Excitation                                                                 !
             !****************************************************************************!
          ELSE IF ( switch .EQ. 1 ) THEN ! Dissociate target species on track and place prods
             CALL cern( o3_prod,o3_dest,null, ev_coords, switch)
             IF ( null .EQ. 1 ) RETURN
             !****************************************************************************!
             ! Ionization                                                                 !
             !****************************************************************************!
          ELSE IF ( switch .EQ. 2 .AND. z+step .NE. 1 .AND. z+step .NE. 2 ) THEN
             IF ( se_box%se_energy .LE. ECUTOFF ) THEN
                CALL base_ionization( o3_prod,o3_dest,ev_coords, null )
                IF ( null .EQ. 1 ) RETURN
             ELSE
                !
                ! Generate secondary electrons/electron track
                !
                !*******************************************************************
                ! Call base_ionization to generate the first-generation secondary electron
                ! NB: the electron should be the second product in the "prods" array
                !*******************************************************************
                null = 0
                CALL base_ionization( o3_prod,o3_dest,ev_coords, null,elec_coords )

                IF ( TRACKPLOT .EQV. .TRUE. ) THEN
                   count_count = count_count + 1
                   WRITE(TRACKPLOT_UNIT_NUM,*) elec_coords(1),',',elec_coords(2),',',elec_coords(3),", electron, Ionization"
                END IF

                !GOTO jumps down to calling next random number
                IF ( null .EQ. 1 ) GOTO 100
                !*******************************************************************
                ! Electron Impact Processes
                !*******************************************************************
                ! Note:
                !  The processes in the following loop correspond to conventional
                ! processes such as electron-impact excitation and ionization. These
                ! types of collisional events are treated semi-classically.
                !*******************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! ELECTRON_IMPACT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                curr     = ev_coords
                next     = elec_coords
                ion_dist = 0
                enull1   = 0
                enull2   = 0
                emfp     = 0
                p        = 0
                de       = 0
                estep    = 0
                eswitch  = 0
                ee_loss  = 0
                !Initialize se_box
                IF ( DEBUG .EQV. .TRUE. ) PRINT *, "Now initializing the se_box"
                CALL se_info_init(se_box)
                IF ( DEBUG .EQV. .TRUE. ) PRINT *, 'se energy is:',se_box%se_energy,' and ECUTOFF is',ECUTOFF
                DO WHILE ( se_box%se_energy .GE. ECUTOFF )
                   IF ( DEBUG .EQV. .TRUE. ) PRINT *, 'Electron box initialized, calling loop'
                   !If the electron no longer has sufficient
                   !energy, exit the loop.
                   IF ( enull1 .NE. 0 .AND. enull2 .NE. 0 ) EXIT

                   !Calculate electron cross_sections
                   CALL esigma_suite(se_box)
                   IF ( DEBUG .EQV. .TRUE. ) PRINT *, "se_box%ineltot=", se_box%se_ineltot
                   IF ( DEBUG .EQV. .TRUE. ) PRINT *, "se_box%se_energy=", se_box%se_energy

                   !Calculate hopping distance
                   emfp  = 1./(RHO*(se_box%se_ineltot+1.0e-17))
                   CALL RANDOM_NUMBER(p)
                   de    = -1.*emfp*LOG(1.-p)
                   de    = de*ESTEPFAC
                   estep = INT(de/C_PR)

                   !Have a minumum hopping distance of 1
                   IF ( estep .EQ. 0 ) estep = 1

                   DO n=1,estep
                      !Each transport hop is like one step
                      prev = curr
                      curr = next
                      CALL transport(prev,curr,next,matrix)

                      IF ( matrix(curr(1),curr(2),curr(3))%sp_num .EQ. O3NUM ) THEN
                         CALL RANDOM_NUMBER(rand1)
                         IF ( rand1 .LE. O3_DIS_BRANCHING ) THEN
                            CALL cern( o3_prod,o3_dest,null, curr, 1 )
                            WRITE(TRACKPLOT_UNIT_NUM,*) curr(1),',',curr(2),',',curr(3),", electron , Excitation"
                         END IF
                      END IF

                      IF ( TRACKPLOT .EQV. .TRUE. ) THEN
                         count_count = count_count + 1
                         WRITE(TRACKPLOT_UNIT_NUM,*) curr(1),',',curr(2),',',curr(3),", electron , Movement"
                      END IF
                      ee_loss = se_box%se_energy*ELASTIC_LOSS
                      se_box%se_energy =  se_box%se_energy - ee_loss
                   END DO

                   !Determine nature of event
                   CALL RANDOM_NUMBER(erand)
                   IF ( erand .GT. 0 .AND. erand .LE. (se_box%se_iontot/se_box%se_ineltot) ) THEN
                      eswitch = 1
                   ELSE
                      eswitch = 0
                   END IF

                   !Carry out impact collision
                   IF ( (matrix(next(1),next(2),next(3))%sp_num .NE. 0) .AND. (eswitch .EQ. 1) ) THEN
                      !Electron impact ionization
                      IF ( DEBUG .EQV. .TRUE. ) PRINT *, 'EII: SE is hopping to site with',matrix(next(1),next(2),next(3))%sp_num
                      CALL base_ionization(o3_prod,o3_dest,next, null )
                      IF ( null .EQ. 1 ) GOTO 100
                      CALL e_ion_select(se_box,e_ion,enull1)
                      ee_loss = e_ion
                      WRITE(TRACKPLOT_UNIT_NUM,*) next(1),',',next(2),',',next(3),", electron , Ionization"
                   ELSE IF ( (matrix(next(1),next(2),next(3))%sp_num .NE. 0) .AND. (eswitch .EQ. 0) ) THEN
                      !Electron impact excitation
                      IF ( DEBUG .EQV. .TRUE. ) PRINT *, 'EIE: SE is hopping to site with',matrix(next(1),next(2),next(3))
                      CALL RANDOM_NUMBER(rand1)

                      IF ( (rand1 .LE. DISPROB) .AND. (matrix(curr(1),curr(2),curr(3))%sp_num .NE. 0) ) THEN
                         CALL cern( o3_prod,o3_dest,null, next, 1)
                      ELSE IF ( ANY( FRAGILE .EQ. matrix(curr(1),curr(2),curr(3))%sp_num)) THEN
                         ! Test for fragile species
                         CALL cern( o3_prod,o3_dest,null, next, 1)
                      END IF

                      CALL e_ex_select(se_box,e_exc,enull2)
                      ee_loss = e_exc
                      WRITE(TRACKPLOT_UNIT_NUM,*) next(1),',',next(2),',',next(3),", electron , Excitation"
                   END IF

                   !Update the secondary electron energy
                   !              PRINT *, 'E_se:',ese_point,' E_loss:',ee_loss
                   se_box%se_energy = se_box%se_energy - ee_loss
                END DO

                !*******************************************************************
                ! Sub-Excitation Processes
                !*******************************************************************
                ! Note:
                !  The processes in the following loop correspond to low-energy, or
                ! sub-excitation processes, in which the electron has lost enough energy
                ! to be unable to excite the target efficiently.
                !*******************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! SUB_EXCITATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                IF ( NSUBEX .NE. 0 ) THEN
                   nn = 0
                   exitcount = 0

                   DO WHILE ( nn .LT. NSUBEX .AND. exitcount .LT. NEXIT) !NSUBEX is the number of sub-excitation collisions
                      IF ( DEBUG .EQV. .TRUE. ) THEN
                         PRINT *, 'nn=',nn
                         PRINT *, 'exitcount=',exitcount
                      END IF
                      exitcount = exitcount + 1
                      prev = curr
                      curr = next
                      CALL transport(prev,curr,next)

                      IF ( TRACKPLOT .EQV. .TRUE. ) THEN
                         count_count = count_count + 1
                         WRITE(TRACKPLOT_UNIT_NUM,*) curr(1),',',curr(2),',',curr(3), ", subex electron, Movement"
                      END IF

                      IF ( matrix(next(1),next(2),next(3))%sp_num .NE. 0 ) THEN
                         IF ( DEBUG .EQV. .TRUE. ) PRINT *, 'matrix in subexloop=',matrix(next(1),next(2),next(3))%sp_num
                         !Carry out dissociate electron attachment
                         !NB: In the model, this is functionally identical to
                         !an ordinary ionization
                         CALL base_ionization(o3_prod,o3_dest,next, null )
                         nn = nn + 1
                         WRITE(TRACKPLOT_UNIT_NUM,*) next(1),',',next(2),',',next(3), ", subex electron, Ionization"
                      END IF
                   END DO
                END IF
                IF ( DEBUG .EQV. .TRUE. ) PRINT *, 'Finishing sub_excitation processes'

                !Manual garbage collection
                CALL se_info_garbage(se_box)
             END IF
          END IF
       ELSE
          IF (ev_coords(1) .EQ. SIZE(matrix,1)/2) THEN
             RETURN
          ELSE
             CONTINUE
          END IF
       END IF

       !********************!
       ! Track Plotting bit !
       !********************!
       dist_trav = dist_trav + dz

       IF ( DEBUG .EQV. .TRUE. ) PRINT *, "z=",z,"of",dimens(1)," which is",(REAL(z)/REAL(dimens(1)))*100,"% of thickness"

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
       step = INT((dz/(THICK/REAL(dimens(1)))))

       ! Make sure the next site is different than the previous one
       IF ( z+step .EQ. z ) THEN
          z = z + 1
       END IF

       ! Make sure the site is greater than the previous one
       IF (step .LE. 0. ) GOTO 100
       IF (z+step .GE. dimens(1) ) THEN
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

  SUBROUTINE krell( in_coords,out_coords,null )
    ! Purose:
    !   This subtroutine takes some ion/bulk interaction site and finds a nearby
    !  empty site to put a second product. The return of the function is a set of
    !  coordinates.
    !
    ! Note:
    !   null = 1 => no empty sites
    !   null = 0 => empty site available
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! KRELL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE
    !*****************
    ! Input and output
    !*****************
    INTEGER            , INTENT(IN) , DIMENSION(3)               :: in_coords !coords of reaction site
    INTEGER            , INTENT(OUT), DIMENSION(3)               :: out_coords !coords for product
    INTEGER(KIND=SHORT)                                          :: null !null error flag
    !****************
    ! Local variables
    !****************
    INTEGER                                                      :: large_count
    INTEGER                                                      :: small_count
    INTEGER                                                      :: lucky !index of selected site, from rand
    INTEGER                                                      :: i,n !counters
    INTEGER                                                      :: i_re, j_re, k_re !in coords
    INTEGER                                                      :: i_re2,j_re2,k_re2 !out coords
    INTEGER                         , DIMENSION(6,3)             :: large_temp
    INTEGER                         , DIMENSION(4,3)             :: small_temp
    INTEGER            , ALLOCATABLE, DIMENSION(:,:)             :: temp_arr !temporary empty site array
    REAL                                                         :: rand !random number

    ! Initialize counters and arrays
    large_count   = 0
    small_count  = 0

    ! Assign i_re, j_re, k_re to in_coords
    i_re = in_coords(1)
    j_re = in_coords(2)
    k_re = in_coords(3)

    large_temp = 0
    small_temp = 0
    DO n=1,6
       IF ( i_re .EQ. 1 .AND. ( n .EQ. 5 .OR. n .EQ. 6 ) ) THEN
          ! If on top layer, stay on top layer
          CONTINUE
       ELSE IF ( i_re .EQ. dimens(1) .AND. n .EQ. 6 ) THEN
          ! Don't hop down if on bottom layer
          CONTINUE
       ELSE
          CALL hopping(i_re,j_re,k_re,i_re2,j_re2,k_re2,n )
          IF ( matrix(i_re2,j_re2,k_re2)%sp_num .EQ. 0 ) THEN
             large_temp(n,1)=i_re2
             large_temp(n,2)=j_re2
             large_temp(n,3)=k_re2
             large_count = large_count + 1
          END IF
       END IF
    END DO

    IF ( large_count .GT. 0 ) THEN
       CONTINUE
    ELSE

       DO n=1,4
          ! Go to a phantom position to hop to nearest neighbors
          SELECT CASE (n)
          CASE (1)
             IF ( (j_re-1 .GT. 0) .AND. (k_re-1 .GT. 0) .AND. (k_re+1 .LE. dimens(3))) THEN
                CALL hopping(i_re,j_re-1,k_re+1,i_re2,j_re2,k_re2,1)
             ELSE
                GOTO 1944
             END IF

          CASE (2)
             IF ( (j_re-1 .GT. 0) .AND. (k_re-1 .GT. 0) .AND. (k_re+1 .LE. dimens(3))) THEN
                CALL hopping(i_re,j_re-1,k_re-1,i_re2,j_re2,k_re2,2)
             ELSE
                GOTO 1944
             END IF

          CASE (3)
             IF ( (j_re+1 .LE. dimens(2)) .AND. (k_re-1 .GT. 0) .AND. (k_re+1 .LE. dimens(3))) THEN
                CALL hopping(i_re,j_re+1,k_re-1,i_re2,j_re2,k_re2,2)
             ELSE
                GOTO 1944
             END IF

          CASE (4)
             IF ( (j_re+1 .LE. dimens(2)) .AND. (k_re-1 .GT. 0) .AND. (k_re+1 .LE. dimens(3))) THEN
                CALL hopping(i_re,j_re+1,k_re+1,i_re2,j_re2,k_re2,1)
             ELSE
                GOTO 1944
             END IF
          END SELECT

          IF ( matrix(i_re2,j_re2,k_re2)%sp_num .EQ. 0 ) THEN
             small_temp(n,1)=i_re2
             small_temp(n,2)=j_re2
             small_temp(n,3)=k_re2
             small_count = small_count + 1
          END IF
1944      CONTINUE
       END DO
    END IF

    ! Determine if there has been a null event
    IF ( large_count .EQ. 0 .AND. small_count .EQ. 0 ) THEN
       null = 1
       !      PRINT *, 'ERROR: No reaction in Krell possible'
       out_coords = 314159
       RETURN
    ELSE
       null = 0
    END IF

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
  END SUBROUTINE krell

  SUBROUTINE action_figure ( temp )
    !
    ! Purpose:
    !  The purpose of this subroutine is to take the
    !  first element of the waiting list, which will
    !  have been determined using the "roll_call"
    !  subroutine, decide which action should be
    !  performed, i.e. desorption or diffusion. In
    !  the case of bulk species, i.e. i .NE. 1,
    !  there is only the possibility of bulk diffusion.
    !
    ! Note:
    !  The subroutine returns an integer value,
    !  called the "flag" that is used by the code
    !  to execute the appropriate action, e.g.
    !  thermal hopping. The values of the flag are:
    !
    !  -- act_type = 1 => thermal hopping
    !  -- act_type = 2 => desorption
    !  -- act_type = 3 => fast reaction
    !
    ! Note:
    !  This subroutine is only called in the case of
    !  the normal motions of mobile species. It is
    !  NOT called for cosmic-ray or photon events,
    !  which are treated separately.
    !
    ! Documentation:
    !  DATE          PROGRAMMER           DESCRIPTION
    !  ========      ==========           ===========
    !  20150415      C. Shingledecker     Original code
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! ACTION_FIGURE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    ! Data dictionary
    TYPE(node), POINTER                                 :: temp
    REAL(KIND=DBL)                                      :: b_1 !thermal surface hopping rate
    REAL(KIND=DBL)                                      :: b_2 !surface desorption rate
    REAL(KIND=DBL)                                      :: comp_val !to determine which action occurs
    REAL(KIND=DBL), INTENT(IN)                          :: rand_num
    REAL                      , DIMENSION(:)  , POINTER :: en_list
    TYPE (wait_info)          , DIMENSION(:)  , POINTER :: wait_list

    ! (1) Decide whether or not the species is on the surface
    IF ( temp%coord1 .EQ. 1 ) THEN
       ! (1a) Species is on the surface
       ! Calculate b-rates to compare
       b_1 = trl_nu * EXP( -1*(  en_list(wait_list(index)%sp_num)*E_SURF / kin_temp  ) )
       b_2 = trl_nu * EXP( -1*(  en_list(wait_list(index)%sp_num)        / kin_temp  ) )
       comp_val = b_1 / (b_1 + b_2)
       ! Decide whether desorption or hopping occurs
       IF ( rand_num .LT. comp_val ) THEN
          ! Diffusion occurs
          temp%act_type = 1
       ELSE
          ! Desorption occurs
          temp%act_type = 2
       END IF
    ELSE
       ! (1b) Species is in the bulk
       ! Only hopping (diffusion) can occur
       temp%act_type = 1
    END IF

    ! If the species is in fast reacting
    IF( ANY(FAST_REACTS .EQ. temp%sp_num) ) THEN
       temp%act_type = 3
    END IF
  END SUBROUTINE action_figure

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
    ! Warning!: As of original code, lateral bonds are not
    !  considered, as described in CH14.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! WAIT_CALC !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    ! Data dictionary: variables passed to the subroutine
    TYPE(node), POINTER :: temp
    REAL(KIND=DBL)                                       :: rand_num  !pseudorandom number
    REAL(KIND=DBL)                                       :: b_1       !surface thermal hopping rate
    REAL(KIND=DBL)                                       :: b_2       !surface desorption rate
    REAL(KIND=DBL)                                       :: b_3       !bulk diffusion rate
    REAL(KIND=DBL)                                       :: b         !total rate, from CH14

    IF ( wait_list(index)%i .EQ. 1 ) THEN
       ! Surface species, separate rates for
       ! desorption and diffusion
       b_1 = trl_nu*EXP( -1*( ( en_list(temp%sp_num)*E_SURF) / kin_temp ) )
       b_2 = trl_nu*EXP( -1*( en_list(temp%sp_num)           / kin_temp ) )
       b = b_1 + b_2
    ELSE
       ! Bulk species, only bulk diffusion
       b_3 = trl_nu*EXP( -1*( en_list(temp%sp_num)*E_BULK    / kin_temp ) )
       b = b_3
    END IF
    CALL RANDOM_NUMBER(rand_num)
    ! Calculate waiting time
    temp%wait_time = (-1*LOG(rand_num) / b) + time

    ! Assign action type for next move
    CALL action_figure(temp)
  END SUBROUTINE wait_calc

  SUBROUTINE counter(o3_prod,o3_dest,numprotons,sp1,sp2)
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! COUNTER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    INTEGER                                          , POINTER :: o3_prod,o3_dest
    INTEGER(KIND=LONG) , INTENT(IN)                            :: numprotons
    INTEGER                                                    :: i,j,k,nn
    INTEGER                                                    :: o_count,o2_count,o3_count,sp3
    INTEGER            , INTENT(IN)                            :: sp1, sp2
    REAL(KIND=DBL)                                             :: volume
    REAL(KIND=DBL)                                             :: denom
    REAL(KIND=DBL)                                             :: area
    REAL(KIND=DBL)                                             :: fluence
    CHARACTER(len=80)                                          :: varfmt

    !    volume = THICK*EDGE*EDGE
    area   = EDGE*EDGE
    denom = THICK*EDGE*EDGE*1.0E20
    sp3 = 1
    o_count = 0
    o2_count = 0
    o3_count = 0
    wrong_count = 0
    ! Method 1 of fluence calculation
    fluence  = CR_FLUX*time ! Note: This is the x-value for the objective function
    ! Method 2 of fluence calculation (only use 1 at a time )
    ! fluence = numprotons/area


    IF ( TEST_WRONG .EQV. .TRUE. ) OPEN(UNIT=1013,FILE="counter_test_wrong_spaces.txt",STATUS='REPLACE')
    DO k = 1,dimens(3)
       DO j = 1,dimens(2)
          DO i = 1,dimens(1)
             IF ( matrix(i,j,k)%sp_num .EQ. O3NUM ) THEN
                o3_count = o3_count + 1
             ELSE IF ( matrix(i,j,k)%sp_num .EQ. O2NUM ) THEN
                o2_count = o2_count + 1
             ELSE IF ( matrix(i,j,k)%sp_num .EQ. ONUM ) THEN
                o_count = o_count + 1
             END IF
          END DO
       END DO
    END DO

    O_ABUNDANCE  = o_count
    O2_ABUNDANCE = o2_count
    O3_ABUNDANCE = o3_count

    IF ( NO_OUTPUT .EQV. .FALSE. ) THEN
       WRITE(AB_UNIT_NUM,*) ALTFLUENCE,',', fluence,',',time,',',o2_count,',',o_count,',',o3_count,',',numprotons
    END IF

    IF ( QUIET .EQV. .FALSE. ) THEN
       varfmt = "(A6,ES10.4,A9,ES10.4)"
       PRINT varfmt, " TIME=",time,"FLUENCE=",fluence
       varfmt = "(A5,ES10.4,A6,ES10.4)"
       PRINT varfmt, " [O]=",o_count/denom," [O3]=",o3_count/denom
       varfmt = "(A16,F10.4,A16,I10)"
       PRINT *, '[O3] PROD/DEST =', (REAL(o3_prod)/REAL(o3_dest))," WAIT LENGTH=",wait_len
       PRINT *, 'O3_prod=',o3_prod, 'O3_dest=',o3_dest
       PRINT *, 'denom=',denom
       PRINT *, '***********************************************************************'
    END IF
  END SUBROUTINE counter

  SUBROUTINE transport(prev,curr,next,matrix)
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! TRANSPORT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    INTEGER            , INTENT(IN)          , DIMENSION(3)     :: prev !coordinates of previous location
    INTEGER            , INTENT(IN)          , DIMENSION(3)     :: curr !current coordinates
    INTEGER            , INTENT(OUT)         , DIMENSION(3)     :: next !coordinates of next position
    INTEGER                         , POINTER, DIMENSION(:,:,:) :: matrix !the solid matrix

    !****************
    ! Local variables
    !****************

    INTEGER                                                      :: n !counters
    INTEGER                                  , DIMENSION(3)      :: dimens !dimensions of ice matrix
    REAL                                                         :: insides
    REAL                                                         :: hopdist,nextdist
    REAL                                                         :: rand,rand2 !random number
    REAL                                                         :: sigma
    LOGICAL                                                      :: carnap

    ! Obtain the dimensions of the matrix
    dimens(1) = SIZE(matrix,1)
    dimens(2) = SIZE(matrix,2)
    dimens(3) = SIZE(matrix,3)

    !    PRINT *,'The coordinates of the event are',in_coords

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
          ELSE IF ( curr(1) .EQ. dimens(1) .AND. n .EQ. 6 ) THEN
             ! Don't hop down if on bottom layer
             CONTINUE
          ELSE
             CALL hopping(curr(1),curr(2),curr(3),next(1),next(2),next(3),n,dimens )
          END IF
       ELSE
          n = 1 + FLOOR(4*rand2)
          ! Go to a phantom position to hop to nearest neighbors
          SELECT CASE (n)
          CASE (1)
             IF ( curr(2)-1 .GT. 0          .AND. &
                  curr(3)-1 .GT. 0          .AND. &
                  curr(3)+1 .LE. dimens(3) ) THEN
                CALL hopping(curr(1),curr(2)-1,curr(3)+1,next(1),next(2),next(3),1,dimens)
             ELSE
                CONTINUE
             END IF
          CASE (2)
             IF ( curr(2)-1 .GT. 0          .AND. &
                  curr(3)-1 .GT. 0          .AND. &
                  curr(3)+1 .LE. dimens(3) ) THEN
                CALL hopping(curr(1),curr(2)-1,curr(3)-1,next(1),next(2),next(3),2,dimens)
             ELSE
                CONTINUE
             END IF
          CASE (3)
             IF ( curr(2)+1 .LE. dimens(2)  .AND. &
                  curr(3)-1 .GT. 0          .AND. &
                  curr(3)+1 .LE. dimens(3) ) THEN
                CALL hopping(curr(1),curr(2)+1,curr(3)-1,next(1),next(2),next(3),2,dimens)
             ELSE
                CONTINUE
             END IF
          CASE (4)
             IF ( curr(2)+1 .LE. dimens(2)  .AND. &
                  curr(3)-1 .GT. 0          .AND. &
                  curr(3)+1 .LE. dimens(3) ) THEN
                CALL hopping(curr(1),curr(2)+1,curr(3)+1,next(1),next(2),next(3),1,dimens)
             ELSE
                CONTINUE
             END IF
          END SELECT
       END IF
       IF ( n .NE. 5 .AND. n .NE. 6 ) THEN
          !  IF ( (curr(1) .NE. next(1) ) .AND. (curr(2) .NE. next(2) ) .AND. (curr(3) .NE. next(3))) THEN
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! PSIGMA_SUITE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! ESIGMA_SUITE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    ! PRINT *, "Sum of allowed sigs is", se_box%se_alwd_extot

    ! !(4) Calculate forbidden excitation cross-sections
    ! !NB: As above, no loop is required, since the subroutine
    ! !    returns an array of values
    ! DO n=1,SIZE(o2_e_ex_fbdn)
    !   ASSOCIATE( e => se_box%se_energy            , &
    !     w => se_box%se_fbdn(n)%wj_fbdn   , &
    !     f => se_box%se_fbdn(n)%fj_fbdn   , &
    !     o => se_box%se_fbdn(n)%omega_fbdn, &
    !     a => se_box%se_fbdn(n)%alpha_fbdn, &
    !     b => se_box%se_fbdn(n)%beta_fbdn    )
    !     se_box%se_fbdnsigs(n) = 0.0 !greendutta(e,f,w,o,a,b)
    !   END ASSOCIATE
    ! END DO
    ! se_box%se_fbdn_extot = SUM(se_box%se_fbdnsigs)

    !(5) The total electron impact excitation is the sum of the
    !    allowed and forbidden transition cross-sections
    se_box%se_extot = se_box%se_alwd_extot !+ se_box%se_fbdn_extot

    !(6) Calculate the total cross-section as the sum of the
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! P_ION_SELECT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    !    rn       = RAND()
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! P_EX_SELECT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    !    rn       = RAND()
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

  SUBROUTINE e_ion_select(se_box,e_loss,null)
    !
    ! Purpose:
    !   This subroutine is to determine the specific ionization state that an
    !  inelastic collision ionizes from.
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! E_ION_SELECT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    !Data dictionary: Calling parameters
    TYPE(se_info)                                  :: se_box
    DOUBLE PRECISION   , INTENT(OUT)                     :: e_loss
    INTEGER            , INTENT(OUT)                     :: null

    !Data dictionary: Local variables
    DOUBLE PRECISION                                     :: e_ion,e_se,e_exc
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
    e_exc = 0D0
    e_loss = 0D0

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
    !    rn       = RAND()
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

    !Draw another pseudo-random number, this time from a Gamma distribution
    !to determine the kinetic energy of the low-energy electron.
    !
    !NB: The input to rgamma, aval, is a global parameter
    e_se = rgamma(AVAL)
    e_loss = e_ion + e_se
    RETURN
  END SUBROUTINE e_ion_select

  SUBROUTINE e_ex_select(se_box,e_exc,null)
    !
    ! Purpose:
    !   This subroutine is to determine the specific excited state that an
    !  inelastic collision results in the target species being promoted to.
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! E_EX_SELECT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

    !Determine which type of transition will occur
    !    rn       = RAND()
    ! CALL RANDOM_NUMBER(rn)
    ! prob = (se_box%se_alwd_extot/se_box%se_ineltot)
    ! IF ( rn .GT. prob ) THEN
    !   ALLOCATE(arr(SIZE(se_box%se_fbdnsigs),2))
    !   arr = 0
    !   arr(:,1) = se_box%se_fbdnsigs
    !   arr(:,2) = se_box%se_fbdn%wj_fbdn
    !   sigtot   = se_box%se_fbdn_extot
    ! ELSE
    ALLOCATE(arr(SIZE(se_box%se_alwdsigs),2))
    arr = 0
    arr(:,1) = se_box%se_alwdsigs
    arr(:,2) = se_box%se_alwd%wj_alwd
    sigtot   = se_box%se_alwd_extot
    ! END IF

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
    !    rn       = RAND()
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! ELASTIC_EVENT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    !    rn = RAND()
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
    !  an energy assigned at the time of calling.
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! SE_INFO_INIT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    TYPE(se_info) :: se_box

    !(1) Initialize scalar values
    se_box%se_iontot     = 0
    se_box%se_extot      = 0
    se_box%se_ineltot    = 0
    se_box%se_alwd_extot = 0
    ! se_box%se_fbdn_extot = 0

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

       ! !Initialize forbidden excitation arrays
       ! ALLOCATE(se_box%se_fbdn(SIZE(o2_e_ex_fbdn)))
       ! se_box%se_fbdn = o2_e_ex_fbdn
       ! ALLOCATE(se_box%se_fbdnsigs(SIZE(o2_e_ex_fbdn)))
       ! se_box%se_fbdnsigs = 0
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! SE_INFO_GARBAGE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    TYPE(se_info) :: se_box

    IF ( ALLOCATED(se_box%se_ionst) .EQV. .TRUE. ) THEN
       DEALLOCATE(se_box%se_ionst)
       DEALLOCATE(se_box%se_ionsigs)
       DEALLOCATE(se_box%se_alwd)
       DEALLOCATE(se_box%se_alwdsigs)
       ! DEALLOCATE(se_box%se_fbdn)
       ! DEALLOCATE(se_box%se_fbdnsigs)
       RETURN
    ELSE
       RETURN
    END IF
  END SUBROUTINE se_info_garbage

  SUBROUTINE init_node(temp_node,x,y,z)
    !
    ! Purpose
    !   This is a subroutine that compares a string value to values
    !  in a list and gives the index of a matching result and an
    !  error if there is no match.
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! LOOKUP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE
    TYPE(node), POINTER :: temp_node
    INTEGER :: x,y,z

    temp_node%wait_time = 0.0
    temp_node%coord1 = x
    temp_node%coord2 = y
    temp_node%coord3 = z
    temp_node%sec_sp_num = 0
    temp_node%act_type = 0
    temp_node%hop_dir = 0
    temp_node%leftRight = 0
    IF ( (MOD(y,2) .EQ. 1) .AND. (MOD(z,2) .EQ. 1) ) THEN
       temp_node%sp_num = -1
    END IF
    NULLIFY(temp_node%before,temp_node%after,temp_node%parent)
  END SUBROUTINE init_node

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! NEW_REACTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE
    INTEGER, INTENT(INOUT), DIMENSION(3) :: i, j, k
    INTEGER               , DIMENSION(3) :: pr_coords
    INTEGER :: r1, r2, pr
    INTEGER :: error
    TYPE(node), POINTER :: root, temp, prevNode, nextNode

    ! Initialize variables


    SELECT CASE (k(3))
    CASE(1)
       r1 = matrix(i(1),i(2),i(3))%sp_num
       r2 = MERGE(matrix(i(1),i(2),i(3))%sec_sp_num,&
            matrix(j(1),j(2),j(3))%sp_num,&
            ALL(ABS(i-j) .EQ. 0))
       k(1) = r1
       k(2) = r2
       pr = REACT_CUBE(k(1),k(2),k(3))
       pr_coords = j
       IF ( pr .EQ. 0 ) THEN
          error = 1
          RETURN
       END IF
       k(3) = 2    
    CASE(2)
       pr = REACT_CUBE(k(1),k(2),k(3))
       IF ( pr .EQ. 0 ) THEN
          temp => matrix(i(1),i(2),i(3))
          CALL delete_node(root,temp,prevNode,nextNode)
          error = 0
          RETURN
       END IF
       pr_coords = i
       k(3) = 3
    CASE(3)
       pr = REACT_CUBE(k(1),k(2),k(3))
       IF ( pr .EQ. 0 ) THEN
          error = 0
          RETURN
       END IF
       CALL thirdman(pr_coords)
    END SELECT




    temp => matrix(pr_coords(1),pr_coords(2),pr_coords(3))
    CALL delete_node(root,temp,prevNode,nextNode)
    IF ( pr .EQ. ELECNUM ) THEN
       temp%sec_sp_num = pr
       temp%sp_num = 0
       temp%wait_time = 0
       temp%act_type = 0
       temp%hop_dir = 0
    ELSE
       temp%sp_num = pr
       temp%sec_sp_num = 0
       CALL wait_calc(temp)
       CALL add_node(root,temp)
    END IF

    CALL new_reaction(i,j,k,root,temp,prevNode,nextNode,error)
  END SUBROUTINE new_reaction

  RECURSIVE SUBROUTINE new_electron()
    IMPLICIT NONE

    curr = i
    next = j
    ion_dist = 0
    enull1 = 0
    enull2 = 0
    emfp = 0
    p = 0
    de = 0
    estep = 0
    eswitch = 0
    ee_loss = 0
    ! Calculate track until the electron's energy is depleted
    DO WHILE ( se_box%se_energy .GE. ECUTOFF )
       ! The electron's mean-free-path is a function of the total cross sections.
       ! Note: here, we have explicitly calculated the inelastic cross section and
       ! have approximated the elastic cross section to be 1.0E-17 cm^2
       CALL RANDOM_NUMBER(p)
       emfp = 1./(RHO*(se_box%se_ineltot+1.0E-17)
       ! Determine the actual distance travelled
       de = -1.*emfp*LOG(1.-p)
       ! Multiply by shortening/lengthening factor
       de = de*ESTEPFAC
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
       END DO

       ! Determine the nature of the inelastic event
       CALL RANDOM_NUMBER(erand)
       IF ( (erand .GT. 0) .AND. (erand .LE. se_box%se_iontot/se_box%se_ineltot) ) THEN
          eswitch = 1
       ELSE
          eswitch = 0
       END IF

       SELECT CASE (eswitch)
       CASE(1)
          IF ( matrix(next(1),next(2),next(3))%sp_num .NE. 0 ) THEN
             ! Calculate energy loss
             CALL e_ion_select(se_box,e_ion,enull1)
             ! Call random number
             CALL RANDOM_NUMBER(p)
             ! Calculate new electron energy based on \DeltaE
             new_e_energy = p*(se_box%se_energy - e_ion)
             ee_loss = e_ion + new_e_energy
             se_box%se_energy = se_box%se_energy - ee_loss
             ! Ionize the species and call new_electron again
          END IF
       CASE(0)       
       END SELECT
    END DO

    ! Once the electron has fallen below the energy threshold, make
    ! it react to form an anion with a surrounding species

  END SUBROUTINE new_electron


END MODULE subroutines
