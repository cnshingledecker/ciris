  MODULE subroutines
    !USE IFPORT
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

    SUBROUTINE lookaroundyou ( react_cube, matrix, coords, null, small_count, &
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
      INTEGER                         , DIMENSION(:,:,:), POINTER  :: matrix
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
      INTEGER                         , DIMENSION(3)               :: dimens



      ! Initialize dimensions
      dimens(1) = SIZE(matrix,1)
      dimens(2) = SIZE(matrix,2)
      dimens(3) = SIZE(matrix,3)

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


      ! Make sure the reactant isn't a zero
  !    PRINT *, "In lookaroundyou, the dimensions of the matrix are:"
  !    PRINT *, dimens

  !    PRINT *, "The value of coords is:"
  !    PRINT *, coords

  !    PRINT *, "The value of matrix in lookaroundyou is:",matrix(i_re,j_re,k_re)

  !    IF ( matrix(i_re,j_re,k_re) .EQ. 0 ) THEN
  !      PRINT *, "***************************************"
  !      PRINT *, "ERROR in lookaroundyou at matrix lookup"
  !      PRINT *, " Product is zero"
  !      PRINT *, "***************************************"
  !    END IF


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
            IF ( matrix(i_re2,j_re2,k_re2) .NE. 0 ) THEN
              IF ( ALL(coords .EQ. new_coords) .EQV. .FALSE.) THEN
              ! Determine if the hopped to species can react with the hopping species
                r1 = matrix(i_re,j_re,k_re)
                r2 = matrix(i_re2,j_re2,k_re2)
                CALL canreact(r1,r2,react_cube,wait_list,null)
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
              IF ( j_re-1 .GT. 0           .AND. &
                   k_re-1 .GT. 0           .AND. &
                   k_re+1 .LE. dimens(3) ) THEN
                CALL hopping(i_re,j_re-1,k_re+1,i_re2,j_re2,k_re2,1,dimens)
              ELSE
                i_re2 = i_re
                j_re2 = 1
                k_re2 = 2
              END IF
            ELSE IF ( n .EQ. 2 ) THEN
              IF ( j_re-1 .GT. 0           .AND. &
                   k_re-1 .GT. 0           .AND. &
                   k_re+1 .LE. dimens(3) ) THEN
                CALL hopping(i_re,j_re-1,k_re-1,i_re2,j_re2,k_re2,2,dimens)
              ELSE
                i_re2 = i_re
                j_re2 = 1
                k_re2 = 2
              END IF
            ELSE IF ( n .EQ. 3 ) THEN
              IF ( j_re+1 .LE. dimens(2)    .AND. &
                   k_re-1 .GT. 0            .AND. &
                   k_re+1 .LE. dimens(3) ) THEN
                CALL hopping(i_re,j_re+1,k_re+1,i_re2,j_re2,k_re2,1,dimens)
              ELSE
                i_re2 = i_re
                j_re2 = 1
                k_re2 = 2
              END IF
            ELSE
              IF ( j_re+1 .LE. dimens(2)    .AND. &
                   k_re-1 .GT. 0            .AND. &
                   k_re+1 .LE. dimens(3) ) THEN
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
            IF ( matrix(i_re2,j_re2,k_re2) .NE. 0 ) THEN
              IF ( ALL(coords .EQ. new_coords) .EQV. .FALSE. ) THEN
                r1 = matrix(i_re,j_re,k_re)
                r2 = matrix(i_re2,j_re2,k_re2)
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

    SUBROUTINE thirdman( prod,i_re,j_re,k_re, null, matrix, en_list, time, &
                         wait_list, wait_len, prod_coords)
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
      INTEGER                                           , POINTER           :: wait_len
      INTEGER            , INTENT(IN)                                       :: i_re,j_re,k_re !coords of original site
      INTEGER            , INTENT(OUT), DIMENSION(3)             , OPTIONAL :: prod_coords !product placement coords
      INTEGER                         , DIMENSION(:,:,:), POINTER           :: matrix !ice-mantle matrix
      TYPE (wait_info)                , DIMENSION(:)    , POINTER           :: wait_list
      REAL                            , DIMENSION(:)    , POINTER           :: en_list
      REAL(KIND=DBL)                                    , POINTER           :: time

      !*****************!
      ! Local variables !
      !*****************!

      INTEGER                                                               :: n,m !counters
      INTEGER                                                               :: i_pr,j_pr,k_pr !product coords
      INTEGER                         , DIMENSION(3)                        :: dimens


      ! Allocate and assign dimens pointer
      dimens(1) = SIZE(matrix,1)
      dimens(2) = SIZE(matrix,2)
      dimens(3) = SIZE(matrix,3)

      DO n=1,6
        IF ( i_re .EQ. 1 .AND. ( n .EQ. 5 .OR. n .EQ. 6 ) ) THEN
          CONTINUE
        ELSE IF ( i_re .EQ. dimens(1) .AND. n .EQ. 6 ) THEN
          ! Don't hop down if on bottom layer
          CONTINUE
        ELSE
          CALL hopping(i_re,j_re,k_re,i_pr,j_pr,k_pr,n,dimens)
          IF ( matrix(i_pr,j_pr,k_pr) .EQ. 0 ) THEN
            GOTO 1985
          END IF
        END IF
      END DO

        DO m=1,4 ! Go to a phantom position
          SELECT CASE (m)
          CASE (1)
            IF ( j_re-1 .GT. 0           .AND. &
                 k_re-1 .GT. 0           .AND. &
                 k_re+1 .LE. dimens(3) ) THEN
              CALL hopping(i_re,j_re-1,k_re+1,i_pr,j_pr,k_pr,1,dimens)
              IF (matrix(i_pr,j_pr,k_pr) .EQ. 0 ) THEN
                GOTO 1985
              ELSE
                CONTINUE
              END IF
            ELSE
              CONTINUE
            END IF
          CASE (2)
            IF ( j_re-1 .GT. 0           .AND. &
                 k_re-1 .GT. 0           .AND. &
                 k_re+1 .LE. dimens(3) ) THEN
              CALL hopping(i_re,j_re-1,k_re-1,i_pr,j_pr,k_pr,2,dimens)
              IF (matrix(i_pr,j_pr,k_pr) .EQ. 0 ) THEN
                GOTO 1985
              ELSE
                CONTINUE
              END IF
            ELSE
              CONTINUE
            END IF
          CASE (3)
            IF ( j_re+1 .LE. dimens(2)   .AND. &
                 k_re-1 .GT. 0           .AND. &
                 k_re+1 .LE. dimens(3) ) THEN
              CALL hopping(i_re,j_re+1,k_re+1,i_pr,j_pr,k_pr,1,dimens)
              IF (matrix(i_pr,j_pr,k_pr) .EQ. 0 ) THEN
                GOTO 1985
              ELSE
                CONTINUE
              END IF
            ELSE
              CONTINUE
            END IF
          CASE (4)
            IF ( j_re+1 .LE. dimens(2)   .AND. &
                 k_re-1 .GT. 0           .AND. &
                 k_re+1 .LE. dimens(3) ) THEN
              CALL hopping(i_re,j_re+1,k_re-1,i_pr,j_pr,k_pr,2,dimens)
              IF (matrix(i_pr,j_pr,k_pr) .EQ. 0 ) THEN
                GOTO 1985
              ELSE
                GOTO 2001
              END IF
            ELSE
              GOTO 2001
            END IF
          END SELECT

          2001 IF ( m .EQ. 4 ) THEN
            null = 1
  !          PRINT *, 'No reaction possible! ERROR!!!'
            RETURN
          END IF
        END DO

      ! Place reactant at chosen site
      1985 IF ( ANY( MOBILE_LIST .EQ. prod ) ) THEN
        wait_len = wait_len + 1
        wait_list(wait_len)%i = i_pr
        wait_list(wait_len)%j = j_pr
        wait_list(wait_len)%k = k_pr
        wait_list(wait_len)%sp_num = prod
        CALL wait_calc(wait_list,wait_len,en_list,time)
        matrix(i_pr,j_pr,k_pr) = wait_len
  !      PRINT *, "Placing ",prod,"at ",i_pr,j_pr,k_pr,"at index ",wait_len,&
  !               "and ",matrix(i_pr,j_pr,k_pr),"should be ",wait_len
      ELSE
        matrix(i_pr,j_pr,k_pr) = -1*prod
  !      PRINT *, "Placing ",prod,"at ",i_pr,j_pr,k_pr,&
  !               "and ",matrix(i_pr,j_pr,k_pr),"should be ",-1*prod
      END IF

      ! Assign output coordinates array
      IF ( PRESENT(prod_coords) ) THEN
        prod_coords(1) = i_pr
        prod_coords(2) = j_pr
        prod_coords(3) = k_pr
      END IF

    END SUBROUTINE thirdman

    SUBROUTINE bresenham( x1,y1,x2,y2,track )
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! BRESENHAM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !   This subroutine uses the bresenham line algorithm results to simulate the
    ! track of a cosmic ray or other form of irradiation in a solid represented by
    ! a 3D crystal lattice structure with both normal and interstitial sites.
    ! Note that here, the axes are changed such that a "slice" from top to bottom
    ! is in the x-y plane
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IMPLICIT NONE

      !*****************
      ! Input and output
      !*****************

      INTEGER, INTENT(IN)                                   :: x1,y1,x2,y2
      INTEGER, INTENT(OUT), ALLOCATABLE, DIMENSION(:,:)     :: track

      !****************
      ! Local variables
      !****************
      INTEGER                                               :: dx, dy, i, e
      INTEGER                                               :: incx,incy,inc1,inc2
      INTEGER                                               :: x,y

      dx = x2 - x1
      dy = y2 - y1

      IF ( dx .LT. 0 ) dx = -dx
      IF ( dy .LT. 0 ) dy = -dy
      incx=1
      IF ( x2 .LT. x1 ) incx = -1
      incy = 1
      IF ( y2 .LT. y1 ) incy = -1
      x = x1
      y = y1

      IF ( dx .GT. dy ) THEN
        ALLOCATE(track(dx+1,2))
  !      PRINT *, x,y
        track(1,1) = x
        track(1,2) = y
        e = 2*dy - dx
        inc1 = 2*(dy-dx)
        inc2 = 2*dy
        DO i=0,dx-1
          IF ( e .GE. 0 ) THEN
            y = y + incy
            e = e + inc1
          ELSE
            e = e + inc2
          END IF
          x = x + incx
  !        PRINT *, x,y
          track(i+2,1) = x
          track(i+2,2) = y
        END DO
      ELSE
        ALLOCATE(track(dy+1,2))
  !      PRINT *, x,y
        track(1,1) = x
        track(1,2) = y
        e = 2*dx - dy
        inc1 = 2*(dx-dy)
        inc2 = 2*dx
        DO i=0,dy-1
          IF ( e .GE. 0 ) THEN
            x = x + incx
            e = e + inc1
          ELSE
            e = e + inc2
          END IF
          y = y + incy
  !        PRINT *, x,y
          track(i+2,1) = x
          track(i+2,2) = y

        END DO
      END IF
    END SUBROUTINE bresenham

    SUBROUTINE cern ( o3_prod,o3_dest,null,en_list, react_cube, matrix, event_num, &
                      event_coords, switch, wait_list, wait_len, time, elec_coords )
    !
    ! Purpose:
    !   This subroutine handles interaction events between ionizing radiation and a
    !  target species in a solid
    !
    ! Note:
    !   The array event_num contains the species number for pseudo reactants,
    !  i.e. excitations, excitations, and electrons. Specifically:
    !
    ! -- event_num(1) = exc_num
    ! -- event_num(2) = ion_num
    ! -- event_num(3) = e_num
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! CERN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IMPLICIT NONE

      !******************!
      ! Input and output !
      !******************!
      INTEGER                                               , POINTER          :: o3_prod,o3_dest
      INTEGER(KIND=SHORT), INTENT(INOUT)                                       :: null
      INTEGER            , INTENT(IN)    , DIMENSION(3)                        :: event_coords !coordinated of collision
      INTEGER            , INTENT(IN)                                          :: switch !type of event
      INTEGER            , INTENT(IN)    , DIMENSION(3)                        :: event_num !pseudo species numbers
      INTEGER(KIND=SHORT)                , DIMENSION(:,:,:), POINTER           :: react_cube !cube of reactions
      INTEGER                            , DIMENSION(:,:,:), POINTER           :: matrix !ice-mantle matrix
      INTEGER                                              , POINTER           :: wait_len
      INTEGER            , INTENT(OUT)   , DIMENSION(3)             , OPTIONAL :: elec_coords
      REAL(KIND=DBL)                                       , POINTER           :: time
      REAL                               , DIMENSION(:)    , POINTER           :: en_list !list of binding energies
      TYPE(wait_info)                    , DIMENSION(:)    , POINTER           :: wait_list

      !*****************!
      ! Local variables !
      !*****************!
      INTEGER                                                                  :: i,j,k
      INTEGER                                                                  :: i_re2,j_re2,k_re2 !reactant coordinates
      INTEGER                                                                  :: i_pr,j_pr,k_pr !product coordinates
      INTEGER                                                                  :: index,index2
      INTEGER                                                                  :: matrix_num
      INTEGER                                                                  :: r1, r2 !reactants 1 and 2
      INTEGER                                                                  :: prods_case
      INTEGER                            , DIMENSION(3)                        :: prods !array of reaction products
      INTEGER                            , DIMENSION(3)                        :: coords !temporary coordinates array
      INTEGER                            , DIMENSION(3)                        :: dimens
      INTEGER :: n
      INTEGER :: original_value
      REAL    :: rnum
      INTEGER                                                                  :: case_num

      ! DEBUG = .TRUE.
      IF ( DEBUG .EQV. .TRUE. ) PRINT *, '*****STARTING Cern*****'

      r1 = 0
      r2 = 0
      i_pr = 0
      j_pr = 0
      k_pr = 0

      ! Assign the coordinates
      i_re2 = event_coords(1)
      j_re2 = event_coords(2)
      k_re2 = event_coords(3)
      original_value = matrix(i_re2,j_re2,k_re2)

      IF ( DEBUG .EQV. .TRUE. ) THEN
        PRINT *, 'Event_coords are:', event_coords
        PRINT *, 'Original value is:',original_value
        PRINT *, 'Switch is:',switch
      END IF


      ! If switch equals 2, then the first reactant, r1, equals an excitation,
      ! otherwise, it equals and ionizing CRP

      SELECT CASE (switch)
      CASE (2) ! Ionization
        r1 = event_num(2)
      CASE (1) ! Excitation
        r1 = event_num(1)
      END SELECT

      index = 0
      matrix_num = matrix(i_re2,j_re2,k_re2)

      IF ( matrix_num .LT. 0 ) THEN
        r2 = ABS(matrix_num)
      ELSE IF ( matrix_num .GT. 0 ) THEN
        r2 = wait_list(matrix_num)%sp_num
        index2 = matrix_num
  !      IF ( index .EQ. 53 ) PRINT *, 'The wait_list at',index,'is:'
  !      IF ( index .EQ. 53 ) PRINT *, wait_list(index)
      END IF

      !****************************************************************************
      !**** ERROR CHECKING ********************************************************
      !****************************************************************************
      IF ( r2 .EQ. 0 ) THEN
  !      matrix(event_num(1),event_num(2),event_num(3)) = 0
  !      PRINT *, 'WHOOPS! FIXING ERROR!'
        PRINT *, 'r1 =',r1
        PRINT *, 'r2 =',r2
        PRINT *, 'original_value =',original_value
        PRINT *, 'event_coords =',event_coords
        dimens(1) = SIZE(matrix,1)
        dimens(2) = SIZE(matrix,2)
        dimens(3) = SIZE(matrix,3)
        OPEN(UNIT=1013,FILE="cern_test_wrong_spaces.txt")
        DO k=1,dimens(3)
          DO j=1,dimens(2)
            DO i=1,dimens(1)
              IF ( matrix(i,j,k) .GT. 0 ) THEN
                IF( ANY(MOBILE_LIST .NE. wait_list(matrix(i,j,k))%sp_num ) ) THEN
                  WRITE(1013,*)  'Matrix=',matrix(i,j,k),'wait_list=',wait_list(matrix(i,j,k))
                  PRINT *, 'We have a wrong space:',matrix(i,j,k),' at',i,j,k
                END IF
              END IF
            END DO
          END DO
        END DO
        CLOSE(1013)
        PRINT *, "BEEP BOOP 1 0"
        PRINT *, "ERROR: wait_len =",wait_len,"< matrix_num =",matrix_num
        CALL EXIT()
      END IF

      prods = 0
      prods = react_cube(r1,r2,:) ! populate the product array

      !****************************************************************************
      !**** BRANCHING RATIOS ******************************************************
      !****************************************************************************
      !Determine if there is dissociation
      IF ( r1 .EQ. 7 .AND. r2 .EQ. event_num(1) .OR. r1 .EQ. event_num(1) .AND. r2 .EQ. 7 ) THEN
    !   PRINT *, "Weve got branching"
!        rnum = RAND()
        CALL RANDOM_NUMBER(rnum)
        IF ( rnum .GT. O3_DIS_BRANCHING ) THEN
          prods = (/ 1, 4, 0 /)
        END IF
      END IF

      !****************************************************************************
      !****ANALYTICS***************************************************************
      !****************************************************************************
      IF ( r1 .EQ. 7 .OR. r2 .EQ. 7 ) o3_dest = o3_dest + 1
      IF ( ANY(prods .EQ. 7 ) ) o3_prod = o3_prod + 1

      !****************************************************************************
      !****ERROR CHECKING**********************************************************
      !****************************************************************************
      IF ( DEBUG .EQV. .TRUE. ) THEN
        PRINT *, "wait_len in cern is ",wait_len
        PRINT *,     r1, r2, prods
        WRITE(777,*) r1,',',r2,',',prods(1),',',prods(2),',',prods(3), ',',wait_len
        IF ( r1 .EQ. 7 .OR. r2 .EQ. 7 ) THEN
          PRINT *, '\/ In cern \/'
          PRINT *, r1,',',r2,',',prods(1),',',prods(2),',',prods(3)
        END IF

        IF ( ANY(prods .EQ. 7 ) ) THEN
          PRINT *, '\/ In cern \/'
          PRINT *, r1,',',r2,',',prods(1),',',prods(2),',',prods(3)
        END IF
      END IF

      IF ( prods(1) .EQ. 0 ) THEN
        PRINT *, 'Uh-oh, we have a problem in Cern!'
        PRINT *, 'r1=',r1,'r2=',r2
        PRINT *, 'prods=',prods
        OPEN(UNIT=1013,FILE="test_wait_list.txt")
        DO n=1,wait_len+2
          WRITE(1013,*) wait_list(n)
        END DO
        CLOSE(1013)

        OPEN(UNIT=1014,FILE="wait_list_flaw.txt")
        DO n=1,wait_len
          IF ( wait_list(n)%sp_num .NE. 20 ) THEN
            IF ( matrix(wait_list(n)%i,wait_list(n)%j,wait_list(n)%k) .NE. n ) THEN
              WRITE(1014,*) matrix(wait_list(n)%i,wait_list(n)%j,wait_list(n)%k),", ",n,",",wait_len
            END IF
          END IF
        END DO
        CLOSE(1014)
        CALL EXIT()
      END IF

      !****************************************************************************
      !**** SANITY CHECK **********************************************************
      !****************************************************************************
      ! If there are two or more products: do a sanity check to see if there are
      ! sufficient empty sites
      !****************************************************************************
      original_value = matrix(i_re2,j_re2,k_re2)
      null = 0
      IF ( prods(2) .NE. 0 ) THEN
      ! NB: Unlike in the reaction subroutine, in Krell, one needs to
      ! find a suitable location for the second product
      ! NB: In the event of an ionization, the electron should ALWAYS be
      ! stored as prods(2), doing otherwise will result in errors
        CALL krell( event_coords, coords, matrix,null )
        IF (null .EQ. 1 ) THEN ! Can't place 2nd product: no good sites
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !TEMPORARY KLUDGE
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! DEBUG = .FALSE.
          RETURN
        ELSE ! Save new coordinates
          i_pr = coords(1)
          j_pr = coords(2)
          k_pr = coords(3)
        END IF
      END IF

      !*****************************************************************************
      ! Save the coords of the electron, if necessary
      !*****************************************************************************
      IF ( prods(2) .EQ. event_num(3) .AND. PRESENT(elec_coords) ) THEN
        elec_coords(1) = i_pr
        elec_coords(2) = j_pr
        elec_coords(3) = k_pr
      END IF

      !*****************************************************************************
      ! Determine the product case
      !*****************************************************************************
      CALL rtype_tree(r1,r2,prods,case_num)

      !*****************************************************************************
      ! Place the products
      !*****************************************************************************
      CALL place_product(i_re2,j_re2,k_re2,i_pr,j_pr,k_pr,index,index2,r1,r2,&
                         prods,case_num,matrix, wait_list,wait_len,en_list,time)


      IF ( DEBUG .EQV. .TRUE. ) PRINT *, '*****ENDING Cern*****'
      ! DEBUG = .FALSE.
    END SUBROUTINE cern

    SUBROUTINE fallout ( o3_prod,o3_dest,react_cube, matrix,  en_list, ionlist, &
                         wait_list, wait_len, time, ev_nums)
    ! Purpose:
    !   To calculate the track of a particle of ionizing radiation through a
    !  crystaline solid.
    !
    ! Note:
    !   The input array, sigmas, contains the proton collision cross-sections.
    !
    ! Note:
    !   The array ev_nums contains the following values for pseudo-reactants:
    !
    ! -- ev_nums(1) = cosmic-ray species number
    ! -- ev_nums(2) = excitation species number
    ! -- ev_nums(3) = electron species number
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
    INTEGER                              , POINTER :: wait_len
    INTEGER            , DIMENSION(:)    , POINTER :: ionlist !list of anionic species
    INTEGER(KIND=SHORT), DIMENSION(:,:,:), POINTER :: react_cube !array of products/reactions
    INTEGER            , DIMENSION(:,:,:), POINTER :: matrix !ice-mantle matrix
    REAL               , DIMENSION(:)    , POINTER :: en_list !list of binding and desorption energies
    REAL(KIND=DBL)                       , POINTER :: time
    TYPE(wait_info)    , DIMENSION(:)    , POINTER :: wait_list
    INTEGER            , DIMENSION(3)              :: ev_nums !values of special pseudoreactants

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
    INTEGER            , DIMENSION(3)              :: dimens !dimensions of matrix
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
      CLOSE(2016)
      OPEN(UNIT=2016,FILE="trackplot.csv", STATUS='REPLACE')
    END IF


    !****************************************************************************!
    ! Preliminary  calculations                                                  !
    !****************************************************************************!

    count_count = 0
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
!    PRINT *, "The initial proton cross-secions are:"
!    DO n=1,3
!      PRINT *, psigmas(n)
!    END DO

    ! Calculate the total cross-section and the mean free path
    ! Note, sigma_i is the inelastic ionization cross-section and
    ! sigma_e is the inelastic excitation cross section
    sigma_tot = SUM(psigmas%cross_section)
    mfp       = 1./(rho*sigma_tot)

  ! Get dimensions of matrix
    dimens(1) = SIZE(matrix,1)
    dimens(2) = SIZE(matrix,2)
    dimens(3) = SIZE(matrix,3)
!   PRINT *, 'Dimens are:',dimens

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
        IF ( x .GT. 300 .AND. x .LT. 600 .AND. y .GT. 300 .AND. y .LT. 600 ) proceed = .TRUE.
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


    main_loop: DO WHILE (z .LE. dimens(1) .AND. ione .GE. 5.0 )
      IF ( DEBUG .EQV. .TRUE. ) PRINT *, 'Now entering loop: z=',z,' and dimens(1)=',dimens(1),' and step=',step
      count_count = count_count + 1
!      IF ( MOD(count_count,1000) .EQ. 0 ) CALL counter(time, AB_UNIT_NUM, matrix, wait_list, 4,7)
      ! Define event coords
      ev_coords(1) = z+step
      ev_coords(2) = y
      ev_coords(3) = x
      IF ( TRACKPLOT .EQV. .TRUE. ) THEN
        count_count = count_count + 1
        WRITE(2016,*) ev_coords(1),',',ev_coords(2),',',ev_coords(3)
      END IF
!      PRINT *, 'The event coords are:',ev_coords

      thinghit = matrix(ev_coords(1),ev_coords(2),ev_coords(3))


      IF ( matrix(z+step,y,x) .NE. 0 ) THEN
        IF ( DEBUG .EQV. .TRUE. ) PRINT *, "The value of the matrix is:",matrix(z+step,y,x)
        ! If the site is occupied, then determine the type of event to occur
!        DO WHILE ( u .EQ. 0.0 .AND. rand .EQ. 0.0 )
!        u = RAND()
!        rand1 = RAND()
        CALL RANDOM_NUMBER(u)
        CALL RANDOM_NUMBER(rand1)
!        END DO


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
          IF ( u .GT. 0.0 .AND. u .LE. (sigma_i + sigma_e)/sigma_tot ) THEN
            IF ( u .GT. 0 .AND. u .LE. sigma_i/(sigma_i + sigma_e) ) THEN
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
        END ASSOCIATE
!        PRINT *, "Ion energy:",ione, "loss:",e_loss,'nature:',nature
!        WRITE(10,*) ione,",",e_loss,",",ione-e_loss,',',nature
        ione = ione - e_loss
        CALL psigma_suite(ione,psigmas,psigij,psigexj)

        !Debugging
!        switch = 2
!        PRINT *, 'The value of the switch is:',switch
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
          CALL cern( o3_prod,o3_dest,null,en_list, react_cube, matrix, ev_nums, ev_coords, switch, &
                     wait_list, wait_len, time )

          IF ( null .EQ. 1 ) RETURN
  !****************************************************************************!
  ! Ionization                                                                 !
  !****************************************************************************!
        ELSE IF ( switch .EQ. 2 .AND. z+step .NE. 1 .AND. z+step .NE. 2 ) THEN
        !  PRINT *, 'SE energy is',se_box%se_energy
          IF ( se_box%se_energy .LE. ECUTOFF ) THEN
          !  PRINT *, 'SE energy is less than or equal to cutoff'
!            PRINT *, 'ev_coords=',ev_coords
            CALL base_ionization( o3_prod,o3_dest,ev_coords, react_cube, matrix, en_list, &
                                  ionlist, wait_list, wait_len, &
                                  time, ev_nums, null )
            IF ( null .EQ. 1 ) RETURN
          ELSE
            ! PRINT *, 'SE energy above ecutoff' !*******************************************************************
            !
            ! Generate secondary electrons/electron track
            !
            !*******************************************************************
            ! Call base_ionization to generate the first-generation secondary electron
            ! NB: the electron should be the second product in the "prods" array
            !*******************************************************************
            null = 0
            CALL base_ionization( o3_prod,o3_dest,ev_coords, react_cube, matrix, en_list, &
                                  ionlist, wait_list, wait_len, &
                                  time, ev_nums, null,elec_coords )

            IF ( TRACKPLOT .EQV. .TRUE. ) THEN
              count_count = count_count + 1
              WRITE(2016,*) elec_coords(1),',',elec_coords(2),',',elec_coords(3)
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
            CALL se_info_init(se_box)
            ! PRINT *, 'se energy is:',se_box%se_energy,' and ECUTOFF is',ECUTOFF
            DO WHILE ( se_box%se_energy .GE. ECUTOFF )
              ! PRINT *, 'Electron box initialized, calling loop'
              !If the electron no longer has sufficient
              !energy, exit the loop.
              IF ( enull1 .NE. 0 .AND. enull2 .NE. 0 ) EXIT

              !Calculate electron cross_sections
              CALL esigma_suite(se_box)

              !Calculate hopping distance
              emfp  = 1./(RHO*(se_box%se_ineltot))
!              p     = RAND()
              CALL RANDOM_NUMBER(p)
              de    = -1.*emfp*LOG(1.-p)
              estep = INT(de/C_PR)

              !Have a minumum hopping distance of 1
              IF ( estep .EQ. 0 ) estep = 1

              ! PRINT *, 'In SE routine, curr=',curr

              DO n=1,estep
                !Each transport hop is like one step
                prev = curr
                curr = next
                CALL transport(prev,curr,next,matrix)

                IF ( TRACKPLOT .EQV. .TRUE. ) THEN
                  count_count = count_count + 1
                  WRITE(2016,*) curr(1),',',curr(2),',',curr(3)
                END IF
              END DO

              !Determine nature of event
!              erand = RAND()
              CALL RANDOM_NUMBER(erand)
              IF ( erand .GT. 0 .AND. erand .LE. (se_box%se_iontot/se_box%se_ineltot) ) THEN
                eswitch = 1
              ELSE
                eswitch = 0
              END IF

              !Carry out impact collision
              IF ( matrix(next(1),next(2),next(3)) .NE. 0 .AND. eswitch .EQ. 1 ) THEN
                !Electron impact ionization
                IF ( DEBUG .EQV. .TRUE. ) PRINT *, 'EII: SE is hopping to site with',matrix(next(1),next(2),next(3))
                CALL base_ionization(o3_prod,o3_dest,next, react_cube, matrix, en_list, &
                                     ionlist, wait_list, wait_len, &
                                     time, ev_nums,null )
                IF ( null .EQ. 1 ) GOTO 100
                CALL e_ion_select(se_box,e_ion,enull1)
                ee_loss = e_ion
              ELSE IF ( matrix(next(1),next(2),next(3)) .NE. 0 .AND. eswitch .EQ. 0 ) THEN
              !Electron impact excitation
                IF ( DEBUG .EQV. .TRUE. ) PRINT *, 'EIE: SE is hopping to site with',matrix(next(1),next(2),next(3))
!                rand1 = RAND()
                CALL RANDOM_NUMBER(rand1)
                IF ( rand1 .LE. DISPROB .AND. matrix(curr(1),curr(2),curr(3)).NE. 0 ) THEN
                  CALL cern( o3_prod,o3_dest,null,en_list, react_cube, matrix, ev_nums, next, 1, &
                             wait_list, wait_len, time )
                ELSE IF ( ANY( FRAGILE .EQ. -1*matrix(curr(1),curr(2),curr(3))  ) ) THEN
                  ! Test for fragile species
                  CALL cern( o3_prod,o3_dest,null,en_list, react_cube, matrix, ev_nums, next, 1, &
                             wait_list, wait_len, time )
                END IF
                CALL e_ex_select(se_box,e_exc,enull2)
                ee_loss = e_exc
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
     !       PRINT *, 'Starting sub_excitation process'
            IF ( NSUBEX .NE. 0 ) THEN
              nn = 0
              exitcount = 0

              DO WHILE ( nn .LT. NSUBEX .AND. exitcount .LT. NEXIT) !NSUBEX is the number of sub-excitation collisions
                IF ( DEBUG .EQV. .TRUE. ) THEN
                  PRINT *, 'nn=',nn
                  PRINT *, 'exitcount=',exitcount
                END IF
                exitcount = exitcount + 1
                !DO jj=1,estep
                  prev = curr
                  curr = next
                  CALL transport(prev,curr,next,matrix)

                  IF ( TRACKPLOT .EQV. .TRUE. ) THEN
                    count_count = count_count + 1
                    WRITE(2016,*) curr(1),',',curr(2),',',curr(3)
                  END IF
                !END DO
                IF ( matrix(next(1),next(2),next(3)) .NE. 0 ) THEN
!		  subexrand = RAND()
!		  IF ( subexrand .GT. SUBEXHITPROB ) THEN
                      IF ( DEBUG .EQV. .TRUE. ) PRINT *, 'matrix in subexloop=',matrix(next(1),next(2),next(3))
                      !Carry out dissociate electron attachment
                      !NB: In the model, this is functionally identical to
                      !an ordinary ionization
                      CALL base_ionization(o3_prod,o3_dest,next, react_cube, matrix, en_list, &
                                           ionlist, wait_list, wait_len, &
                                           time, ev_nums,null )
!		  END IF
                  nn = nn + 1
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
!      PRINT *,'Before writing to file, the coordinates are:'
!      PRINT *, 'i_pr =',i_pr,'This should be a large number'
!      PRINT *, 'j_pr =',j_pr
!      PRINT *, 'k_pr =',k_pr

!      IF ( switch .EQ. 2 ) THEN
!        PRINT *, 'Writing to file:'
!        WRITE (10,*) z+step, ', 100'
!      END IF

      ! Increment z for next cycle
      z = z + step


      !********************!
      ! GOTO jumps to here !
      !********************!
      ! Call a random number between [0,1)
      100 CALL RANDOM_NUMBER(p)

      ! Make sure the random number does not equal 1
      IF ( p .EQ. 1.0 ) THEN
        DO
          IF ( p .NE. 1.0 ) EXIT
          CALL RANDOM_NUMBER(p)
        END DO
      END IF

      ! Determine the travel distance
      dz = -mfp*LOG(1-p)
!      PRINT *, 'The travel distance is:',dz,'m'

      ! Determine whether or not the site is occupied by dividing the
      ! Delta z by the height of the crystal cube, i.e. \Delta ml =
      ! \Delta z(m) * (1ml/c(m))
      step = INT(STEPFAC*(dz/c_pr))

      ! Make sure the next site is different than the previous one
      IF ( z+step .EQ. z ) THEN
        z = z + 1
      END IF

      ! Make sure the site is greater than the previous one
      IF (step .LE. 0. ) GOTO 100
      IF (z+step .GE. dimens(1) ) THEN
        IF ( TRACKPLOT .EQV. .TRUE. ) THEN
          IF ( count_count .GT. TRACKMAX ) CALL EXIT()
        ELSE
          RETURN
        END IF
      END IF

      IF ( TRACKPLOT .EQV. .TRUE. ) THEN
        PRINT *, 'Count_count=',count_count
        !IF ( count_count .GT. 80000 ) CALL EXIT()
      END IF
      IF ( DEBUG .EQV. .TRUE. ) PRINT *, 'Now at end of main loop in Fallout'
    END DO main_loop

    IF ( DEBUG .EQV. .TRUE. ) PRINT *, '*****Ending Fallout*****'
!    CLOSE(10)
  END SUBROUTINE fallout

  SUBROUTINE krell( in_coords,out_coords,matrix,null )
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
    INTEGER                         , DIMENSION(:,:,:), POINTER  :: matrix !ice-mantle matrix
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
    INTEGER                         , DIMENSION(3)               :: dimens !dimensions of ice matrix
    INTEGER                         , DIMENSION(6,3)             :: large_temp
    INTEGER                         , DIMENSION(4,3)             :: small_temp
    INTEGER            , ALLOCATABLE, DIMENSION(:,:)             :: temp_arr !temporary empty site array
    REAL                                                         :: rand !random number

!    PRINT *, 'In krell, in_coords=',in_coords

    ! Obtain the dimensions of the matrix
    dimens(1) = SIZE(matrix,1)
    dimens(2) = SIZE(matrix,2)
    dimens(3) = SIZE(matrix,3)

!    PRINT *, 'In krell, dimens=',dimens

    ! Initialize counters and arrays
    large_count   = 0
    small_count  = 0

!    PRINT *,'The coordinates of the event are',in_coords

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
        CALL hopping(i_re,j_re,k_re,i_re2,j_re2,k_re2,n,dimens )
        IF ( matrix(i_re2,j_re2,k_re2) .EQ. 0 ) THEN
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
!      PRINT *, 'Phantom hopping'
      DO n=1,4
      ! Go to a phantom position to hop to nearest neighbors
        SELECT CASE (n)
        CASE (1)
          IF ( j_re-1 .GT. 0          .AND. &
               k_re-1 .GT. 0          .AND. &
               k_re+1 .LE. dimens(3) ) THEN
            CALL hopping(i_re,j_re-1,k_re+1,i_re2,j_re2,k_re2,1,dimens)
          ELSE
            GOTO 1944
          END IF

        CASE (2)
          IF ( j_re-1 .GT. 0          .AND. &
               k_re-1 .GT. 0          .AND. &
               k_re+1 .LE. dimens(3) ) THEN
            CALL hopping(i_re,j_re-1,k_re-1,i_re2,j_re2,k_re2,2,dimens)
          ELSE
            GOTO 1944
          END IF

        CASE (3)
          IF ( j_re+1 .LE. dimens(2)  .AND. &
               k_re-1 .GT. 0          .AND. &
               k_re+1 .LE. dimens(3) ) THEN
            CALL hopping(i_re,j_re+1,k_re-1,i_re2,j_re2,k_re2,2,dimens)
          ELSE
            GOTO 1944
          END IF

        CASE (4)
          IF ( j_re+1 .LE. dimens(2)  .AND. &
               k_re-1 .GT. 0          .AND. &
               k_re+1 .LE. dimens(3) ) THEN
            CALL hopping(i_re,j_re+1,k_re+1,i_re2,j_re2,k_re2,1,dimens)
          ELSE
            GOTO 1944
          END IF
        END SELECT

        IF ( matrix(i_re2,j_re2,k_re2) .EQ. 0 ) THEN
          small_temp(n,1)=i_re2
          small_temp(n,2)=j_re2
          small_temp(n,3)=k_re2
          small_count = small_count + 1
        END IF
        1944 CONTINUE
      END DO
    END IF

    ! Determine if there has been a null event
!    PRINT *, 'Determining null event'
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
!      PRINT *, 'The contents of int_arr are:'
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

!    PRINT *, 'The contents of temp_arr are:'
!    DO n=1,i-1
!      PRINT *, temp_arr(n,1), temp_arr(n,2), temp_arr(n,3)
!    END DO

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

!    PRINT *, 'The coordinates of the second site are:'
!    PRINT *, out_coords
  END SUBROUTINE krell

  SUBROUTINE ioncount(species_file,anion_num,cation_num)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! IONCOUNT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subroutine counts the number of ions in a species list file. It can be
  ! set to count either cations or anions
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE
    INTEGER                                        :: n ! Counter
    INTEGER          , INTENT(OUT)                 :: anion_num
    INTEGER          , INTENT(OUT)                 :: cation_num
    INTEGER                                        :: ierror1
    INTEGER                                        :: headerlines
    INTEGER(KIND=SHORT)                            :: nlines1
    REAL             , ALLOCATABLE, DIMENSION(:)   :: energy_array
    CHARACTER(len=10), ALLOCATABLE, DIMENSION(:)   :: speciesList
    CHARACTER(len=80), INTENT(in)                  :: species_file

    ! Open files
    OPEN (UNIT=1,FILE=species_file,STATUS='OLD',ACTION='READ',IOSTAT=ierror1)

    IF ( ierror1 .EQ. 0 ) THEN
  !    PRINT *, 'The files have been opened.'

      ! Count the number of lines in the files
      CALL linecount(1,ierror1,nlines1,headerlines)
      ALLOCATE( speciesList(nlines1-headerlines) )
      ALLOCATE( energy_array(nlines1-headerlines) )

      ! Read the contents of the species file and create the speciesList and
      ! energy_list
      anion_num = 0
      cation_num = 0
      DO n=1,nlines1
        READ(1,*,IOSTAT=ierror1) speciesList(n), energy_array(n)
  !      PRINT *, speciesList(n)

  !      tempName = TRIM(speciesList(n))
  !      PRINT *, tempName
  !      IF ( tempName(LEN(TRIM(tempName)):LEN(TRIM(tempName))) == '-' ) THEN
  !        PRINT *, TRIM(tempName), 'is an anion'
  !        anion_num = anion_num + 1
  !      ELSE IF ( tempName(LEN(TRIM(tempName)):LEN(TRIM(tempName))) == '+' ) THEN
  !        PRINT *, TRIM(tempName), 'is an cation'
  !        cation_num = cation_num + 1
  !      END IF

        IF ( ierror1 .NE. 0 ) STOP "Error reading species file."
      END DO
    ELSE
      PRINT *, 'Unable to open file...'
      CALL EXIT()
    END IF

  END SUBROUTINE ioncount

  SUBROUTINE roll_call ( wait_list, time, wait_len, mindex )
  !
  ! Purpose:
  !  The purpose of this subroutine is to take the
  !  waiting list, sort it, and "skim" off the first
  !  row, returning all of those values but the first,
  !  i.e. the old waiting time (which will have expired
  !  since it is at the top of the list.
  !
  ! Note:
  !  The types in wait_list are real, but the coords and
  !  species number are integers
  !
  ! Documentation:
  !  DATE          PROGRAMMER           DESCRIPTION
  !  ========      ==========           ===========
  !  20150415      C. Shingledecker     Original code
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! ROLL_CALL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    ! Data dictionary: variables passed to and by the subroutine
    INTEGER            , INTENT(OUT)                          :: mindex !index of smallest waiting time row
    INTEGER                                       , POINTER   :: wait_len !length of non-zero entries in wait_list
    REAL(KIND=DBL)                                , POINTER   :: time !total elapsed time
    TYPE(wait_info)                 , DIMENSION(:), POINTER   :: wait_list !waiting list of mobile species/events

    ! Data dictionary: local variables
    INTEGER                         , DIMENSION(3)            :: coords
    INTEGER                                                   :: species
    INTEGER                         , DIMENSION(1)            :: ind

    mindex = 0
  !  CALL minmod(wait_list,mindex,wait_len)
    ind =  MINLOC(wait_list(1:wait_len)%wait_time)
    mindex = ind(1)


    ! (2) Read the first element and save to output
    coords(1) = wait_list(mindex)%i
    coords(2) = wait_list(mindex)%j
    coords(3) = wait_list(mindex)%k

  !  PRINT *, "The coords are:",coords
    species   = wait_list(mindex)%sp_num
  !  PRINT *, "The species is:",species

    ! (3) Increment total simulation time
    time = wait_list(mindex)%wait_time
  END SUBROUTINE roll_call

  SUBROUTINE action_figure ( wait_list, index, rand_num, en_list )
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
    INTEGER       , INTENT(IN)                          :: index !index of species to
    REAL(KIND=DBL)                                      :: b_1 !thermal surface hopping rate
    REAL(KIND=DBL)                                      :: b_2 !surface desorption rate
    REAL(KIND=DBL)                                      :: comp_val !to determine which action occurs
    REAL(KIND=DBL), INTENT(IN)                          :: rand_num
    REAL                      , DIMENSION(:)  , POINTER :: en_list
    TYPE (wait_info)          , DIMENSION(:)  , POINTER :: wait_list

    ! (1) Decide whether or not the species is on the surface
    IF ( wait_list(index)%i .EQ. 1 ) THEN
    ! (1a) Species is on the surface
      ! Calculate b-rates to compare
      b_1 = trl_nu * EXP( - (  en_list(wait_list(index)%sp_num)*E_SURF / kin_temp  ) )
      b_2 = trl_nu * EXP( - (  en_list(wait_list(index)%sp_num)        / kin_temp  ) )
      comp_val = b_1 / (b_1 + b_2)
      ! Decide whether desorption or hopping occurs
      IF ( rand_num .LT. comp_val ) THEN
        ! Diffusion occurs
        wait_list(index)%act_type = 1
      ELSE
        ! Desorption occurs
        wait_list(index)%act_type = 2
      END IF
    ELSE
    ! (1b) Species is in the bulk
      ! Only hopping (diffusion) can occur
      wait_list(index)%act_type = 1
    END IF

    ! If the species is in fast reacting
    IF( ANY(FAST_REACTS .EQ. wait_list(index)%sp_num) ) THEN
      wait_list(index)%act_type = 3
    END IF
  END SUBROUTINE action_figure

  SUBROUTINE meta_hop ( o3_prod,o3_dest,mindex, react_ptr, matrix_ptr, en_list, wait_list, &
                        result, wait_len, time )
  !
  ! Purpose:
  !  The purpose of this subroutine is to serve as
  !  a central calling function for the subroutines
  !  involved in a reaction.This subroutine is called
  !  any time there is a reaction in the bulk or on
  !  the surface.
  !
  ! Note:
  !  The subroutine checks whether or not the species is
  !  on the surface. If the species is on the surface
  !  i.e. coords(1) = 1, then the species can
  !  only hop in one of four possible directions.
  !  On the other hand, bulk species can hop in
  !  one of six possible directions since they
  !  can also hop up or down on an interstitial
  !  site.
  !
  ! Note:
  !  To convert a real random number in interval [0,1) to
  !  an integer in the interval [n,m], use the formula:
  !
  !  j = n + FLOOR( (m+1-n) * rand_num )
  !
  !  For intervals where n = 0, this simplifies to
  !
  !  j = FLOOR( (m+1) * rand_num )
  !
  ! Documentation:
  !  DATE          PROGRAMMER           DESCRIPTION
  !  ========      ==========           ===========
  !  20150416      C. Shingledecker     Original code
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! META_HOP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    ! Data dictionary: variables passed to and by the subroutine
    INTEGER                                           , POINTER :: o3_prod,o3_dest
    INTEGER            , INTENT(IN)                             :: mindex !index of the moving species in wait_list
    INTEGER                                           , POINTER :: wait_len !length of non-zero elements of wait_list
    INTEGER(KIND=SHORT)             , DIMENSION(:,:,:), POINTER :: react_ptr !pointer to the "cube" of reactions
    INTEGER                         , DIMENSION(:,:,:), POINTER :: matrix_ptr !pointer to the matrix
    REAL                            , DIMENSION(:)    , POINTER :: en_list !energy list
    REAL(KIND=DBL)                                    , POINTER :: time !total sim time
    TYPE(wait_info)                 , DIMENSION(:)    , POINTER :: wait_list

    ! Data dictionary: local variables
    INTEGER                                                     :: index2 !index of second reactant in wait_list
    INTEGER                                                     :: species !species moving
    INTEGER                                                     :: species2 !potential reactant at new site
    INTEGER                                                     :: hop_dir !governs hop direction
    INTEGER                                                     :: i, j, k !coordinates of hopping species
    INTEGER                                                     :: i_hop, j_hop, k_hop !coords of new location
    INTEGER                                                     :: result !indicates whether hopping or reaction occurs
    INTEGER                         , DIMENSION(3)              :: dimens !dimension of matrix
    REAL                            , DIMENSION(3,5)            :: prods !product array
    REAL                                                        :: rand_num !random number

  !  OPEN(UNIT=1015,FILE="mindex_list.txt",POSITION="append")
!    PRINT *, 'In meta_hop, wait_list at index is:'
 !   PRINT *, wait_list(mindex)

    ! Initialize the species
    species = wait_list(mindex)%sp_num
!    PRINT *, 'species in meta_hop is',species

    ! Initialize prods array to 0
    prods = 0

    ! Allocate the dimensions of the matrix pointer
    dimens(1) = SIZE(matrix_ptr,1)
    dimens(2) = SIZE(matrix_ptr,2)
    dimens(3) = SIZE(matrix_ptr,3)


    ! Call a random number
    CALL RANDOM_NUMBER(rand_num)

    ! Check whether species is on the surface or not,
    ! to determine direction of hopping, then save
    ! the coordinates of the new site
    IF ( wait_list(mindex)%i .EQ. 1 ) THEN
      ! Species is on surface
!      PRINT *, "Surface species"
      ! integer in interval [1,4]
      hop_dir = 1 + FLOOR(4*rand_num)
  !    PRINT *, "The surface hopping direction is:", hop_dir
    ELSE IF ( wait_list(mindex)%i .EQ. dimens(1) ) THEN
!      PRINT *, "Species on bottom"
      ! integer in interval [1,5]
      hop_dir = 1 + FLOOR(5*rand_num)
    ELSE
!      PRINT *, 'Species in bulk'
      ! Species is in bulk
      ! Convert randon number to one in interval [1,6]
      hop_dir = 1 + FLOOR(6*rand_num)
    END IF

    i = wait_list(mindex)%i
    j = wait_list(mindex)%j
    k = wait_list(mindex)%k
    i_hop = 0
    j_hop = 0
    k_hop = 0
    IF ( DEBUG .EQV. .TRUE. ) PRINT *, "Initially, i, j, k are ",i,j,k    ! Call hopping to new location
    IF ( DEBUG .EQV. .TRUE. ) PRINT *, 'Initially, at such points, matrix=',matrix_ptr(i,j,k)
      CALL hopping(i,j,k, i_hop, j_hop, k_hop, hop_dir, dimens)
  !    PRINT *, "Hop dir is ",hop_dir
!      PRINT *, "i_hop,j_hop,k_hop are ",i_hop,j_hop,k_hop

      index2 = matrix_ptr(i_hop,j_hop,k_hop)
!      PRINT *, "Index is ",mindex
!      PRINT *, "Index2 is ",index2
  !    PRINT *, "Wait_len is ",wait_len
      IF ( DEBUG .EQV. .TRUE. ) THEN
        PRINT *, 'In Meta_hop, original coords are:',i,j,k
        PRINT *, 'In Meta_hop, new coords are:',i_hop,j_hop,k_hop
      END IF
      IF ( index2 .GT. 0 ) THEN
        species2 = wait_list(index2)%sp_num
      ELSE IF ( index2 .LT. 0 ) THEN
        species2 = ABS(index2)
      END IF

!      PRINT *, "Species 1 is ",species,"and species 2 is ",species2
      ! Look at location, if it is empty, move there
      ! else, check whether the occupant of the new
      ! site is a reaction partner. If it is, react
      ! else, stay at original location and "end turn"
      IF ( index2 .EQ. 0 ) THEN ! Site is empty, move to site
!        PRINT *, "Hopping to new site"


        ! Remove reactant from old site and set to empty (i.e. 0)
        matrix_ptr(i,j,k) = 0

        ! Update coordinate info
        wait_list(mindex)%i = i_hop
        wait_list(mindex)%j = j_hop
        wait_list(mindex)%k = k_hop

        matrix_ptr(i_hop,j_hop,k_hop) = mindex

       ! Calculate a new waiting time
        CALL wait_calc( wait_list, mindex, en_list, time )

        ! Event 1
        result = 1
  !      WRITE(1015,*) mindex,"hops"
      ELSE ! Site is occupied, check whether the species can react
        IF ( react_ptr(species, species2, 1) .EQ. 0 ) THEN ! Null event, species does not move
  !        PRINT *, "No reaction possible"
          ! Re-add species to wait_list with new coords and waiting time
          CALL wait_calc ( wait_list, mindex, en_list, time )
          result = 0
  !        WRITE(1015,*) mindex,"no react"
        ELSE ! Have the two species react and place products
  !        PRINT *, "Reaction possible"
  !        PRINT *, "Right before reaction, mindex is ",mindex
  !        PRINT *, "Right before reaction, i, j, k are ",i,j,k
!          PRINT *, "The value of the matrix at that point is ",matrix_ptr(i,j,k),'which should be',mindex
!          PRINT *, 'Calling reaction in meta_hop'
          CALL reaction( o3_prod,o3_dest,react_ptr, en_list, matrix_ptr, wait_list, wait_len, &
                          time, i, j, k, i_hop, j_hop, k_hop  )
          result = 2
  !        WRITE(1015,*) mindex,"reacts"
        END IF
      END IF
  !    PRINT *, "Ending meta_hop"
  !    CLOSE(1015)
  END SUBROUTINE meta_hop

  SUBROUTINE reaction( o3_prod,o3_dest,qube, en_list, matrix,  wait_list, wait_len, time, i_re, j_re, &
                       k_re, i_re2, j_re2, k_re2, ion_coords, anion_list)
  !
  ! Purpose:
  !   The purpose of this subroutine is to take some set of target coordinates,
  !  representing the location in the matrix of a second reactant, and place
  !  the products. If there is only one product, it is placed at the second site.
  !  If there are two products, they are placed on the old and new sites. In the
  !  case where there are three, a third location is selected at random.
  !
  ! Documentation:
  !  DATE          PROGRAMMER           DESCRIPTION
  !  ========      ==========           ===========
  !  20150513      C. Shingledecker     Original code
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! REACTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    INTEGER                                           , POINTER           :: o3_prod,o3_dest
    INTEGER            , INTENT(IN)                                       :: i_re,j_re,k_re
    INTEGER            , INTENT(IN)                                       :: i_re2,j_re2,k_re2
    INTEGER(KIND=SHORT)                                                   :: null
    INTEGER                                                               :: r1,r2
    INTEGER                                                               :: case_num
    INTEGER                                                               :: index,index2
    INTEGER                                           , POINTER           :: wait_len
    INTEGER                         , DIMENSION(:,:,:), POINTER           :: matrix
    INTEGER                         , DIMENSION(3)                        :: third_coords
    INTEGER            , INTENT(OUT), DIMENSION(3)             , OPTIONAL :: ion_coords
    INTEGER                         , DIMENSION(:)    , POINTER, OPTIONAL :: anion_list
    INTEGER(KIND=SHORT)             , DIMENSION(:,:,:), POINTER           :: qube
    INTEGER                         , DIMENSION(3)                        :: prods
    REAL(KIND=DBL)                                    , POINTER           :: time
    REAL                            , DIMENSION(:)    , POINTER           :: en_list
    TYPE(wait_info)                 , DIMENSION(:)    , POINTER           :: wait_list
    INTEGER                                                               :: n
    REAL                                                                  :: rnum
    INTEGER                                                               :: i, j, k
    INTEGER                         , DIMENSION(3,3)                      :: prod_coords
    INTEGER                                                               :: p1,p2
    INTEGER                                                               :: dummy_protons
    dummy_protons = 1


    ! Find the species numbers
    ! NB: if the number in the matrix is negative, then the species number
    ! is just the absolute value, but if it is positive, then it corresponds
    ! to an entry in the wait_list, i.e. that the species is mobile, and its
    ! species number can be found in the wait_list derived type structure.

    r1 = 0
    r2 = 0

    index  = matrix(i_re,j_re,k_re)
    index2 = matrix(i_re2,j_re2,k_re2)

    IF ( index .LT. 0 ) THEN ! First reactant
      r1 = ABS(index)
    ELSE IF ( index .GT. 0 ) THEN
      r1 = wait_list(index)%sp_num
    END IF

    IF ( index2 .LT. 0 ) THEN ! Second reactant
      r2 = ABS(index2)
    ELSE IF ( index2 .GT. 0 ) THEN
      r2 = wait_list(index2)%sp_num
    END IF

    prods = qube(r1,r2,:)

    !****************************************************************************
    !****BRANCHING_RATIOS********************************************************
    !****************************************************************************
    !VERY TEMPORARY FIX TO SIMULATE BRANCHING RATIOS: FIX!!!
    !****************************************************************************
    !****************************************************************************
!    IF ( r1 .EQ. 4 .AND. r2 .EQ. 1 .OR. r1 .EQ. 1 .AND. r2 .EQ. 4 ) THEN
!      rnum = RAND()
!      CALL RANDOM_NUMBER(rnum)
!      IF ( rnum .LE. O_O2_BRANCHING ) THEN
!        prods = (/ 1, 4, 0 /)
!      END IF
    IF ( r1 .EQ. 2 .AND. r2 .EQ. 3 .OR. r1 .EQ. 3 .AND. r2 .EQ. 2 ) THEN
!      rnum = RAND()
      CALL RANDOM_NUMBER(rnum)
      IF ( rnum .LE. O2_ION_BRANCHING ) THEN
        prods = (/ 7, 4, 0 /)
      END IF
    ELSE IF ( r1 .EQ. 8 .AND. r2 .EQ. 6 .OR. r1 .EQ. 6 .AND. r2 .EQ. 8 ) THEN
!      rnum = RAND()
      CALL RANDOM_NUMBER(rnum)
      IF ( rnum .LE. O3_O_ION_BRANCHING ) THEN
        prods = (/ 1, 4, 4 /)
      END IF
    ELSE IF ( r1 .EQ. 9 .AND. r2 .EQ. 5 .OR. r1 .EQ. 5 .AND. r2 .EQ. 9 ) THEN
!      rnum = RAND()
      CALL RANDOM_NUMBER(rnum)
      IF ( rnum .LE. O3_O_ION_BRANCHING ) THEN
        prods = (/ 1, 4, 4 /)
      END IF
    ELSE IF ( r1 .EQ. 8 .AND. r2 .EQ. 3 .OR. r1 .EQ. 3 .AND. r2 .EQ. 8 ) THEN
!      rnum = RAND()
      CALL RANDOM_NUMBER(rnum)
      IF ( rnum .LE. O3_O2_ION_BRANCHING ) THEN
        prods = (/ 7, 4, 4 /)
      END IF
    ELSE IF ( r1 .EQ. 9 .AND. r2 .EQ. 2 .OR. r1 .EQ. 2 .AND. r2 .EQ. 9 ) THEN
!      rnum = RAND()
      CALL RANDOM_NUMBER(rnum)
      IF ( rnum .LE. O3_O2_ION_BRANCHING ) THEN
        prods = (/ 7, 4, 4 /)
      END IF
    END IF

    !****************************************************************************
    !****ANALYTICS***************************************************************
    !****************************************************************************
    IF ( r1 .EQ. 7 .OR. r2 .EQ. 7 ) o3_dest = o3_dest + 1
    IF ( ANY(prods .EQ. 7 ) ) o3_prod = o3_prod + 1
    !****************************************************************************

    !****************************************************************************
    !****ERROR CHECKING**********************************************************
    !****************************************************************************

    ! Test to see if two reacting coordinates are the same
    IF ( i_re .EQ. i_re2 .AND. &
         j_re .EQ. j_re2 .AND. &
         k_re .EQ. k_re2 ) THEN
!      PRINT *, 'In reaction: i_re = i_re2...'
      RETURN
    END IF

    ! Test to see if second reactant is 0
    IF ( r2 .EQ. 0 ) THEN
      OPEN(UNIT=1013,FILE="test_wrong_spaces.txt")
      DO k=1,SIZE(matrix,3)
        DO j=1,SIZE(matrix,2)
          DO i=1,SIZE(matrix,1)
            IF ( ABS(matrix(i,j,k)) .GT. wait_len ) THEN
!              WRITE(1013,*)  matrix(i,j,k)
            END IF
          END DO
        END DO
      END DO
      CLOSE(1013)
      PRINT *, "ERROR: reactant 2 is 0"
      CALL EXIT()
    END IF

    !Print debug info to file if second reactant is O3 and debug set to on
    IF ( DEBUG .EQV. .TRUE. ) THEN
      PRINT *, 'In reaction: the coordinates are:',i_re,j_re,k_re,'and ',i_re2,j_re2,k_re2
      PRINT *, "wait_len in reaction is ",wait_len
      PRINT *,     r1, r2, prods
      WRITE(777,*) r1,',',r2,',',prods(1),',',prods(2),',',prods(3), ',', wait_len

      IF ( r1 .EQ. 7 .OR. r2 .EQ. 7 ) THEN
        PRINT *, '\/ In reaction \/'
        PRINT *, r1,',',r2,',',prods(1),',',prods(2),',',prods(3)
      END IF

      IF ( ANY(prods .EQ. 7 ) ) THEN
        PRINT *, '\/ In reaction \/'
        PRINT *, r1,',',r2,',',prods(1),',',prods(2),',',prods(3)
      END IF
    END IF

    ! Test to see if products are 0
    IF ( prods(1) .EQ. 0 ) THEN
      PRINT *, 'Uh-oh, we have a problem in Reaction!'
      PRINT *, 'r1=',r1,'r2=',r2
      PRINT *, 'prods=',prods
      PRINT *, 'Index of r1= ',index
      PRINT *, 'Index of r2= ',index2
      PRINT *, 'Coords are ',i_re,j_re,k_re
      PRINT *, 'Wait_len =',wait_len
      OPEN(UNIT=1013,FILE="test_wait_list.txt")
      DO n=1,wait_len+2
        WRITE(1013,*) wait_list(n)
      END DO
      CLOSE(1013)

      OPEN(UNIT=1014,FILE="wait_list_flaw.txt")
      DO n=1,wait_len
        IF ( wait_list(n)%sp_num .NE. 20 ) THEN
          IF ( matrix(wait_list(n)%i,wait_list(n)%j,wait_list(n)%k) .NE. n ) THEN
            PRINT *, 'Wrong space at line',n
            PRINT *, wait_list(n)
            WRITE(1014,*) 'Matrix val=',matrix(wait_list(n)%i,wait_list(n)%j,wait_list(n)%k),", SP_NUM=" &
                         ,wait_list(n)%sp_num,',INDEX=',n,", WAIT_LEN=",wait_len,',I=',wait_list(n)%i, &
                          ',J=', wait_list(n)%j,', K=',wait_list(n)%k
            WRITE(1014,*) '******************************************************'
          END IF
        END IF
      END DO
      CLOSE(1014)
      PRINT *,'Quiting...'
      CALL EXIT()
    END IF

    IF (DEBUG .EQV. .TRUE. ) THEN
      PRINT *, 'Initially matrix r1=',matrix(i_re,j_re,k_re)
      PRINT *, 'Initially matrix r2=',matrix(i_re2,j_re2,k_re2)
    END IF


    !*****************************************************************************
    !****Determine Products Case**************************************************
    !*****************************************************************************
    CALL rtype_tree(r1,r2,prods,case_num)

    !*****************************************************************************
    !****Place Products***********************************************************
    !*****************************************************************************
    prod_coords = 0
    CALL place_product(i_re,j_re,k_re,i_re2,j_re2,k_re2,index,index2,r1,r2,&
                       prods,case_num,matrix, wait_list,wait_len,en_list,time, prod_coords)


    !*****************************************************************************
    !****Save ion coordinates, if one of the products is an ion
    !*****************************************************************************
    IF ( PRESENT(ion_coords) ) THEN
      ion_coords = 0
      DO n = 1,3
        IF (prod_coords(n,1) .NE. 0 ) THEN
          IF ( ANY(anion_list .EQ. -1*matrix(prod_coords(n,1),prod_coords(n,2),prod_coords(n,3)) ) ) THEN
            ion_coords(1) = prod_coords(n,1)
            ion_coords(2) = prod_coords(n,2)
            ion_coords(3) = prod_coords(n,3)
          END IF
        END IF
      END DO
    END IF

   IF ( DEBUG .EQV. .TRUE. ) THEN
     CALL counter(o3_prod,o3_dest,dummy_protons,time,AB_UNIT_NUM,matrix,wait_list,4,7,wait_len)
     PRINT *, 'Afterwards matrix r1=',matrix(i_re,j_re,k_re)
     IF ( matrix(i_re,j_re,k_re) .GT. 0 ) PRINT *, &
     'Wait_list at ^ is',wait_list(matrix(i_re,j_re,k_re))
     PRINT *, 'Afterwards matrix r2=',matrix(i_re2,j_re2,k_re2)
     IF (matrix(i_re2,j_re2,k_re2) .GT. 0 ) PRINT *, &
     'Wait_list at ^ is',wait_list(matrix(i_re2,j_re2,k_re2))
     DO n=1,3
        PRINT *,prod_coords(n,:)
     END DO
   END IF

  END SUBROUTINE reaction

  SUBROUTINE canreact ( reactant1, reactant2, qube, wait_list, null )
  ! --null .EQ. 0 => reaction can occur
  ! --null .EQ. 1 => no reaction can occur
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! CANREACT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    INTEGER                                                     :: species1
    INTEGER                                                     :: species2
    INTEGER            , INTENT(IN)                             :: reactant1
    INTEGER            , INTENT(IN)                             :: reactant2
    INTEGER(KIND=SHORT), INTENT(OUT)                            :: null
    INTEGER(KIND=SHORT)             , DIMENSION(:,:,:), POINTER :: qube
    TYPE(wait_info)                 , DIMENSION(:)    , POINTER :: wait_list

    species1 = 0
    species2 = 0


    ! Find the first species number
    IF ( reactant1 .LT. 0 ) THEN
      species1 = ABS(reactant1)
    ELSE IF ( reactant1 .GT. 0 ) THEN
      species1 = wait_list(reactant1)%sp_num
    END IF

    ! Find the second species number
    IF ( reactant2 .LT. 0 ) THEN
      species2 = ABS(reactant2)
    ELSE IF ( reactant2 .GT. 0 ) THEN
      species2 = wait_list(reactant2)%sp_num
    END IF

    ! Determine if species can react
    IF ( qube(species1,species2,1) .EQ. 0 ) THEN
      null = 1 ! No reaction
    ELSE
      null = 0 ! Reaction
    END IF
  END SUBROUTINE canreact

  SUBROUTINE reactant_remove( wait_list, index, matrix_rr, wait_len )
  !
  ! Purpose:
  !  The purpose of this subroutine is to take the
  !  coordinates of a mobile reactant and remove it
  !  from the waiting list of mobile species.
  !
  !  Documentation:
  !  DATE          PROGRAMMER           DESCRIPTION
  !  ========      ==========           ===========
  !  20150414      C. Shingledecker     Original code
  !  20150512      C. Shingledecker     Changed to handle TYPEs
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! REACTANT_REMOVE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    ! Data dictionary: variables passed to subroutine
    INTEGER                                          , POINTER :: wait_len !number of non-zero entries in wait_list
    INTEGER         , INTENT(IN)                               :: index !index of row to be removed
    INTEGER                        , DIMENSION(:,:,:), POINTER :: matrix_rr
    TYPE (wait_info)               , DIMENSION(:)    , POINTER :: wait_list

    OPEN(UNIT=1015,FILE="remove_row_list.txt",POSITION="append")
    IF ( index .NE. wait_len ) THEN
      ! Copy information in last entry to index
      wait_list(index)%wait_time    = wait_list(wait_len)%wait_time
      wait_list(index)%i            = wait_list(wait_len)%i
      wait_list(index)%j            = wait_list(wait_len)%j
      wait_list(index)%k            = wait_list(wait_len)%k
      wait_list(index)%sp_num       = wait_list(wait_len)%sp_num
      wait_list(index)%act_type     = wait_list(wait_len)%act_type

      ! Update matrix
      matrix_rr(wait_list(index)%i,wait_list(index)%j,wait_list(index)%k) = index

      ! Make info in last entry equal to 0
      wait_list(wait_len)%wait_time = 0
      wait_list(wait_len)%i         = 0
      wait_list(wait_len)%j         = 0
      wait_list(wait_len)%k         = 0
      wait_list(wait_len)%sp_num    = 0
      wait_list(wait_len)%act_type  = 0

      ! Update number of non-zero species by -1
      wait_len = wait_len - 1
      WRITE(1015,*) 'Index=',index,"wait_len=",wait_len,"matrix=",&
                     matrix_rr(wait_list(index)%i,wait_list(index)%j,wait_list(index)%k), &
                    'coords=', wait_list(index)%i,wait_list(index)%j,wait_list(index)%k
    ELSE IF ( index .EQ. wait_len ) THEN
      ! Make info in last entry equal to 0
      WRITE(1015,*) 'Index=',index,"wait_len=",wait_len,"matrix=", &
                     matrix_rr(wait_list(index)%i,wait_list(index)%j,wait_list(index)%k), &
                    'coords=', wait_list(index)%i,wait_list(index)%j,wait_list(index)%k

      wait_list(wait_len)%wait_time = 0
      wait_list(wait_len)%i         = 0
      wait_list(wait_len)%j         = 0
      wait_list(wait_len)%k         = 0
      wait_list(wait_len)%sp_num    = 0
      wait_list(wait_len)%act_type  = 0

      wait_len = wait_len - 1
    END IF
    CLOSE(1015)

  END SUBROUTINE reactant_remove

  SUBROUTINE wait_calc ( wait_list, index, en_list, time )
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
    INTEGER       , INTENT(IN)                           :: index     !index of species
    REAL(KIND=DBL)                             , POINTER :: time      !total simulation time
    REAL(KIND=DBL)                                       :: rand_num  !pseudorandom number
    REAL(KIND=DBL)                                       :: b_1       !surface thermal hopping rate
    REAL(KIND=DBL)                                       :: b_2       !surface desorption rate
    REAL(KIND=DBL)                                       :: b_3       !bulk diffusion rate
    REAL(KIND=DBL)                                       :: b         !total rate, from CH14
    REAL                       , DIMENSION(:)  , POINTER :: en_list    !binding/diffusion energy list
    TYPE (wait_info)           , DIMENSION(:)  , POINTER :: wait_list !list of mobile species


!    IF ( ANY(FAST_REACTS .EQ. wait_list(index)%sp_num) ) THEN
!      wait_list(index)%wait_time = 1.0D-14 + time
!    ELSE
      IF ( wait_list(index)%i .EQ. 1 ) THEN
        ! Surface species, separate rates for
        ! desorption and diffusion
        b_1 = trl_nu*EXP( - ( ( en_list(wait_list(index)%sp_num)*E_SURF) / kin_temp ) )
        b_2 = trl_nu*EXP( - ( en_list(wait_list(index)%sp_num)           / kin_temp ) )
        b = b_1 + b_2
      ELSE
        ! Bulk species, only bulk diffusion
        b_3 = trl_nu*EXP( - ( en_list(wait_list(index)%sp_num)*E_BULK    / kin_temp ) )
        b = b_3
      END IF
      CALL RANDOM_NUMBER(rand_num)
      ! Calculate waiting time
      wait_list(index)%wait_time = (-LOG(rand_num) / b) + time
!    END IF

    ! Assign action type for next move
    CALL action_figure(wait_list,index,rand_num,en_list)
  END SUBROUTINE wait_calc

  SUBROUTINE minmod ( wait_list, mindex, wait_len )
  !
  ! Purpose:
  !  The purpose of this subroutine is to find the
  ! minimum waiting time in the wait_list structure
  !
  ! Note:
  !   This subroutine is required since MINLOC will
  !  just return a zero value AND go through the
  !  entire structure, which is not desired.
  !
  !  Documentation:
  !  DATE          PROGRAMMER           DESCRIPTION
  !  ========      ==========           ===========
  !  20150512      C. Shingledecker     Original code
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! MINMOD !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    ! Data dictionary: variables passed to subroutine
    INTEGER                                              :: i !counter
    INTEGER                                    , POINTER :: wait_len !number of non-zero entries in wait_list
    INTEGER         , INTENT(OUT)                        :: mindex !index of row with minimum time
    REAL(KIND=DBL)                                       :: time_temp
    TYPE (wait_info)             , DIMENSION(:), POINTER :: wait_list

    ! find minimum time
    mindex = 1
    time_temp = wait_list(1)%wait_time
    DO i=2,wait_len
      IF (wait_list(i)%wait_time .LT. time_temp ) THEN
        mindex = i
        time_temp = wait_list(i)%wait_time
      END IF
    END DO
  END SUBROUTINE minmod

  SUBROUTINE counter(o3_prod,o3_dest,numprotons,time, unit_num, matrix,wait_list,sp1,sp2,wait_len)
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
    INTEGER            , INTENT(IN)                            :: unit_num
    INTEGER                                                    :: i,j,k,nn
    INTEGER                                                    :: sp1_count, sp2_count
    INTEGER            , INTENT(IN)                            :: sp1, sp2
    INTEGER                        , DIMENSION(3)              :: dimens
    INTEGER                        , DIMENSION(:,:,:), POINTER :: matrix
    INTEGER                              , POINTER :: wait_len
    REAL(KIND=DBL)                                   , POINTER :: time
    REAL(KIND=DBL)                                             :: volume
    REAL(KIND=DBL)                                             :: denom
    REAL(KIND=DBL)                                             :: area
    REAL(KIND=DBL)                                             :: fluence
    TYPE(wait_info)                , DIMENSION(:)    , POINTER :: wait_list
    CHARACTER(len=80)                                          :: varfmt
    INTEGER                                                    :: wrong_count


!    volume = THICK*EDGE*EDGE
    area   = EDGE*EDGE
    denom = THICK*EDGE*EDGE*1.0E20
    sp1_count = 0
    sp2_count = 0
    wrong_count = 0
    ! Method 1 of fluence calculation
    fluence  = CR_FLUX*time ! Note: This is the x-value for the objective function
    ! Method 2 of fluence calculation (only use 1 at a time )
    ! fluence = numprotons/area


    dimens(1) = SIZE(matrix,1)
    dimens(2) = SIZE(matrix,2)
    dimens(3) = SIZE(matrix,3)


    IF ( TEST_WRONG .EQV. .TRUE. ) OPEN(UNIT=1013,FILE="counter_test_wrong_spaces.txt",STATUS='REPLACE')
    DO k = 1,dimens(3)
      DO j = 1,dimens(2)
        DO i = 1,dimens(1)
          IF ( matrix(i,j,k) .NE. 0 ) THEN
            IF ( matrix(i,j,k) .LT. 0 ) THEN
              IF ( ABS(matrix(i,j,k)) .EQ. sp2 ) THEN
                sp2_count = sp2_count + 1
              END IF
            ELSE IF ( matrix(i,j,k) .GT. 0 ) THEN
              IF ( wait_list(matrix(i,j,k))%sp_num .EQ. sp1 ) THEN
                sp1_count = sp1_count + 1
              ELSE IF ( wait_list(matrix(i,j,k))%sp_num .EQ. sp2 ) THEN
                sp2_count = sp2_count + 1
              END IF

              IF ( TEST_WRONG .EQV. .TRUE. ) THEN
                IF( ANY(MOBILE_LIST .EQ. wait_list(matrix(i,j,k))%sp_num ) ) THEN
                  CONTINUE
                ELSE
                  wrong_count = wrong_count + 1
                  WRITE(1013,*)  'Matrix=',matrix(i,j,k),'wait_list=',wait_list(matrix(i,j,k))
                  PRINT *, 'We have a wrong space:',matrix(i,j,k),' at',i,j,k
                END IF
              END IF
            END IF
          END IF
        END DO
      END DO
    END DO
!    END IF

    IF ( TEST_WRONG .EQV. .TRUE. ) THEN
      OPEN(UNIT=1014,FILE="counter_test_wait_list.txt",STATUS='REPLACE')
      DO nn = 1,wait_len
        IF ( wait_list(nn)%sp_num .NE. 20 ) THEN
          IF( matrix(wait_list(nn)%i,wait_list(nn)%j,wait_list(nn)%k) .NE. nn ) THEN
            wrong_count = wrong_count + 1
            WRITE(1014,*)  'Matrix=',matrix(wait_list(nn)%i,wait_list(nn)%j,wait_list(nn)%k),'wait_list=',wait_list(nn)
            PRINT *, 'We have a list flaw at:',wait_list(nn)%i,wait_list(nn)%j,wait_list(nn)%k
          END IF
        END IF
      END DO
    END IF

    IF ( TEST_WRONG .EQV. .TRUE. ) THEN
      CLOSE(1013)
      CLOSE(1014)
      IF ( wrong_count .GT. 0 ) THEN
        PRINT *, 'Wrong spaces found in Counter!! Exiting!!'
        CALL EXIT()
      END IF
    END IF

 !   sp2_count = REAL(o3_prod - o3_dest)
    IF ( NO_OUTPUT .EQV. .FALSE. ) THEN
      WRITE(unit_num,*) time,',', fluence,',',sp1_count,',',sp2_count,',',numprotons,',' &
                        ,(REAL(o3_prod)/REAL(o3_dest))
    END IF

    IF ( QUIET .EQV. .FALSE. ) THEN
      varfmt = "(A6,ES10.4,A9,ES10.4)"
      PRINT varfmt, " TIME=",time,"FLUENCE=",fluence
      varfmt = "(A5,ES10.4,A6,ES10.4)"
      PRINT varfmt, " [O]=",sp1_count/denom," [O3]=",sp2_count/denom
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

  SUBROUTINE base_ionization( o3_prod,o3_dest,ev_coords,react_cube, matrix,  en_list, ionlist, &
                              wait_list, wait_len, time, ev_nums,null,elec_out )
    ! Purpose:
    !   To calculate the track of a particle of ionizing radiation through a
    !  crystaline solid.
    !
    ! Note:
    !   The input array, sigmas, contains the proton collision cross-sections.
    !  The contents of the array are:
    !
    ! Note:
    !   The array ev_nums contains the following values for pseudo-reactants:
    !
    ! -- ev_nums(1) = cosmic-ray species number
    ! -- ev_nums(2) = excitation species number
    ! -- ev_nums(3) = electron species number
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! BASE_IONIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    !******************!
    ! Input and output !
    !******************!
    INTEGER                              , POINTER :: o3_prod,o3_dest
    INTEGER                              , POINTER :: wait_len
    INTEGER            , DIMENSION(:)    , POINTER :: ionlist !list of anionic species
    INTEGER            , DIMENSION(:,:,:), POINTER :: matrix !ice-mantle matrix
    INTEGER            , DIMENSION(3)              :: ev_coords !coordiantes of collision
    INTEGER(KIND=SHORT), DIMENSION(:,:,:), POINTER :: react_cube !array of products/reactions
    REAL               , DIMENSION(:)    , POINTER :: en_list !list of binding and desorption energies
    REAL(KIND=DBL)                       , POINTER :: time
    TYPE(wait_info)    , DIMENSION(:)    , POINTER :: wait_list
    INTEGER            , DIMENSION(3)              :: ev_nums !values of special pseudoreactants
    INTEGER            , INTENT(OUT)   , DIMENSION(3), OPTIONAL :: elec_out

    !*****************!
    ! Local variables !
    !*****************!
    INTEGER(KIND=SHORT), INTENT(OUT)               :: null !null flag
    INTEGER                                        :: breakout
    INTEGER                                        :: large_count
    INTEGER                                        :: small_count
    INTEGER                                        :: switch !variable that determines the nature of collisions
    INTEGER, DIMENSION(3)                          :: elec_coords
    INTEGER, DIMENSION(3)                          :: ion_coords
    INTEGER, DIMENSION(3)                          :: coords
    INTEGER, DIMENSION(6,4)                        :: large_temp
    INTEGER, DIMENSION(4,4)                        :: small_temp



  !  PRINT *, 'Now in base_ionization'
    ! Switch = 2 => ionization
    switch = 2

    ! Call Cern to generate electron
    ! NB: the electron should be the second product in the "prods" array
    ! in Cern.
  !  PRINT *, 'At the beginning of base_ionization, ev_coords are:',ev_coords
    null = 0
    CALL cern( o3_prod,o3_dest,null,en_list, react_cube, matrix, ev_nums, ev_coords, switch, &
               wait_list, wait_len, time, elec_coords )
  !  PRINT *, 'After cern in base_ionization, elec_coords=',elec_coords
    IF ( PRESENT(elec_out) ) THEN
      elec_out = elec_coords
    END IF

    IF ( null .EQ. 1 ) RETURN

    ! Find a potential reaction partners for electron
    ! NB: Pass elec_coords to lookaroundyou
    IF ( DEBUG .EQV. .TRUE. ) THEN
      PRINT *, 'Calling lookaroundyou'
      PRINT *, 'In base_ionization, elec_coords=',elec_coords
    END IF
    CALL lookaroundyou( react_cube, matrix, elec_coords, null, small_count, &
                        large_count, small_temp, large_temp, wait_list )

    IF ( DEBUG .EQV. .TRUE. ) PRINT *, "Calling reaction for electron"
    IF ( null .EQ. 1 ) THEN ! No other reactants: have electron and cation react
      IF ( DEBUG .EQV. .TRUE. ) PRINT *, "Making electron and initial ion react"
      CALL reaction( o3_prod,o3_dest,react_cube, en_list, matrix, wait_list, wait_len, time, &
                     ev_coords(1), ev_coords(2), ev_coords(3), elec_coords(1), &
                     elec_coords(2), elec_coords(3) )
    ELSE
     ! Choose one at random
      IF ( DEBUG .EQV. .TRUE. ) PRINT *, "Calling solarlottery in base_ionization"
      breakout = 0
      DO
        IF ( DEBUG .EQV. .TRUE. ) PRINT *, 'Actually calling solarlottery in base_ionization'
        CALL solarlottery( small_count,large_count,small_temp,large_temp,coords )
        IF ( matrix(coords(1),coords(2),coords(3)) .NE.  matrix(ev_coords(1),ev_coords(2),ev_coords(3)) ) EXIT
        breakout = breakout + 1
        IF ( DEBUG .EQV. .TRUE. ) PRINT *, 'breakout is:',breakout
        IF ( breakout .GE. 10 ) THEN
          IF ( DEBUG .EQV. .TRUE. ) PRINT *, 'Calling reaction in base_ionization due to breakout'
          CALL reaction( o3_prod,o3_dest,react_cube, en_list, matrix, wait_list, wait_len, time, &
                         ev_coords(1), ev_coords(2), ev_coords(3), elec_coords(1), &
                         elec_coords(2), elec_coords(3) )
          RETURN
        END IF
      END DO
      ! Make the electron that has just formed react
      ! NB: Pass elec_coords and coords chosen by solarlottery
      IF ( DEBUG .EQV. .TRUE. ) PRINT *, "Calling reaction in BI to make electron create anion"
      CALL reaction( o3_prod,o3_dest,react_cube, en_list, matrix, wait_list, wait_len, time, &
                     elec_coords(1), elec_coords(2), elec_coords(3), &
                     coords(1), coords(2), coords(3), &
                     ion_coords, ionlist )
      ! Recombine the ions
      IF ( DEBUG .EQV. .TRUE. ) THEN
        PRINT *, "Recombining ions"
        PRINT *, "ion_coords are: ",ion_coords
        PRINT *, "At ion_coords, matrix=",matrix(ion_coords(1),ion_coords(2),ion_coords(3))
        PRINT *, "ev_coords are: ",ev_coords
        PRINT *, "At ev_coords, matrix=",matrix(ev_coords(1),ev_coords(2),ev_coords(3))
        PRINT *, 'Calling reaction in base_ionization'
      END IF
      CALL reaction( o3_prod,o3_dest,react_cube, en_list, matrix, wait_list, wait_len, time, &
                     ev_coords(1), ev_coords(2), ev_coords(3), &
                     ion_coords(1), ion_coords(2), ion_coords(3) )
      IF ( DEBUG .EQV. .TRUE. ) THEN
        PRINT *, "!!!IONS RECOMBINED!!!"
        PRINT *, "NOW:ion_coords are: ",ion_coords
        PRINT *, "NOW:At ion_coords, matrix=",matrix(ion_coords(1),ion_coords(2),ion_coords(3))
        PRINT *, "NOW:ev_coords are: ",ev_coords
        PRINT *, "NOW:At ev_coords, matrix=",matrix(ev_coords(1),ev_coords(2),ev_coords(3))
      END IF
    END IF

    IF ( DEBUG .EQV. .TRUE. ) PRINT *, '****Ending Base_Ionization****'

  END SUBROUTINE base_ionization

  SUBROUTINE make_react(o3_prod,o3_dest,qube, en_list, matrix, wait_list, wait_len, time,spec_index )
    ! Purpose:
    !   The purpose of this subroutine is to make certain species react as soon
    !  as they are formed. This would happen in the case of a species formed in
    !  the bulk that reacts rapidly with the matrix species, as in the cas of:
    !
    !    O + O2 -> O3 + energy(absorbed by the solid)
    !
    !  where here, O2 makes up the solid and O is the newly formed species.
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! MAKE_REACT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    INTEGER                                          , POINTER :: o3_prod,o3_dest
    INTEGER                                                    :: large_count
    INTEGER                                                    :: small_count
    INTEGER                                                    :: breakout
    INTEGER(KIND=SHORT)                                        :: null
    INTEGER                        , DIMENSION(3)              :: ev_coords !location of fast-reacting species
    INTEGER                        , DIMENSION(3)              :: coords
    INTEGER                        , DIMENSION(6,4)            :: large_temp
    INTEGER                        , DIMENSION(4,4)            :: small_temp
    INTEGER                                          , POINTER :: wait_len
    INTEGER                        , DIMENSION(:,:,:), POINTER :: matrix
    INTEGER(KIND=SHORT)            , DIMENSION(:,:,:), POINTER :: qube !array of products/reactions
    REAL                           , DIMENSION(:)    , POINTER :: en_list !list of binding and desorption energies
    REAL(KIND=DBL)                                   , POINTER :: time
    TYPE(wait_info)                , DIMENSION(:)    , POINTER :: wait_list !list of mobile species
    INTEGER                                                    :: spec_index

    ev_coords(1) = wait_list(spec_index)%i
    ev_coords(2) = wait_list(spec_index)%j
    ev_coords(3) = wait_list(spec_index)%k


    ! Check to see if there are any potential second reactants
!    PRINT *, 'Calling lookaroundyou in make_react'
    CALL lookaroundyou( qube, matrix, ev_coords, null, small_count, &
                        large_count, small_temp, large_temp, wait_list )
    IF ( null .EQ. 1 ) THEN
      wait_list(spec_index)%wait_time = 10.00
      RETURN ! If there aren't, keep fast-reacting species to the wait_list
    ELSE IF ( null .EQ. 0 ) THEN
      DO ! If there are, choose one at random
        CALL solarlottery( small_count,large_count,small_temp,large_temp,coords )
        IF ( matrix(coords(1),coords(2),coords(3)) .NE.  matrix(ev_coords(1),ev_coords(2),ev_coords(3)) ) EXIT
        breakout = breakout + 1
        IF ( breakout .GE. 10 ) RETURN
      END DO
      ! Make the species react
!      PRINT *, 'Calling reaction in make_react'
      CALL reaction( o3_prod,o3_dest,qube, en_list, matrix, wait_list, wait_len, time, &
                     ev_coords(1), ev_coords(2), ev_coords(3), &
                     coords(1), coords(2), coords(3))
    END IF
  END SUBROUTINE make_react

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
!        PRINT *, '#',n,'For energy:',energy
!        PRINT *, 'i=',i
!        PRINT *, 'k=',k
!        PRINT *, 'kb=',kb
!        PRINT *, 'j=',j
!        PRINT *, 'jb=',jb
!        PRINT *, 'jc=',jc
!        PRINT *, 'gs=',gs
!        PRINT *, 'gb=',gb
!        PRINT *, 'ts=',ts
!        PRINT *, 'ta=',ta
!        PRINT *, 'tb=',tb
        IF ( en .LT. i ) THEN
          sig = 0D0
        ELSE
          !i. Calculate the A(E) value from Green & Sawada
          ae = a_gs(en,k,kb,j,jb,jc)
!          PRINT *, '#',n,' ae=',ae
          !ii. Calculate the \Gamma(E) factor
          ge = gamma_gs(en,gs,gb)
!          PRINT *, '#',n,' ge=',ge
          !iii. Calculate the T_0 value
          tnaught = t_0_gs(en,ta,tb,ts)
!          PRINT *, '#',n,' tnaught=',tnaught
          !iv. Calculate the Tmas value
          tmax = t_max_gs(en,i)
!          PRINT *, '#',n,' tmax=',tmax
          !v. Calculate the cross-section for the state
          sig = green_sawada(ae,ge,tmax,tnaught)
        END IF
!        PRINT *, 'the',n,' value of sig is:',sig
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
      END ASSOCIATE
    END DO
    se_box%se_alwd_extot = SUM(se_box%se_alwdsigs)

    !(4) Calculate forbidden excitation cross-sections
    !NB: As above, no loop is required, since the subroutine
    !    returns an array of values
    DO n=1,SIZE(o2_e_ex_fbdn)
      ASSOCIATE( e => se_box%se_energy            , &
                 w => se_box%se_fbdn(n)%wj_fbdn   , &
                 f => se_box%se_fbdn(n)%fj_fbdn   , &
                 o => se_box%se_fbdn(n)%omega_fbdn, &
                 a => se_box%se_fbdn(n)%alpha_fbdn, &
                 b => se_box%se_fbdn(n)%beta_fbdn    )
        se_box%se_fbdnsigs(n) = greendutta(e,f,w,o,a,b)
      END ASSOCIATE
    END DO
    se_box%se_fbdn_extot = SUM(se_box%se_fbdnsigs)

    !(5) The total electron impact excitation is the sum of the
    !    allowed and forbidden transition cross-sections
    se_box%se_extot = se_box%se_alwd_extot + se_box%se_fbdn_extot

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
    CALL RANDOM_NUMBER(rn)
    prob = (se_box%se_alwd_extot/se_box%se_ineltot)
    IF ( rn .GT. prob ) THEN
      ALLOCATE(arr(SIZE(se_box%se_fbdnsigs),2))
      arr = 0
      arr(:,1) = se_box%se_fbdnsigs
      arr(:,2) = se_box%se_fbdn%wj_fbdn
      sigtot   = se_box%se_fbdn_extot
    ELSE
      ALLOCATE(arr(SIZE(se_box%se_alwdsigs),2))
      arr = 0
      arr(:,1) = se_box%se_alwdsigs
      arr(:,2) = se_box%se_alwd%wj_alwd
      sigtot   = se_box%se_alwd_extot
    END IF

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
    se_box%se_fbdn_extot = 0

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

      !Initialize forbidden excitation arrays
      ALLOCATE(se_box%se_fbdn(SIZE(o2_e_ex_fbdn)))
      se_box%se_fbdn = o2_e_ex_fbdn
      ALLOCATE(se_box%se_fbdnsigs(SIZE(o2_e_ex_fbdn)))
      se_box%se_fbdnsigs = 0
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
      DEALLOCATE(se_box%se_fbdn)
      DEALLOCATE(se_box%se_fbdnsigs)
      RETURN
    ELSE
      RETURN
    END IF
  END SUBROUTINE se_info_garbage

  SUBROUTINE fast_reaction(o3_prod,o3_dest,spec_index,wait_len,matrix,react_cube,en_list,time,wait_list)
  !
  ! Purpose:
  !   This subroutine is designed for those species that react quickly with the
  !   surrounding matrix. Such species are specified in the parameters.f03 file
  !   and are given a special react type when they are formed.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! FAST_REACTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE
    !******************!
    ! Input and output !
    !******************!
    INTEGER                              , POINTER :: o3_prod,o3_dest
    INTEGER                                        :: spec_index !index of fast-reacting species
    INTEGER                              , POINTER :: wait_len
    INTEGER            , DIMENSION(:,:,:), POINTER :: matrix !ice-mantle matrix
    INTEGER(KIND=SHORT), DIMENSION(:,:,:), POINTER :: react_cube !array of products/reactions
    REAL               , DIMENSION(:)    , POINTER :: en_list !list of binding and desorption energies
    REAL(KIND=DBL)                       , POINTER :: time
    TYPE(wait_info)    , DIMENSION(:)    , POINTER :: wait_list

    !*****************!
    ! Local variables !
    !*****************!
    INTEGER(KIND=SHORT)                            :: null_1,null_2 !null flag
    INTEGER                                        :: large_count
    INTEGER                                        :: small_count
    INTEGER, DIMENSION(3)                          :: coords
    INTEGER, DIMENSION(3)                          :: temp_coords
    INTEGER, DIMENSION(6,4)                        :: large_temp
    INTEGER, DIMENSION(4,4)                        :: small_temp
    INTEGER                                        :: result
    REAL(KIND=DBL)                                 :: rand_num,b

    IF ( DEBUG .EQV. .TRUE. ) PRINT *, '*****CALLING FAST REACTION*****'

    !Step i. Initialize temporary variables
    null_1      = 0
    null_2      = 0
    large_count = 0
    small_count = 0
    coords      = 0
    temp_coords = 0
    large_temp  = 0
    small_temp  = 0

    !Step ii. Get coords from struct
    coords(1) = wait_list(spec_index)%i
    coords(2) = wait_list(spec_index)%j
    coords(3) = wait_list(spec_index)%k

    !Step 1. Look around to see if there are any potential co-reactants
    CALL lookaroundyou( react_cube, matrix, coords, null_1, small_count, &
                        large_count, small_temp, large_temp, wait_list )

    SELECT CASE (null_1)
    CASE(0)
    IF ( DEBUG .EQV. .TRUE. ) PRINT *, '*****CASE (0) IN FAST REACTION*****'
    !Step 2. If there are potential co-reactants, randomly choose one and react.
    CALL solarlottery( small_count,large_count,small_temp,large_temp,temp_coords )
      IF ( coords(1) .EQ. temp_coords(1) .AND. &
           coords(2) .EQ. temp_coords(2) .AND. &
           coords(3) .EQ. temp_coords(3) )       THEN
!        PRINT *, 'coords =',coords,'= temp_coords= ', temp_coords
        b = TRL_NU*EXP( - ( en_list(wait_list(spec_index)%sp_num)*E_BULK    / KIN_TEMP ) )
        CALL RANDOM_NUMBER(rand_num)
        wait_list(spec_index)%wait_time = (-LOG(rand_num) / b) + time
        wait_list(spec_index)%act_type = 1
        RETURN
      ELSE
        IF ( matrix(temp_coords(1),temp_coords(2),temp_coords(3)) .NE. 0  ) THEN
          CALL reaction( o3_prod,o3_dest,react_cube, en_list, matrix, wait_list, wait_len, time, &
                         coords(1), coords(2), coords(3), &
                         temp_coords(1),temp_coords(2),temp_coords(3) )
        ELSE
        PRINT *, 'ERROR: matrix=0'
        CALL wait_calc(wait_list,spec_index,en_list,time)
        wait_list(spec_index)%act_type = 1
        END IF
      END IF
    CASE(1)
    !Step 3. If there are no co-reactants, look for empty species to hop to
      b = TRL_NU*EXP( - ( en_list(wait_list(spec_index)%sp_num)*E_BULK    / KIN_TEMP ) )
      CALL RANDOM_NUMBER(rand_num)
      wait_list(spec_index)%wait_time = (-LOG(rand_num) / b) + time
      wait_list(spec_index)%act_type = 1
      IF ( DEBUG .EQV. .TRUE. ) PRINT *, '*****CASE (1) IN FAST REACTION*****'
!      CALL meta_hop( o3_prod,o3_dest,spec_index,react_cube,matrix,en_list,wait_list,result,wait_len,time )
    END SELECT
    IF ( DEBUG .EQV. .TRUE. ) PRINT *, '*****ENDING FAST REACTION*****'
  END SUBROUTINE fast_reaction

  SUBROUTINE find_cr ( wait_list, cr_index, wait_len, cr_num, en_list, time )
  !
  ! Purpose:
  !  The purpose of this subroutine is to find the
  ! next cosmic-ray time in the wait list
  !
  !  Documentation:
  !  DATE          PROGRAMMER           DESCRIPTION
  !  ========      ==========           ===========
  !  20160209      C. Shingledecker     Original code
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! FIND_CR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    ! Data dictionary: variables passed to subroutine
    INTEGER                                              :: i !counter
    INTEGER                                    , POINTER :: wait_len !number of non-zero entries in wait_list
    INTEGER         , INTENT(OUT)                        :: cr_index !index of row with minimum time
    INTEGER                                              :: cr_num
    TYPE (wait_info)             , DIMENSION(:), POINTER :: wait_list
    REAL                         , DIMENSION(:), POINTER :: en_list !list of binding and desorption energies
    REAL(KIND=DBL)                             , POINTER :: time
    INTEGER                                              :: n

    ! find minimum time
    cr_index = 1
    DO i=2,wait_len
      IF (wait_list(i)%sp_num .EQ. cr_num ) THEN
        cr_index = i
        RETURN
      END IF
    END DO

    time = wait_list(cr_index)%wait_time

    DO n=1,wait_len
        IF ( n .NE. cr_index ) THEN
            CALL wait_calc(wait_list,n,en_list,time)
        END IF
    END DO
  END SUBROUTINE find_cr

  SUBROUTINE rtype_tree(r1,r2,prods,rtype)
  !
  ! Purpose:
  !  The purpose of this subroutine is to find the
  ! type of the reaction based on the reactants and
  ! products.
  !
  !  Documentation:
  !  DATE          PROGRAMMER           DESCRIPTION
  !  ========      ==========           ===========
  !  20160215      C. Shingledecker     Original code
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! RTYPE_TREE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: r1,r2 ! Reactants 1 & 2
    INTEGER, INTENT(IN), DIMENSION(3) :: prods ! Array of products
    INTEGER, INTENT(OUT) :: rtype ! The type of the reaction

    LOGICAL :: r1_sp
    LOGICAL :: r1_m
    LOGICAL :: r2_m
    LOGICAL :: p1_m
    LOGICAL :: p2_m
    LOGICAL :: p2_0

    ! Initialize logical variables to false
    r1_sp = .FALSE.
    r1_m  = .FALSE.
    r2_m  = .FALSE.
    p1_m  = .FALSE.
    p2_m  = .FALSE.
    p2_0  = .FALSE.

    ! Change logical variables, as appropriate
    IF ( ANY(SPECIAL_LIST .EQ. r1 ) )        r1_sp = .TRUE.
    IF ( ANY(MOBILE_LIST .EQ. r1 ) )         r1_m  = .TRUE. ! Check is r1 is mobile
    IF ( ANY(MOBILE_LIST .EQ. r2 ) )         r2_m  = .TRUE. ! Check if r2 is mobile
    IF ( ANY(MOBILE_LIST .EQ. prods(1) ) )   p1_m  = .TRUE. ! Check is p1 is mobile
    IF ( prods(2) .EQ. 0 ) THEN                            ! Check is p2 is present
      p2_0 = .TRUE.  ! No p2
    ELSE
      p2_0 = .FALSE. ! p2 exists
      IF ( ANY(MOBILE_LIST .EQ. prods(2) ) ) p2_m  = .TRUE. ! Check is p2 is mobile
    END IF

    ! Initialize type
    rtype = 0

    ! Assign type
    IF ( r1_sp .EQV. .FALSE.) THEN
      IF ( r1_m .EQV. .TRUE. ) THEN         ! r1 = m
        IF ( r2_m .EQV. .TRUE. ) THEN       ! r1 = m, r2 = m
          IF ( p1_m .EQV. .TRUE. ) THEN     ! r1 = m, r2 = m, p1 = m
            IF ( p2_0 .EQV. .TRUE. ) THEN   ! r1 = m, r2 = m, p1 = m, p2 = 0 (1)
              rtype = 1
            ELSE
              IF ( p2_m .EQV. .TRUE. ) THEN ! r1 = m, r2 = m, p1 = m, p2 = m (9)
                rtype = 9
              ELSE                          ! r1 = m, r2 = m, p1 = m, p2 = i (10)
                rtype = 10
              END IF
            END IF
          ELSE                              ! r1 = m, r2 = m, p1 = i
            IF ( p2_0 .EQV. .TRUE. ) THEN   ! r1 = m ,r2 = m, p1 = i,  p2 = 0 (2)
              rtype = 2
            ELSE
              IF ( p2_m .EQV. .TRUE. ) THEN
                rtype = 11                  ! r1 = m, r2 = m, p1 = i, p2 = m (11)
              ELSE
                rtype = 12                  ! r1 = m, r2 = m, p1 = i, p2 = i (12)
              END IF
            END IF
          END IF
        ELSE                                ! r1 = m, r2 = i
          IF ( p1_m .EQV. .TRUE. ) THEN     ! r1 = m ,r2 = i, p1 = m
            IF ( p2_0 .EQV. .TRUE. ) THEN   ! r1 = m, r2 = i, p1 = m, p2 = 0 (3)
              rtype = 3
            ELSE
              IF ( p2_m .EQV. .TRUE. ) THEN ! r1 = m, r2 = i, p1 = m, p2 = m (13)
                rtype = 13
              ELSE                          ! r1 = m, r2 = i, p1 = m, p2 = i (14)
                rtype = 14
              END IF
            END IF
          ELSE                              ! r1 = m, r2 = i, p1 = i
            IF ( p2_0 .EQV. .TRUE. ) THEN   ! r1 = m, r2 = i, p1 = i, p2 = 0 (4)
              rtype = 4
            ELSE
              IF ( p2_m .EQV. .TRUE. ) THEN ! r1 = m, r2 = i, p1 = i, p2 = m (15)
                rtype = 15
              ELSE                          ! r1 = m, r2 = i, p1 = i, p2 = i (16)
                rtype = 16
              END IF
            END IF
          END IF
        END IF
      ELSE                                  ! r1 = i
        IF ( r2_m .EQV. .TRUE. ) THEN       ! r1 = i, r2 = m
          IF ( p1_m .EQV. .TRUE. ) THEN     ! r1 = i, r2 = m, p1 = m
            IF ( p2_0 .EQV. .TRUE. ) THEN   ! r1 = i, r2 = m, p1 = m, p2 = 0 (5)
              rtype = 5
            ELSE
              IF ( p2_m .EQV. .TRUE. ) THEN ! r1 = i, r2 = m, p1 = m, p2 = m (17)
                rtype = 17
              ELSE                          ! r1 = i, r2 = m, p1 = m, p2 = i (18)
                rtype = 18
              END IF
            END IF
          ELSE                              ! r1 = i, r2 = m, p1 = i
            IF ( p2_0 .EQV. .TRUE. ) THEN   ! r1 = i, r2 = m, p1 = i, p2 = 0 (6)
              rtype = 6
            ELSE
              IF ( p2_m .EQV. .TRUE. ) THEN ! r1 = i, r2 = m, p1 = i, p2 = m (19)
                rtype = 19
              ELSE                          ! r1 = i, r2 = m, p1 = i, p2 = i (20)
                rtype = 20
              END IF
            END IF
          END IF
        ELSE                                ! r1 = i, r2 = i
          IF ( p1_m .EQV. .TRUE. ) THEN     ! r1 = i, r2 = i, p1 = m
            IF ( p2_0 .EQV. .TRUE. ) THEN   ! r1 = i, r2 = i, p1 = m, p2 = 0 (7)
              rtype = 7
            ELSE
              IF ( p2_m .EQV. .TRUE. ) THEN ! r1 = i, r2 = i, p1 = m, p2 = m (21)
                rtype = 21
              ELSE                          ! r1 = i, r2 = i, p1 = m, p2 = i (22)
                rtype = 22
              END IF
            END IF
          ELSE                              ! r1 = i, r2 = i, p1 = i
            IF ( p2_0 .EQV. .TRUE. ) THEN   ! r1 = i, r2 = i, p1 = i, p2 = 0 (8)
              rtype = 8
            ELSE
              IF ( p2_m .EQV. .TRUE. ) THEN ! r1 = i, r2 = i, p1 = i, p2 = m (23)
                rtype = 23
              ELSE                          ! r1 = i, r2 = i, p1 = i, p2 = i (24)
                rtype = 24
              END IF
            END IF
          END IF
        END IF
      END IF
    ELSE                                    ! r1 = sp
      IF ( r2_m .EQV. .TRUE. ) THEN         ! r1 = sp, r2 = m
        IF ( p1_m .EQV. .TRUE. ) THEN       ! r1 = sp, r2 = m, p1 = m
          IF ( p2_0 .EQV. .TRUE. ) THEN     ! r1 = sp, r2 = m, p1 = m, p2 = 0
            rtype = 25
          ELSE
            IF ( p2_m .EQV. .TRUE. ) THEN   ! r1 = sp, r2 = m, p1 = m, p2 = m
              rtype = 26
            ELSE                            ! r1 = sp, r2 = m, p1 = m, p2 = i
              rtype = 27
            END IF
          END IF
        ELSE                                ! r1 = sp, r2 = m, p1 = i
          IF ( p2_0 .EQV. .TRUE. ) THEN     ! r1 = sp, r2 = m, p1 = i, p2 = 0
            rtype = 28
          ELSE
            IF ( p2_m .EQV. .TRUE. ) THEN   ! r1 = sp, r2 = m, p1 = i, p2 = m
              rtype = 29
            ELSE                            ! r1 = sp, r2 = m, p1 = i, p2 = i
              rtype = 30
            END IF
          END IF
        END IF
      ELSE                                  ! r1 = sp, r2 = i
        IF ( p1_m .EQV. .TRUE. ) THEN       ! r1 = sp, r2 = i, p1 = m
          IF ( p2_0 .EQV. .TRUE. ) THEN     ! r1 = sp, r2 = i, p1 = m, p2 = 0
            rtype = 31
          ELSE
            IF ( p2_m .EQV. .TRUE. ) THEN   ! r1 = sp, r2 = i, p1 = m, p2 = m
              rtype = 32
            ELSE                            ! r1 = sp, r2 = i, p1 = m, p2 = i
              rtype = 33
            END IF
          END IF
        ELSE                                ! r1 = sp, r2 = i, p1 = i
          IF ( p2_0 .EQV. .TRUE. ) THEN     ! r1 = sp, r2 = i, p1 = i, p2 = 0
            rtype = 34
          ELSE
            IF ( p2_m .EQV. .TRUE. ) THEN   ! r1 = sp, r2 = i, p1 = i, p2 = m
              rtype = 35
            ELSE                            ! r1 = sp, r2 = i, p1 = i, p2 = i
              rtype = 36
            END IF
          END IF
        END IF
      END IF
    END IF
  END SUBROUTINE rtype_tree

  SUBROUTINE place_product(i_re,j_re,k_re,i_re2,j_re2,k_re2,index,index2,r1,r2,&
                           prods,casetype,matrix, wait_list,wait_len,en_list,time,&
                           prod_coords)
  !
  ! Purpose:
  !  The purpose of this subroutine is to place the products of
  !  a reaction in the matrix.
  !
  !  Documentation:
  !  DATE          PROGRAMMER           DESCRIPTION
  !  ========      ==========           ===========
  !  20160216      C. Shingledecker     Original code
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! PLACE_PRODUCT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    INTEGER, INTENT(IN)                                       :: i_re,j_re,k_re
    INTEGER, INTENT(IN)                                       :: i_re2,j_re2,k_re2
    INTEGER, INTENT(INOUT)                                    :: index,index2
    INTEGER, INTENT(IN)                                       :: r1,r2 ! Reactants 1 & 2
    INTEGER, INTENT(IN) , DIMENSION(3)                        :: prods ! Array of products
    INTEGER, INTENT(IN)                                       :: casetype ! The type of the reaction
    INTEGER                               , POINTER           :: wait_len
    INTEGER             , DIMENSION(:,:,:), POINTER           :: matrix
    REAL(KIND=DBL)                        , POINTER           :: time
    REAL                , DIMENSION(:)    , POINTER           :: en_list !list of binding energies
    TYPE(wait_info)     , DIMENSION(:)    , POINTER           :: wait_list
    INTEGER                                                   :: p1,p2
    INTEGER             , DIMENSION(3)                        :: third_coords
    INTEGER, INTENT(OUT), DIMENSION(3,3)           , OPTIONAL :: prod_coords
    INTEGER(KIND=SHORT)                                       :: null


    third_coords = 0
    p1 = prods(1)
    p2 = prods(2)

 !   IF ( i_re .EQ. 340 .AND. j_re .EQ. 71 .AND. k_re .EQ. 141 ) DEBUG = .TRUE.
 !   IF ( i_re2 .EQ. 340 .AND. j_re2 .EQ. 71 .AND. k_re2 .EQ. 141 ) DEBUG = .TRUE.

    IF ( O3_ANALYTICS .EQV. .TRUE. ) THEN
    	IF ( r1 .EQ. 7 .OR. r2 .EQ. 7 .OR. p1 .EQ. 7 .OR. p2 .EQ. 7 ) THEN
	    WRITE(O3_NUM,*) r1,',',r2,',',prods(1),',',prods(2),',',prods(3),',',time*CR_FLUX
	END IF
    END IF

    IF ( DEBUG .EQV. .TRUE.) PRINT *, 'In Place_Product, case=',casetype

    SELECT CASE (casetype)
    CASE (1)
      !***************************************************************************
      ! r1 = m, r2 = m, p1 = m, p2 = 0
      !***************************************************************************
      ! Place p1 at r2 site (since p1 = m, leave index2 at site as-is)
      ! Update sp_num
      ! Update wait_time in wait_list
      ! Remove r1 from list
      ! Make r1 site empty
      !***************************************************************************
      wait_list(index2)%sp_num = p1
      CALL wait_calc( wait_list,index2,en_list,time )
      CALL reactant_remove( wait_list,index,matrix,wait_len )
      matrix(i_re,j_re,k_re)    = 0
      matrix(i_re2,j_re2,k_re2) = index2
    CASE (2)
      !***************************************************************************
      ! r1 = m, r2 = m, p1 = i, p2 = 0
      !***************************************************************************
      ! Remove r1 from list
      ! Update index2
      ! Remove r2 from list
      ! Place product at r2 coords
      ! Make r1 site empty
      !***************************************************************************
      CALL reactant_remove(wait_list,index,matrix,wait_len)
      index2 = matrix(i_re2,j_re2,k_re2)
      CALL reactant_remove(wait_list,index2,matrix,wait_len)
      matrix(i_re2,j_re2,k_re2) = -1*p1
      matrix(i_re,j_re,k_re)    = 0
    CASE (3)
      !***************************************************************************
      ! r1 = m, r2 = i, p1 = m, p2 = 0
      !***************************************************************************
      ! Place p1 at r2 site (since p1 = m, leave index2 at site as-is)
      ! Update sp_num
      ! Update wait_time in wait_list
      ! Make r1 site empty
      ! Update coordinates in wait_list at index
      !***************************************************************************
      wait_list(index)%sp_num = p1
      CALL wait_calc(wait_list,index,en_list,time)
      matrix(i_re2,j_re2,k_re2) = index
      matrix(i_re,j_re,k_re)    = 0
      wait_list(index)%i = i_re2
      wait_list(index)%j = j_re2
      wait_list(index)%k = k_re2
    CASE (4)
      !***************************************************************************
      ! r1 = m, r2 = i, p1 = i, p2 = 0
      !***************************************************************************
      ! Remove r1 from list
      ! Place p1 at r2
      ! Make r1 site empty
      !***************************************************************************
      CALL reactant_remove(wait_list,index,matrix,wait_len)
      matrix(i_re2,j_re2,k_re2) = -1*p1
      matrix(i_re,j_re,k_re)    = 0
    CASE (5)
      !***************************************************************************
      ! r1 = i, r2 = m, p1 = m, p2 = 0
      !***************************************************************************
      ! Update sp_num at index2 to p1
      ! Call wait_calc for p1 at index2
      ! Make r1 site empty
      !***************************************************************************
      wait_list(index2)%sp_num  = p1
      CALL wait_calc(wait_list,index2,en_list,time)
      matrix(i_re,j_re,k_re)    = 0
      matrix(i_re2,j_re2,k_re2) = index2
    CASE (6)
      !***************************************************************************
      ! r1 = i, r2 = m, p1 = i, p2 = 0
      !***************************************************************************
      ! Remove r2 from list
      ! Place p1 at r2 site
      ! Make r1 site empty
      !***************************************************************************
      CALL reactant_remove(wait_list,index2,matrix,wait_len)
      matrix(i_re2,j_re2,k_re2) = -1*p1
      matrix(i_re,j_re,k_re)    = 0
    CASE (7)
      !***************************************************************************
      ! r1 = i, r2 = i, p1 = m, p2 = 0
      !***************************************************************************
      ! Increase wait_len by 1
      ! Add p1 to wait_list at wait_len
      ! Call wait_calc for p1
      ! Make r2 site wait_len
      ! Make r1 site empty
      !***************************************************************************
      wait_len = wait_len + 1
      wait_list(wait_len)%sp_num = p1
      wait_list(wait_len)%i = i_re2
      wait_list(wait_len)%j = j_re2
      wait_list(wait_len)%k = k_re2
      CALL wait_calc(wait_list,wait_len,en_list,time)
      matrix(i_re2,j_re2,k_re2) = wait_len
      matrix(i_re,j_re,k_re)    = 0
    CASE (8)
      !***************************************************************************
      ! r1 = i, r2 = i, p1 = i, p2 = 0
      !***************************************************************************
      ! Place p1 at r2 site
      ! Make r1 site empty
      !***************************************************************************
      matrix(i_re2,j_re2,k_re2) = -1*p1
      matrix(i_re,j_re,k_re)    = 0
    CASE (9)
      !***************************************************************************
      ! r1 = m, r2 = m, p1 = m, p2 = m
      !***************************************************************************
      ! Make index2 spe_num = p2
      ! Call wait_calc for index2
      ! Make index sp_num = p1
      ! Call wait_calc for index
      !***************************************************************************
      wait_list(index)%sp_num = p1
      CALL wait_calc(wait_list,index,en_list,time)
      wait_list(index2)%sp_num = p2
      CALL wait_calc(wait_list,index2,en_list,time)
      matrix(i_re,j_re,k_re)    = index
      matrix(i_re2,j_re2,k_re2) = index2
    CASE (10)
      !***************************************************************************
      ! r1 = m, r2 = m, p1 = m, p2 = i
      !***************************************************************************
      ! Relace index sp_num = p1
      ! Call wait_calc for index p1
      ! Remove r2 from list
      ! Place p2 at r2 site
      !***************************************************************************
      wait_list(index)%sp_num = p1
      CALL wait_calc(wait_list,index,en_list,time)
      CALL reactant_remove(wait_list,index2,matrix,wait_len)
      matrix(i_re2,j_re2,k_re2) = -1*p2
      matrix(i_re,j_re,k_re)    = index
    CASE (11)
      !***************************************************************************
      ! r1 = m, r2 = m, p1 = i, p2 = m
      !***************************************************************************
      ! Replace index sp_num = p2
      ! Call wait_calc for index p2
      ! Remove index2 from list
      ! Place p1 at r2 site
      !***************************************************************************
      wait_list(index)%sp_num = p2
      CALL wait_calc(wait_list,index,en_list,time)
      matrix(i_re2,j_re2,k_re2) = -1*p1
      matrix(i_re,j_re,k_re)    = index
      CALL reactant_remove(wait_list,index2,matrix,wait_len)
    CASE (12)
      !***************************************************************************
      ! r1 = m, r2 = m, p1 = i, p2 = i
      !***************************************************************************
      ! Remove r1 from list
      ! Update index2
      ! Remove r2 from list
      ! Place p1 at r1 site
      ! Place p2 at r2 site
      !***************************************************************************
      CALL reactant_remove(wait_list,index,matrix,wait_len)
      index2 = matrix(i_re2,j_re2,k_re2)
      CALL reactant_remove(wait_list,index2,matrix,wait_len)
      matrix(i_re,j_re,k_re)    = -1*p1
      matrix(i_re2,j_re2,k_re2) = -1*p2
    CASE (13)
      !***************************************************************************
      ! r1 = m, r2 = i, p1 = m, p2 = m
      !***************************************************************************
      ! Change index sp_num = p1
      ! Call wait_calc for p1
      ! Increase wait_len +1
      ! Call wait_calc for (wait_len)
      ! Place wait_len at r2
      !***************************************************************************
      wait_list(index)%sp_num = p1
      CALL wait_calc(wait_list,index,en_list,time)
      wait_len = wait_len + 1
      wait_list(wait_len)%sp_num = p2
      wait_list(wait_len)%i = i_re2
      wait_list(wait_len)%j = j_re2
      wait_list(wait_len)%k = k_re2
      CALL wait_calc(wait_list,wait_len,en_list,time)
      matrix(i_re,j_re,k_re)    = index
      matrix(i_re2,j_re2,k_re2) = wait_len
    CASE (14)
      !***************************************************************************
      ! r1 = m, r2 = i, p1 = m, p2 = i
      !***************************************************************************
      ! Change index sp_num = p1
      ! Call wait_calc for p1
      ! Place p2 at r2
      !***************************************************************************
      wait_list(index)%sp_num    = p1
      CALL wait_calc(wait_list,index,en_list,time)
      matrix(i_re,j_re,k_re)     = index
      matrix(i_re2,j_re2,k_re2)  = -1*p2
    CASE (15)
      !***************************************************************************
      ! r1 = m, r2 = i, p1 = i, p2 = m
      !***************************************************************************
      ! Change index sp_num = p2
      ! Call wait_calc for p2
      ! Place p1 at r2
      !***************************************************************************
      wait_list(index)%sp_num    = p2
      CALL wait_calc(wait_list,index,en_list,time)
      matrix(i_re,j_re,k_re)     = index
      matrix(i_re2,j_re2,k_re2)  = -1*p1
    CASE (16)
      !***************************************************************************
      ! r1 = m, r2 = i, p1 = i, p2 = i
      !***************************************************************************
      ! Remove r1 from list
      ! Place p1 at r1
      ! Place p2 at r2
      !***************************************************************************
      CALL reactant_remove(wait_list,index,matrix,wait_len)
      matrix(i_re,j_re,k_re)     = -1*p1
      matrix(i_re2,j_re2,k_re2)  = -1*p2
    CASE (17)
      !***************************************************************************
      ! r1 = i, r2 = m, p1 = m, p2 = m
      !***************************************************************************
      ! Change index2 sp_num = p1
      ! Call wait_calc for p1
      ! Increase wait_list +1
      ! Populate sp_num and coords for p2
      ! Call wait_calc for p2
      ! Matrix r1 site = wait_len
      !***************************************************************************
      wait_list(index2)%sp_num   = p1
      CALL wait_calc(wait_list,index2,en_list,time)
      wait_len                   = wait_len + 1
      wait_list(wait_len)%sp_num = p2
      wait_list(wait_len)%i      = i_re
      wait_list(wait_len)%j      = j_re
      wait_list(wait_len)%k      = k_re
      CALL wait_calc(wait_list,wait_len,en_list,time)
      matrix(i_re,j_re,k_re)     = wait_len
      matrix(i_re2,j_re2,k_re2)  = index2
    CASE (18)
      !***************************************************************************
      ! r1 = i, r2 = m, p1 = m, p2 = i
      !***************************************************************************
      ! Change index2 sp_num = p1
      ! Call wait_calc for p1
      ! Save p2 at r1 site
      !***************************************************************************
      wait_list(index2)%sp_num   = p1
      CALL wait_calc(wait_list,index2,en_list,time)
      matrix(i_re,j_re,k_re)     = -1*p2
      matrix(i_re2,j_re2,k_re2)  = index2
    CASE (19)
      !***************************************************************************
      ! r1 = i, r2 = m, p1 = i, p2 = m
      !***************************************************************************
      ! Change index2 sp_num = p2
      ! Call wait_calc for p2
      ! Place p1 at r1
      !***************************************************************************
      wait_list(index2)%sp_num   = p2
      CALL wait_calc(wait_list,index2,en_list,time)
      matrix(i_re,j_re,k_re)     = -1*p1
      matrix(i_re2,j_re2,k_re2)  = index2
    CASE (20)
      !***************************************************************************
      ! r1 = i, r2 = m, p1 = i, p2 = i
      !***************************************************************************
      ! Remove r2 from list
      ! Place p1 at r1
      ! Place p2 at r2
      !***************************************************************************
      CALL reactant_remove(wait_list,index2,matrix,wait_len)
      matrix(i_re,j_re,k_re)     = -1*p1
      matrix(i_re2,j_re2,k_re2)  = -1*p2
    CASE (21)
      !***************************************************************************
      ! r1 = i, r2 = i, p1 = m, p2 = m
      !***************************************************************************
      ! Increase wait_len +1
      ! Populate wait_len with p1
      ! Call wait_calc for p1
      ! Write wait_len at r1
      ! Increase wait_len +1
      ! Populate wait_len with p2
      ! Call wait_calc for p2
      ! Write wait_len at r2
      !***************************************************************************
      wait_len                   = wait_len + 1
      wait_list(wait_len)%sp_num = p1
      wait_list(wait_len)%i      = i_re
      wait_list(wait_len)%j      = j_re
      wait_list(wait_len)%k      = k_re
      CALL wait_calc(wait_list,wait_len,en_list,time)
      matrix(i_re,j_re,k_re)     = wait_len
      wait_len                   = wait_len + 1
      wait_list(wait_len)%sp_num = p2
      wait_list(wait_len)%i      = i_re2
      wait_list(wait_len)%j      = j_re2
      wait_list(wait_len)%k      = k_re2
      CALL wait_calc(wait_list,wait_len,en_list,time)
      matrix(i_re2,j_re2,k_re2)  = wait_len
    CASE (22)
      !***************************************************************************
      ! r1 = i, r2 = i, p1 = m, p2 = i
      !***************************************************************************
      ! Increase wait_len +1
      ! Populate wait_len
      ! Call wait_calc for wait_len
      ! Place wait_len at r1
      ! Place p2 at r2
      !***************************************************************************
      wait_len                   = wait_len + 1
      wait_list(wait_len)%sp_num = p1
      wait_list(wait_len)%i      = i_re
      wait_list(wait_len)%j      = j_re
      wait_list(wait_len)%k      = k_re
      CALL wait_calc(wait_list,wait_len,en_list,time)
      matrix(i_re,j_re,k_re)     = wait_len
      matrix(i_re2,j_re2,k_re2)  = -1*p2
    CASE (23)
      !***************************************************************************
      ! r1 = i, r2 = i, p1 = i, p2 = m
      !***************************************************************************
      ! Increase wait_len + 1
      ! Populate wait_len
      ! Call wait_calc for wait_len
      ! Place wait_calc at r1
      ! Place p1 at r2
      !***************************************************************************
      wait_len                   = wait_len + 1
      wait_list(wait_len)%sp_num = p2
      wait_list(wait_len)%i      = i_re
      wait_list(wait_len)%j      = j_re
      wait_list(wait_len)%k      = k_re
      CALL wait_calc(wait_list,wait_len,en_list,time)
      matrix(i_re,j_re,k_re)     = wait_len
      matrix(i_re2,j_re2,k_re2)  = -1*p1
    CASE (24)
      !***************************************************************************
      ! r1 = i, r2 = i, p1 = i, p2 = i
      !***************************************************************************
      ! Place p1 at r1
      ! Place p2 at r2
      !***************************************************************************
      matrix(i_re,j_re,k_re)     = -1*p1
      matrix(i_re2,j_re2,k_re2)  = -1*p2
    CASE (25)
      !***************************************************************************
      ! r1 = sp, r2 = m, p1 = m, p2 = 0
      !***************************************************************************
      ! Net 0 mobile
      ! Populate wait_list r2 with p1
      ! Call wait_calc for p1
      ! NOTE: r2 is at i_re...
      wait_list(index2)%sp_num = p1
      CALL wait_calc(wait_list,index2,en_list,time)
      matrix(i_re,j_re,k_re) = index2
    CASE (26)
      !***************************************************************************
      ! r1 = sp, r2 = m, p1 = m, p2 = m
      !***************************************************************************
      ! Net +1 mobile
      ! Populate r2 with p1
      ! Call wait_calc for p1
      ! Increase wait_len +1
      ! Populate wait_len with p2
      ! Call wait_calc for p2
      ! Place wait_len at new coords
      ! NOTE: r2 is at i_re...
      wait_list(index2)%sp_num = p1
      wait_list(index2)%i = i_re
      wait_list(index2)%j = j_re
      wait_list(index2)%k = k_re
      CALL wait_calc(wait_list,index2,en_list,time)
      wait_len = wait_len + 1
      wait_list(wait_len)%sp_num = p2
      wait_list(wait_len)%i = i_re2
      wait_list(wait_len)%j = j_re2
      wait_list(wait_len)%k = k_re2
      CALL wait_calc(wait_list,wait_len,en_list,time)
      matrix(i_re,j_re,k_re)    = index2
      matrix(i_re2,j_re2,k_re2) = wait_len
    CASE (27)
      !***************************************************************************
      ! r1 = sp, r2 = m, p1 = m, p2 = i
      !***************************************************************************
      ! Net 0 mobile
      ! Populate r2 with p1
      ! Call wait_calc for p1
      ! Place p2 at new coords
      ! NOTE: r2 is at i_re...
      wait_list(index2)%sp_num = p1
      CALL wait_calc(wait_list,index2,en_list,time)
      matrix(i_re,j_re,k_re)    = index2
      matrix(i_re2,j_re2,k_re2) = -1*p2
    CASE (28)
      !***************************************************************************
      ! r1 = sp, r2 = m, p1 = i, p2 = 0
      !***************************************************************************
      ! Net -1 mobile
      ! Remove r2 from list
      ! Place p1 at r2
      ! NOTE: r2 is at i_re...
      CALL reactant_remove(wait_list,index2,matrix,wait_len)
      matrix(i_re,j_re,k_re) = -1*p1
    CASE (29)
      !***************************************************************************
      ! r1 = sp, r2 = m, p1 = i, p2 = m
      !***************************************************************************
      ! NOTE: r2 is at i_re...
      ! Net 0 mobile
      ! Populate r2 with p2
      ! Call wait_calc for p2
      ! Place p1 at i_re2..
      wait_list(index2)%sp_num = p2
      CALL wait_calc(wait_list,index2,en_list,time)
      matrix(i_re,j_re,k_re)    = index2
      matrix(i_re2,j_re2,k_re2) = -1*p1
    CASE (30)
      !***************************************************************************
      ! r1 = sp, r2 = m, p1 = i, p2 = i
      !***************************************************************************
      ! NOTE: r2 is at i_re...
      ! NOTE: r2 has index2
      ! Net -1 mobile
      ! Remove r2 from list
      ! Place p1 at i_re..
      ! Place p2 at i_re2..
      CALL reactant_remove(wait_list,index2,matrix,wait_len)
      matrix(i_re,j_re,k_re) = -1*p1
      matrix(i_re2,j_re2,k_re2) = -1*p2
    CASE (31)
      !***************************************************************************
      ! r1 = sp, r2 = i, p1 = m, p2 = 0
      !***************************************************************************
      ! NOTE: r2 is at i_re...
      ! NOTE: r2 has index2
      ! Net +1 mobile
      ! Increase wait_len +1
      ! Populate wait_len with i_re...
      ! Call wait_calc for p1
      ! Place p1 at i_re...
      wait_len = wait_len + 1
      wait_list(wait_len)%sp_num = p1
      wait_list(wait_len)%i = i_re
      wait_list(wait_len)%j = j_re
      wait_list(wait_len)%k = k_re
      CALL wait_calc(wait_list,index2,en_list,time)
      matrix(i_re,j_re,k_re) = wait_len
    CASE (32)
      !***************************************************************************
      ! r1 = sp, r2 = i, p1 = m, p2 = m
      !***************************************************************************
      ! NOTE: r2 is at i_re...
      ! NOTE: r2 has index2
      ! Net +2 mobile
      ! Increase wait_len + 1
      ! Populate wait_len with p1 at i_re...
      ! Call wait_calc for p1
      ! Place wait_len at i_re...
      ! Increase wait_len +1
      ! Populate wait_len with p2 at i_re2...
      ! Call wait_calc for p2
      ! Place wait_len at i_re2...
      wait_len = wait_len + 1
      wait_list(wait_len)%sp_num = p1
      wait_list(wait_len)%i = i_re
      wait_list(wait_len)%j = j_re
      wait_list(wait_len)%k = k_re
      CALL wait_calc(wait_list,wait_len,en_list,time)
      matrix(i_re,j_re,k_re) = wait_len
      wait_len = wait_len + 1
      wait_list(wait_len)%sp_num = p2
      wait_list(wait_len)%i = i_re2
      wait_list(wait_len)%j = j_re2
      wait_list(wait_len)%k = k_re2
      CALL wait_calc(wait_list,wait_len,en_list,time)
      matrix(i_re2,j_re2,k_re2) = wait_len
    CASE (33)
      !***************************************************************************
      ! r1 = sp, r2 = i, p1 = m, p2 = i
      !***************************************************************************
      ! NOTE: r2 is at i_re...
      ! NOTE: r2 has index2
      ! Net +1 mobile
      ! Increase wait_len +1
      ! Populate wait_len with p1 at i_re...
      ! Call wait_calc for wait_len
      ! Place wait_len at i_re...
      ! Place p2 at i_re2...
      wait_len = wait_len +1
      wait_list(wait_len)%sp_num = p1
      wait_list(wait_len)%i = i_re
      wait_list(wait_len)%j = j_re
      wait_list(wait_len)%k = k_re
      CALL wait_calc(wait_list,wait_len,en_list,time)
      matrix(i_re,j_re,k_re)    = wait_len
      matrix(i_re2,j_re2,k_re2) = -1*p2
    CASE (34)
      !***************************************************************************
      ! r1 = sp, r2 = i, p1 = i, p2 = 0
      !***************************************************************************
      ! NOTE: r2 is at i_re...
      ! NOTE: r2 has index2
      ! Net 0 mobile
      ! Place p1 at r2
      matrix(i_re,j_re,k_re) = -1*p1
    CASE (35)
      !***************************************************************************
      ! r1 = sp, r2 = i, p1 = i, p2 = m
      !***************************************************************************
      ! NOTE: r2 is at i_re...
      ! NOTE: r2 has index2
      ! Net +1 mobile
      ! Increase wait_len +1
      ! Populate wait_len with p2 at i_re...
      ! Call wait_calc for p2
      ! Place wait_len at i_re...
      ! Place p1 at i_re2...
      wait_len = wait_len + 1
      wait_list(wait_len)%sp_num = p2
      wait_list(wait_len)%i      = i_re
      wait_list(wait_len)%j      = j_re
      wait_list(wait_len)%k      = k_re
      CALL wait_calc(wait_list,wait_len,en_list,time)
      matrix(i_re,j_re,k_re)     = wait_len
      matrix(i_re2,j_re2,k_re2)  = -1*p1
    CASE (36)
      !***************************************************************************
      ! r1 = sp, r2 = i, p1 = i, p2 = i
      !***************************************************************************
      ! NOTE: r2 is at i_re...
      ! NOTE: r2 has index2
      ! Net 0 mobile
      ! Place p1 at i_re...
      ! Place p2 at i_re2...
      matrix(i_re,j_re,k_re)     = -1*p1
      matrix(i_re2,j_re2,k_re2)  = -1*p2
    END SELECT

    !***************************************************************************
    ! Save product coordinates
    !***************************************************************************
    IF ( PRESENT(prod_coords) ) THEN
      ! Coordinates of first product
      prod_coords(1,1) = i_re2
      prod_coords(1,2) = j_re2
      prod_coords(1,3) = k_re2
      ! Coordinates of second product
      IF ( casetype .GT. 8 ) THEN
        prod_coords(2,1) = i_re
        prod_coords(2,2) = j_re
        prod_coords(2,3) = k_re
      END IF
    END IF

    !***************************************************************************
    ! Place 3rd product if necessary
    !***************************************************************************
    ! Find an empty location and save it to prod_coords
    !***************************************************************************
    IF ( prods(3) .NE. 0 ) THEN
      ! Make sure to pass back the coords so it can be checked for ion
      CALL thirdman(prods(3),i_re, j_re, k_re, null, matrix, en_list, time, &
                    wait_list, wait_len, third_coords)
!      IF ( third_coords(1) .EQ. 340 .AND. third_coords(2) .EQ. 71 .AND. third_coords(3) .EQ. 141 ) DEBUG = .TRUE.
      IF ( PRESENT(prod_coords) ) THEN
        ! Save product coordinates
        prod_coords(3,1) = third_coords(1)
        prod_coords(3,2) = third_coords(2)
        prod_coords(3,3) = third_coords(3)
      END IF
    END IF
  END SUBROUTINE place_product
END MODULE subroutines
