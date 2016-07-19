MODULE gp
  USE parameters
  USE typedefs
CONTAINS

  SUBROUTINE fitness(unfit,o3_prod,o3_dest,fluence,total_fitness,dimens,matrix,wait_list)
  ! Purpose:
  !    This subroutine measures the fitness of the current simulation.
  !  the resulting fitness is added to the total fitness thusfar.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! FITNESS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    LOGICAL                                                    :: unfit
    INTEGER                                          , POINTER :: o3_prod,o3_dest
    REAL(KIND=DBL)                                             :: fluence
    REAL(KIND=DBL)   , INTENT(OUT)                             :: total_fitness
    INTEGER                                                    :: err
    CHARACTER(LEN=80), PARAMETER                               :: FITNESS_FILE = 'fitness_results'
    CHARACTER(LEN=80)                                          :: varfmt
    REAL(KIND=DBL)                                             :: denom
    REAL(KIND=DBL)                                             :: objective
    REAL(KIND=DBL)                                             :: model
    REAL(KIND=DBL)                                             :: fit
    REAL(KIND=DBL)                                             :: part1, part2
    INTEGER                                                    :: ozone_count
    TYPE(wait_info)                , DIMENSION(:)    , POINTER :: wait_list
    INTEGER                        , DIMENSION(:,:,:), POINTER :: matrix
    INTEGER, DIMENSION(3)                                      :: dimens
    INTEGER                                                    :: i, j, k

    ozone_count = 0
    fit = 0

    DO k = 1,dimens(3)
      DO j = 1,dimens(2)
        DO i = 1,dimens(1)
          IF ( matrix(i,j,k) .NE. 0 ) THEN
            IF ( matrix(i,j,k) .LT. 0 ) THEN
              IF ( ABS(matrix(i,j,k)) .EQ. 7 ) THEN
                ozone_count = ozone_count + 1
              END IF
            ELSE IF ( matrix(i,j,k) .GT. 0 ) THEN
              IF  ( wait_list(matrix(i,j,k))%sp_num .EQ. 7 ) THEN
                ozone_count = ozone_count + 1
              END IF
            END IF
          END IF
        END DO
      END DO
    END DO


    denom          = THICK*EDGE*EDGE*1.0E20    ! volume * 1E20
    objective      = (4*(fluence**0.8))/(1E13**0.8+fluence**0.8) ! hard-ncoded expected value (objective) function
    model          = REAL(ozone_count)/denom
    !NB: Original method for calculating fitness
    !This method depends on the number of times counter is called
!    fit = (objective - model) ** 2

    !NB: Alternate method for calculating fitness
    part1 = (((objective + ABS(objective - model)) / objective)*(-100.0)) + 100.0
    part2 = (((LOG10(FLUENCE_TOTAL) + ABS(LOG10(FLUENCE_TOTAL/fluence))) / LOG10(FLUENCE_TOTAL))*(-100.0)) + 100.0
    fit = ABS(part1) + ABS(part2)

    IF ( ISNAN(fit) .EQV. .FALSE. ) total_fitness  = total_fitness + fit
    IF ( QUIET .EQV. .FALSE. ) THEN
      varfmt = "(A7,ES10.4,A5,F10.4)"
      WRITE (*,varfmt), 'F_obj( ', fluence, ' ) = ', objective
      varfmt = "(A9,ES10.4,A5,F10.4)"
      WRITE (*,varfmt), 'F_model( ', fluence, ' ) = ', model
  !    varfmt = "A15,F10.4)"
      PRINT *, 'Total fitness: ', total_fitness
      PRINT *, '***********************************************************************'
    END IF

    ! if current solution's total_fitness score is too big, save time and end the simulation
    IF (total_fitness > FITNESS_THRESHOLD) THEN
        unfit = .TRUE.
    END IF

    ! save results to a file
!    OPEN(UNIT=201, FILE=FITNESS_FILE, ACCESS='APPEND', ACTION='WRITE', IOSTAT=err)
    OPEN(UNIT=201, FILE=FITNESS_FILE, STATUS='REPLACE', ACTION='WRITE', IOSTAT=err)
    IF (err .NE. 0) THEN
        PRINT *, "ERROR: Failed to open fitness_results file for writing"
    ELSE
        WRITE(201, *) fluence, objective, model, total_fitness
        CLOSE(201)
    END IF
  END SUBROUTINE fitness

  SUBROUTINE store_rand()
  ! Purpose:
  !    This subroutine stores the currently used random seed in a file.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! STORE_RAND !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    INTEGER :: seed_size
    INTEGER, ALLOCATABLE :: seed(:)

    INTEGER :: err
    CHARACTER(LEN=80), PARAMETER  :: SEED_FILE       = 'seed'

    ! get seed and save it to a file
    CALL RANDOM_SEED(size=seed_size)
    ALLOCATE(seed(seed_size))
    CALL RANDOM_SEED(get=seed)

    OPEN(UNIT=201, FILE=SEED_FILE, ACCESS='APPEND', ACTION='WRITE', IOSTAT=err)
    IF (err .NE. 0) THEN
        PRINT *, "ERROR: Failed to open seed file for writing"
    ELSE
        WRITE(201, *) seed
        CLOSE(201)
    END IF
    DEALLOCATE(seed)
  END SUBROUTINE store_rand

END MODULE gp
