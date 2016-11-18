MODULE gp
  USE parameters
  USE typedefs
  USE subroutines
CONTAINS
  SUBROUTINE fitness(unfit,fluence,total_fitness)
    ! Purpose:
    !    This subroutine measures the fitness of the current simulation.
    !  the resulting fitness is added to the total fitness thusfar.
    !
    !! FITNESS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE
    LOGICAL                                                    :: unfit
    DOUBLE PRECISION                                             :: fluence
    DOUBLE PRECISION   , INTENT(OUT)                             :: total_fitness
    INTEGER                                                    :: err
    CHARACTER(LEN=80), PARAMETER                               :: FITNESS_FILE = 'fitness_results'
    CHARACTER(LEN=80)                                          :: varfmt
    DOUBLE PRECISION                                             :: denom
    DOUBLE PRECISION                                             :: objective
    DOUBLE PRECISION                                             :: model
    DOUBLE PRECISION                                             :: fit
    DOUBLE PRECISION                                             :: part1, part2

    CALL counter()

    fit       = 0
    denom     = THICK*EDGE*EDGE*1.0E20    ! volume * 1E20

    ! hard-ncoded expected value (objective) function
    objective = (4*(fluence**0.8))/(1E13**0.8+fluence**0.8) 
    model     = REAL(O3_ABUNDANCE)/denom

    !NB: Alternate method for calculating fitness
    part1 = (((objective + ABS(objective - model)) / objective)*(-100.0)) + 100.0
    part2 = (((LOG10(FLUENCE_TOTAL) + ABS(LOG10(FLUENCE_TOTAL/fluence))) / &
         LOG10(FLUENCE_TOTAL))*(-100.0)) + 100.0
    fit = ABS(part1) + ABS(part2)

    IF ( ISNAN(fit) .EQV. .FALSE. ) total_fitness  = total_fitness + fit
    IF ( QUIET .EQV. .FALSE. ) THEN
       varfmt = "(A7,ES10.4,A5,F10.4)"
       WRITE (*,varfmt) 'F_obj( ', fluence, ' ) = ', objective
       varfmt = "(A9,ES10.4,A5,F10.4)"
       WRITE (*,varfmt) 'F_model( ', fluence, ' ) = ', model
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
    varfmt = "(A9,ES10.4,A12,ES10.4)"
    WRITE(*,varfmt) "FITNESS=",total_fitness,"at FLUENCE=",fluence
    OPEN(UNIT=201, FILE=FITNESS_FILE, STATUS='REPLACE', ACTION='WRITE', IOSTAT=err)
    IF (err .NE. 0) THEN
       PRINT *, "ERROR: Failed to open fitness_results file for writing"
       CALL EXIT()
    ELSE
       WRITE(201, *) fluence, objective, model, total_fitness
       CLOSE(201)
    END IF
  END SUBROUTINE fitness

  SUBROUTINE store_rand()
    ! Purpose:
    !    This subroutine stores the currently used random seed in a file.
    !
    !! STORE_RAND !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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


