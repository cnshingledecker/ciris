MODULE gp
  USE parameters
CONTAINS

  SUBROUTINE fitness(unfit,o3_prod,o3_dest,fluence,total_fitness)
  ! Purpose:
  !    This subroutine measures the fitness of the current simulation.
  !  the resulting fitness is added to the total fitness thusfar.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! FITNESS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    LOGICAL                       :: unfit
    INTEGER,              POINTER :: o3_prod,o3_dest
    REAL(KIND=DBL)                :: fluence
    REAL(KIND=DBL),   INTENT(OUT) :: total_fitness
    INTEGER :: err
    CHARACTER(LEN=80), PARAMETER  :: FITNESS_FILE       = 'fitness_results'
    REAL(KIND=DBL)                :: denom
    REAL(KIND=DBL)                :: objective
    REAL(KIND=DBL)                :: model
    REAL(KIND=DBL)                :: fit

    denom          = THICK*EDGE*EDGE*1E20    ! volume * 1E20
    objective      = (4*(fluence**0.8))/(1E13**0.8+fluence**0.8) ! hard-coded expected value (objective) function
    model          = REAL(o3_prod-o3_dest)/denom
    fit = (objective - model) ** 2
    total_fitness  = total_fitness + fit

!    PRINT *, 'F_obj(', fluence, ') = ', objective
!    PRINT *, 'F_model(', fluence, ') = ', model
!    PRINT *, 'Total fitness: ', total_fitness
!    PRINT *, '***********************************************************************'

    ! if current solution's total_fitness score is too big, save time and end the simulation
    IF (total_fitness > FITNESS_THRESHOLD) THEN
        unfit = .TRUE.
    END IF

    ! save results to a file
    OPEN(UNIT=201, FILE=FITNESS_FILE, ACCESS='APPEND', ACTION='WRITE', IOSTAT=err)
    IF (err .NE. 0) THEN
        PRINT *, "ERROR: Failed to open fitness_results file for reading"
    ELSE
        WRITE(201, *) fluence, objective, model, total_fitness
        CLOSE(201)
    END IF
  END SUBROUTINE fitness

END MODULE gp
