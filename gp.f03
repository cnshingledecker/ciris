MODULE gp
  USE parameters
CONTAINS

  SUBROUTINE fitness(o3_prod,o3_dest,time,total_fitness)
  ! Purpose:
  !    This subroutine measures the fitness of the current simulation.
  !  the resulting fitness is added to the total fitness thusfar.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! FITNESS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    INTEGER                                          , POINTER :: o3_prod,o3_dest
    REAL(KIND=DBL)                                   , POINTER :: time
    REAL(KIND=DBL)                                             :: fluence

    fluence = CR_FLUX*time

    PRINT *, 'F_obj(', fluence, ') = ', (4*(fluence**0.8))/(1E13**0.8+fluence**0.8) ! hard-coded expected value (objective) function
    PRINT *, 'F_model(', fluence, ') = ', REAL(o3_prod-o3_dest)/denom
    PRINT *, '***********************************************************************'
  END SUBROUTINE fitness

END MODULE gp
