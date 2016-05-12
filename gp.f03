MODULE gp
  USE parameters
  USE typedefs
CONTAINS

  SUBROUTINE fitness(o3_prod,o3_dest,fluence,dimens,matrix,wait_list)!,total_fitness)
  ! Purpose:
  !    This subroutine measures the fitness of the current simulation.
  !  the resulting fitness is added to the total fitness thusfar.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! FITNESS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    INTEGER                                          , POINTER :: o3_prod,o3_dest
    INTEGER                        , DIMENSION(:,:,:), POINTER :: matrix
    INTEGER, DIMENSION(3)                                      :: dimens
    REAL(KIND=DBL)                                             :: denom
    REAL(KIND=DBL)                                             :: fluence
    !INTEGER                                          , POINTER :: total_fitness
    INTEGER                                                    :: i,j,k
    INTEGER                                                    :: ozone_count
    TYPE(wait_info)                , DIMENSION(:)    , POINTER :: wait_list

    ! volume * 1E20
    ozone_count = 0
    denom = THICK*EDGE*EDGE*1E20
    ! total_fitness = total_fitness + some_comparative_diff

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



    PRINT *, 'F_obj(', fluence, ') = ', (4*(fluence**0.8))/(1E13**0.8+fluence**0.8) ! hard-coded expected value (objective) function
    PRINT *, 'F_model(', fluence, ') = ', REAL(ozone_count)/denom
    PRINT *, '***********************************************************************'
  END SUBROUTINE fitness

END MODULE gp
