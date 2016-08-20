PROGRAM sigma_test
! This is a program to generate a set of cross sections over some
! user defined energy range suitable for plotting.
! This test suite used the PSIGMA_SUITE and ESIGMA_SUITE subroutines
! of the CIRISS model.
  USE parameters
  USE specdata
  USE functiondefs
  USE typedefs
  USE subroutines
  USE mc_toolbox
  IMPLICIT NONE
  DOUBLE PRECISION   , DIMENSION(:)    , POINTER :: psigij,psigexj
  TYPE(SIGMA_BOX)    , DIMENSION(:)    , POINTER :: psigmas
  TYPE(SE_INFO)                                  :: se_box
  DOUBLE PRECISION   , DIMENSION(:), ALLOCATABLE, TARGET :: psigij_target,psigexj_target
  TYPE(SIGMA_BOX)    , DIMENSION(:), ALLOCATABLE, TARGET :: psigmas_target
  DOUBLE PRECISION                     , POINTER :: ione
  DOUBLE PRECISION                     , TARGET  :: energy_target

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
  PRINT *, "Hello, world"

END PROGRAM sigma_test
