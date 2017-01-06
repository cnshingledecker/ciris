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
  USE bsimple
  IMPLICIT NONE
  INTEGER            , PARAMETER                          :: final_energy=1e6
  DOUBLE PRECISION   , DIMENSION(:)             , POINTER :: psigij,psigexj
  TYPE(SIGMA_BOX)    , DIMENSION(:)             , POINTER :: psigmas
  TYPE(SIGMA_BOX)    , POINTER :: o_psigmas(:)
  TYPE(SIGMA_BOX)    , POINTER :: o3_psigmas(:)
  TYPE(SE_INFO)                                           :: se_box
  DOUBLE PRECISION   , DIMENSION(:), ALLOCATABLE, TARGET  :: psigij_target,psigexj_target
  TYPE(SIGMA_BOX)    , DIMENSION(:), ALLOCATABLE, TARGET  :: psigmas_target
  DOUBLE PRECISION, ALLOCATABLE, TARGET :: o_psigij_target(:),o_psigexj_target(:)
  TYPE(SIGMA_BOX),  ALLOCATABLE, TARGET :: o_psigmas_target(:)
  DOUBLE PRECISION, ALLOCATABLE, TARGET :: o3_psigij_target(:),o3_psigexj_target(:)
  TYPE(SIGMA_BOX),  ALLOCATABLE, TARGET :: o3_psigmas_target(:)
  DOUBLE PRECISION   , POINTER :: o_psigij(:),o_psigexj(:)
  DOUBLE PRECISION   , POINTER :: o3_psigij(:),o3_psigexj(:)
  DOUBLE PRECISION                              , POINTER :: ione
  DOUBLE PRECISION                              , TARGET  :: energy_target
  INTEGER                                                 :: i,j
  DOUBLE PRECISION   , DIMENSION(3)                       :: final_psigs
  DOUBLE PRECISION   , DIMENSION(4)                       :: final_esigs
  DOUBLE PRECISION                                        :: final_esig
  CHARACTER(len=20)  , DIMENSION(4)                       :: se_coll
  DOUBLE PRECISION                                        :: temp_energy

  OPEN(UNIT=1,FILE="proton_sigmas.csv",POSITION='APPEND', STATUS='REPLACE')
  OPEN(UNIT=2,FILE="electron_sigmas.csv",POSITION='APPEND', STATUS='REPLACE')

  ! For O2
  ALLOCATE(psigmas_target(3))
  psigmas => psigmas_target
  ALLOCATE(psigij_target(SIZE(o2_p_ion)))
  psigij => psigij_target
  ALLOCATE(psigexj_target(SIZE(o2_p_ex)))
  psigexj => psigexj_target
  psigmas%cross_section = 0D0
  psigij  = 0D0
  psigexj = 0D0

  ! For O
  ALLOCATE(o_psigmas_target(3))
  o_psigmas => o_psigmas_target
  ALLOCATE(o_psigij_target(SIZE(o_p_ion)))
  o_psigij => o_psigij_target
  ALLOCATE(o_psigexj_target(SIZE(o_p_ex)))
  o_psigexj => o_psigexj_target
  o_psigmas%cross_section = 0D0
  o_psigij  = 0D0
  o_psigexj = 0D0

  ! For O3
  ALLOCATE(o3_psigmas_target(3))
  o3_psigmas => o3_psigmas_target
  ALLOCATE(o3_psigij_target(SIZE(o3_p_ion)))
  o3_psigij => o3_psigij_target
  ALLOCATE(o3_psigexj_target(SIZE(o3_p_ex)))
  o3_psigexj => o3_psigexj_target
  o3_psigmas%cross_section = 0D0
  o3_psigij  = 0D0
  o3_psigexj = 0D0


  ione => energy_target
  ione = EINIT
  temp_energy = 1
  DO WHILE ( temp_energy .LT. final_energy )
    energy_target = temp_energy
    CALL psigma_suite(ione,&
         psigmas,psigij,psigexj,&
         o_psigmas,o_psigij,o_psigexj,&
         o3_psigmas,o3_psigij,o3_psigexj)
    ! For O2
    final_psigs = psigmas%cross_section
    DO j=1,3
      IF ( ISNAN(final_psigs(j) ) ) final_psigs(j) = 0.0
      WRITE(1,*) ione,',',final_psigs(j),',',psigmas(j)%description, ", O2"
    END DO

    ! For O
    final_psigs = o_psigmas%cross_section
    DO j=1,3
       IF ( ISNAN(final_psigs(j) ) ) final_psigs(j) = 0.0
       WRITE(1,*) ione,',',final_psigs(j),',',o_psigmas(j)%description, ", O"
    END DO

    ! For O3
    final_psigs = o3_psigmas%cross_section
    DO j=1,3
       IF ( ISNAN(final_psigs(j) ) ) final_psigs(j) = 0.0
       WRITE(1,*) ione,',',final_psigs(j),',',o3_psigmas(j)%description, ", O3"
    END DO

    temp_energy = temp_energy + temp_energy*0.1
  END DO

  !Now for electrons
  se_coll(1) = ", Ionization"
  se_coll(2) = ", Excitation"
  se_coll(3) = ", Allowed"
  se_coll(4) = ", Forbidden"
  temp_energy = 1
  DO WHILE ( temp_energy .LT. final_energy )
    !Initialize se_box
    se_box%se_energy = temp_energy
    CALL se_info_init(se_box)
    !Calculate electron cross_sections
    CALL esigma_suite(se_box)
    ! For O2
    final_esigs(1) = se_box%se_iontot
    final_esigs(2) = se_box%se_extot
    final_esigs(3) = se_box%se_alwd_extot
    final_esigs(4) = se_box%se_fbdn_extot
    DO j=1,4
      IF ( ISNAN(final_esigs(j) ) ) final_esigs(j) = 0.0
      final_esig = final_esigs(j)
      WRITE(2,*) se_box%se_energy,',',final_esig,se_coll(j), ", O2"
    END DO
    ! For O
    final_esigs = 0
    final_esigs(1) = se_box%se_o_iontot
    final_esigs(2) = se_box%se_o_extot
    final_esigs(3) = se_box%se_o_alwd_extot
    final_esigs(4) = se_box%se_o_fbdn_extot
    DO j=1,4
       IF ( ISNAN(final_esigs(j) ) ) final_esigs(j) = 0.0
       final_esig = final_esigs(j)
       WRITE(2,*) se_box%se_energy,',',final_esig,se_coll(j), ", O"
    END DO
    ! For O3
    final_esigs = 0
    final_esigs(1) = se_box%se_o3_iontot
    final_esigs(2) = se_box%se_o3_extot
    DO j=1,2
       IF ( ISNAN(final_esigs(j) ) ) final_esigs(j) = 0.0
       final_esig = final_esigs(j)
       WRITE(2,*) se_box%se_energy,',',final_esig,se_coll(j), ", O3"
    END DO
    temp_energy = temp_energy + temp_energy*0.1
  END DO

!  PRINT *, se_box%se_ionst
  PRINT *, se_box%se_ionsigs
END PROGRAM sigma_test
