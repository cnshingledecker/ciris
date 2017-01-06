MODULE typedefs

  TYPE :: wait_info
     ! Purpose:
     !   This derived data type is designed to contain
     !  the information related to species waiting times
     DOUBLE PRECISION    :: wait_time !waiting time
     INTEGER             :: i,j,k     !coordinates in matrix
     INTEGER             :: sp_num    !species identifier
     INTEGER             :: act_type  !1 => hopping, 0 => desorption
  END TYPE wait_info


  TYPE :: rate_info
     ! Purpose:
     !   This derived data type is designed to contain
     !  the information related to the number of species produced per Δt
     INTEGER             :: r1 ! Reactant 1
     INTEGER             :: r2 ! Reactant 2
     INTEGER             :: count ! The number times this reaction has occured
  END TYPE rate_info

  TYPE :: sigma_box
     ! Purpose: A structure (or "box") to contain the cross-sections that will be used to
     ! calculate track parameters and energy transfers
     DOUBLE PRECISION    :: cross_section
     CHARACTER(len=20)   :: description
  END TYPE sigma_box

  TYPE :: ionstate
     ! Purpose: This type defines the parameters needed to calculate the
     ! cross-sections for ionization associated with specific properties of
     ! a species (i.e. states).
     !
     ! Note: For an example of the values given here, see Jackman, Garvey, &
     ! Green 1977
     !
     ! Note: See the above reference for more detail on the formula to use with
     ! these values
     CHARACTER(len=20) :: termsym_ion !term symbol for the state in LaTeX format
     DOUBLE PRECISION  :: i_energy !Ionization energy in eV
     DOUBLE PRECISION  :: k_ion
     DOUBLE PRECISION  :: kb_ion
     DOUBLE PRECISION  :: j_ion !eV?
     DOUBLE PRECISION  :: jb_ion
     DOUBLE PRECISION  :: jc_ion
     DOUBLE PRECISION  :: gams_ion
     DOUBLE PRECISION  :: gamb_ion
     DOUBLE PRECISION  :: ts_ion
     DOUBLE PRECISION  :: ta_ion
     DOUBLE PRECISION  :: tb_ion
  END TYPE ionstate

  TYPE :: alwd_exstate
     ! Purpose: This type stores data for the allowed, discrete transitions
     ! of a species (i.e. excitations). The variable names come from Porter,
     ! Jackman, Green 1976 and the formula with which these are designed to be
     ! used can be found in Jackman, Garvey, and Green 1977.
     !
     ! For more detail on the physical meaning of these values, see the above
     ! reference.
     CHARACTER(len=20) :: termsym_alwd !term symbol for the state in LaTeX format
     DOUBLE PRECISION  :: wj_alwd !Excitation energy in eV
     DOUBLE PRECISION  :: fj_alwd !oscillator strength
     DOUBLE PRECISION  :: cj_alwd !determined from the oscillator strength
     DOUBLE PRECISION  :: alpha_alwd !fitting parameter
     DOUBLE PRECISION  :: beta_alwd !fitting parameter
  END TYPE alwd_exstate

  TYPE :: fbdn_exstate
     !
     ! Purpose: This type stores parameters for use with the Green and Dutta 1967
     ! formula for calculating the cross-section for forbidden excitation
     ! transitions.
     CHARACTER(len=20) :: termsym_fbdn !term symbol for the state in LaTeX format
     DOUBLE PRECISION  :: wj_fbdn !excitation energy in eV
     DOUBLE PRECISION  :: fj_fbdn !oscillator strength
     DOUBLE PRECISION  :: omega_fbdn !fitting parameter
     DOUBLE PRECISION  :: alpha_fbdn !fitting parameter
     DOUBLE PRECISION  :: beta_fbdn !fitting parameter
  END TYPE fbdn_exstate

  TYPE :: epg_exstate
     !
     ! Purpose: This type contains parameters for calculating the inelastic
     ! excitation cross-sections for protons with the Green-McNeal formula
     ! as described in Edgar, Porter, and Green 1975.
     CHARACTER(len=20) :: termsym_epgex !term symbol for process in LaTeX format
     DOUBLE PRECISION  :: a_epgex
     DOUBLE PRECISION  :: j_epgex !eV, though given in keV in EPG75
     DOUBLE PRECISION  :: nu_epgex
     DOUBLE PRECISION  :: omega_epgex
     DOUBLE PRECISION  :: w_epgex !excitation energy in eV
  END TYPE epg_exstate

  TYPE :: epg_ionstate
     !
     ! Purpose: This type contains data for calculating the inelastic
     ! proton excitation cross-section using the Green-McNeal formalism.
     CHARACTER(len=20) :: termsym_epgion !term symbol for process in LaTeX format
     DOUBLE PRECISION  :: a_epgion
     DOUBLE PRECISION  :: j_epgion !eV, though given in keV in EPG75
     DOUBLE PRECISION  :: nu_epgion
     DOUBLE PRECISION  :: omega_epgion
     DOUBLE PRECISION  :: i_epgion !ionization energy in eV
  END TYPE epg_ionstate

  TYPE :: se_info
     ! Purpose:
     !  The purpose of this type is to contain the information relevant to a
     ! secondary electron, viz. allowed and forbidden excitation cross-section
     ! info and ionization cross-section info. The energy is also contained in this
     ! type for convenience.
     !
     DOUBLE PRECISION   :: se_energy !The secondary electron energy
     DOUBLE PRECISION   :: se_iontot !Total SE ionization cross-section
     DOUBLE PRECISION   :: se_extot !Total SE excitation cross-section
     DOUBLE PRECISION   :: se_ineltot !Total SE inelastic cross-section
     DOUBLE PRECISION   :: se_alwd_extot !Total allowed excitation cross-section
     DOUBLE PRECISION   :: se_fbdn_extot !Total forbidden excitation cross_section
     DOUBLE PRECISION   :: se_o_iontot !Total SE ionization cross-section
     DOUBLE PRECISION   :: se_o_extot !Total SE excitation cross-section
     DOUBLE PRECISION   :: se_o_ineltot !Total SE inelastic cross-section
     DOUBLE PRECISION   :: se_o_alwd_extot !Total allowed excitation cross-section
     DOUBLE PRECISION   :: se_o_fbdn_extot !Total forbidden excitation cross_section
     DOUBLE PRECISION   :: se_o3_extot ! Total excitation cross-section for ozone
     DOUBLE PRECISION   :: se_o3_iontot ! Total ionization cross-section for ozone
     INTEGER, DIMENSION(3) :: parent_coords ! The coordinates at which the electron formed/site of cation
     TYPE(ionstate), ALLOCATABLE, DIMENSION(:) :: se_ionst !Information regarding the ionization states of the target
     TYPE(alwd_exstate), ALLOCATABLE, DIMENSION(:) :: se_alwd !Information on the allowed transitions of the target
     TYPE(fbdn_exstate), ALLOCATABLE, DIMENSION(:) :: se_fbdn !Information on the forbidden transitions of the target
     TYPE(ionstate), ALLOCATABLE, DIMENSION(:) :: se_o_ionst !Information regarding the ionization states of the target
     TYPE(alwd_exstate), ALLOCATABLE, DIMENSION(:) :: se_o_alwd !Information on the allowed transitions of the target
     TYPE(fbdn_exstate), ALLOCATABLE, DIMENSION(:) :: se_o_fbdn !Information on the forbidden transitions of the target
     DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: se_ionsigs !Ionization cross-sections
     DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: se_alwdsigs !Allowed excitation cross-sections
     DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: se_fbdnsigs !Forbidden exc. cross-sections
     DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: se_o_ionsigs !Ionization cross-sections
     DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: se_o_alwdsigs !Allowed excitation cross-sections
     DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: se_o_fbdnsigs !Forbidden exc. cross-sections

  END TYPE se_info

  TYPE :: node
     ! Purpose:
     !   This derived data type is designed to contain
     !  the information related to the ice matrix
     LOGICAL             :: normal
     LOGICAL             :: interstitial
     DOUBLE PRECISION    :: wait_time !waiting time
     INTEGER             :: coord1     !coordinates in matrix
     INTEGER             :: coord2
     INTEGER             :: coord3
     INTEGER             :: sp_num    !species identifier
     INTEGER             :: sec_sp_num ! Secondary species at site
     INTEGER             :: act_type  !1 => hopping, 0 => desorption
     INTEGER             :: hop_dir   ! Direction of hopping
     TYPE (node), POINTER :: before
     TYPE (node), POINTER :: after
     TYPE (node), POINTER :: parent
     INTEGER                   :: leftRight
  END TYPE node

  TYPE :: species
     CHARACTER(len=10) :: name
     INTEGER :: id
     integer :: atoms
     integer :: charge
     DOUBLE PRECISION :: e_d
     LOGICAL :: ion
     logical :: exc
     LOGICAL :: special
     TYPE(species), POINTER :: next
  END TYPE species

  TYPE :: reaction
     CHARACTER(len=10) :: r1
     CHARACTER(len=10) :: r2
     CHARACTER(len=10) :: p1
     CHARACTER(len=10) :: p2
     CHARACTER(len=10) :: p3
     INTEGER :: nr1
     INTEGER :: nr2
     INTEGER :: np1
     INTEGER :: np2
     INTEGER :: np3
     integer :: id
     integer :: count
     DOUBLE PRECISION :: arrh_alpha
     DOUBLE PRECISION :: arrh_beta
     DOUBLE PRECISION :: arrh_gamma
     INTEGER :: rtype
     TYPE(reaction), POINTER :: next
  END TYPE reaction

  TYPE :: reaction_info
     INTEGER :: nr1
     INTEGER :: nr2
     INTEGER :: np1
     INTEGER :: np2
     INTEGER :: np3
     DOUBLE PRECISION :: arrh_alpha
     DOUBLE PRECISION :: arrh_beta
     DOUBLE PRECISION :: arrh_gamma
     INTEGER :: rtype
     integer :: id
  END TYPE reaction_info

END MODULE typedefs
