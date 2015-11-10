MODULE typedefs


TYPE :: wait_info
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Purpose: 
  !   This derived data type is designed to contain
  !  the information related to species waiting times
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE PRECISION    :: wait_time !waiting time
  INTEGER             :: i,j,k     !coordinates in matrix
  INTEGER             :: sp_num    !species identifier
  INTEGER             :: act_type  !1 => hopping, 0 => desorption
END TYPE wait_info

TYPE :: sigma_box
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Purpose: A structure (or "box") to contain the cross-sections that will be used to 
  ! calculate track parameters and energy transfers
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
  DOUBLE PRECISION  :: j_alwd
  DOUBLE PRECISION  :: nu_alwd
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
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Purpose: 
  !  The purpose of this type is to contain the information relevant to a
  ! secondary electron, viz. allowed and forbidden excitation cross-section 
  ! info and ionization cross-section info. The energy is also contained in this
  ! type for convenience. 
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    TYPE(ionstate), ALLOCATABLE, DIMENSION(:) :: se_ionst !Information regarding the ionization states of the target
    TYPE(alwd_exstate), ALLOCATABLE, DIMENSION(:) :: se_alwd !Information on the allowed transitions of the target
    TYPE(fbdn_exstate), ALLOCATABLE, DIMENSION(:) :: se_fbdn !Information on the forbidden transitions of the target
    DOUBLE PRECISION   :: se_energy !The secondary electron energy
    DOUBLE PRECISION   :: se_iontot !Total SE ionization cross-section
    DOUBLE PRECISION   :: se_extot !Total SE excitation cross-section
    DOUBLE PRECISION   :: se_ineltot !Total SE inelastic cross-section
    DOUBLE PRECISION   :: se_alwd_extot !Total allowed excitation cross-section
    DOUBLE PRECISION   :: se_fbdn_extot !Total forbidden excitation cross_section
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: se_ionsigs !Ionization cross-sections
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: se_alwdsigs !Allowed excitation cross-sections
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: se_fbdnsigs !Forbidden exc. cross-sections
END TYPE se_info


END MODULE typedefs
