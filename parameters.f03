MODULE parameters
  SAVE

  !******************************************************************************
  ! Input file names
  !******************************************************************************
  CHARACTER(LEN=80), PARAMETER :: SPECIES_FILE   = 'species.dat'   ! Name of species file
  CHARACTER(LEN=80), PARAMETER :: REACTIONS_FILE = 'reactions.dat' ! Name of reactions file

  !******************************************************************************
  ! Initial Ion Energy
  !******************************************************************************
  DOUBLE PRECISION, PARAMETER :: EINIT       = 100D3                      ! Initial ion energy in eV

  !******************************************************************************
  ! Precision Parameters
  !******************************************************************************
  INTEGER         , PARAMETER :: SHORT       = SELECTED_INT_KIND(3)
  INTEGER         , PARAMETER :: LONG        = SELECTED_INT_KIND(9)
  INTEGER         , PARAMETER :: SGL         = SELECTED_REAL_KIND(p=6,r=37)
  INTEGER         , PARAMETER :: DBL         = SELECTED_REAL_KIND(p=13,r=200)

  !******************************************************************************
  ! Matrix/Crystal Structure Parameters
  !******************************************************************************
  REAL(KIND=DBL)  , PARAMETER :: CDIM        = 3.414e-8                    ! Crystal dimension in cm
  REAL(KIND=DBL)  , PARAMETER :: BDIM        = 6.668e-8                    !    "
  REAL(KIND=DBL)  , PARAMETER :: ADIM        = 9.225e-8                    !    "
  REAL(KIND=DBL)  , PARAMETER :: BETACRYS    = 85.05                       ! Beta parameter in deg
  REAL(KIND=DBL)  , PARAMETER :: C_PR        = CDIM*COS(90-BETACRYS)       ! Actual height of the crystal cube
  REAL(KIND=DBL)  , PARAMETER :: RHO         = 1.313E22                    !4.78E27  ! Crystal density in cm^-3
  DOUBLE PRECISION, PARAMETER :: RHO2        = 0.01313                     !0.0286   !in Angstrom^-3

  !******************************************************************************
  ! Physical Conditions
  !******************************************************************************
  REAL(KIND=DBL)  , PARAMETER :: THICK       = 5.0e-6                      !1.0e-5 ! Thickness of the ice in cm
  REAL(KIND=DBL)  , PARAMETER :: EDGE        = 5.0e-6                      !1.0e-7 ! The edge of the crystal in cm
  REAL(KIND=DBL)  , PARAMETER :: KIN_TEMP    = 5.0D0                       ! Kinetic temperature in Kelvin
  REAL(KIND=DBL)  , PARAMETER :: AREA        = EDGE*EDGE                   ! Area of irradiated surface in cm
  REAL(KIND=DBL)  , PARAMETER :: CR_FLUX     = 1.0D11                      ! Proton/Cosmic-ray flux in n(H+) cm^-2 s^-1
  REAL(KIND=DBL)  , PARAMETER :: CR_RATE     = CR_FLUX*AREA                ! Rate of proton arrival
  REAL(KIND=DBL)  , PARAMETER :: NELEM       = 3.0*RHO*(THICK*EDGE*EDGE)   !Total matrix elements
  REAL(KIND=DBL)  , PARAMETER :: TER         = THICK/EDGE                  !Thick to edge ratio
  INTEGER         , PARAMETER :: NEDGE       = FLOOR((NELEM/TER)**(1./3.)) !Edge elements
  INTEGER         , PARAMETER :: NTHICK      = FLOOR(NEDGE*TER)            !Thickness elements

  !******************************************************************************
  ! Diffusion Energy Fractions
  !******************************************************************************
  REAL            , PARAMETER :: E_SURF      = 0.5                         ! Surface diffusion energy fraction
  REAL            , PARAMETER :: E_BULK      = 0.7                         ! Bulk diffusion energy fraction


  !******************************************************************************
  ! Kinetic Parameters
  !******************************************************************************
  DOUBLE PRECISION, PARAMETER :: TRL_NU      = 2.6E8                        ! Trial frequency, for the rates, in 1/s
  DOUBLE PRECISION, PARAMETER :: DISPROB     = 0.0                         ! Probability of excitative dissociation
  DOUBLE PRECISION, PARAMETER :: ZP          = 1.D0                        ! Proton number
  DOUBLE PRECISION, PARAMETER :: ZO1         = 8.D0                        ! Atomic oxygen proton number
  DOUBLE PRECISION, PARAMETER :: ZO2         = 16.D0                       ! Molecular oxygen proton number
  DOUBLE PRECISION, PARAMETER :: ENERG       = 100*1E3                     ! Ion energy in eV
  DOUBLE PRECISION, PARAMETER :: MP          = 1                           ! Ion mass in amu
  DOUBLE PRECISION, PARAMETER :: MO2         = 16                          ! Target mass in amu
  DOUBLE PRECISION, PARAMETER :: A0          = 0.529177                    ! Bohr radius in Angstroms
  DOUBLE PRECISION, PARAMETER :: ECHARG2     = 14.39                       ! Fundamental charge**2 in eV*Angstroms
  DOUBLE PRECISION, PARAMETER :: Q0          = 6.513E-14                   ! (eV*cm)**2
  DOUBLE PRECISION, PARAMETER :: PI          = 4.D0*DATAN(1.D0)            ! Pi
  DOUBLE PRECISION, PARAMETER :: EBASE       = EXP(1.D0)                   ! Natural base


  !******************************************************************************
  ! Model Parameters
  !******************************************************************************
  INTEGER         , PARAMETER :: IONS        = 6                           ! Number of anions in species list
  INTEGER         , PARAMETER :: TIME_COUNTS = 2                           ! Times the model will check abundances
  INTEGER         , PARAMETER :: NSUBEX      = 5                           ! Number of sub-excitation interactions
  INTEGER         , PARAMETER :: NEXIT       = 10*NSUBEX                   ! Max sub-ex loop iters
  REAL            , PARAMETER :: STEPFAC     = 0.01                           ! Determines freq. between colls.
  REAL(KIND=DBL)  , PARAMETER :: TIME_TOTAL  = 1D5                         ! Total time in s
  REAL(KIND=DBL)  , PARAMETER :: AVAL        = 15                          ! Parameter for Gamma distribution
  REAL(KIND=DBL)  , PARAMETER :: ECUTOFF     = 4.5D0                       ! Secondary cutoff energy in eV
  REAL(KIND=DBL)  , PARAMETER :: PCUTOFF     = 5.0D0                       ! Primary ion cutoff energy in eV

  !******************************************************************************
  ! Branching Ratios
  !******************************************************************************
  REAL            , PARAMETER :: O_O2_BRANCHING     = 0.3 ! For O + O2
  REAL            , PARAMETER :: O2_ION_BRANCHING   = 0.0 ! For O2- + O2+
  REAL            , PARAMETER :: O_O2_ION_BRANCHING = 0.0 ! For O+ + O2- and O- + O2+

  !******************************************************************************
  ! Output File Unit Numbers
  !******************************************************************************
  INTEGER         , PARAMETER :: AB_UNIT_NUM = 1009                        ! Abundance output file


  !******************************************************************************
  ! Array Parameters
  !******************************************************************************
  INTEGER, DIMENSION(2)       :: FAST_REACTS  = (/ 4,7 /)                    ! 4 ! Species that react upon formation
  INTEGER, DIMENSION(1)       :: FRAGILE      = (/ 1 /)                    ! 7 ! Species that dissociate easily
  INTEGER, DIMENSION(2)       :: MOBILE_LIST  = (/ 4,7 /)
  INTEGER, DIMENSION(2)       :: SPECIAL_LIST = (/ 20, 21 /)
  INTEGER, PARAMETER          :: TIME_FREQ    = 1000000

  LOGICAL, PARAMETER :: SECELEC    = .TRUE.
  LOGICAL, PARAMETER :: DEBUG      = .FALSE.
  LOGICAL, PARAMETER :: TEST_WRONG = .FALSE.

END MODULE parameters
