MODULE parameters
  USE typedefs
  SAVE
  !******************************************************************************
  ! Misc. Global Variables
  !******************************************************************************
  TYPE(node)       , ALLOCATABLE, TARGET :: MATRIX(:,:,:)
  INTEGER                                :: DIMENS(3)
  INTEGER                                :: PRODS(3)
  DOUBLE PRECISION , ALLOCATABLE         :: EN_LIST(:)
  CHARACTER(LEN=10), ALLOCATABLE         :: SP_LIST(:)
  INTEGER          , ALLOCATABLE         :: IONLIST(:)
  INTEGER          , ALLOCATABLE         :: SP_PROD_DEST(:,:)
  integer          , allocatable         :: fast_reacts(:)
  INTEGER                                :: NUM_SPECIES = 0
  INTEGER                                :: NUM_REACTS  = 0
  INTEGER                                :: NUMPROTONS  = 0
  DOUBLE PRECISION                       :: TIME        = 0
  TYPE(reaction), POINTER                :: RE_HEAD
  integer :: initialatoms = 0
  REAL                                   :: GVALUE(3) = 0
  REAL                                   :: CUIRCT(3) = 0
  REAL                                   :: CUPREL = 0

  !******************************************************************************
  ! Input file names
  !******************************************************************************
  CHARACTER(LEN=80), PARAMETER :: SPECIES_FILE   = 'species.dat'   ! Name of species file
  CHARACTER(LEN=80), PARAMETER :: REACTIONS_FILE = 'reactions.dat' ! Name of reactions file
  CHARACTER(LEN=80), PARAMETER :: PARAMS_FILE    = 'params.dat'    ! Name of constants file (for GP)
  CHARACTER(LEN=80), PARAMETER :: GEMINACY_FILE  = 'geminacies.wsv'

  !******************************************************************************
  ! Initial Ion Energy
  !******************************************************************************
  DOUBLE PRECISION, PARAMETER :: EINIT    = 100D3                      ! Initial ion energy in eV

  !******************************************************************************
  ! Matrix/Crystal Structure Parameters
  !******************************************************************************
  DOUBLE PRECISION, PARAMETER :: CDIM     = 3.414e-8                    ! Crystal dimension in cm
  DOUBLE PRECISION, PARAMETER :: BDIM     = 6.668e-8                    !    "
  DOUBLE PRECISION, PARAMETER :: ADIM     = 9.225e-8                    !    "
  DOUBLE PRECISION, PARAMETER :: BETACRYS = 85.05                       ! Beta parameter in deg
  DOUBLE PRECISION, PARAMETER :: C_PR     = CDIM*COS(90-BETACRYS)       ! Actual height of the crystal cube
  DOUBLE PRECISION, PARAMETER :: RHO      = 1.313E22                    !4.78E27  ! Crystal density in cm^-3
  DOUBLE PRECISION, PARAMETER :: RHO2     = 0.01313                     !0.0286   !in Angstrom^-3

  !******************************************************************************
  ! Physical Conditions
  !******************************************************************************
  INTEGER         , PARAMETER :: FIX1     = 150
  INTEGER         , PARAMETER :: FIX2     = 150
  INTEGER         , PARAMETER :: FIX3     = 150
  DOUBLE PRECISION, PARAMETER :: THICK    = 1.0e-5                      !1.0e-5 ! Thickness of the ice in cm
  DOUBLE PRECISION, PARAMETER :: EDGE     = 3.5e-6                      !1.0e-7 ! The edge of the crystal in cm
  DOUBLE PRECISION, PARAMETER :: VOLUME   = THICK*EDGE*EDGE             ! Volume of ice chunk
  DOUBLE PRECISION, PARAMETER :: KIN_TEMP = 5.0D0                       ! Kinetic temperature in Kelvin
  DOUBLE PRECISION, PARAMETER :: AREA     = EDGE*EDGE                   ! Area of irradiated surface in cm
  DOUBLE PRECISION, PARAMETER :: CR_FLUX  = 1.0D11                      ! Proton/Cosmic-ray flux in n(H+) cm^-2 s^-1
  DOUBLE PRECISION, PARAMETER :: CR_RATE  = CR_FLUX*AREA                ! Rate of proton arrival
  DOUBLE PRECISION, PARAMETER :: NELEM    = 3.0*RHO*(THICK*EDGE*EDGE)   !Total matrix elements
  DOUBLE PRECISION, PARAMETER :: TER      = THICK/EDGE                  !Thick to edge ratio
  INTEGER         , PARAMETER :: NEDGE    = FLOOR((NELEM/TER)**(1./3.)) !Edge elements
  INTEGER         , PARAMETER :: NTHICK   = FLOOR(NEDGE*TER)            !Thickness elements

  !******************************************************************************
  ! Diffusion Energy Fractions
  !******************************************************************************
  REAL            , PARAMETER :: E_SURF      = 0.5                         ! Surface diffusion energy fraction
  REAL            , PARAMETER :: E_BULK      = 0.7                         ! Bulk diffusion energy fraction

  !******************************************************************************
  ! Kinetic Parameters
  !******************************************************************************
  DOUBLE PRECISION, PARAMETER :: SHORT_TIME  = 100.0            ! K
  DOUBLE PRECISION, PARAMETER :: ZP          = 1.D0             ! Proton number
  DOUBLE PRECISION, PARAMETER :: ZO1         = 8.D0             ! Atomic oxygen proton number
  DOUBLE PRECISION, PARAMETER :: ZO2         = 16.D0            ! Molecular oxygen proton number
  DOUBLE PRECISION, PARAMETER :: ZO3         = 24.D0            ! Molecular oxygen proton number
  DOUBLE PRECISION, PARAMETER :: ENERG       = 100*1E3          ! Ion energy in eV
  DOUBLE PRECISION, PARAMETER :: MP          = 1                ! Ion mass in amu
  DOUBLE PRECISION, PARAMETER :: MO3         = 24               ! Target mass in amu
  DOUBLE PRECISION, PARAMETER :: MO2         = 16               ! Target mass in amu
  DOUBLE PRECISION, PARAMETER :: MO1         = 8               ! Target mass in amu
  DOUBLE PRECISION, PARAMETER :: A0          = 0.529177         ! Bohr radius in Angstroms
  DOUBLE PRECISION, PARAMETER :: ECHARG2     = 14.39            ! Fundamental charge**2 in eV*Angstroms
  DOUBLE PRECISION, PARAMETER :: Q0          = 6.513E-14        ! (eV*cm)**2
  DOUBLE PRECISION, PARAMETER :: PI          = 4.D0*DATAN(1.D0) ! Pi
  DOUBLE PRECISION, PARAMETER :: EBASE       = EXP(1.D0)        ! Natural base

  !******************************************************************************
  ! Model Parameters
  !******************************************************************************
  INTEGER         , PARAMETER :: IONS              = 6      ! Number of anions in species list
  INTEGER         , PARAMETER :: TIME_COUNTS       = 2      ! Times the model will check abundances
  DOUBLE PRECISION            :: FLUENCE           = 0.d0
  DOUBLE PRECISION, PARAMETER :: TIME_TOTAL        = 1D5    ! Total time in s
  DOUBLE PRECISION, PARAMETER :: ECUTOFF           = 0.98D0 ! Secondary cutoff energy in eV
  DOUBLE PRECISION, PARAMETER :: PCUTOFF           = 4.0D0  ! Primary ion cutoff energy in eV
  DOUBLE PRECISION, PARAMETER :: FLUENCE_TOTAL     = 5.0D17
  DOUBLE PRECISION, PARAMETER :: SUBEXHITPROB      = 0.5
  DOUBLE PRECISION, PARAMETER :: FITNESS_THRESHOLD = 1E20   ! Terminate if FITNESS > this
  DOUBLE PRECISION, PARAMETER :: ELASTIC_LOSS      = 0.0    ! Percent of Etot lost per electron hop
  DOUBLE PRECISION, PARAMETER :: TRL_NU            = 1.0E12 ! Trial frequency, for the rates, in 1/s
  INTEGER         , PARAMETER :: STEPFAC           = 1  ! Determines freq. between colls. for protons


  !******************************************************************************
  ! Output File Unit Numbers
  !******************************************************************************
  INTEGER         , PARAMETER :: SPECIES_LUN   = 1
  INTEGER         , PARAMETER :: REACTIONS_LUN = 2
  INTEGER         , PARAMETER :: AB_LUN        = 3 ! Abundance output file
  INTEGER         , PARAMETER :: RATE_LUN      = 4 ! Rate output file
  INTEGER         , PARAMETER :: TRACKPLOT_LUN = 5 ! Track file
  INTEGER         , PARAMETER :: TRACK_AN_LUN  = 7777
  INTEGER         , PARAMETER :: GEMINACY_LUN  = 8888
  !******************************************************************************
  ! Analytics Parameters
  !******************************************************************************
  INTEGER         , PARAMETER :: TRACKMIN     = 5000
  INTEGER         , PARAMETER :: TRACKMAX     = 1000000
  INTEGER                     :: COUNT_COUNT  = 0.0
  INTEGER                     :: FIND_EMPTY_COUNT = 0
  INTEGER                     :: FINDMAX     = 1000
  INTEGER                     :: BI_CALLS     = 0
  INTEGER                     :: O_ABUNDANCE  = 0
  INTEGER                     :: O2_ABUNDANCE = 0
  INTEGER                     :: O3_ABUNDANCE = 0
  INTEGER                     :: N_OHOP       = 0
  INTEGER                     :: N_O3HOP      = 0
  INTEGER                     :: PROTON_COLL  = 0
  DOUBLE PRECISION            :: DELTA_TIME   = 0.d0
  DOUBLE PRECISION            :: PROTON_ELOSS = 0.d0
  DOUBLE PRECISION            :: TOTAL_PROTON_ELOSS = 0.d0
  TYPE(rate_info)             :: RATEINFO(8)

  !******************************************************************************
  ! Reference Parameters
  !******************************************************************************
  INTEGER                     :: CRPNUM  = 0 ! Index of primary ion in code
  INTEGER                     :: EXCNUM  = 0 ! Index of excitation in code
  INTEGER                     :: ELECNUM = 0 ! Index of electron in code
  INTEGER                     :: MNUM    = 0 ! Index of electron in code
  INTEGER                     :: O3NUM   = 7 ! Index of ozone in code
  INTEGER                     :: O2NUM   = 1 ! Index of molecular oxygen in code
  INTEGER                     :: ONUM    = 4 ! Index of atomic oxygen in the code
  INTEGER                     :: OSTARNUM = 11
  INTEGER                     :: O2STARNUM = 10
  INTEGER                     :: O3STARNUM = 12


  !******************************************************************************
  ! Array Variables
  !******************************************************************************
  INTEGER                     :: SPECIAL_LIST(4) = 0
  INTEGER                     :: TIME_FREQ       = 1 !1000000

  !******************************************************************************
  ! Switches
  !******************************************************************************
  LOGICAL         , PARAMETER :: FIXED_SIZE   = .FALSE.
  LOGICAL         , PARAMETER :: NO_OUTPUT    = .FALSE.
  LOGICAL         , PARAMETER :: QUIET        = .FALSE.
  LOGICAL         , PARAMETER :: SECELEC      = .TRUE.
  LOGICAL         , PARAMETER :: DEBUG        = .FALSE.
  LOGICAL         , PARAMETER :: TRACKPLOT    = .FALSE.
  LOGICAL         , PARAMETER :: O3_ANALYTICS = .TRUE.
  LOGICAL         , PARAMETER :: CALC_RATES   = .FALSE.
  LOGICAL         , PARAMETER :: FIX_FREQ     = .FALSE.
  LOGICAL                     :: ISBARRIER    = .FALSE.
  LOGICAL         , PARAMETER :: TRACK_ANALYTICS = .TRUE.
  LOGICAL         , PARAMETER :: FIXED_STEP   = .FALSE.

  !******************************************************************************
  ! Fitting parameters
  !******************************************************************************
  DOUBLE PRECISION :: AVAL         = 33.0d0 !21.0d0     ! Parameter for Gamma distribution
END MODULE parameters
