MODULE parameters
  USE typedefs
  SAVE
  !******************************************************************************
  ! Misc. Global Variables
  !******************************************************************************
  TYPE(node)       , ALLOCATABLE, TARGET :: MATRIX(:,:,:)
  INTEGER                                :: DIMENS(3)
  DOUBLE PRECISION                       :: TIME
  DOUBLE PRECISION , ALLOCATABLE         :: EN_LIST(:)
  CHARACTER(LEN=10), ALLOCATABLE         :: SP_LIST(:)
  INTEGER          , ALLOCATABLE         :: REACT_CUBE(:,:,:)
  INTEGER          , ALLOCATABLE         :: IONLIST(:)
  INTEGER                                :: NUM_SPECIES
  INTEGER                                :: NUM_REACTS

  !******************************************************************************
  ! Input file names
  !******************************************************************************
  CHARACTER(LEN=80), PARAMETER :: SPECIES_FILE   = 'species.dat'   ! Name of species file
  CHARACTER(LEN=80), PARAMETER :: REACTIONS_FILE = 'reactions.dat' ! Name of reactions file
  CHARACTER(LEN=80), PARAMETER :: PARAMS_FILE    = 'params.dat'    ! Name of constants file (for GP)

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
  DOUBLE PRECISION, PARAMETER :: ENERG       = 100*1E3          ! Ion energy in eV
  DOUBLE PRECISION, PARAMETER :: MP          = 1                ! Ion mass in amu
  DOUBLE PRECISION, PARAMETER :: MO2         = 16               ! Target mass in amu
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
  DOUBLE PRECISION            :: FLUENCE        = 0.d0
  DOUBLE PRECISION, PARAMETER :: TIME_TOTAL        = 1D5    ! Total time in s
  DOUBLE PRECISION, PARAMETER :: ECUTOFF           = 9.0D0  ! Secondary cutoff energy in eV
  DOUBLE PRECISION, PARAMETER :: PCUTOFF           = 5.0D0  ! Primary ion cutoff energy in eV
  DOUBLE PRECISION, PARAMETER :: FLUENCE_TOTAL     = 1.0D16
  DOUBLE PRECISION, PARAMETER :: SUBEXHITPROB      = 0.5
  DOUBLE PRECISION, PARAMETER :: FITNESS_THRESHOLD = 1E20   ! Terminate if FITNESS > this
  DOUBLE PRECISION, PARAMETER :: ELASTIC_LOSS      = 0.0    ! Percent of Etot lost per electron hop
  DOUBLE PRECISION, PARAMETER :: TRL_NU            = 1.0E12 ! Trial frequency, for the rates, in 1/s


  !******************************************************************************
  ! Output File Unit Numbers
  !******************************************************************************
  INTEGER         , PARAMETER :: AB_UNIT_NUM        = 1009                        ! Abundance output file
  INTEGER         , PARAMETER :: RATE_UNIT_NUM      = 1946                        ! Rate output file
  INTEGER         , PARAMETER :: TRACKPLOT_UNIT_NUM = 2016
  INTEGER         , PARAMETER :: REACTIONS_UNIT_NUM = 777

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
  DOUBLE PRECISION            :: DELTA_TIME   = 0.d0
  DOUBLE PRECISION            :: PROTON_ELOSS = 0.d0
  TYPE(rate_info)             :: RATEINFO(8)

  !******************************************************************************
  ! Reference Parameters
  !******************************************************************************
  INTEGER                     :: CRPNUM  = 0 ! Index of primary ion in code
  INTEGER                     :: EXCNUM  = 0 ! Index of excitation in code
  INTEGER                     :: ELECNUM = 0 ! Index of electron in code
  INTEGER                     :: O3NUM   = 7 ! Index of ozone in code
  INTEGER                     :: O2NUM   = 1 ! Index of molecular oxygen in code
  INTEGER                     :: ONUM    = 4 ! Index of atomic oxygen in the code


  !******************************************************************************
  ! Array Variables
  !******************************************************************************
  INTEGER                     :: SPECIAL_LIST(3) = 0
  INTEGER                     :: TIME_FREQ       = 10000  !1000000

  !******************************************************************************
  ! Switches
  !******************************************************************************
  LOGICAL         , PARAMETER :: FIXED_SIZE   = .FALSE.
  LOGICAL         , PARAMETER :: NO_OUTPUT    = .FALSE.
  LOGICAL         , PARAMETER :: QUIET        = .FALSE.
  LOGICAL         , PARAMETER :: SECELEC      = .TRUE.
  LOGICAL         , PARAMETER :: DEBUG        = .FALSE.
  LOGICAL         , PARAMETER :: TRACKPLOT    = .FALSE.
  LOGICAL         , PARAMETER :: O3_ANALYTICS = .FALSE.
  LOGICAL         , PARAMETER :: CALC_RATES   = .FALSE.

  !******************************************************************************
  ! Fitting parameters
  !******************************************************************************
  DOUBLE PRECISION :: O2_DISPROB   = 0      ! O2 dissociation probability
  DOUBLE PRECISION :: O3_DISPROB   = 0      ! O3 dissociation probability
  DOUBLE PRECISION :: AVAL         = 33     ! Parameter for Gamma distribution
  DOUBLE PRECISION :: STEPFAC      = 0.1    ! Determines freq. between colls. for protons

CONTAINS
  SUBROUTINE initconstants ()
    INTEGER :: err
    CHARACTER(LEN=32) :: var
    CHARACTER(LEN=32) :: val

    ! Open file for reading
    OPEN(UNIT=200, FILE=PARAMS_FILE, STATUS='OLD', ACTION='READ', IOSTAT=err)
    IF (err .NE. 0) THEN
       PRINT *, "ERROR: Failed to open params.dat file for reading"
       CALL EXIT(-1)
    END IF

    ! Read in the file
    DO
       READ(200,*,IOSTAT=err) var, val
       IF ( err .NE. 0 ) EXIT
       ! Store value
       SELECT CASE (var)
       CASE ("O")
          CONTINUE
       CASE ("O3")
          CONTINUE
       CASE ("O2_DISPROB")
          READ(val, *) O2_DISPROB
       CASE ("O3_DISPROB")
          READ(val, *) O3_DISPROB
       CASE ("STEPFAC")
          READ(val, *) STEPFAC
       CASE ("AVAL")
          READ(val, *) AVAL
       CASE DEFAULT
          PRINT *, "WARNING: Unexpect variable name '", var, "'. Ignoring..."
       END SELECT
    END DO

    ! NEXIT seems to be the only variable that depended on one of these...
    NEXIT = 10*NSUBEX
    CLOSE(200)
  END SUBROUTINE initconstants

END MODULE parameters
