MODULE parameters
  IMPLICIT NONE
  SAVE

  !******************************************************************************
  ! Initial Ion Energy 
  !******************************************************************************
  DOUBLE PRECISION, PARAMETER :: EINIT       = 100D3 ! Initial ion energy in eV

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
  REAL            , PARAMETER :: CDIM        = 3.414e-10      ! Crystal dimension in m
  REAL            , PARAMETER :: BDIM        = 6.668e-10      !    "
  REAL            , PARAMETER :: ADIM        = 9.225e-10      !    "
  REAL            , PARAMETER :: BETACRYS    = 85.05          ! Beta parameter in deg
  REAL            , PARAMETER :: C_PR        = CDIM*COS(90-BETACRYS) ! Actual height of the crystal cube
  REAL            , PARAMETER :: RHO         = 4.78E27        ! Crystal density in m^-3
  DOUBLE PRECISION, PARAMETER :: RHO2        =0.0286 !in Angstrom^-3

  !******************************************************************************
  ! Physical Conditions 
  !******************************************************************************
  REAL            , PARAMETER :: THICK       = 1.0e-7 !1.0e-5 ! Thickness of the ice in m
  REAL            , PARAMETER :: EDGE        = 1.0e-7 !1.0e-7 ! The edge of the crystal in m 
  REAL            , PARAMETER :: KIN_TEMP    = 5              ! Kinetic temperature in Kelvin
  REAL            , PARAMETER :: CR_FLUX     = 3.8E12         ! Proton/Cosmic-ray flux in n(H+) cm^-2 s^-1
  REAL            , PARAMETER :: AREA        = 1E-12          ! Area of irradiated surface in cm 

  !******************************************************************************
  ! Diffusion Energy Fractions 
  !******************************************************************************
  REAL          , PARAMETER :: E_SURF      = 0.5    ! Surface diffusion energy fraction
  REAL          , PARAMETER :: E_BULK      = 0.7    ! Bulk diffusion energy fraction 
  

  !******************************************************************************
  ! Cross-sections 
  !******************************************************************************
  ! NB: These are now obsolete, as they are calculated at the beginning of the model
  ! and after every energy loss event. 
!  REAL          , PARAMETER :: sigma_i     = 1.68916E-16 ! Ionization proton cross-section in m^2
!  REAL          , PARAMETER :: sigma_e     = 1.06891E-19 ! Excitation proton cross-section in m^2
!  REAL          , PARAMETER :: sigma_el    = 7.258E-20   ! Elastic proton cross-section in m^2

  !******************************************************************************
  ! Kinetic Parameters 
  !******************************************************************************
  REAL          , PARAMETER :: TRL_NU      = 1E10     ! Trial frequency, for the rates, in 1/s
  REAL          , PARAMETER :: DISPROB     = 0.0      ! Probability of excitative dissociation

  !******************************************************************************
  ! Model Parameters 
  !******************************************************************************
  INTEGER       , PARAMETER :: IONS        = 6        ! Number of anions in species list
  INTEGER       , PARAMETER :: TIME_COUNTS = 2        ! Number of times the model will check abundances
  INTEGER       , PARAMETER :: NSGSE       = 10       ! Number of second-generation secondary electrons
  REAL(KIND=DBL), PARAMETER :: TIME_TOTAL  = 1D2      ! Total time in s
  REAL(KIND=DBL), PARAMETER :: AVAL        = 13.0     ! Parameter for Gamma distribution
  REAL(KIND=DBL), PARAMETER :: ECUTOFF     = 3D0      ! Secondary cutoff energy in eV
  REAL(KIND=DBL), PARAMETER :: PCUTOFF     = 5D0      ! Primary ion cutoff energy in eV
  
  !******************************************************************************
  ! Array Parameters 
  !******************************************************************************
  INTEGER       , PARAMETER, DIMENSION(1) ::  FAST_REACTS = (/ 4 /) 

  !******************************************************************************
  ! Output File Unit Numbers 
  !******************************************************************************
  INTEGER       , PARAMETER :: AB_UNIT_NUM = 1009            ! Abundance output file

END MODULE parameters 
