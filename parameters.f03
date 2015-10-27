MODULE parameters
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
  REAL(KIND=DBL)  , PARAMETER :: CDIM        = 3.414e-8      ! Crystal dimension in cm
  REAL(KIND=DBL)  , PARAMETER :: BDIM        = 6.668e-8      !    "
  REAL(KIND=DBL)  , PARAMETER :: ADIM        = 9.225e-8      !    "
  REAL(KIND=DBL)  , PARAMETER :: BETACRYS    = 85.05          ! Beta parameter in deg
  REAL(KIND=DBL)  , PARAMETER :: C_PR        = CDIM*COS(90-BETACRYS) ! Actual height of the crystal cube
  REAL(KIND=DBL)  , PARAMETER :: RHO         = 1.313E22 !4.78E27  ! Crystal density in cm^-3
  DOUBLE PRECISION, PARAMETER :: RHO2        = 0.01313  !0.0286   !in Angstrom^-3

  !******************************************************************************
  ! Physical Conditions 
  !******************************************************************************
  REAL(KIND=DBL)  , PARAMETER :: THICK       = 1.0e-5 !1.0e-5 ! Thickness of the ice in cm
  REAL(KIND=DBL)  , PARAMETER :: EDGE        = 1.0e-5 !1.0e-7 ! The edge of the crystal in cm 
  REAL(KIND=DBL)  , PARAMETER :: KIN_TEMP    = 5.0D0          ! Kinetic temperature in Kelvin
  REAL(KIND=DBL)  , PARAMETER :: AREA        = EDGE*EDGE      ! Area of irradiated surface in cm 
  REAL(KIND=DBL)  , PARAMETER :: CR_FLUX     = 1.0D11         ! Proton/Cosmic-ray flux in n(H+) cm^-2 s^-1
  REAL(KIND=DBL)  , PARAMETER :: CR_RATE     = CR_FLUX*AREA   ! Rate of proton arrival
  REAL(KIND=DBL)  , PARAMETER :: NELEM       = 3.0*RHO*(THICK*EDGE*EDGE) !Total matrix elements
  REAL(KIND=DBL)  , PARAMETER :: TER         = THICK/EDGE !Thick to edge ratio
  REAL(KIND=DBL)  , PARAMETER :: NEDGE       = (NELEM/TER)**(1./3.) !Edge elements
  REAL(KIND=DBL)  , PARAMETER :: NTHICK      = NEDGE*TER !Thickness elements

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
  DOUBLE PRECISION, PARAMETER :: TRL_NU      = 1E10     ! Trial frequency, for the rates, in 1/s
  DOUBLE PRECISION, PARAMETER :: DISPROB     = 1.0      ! Probability of excitative dissociation
  DOUBLE PRECISION, PARAMETER :: ZP=1.
  DOUBLE PRECISION, PARAMETER :: ZO1=8.
  DOUBLE PRECISION, PARAMETER :: ZO2=16.
  DOUBLE PRECISION, PARAMETER :: ENERG =100*1E3 !1.602E-14 !in eV 
  DOUBLE PRECISION, PARAMETER :: MP=1 !Ion mass in amu 
  DOUBLE PRECISION, PARAMETER :: MO2=16 !Target mass in amu
  DOUBLE PRECISION, PARAMETER :: A0=0.529177 !Bohr radius in Angstroms 
  DOUBLE PRECISION, PARAMETER :: ECHARG2 = 14.39 !Square of the electron charge in eV*Angstroms 
  DOUBLE PRECISION, PARAMETER :: Q0=6.513E-14 ! (eV*cm)**2
  DOUBLE PRECISION, PARAMETER :: PI=4.D0*DATAN(1.D0)
  DOUBLE PRECISION, PARAMETER :: EBASE=EXP(1.D0)


  !******************************************************************************
  ! Model Parameters 
  !******************************************************************************
  INTEGER       , PARAMETER :: IONS        = 6     ! Number of anions in species list
  INTEGER       , PARAMETER :: TIME_COUNTS = 2     ! Times the model will check abundances
  INTEGER       , PARAMETER :: NSUBEX      = 1     ! Number of sub-excitation interactions 
  REAL(KIND=DBL), PARAMETER :: TIME_TOTAL  = 1D5   ! Total time in s
  REAL(KIND=DBL), PARAMETER :: AVAL        = 13.0  ! Parameter for Gamma distribution
  REAL(KIND=DBL), PARAMETER :: ECUTOFF     = 4.5D0 ! Secondary cutoff energy in eV
  REAL(KIND=DBL), PARAMETER :: PCUTOFF     = 5D0   ! Primary ion cutoff energy in eV
  
  !******************************************************************************
  ! Array Parameters 
  !******************************************************************************
  INTEGER       , PARAMETER, DIMENSION(1) ::  FAST_REACTS = (/ 4 /) 

  !******************************************************************************
  ! Output File Unit Numbers 
  !******************************************************************************
  INTEGER       , PARAMETER :: AB_UNIT_NUM = 1009            ! Abundance output file

END MODULE parameters 
