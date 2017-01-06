MODULE functiondefs
  USE parameters
  USE typedefs
CONTAINS
  FUNCTION arrhenius(a,b,c)
    ! This function calculates rate-coefficients based on the 
    ! Arrhenius formula
    DOUBLE PRECISION :: arrhenius
    DOUBLE PRECISION :: a, b, c

    arrhenius = a*((KIN_TEMP/300.0d0)**b)*EXP(c/KIN_TEMP)
    RETURN
  END FUNCTION arrhenius

  FUNCTION green_mcneal(energy,a,j,nu,omega,z,i)
    !
    !  Purpose:
    !    To calculate the Green-McNeal scaled proton cross-sections
    !  as described in Miller & Green 1971.
    !
    !  Note:
    !    This formula can be used for both excitation and ionization
    !  cross-sections, given the appropriate input.
    !
    !! GREEN_MCNEAL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    ! Data dictionary: declare calling parameters
    DOUBLE PRECISION             :: green_mcneal
    DOUBLE PRECISION, INTENT(IN) :: energy
    DOUBLE PRECISION, INTENT(IN) :: a
    DOUBLE PRECISION, INTENT(IN) :: j !units of eV
    DOUBLE PRECISION, INTENT(IN) :: nu
    DOUBLE PRECISION, INTENT(IN) :: omega
    DOUBLE PRECISION, INTENT(IN) :: z !atomic number of the TARGET!
    DOUBLE PRECISION, INTENT(IN) :: i !ionization energy in eV

    ! Data dictionary: declare local vals
    DOUBLE PRECISION             :: numerator
    DOUBLE PRECISION             :: denominator

    ! Initialize values
    green_mcneal = 0
    numerator    = 0
    denominator  = 0

    ! Perform calculation
    numerator    = ((z*a)**omega)*((energy)**nu)
    numerator    = numerator*1.E-16
    denominator  = (j**(omega + nu)) + (energy**(omega+nu))
    green_mcneal = numerator/denominator
    RETURN
  END FUNCTION green_mcneal

  FUNCTION a_gs(energy,k,k_b,j,j_b,j_c)
    !
    !  Purpose:
    !    To calculate the A(E) prefactor for the Green-Sawada (1973) electron
    !  impact cross-sections.
    !
    !! A_GS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    ! Data dictionary: declare calling parameters
    DOUBLE PRECISION             :: a_gs
    DOUBLE PRECISION, INTENT(IN) :: energy
    DOUBLE PRECISION, INTENT(IN) :: k
    DOUBLE PRECISION, INTENT(IN) :: k_b
    DOUBLE PRECISION, INTENT(IN) :: j
    DOUBLE PRECISION, INTENT(IN) :: j_b
    DOUBLE PRECISION, INTENT(IN) :: j_c

    ! Data dictionary: declare local vals
    DOUBLE PRECISION             :: factor1
    DOUBLE PRECISION             :: factor2

    ! Initialize values
    a_gs = 0
    factor1 = 0
    factor2 = 0

    ! Carry out computation
    factor1 = (k/energy + k_b)
    factor2 = energy/j + j_b + j_c/energy
    a_gs    = factor1*DLOG(factor2)
    RETURN
  END FUNCTION a_gs

  FUNCTION gamma_gs(energy,gamma_s,gamma_b)
    !
    !  Purpose:
    !    To calculate the Gamma(E) for the Green-Sawada (1973) electron
    !  impact cross-sections.
    !
    !! GAMMA_GS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    ! Data dictionary: declare calling parameters
    DOUBLE PRECISION             :: gamma_gs
    DOUBLE PRECISION, INTENT(IN) :: energy
    DOUBLE PRECISION, INTENT(IN) :: gamma_s
    DOUBLE PRECISION, INTENT(IN) :: gamma_b

    ! Data dictionary: declare local vals
    DOUBLE PRECISION             :: numerator
    DOUBLE PRECISION             :: denominator

    ! Initialize variables
    gamma_gs = 0
    numerator = 0
    denominator = 0

    ! Perform calculation
    numerator   = gamma_s*energy
    denominator = energy + gamma_b
    gamma_gs    = numerator/denominator
    RETURN
  END FUNCTION gamma_gs

  FUNCTION t_0_gs(energy,t_a,t_b,t_s)
    !
    !  Purpose:
    !    To calculate T_0 for the Green-Sawada (1973) electron
    !  impact cross-sections.
    !
    !! T_0_GS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    ! Data dictionary: declare calling parameters
    DOUBLE PRECISION             :: t_0_gs
    DOUBLE PRECISION, INTENT(IN) :: energy
    DOUBLE PRECISION, INTENT(IN) :: t_a
    DOUBLE PRECISION, INTENT(IN) :: t_b
    DOUBLE PRECISION, INTENT(IN) :: t_s

    ! Data dictionary: declare local vals
    DOUBLE PRECISION             :: bracket

    ! Initialize values
    t_0_gs = 0
    bracket = 0

    ! Perform calculation
    bracket = t_a/(energy + t_b)
    t_0_gs  = t_s - bracket
    RETURN
  END FUNCTION t_0_gs

  FUNCTION t_max_gs(energy,i)
    !
    !  Purpose:
    !    To calculate T_max for the Green-Sawada (1973) electron
    !  impact cross-sections.
    !
    !! T_MAX_GS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    ! Data dictionary: declare calling parameters
    DOUBLE PRECISION             :: t_max_gs
    DOUBLE PRECISION, INTENT(IN) :: energy
    DOUBLE PRECISION, INTENT(IN) :: i

    ! Initialize values
    t_max_gs = 0

    ! Perform calculation
    t_max_gs = 0.5*(energy - i)
    RETURN
  END FUNCTION t_max_gs

  FUNCTION green_sawada(a,gamma_fac,t_max,t_0)
    !
    !  Purpose:
    !    To calculate the Green-Sawada (1973) electron
    !  impact cross-sections.
    !
    !! GREEN_SAWADA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    ! Data dictionary: declare calling parameters
    DOUBLE PRECISION             :: green_sawada
    DOUBLE PRECISION, INTENT(IN) :: a
    DOUBLE PRECISION, INTENT(IN) :: gamma_fac
    DOUBLE PRECISION, INTENT(IN) :: t_max
    DOUBLE PRECISION, INTENT(IN) :: t_0

    ! Data dictionary: declare local vals
    DOUBLE PRECISION             :: bracket
    DOUBLE PRECISION             :: paren
    DOUBLE PRECISION             :: insides

    ! Initialize values
    green_sawada = 0
    bracket = 0
    paren = 0
    insides = 0

    ! Perform calculation
    bracket = (t_max - t_0)/gamma_fac
    ! PRINT *, "Inside atan1=",bracket
    paren   = t_0/gamma_fac
    ! PRINT *, "Inside atan2=", paren
    insides = ATAN(bracket) + ATAN(paren)
    ! PRINT *, "insides are ", insides
    green_sawada = 1E-16*a*gamma_fac*insides
    RETURN
  END FUNCTION green_sawada

  FUNCTION ne_mg(gamma_fac,t_max,t_0)
    !
    !  Purpose:
    !   To calculate the scaling factor N(E) for the inelastic H+ stopping
    !  cross-section as taken from Miller & Green 1971.
    !
    !! NE_MG !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    ! Data dictionary: declare calling parameters
    REAL             :: ne_mg
    REAL, INTENT(IN) :: gamma_fac
    REAL, INTENT(IN) :: t_0
    REAL, INTENT(IN) :: t_max

    ! Data dictionary: declare local vals
    REAL             :: part1
    REAL             :: part2

    ! Initialize values
    ne_mg = 0
    part1 = 0
    part2 = 0

    ! Perform calculation
    part1 = ATAN((t_max - t_0)/gamma_fac)
    part2 = ATAN(t_0/gamma_fac)
    ne_mg = gamma_fac*(part1 + part2)
    RETURN
  END FUNCTION ne_mg

  SUBROUTINE magic(eps,b,c2,s2,theta)
    !
    !  The famous Magic formula described in Biersack and Haggmark 1980
    !
    !
    !  Purpose:
    !   To calculate the center-of-mass scattering angle for a given potential
    !  using the formalism developed by Biersack and Haggmark. This algorithm
    !  is the same as is used in TRIM and SRIM by Ziegler.
    !
    !! MAGIC !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    !Data dictionary: Calling Parameters and Output
    DOUBLE PRECISION, INTENT(IN)  :: eps
    DOUBLE PRECISION, INTENT(IN)  :: b
    DOUBLE PRECISION, INTENT(OUT) :: c2,s2
    DOUBLE PRECISION, INTENT(OUT) :: theta

    !Data dictionary: Local Variables
    DOUBLE PRECISION              :: r, rr
    DOUBLE PRECISION              :: ex1,ex2,ex3,ex4
    DOUBLE PRECISION              :: v,v1
    DOUBLE PRECISION              :: fr,fr1
    DOUBLE PRECISION              :: q
    DOUBLE PRECISION              :: roc
    DOUBLE PRECISION              :: sqe
    DOUBLE PRECISION              :: cc,aa,ff
    DOUBLE PRECISION              :: delta
    DOUBLE PRECISION              :: co


    ! Initialize values
    c2    = 0
    s2    = 0
    theta = 0
    r     = 0
    rr    = 0
    ex1   = 0
    ex2   = 0
    ex3   = 0
    ex4   = 0
    v     = 0
    v1    = 0
    fr    = 0
    fr1   = 0
    q     = 0
    roc   = 0
    sqe   = 0
    cc    = 0
    aa    = 0
    ff    = 0
    delta = 0
    co    = 0

    ! Perform calculation
    r     = b
    rr    = -2.7*DLOG(eps*b)
    IF ( rr .LT. b ) GOTO 1980
    rr    = -2.7*DLOG(eps*b)
    IF ( rr .LT. b ) GOTO 1980
    r     = rr
1980 ex1   =  0.18175*EXP(-3.1998*r)
    ex2   =  0.50986*EXP(-0.94229*r)
    ex3   =  0.28022*EXP(-0.4029*r)
    ex4   = 0.028171*EXP(-0.20162*r)
    v     = (ex1 + ex2 + ex3 + ex4)/r
    v1    = -(v+3.1998*ex1+0.94229*ex2+0.4029*ex3+0.20162*ex4)/r
    fr    = b*b/r+v*r/eps-r
    fr1   = -b*b/(r*r)+(v+v1*r)/eps-1.0
    q     = fr/fr1
    r     = r-q
    IF ( ABS(q/r) .GT. 0.001) GOTO 1980
    roc   = -2.0*(eps-v)/v1
    sqe   = SQRT(eps)
    cc    = (0.011615+sqe)/(0.0071222+sqe)
    aa    = 2.0*eps*(1.0+(0.99229/sqe))*b**cc
    ff    = (SQRT(aa**2+1.0)-aa)*((9.3066+eps)/(14.813+eps))
    delta = (r-b)*aa*ff/(ff+1.0)
    co    = (b+delta+roc)/(r+roc)
    c2    = co*co
    s2    = 1.0-c2
    theta = 2.0*ACOS(co)
    RETURN
  END SUBROUTINE magic

  FUNCTION b_magic(rn,a,rho)
    !
    !  Purpose:
    !   To calculate the reduced impact parameter based on a uniform random
    !  number and a screening length.
    !
    !! B_MAGIC !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    !Data dictionary: Calling parameters
    DOUBLE PRECISION, INTENT(IN) :: rn
    DOUBLE PRECISION, INTENT(IN) :: rho
    DOUBLE PRECISION, INTENT(IN) :: a
    DOUBLE PRECISION             :: b_magic

    !Data dictionary: Local variables
    DOUBLE PRECISION             :: p


    ! Initialize values
    b_magic = 0
    p       = 0

    ! Perform calculation
    p = SQRT(rn/(3.14159*(rho**(2./3.))))
    b_magic = p/a
    RETURN
  END FUNCTION b_magic

  FUNCTION au(z1,z2)
    !
    !  Purpose:
    !   To calculate the screening length for a two-particle interaction.
    !
    !! AU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    !Data dictionary: Calling parameters
    DOUBLE PRECISION, INTENT(IN) :: z1,z2
    DOUBLE PRECISION             :: au


    ! Initialize values
    au = 0

    ! Perform calculation
    au = (0.8853*A0)/(z1**0.23 + z2**0.23)
    RETURN
  END FUNCTION au

  FUNCTION eps(en,z1,z2,m1,m2,au)
    !
    !  Purpose:
    !   To calculate the Lindhard-Scharff-Sigmund reduced energy.
    !
    !! EPS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Note: All units must be Gaussian-CGS
    IMPLICIT NONE

    !Data dictionary: Calling parameters
    DOUBLE PRECISION, POINTER    :: en
    DOUBLE PRECISION, INTENT(IN) :: z1,z2
    DOUBLE PRECISION, INTENT(IN) :: m1,m2
    DOUBLE PRECISION, INTENT(IN) :: au
    DOUBLE PRECISION             :: eps

    !Data dictionary: Local variables
    DOUBLE PRECISION             :: fac1,fac2,fac3

    ! Initialize values
    eps = 0
    fac1 = 0
    fac2 = 0
    fac3 = 0

    ! Perform calculation
    fac1 = au/ECHARG2
    fac2 = m2/(m1+m2)
    fac3 = 1./(z1+z2)
    eps = en*fac1*fac2*fac3
  END FUNCTION eps

  FUNCTION mass_fac(m1,m2)
    !
    !  Purpose:
    !   To calculate the mass factor, denoted /gamma in most equations
    !
    !! MASS_FAC !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IMPLICIT NONE

    !Data dictionary: Calling parameters
    DOUBLE PRECISION, INTENT(IN) :: m1,m2
    DOUBLE PRECISION             :: mass_fac


    ! Initialize values
    mass_fac = 0

    ! Perform calculation
    mass_fac = (4.*m1*m2)/((m1+m2)**2)
    RETURN
  END FUNCTION mass_fac

  FUNCTION t_coll(e,mf,s2)
    !
    !  Purpose:
    !   To calculate the energy transferred in an elastic, i.e. nuclear,
    !  collision between a target species and an incoming ion. This value
    !  is calculated based on the output of the Magic Formula, from which
    !  one can calculate s2, which is the sin^2(theta/2)
    !
    !! T_COLL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IMPLICIT NONE

    !Data dictionary: Calling parameters
    DOUBLE PRECISION, POINTER    :: e
    DOUBLE PRECISION, INTENT(IN) :: mf
    DOUBLE PRECISION, INTENT(IN) :: s2
    DOUBLE PRECISION             :: t_coll


    ! Initialize values
    t_coll = 0

    ! Perform calculation
    t_coll = mf*e*s2
    RETURN
  END FUNCTION t_coll

  FUNCTION lab_theta(cmtheta,m1,m2)
    !
    !  Purpose:
    !   To calculate the laboratory frame of reference scattering angle based
    !  on the output of the Magic Formula, from which the model calculated the
    !  center-of-mass scattering angle.
    !
    !! LAB_THETA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    !Data dictionary: Calling parameters
    DOUBLE PRECISION, INTENT(IN) :: cmtheta
    DOUBLE PRECISION, INTENT(IN) :: m1,m2
    DOUBLE PRECISION             :: lab_theta

    !Data dictionary: Local variables
    DOUBLE PRECISION             :: insides


    ! Initialize values
    lab_theta = 0
    insides   = 0

    ! Perform calculation
    insides = SIN(cmtheta)/(COS(cmtheta)+(m1/m2))
    lab_theta = ATAN(insides)
    RETURN
  END FUNCTION lab_theta

  FUNCTION sneps(eps)
    !
    !  Purpose:
    !   To calculate the reduced energy elastic stopping cross-section of
    !  nuclear collisions between a target and incoming ion.
    !
    !! SNEPS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    !Data dictionart: Calling Parameters
    DOUBLE PRECISION, INTENT(IN) :: eps
    DOUBLE PRECISION             :: sneps

    !Data dictionary: Local variables
    DOUBLE PRECISION             :: num
    DOUBLE PRECISION             :: den1,den2,den3,den

    ! Initialize values
    sneps = 0
    num   = 0
    den1  = 0
    den2  = 0
    den3  = 0
    den   = 0

    ! Perform calculation
    IF ( eps .LE. 30. ) THEN
       num = DLOG(1.+1.1383*eps)
       den1 = eps
       den2 = 0.01321*(eps**0.21226)
       den3 = 0.19593*(eps**0.5)
       den  = 2*(den1+den2+den3)
       sneps = num/den
    ELSE
       sneps = DLOG(eps)/(2*eps)
    END IF
    RETURN
  END FUNCTION sneps

  FUNCTION sne(sneps,z1,z2,m1,m2)
    !
    !  Purpose:
    !   To calculate the elastic stopping cross-section, S(E,T) of nuclear
    !  collisions between a target and incoming ion.
    !
    !! SNE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IMPLICIT NONE

    !Data dicionary: Calling parameters
    DOUBLE PRECISION, INTENT(IN) :: sneps
    DOUBLE PRECISION, INTENT(IN) :: z1,z2
    DOUBLE PRECISION, INTENT(IN) :: m1,m2
    DOUBLE PRECISION             :: sne

    !Data dictionary: Local variables
    DOUBLE PRECISION             :: num,den
    DOUBLE PRECISION             :: den1,den2

    ! Initialize values
    sne  = 0
    num  = 0
    den  = 0
    den1 = 0
    den2 = 0

    num = (8.462E-15)*z1*z2*m1*sneps
    den1 = m1 + m2
    den2 = z1**0.23 + z2**0.23
    den = den1*den2

    sne = num/den
    RETURN
  END FUNCTION sne

  FUNCTION pelsig(energy,sne,mf)
    !
    !  Purpose:
    !   To calculate the elastic proton cross-section, sigma(E) of nuclear
    !  collisions between a target and incoming ion.
    !
    !  pelsig => Proton ELastic SIGma
    !
    !! PELSIG !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    !Data dictionary: Calling parameters
    DOUBLE PRECISION, INTENT(IN) :: energy
    DOUBLE PRECISION, INTENT(IN) :: sne
    DOUBLE PRECISION, INTENT(IN) :: mf
    DOUBLE PRECISION             :: pelsig

    ! Initialize values
    pelsig = 0

    ! Perform calculation
    pelsig = (2*sne)/(mf*energy)
    RETURN
  END FUNCTION pelsig

  FUNCTION greendutta(energy,w,f,o,a,b)
    !
    !  Purpose:
    !   To calculate the Green & Dutta (1967) excitation cross-sections
    !  for forbidden transitions caused by electron impact.
    !
    !  NB:
    !    For further explanation of the formula used, see Jackman, Garvey, &
    !  Green 1977.
    !
    !  OUTPUT:
    !    The output of this function is an array of cross-sections, the number
    !  of which equals the number of forbidden states included in the struct.
    !
    !! GREENDUTTA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    ! Data dictionary: Calling parameters
    DOUBLE PRECISION, INTENT(IN) :: energy !eV
    DOUBLE PRECISION, INTENT(IN) :: w,f,o,a,b
    DOUBLE PRECISION             :: greendutta

    ! Local variables
    DOUBLE PRECISION             :: fac1,fac2,fac3

    ! Initialize values
    greendutta = 0
    fac1       = 0
    fac2       = 0
    fac3       = 0

    ! Perform calculation
    IF ( energy .LT. w ) THEN
       greendutta  = 0D0
    ELSE
       fac1       = (Q0*f)/(w*w)
       fac2       = (1.-(w/energy)**a)**b
       fac3       = (w/energy)**o
       greendutta = fac1*fac2*fac3
    END IF
    RETURN
  END FUNCTION greendutta

  FUNCTION pjgsigma(energy,f,w,c,a,b)
    !
    !  Purpose:
    !   To calculate the Porter, Jackman, Green (1976) excitation cross-section
    !  for allowed transitions caused by electron impact.
    !
    !  NB:
    !    For further explanation of the formula used, see Jackman, Garvey, &
    !  Green 1977.
    !
    !  OUTPUT:
    !    The output of this function is an array of cross-sections, the number
    !  of which equals the number of allowed states included in the struct.
    !
    !! PJGSIGMA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    ! Calling parameters
    DOUBLE PRECISION, INTENT(IN)                                :: energy
    DOUBLE PRECISION, INTENT(IN)                                :: f,w,c,a,b
    DOUBLE PRECISION                                            :: pjgsigma

    ! Local variables
    DOUBLE PRECISION                                            :: num, den, insides

    ! Initialize values
    pjgsigma = 0
    num      = 0
    den      = 0
    insides  = 0

    ! Perform calculation
    ! PRINT *, "f=", f
    ! PRINT *, "w=", w
    ! PRINT *, "c=", c
    ! PRINT *, "a=", a
    ! PRINT *, "b=", b
    IF ( energy .LT. w ) THEN
       pjgsigma = 0D0
    ELSE
       num      = Q0*f*(1.-(w/energy)**a)**b
       ! PRINT *, " num=", num
       den      = energy*w
       ! PRINT *, "den=", den
       insides  = (4.*energy*c)/w + EBASE
       ! PRINT *, "insides=", insides
       pjgsigma = (num/den)*DLOG(insides)
       ! PRINT *, "pjgsigma=", pjgsigma
    END IF
    RETURN
  END FUNCTION pjgsigma
END MODULE functiondefs
