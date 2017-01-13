MODULE branchmod
  USE parameters
  USE typedefs
  USE functiondefs
CONTAINS
  FUNCTION branching(reactant1,reactant2,pr)
    ! On first call, ultimately populate global PRODS array with three values,
    ! chosen based on the reaction selected. If there is more than one pathway,
    ! choose one based on it's relative rate-coefficient, if available, or equal
    ! weight statistical probability if not.
    USE typedefs
    IMPLICIT NONE
    INTEGER :: n
    INTEGER :: branching
    INTEGER :: reactant1,reactant2,pr ! Indices of the reaction look-up-table
    integer :: tempr1, tempr2
    INTEGER :: prods(3)
    INTEGER :: reactionid(5)
    REAL    :: rnum
    CHARACTER(80) :: varfmt
    TYPE(reaction), POINTER :: re_temp
    TYPE(reaction_info) :: pathways(5)
    DOUBLE PRECISION :: vals(4)
    DOUBLE PRECISION :: prob1, prob2
    DOUBLE PRECISION :: barrier
    double precision :: rand
    double precision :: e_d

    reactionid = 0
    varfmt = "(5I10,ES15.4)"

    firstcall: IF ( pr .EQ. 1 ) THEN
       ! Initialize pathways data structure
       ISBARRIER = .FALSE.
       PRODS     = 0
       vals      = 0
       prob1     = 0
       prob2     = 0
       barrier   = 0
       initpaths: DO n=1,SIZE(pathways,1)
          pathways(n)%nr1 = 0
          pathways(n)%nr2 = 0
          pathways(n)%np1 = 0
          pathways(n)%np2 = 0
          pathways(n)%np3 = 0
          pathways(n)%arrh_alpha = 0
          pathways(n)%arrh_beta = 0
          pathways(n)%arrh_gamma = 0
          pathways(n)%rtype = 0
          pathways(n)%id = 0
       END DO initpaths

       ! 1. Go through list of reactions and find all between reactant1 and r2
       re_temp => RE_HEAD
       n = 0
       SP_PROD_DEST(reactant1,2) = SP_PROD_DEST(reactant1,2) + 1
       SP_PROD_DEST(reactant2,2) = SP_PROD_DEST(reactant2,2) + 1
       getpaths: DO
          tempr1 = re_temp%nr1
          tempr2 = re_temp%nr2
          IF ( ((tempr1 .EQ. reactant1) .AND. (tempr2 .EQ. reactant2)) .or. &
               ((tempr1 .EQ. reactant2) .AND. (tempr2 .EQ. reactant1))) then
             n = n + 1
             pathways(n)%nr1 = re_temp%nr1
             pathways(n)%nr2 = re_temp%nr2
             pathways(n)%np1 = re_temp%np1
             pathways(n)%np2 = re_temp%np2
             pathways(n)%np3 = re_temp%np3
             pathways(n)%arrh_alpha = re_temp%arrh_alpha
             pathways(n)%arrh_beta = re_temp%arrh_beta
             pathways(n)%arrh_gamma = re_temp%arrh_gamma
             pathways(n)%rtype = re_temp%rtype
             pathways(n)%id = re_temp%id
          END IF
          IF ( .NOT. ASSOCIATED(re_temp%next)) EXIT
          re_temp => re_temp%next
       END DO getpaths

       onepath: IF ( n .EQ. 1 ) THEN
          ! 2. If there is just one possible reaction, proceed, else, calculate rate coeff.
          PRODS(1) = pathways(1)%np1
          PRODS(2) = pathways(1)%np2
          PRODS(3) = pathways(1)%np3
          barrier  = pathways(1)%arrh_gamma
       ElSE
          IF ( pathways(1)%rtype .NE. pathways(2)%rtype ) THEN
             ! If pathways rtype is not the same, we can't compare them directly
             PRINT *,"Product pathways are not of the same type"
             CALL EXIT()
          END IF

          ! 3. Based on reaction type, calculate rate coeff. and store that in data-structure
          SELECT CASE (pathways(1)%rtype)
          CASE(2)
             DO n=1,SIZE(vals,1)
                vals(n) = arrhenius(&
                     pathways(n)%arrh_alpha,&
                     pathways(n)%arrh_beta,&
                     pathways(n)%arrh_gamma)
             END DO
          CASE DEFAULT
             DO n=1,SIZE(vals,1)
                vals(n) = pathways(n)%arrh_alpha
             END DO
          END SELECT

          CALL selectbranch(vals,n)

          PRODS(1) = pathways(n)%np1
          PRODS(2) = pathways(n)%np2
          PRODS(3) = pathways(n)%np3
          barrier  = pathways(n)%arrh_gamma
       END IF onepath

       ! 6. If there is a barrier to the selected pathway, engage the competition mechanism
       !    to determine if the reaction proceeds or the hopping species diffuses away
       competition: IF ( barrier .GT. 0.0d0) THEN
          e_d = E_BULK*EN_LIST(reactant1)
          prob1 = e_d/(e_d+barrier)
          prob2 = barrier/(e_d+barrier)
          CALL RANDOM_NUMBER(rand)
          IF ( rand .LE. prob1 ) THEN
             PRODS(1) = 0
             PRODS(2) = 0
             PRODS(3) = 0
             n = 0
             ISBARRIER = .TRUE.
          END IF
       END IF competition

       ! On first call: bump the count for that particular reaction for later
       ! time-dependent analysis
       IF ( n .gt. 0 ) THEN
          if ( pathways(n)%id .GT. 0 ) call bumpreaction(pathways(n)%id,RE_HEAD)
       END IF
       IF ( debug .eqv. .true. ) print *, "In branchmod:",pathways(n)
    END IF firstcall

    branching = PRODS(pr)
    SP_PROD_DEST(branching,1) = SP_PROD_DEST(branching,1) + 1
    RETURN
  END FUNCTION branching

  RECURSIVE SUBROUTINE bumpreaction(id,node)
    !
    ! Purpose:
    !   This is a subroutine that takes some reaction id number and increments its
    !  counter for later analytics.
    !! BUMPREACTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE
    INTEGER            :: id
    type(reaction) :: node

    IF ( id .EQ. node%id ) THEN
       node%count = node%count + 1
       RETURN
    ELSE
       IF ( ASSOCIATED(node%next)) THEN
          CALL bumpreaction(id,node%next)
       ELSE
          PRINT *, id
          PRINT *, node%id
          PRINT *, "ERROR! No match in bumpreaction!"
          CALL EXIT()
       END IF
    END IF
    RETURN
  END SUBROUTINE bumpreaction

  SUBROUTINE selectbranch(kvals,ndex)
    ! This subroutine randomly selects a branching pathway
    ! based on the relative magnitude of the rate coefficient
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: kvals(4)
    INTEGER, INTENT(OUT) :: ndex
    INTEGER :: n
    DOUBLE PRECISION :: sumval
    DOUBLE PRECISION :: tempkvals(4)
    DOUBLE PRECISION :: rand

    tempkvals = kvals
    sumval = SUM(tempkvals)
    tempkvals = tempkvals/sumval

    CALL RANDOM_NUMBER(rand)

    DO n=1,SIZE(kvals,1)
       IF ( n .EQ. 1 ) THEN
          tempkvals(n) = tempkvals(n)
       ELSE
          tempkvals(n) = (tempkvals(n) + tempkvals(n-1))
       END IF
       ! Check to see if the  index is the lucky number
       IF ( rand .LE. tempkvals(n)) THEN
          EXIT
       END IF
    END DO

    ndex = n
    IF ( DEBUG .EQV. .TRUE. ) PRINT *, "The lucky number is ",ndex

  END SUBROUTINE selectbranch
END MODULE branchmod
