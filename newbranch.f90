MODULE branchmod
  USE parameters
  USE functiondefs
  USE typedefs
CONTAINS
  FUNCTION branching(r1,r2,pr)
    ! On first call, ultimately populate global PRODS array with three values,
    ! chosen based on the reaction selected. If there is more than one pathway,
    ! choose one based on it's relative rate-coefficient, if available, or equal
    ! weight statistical probability if not.
    IMPLICIT NONE
    INTEGER :: branching
    INTEGER :: r1,r2,pr ! Indices of the reaction look-up-table
    INTEGER :: prods(3)
    INTEGER :: reaction(5)
    REAL    :: rnum
    CHARACTER(80) :: varfmt
    TYPE(reaction), POINTER :: re_temp
    TYPE(reaction_info) :: pathways(2)
    DOUBLE PRECISION :: val1, val2
    DOUBLE PRECISION :: prob1, prob2
    DOUBLE PRECISION :: barrier

    reaction = 0
    varfmt = "(5I10,ES15.4)"

    firstcall: IF ( pr .EQ. 1 ) THEN
       ! Initialize pathways data structure
       PRODS = 0
       val1    = 0
       val2    = 0
       prob1   = 0
       prob2   = 0
       barrier = 0
       initpaths: DO n=1,2
          pathways(n)%nr1 = 0
          pathways(n)%nr2 = 0
          pathways(n)%np1 = 0
          pathways(n)%np2 = 0
          pathways(n)%np3 = 0
          pathways(n)%arrh_alpha = 0
          pathways(n)%arrh_beta = 0
          pathways(n)%arrh_gamma = 0
          pathways(n)%rtype = 0
       END DO initpaths

       ! 1. Go through list of reactions and find all between r1 and r2
       re_temp => RE_HEAD
       n = 0
       getpaths: DO
          IF ( .NOT. ASSOCIATED(re_temp%next)) EXIT
          IF ( ((re_temp%nr1 .EQ. r1) .AND. (re_temp%nr2 .EQ. r2)) .OR. &
               ((re_temp%nr1 .EQ. r2) .AND. (re_temp%nr2 .EQ. r1)) ) THEN
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
          END IF
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
          CASE(1)
             PRINT *, "Case",pathways(1)%rtype,"is not implemented!"
             CALL EXIT()
          CASE(2)
             val1 = arrhenius(&
                  pathways(1)%arrh_alpha,&
                  pathways(1)%arrh_beta,&
                  pathways(1)%arrh_gamma)
             val2 = arrhenius(&
                  pathways(2)%arrh_alpha,&
                  pathways(2)%arrh_beta,&
                  pathways(2)%arrh_gamma)
          CASE(3)
             PRINT *, "Case",pathways(1)%rtype,"is not implemented!"
             CALL EXIT()
          CASE(4)
             PRINT *, "Case",pathways(1)%rtype,"is not implemented!"
             CALL EXIT()
          CASE(5)
             PRINT *, "Case",pathways(1)%rtype,"is not implemented!"
             CALL EXIT()
          CASE(6)
             val1 = pathways(1)%arrh_alpha
             val2 = pathways(2)%arrh_alpha
          CASE(7)
             PRINT *, "Case",pathways(1)%rtype,"is not implemented!"
             CALL EXIT()
          CASE(8)
             PRINT *, "Case",pathways(1)%rtype,"is not implemented!"
             CALL EXIT()
          END SELECT

          ! 4. Call a random number and select pathway based on result
          prob1 = val1/(val1+val2)
          prob2 = val2/(val1+val2)
          CALL RANDOM_NUMBER(rand)

          ! 5. Populate PRODS array with result
          IF ( rand .LE. prob1 ) THEN
             PRODS(1) = pathways(1)%p1
             PRODS(2) = pathways(1)%p2
             PRODS(3) = pathways(1)%p3
             barrier  = pathways(1)%arrh_gamma
          ELSE
             PRODS(1) = pathways(2)%p1
             PRODS(2) = pathways(2)%p2
             PRODS(3) = pathways(2)%p3
             barrier  = pathways(2)
          END IF
       END IF onepath

       ! 6. If there is a barrier to the selected pathway, engage the competition mechanism
       !    to determine if the reaction proceeds or the hopping species diffuses away
       e_d = E_BULK*EN_LIST(r1)
       prob1 = e_d/(e_d+barrier)
       prob2 = barrier/(e_d+barrier)
       CALL RANDOM_NUMBER(rand)
       IF ( rand .LE. prob1 ) THEN
          PRODS(1) = 0
          PRODS(2) = 0
          PRODS(3) = 0
       END IF
    END IF firstcall

    reaction(1) = r1
    reaction(2) = r2
    reaction(3:5) = PRODS
    ! If checking reactions, print out reactants and products
    o3check: IF ( (O3_ANALYTICS .EQV. .TRUE.) ) THEN
       IF ( ANY(reaction .EQ. O3NUM) ) THEN
          IF ( SUM(PRODS) .GT. 0 ) THEN
             OPEN(FILE="ozone_reactions.wsv", &
                  UNIT=REACTIONS_UNIT_NUM, &
                  STATUS="UNKNOWN", &
                  POSITION="APPEND")
             IF ( r1 .LE. r2 ) THEN
                WRITE(REACTIONS_UNIT_NUM,varfmt) r1,r2,PRODS(1),PRODS(2),PRODS(3),FLUENCE
             ELSE
                WRITE(REACTIONS_UNIT_NUM,varfmt) r2,r1,PRODS(1),PRODS(2),PRODS(3),FLUENCE
             END IF
             CLOSE(REACTIONS_UNIT_NUM)
          END IF
       END IF
    END IF o3check

    branching = PRODS(pr)
    RETURN
  END FUNCTION branching
END MODULE branchmod
