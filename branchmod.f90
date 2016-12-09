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
    TYPE(reaction_info) :: pathways(2)
    DOUBLE PRECISION :: val1, val2
    DOUBLE PRECISION :: prob1, prob2
    DOUBLE PRECISION :: barrier
    double precision :: rand
    double precision :: e_d

    reactionid = 0
    varfmt = "(5I10,ES15.4)"

    firstcall: IF ( pr .EQ. 1 ) THEN
       if ( ((reactant1 .eq. 9) .and. (reactant2 .eq. 2)) .or. &
            ((reactant1 .eq. 2) .and. (reactant2 .eq. 9))) then
          continue
       end if
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
          pathways(n)%id = 0
       END DO initpaths

       ! 1. Go through list of reactions and find all between reactant1 and r2
       re_temp => RE_HEAD
       n = 0
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
             val1 = pathways(1)%arrh_alpha
             val2 = pathways(2)%arrh_alpha
          END SELECT

          ! 4. Call a random number and select pathway based on result
          prob1 = val1/(val1+val2)
          prob2 = val2/(val1+val2)
          CALL RANDOM_NUMBER(rand)

          ! 5. Populate PRODS array with result
          n = 0
          IF ( rand .LE. prob1 ) THEN
             PRODS(1) = pathways(1)%np1
             PRODS(2) = pathways(1)%np2
             PRODS(3) = pathways(1)%np3
             barrier  = pathways(1)%arrh_gamma
             n = pathways(1)%id
          ELSE
             PRODS(1) = pathways(2)%np1
             PRODS(2) = pathways(2)%np2
             PRODS(3) = pathways(2)%np3
             barrier  = pathways(2)%arrh_gamma
             n = pathways(2)%id
          END IF
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
          END IF
       END IF competition

       if ( n .gt. 0 ) call bumpreaction(n,RE_HEAD)
    END IF firstcall

    branching = PRODS(pr)
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
          id = -1
          PRINT *, "ERROR! No match!"
          CALL EXIT()
       END IF
    END IF
    RETURN
  END SUBROUTINE bumpreaction
END MODULE branchmod
