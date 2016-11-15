MODULE branchmod
  USE parameters
CONTAINS
  FUNCTION branching(r1,r2,pr)
    IMPLICIT NONE
    INTEGER :: branching
    INTEGER :: r1,r2,pr ! Indices of the reaction look-up-table
    INTEGER :: prods(3)
    INTEGER :: reaction(5)
    REAL    :: rnum

    firstcall: IF ( pr .EQ. 1 ) THEN
       reaction = 0
       CALL RANDOM_NUMBER(rnum)
       prods = REACT_CUBE(r1,r2,:)

       !****************************************************************************
       !**** BRANCHING RATIOS ******************************************************
       !****************************************************************************
       isdis: IF ( (r1 .EQ. EXCNUM) .OR. (r2 .EQ. EXCNUM) ) THEN
          IF ( (r1 .EQ. O2NUM) .OR. (r2 .EQ. O2NUM) ) THEN
             IF ( rnum .LE. O2_DISPROB ) THEN
                prods = (/ ONUM, ONUM, 0 /)
             ELSE
                prods = (/ 0, 0, 0 /)
             END IF
          ELSE IF ( (r1 .EQ. O3NUM) .OR. (r2 .EQ. O3NUM) ) THEN
             IF ( rnum .LE. O3_DISPROB) THEN
                prods = (/ O2NUM, ONUM, 0 /)
             ELSE
                prods = (/ 0, 0, 0 /)
             END IF
          END IF
       ELSE
          branchcond: IF ( r1 .EQ. ONUM .AND. r2 .EQ. O2NUM .OR. r1 .EQ. O2NUM .AND. r2 .EQ. ONUM ) THEN
             IF ( rnum .LE. O3_DISPROB ) THEN
                ! O + O2 -> O3* -> O + O2
                prods = (/ O2NUM, ONUM, 0 /)
             ELSE
                ! O + O2 -> O3
                prods = (/O3NUM, 0, 0 /)
             END IF
          ELSE IF ( r1 .EQ. 2 .AND. r2 .EQ. 3 .OR. r1 .EQ. 3 .AND. r2 .EQ. 2 ) THEN
             IF ( rnum .LE. O2_DISPROB ) THEN
                ! O2+ + O2- -> O2* + O2 -> O + O + O2
                prods = (/ ONUM, ONUM, O2NUM /)
             ELSE
                ! O2+ + O2- -> O2 + O2
                prods = (/ O2NUM, O2NUM, 0 /)
             END IF
          ELSE IF ( r1 .EQ. 8 .AND. r2 .EQ. 6 .OR. r1 .EQ. 6 .AND. r2 .EQ. 8 ) THEN
             IF ( rnum .LE. O2_DISPROB ) THEN
                ! O- + O3+ -> O2* + O2 -> O + O + O2
                prods = (/ ONUM, ONUM, O2NUM /)
             ELSE
                ! O- + O3+ -> O2 + O2
                prods = (/ O2NUM, O2NUM, 0 /)
             END IF
          ELSE IF ( r1 .EQ. 9 .AND. r2 .EQ. 5 .OR. r1 .EQ. 5 .AND. r2 .EQ. 9 ) THEN
             IF ( rnum .LE. O2_DISPROB ) THEN
                ! O3- + O+ -> O2* + O2 -> O + O + O2
                prods = (/ ONUM, ONUM, O2NUM /)
             ELSE
                ! O3- + O+ -> O2 + O2
                prods = (/ O2NUM, O2NUM, 0 /)
             END IF
          ELSE IF ( r1 .EQ. 5 .AND. r2 .EQ. 6 .OR. r1 .EQ. 6 .AND. r2 .EQ. 5 ) THEN
             IF ( rnum .LE. O2_DISPROB ) THEN
                ! O+ + O- -> O2* -> O + O
                prods = (/ ONUM, ONUM, 0 /)
             ELSE
                ! O+ + O- -> O2
                prods = (/ O2NUM, 0, 0 /)
             END IF
          ELSE IF ( r1 .EQ. 8 .AND. r2 .EQ. 3 .OR. r1 .EQ. 3 .AND. r2 .EQ. 8 ) THEN
             IF ( rnum .LE. O3_DISPROB ) THEN
                ! O3+ + O2- -> O3* + O2 -> O + O2 + O2
                prods = (/ ONUM, O2NUM, O2NUM /)
             ELSE
                ! O3+ + O2- -> O3 + O + O
                prods = (/ O3NUM, O2NUM, 0 /)
             END IF
          ELSE IF ( r1 .EQ. 9 .AND. r2 .EQ. 2 .OR. r1 .EQ. 2 .AND. r2 .EQ. 9 ) THEN
             IF ( rnum .LE. O3_DISPROB ) THEN
                ! O3- + O2+ -> O3* + O2 -> O + O2 + O2
                prods = (/ ONUM, O2NUM, O2NUM /)
             ELSE
                ! O3- + O2+ -> O3 + O + O
                prods = (/ O3NUM, O2NUM, 0 /)
             END IF
          END IF branchcond

          ! If none of the conditionals proc, the react cube should get the original values
          REACT_CUBE(r1,r2,:) = prods
          reaction(1) = r1
          reaction(2) = r2
          reaction(3:5) = prods

          ! If checking reactions, print out reactants and products
          o3check: IF ( (O3_ANALYTICS .EQV. .TRUE.) ) THEN
             IF ( ANY(reaction .EQ. O3NUM) ) THEN
                IF ( SUM(prods) .GT. 0 ) THEN
                   OPEN(FILE="ozone_reactions.csv", UNIT=REACTIONS_UNIT_NUM, STATUS="UNKNOWN", POSITION="APPEND")
                   WRITE(REACTIONS_UNIT_NUM,*) r1,',',r2,',',prods(1),',',prods(2),',',prods(3),',',TIME*CR_FLUX
                   CLOSE(REACTIONS_UNIT_NUM)
                END IF
             END IF
          END IF o3check
       END IF isdis
    END IF firstcall

    branching = REACT_CUBE(r1,r2,pr)

    dbug: IF ( DEBUG .EQV. .TRUE. ) THEN
       PRINT *, REACT_CUBE(r1,r2,:)
       PRINT *, "BRANCHING=",branching
    END IF dbug
    RETURN
  END FUNCTION branching
END MODULE branchmod
