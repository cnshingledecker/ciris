MODULE branchmod
  USE parameters
CONTAINS
  FUNCTION branching(r1,r2,pr)
    IMPLICIT NONE
    INTEGER :: branching
    INTEGER :: r1,r2,pr ! Indices of the reaction look-up-table
    INTEGER :: prods(3)
    REAL    :: rnum

    IF ( pr .EQ. 1 ) THEN   
       CALL RANDOM_NUMBER(rnum)
       prods = REACT_CUBE(r1,r2,:)

       !****************************************************************************
       !**** BRANCHING RATIOS ******************************************************
       !****************************************************************************
       IF ( r1 .EQ. ONUM .AND. r2 .EQ. O2NUM .OR. r1 .EQ. O2NUM .AND. r2 .EQ. ONUM ) THEN
          ! O + O2 -> O + O2
          IF ( rnum .LE. O_O2_BRANCHING ) THEN
             prods = (/ ONUM, O2NUM, 0 /)
          END IF
       ELSE IF ( r1 .EQ. 2 .AND. r2 .EQ. 3 .OR. r1 .EQ. 3 .AND. r2 .EQ. 2 ) THEN
          ! O2+ + O2- -> O + O + O2
          IF ( rnum .LE. O2_ION_BRANCHING ) THEN
             prods = (/ ONUM, ONUM, O2NUM /)
          END IF
       ELSE IF ( r1 .EQ. 8 .AND. r2 .EQ. 6 .OR. r1 .EQ. 6 .AND. r2 .EQ. 8 ) THEN
          ! O- + O3+ -> O + O + O2
          IF ( rnum .LE. O3_O_ION_BRANCHING ) THEN
             prods = (/ ONUM, ONUM, O2NUM /)
          END IF
       ELSE IF ( r1 .EQ. 5 .AND. r2 .EQ. 6 .OR. r1 .EQ. 6 .AND. r2 .EQ. 5 ) THEN
          ! O+ + O- -> O2
          IF ( rnum .LE. O_ION_BRANCHING ) THEN
             prods = (/ O2NUM, 0, 0 /)
          END IF
       ELSE IF ( r1 .EQ. 9 .AND. r2 .EQ. 5 .OR. r1 .EQ. 5 .AND. r2 .EQ. 9 ) THEN
          ! O3- + O+ -> O + O + O2
          IF ( rnum .LE. O3_O_ION_BRANCHING ) THEN
             prods = (/ ONUM, ONUM, O2NUM /)
          END IF
       ELSE IF ( r1 .EQ. 8 .AND. r2 .EQ. 3 .OR. r1 .EQ. 3 .AND. r2 .EQ. 8 ) THEN
          ! O3+ + O2- -> O3 + O2
          IF ( rnum .LE. O3_O2_ION_BRANCHING ) THEN
             prods = (/ O3NUM, O2NUM, 0 /)
          END IF
       ELSE IF ( r1 .EQ. 9 .AND. r2 .EQ. 2 .OR. r1 .EQ. 2 .AND. r2 .EQ. 9 ) THEN
          ! O3- + O2+ -> O3 + O2
          IF ( rnum .LE. O3_O2_ION_BRANCHING ) THEN
             prods = (/ O3NUM, O2NUM, 0 /)
          END IF
       END IF

       ! If none of the conditionals proc, the react cube should get the original values
       REACT_CUBE(r1,r2,:) = prods
    END IF

    branching = REACT_CUBE(r1,r2,pr)

    IF ( DEBUG .EQV. .TRUE. ) PRINT *, "BRANCHING=",branching
    IF ( branching .EQ. 0 ) THEN
       ! Restore original values to react cube
       !****************************************************************************
       !**** ORIGINAL VALUES *******************************************************
       !****************************************************************************
       IF ( r1 .EQ. ONUM .AND. r2 .EQ. O2NUM .OR. r1 .EQ. O2NUM .AND. r2 .EQ. ONUM ) THEN
          ! O + O2 -> O3
          prods = (/ O3NUM, 0, 0 /)
       ELSE IF ( r1 .EQ. 2 .AND. r2 .EQ. 3 .OR. r1 .EQ. 3 .AND. r2 .EQ. 2 ) THEN
          ! O2+ + O2- -> O2 + O2
          prods = (/ O2NUM, O2NUM, 0 /)
       ELSE IF ( r1 .EQ. 8 .AND. r2 .EQ. 6 .OR. r1 .EQ. 6 .AND. r2 .EQ. 8 ) THEN
          ! O- + O3+ -> O2 + O2
          prods = (/ O2NUM, O2NUM, 0 /)
       ELSE IF ( r1 .EQ. 5 .AND. r2 .EQ. 6 .OR. r1 .EQ. 6 .AND. r2 .EQ. 5 ) THEN
          ! O+ + O- -> O + O
          prods = (/ ONUM, ONUM, 0 /)
       ELSE IF ( r1 .EQ. 9 .AND. r2 .EQ. 5 .OR. r1 .EQ. 5 .AND. r2 .EQ. 9 ) THEN
          ! O3- + O+ -> O2 + O2
          prods = (/ O2NUM, O2NUM, 0 /)
       ELSE IF ( r1 .EQ. 8 .AND. r2 .EQ. 3 .OR. r1 .EQ. 3 .AND. r2 .EQ. 8 ) THEN
          ! O3+ + O2- -> O3 + O + O
          prods = (/ O3NUM, ONUM, ONUM /)
       ELSE IF ( r1 .EQ. 9 .AND. r2 .EQ. 2 .OR. r1 .EQ. 2 .AND. r2 .EQ. 9 ) THEN
          ! O3- + O2+ -> O3 + O + O
          prods = (/ O3NUM, ONUM, ONUM /)
       END IF
    END IF

    RETURN
  END FUNCTION branching
END MODULE branchmod
