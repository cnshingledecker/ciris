  SUBROUTINE lookup( string, lines, array, n_line)
  !
  ! Purpose:
  !   This is a subroutine that compares a string value to values
  !  in a list and gives the index of a matching result and an
  !  error if there is no match.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! LOOKUP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IMPLICIT NONE

    !***************** 
    ! Input and output
    !*****************

    CHARACTER(*)     , INTENT(IN)                     :: string
    INTEGER          , INTENT(IN)                     :: lines
    CHARACTER(len=10)             , DIMENSION(lines)  :: array
    INTEGER          , INTENT(OUT)                    :: n_line

    !****************
    ! Local variables
    !****************

    INTEGER                                          :: i
    CHARACTER(len=10), DIMENSION(1)                  :: string_arr

    ! Go through the array and compare the supplied string with the 
    ! strings in the array
    n_line = 0
    string_arr = (/ string /)
    DO i=1,lines
      IF ( TRIM(string_arr(1)) .EQ. TRIM(array(i)) )THEN
        n_line=i
      END IF
    END DO
  END SUBROUTINE lookup

