MODULE newmod
  IMPLICIT NONE

  TYPE species
     CHARACTER(len=10) :: name
     INTEGER :: id
     DOUBLE PRECISION :: e_d
     LOGICAL :: ion
     LOGICAL :: special
     TYPE(species), POINTER :: next
  END TYPE species

  TYPE reaction
     CHARACTER(len=10) :: r1
     CHARACTER(len=10) :: r2
     CHARACTER(len=10) :: p1
     CHARACTER(len=10) :: p2
     CHARACTER(len=10) :: p3
     INTEGER :: nr1
     INTEGER :: nr2
     INTEGER :: np1
     INTEGER :: np2
     INTEGER :: np3
     DOUBLE PRECISION :: arrh_alpha
     DOUBLE PRECISION :: arrh_beta
     DOUBLE PRECISION :: arrh_gamma
     INTEGER :: rtype
     TYPE(reaction), POINTER :: next
  END TYPE reaction

CONTAINS
  RECURSIVE SUBROUTINE add_species(root,temp)
    IMPLICIT NONE
    TYPE(species), POINTER :: root, temp

    IF ( .NOT. ASSOCIATED(root) ) THEN
       root => temp
    ElSE
       IF ( .NOT. ASSOCIATED(root%next)) THEN
          root%next => temp
       ELSE
          CALL add_species(root%next,temp)
       END IF
    END IF
  END SUBROUTINE add_species

  RECURSIVE SUBROUTINE add_reaction(root,temp)
    IMPLICIT NONE
    TYPE(reaction), POINTER :: root, temp

    IF ( .NOT. ASSOCIATED(root) ) THEN
       root => temp
    ElSE
       IF ( .NOT. ASSOCIATED(root%next)) THEN
          root%next => temp
       ELSE
          CALL add_reaction(root%next,temp)
       END IF
    END IF
  END SUBROUTINE add_reaction

  RECURSIVE SUBROUTINE lookup(name,node,id)
    !
    ! Purpose:
    !   This is a subroutine that compares a string value to values
    !  in a list and gives the index of a matching result and an
    !  error if there is no match.
    !
    !! LOOKUP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE
    INTEGER            :: id
    CHARACTER(len=10)  :: name
    TYPE(species) :: node

    IF ( TRIM(name) .EQ. TRIM(node%name) ) THEN
       id = node%id
       RETURN
    ELSE
       IF ( ASSOCIATED(node%next)) THEN
          CALL lookup(name,node%next,id)
       ELSE
          id = -1
          PRINT *, "ERROR! No match!"
          CALL EXIT()
       END IF
    END IF
    RETURN
  END SUBROUTINE lookup
END MODULE newmod
