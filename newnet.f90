PROGRAM newnet
  IMPLICIT NONE

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

  RECURSIVE FUNCTION lookup(name,node)
    !
    ! Purpose:
    !   This is a subroutine that compares a string value to values
    !  in a list and gives the index of a matching result and an
    !  error if there is no match.
    !
    !! LOOKUP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE
    INTEGER            :: lookup
    CHARACTER(len=10)  :: name
    TYPE(species) :: node

    IF ( TRIM(name) .EQ. TRIM(node%name) ) THEN
       lookup = node%id
       RETURN
    ELSE
       IF ( ASSOCIATED(node%next)) THEN
          CALL lookup(name,node%next)
       ELSE
          lookup = -1
          PRINT *, "ERROR! No match!"
          CALL EXIT()
       END IF
    END IF
    RETURN
  END FUNCTION lookup


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

  ! Input and output
  TYPE(species) , POINTER          :: sp_root, sp_temp, sp_temp2
  TYPE(reaction), POINTER          :: re_root, re_temp
  INTEGER                          :: ierror1, ierror2
  INTEGER                          :: n
  INTEGER                          :: i_count
  CHARACTER(LEN=10)                :: tempName
  CHARACTER(LEN=10)                :: line

  INTEGER :: CRPNUM, EXCNUM, ELECNUM, SPECIAL_LIST(3)
  INTEGER, ALLOCATABLE :: IONLIST(:)

  ! Open files
  OPEN(UNIT=1,FILE="species.dat"  ,STATUS='OLD',ACTION='READ',IOSTAT=ierror1)
  OPEN(UNIT=2,FILE="reactions.dat",STATUS='OLD',ACTION='READ',IOSTAT=ierror2)
  ALLOCATE(sp_root,sp_temp,sp_temp2)
  ALLOCATE(re_root,re_temp)

  fileopen: IF ( ierror1 .EQ. 0 .AND. ierror2 .EQ. 0 ) THEN

     ! Initialize local values
     i_count = 0
     line    = '0'

     ! Read the contents of the species file and create the SP_LIST and energy_list
     countspecies: DO
        READ (1,*,IOSTAT=ierror1) line
        IF ( ierror1 .NE. 0 ) EXIT
        IF ( line(1:1) .NE. "!" ) THEN
           NUM_SPECIES = NUM_SPECIES + 1
           BACKSPACE(UNIT=1,IOSTAT=ierror1)
           ! Re-read line and save values to appropriate data-structure elements
           READ(1,*,IOSTAT=ierror1) sp_temp%name, sp_temp%e_d
           sp_temp%ion = .FALSE.
           sp_temp%species = .FALSE.
           sp_temp%id = NUM_SPECIES
           CALL add_species(sp_root,sp_temp)
           tempName = TRIM(sp_temp%name)
           ! Tally up the total number of anions, cations, and special species
           specialion: IF ( (tempName(LEN(TRIM(tempName)):LEN(TRIM(tempName))) .EQ. '-') &
                .OR. (tempName(LEN(TRIM(tempName)):LEN(TRIM(tempName))) .EQ. '+')) THEN
              IF ( tempName(1:1) .NE. '*') THEN
                 sp_temp%ion = .TRUE.
                 i_count = i_count + 1
              END IF
           ELSE IF  (tempName = "CRP") THEN
              CRPNUM = NUM_SPECIES
              sp_temp%special = .TRUE.
           ELSE IF (tempName .EQ. "e") THEN
              ELECNUM = NUM_SPECIES
              sp_temp%special = .TRUE.
           ELSE IF (tempName .EQ. "*") THEN
              EXCNUM = NUM_SPECIES
              sp_temp%special = .TRUE.
           END IF specialion
        ELSE
           CONTINUE
        END IF
     END DO countspecies

     SPECIAL_LIST = (/ CRPNUM, EXCNUM, ELECNUM /)
     ALLOCATE(IONLIST(i_count))


     sp_temp => sp_root
     n = 0
     popions: DO
        IF ( sp_temp%ion .EQV. .TRUE. ) THEN
           n = n + 1
           IONLIST(n) = id
           IF ( ASSOCIATED(sp_temp%next) ) THEN
              sp_temp2 => sp_temp%next
              sp_temp => sp_temp2
           ELSE
              PRINT *, "ERROR! END OF THE LINE FOR SPECIES LIST"
              CALL EXIT()
           END IF
        END IF
        IF ( n .EQ. i_count ) EXIT
     END DO popions

     ! Read the contents of the reactions file and make the reactionCube
     countreactions: DO
        READ(2,*,IOSTAT=ierror2) line
        ! Exit upon reaching the end of the file
        IF ( ierror2 .NE. 0 ) EXIT
        IF ( line(1:1) .NE. "!" ) THEN
           NUM_REACTS = NUM_REACTS + 1
           ! Rewind one line and re-read
           BACKSPACE (UNIT=2,IOSTAT=ierror2)
           ! Read in variables
           READ(2,*,IOSTAT=ierror2) &
                re_temp%r1, &
                re_temp%r2, &
                re_temp%p1, &
                re_temp%p2, &
                re_temp%p3, &
                re_temp%arrh_alpha, &
                re_temp%arrh_beta , &
                re_temp%arrh_gamma, &
                re_temp%rtype
           ! Get integer identities of species
           re_temp%nr1 =  lookup(re_temp%r1,sp_root)
           re_temp%nr2 =  lookup(re_temp%r2,sp_root)
           re_temp%np1 =  lookup(re_temp%p1,sp_root)
           IF ( TRIM(re_temp%p2) .NE. "0" ) re_temp%np2 = lookup(re_temp%p2,sp_root)
           IF ( TRIM(re_temp%p3) .NE. "0" ) re_temp%np3 = lookup(re_temp%p3,sp_root)
        ELSE
           CONTINUE
        END IF
     END DO countreactions
  ELSE fileopen
     IF ( ierror1 .NE. 0 .AND. ierror2 .NE. 0 ) THEN
        PRINT *, 'Could not open any of the input files'
     ELSE IF ( ierror1 .NE. 0 .AND. ierror2 .EQ. 0 ) THEN
        PRINT *, 'Could not open: ', species_file
     ELSE
        PRINT *, 'Could not open: ', reactions_file
     END IF
  END IF fileopen
END PROGRAM newnet

