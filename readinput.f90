MODULE readinput
  USE typedefs
  USE parameters
  USE subroutines
  IMPLICIT NONE
CONTAINS
  SUBROUTINE buildnetwork()
    IMPLICIT NONE
    ! Input and output
    TYPE(species) , POINTER          :: sp_head, sp_tail, sp_temp
    TYPE(reaction), POINTER          :: re_tail
    INTEGER                          :: ierror1, ierror2
    INTEGER                          :: n
    INTEGER                          :: i_count
    INTEGER                          :: rtype
    CHARACTER(LEN=10)                :: tempName,r1,r2,p1,p2,p3
    CHARACTER(LEN=10)                :: line
    DOUBLE PRECISION                 :: e_d
    DOUBLE PRECISION                 :: arrh_alpha,arrh_beta,arrh_gamma

    ! Open files
    OPEN(UNIT=1,FILE=SPECIES_FILE  ,STATUS='OLD',ACTION='READ',IOSTAT=ierror1)
    OPEN(UNIT=2,FILE=REACTIONS_FILE,STATUS='OLD',ACTION='READ',IOSTAT=ierror2)

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
             READ(1,*,IOSTAT=ierror1) tempName,e_d
             addspecies: IF ( .NOT. ASSOCIATED(sp_head)) THEN
                ALLOCATE(sp_head)
                sp_tail => sp_head
                NULLIFY(sp_tail%next)
             ELSE
                ALLOCATE(sp_tail%next)
                sp_tail => sp_tail%next
                NULLIFY(sp_tail%next)
             END IF addspecies
             sp_tail%name = tempName
             sp_tail%e_d = e_d
             sp_tail%ion = .FALSE.
             sp_tail%special = .FALSE.
             sp_tail%id = NUM_SPECIES
             tempName = TRIM(sp_tail%name)
             ! Tally up the total number of anions, cations, and special species
             specialion: IF ( (tempName(LEN(TRIM(tempName)):LEN(TRIM(tempName))) .EQ. '-') &
                  .OR. (tempName(LEN(TRIM(tempName)):LEN(TRIM(tempName))) .EQ. '+')) THEN
                IF ( tempName(1:1) .NE. '*') THEN
                   sp_tail%ion = .TRUE.
                   i_count = i_count + 1
                END IF
             ELSE IF  (tempName .EQ. "CRP") THEN
                CRPNUM = NUM_SPECIES
                sp_tail%special = .TRUE.
             ELSE IF (tempName .EQ. "e") THEN
                ELECNUM = NUM_SPECIES
                sp_tail%special = .TRUE.
             ELSE IF (tempName .EQ. "*") THEN
                EXCNUM = NUM_SPECIES
                sp_tail%special = .TRUE.
             END IF specialion
          ELSE
             CONTINUE
          END IF
       END DO countspecies

       SPECIAL_LIST = (/ CRPNUM, EXCNUM, ELECNUM /)
       ALLOCATE(IONLIST(i_count))


       sp_temp => sp_head
       n = 0
       popions: DO
          IF ( .NOT. ASSOCIATED(sp_temp%next)) EXIT
          IF ( sp_temp%ion .EQV. .TRUE. ) THEN
             n = n + 1
             IONLIST(n) = sp_temp%id
          END IF
          sp_temp => sp_temp%next
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
             READ(2,*,IOSTAT=ierror2) r1,r2,p1,p2,p3,arrh_alpha,arrh_beta,arrh_gamma,rtype
             addreact: IF ( .NOT. ASSOCIATED(RE_HEAD)) THEN
                ALLOCATE(RE_HEAD)
                re_tail => RE_HEAD
                NULLIFY(re_tail%next)
             ELSE
                ALLOCATE(re_tail%next)
                re_tail => re_tail%next
             END IF addreact
             re_tail%r1 = r1
             re_tail%r2 = r2
             re_tail%p1 = p1
             re_tail%p2 = p2
             re_tail%p3 = p3
             re_tail%arrh_alpha = arrh_alpha
             re_tail%arrh_beta = arrh_beta
             re_tail%arrh_gamma = arrh_gamma
             re_tail%rtype = rtype
             ! Get integer identities of species
             CALL lookup(re_tail%r1,sp_head,re_tail%nr1)
             CALL lookup(re_tail%r2,sp_head,re_tail%nr2)
             CALL lookup(re_tail%p1,sp_head,re_tail%np1)
             IF ( TRIM(re_tail%p2) .NE. "0" ) CALL lookup(re_tail%p2,sp_head,re_tail%np2)
             IF ( TRIM(re_tail%p3) .NE. "0" ) CALL lookup(re_tail%p3,sp_head,re_tail%np3)
          ELSE
             CONTINUE
          END IF
       END DO countreactions
    ELSE fileopen
       IF ( ierror1 .NE. 0 .AND. ierror2 .NE. 0 ) THEN
          PRINT *, 'Could not open any of the input files'
       ELSE IF ( ierror1 .NE. 0 .AND. ierror2 .EQ. 0 ) THEN
          PRINT *, 'Could not open: species_file'
       ELSE
          PRINT *, 'Could not open: reactions_file'
       END IF
    END IF fileopen

    PRINT *, "There are",NUM_SPECIES,"species"
    PRINT *, "There are",NUM_REACTS,"reactions"
    PRINT *, "There are",i_count,"ions"

    ! Go through species linked list and build global array, while de-allocating
    ! linked list memory
    ALLOCATE(EN_LIST(NUM_SPECIES))
    ALLOCATE(SP_LIST(NUM_SPECIES))
    DO
       sp_temp => sp_head
       n = sp_temp%id
       EN_LIST(n) = sp_temp%e_d
       SP_LIST(n) = sp_temp%name
       if ( .not. associated(sp_temp%next)) exit
       sp_tail => sp_temp%next
       NULLIFY(sp_head)
       sp_head => sp_tail
    end do
  END SUBROUTINE buildnetwork
END MODULE readinput
