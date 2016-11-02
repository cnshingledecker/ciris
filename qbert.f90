!-*- f90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  System        :
!  Module        :
!  Object Name   : qbert
!  Revision      : 1
!  Date          : $Date$
!  Author        : Christopher N. Shingledecker
!  Created By    : Christopher N. Shingledecker
!  Created       : Thu Jul 10 10:52:59 2014
!  Last Modified : <140717.1648>
!
!  Description: This is a subroutine to "merge" species and reaction files into a
!               three dimensional array that the main program can use to determine
!               whether two species react or not. It also generates an array that
!               correlates a species' name with its number. Finally, it creates a
!               two dimensional array that lists the Ed,Eb, and Eb2 (Chang & Herbst 2014)
!               values for each surface species.
!
!  Notes:
!
!  History:
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Copyright (c) 2014 Christopher N. Shingledecker.
!
!  All Rights Reserved.
!
! This  document  may  not, in  whole  or in  part, be  copied,  photocopied,
! reproduced,  translated,  or  reduced to any  electronic  medium or machine
! readable form without prior written consent from Christopher N. Shingledecker.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE qbert()
  USE subroutines
  USE parameters
  IMPLICIT NONE
  ! Input and output
  INTEGER                          :: ierror1, ierror2
  INTEGER                          :: n, k, n1, n3
  INTEGER                          :: prod, i, j
  INTEGER                          :: i_count
  INTEGER                          :: ion_num
  INTEGER                          :: spec_header_num, reac_header_num
  INTEGER                          :: species_count
  INTEGER            , ALLOCATABLE :: ion_temp(:)
  CHARACTER(LEN=10)                :: tempName
  CHARACTER(LEN=10)                :: line
  CHARACTER(LEN=10)  , ALLOCATABLE :: r_array(:,:)

  ! Open files
  OPEN (UNIT=1,FILE=SPECIES_FILE  ,STATUS='OLD',ACTION='READ',IOSTAT=ierror1)
  OPEN (UNIT=2,FILE=REACTIONS_FILE,STATUS='OLD',ACTION='READ',IOSTAT=ierror2)

  fileopen: IF ( ierror1 .EQ. 0 .AND. ierror2 .EQ. 0 ) THEN
     ! Read the species and reactions file header line numbers
     CALL linecount(1,ierror1,NUM_SPECIES,spec_header_num)
     CALL linecount(2,ierror2,NUM_REACTS,reac_header_num)

     ALLOCATE( REACT_CUBE(NUM_SPECIES,NUM_SPECIES,3), &
          SP_LIST(NUM_SPECIES),&
          EN_LIST(NUM_SPECIES), &
          IONLIST(ions) )
     ALLOCATE( ion_temp(NUM_SPECIES) )
     ALLOCATE( r_array(NUM_REACTS,5) )

     ! Initialize global values
     REACT_CUBE =  0
     EN_LIST    =  0
     IONLIST    =  0
     SP_LIST    = "0"

     ! Initialize local values
     ion_temp      =  0
     ion_num       =  0
     species_count =  0
     line          = '0'
     r_array       = '0'

     ! Read the contents of the species file and create the SP_LIST and energy_list
     DO n=1,NUM_SPECIES
        READ (1,*,IOSTAT=ierror1) line
        IF ( line(1:1) .NE. "!" ) THEN
           species_count = species_count + 1
           BACKSPACE (UNIT=1,IOSTAT=ierror1)
           READ(1,*,IOSTAT=ierror1) SP_LIST(species_count), EN_LIST(species_count)
           tempName = TRIM(SP_LIST(species_count))
           SP_LIST(species_count) = tempName
           ! Tally up the total number of anions and cations
           IF ( (tempName(LEN(TRIM(tempName)):LEN(TRIM(tempName))) .EQ. '-') .OR. &
                (tempName(LEN(TRIM(tempName)):LEN(TRIM(tempName))) .EQ. '+')) THEN
              IF ( tempName(1:1) .NE. '*') THEN 
                ion_num = ion_num + 2
                ion_temp(species_count) = species_count
              END IF
           END IF
        ELSE
           CONTINUE
        END IF
        IF ( ierror1 .NE. 0 ) STOP "Error reading species file."
     END DO

     ! Lookup to numbers of CRP and electron in the listj
     CALL lookup( "CRP", NUM_SPECIES, SP_LIST,  CRPNUM)
     CALL lookup( '*'  , NUM_SPECIES, SP_LIST,  EXCNUM )
     CALL lookup( 'e'  , NUM_SPECIES, SP_LIST,  ELECNUM )
     SPECIAL_LIST = (/ CRPNUM, EXCNUM, ELECNUM /)     

     i_count = 1
     DO n=1,NUM_SPECIES
        IF (ion_temp(n) .NE. 0 ) THEN
           IONLIST(i_count) = ion_temp(n)
           i_count = i_count + 1
        ELSE
           CONTINUE
        END IF
     END DO

     ! Read the contents of the reactions file and make the reactionCube
     qubemake: DO n1=1,NUM_REACTS
        READ(2,*,IOSTAT=ierror2) line
        IF ( line(1:1) .NE. "!" ) THEN
           BACKSPACE (UNIT=2,IOSTAT=ierror2)
           READ(2,*,IOSTAT=ierror2) r_array(n1,1), r_array(n1,2), r_array(n1,3), &
                r_array(n1,4), r_array(n1,5)
           CALL lookup(TRIM(r_array(n1,1)),NUM_SPECIES,SP_LIST,i)
           CALL lookup(TRIM(r_array(n1,2)),NUM_SPECIES,SP_LIST,j)
           k=0
           DO n3=3,5
              k = k + 1
              IF ( r_array(n1,n3) .NE. '0' ) THEN
                 CALL lookup(r_array(n1,n3),NUM_SPECIES,SP_LIST,prod)
                  REACT_CUBE(i,j,k) = prod
                 REACT_CUBE(j,i,k) = prod
              ELSE
                 CONTINUE
              END IF
           END DO
        ELSE
           CONTINUE
        END IF
     END DO qubemake
  ELSE fileopen
     IF ( ierror1 .NE. 0 .AND. ierror2 .NE. 0 ) THEN
        PRINT *, 'Could not open any of the input files'
     ELSE IF ( ierror1 .NE. 0 .AND. ierror2 .EQ. 0 ) THEN
        PRINT *, 'Could not open: ', species_file
     ELSE
        PRINT *, 'Could not open: ', reactions_file
     END IF
  END IF fileopen
END SUBROUTINE qbert

SUBROUTINE sort_energies( prods, A_en, B_en, C_en )
  USE parameters
  IMPLICIT NONE
  INTEGER, INTENT(INOUT), DIMENSION(3) :: prods
  INTEGER               , DIMENSION(3) :: temp
  REAL               , INTENT(IN)                  :: A_en, B_en, C_en

  IF ( A_en .GT. B_en ) THEN
     IF ( A_en .GT. C_en ) THEN
        IF ( B_en .GT. C_en ) THEN
           ! A > B > C
           temp(1) = prods(1)
           temp(2) = prods(2)
           temp(3) = prods(3)
        ELSE
           ! A > C > B
           temp(1) = prods(1)
           temp(2) = prods(3)
           temp(3) = prods(2)
        END IF
     ELSE
        ! C > A > B
        temp(1) = prods(3)
        temp(2) = prods(1)
        temp(3) = prods(2)
     END IF
  ELSE
     IF ( B_en .GT. C_en ) THEN
        IF ( A_en .GT. C_en ) THEN
           ! B > A > C
           temp(1) = prods(2)
           temp(2) = prods(1)
           temp(3) = prods(3)
        ELSE
           ! B > C > A
           temp(1) = prods(2)
           temp(2) = prods(3)
           temp(3) = prods(1)
        END IF
     ELSE
        ! C > B > A
        temp(1) = prods(3)
        temp(2) = prods(2)
        temp(3) = prods(1)
     END IF
  END IF

  prods = temp

END SUBROUTINE sort_energies
