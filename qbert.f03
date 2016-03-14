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
SUBROUTINE qbert(num_species,num_reacts,qube,energy_array,speciesList,anions)
  USE subroutines
  USE parameters

  IMPLICIT NONE

  ! Input and output

  INTEGER(KIND=SHORT)              , DIMENSION(num_species,num_species,3) :: qube
  INTEGER                          , DIMENSION(IONS)                      :: anions
  REAL                             , DIMENSION(num_species)               :: energy_array ! Pointer to en_list
  CHARACTER(len=10)                , DIMENSION(num_species)               :: speciesList  ! Pointer to species list
  INTEGER(KIND=SHORT)                                                     :: nlines1, nlines2
  INTEGER                                                                 :: ierror1, ierror2
  INTEGER                                                                 :: n, k, n1, n3
  INTEGER                                                                 :: prod, i, j
  CHARACTER(len=10)                                                       :: tempName
  CHARACTER(len=10)   , ALLOCATABLE, DIMENSION(:,:)                       :: r_array
  CHARACTER(LEN=10)                                                       :: line
  INTEGER                                                                 :: i_count
  INTEGER                                                                 :: anion_num
  INTEGER             , ALLOCATABLE, DIMENSION(:)                         :: anion_temp
  REAL                                                                    :: Aqbert, Bqbert, Cqbert
  INTEGER(KIND=SHORT)              , DIMENSION(3)                         :: temp_prods
  INTEGER                                                                 :: spec_header_num, reac_header_num
  INTEGER                                                                 :: species_count
  INTEGER                                                                 :: num_species,num_reacts


  ! Open files
  OPEN (UNIT=1,FILE=SPECIES_FILE  ,STATUS='OLD',ACTION='READ',IOSTAT=ierror1)
  OPEN (UNIT=2,FILE=REACTIONS_FILE,STATUS='OLD',ACTION='READ',IOSTAT=ierror2)



  fileopen: IF ( ierror1 .EQ. 0 .AND. ierror2 .EQ. 0 ) THEN
    ! Read the species and reactions file header line numbers
    CALL linecount(1,ierror1,nlines1,spec_header_num)
    CALL linecount(2,ierror2,nlines2,reac_header_num)

    ALLOCATE( anion_temp(num_species) )
    ALLOCATE( r_array(num_reacts,5) )
    ! Initialize values
    anion_temp   =  0
    qube         =  0
    Aqbert       =  0
    Bqbert       =  0
    Cqbert       =  0
    anion_num    =  0
    speciesList  = '0'
    species_count = 0
    line = '0' 

    ! Read the contents of the species file and create the speciesList and energy_list
    DO n=1,nlines1
      READ (1,*,IOSTAT=ierror1) line 
      IF ( line(1:1) .NE. "!" ) THEN
        species_count = species_count + 1
        BACKSPACE (UNIT=1,IOSTAT=ierror1)
        READ(1,*,IOSTAT=ierror1) speciesList(species_count), energy_array(species_count)
        tempName = TRIM(speciesList(species_count))
        IF ( tempName(LEN(TRIM(tempName)):LEN(TRIM(tempName))) .EQ. '-' ) THEN
          anion_num = anion_num + 1
          anion_temp(species_count) = species_count
        END IF
      ELSE
        CONTINUE
      END IF
      IF ( ierror1 .NE. 0 ) STOP "Error reading species file."
    END DO

    i_count = 1
    DO n=1,num_species
      IF (anion_temp(n) .NE. 0 ) THEN
        anions(i_count) = anion_temp(n)
        i_count = i_count + 1
      ELSE
        CONTINUE
      END IF
    END DO

!    PRINT *, 'Now starting qubemake'
    ! Read the contents of the reactions file and make the reactionCube
    qubemake: DO n1=1,nlines2
      READ(2,*,IOSTAT=ierror2) line
      IF ( line(1:1) .NE. "!" ) THEN
        BACKSPACE (UNIT=2,IOSTAT=ierror2)
        READ(2,*,IOSTAT=ierror2) r_array(n1,1), r_array(n1,2), r_array(n1,3), &
                                 r_array(n1,4), r_array(n1,5)
!        PRINT *, r_array(n1,:)
        CALL lookup(r_array(n1,1),num_species,speciesList,i)
        CALL lookup(r_array(n1,2),num_species,speciesList,j)
        k=0
        DO n3=3,5
          k = k + 1
!          PRINT *, 'r_array(n1,n3)=',r_array(n1,n3)
          IF ( r_array(n1,n3) .NE. '0' ) THEN
            CALL lookup(r_array(n1,n3),num_species,speciesList,prod)
!            PRINT *, 'Prod of',r_array(n1,:),'is',prod
            qube(i,j,k) = prod
            qube(j,i,k) = prod
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

  ! Order the products of the reaction by binding energies in descending order
!  DO j=1,nlines1
!    DO i=1,nlines1
!      temp_prods = qube(i,j,:)
!      IF (temp_prods(1) .NE. 0 ) THEN
!        Aqbert = energy_array(temp_prods(1))
!      END IF
!
!      IF (temp_prods(2) .NE. 0 ) THEN
!        Bqbert = energy_array(temp_prods(2))
!      END IF
!
!      IF (temp_prods(3) .NE. 0 ) THEN
!        Cqbert = energy_array(temp_prods(3))
!      END IF
!
!      CALL sort_energies(temp_prods, Aqbert, Bqbert, Cqbert)
!      qube(i,j,:) = temp_prods
!      Aqbert = 0
!      Bqbert = 0
!      Cqbert = 0
!    END DO
!  END DO
END SUBROUTINE qbert

SUBROUTINE sort_energies( prods, A_en, B_en, C_en )
   USE parameters
   IMPLICIT NONE

   INTEGER(KIND=SHORT), INTENT(INOUT), DIMENSION(3) :: prods
   INTEGER(KIND=SHORT)               , DIMENSION(3) :: temp
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
