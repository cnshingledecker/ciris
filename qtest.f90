!-*- f90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  System        :
!  Module        :
!  Object Name   : $RCSfile$
!  Revision      : $Revision$
!  Date          : $Date$
!  Author        : $Author$
!  Created By    : Christopher N. Shingledecker
!  Created       : Fri Jul 11 12:33:07 2014
!  Last Modified : <140717.1657>
!
!  Description: This is a test driver for the subroutine, qbert
!
!  Notes
!
!  History
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
PROGRAM qtest
USE chaco_data
IMPLICIT NONE

CHARACTER(len=80) :: species_file
CHARACTER(len=80) :: reactions_file
CHARACTER(len=10), ALLOCATABLE, DIMENSION(:) :: sp_list
INTEGER, ALLOCATABLE, DIMENSION(:,:,:)  :: qube
REAL, ALLOCATABLE, DIMENSION(:,:)  :: en_list
INTEGER :: err1, err2, unit1, unit2
INTEGER :: lines1,lines2, n, i,j

species_file = "species.txt"
reactions_file = "reactions.txt"
OPEN (UNIT=1, FILE=species_file, STATUS='OLD', ACTION='READ', IOSTAT=err1)
OPEN (UNIT=2, FILE=reactions_file, STATUS='OLD', ACTION='READ', IOSTAT=err2)

CALL linecount(1,err1,lines1)
CALL linecount(2,err2,lines2)

ALLOCATE(qube(lines1,lines1,3),sp_list(lines1),en_list(lines1,3))

CALL qbert(species_file, reactions_file, lines1, lines2, qube, sp_list, en_list)
!PRINT *, SIZE(sp_list)

DO j=1,lines1
  DO i=1,lines1
    PRINT *, qube(i,j,:)
!    PRINT *, ' '
  END DO
  PRINT *, ' '
END DO
  

END PROGRAM qtest
