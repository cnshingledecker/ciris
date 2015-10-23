!-*- f90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  System        :
!  Module        :
!  Object Name   : $RCSfile$
!  Revision      : $Revision$
!  Date          : $Date$
!  Author        : $Author$
!  Created By    : Christopher N. Shingledecker
!  Created       : Thu Jul 10 17:35:24 2014
!  Last Modified : <140711.1042>
!
!  Description
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
PROGRAM lookuptest
IMPLICIT NONE

INTEGER :: n
CHARACTER(len=10), DIMENSION(3) :: data
CHARACTER(len=10), DIMENSION(1) :: string_array
CHARACTER(len=10)               :: string

data = (/ 'a','b','c' /)
string = 'c'
CALL lookup(string,3,data,n)

PRINT *, n


END PROGRAM lookuptest
