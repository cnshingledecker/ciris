MODULE specdata
  USE typedefs
  !******************************************************************************
  ! Spectroscopic Parameters
  !******************************************************************************

  ! Excitation parameters for molecular oxygen (H+)
  ! Values taken from Edgar, Porter, and Green 1974
  TYPE(EPG_EXSTATE), PARAMETER, DIMENSION(5) :: o2_p_ex=(/    &
    EPG_EXSTATE("a\;^1\Delta_g",0.092,2.50E3,0.5,3.0,0.98)  , &
    EPG_EXSTATE("b\;^1\Sigma_g^+",0.109,4.19E3,0.5,3.0,1.64), &
    EPG_EXSTATE("A\;^3\Sigma_u^+",0.57,17.6E3,0.5,0.9,4.5)  , &
    EPG_EXSTATE("B\;^3\Sigma_u^-",4.73,51.7E3,0.5,0.75,8.4) , &
    EPG_EXSTATE("9.9 eV state",0.83,80.7E3,0.5,0.85,9.9)      &
  /)

  ! Ionization parameters for molecular oxygen (H+)
  ! Values taken from Edgar, Porter, and Green 1974
  TYPE(EPG_IONSTATE), PARAMETER, DIMENSION(1) :: o2_p_ion=(/            &
    EPG_IONSTATE("X\;^2\Pi_g",9.56D0,45.8D3,0.61D0,1.26D0,12.1D0) /) !   , &
!    EPG_IONSTATE("a\;^4\Pi_u",5.38D0,127.1D3,1.1D0,0.58D0,16.1D0)     , &
!    EPG_IONSTATE("A\;^2\Pi_g",5.25D0,170.3D3,1.1D0,0.60D0,16.9D0)     , &
!    EPG_IONSTATE("b\;^4\Sigma_g^-",3.19D0,142.4D3,1.1D0,0.58D0,18.2D0), &
!    EPG_IONSTATE("B\;State",0.94D0,140.8D3,1.1D0,0.58D0,23.0D0)       , &
!    EPG_IONSTATE("O^+\;^4S",15.6D0,56.0D3,1.2D0,0.72D0,18.0D0)        , &
!    EPG_IONSTATE("O^+\;^2D",8.75D0,57.0D3,1.2D0,0.72D0,22.0D0)          &
!  /)

  ! Excitation parameters for molecular oxygen by electrons (e-)
  ! Values taken from Porter, Jackman, Green 1976
  !****************************************************************************!
  !**********************Allowed Transitions***********************************!
  !****************************************************************************!
  TYPE(ALWD_EXSTATE), PARAMETER, DIMENSION(2)   :: o2_e_ex_alwd=(/     &
    ALWD_EXSTATE("B\;^3\Sigma_u^-",8.4,0.254,0.037,1.19,2.31,3.35,0.5),&
    ALWD_EXSTATE("9.9 eV state",9.9,0.0285,0.622,1.38,3.44,4.44,0.5)   &
  /)
  !****************************************************************************!
  !***********************Forbidden Transitions********************************!
  !****************************************************************************!
  TYPE(FBDN_EXSTATE), PARAMETER, DIMENSION(3)   :: o2_e_ex_fbdn=(/ &
    FBDN_EXSTATE("b\;^1\Sigma_g^+",1.64,0.0005,3.0,3.0,1.0)      , &
    FBDN_EXSTATE("a\;^1\Delta_g",0.98,0.0005,3.0,3.0,1.0)        , &
    FBDN_EXSTATE("A\;^3\Sigma_u^+",4.5,0.021,0.9,1.0,1.0)          &
  /)

  ! Ionization parameters for molecular oxygen by electrons (e-)
  ! Values taken from Jackman Garvey, Green 1977
  TYPE(IONSTATE), PARAMETER, DIMENSION(1)   :: o2_e_ion=(/                                   &
    IONSTATE("X\;^2\Pi_g",12.1,0.475,0.0,3.760,0.0,0.0,18.50,12.10,1.860,1000.0,24.20)      &
    ! IONSTATE("a\;^4\Pi_u",16.1,1.129,0.0,3.760,0.0,0.0,18.50,16.10,1.860,1000.0,32.20)     , &
    ! IONSTATE("A\;^2\Pi_u",16.9,1.129,0.0,3.760,0.0,0.0,18.50,16.90,1.860,1000.0,33.80)     , &
    ! IONSTATE("b\;^4\Sigma_g^-",18.2,1.010,0.0,3.760,0.0,0.0,18.50,18.20,1.860,1000.0,36.40), &
    ! IONSTATE("B\;^2\Sigma_g^-",20.0,0.653,0.0,3.760,0.0,0.0,18.50,20.30,1.860,1000.0,40.60), &
    ! IONSTATE("c\;^4\Sigma_u^-",23.0,0.950,0.0,3.760,0.0,0.0,18.50,23.00,1.860,1000.0,46.00), &
    ! IONSTATE("37-eV state",37.0,0.594,0.0,3.760,0.0,0.0,18.50,37.00,1.860,1000.0,74.00)      &
  /)
END MODULE specdata
