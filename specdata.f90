MODULE specdata
  USE typedefs
  !******************************************************************************
  ! Spectroscopic Parameters
  !******************************************************************************

  ! Excitation parameters for molecular oxygen (H+)
  ! Values taken from Edgar, Porter, and Green 1974
  TYPE(EPG_EXSTATE), PARAMETER, DIMENSION(5) :: o2_p_ex=(/    &
       EPG_EXSTATE("a\;^1\Delta_g"  ,0.092,2.50,0.5,3.0 ,0.98), &
       EPG_EXSTATE("b\;^1\Sigma_g^+",0.109,4.19,0.5,3.0 ,1.64), &
       EPG_EXSTATE("A\;^3\Sigma_u^+",0.57 ,17.6,0.5,0.9 ,4.5) , &
       EPG_EXSTATE("B\;^3\Sigma_u^-",4.73 ,51.7,0.5,0.75,8.4) , &
       EPG_EXSTATE("9.9 eV state"   ,0.83 ,80.7,0.5,0.85,9.9)   &
       /)

  ! Excitation parameters for atomic oxygen
  ! Values taken from Edgar, Porter, and Green 1974
  TYPE(EPG_EXSTATE), PARAMETER, DIMENSION(4) :: o_p_ex = (/ &
       EPG_EXSTATE("^1D",0.51 ,5.4 ,1.0,1.0,1.85), &
       EPG_EXSTATE("^1S",0.075,8.7 ,1.0,1.0,4.18), &
       EPG_EXSTATE("^3S",0.38 ,37.2,1.0,1.0,9.53),&
       EPG_EXSTATE("^5S",0.55 ,16.1,1.0,1.0,9.20) &
       /)

  ! Excitation parameters for ozone
  ! Values are taken from Shingledecker et al. 2017
  TYPE(EPG_EXSTATE), PARAMETER, DIMENSION(1) :: o3_p_ex = (/ &
       EPG_EXSTATE("^1B_2", 1.0, 20.00, 1.0, 2.5, 4.9) /)

  ! Ionization parameters for molecular oxygen (H+)
  ! Values taken from Edgar, Porter, and Green 1974
  TYPE(EPG_IONSTATE), PARAMETER, DIMENSION(7) :: o2_p_ion=(/            &
       EPG_IONSTATE("X\;^2\Pi_g",9.56D0,45.8,0.61D0,1.26D0,12.1D0)     , &
       EPG_IONSTATE("a\;^4\Pi_u",5.38D0,127.1,1.1D0,0.58D0,16.1D0) , &
       EPG_IONSTATE("A\;^2\Pi_g",5.25D0,170.3,1.1D0,0.60D0,16.9D0) , &
       EPG_IONSTATE("b\;^4\Sigma_g^-",3.19D0,142.4,1.1D0,0.58D0,18.2D0), &
       EPG_IONSTATE("B\;State",0.94D0,140.8,1.1D0,0.58D0,23.0D0) , &
       EPG_IONSTATE("O^+\;^4S",15.6D0,56.0,1.2D0,0.72D0,18.0D0) , &
       EPG_IONSTATE("O^+\;^2D",8.75D0,57.0,1.2D0,0.72D0,22.0D0) &
       /)

  ! Ionization parameters for atomic oxygen (H+)
  ! Values taken from Edgar, Porter, and Green 1974
  TYPE(EPG_IONSTATE), PARAMETER, DIMENSION(3) :: o_p_ion = (/ &
       EPG_IONSTATE("OII\;^4S",20.4,61.5,0.82,0.75,13.6), &
       EPG_IONSTATE("O!!\;^2D",25.4,69.5,0.82,0.75,16.9), &
       EPG_IONSTATE("O11\;^4P",10.8,72.8,0.82,0.75,18.5) /)

  ! Ionization parameters for ozone
  ! Values are taken from Shingledecker et al. 2017
  TYPE(EPG_IONSTATE), PARAMETER, DIMENSION(1) :: o3_p_ion = (/ &
       EPG_IONSTATE("^1A_1", 40.00, 105.00, 1.00, 0.75, 12.43 ) /)

  ! Excitation parameters for molecular oxygen by electrons (e-)
  ! Values taken from Porter, Jackman, Green 1976
  !****************************************************************************!
  !**********************Allowed Transitions***********************************!
  !****************************************************************************!
  TYPE(ALWD_EXSTATE), PARAMETER, DIMENSION(2)   :: o2_e_ex_alwd=(/     &
       ALWD_EXSTATE("B\;^3\Sigma_u^-",8.4,0.254,0.037,1.19,2.31),&
       ALWD_EXSTATE("9.9 eV state",9.9,0.0285,0.622,1.38,3.44)   &
       /)

  ! Excitation parameters for molecular oxygen by electrons (e-)
  ! Values taken from Porter, Jackman, Green 1976
  !****************************************************************************!
  !***********************Forbidden Transitions********************************!
  !****************************************************************************!
  TYPE(FBDN_EXSTATE), PARAMETER, DIMENSION(3)   :: o2_e_ex_fbdn=(/ &
       FBDN_EXSTATE("b\;^1\Sigma_g^+",1.64,0.0005,3.0,3.0,1.0)      , &
       FBDN_EXSTATE("a\;^1\Delta_g",0.98,0.0005,3.0,3.0,1.0)        , &
       FBDN_EXSTATE("A\;^3\Sigma_u^+",4.5,0.021,0.9,1.0,1.0)          &
       /)

  ! Excitation parameters for atomic oxygen by electrons (e-)
  ! Values taken from Jackman, Garvey, & Green 1977
  !****************************************************************************!
  !**********************Allowed Transitions***********************************!
  !****************************************************************************!
  TYPE(ALWD_EXSTATE), PARAMETER, DIMENSION(2)   :: o_e_ex_alwd=(/     &
       ALWD_EXSTATE("^3S^0",9.53,0.0560,0.320,0.86,1.440),   &
       ALWD_EXSTATE("^3D^0",12.10,0.0310,0.610,1.26,0.490)  &
       /)

  ! Excitation parameters for atomic oxygen by electrons (e-)
  ! Values taken from Jackman, Garvey, & Green 1977
  !****************************************************************************!
  !***********************Forbidden Transitions********************************!
  !****************************************************************************!
  TYPE(FBDN_EXSTATE), PARAMETER, DIMENSION(2)   :: o_e_ex_fbdn=(/ &
       FBDN_EXSTATE("^1D",1.85,0.0100,1.000,1.00,2.00)      , &
       FBDN_EXSTATE("^1S",4.18,0.0042,1.00,0.50,1.000)      &
       /)

  ! Ionization parameters for molecular oxygen by electrons (e-)
  ! Values taken from Jackman Garvey, Green 1977
  TYPE(IONSTATE), PARAMETER, DIMENSION(1)   :: o2_e_ion=(/                                   &
       IONSTATE("X\;^2\Pi_g",12.1,0.475,0.0,3.760,0.0,0.0,18.50,12.10,1.860,1000.0,24.20) /)

  ! Ionization parameters for atomic oxygen by electrons (e-)
  ! Values taken from Jackman Garvey, Green 1977
  TYPE(IONSTATE), PARAMETER, DIMENSION(1)   :: o_e_ion=(/                                   &
       IONSTATE("^4S^0",13.6,1.03,0.0,1.81,0.0,0.0,13.0,-0.815,6.41,3450.0,162.0) /)

END MODULE specdata
