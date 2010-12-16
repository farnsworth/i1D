!it is written to work with even number of fermions -> add for odd number of fermions
MODULE system
  !
  INTEGER, PARAMETER :: dp = 8
  !
  !
  ! system parameters
  INTEGER :: L = 8
  REAL (kind=8) :: h = 1.0_dp,gamma = 1.0_dp
  REAL (kind=8) :: h0 = 1.0_dp,gamma0 = 1.0_dp
  !
  REAL (kind=8) :: pi = 3.141592653589793_dp
  !
  CHARACTER(LEN=30) :: datadir = 'data/'
  !
  !
CONTAINS
  !
  !
  FUNCTION sigmaz(k_in)
    !
    IMPLICIT NONE
    !
    REAL(kind=8) :: k_in
    REAL(kind=8) :: num_tmp,den_tmp
    REAL(kind=8) :: sigmaz
    !
    num_tmp  =  dcos(k_in) - h
    den_tmp  =  dsqrt( num_tmp*num_tmp + gamma*gamma*dsin(k_in)*dsin(k_in) )
    sigmaz   =  2.0_dp*num_tmp/den_tmp
    !
  END FUNCTION sigmaz
  !
  !
  !
  FUNCTION kinks(k_in)
    !
    IMPLICIT NONE
    !
    REAL(kind=8) :: k_in, kinks, e
    !
    e = energy(k_in, h, gamma)
    !
    kinks = 2.0d0 * ( 1.0d0 - h * dcos(k_in) ) / e
    !
  END FUNCTION kinks
  !
  !
  !
  FUNCTION energy(k_in, h_in, gamma_in)
    !
    IMPLICIT NONE
    !
    REAL(kind=8) :: k_in
    REAL(kind=8) :: h_in
    REAL(kind=8) :: gamma_in
    !
    REAL(kind=8) :: energy
    !
    energy = 2.0_dp*dsqrt( (dcos(k_in)-h_in)*(dcos(k_in)-h_in) + gamma_in*gamma_in*dsin(k_in)*dsin(k_in) )
    !
  END FUNCTION energy
  !
  !
  !
  FUNCTION Egs()
    !
    IMPLICIT NONE
    !
    INTEGER :: i
    REAL(kind=8) :: Egs,k
    !
    Egs = 0.0d0
    !
    DO i=1,L/2
       !
       k = (pi*dble(2*i-1))/dble(L)
       Egs = Egs + energy( k, h, gamma )
       !
    END DO
    !
    Egs = -Egs/dble(L)
    !
  END FUNCTION Egs
  !
  !
  !
  FUNCTION mgs()
    !
    IMPLICIT NONE
    !
    INTEGER :: i
    REAL(kind=8) :: mgs, k
    !
    mgs = 0.0d0
    !
    DO i=1,L/2
       !
       k = (pi*dble(2*i-1))/dble(L)
       mgs = mgs + sigmaz(k)
       !
    END DO
    !
    mgs = -mgs/dble(L)
    !
  END FUNCTION mgs
  !
  !
  !
  FUNCTION kinksgs()
    !
    IMPLICIT NONE
    !
    INTEGER :: i
    REAL(kind=8) :: kinksgs, k
    !
    kinksgs = 0.0d0
    DO i=1,L/2
       !
       k = (pi*dble(2*i-1))/dble(L)
       kinksgs = kinksgs + kinks(k) 
       !
    END DO
    !
    kinksgs = 0.5d0 - kinksgs/dble(L)
    !
  END FUNCTION kinksgs
  !
  !
  !
  SUBROUTINE error ( from_in, string_in )
    !
    IMPLICIT NONE
    CHARACTER( len = * ), INTENT(in) :: string_in
    CHARACTER( len = * ), INTENT(in) :: from_in
    !
    WRITE(UNIT=*,FMT=*), " "
    WRITE(UNIT=*,FMT='(a)'), "--------------- ERROR ---------------"
    WRITE(UNIT=*,FMT=*), " "
    WRITE(UNIT=*,FMT='(3a)'),"Error from ", from_in," subprogram"
    WRITE(UNIT=*,FMT='(2a)'),"Message: ",string_in
    WRITE(UNIT=*,FMT=*), " "
    !
  END SUBROUTINE error
  !
  !
  !
  SUBROUTINE warning ( from_in, string_in )
    !
    IMPLICIT NONE
    CHARACTER( len = * ), INTENT(in) :: string_in
    CHARACTER( len = * ), INTENT(in) :: from_in
    !
    WRITE(UNIT=*,FMT='(a)'), "--------------- WARNING ---------------"
    WRITE(UNIT=*,FMT='(4a)'),"From ", from_in,": ",string_in
    WRITE(UNIT=*,FMT=*), " "
    !
  END SUBROUTINE warning
  !
  !
END MODULE system
