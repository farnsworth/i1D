!it is written to work with even number of fermions -> add for odd number of fermions
MODULE system
  !
  IMPLICIT NONE
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
  CHARACTER(LEN=500) :: datadir = 'data/'
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
  ! ... at given eigenstate it gives the correlation:
  !        xxcorrelation : for any eigenstate
  !        xxcorrelation_red : reduced space (only state with simmetry k,-k)
  !
  ! ... i must use size because f2py has some problem
  ! ... in using an array of size L (variable inside system)
  FUNCTION xxcorrelation(d, state, size_in )
    !
    INTEGER, INTENT(IN) :: d
    INTEGER, INTENT(IN) :: size_in
    LOGICAL, INTENT(IN), DIMENSION(size_in) :: state
    COMPLEX(kind = dp), DIMENSION ((2*d),(2*d)) :: matrix
    INTEGER :: j,k,delta,info,icol,irow,sign
    INTEGER,DIMENSION(2*d) :: pivot
    COMPLEX(kind = dp) :: det
    !
    REAL(kind= 8) :: xxcorrelation
    !
    !
    matrix = (0.0_dp,0.0_dp)
    det = (1.0_dp,0.0_dp)
    !
    ! to optimize this definition i can compute before vectors of BiAj, AiAj
    ! k is the row, l the column 
    DO k=1,d
       ! diagonal block
       irow = 2*k
       matrix( irow-1, irow   ) =  BiAj(1,state,size_in)
       matrix( irow  , irow-1 ) = -BiAj(1,state,size_in)
       !
       IF (k < d) THEN
          !
          DO j=k+1, d
             delta = j-k
             icol = 2*j
             matrix(irow-1,icol-1) = BiBj(delta,state,size_in) 
             matrix(irow-1,icol) = BiAj(delta+1,state,size_in)
             matrix(irow,icol-1) = -BiAj(-delta+1,state,size_in)
             matrix(irow,icol) = AiAj(delta,state,size_in)
! it is a skew-symmetric matrix
             matrix(icol-1,irow-1) = - matrix(irow-1,icol-1) 
             matrix(icol,irow-1) = - matrix(irow-1,icol)
             matrix(icol-1,irow) = - matrix(irow,icol-1)
             matrix(icol,irow) = - matrix(irow,icol)
          ENDDO
          !
       ENDIF
       !
    ENDDO
    !
    ! calculation of the determinant
    info = 1
    !
    pivot = 0
    !
    CALL zgetrf(2*d,2*d,matrix,2*d,pivot,info)
    !
    !print*,"info",info,"pivot",pivot
    !
    ! compute the determinant
    !
    DO k=1,2*d
       det = det*matrix(k,k)
    ENDDO
    !
    ! compute the sign of the determiant
    !
    sign = 0
    DO k=1,2*d
       IF ( pivot(k) /= k) sign = sign +1
    END DO
    !
    IF (mod(sign,2) == 1) THEN
       det = -det
    END IF
    !
    !
    xxcorrelation = dsqrt(dble(det))
    !
    !
  END FUNCTION xxcorrelation
  !
  !
  ! ... does a calculation of the correlation for states with 
  ! ... symmetric occupation k -k
  !
  FUNCTION xxcorrelation_red(d, state, size_in )
    !
    INTEGER, INTENT(IN) :: d
    INTEGER, INTENT(IN) :: size_in
    LOGICAL, INTENT(IN), DIMENSION(size_in) :: state
    REAL(kind = dp), DIMENSION ( d, d ) :: matrix
    REAL(kind=dp), DIMENSION (-d+2:d ) :: vectorBiAj
    INTEGER :: j,k,info,sign
    INTEGER,DIMENSION(d) :: pivot
    REAL(kind = dp) :: det
    !
    REAL(kind= 8) :: xxcorrelation_red
    !
    !
    matrix = 0.0_dp
    det = 1.0_dp
    !
    DO k=-d+2,d
       !
       vectorBiAj(k) = BiAj_red(k,state,size_in)
       !
    ENDDO
    !
    DO k=1,d
       !
       DO j=1,d
          !
          matrix( k , j ) =  vectorBiAj( j - k + 1 )
          !
       ENDDO
       !
    ENDDO
    !
    ! calculation of the determinant
    info = 1
    !
    pivot = 0
    !
    !
    CALL dgetrf( d , d ,matrix, d , pivot,info)
    !
    !
    !print*,"info",info,"pivot",pivot
    !
    ! compute the determinant
    !
    DO k=1,d
       det = det*matrix(k,k)
    ENDDO
    !
    ! compute the sign of the determiant
    !
    sign = 0
    DO k=1, d
       IF ( pivot(k) /= k) sign = sign +1
    END DO
    !
    IF (mod(sign,2) == 1) THEN
       det = -det
    END IF
    !
    xxcorrelation_red = det
    !
  END FUNCTION xxcorrelation_red
  !
  !
  !
  FUNCTION BiAj(d,state,size_in)
    !
    INTEGER, INTENT(IN) :: size_in
    INTEGER, INTENT(IN) :: d
    LOGICAL, INTENT(IN), DIMENSION(size_in) :: state
    COMPLEX (kind=8) :: BiAj
    !
    INTEGER :: i_tmp,term
    REAL(kind=dp) :: k,e,res
    !
    res = 0.0_dp
    !
    DO i_tmp=1,L/2
       !
       k = (pi*dble(2*i_tmp-1))/dble(L)
       e = energy(k,h,gamma)
       term = 1
       !
       IF (state(L/2 - i_tmp + 1)) term = term - 1
       IF (state(L/2 + i_tmp)) term = term - 1
       !
       SELECT CASE (term)
       CASE(-1) 
          res = res + ( - dcos( k*dble(d-1) ) + h * dcos(k*dble(d) ) )/e
       CASE(1) 
          res = res + ( dcos( k*dble(d-1) ) - h * dcos(k*dble(d)) )/e
       END SELECT
       !
    ENDDO
    !
    res = res*4.0_dp/dble(L)
    !
    BiAj = res*( 1.0_dp, 0.0_dp )
    !
  END FUNCTION BiAj
  !
  !
  FUNCTION AiAj(d,state,size_in)
    !
    INTEGER, INTENT(IN) :: d
    INTEGER, INTENT(IN) :: size_in
    LOGICAL, INTENT(IN), DIMENSION(size_in) :: state
    COMPLEX (kind=8) :: AiAj
    !
    INTEGER :: i_tmp,term,sign
    REAL(kind=dp) :: k,res
    !
    res = 0.0_dp
    !
    !
    IF ( mod(d,L) /= 0 ) THEN
       DO i_tmp=1,L/2
          !
          k = (pi*dble(2*i_tmp-1))/dble(L)
          !
          term = 0
          !
          IF (state(L/2 - i_tmp + 1)) term = term - 1
          IF (state(L/2 + i_tmp)) term = term + 1
          !
          SELECT CASE (term)
          CASE(-1) 
             res = res - dsin(k*dble(d))
          CASE(1) 
             res = res + dsin(k*dble(d))
          END SELECT
          !
       ENDDO
       !
       res = 2.0_dp*res/dble(L)
       AiAj = res*(0.0_dp,1.0_dp)
       !
    ELSE
       ! ... the sign is due to antiperiodic boundary condition
       sign = (-1)**(d/L)
       AiAj = dble(sign)*(1.0_dp,0.0_dp)
    ENDIF
    !
  END FUNCTION AiAj
  !
  !
  FUNCTION BiBj(d,state,size_in)
    !
    INTEGER, INTENT(IN) :: d
    INTEGER, INTENT(IN) :: size_in
    LOGICAL, INTENT(IN), DIMENSION(size_in) :: state
    COMPLEX (kind=8) :: BiBj 
    !
    BiBj = -AiAj(d,state,size_in)
    !
  END FUNCTION BiBj
  !
  !
  FUNCTION BiAj_red(d,state,size_in)
    !
    INTEGER, INTENT(IN) :: size_in
    INTEGER, INTENT(IN) :: d
    LOGICAL, INTENT(IN), DIMENSION(size_in) :: state
    COMPLEX (kind=8) :: BiAj_red
    !
    INTEGER :: i_tmp
    REAL(kind=dp) :: k,e,res
    !
    res = 0.0_dp
    !
    DO i_tmp=1,L/2
       !
       k = (pi*dble(2*i_tmp-1))/dble(L)
       e = energy(k,h,gamma)
       !
       ! non riesco a capire perche' il segno e' questo
       !
       IF (state(i_tmp)) THEN
          res = res - ( dcos( k*dble(d-1) ) - h * dcos(k*dble(d) ) )/e
       ELSE
          res = res + ( dcos( k*dble(d-1) ) - h * dcos(k*dble(d)) )/e
       ENDIF
       !
    ENDDO
    !
    BiAj_red = res*4.0_dp/dble(L)
    !
  END FUNCTION BiAj_red
  !
  !
  !
  FUNCTION state_ex_energy( state, size_in )
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: size_in
    LOGICAL, INTENT(IN), DIMENSION(size_in) :: state
    !
    INTEGER :: i
    REAL(kind= 8) :: state_ex_energy,k
    !
    !
    state_ex_energy = 0.0_dp
    !
    DO i=1,size_in
       !
       k = (pi*dble(2*i-1))/dble(L)
       !
       IF (state(i)) THEN
          state_ex_energy = state_ex_energy + 2.0_dp * energy(k,h,gamma)
       ENDIF
       !
    ENDDO
    !
    state_ex_energy = state_ex_energy/dble(L)
    !
    !
  END FUNCTION state_ex_energy
  !
  !
END MODULE system
