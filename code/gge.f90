MODULE gge
  !
  USE system
  !
CONTAINS
  !
  FUNCTION occupation(k)

    IMPLICIT NONE
    !
    REAL(kind=8) :: k,occupation,temp,e,e0
    !
    temp = 1.0_dp - (h0+h)*dcos(k)+h*h0
    e = energy(k,h,gamma)
    e0 = energy(k,h0,gamma)
    !
    occupation = 0.5d0 - 2.0d0*temp/(e*e0)
    !
  ENDFUNCTION occupation
  !
  FUNCTION gge_BiAj (d)
    !                                                                     
    INTEGER, INTENT(IN) :: d
    REAL (kind = 8) :: gge_BiAj, term, k,ek
    !                                                                     
    gge_BiAj = 0.0_dp
    !                                                                     
    DO i_tmp=1,L/2
       !                                                                  
       k = (pi*dble(2*i_tmp-1))/dble(L)
       ek = energy(k,h,gamma)
       !                                                                  
       term = (dcos(k*dble(d-1))-h*dcos(k*dble(d))) / ek
       !                                                                  
       gge_BiAj = gge_BiAj + term * ( 1.0d0 - 2.0d0*occupation(k) )
       !                                                                  
    ENDDO
    !                                                                     
    gge_BiAj = gge_BiAj*4.0d0/dble(L)
    !                                                                     
    !                                                                     
  END FUNCTION gge_BiAj

  !
  FUNCTION gge_xxcorrealtion(d)
    !
    INTEGER, INTENT(IN) :: d
    REAL(kind = dp), DIMENSION ( d, d ) :: matrix
    REAL(kind=dp), DIMENSION (-d+2:d ) :: vectorBiAj
    INTEGER :: j,k,info,sign
    INTEGER,DIMENSION(d) :: pivot
    REAL(kind = dp) :: det
    !
    REAL(kind= 8) :: gge_xxcorrelation
    !
    !
    matrix = 0.0_dp
    det = 1.0_dp
    !
    DO k=-d+2,d
       !
       vectorBiAj(k) = gge_BiAj(k)
       !
    ENDDO
    !
    DO k=1,d
       !
       DO j=1,d
          matrix( k , j ) =  vectorBiAj( j - k + 1 )
       ENDDO
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
    gge_xxcorrelation = det
    !
  END FUNCTION gge_xxcorrealtion
  !
END MODULE gge
