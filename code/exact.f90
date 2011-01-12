MODULE exact
  !
  USE system
  !
  ! 
  REAL(kind = 8), DIMENSION (:,:), ALLOCATABLE :: spectrum
  INTEGER, DIMENSION (:,:), ALLOCATABLE :: array
  REAL(kind = 8), DIMENSION (:,:), ALLOCATABLE :: real_array
  REAL(kind = 8), DIMENSION (:), ALLOCATABLE :: obs
  INTEGER, DIMENSION(:), ALLOCATABLE :: par
  REAL(kind = 8) :: emax_plot, obsmax_plot 
  !
  !
CONTAINS
  !
  !
  SUBROUTINE spectrum_calc()
    IMPLICIT NONE
    
    INTEGER :: i,j,m
    REAL(kind = 8) :: econf,mconf,k
    REAL(kind = 8), DIMENSION(:,:),ALLOCATABLE :: evect,mvect
    REAL(kind = 8), DIMENSION(2) :: m0,e0
    
    IF ((mod(L,2) /= 0).or.(L>30)) THEN
       CALL error('exact','L is not good')
    END IF
    
    !the second component is for the odd subspace  
    ALLOCATE(evect(L,2),mvect(L,2))
    
    m0 = 0.0_dp
    e0 = 0.0_dp
    
    !k points for an even number of fermions
    DO i=1,L
       k = (pi*dble(2*i-L-1))/dble(L)

       evect(i,1) = energy( k, h, gamma)
       mvect(i,1) = sigmaz( k )
       
       m0(1) = m0(1) - 0.5_dp*mvect(i,1)
       e0(1) = e0(1) - 0.5_dp*evect(i,1)   
    END DO
    
    
    !k points for an odd number of fermions
    DO i=1,L
       k = -pi+dble(i)*(2.0_dp*pi)/dble(L)
       
      IF (i==L/2) THEN
          evect(i,2) = 2.0_dp*( -1.0_dp + h )
          mvect(i,2) =  -2.0_dp
          e0(2) = e0(2) - h
          m0(2) = m0(2) + 1.0_dp
       ELSE IF (i==L) THEN
          evect(i,2) =  2.0_dp*( 1.0_dp + h )
          mvect(i,2) =  -2.0_dp
          e0(2) = e0(2) - h
          m0(2) = m0(2) + 1.0_dp
       ELSE
          evect(i,2) = energy( k, h, gamma )
          mvect(i,2) = sigmaz( k )
          m0(2) = m0(2) - 0.5_dp*mvect(i,2)
          e0(2) =  e0(2) - 0.5_dp*evect(i,2)
       END IF
       
    END DO

    
    
    print*,"GROUND STATE - EVEN SUBSPACE (0 quasiparticles)"
    print*,"Energy =",e0(1)/dble(L)
    print*,"Local magnetization =", m0(1)/L
    print*
    print*
    print*,"GROUND STATE - ODD SUBSPACE (1 quasiparticle with k=0)"
    print*,"Energy =",(e0(2)+evect(L/2,2))/dble(L)
    print*,"Local magnetization =", (m0(2)+mvect(L/2,2))/L
    print*
    ! (8 byte per ogni dobule * 2 dble + 4 byte per ogni intero)
    print '(x,a,e10.2,a)',"estimate use of RAM:",(20.0_dp/1048576.0_dp)*dble(2**L)," Mb"
    
   
    IF (ALLOCATED(spectrum)) DEALLOCATE(spectrum)
    IF (ALLOCATED(par)) DEALLOCATE(par)
    
    ALLOCATE(spectrum(2**L,2))
    ALLOCATE(par(2**L))
   
    
    DO i=0,2**L-1
       
       m = 1 + mod(npart(i),2)
       
       par(i+1) = m - 1
       
       econf = e0(m)
       mconf = m0(m)
       
       DO j=1,L
          IF (btest(i,j-1)) THEN
             econf = econf + evect(j,m)
             mconf = mconf + mvect(j,m)
          END IF
       END DO
       
       spectrum(i+1,1) = econf/dble(L)
       spectrum(i+1,2) = mconf/dble(L)
       
    END DO

    DEALLOCATE(evect,mvect)
    
    
  CONTAINS
    
    FUNCTION npart(i_in)
      
      INTEGER, INTENT(IN) :: i_in
      INTEGER :: npart
      INTEGER :: i_tmp
      
      npart = 0
      DO i_tmp=0,L-1
         IF (btest(i_in,i_tmp)) npart = npart + 1
      END DO
      
    END FUNCTION npart
    
  END SUBROUTINE spectrum_calc
  !
  !
  !
  ! ... Magnetiation for all alpha states
  !
  SUBROUTINE magalpha_calc ()
    !
    IMPLICIT NONE
    !
    REAL(kind = 8), DIMENSION (:), ALLOCATABLE :: mvect
    REAL(kind = 8) :: mgs, k
    INTEGER :: i,j
    !
    ALLOCATE( mvect(L/2) )
    !
    mgs = 0.0_dp
    !
    DO i=1,L/2
       !
       k = (pi*dble(2*i-1))/dble(L)
       mvect(i) = sigmaz( k )/dble(L)
       mgs = mgs - mvect(i)
       !
    END DO
    !
    IF (allocated(obs)) deallocate(obs)
    ALLOCATE( obs( 2**(L/2) ))
    !
    obs(:) = mgs
    !
    DO i=0,2**(L/2)-1
       !
       DO j=1,L/2
          IF ( btest(i,j-1) ) THEN
             obs(i+1) = obs(i+1) + 2.0_dp*mvect(j)
          END IF
       END DO
       !
    END DO
    !
    PRINT '(x,a,e10.2,a)','estimate use of RAM:',(8.0_dp/1048576.0_dp)*dble(2**(L/2))," Mb"
    !
    DEALLOCATE(mvect)
    !
  END SUBROUTINE magalpha_calc
  !
  !
  ! ... Energy for all alpha states
  !
  SUBROUTINE ealpha_calc ()
    !
    IMPLICIT NONE
    !
    REAL(kind = 8), DIMENSION (:), ALLOCATABLE :: evect
    REAL(kind = 8) :: egs, k
    INTEGER :: i,j
    !
    ALLOCATE( evect(L/2) )
    !
    !k points for an even number of fermions
    egs = 0.0_dp
    !
    DO i=1,L/2
       !
       k = (pi*dble(2*i-1))/dble(L)
       evect(i) = energy( k, h, gamma )/dble(L)
       egs = egs - evect(i) !twice because we are doing a summation only over k>0
       !
    END DO
    !
    IF (allocated(obs)) deallocate(obs)
    ALLOCATE( obs( 2**(L/2) ))
    !
    obs(:) = egs
    !
    DO i=0,2**(L/2)-1
       !
       DO j=1,L/2
          IF ( btest(i,j-1) ) THEN
             obs(i+1) = obs(i+1) + 2.0_dp*evect(j)
          END IF
       END DO
       !
    END DO
    !
    PRINT '(x,a,e10.2,a)','estimate use of RAM:',(8.0_dp/1048576.0_dp)*dble(2**(L/2))," Mb"
    !
    DEALLOCATE(evect)
    !
  END SUBROUTINE ealpha_calc
  !
  !
  !
  SUBROUTINE mag_array_calc( ebin, mbin )
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ebin, mbin
    !
    REAL(kind = 8), DIMENSION (:), ALLOCATABLE :: evect, mvect
    REAL(kind = 8) :: egs, mmin, mgs_loc, k, deltam, deltae
    INTEGER :: i
    !
    REAL(kind = 8) :: malpha, ealpha
    INTEGER :: j, iebin, imbin
    !
    IF (allocated(array)) deallocate(array)
    allocate(array(ebin,mbin))
    !
    !
    allocate( evect(L/2), mvect(L/2) )
    !
    !
    !k points for an even number of fermions
    egs = 0.0_dp
    mgs_loc = 0.0_dp
    mmin = 0.0_dp
    !
    DO i=1,L/2
       !
       k = (pi*dble(2*i-1))/dble(L)
       !
       evect(i)  = energy( k, h, gamma )/dble(L)
       egs = egs - evect(i)
       !
       mvect(i) = sigmaz(k)/dble(L)
       !
       mgs_loc = mgs_loc - mvect(i)
       mmin = mmin - dabs( mvect(i) ) 
       !
    END DO
    !
    !
    deltam = 2.0_dp*dabs(mmin)/dble(mbin)
    deltae = 2.0_dp*dabs(egs)/dble(ebin)
    !
    emax_plot = dabs(egs)
    obsmax_plot = dabs(mmin)
    !
    array = 0
    !
    DO i=0,2**(L/2)-1
       ealpha = egs
       malpha  = mgs_loc
       !
       DO j=1,L/2
          IF ( btest(i,j-1) ) THEN
             malpha = malpha + 2.0_dp*mvect(j)
             ealpha = ealpha + 2.0_dp*evect(j)
          END IF
       END DO
       !
       !
       imbin = int((malpha-mmin)/deltam)+1
       iebin = int((ealpha-egs)/deltae)+1
       !
!       IF (i == 0) THEN
!          print*,"malpha",malpha,"imbin",imbin,"iebin",iebin
!       ENDIF
       !
       !
       IF ( (imbin <= mbin).and.( imbin > 0 ).and.(iebin <= ebin) .and. (iebin > 0) ) THEN
          array(iebin, imbin) = array(iebin, imbin) + 1
       ELSE
          IF (iebin > ebin) array(ebin, imbin) = array(ebin, imbin) + 1
          IF (imbin > mbin) array(iebin, mbin) = array(iebin, mbin) + 1
       END IF
       !
    END DO
    !
    !
    DEALLOCATE(mvect, evect)
    !
    !
    RETURN
    !
  END SUBROUTINE mag_array_calc
  !
  !
  !
  SUBROUTINE kinks_array_calc( ebin, kinkbin )
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ebin, kinkbin
    !
    REAL(kind = 8), DIMENSION (:), ALLOCATABLE :: evect, kinkvect
    REAL(kind = 8) :: egs, kinkmin, kinkgs_loc, k, deltakink, deltae
    INTEGER :: i
    !
    REAL(kind = 8) :: kinkalpha, ealpha
    INTEGER :: j, iebin, ikinkbin
    !
    IF ( allocated(array) ) deallocate(array)
    allocate( array(ebin,kinkbin) )
    !
    !
    allocate( evect(L/2), kinkvect(L/2) )
    !
    !
    !k points for an even number of fermions
    egs = 0.0_dp
    kinkgs_loc = 0.0_dp
    kinkmin = 0.0_dp
    !
    DO i=1,L/2
       !
       k = (pi*dble(2*i-1))/dble(L)
       !
       evect(i)  = energy( k, h, gamma )/dble(L)
       egs = egs - evect(i)
       !
       kinkvect(i) = kinks(k)/dble(L)
       !
       kinkgs_loc = kinkgs_loc - kinkvect(i)
       kinkmin = kinkmin - dabs( kinkvect(i) ) 
       !
    END DO
    !
    print*,"kinkgs_loc", 0.5 + kinkgs_loc,"kinkgs",kinksgs()
    !
    deltakink = 2.0_dp*dabs(kinkmin)/dble(kinkbin)
    deltae = 2.0_dp*dabs(egs)/dble(ebin)
    !
    emax_plot = dabs(egs)
    obsmax_plot = dabs( kinkmin )
    !
    array = 0
    !
    DO i=0,2**(L/2)-1
       ealpha = egs
       kinkalpha  = kinkgs_loc
       !
       DO j=1,L/2
          IF ( btest(i,j-1) ) THEN
             kinkalpha = kinkalpha + 2.0_dp*kinkvect(j)
             ealpha = ealpha + 2.0_dp*evect(j)
          END IF
       END DO
       !
       !
       ikinkbin = int((kinkalpha-kinkmin)/deltakink)+1
       iebin = int((ealpha-egs)/deltae)+1
       !
       !
       IF ( (ikinkbin <= kinkbin).and.( ikinkbin > 0 ).and.(iebin <= ebin) .and. (iebin > 0) ) THEN
          array(iebin, ikinkbin) = array(iebin, ikinkbin) + 1
       ELSE
          IF (iebin > ebin) array(ebin, ikinkbin) = array(ebin, ikinkbin) + 1
          IF (ikinkbin > kinkbin) array(iebin, kinkbin) = array(iebin, kinkbin) + 1
       END IF
       !
    END DO
    !
    DEALLOCATE(kinkvect, evect)
    !
    !
    RETURN
    !
  END SUBROUTINE kinks_array_calc
  !
  !
  !
  SUBROUTINE mag_array_calc_allspace( ebin, mbin )
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ebin, mbin
    !
    REAL(kind = 8), DIMENSION (:), ALLOCATABLE :: evect, mvect
    REAL(kind = 8) :: egs, mmin, mgs_loc, k, deltam, deltae
    INTEGER :: i
    !
    REAL(kind = 8) :: malpha, ealpha
    INTEGER :: j, iebin, imbin
    !
    IF (allocated(array)) deallocate(array)
    allocate(array(ebin,mbin))
    !
    IF (allocated(real_array)) deallocate(real_array)
    allocate(real_array(ebin,mbin))
    !
    allocate( evect(L), mvect(L) )
    !
    !
    !k points for an even number of fermions
    egs = 0.0_dp
    mgs_loc = 0.0_dp
    mmin = 0.0_dp
    !
    DO i=1,L
       !
       k = (pi*dble(2*i-1-L))/dble(L)
       !
       evect(i)  = energy( k, h, gamma )/dble(L)
       egs = egs - 0.5d0*evect(i)
       !
       mvect(i) = sigmaz(k)/dble(L)
       !
       mgs_loc = mgs_loc - 0.5d0*mvect(i)
       mmin = mmin - 0.5d0*dabs( mvect(i) ) 
       !
    END DO
    !
    !
    deltam = 2.0_dp*dabs(mmin)/dble(mbin)
    deltae = 2.0_dp*dabs(egs)/dble(ebin)
    !
    emax_plot = dabs(egs)
    obsmax_plot = dabs(mmin)
    !
    array = 0
    !
    DO i=0,2**(L)-1
       ealpha = egs
       malpha  = mgs_loc
       !
       DO j=1,L
          IF ( btest(i,j-1) ) THEN
             malpha = malpha + mvect(j)
             ealpha = ealpha + evect(j)
          END IF
       END DO
       !
       imbin = int((malpha-mmin)/deltam)+1
       iebin = int((ealpha-egs)/deltae)+1
       !
       IF ( (imbin <= mbin).and.( imbin > 0 ).and.(iebin <= ebin) .and. (iebin > 0) ) THEN
          array(iebin, imbin) = array(iebin, imbin) + 1
!       ELSE
!          IF ((iebin > ebin).and.(imbin <= mbin )) array(ebin, imbin) = array(ebin, imbin) + 1
!          IF ((imbin > mbin).and.(iebin <= ebin)) array(iebin, mbin) = array(iebin, mbin) + 1
       END IF
       !
    END DO
    !
    DEALLOCATE(mvect, evect)
    !
    !
    real_array = dble( array ) / ( deltam*deltae*dble(2**L) )
    !
    !
    RETURN
    !
  END SUBROUTINE mag_array_calc_allspace
  !
  !
  ! ... Correlation function in the x direction for all alpha states
  !
  SUBROUTINE corralpha_calc (d)
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: d
    INTEGER :: i,j
    LOGICAL, DIMENSION(L) :: state
    REAL(kind = dp), DIMENSION (L) :: evect
    REAL(kind = dp) :: ener,k
    !
    IF (allocated(real_array)) deallocate(real_array)
    ALLOCATE( real_array( 2, 2**L ) )
    !
    real_array = 0.0_dp
    !
    ener = 0.0_dp
    !
    DO i=1,L/2
       !
       k = (pi*dble(2*i-1))/dble(L)
       !
       evect(L/2 + i) = energy(k,h,gamma)
       evect(L/2 - i + 1) = evect(L/2 + i)
       !
       ener = ener - evect(L/2 + i)
       !
    ENDDO
    !
    ener = ener/dble(L)
    evect = evect/dble(L)
    !
    real_array( 1 , : ) = ener
    !
    DO i=0,2**L-1
       !
       state = .false.
       !
       DO j=1,L
          IF ( btest(i,j-1) ) THEN
             state(j) = .true.
             real_array(1,i+1) = real_array(1,i+1) + evect(j)
          END IF
       END DO
       !
       ! 
       real_array(2,i+1) = xxcorrelation(d,state,L)
       !
    END DO
    !
    !
  END SUBROUTINE corralpha_calc
  !
  !
  !
  ! ... Correlation function in the x direction for all alpha states
  ! ... with the simmetry k,-k
  !
  SUBROUTINE corralpha_calc_simp (d)
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: d
    INTEGER :: i,j
    LOGICAL, DIMENSION(L/2) :: state
    REAL(kind = dp), DIMENSION (L/2) :: evect
    REAL(kind = dp) :: ener,k
    !
    IF (allocated(real_array)) deallocate(real_array)
    ALLOCATE( real_array( 2, 2**(L/2) ) )
    !
    real_array = 0.0_dp
    !
    ener = 0.0_dp
    !
    DO i=1,L/2
       !
       k = (pi*dble(2*i-1))/dble(L)
       !
       evect( i ) = energy(k,h,gamma)
       !
       ener = ener - evect( i )
       !
    ENDDO
    !
    ener = ener/dble(L)
    evect = evect/dble(L)
    !
    real_array( 1 , : ) = ener
    !
    DO i=0,2**(L/2)-1
       !
       state = .false.
       !
       DO j=1,L/2
          !
          IF ( btest(i,j-1) ) THEN
             state(j) = .true.
             real_array(1,i+1) = real_array(1,i+1) + 2.0*evect(j)
          END IF
          !
       END DO
       !
       ! 
       real_array(2,i+1) = xxcorrelation_simp(d,state,L/2)
       !
    END DO
    !
    !
  END SUBROUTINE corralpha_calc_simp
  !
  !
  ! ... i must use size because f2py has some problem
  ! ... in using an array of size L (variable inside system)
  FUNCTION xxcorrelation(d, state, size )
    !
    INTEGER, INTENT(IN) :: d
    INTEGER, INTENT(IN) :: size
    LOGICAL, INTENT(IN), DIMENSION(size) :: state
    COMPLEX(kind = dp), DIMENSION ((4*d-2),(4*d-2)) :: matrix
    INTEGER :: j,k,delta,info,icol,irow,sign
    INTEGER,DIMENSION(4*d-2) :: pivot
    COMPLEX(kind = dp) :: det
    !
    REAL(kind= dp) :: xxcorrelation
    !
    !
    matrix = (0.0_dp,0.0_dp)
    det = (1.0_dp,0.0_dp)
    !
    ! to optimize this definition i can compute before vectors of BiAj, AiAj
    ! k is the row, l the column 
    DO k=1,2*d-1
       ! diagonal block
       irow = 2*k
       matrix( irow-1, irow   ) =  BiAj(1,state,size)
       matrix( irow  , irow-1 ) = -BiAj(1,state,size)
       !
       IF (k<2*d-1) THEN
          !
          DO j=k+1,2*d-1
             delta = j-k
             icol = 2*j
             matrix(irow-1,icol-1) = BiBj(delta,state,size) 
             matrix(irow-1,icol) = BiAj(delta+1,state,size)
             matrix(irow,icol-1) = -BiAj(-delta+1,state,size)
             matrix(irow,icol) = AiAj(delta,state,size)
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
    CALL zgetrf(4*d-2,4*d-2,matrix,4*d-2,pivot,info)
    !
    !print*,"info",info,"pivot",pivot
    !
    ! compute the determinant
    !
    DO k=1,4*d-2
       det = det*matrix(k,k)
    ENDDO
    !
    ! compute the sign of the determiant
    !
    sign = 0
    DO k=1,4*d-2
       IF ( pivot(k) /= k) sign = sign +1
    END DO
    !
    IF (mod(sign,2) == 1) THEN
       det = -det
    END IF
    !
    !print*,"determinant",det
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
  FUNCTION xxcorrelation_simp(d, state, size )
    !
    INTEGER, INTENT(IN) :: d
    INTEGER, INTENT(IN) :: size
    LOGICAL, INTENT(IN), DIMENSION(size) :: state
    REAL(kind = dp), DIMENSION ( d, d ) :: matrix
    REAL(kind=dp), DIMENSION (-d+2:d ) :: vectorBiAj
    INTEGER :: j,k,info,sign
    INTEGER,DIMENSION(d) :: pivot
    REAL(kind = dp) :: det
    !
    REAL(kind= dp) :: xxcorrelation_simp
    !
    !
    matrix = 0.0_dp
    det = 1.0_dp
    !
    DO k=-d+2,d
       !
       vectorBiAj(k) = BiAj_simp(k,state,size)
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
    !print*,"matrix",matrix
    !
    CALL dgetrf( d , d ,matrix, d , pivot,info)
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
    !print*,"determinant",det
    !
    xxcorrelation_simp = det
    !
  END FUNCTION xxcorrelation_simp
  !
  !
  FUNCTION BiAj_simp(d,state,size)
    !
    INTEGER, INTENT(IN) :: size
    INTEGER, INTENT(IN) :: d
    LOGICAL, INTENT(IN), DIMENSION(size) :: state
    COMPLEX (kind=dp) :: BiAj_simp
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
          res = res - ( - dcos( k*dble(d-1) ) + h * dcos(k*dble(d)) )/e
       ENDIF
       !
    ENDDO
    !
    BiAj_simp = res*4.0_dp/dble(L)
    !
  END FUNCTION BiAj_simp
  !
  !
  !
  FUNCTION BiAj(d,state,size)
    !
    INTEGER, INTENT(IN) :: size
    INTEGER, INTENT(IN) :: d
    LOGICAL, INTENT(IN), DIMENSION(size) :: state
    COMPLEX (kind=dp) :: BiAj
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
  FUNCTION AiAj(d,state,size)
    !
    INTEGER, INTENT(IN) :: d
    INTEGER, INTENT(IN) :: size
    LOGICAL, INTENT(IN), DIMENSION(size) :: state
    COMPLEX (kind=dp) :: AiAj
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
  FUNCTION BiBj(d,state,size)
    !
    INTEGER, INTENT(IN) :: d
    INTEGER, INTENT(IN) :: size
    LOGICAL, INTENT(IN), DIMENSION(size) :: state
    COMPLEX (kind=dp) :: BiBj 
    !
    BiBj = -AiAj(d,state,size)
    !
  END FUNCTION BiBj
  !
  !
END MODULE exact
