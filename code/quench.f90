MODULE quench
  !
  USE system
  !
  ! 
  REAL(kind = 8), DIMENSION (:), ALLOCATABLE :: calpha2
  !
  REAL(kind = 8) :: mean_en, sigma_en,sigma_obs
  REAL(kind = 8), DIMENSION (:,:), ALLOCATABLE :: dist
  !
  !
CONTAINS
  !
  ! ... |c_\alpha|^2 for all alpha state
  !
  SUBROUTINE calpha2_calc()
    !
    IMPLICIT NONE
    !
    REAL(kind = 8) :: e0,e,control
    REAL(kind = 8), DIMENSION (:), ALLOCATABLE :: coeff
    REAL(kind = 8) :: cgs, k
    INTEGER :: i,j
    !
    ALLOCATE( coeff(L/2) )
    !
    !k points for an even number of fermions
    cgs = 1.0_dp
    !
    DO i=1,L/2
       !
       k = (pi*dble(2*i-1))/dble(L)
       !
       e  = energy( k, h, gamma )
       e0 = energy( k, h0, gamma0 )
       !
       cgs = cgs * (-4.0_dp*(h0-h)**2 + ( e0 + e )**2  )
       cgs = 0.25_dp*cgs/( e * e0 )
       !
       coeff(i) = ( 4.0_dp*(h0-h)**2-(e0-e)**2 )/( -4.0_dp*(h0-h)**2 + (e0 + e)**2 )
       !
    END DO
    !
    !
    IF (allocated(calpha2)) deallocate(calpha2)
    ALLOCATE( calpha2( 2**(L/2)))
    !
    calpha2(:) = cgs
    control = 0.0_dp
    !
    DO i=0,2**(L/2)-1
       !
       DO j=1,L/2
          IF ( btest(i,j-1) ) THEN
             !
             calpha2(i+1) = calpha2(i+1)*coeff(j)
             !
          END IF
       END DO
       !
       control = control + calpha2(i+1)
       !
    END DO
    !
    PRINT*,"Norm of the ground state: ",control
    PRINT '(x,a,e10.2,a)','estimate use of RAM:',(16.0_dp/1048576.0_dp)*dble(2**(L/2))," Mb"

    DEALLOCATE(coeff)

  END SUBROUTINE calpha2_calc


! instead to calculate and save all the c_aplha and energy
! i calculate the mean value doing a calculation on the fly:
! the result is that a don't use a lot of memory and i don't
! have physical limitation ( apart the time to generate all
! the configurations)
  SUBROUTINE en_dist_calc( nbin )
    IMPLICIT NONE
    !
    REAL(kind = 8) :: e,e0
    REAL(kind = 8) :: cgs, egs, k, mean_en2, delta,ea,ca2
    INTEGER :: i,j,ibin
    INTEGER, INTENT(in) :: nbin
    REAL(kind = 8), DIMENSION (:), ALLOCATABLE :: coeff,evect
    !
    ALLOCATE( coeff(L/2), evect(L/2) )
    !
    IF (allocated(dist)) deallocate(dist)
    ALLOCATE( dist(nbin,2) )
    !
    !k points for an even number of fermions
    cgs = 1.0_dp
    egs = 0.0_dp
    !
    DO i=1,L/2
       !
       k = (pi*dble(2*i-1))/dble(L)
       !
       e  = energy( k, h, gamma )
       e0 = energy( k, h0, gamma0 )
       !
       cgs = cgs * (-4.0_dp*(h0-h)**2 + ( e0 + e )**2  )
       cgs = 0.25_dp*cgs/( e * e0 )
       egs = egs - e/dble(L) !twice because we are doing a summation only over k>0
       !
       coeff(i) = ( 4.0_dp*(h0-h)**2-(e0-e)**2 )/( -4.0_dp*(h0-h)**2 + (e0 + e)**2 )
       evect(i) = e/dble(L)
       !
    END DO
    !
    !
    delta = 2.0*abs(egs)/dble(nbin)

    DO i=1,nbin
       dist( i ,1) = egs + delta * ( dble(i-1) + 0.5 )
    END DO
    !
    dist(:,2) = 0.0d0
    !
    mean_en = 0.0_dp
    mean_en2 = 0.0_dp
    sigma_en = 0.0_dp
    !
    DO i=0,2**(L/2)-1
       !
       ! c^2 of the state alpha
       ! e of the state alpha
       ca2 = cgs
       ea  = egs
       !
       DO j=1,L/2
          IF ( btest(i,j-1) ) THEN
             ea = ea + 2.0_dp*evect(j)
             ca2 = ca2*coeff(j)
          END IF
       END DO
       !
       mean_en = mean_en + ca2*ea
       mean_en2 = mean_en2 + ca2*ea*ea
       !
       ibin = int((ea-egs)/delta)+1
       if ( (ibin <= nbin).and.( ibin > 0 ) ) then
          dist(ibin,2) = dist(ibin,2) + ca2
       end if
       !
    END DO
    !
    sigma_en = sqrt( mean_en2 - mean_en*mean_en)
    dist(:,2) = dist(:,2)/delta
    !
    !
    DEALLOCATE(evect,coeff)

  END SUBROUTINE en_dist_calc
  !
  !
  !
  SUBROUTINE mag_dist_calc( nbin )
    !
    IMPLICIT NONE
    !
    REAL(kind = 8) :: e,e0,mmin
    REAL(kind = 8) :: cgs, mgs, k, delta, ma,ca2, mean_mag,mean_mag2,sigma_mag
    INTEGER :: i,j,ibin
    INTEGER, INTENT(in) :: nbin
    REAL(kind = 8), DIMENSION (:), ALLOCATABLE :: coeff,mvect
    !
    ALLOCATE( coeff(L/2), mvect(L/2) )
    !
    IF (allocated(dist)) deallocate(dist)
    ALLOCATE( dist(nbin,2) )
    !
    !k points for an even number of fermions
    cgs = 1.0_dp
    mgs = 0.0_dp
    mmin = 0.0_dp
    !
    DO i=1,L/2
       !
       k = (pi*dble(2*i-1))/dble(L)
       !
       e  = energy( k, h, gamma )
       e0 = energy( k, h0, gamma0 )
       !
       cgs = cgs * (-4.0_dp*(h0-h)**2 + ( e0 + e )**2  )
       cgs = 0.25_dp*cgs/( e * e0 )
       mvect(i) = sigmaz(k)/dble(L)
       !
       mgs = mgs - mvect(i)
       !
       mmin = mmin - dabs( mvect(i) ) 
       !
       coeff(i) = ( 4.0_dp*(h0-h)**2-(e0-e)**2 )/( -4.0_dp*(h0-h)**2 + (e0 + e)**2 )
       !
    END DO
    !
    !
    !
    !
    delta = 2.0*dabs(mmin)/dble(nbin)
    !
    DO i=1,nbin
       dist( i ,1) = mmin + delta * ( dble(i-1) + 0.5 )
    END DO
    !
    dist(:,2) = 0.0d0
    !
    mean_mag = 0.0_dp
    mean_mag2 = 0.0_dp
    sigma_mag = 0.0_dp
    !
    DO i=0,2**(L/2)-1
       !
       ! c^2 of the state alpha
       ! e of the state alpha
       ca2 = cgs
       ma  = mgs
       !
       DO j=1,L/2
          IF ( btest(i,j-1) ) THEN
             ma = ma + 2.0_dp*mvect(j)
             ca2 = ca2*coeff(j)
          END IF
       END DO
       !
       mean_mag = mean_mag + ca2*ma
       mean_mag2 = mean_mag2 + ca2*ma*ma
       !
       ibin = int((ma-mmin)/delta)+1
       ! ... modifica ...
       if (ibin > nbin) ibin = nbin
       ! ....
       if ( (ibin <= nbin).and.( ibin > 0 ) ) then
          dist(ibin,2) = dist(ibin,2) + ca2
       end if
       !
    END DO
    !
    sigma_mag = sqrt( mean_mag2 - mean_mag*mean_mag)
    sigma_obs = sigma_mag
    dist(:,2) = dist(:,2)/delta
    !
    !
    DEALLOCATE(mvect,coeff)
    !
  END SUBROUTINE mag_dist_calc
  !
  !
  !
  SUBROUTINE kinks_dist_calc( nbin )
    !
    IMPLICIT NONE
    !
    REAL(kind = 8) :: e, e0, kinksmin
    REAL(kind = 8) :: cgs, kinksgs, k, delta, kinksa,ca2,mean_kinks
    INTEGER :: i,j,ibin
    INTEGER, INTENT(in) :: nbin
    REAL(kind = 8), DIMENSION (:), ALLOCATABLE :: coeff, kinksvect
    !
    ALLOCATE( coeff(L/2), kinksvect(L/2) )
    !
    IF (allocated(dist)) deallocate(dist)
    ALLOCATE( dist(nbin,2) )
    !
    !k points for an even number of fermions
    cgs = 1.0_dp
    kinksgs = 0.0_dp
    kinksmin = 0.0_dp
    !
    DO i=1,L/2
       !
       k = (pi*dble(2*i-1))/dble(L)
       !
       e  = energy( k, h, gamma )
       e0 = energy( k, h0, gamma0 )
       !
       cgs = cgs * (-4.0_dp*(h0-h)**2 + ( e0 + e )**2  )
       cgs = 0.25_dp*cgs/( e * e0 )
       kinksvect(i) = kinks(k)/dble(L)
       !
       kinksgs = kinksgs - kinksvect(i)
       !
       kinksmin = kinksmin - dabs( kinksvect(i) ) 
       !
       coeff(i) = ( 4.0_dp*(h0-h)**2-(e0-e)**2 )/( -4.0_dp*(h0-h)**2 + (e0 + e)**2 )
       !
    END DO
    !
    !
    delta = 2.0*dabs(kinksmin)/dble(nbin)
    !
    DO i=1,nbin
       dist( i ,1) = 0.5d0 + kinksmin + delta * ( dble(i-1) + 0.5 )
    END DO
    !
    dist(:,2) = 0.0d0
    !
    mean_kinks = 0.0_dp
    !
    DO i=0,2**(L/2)-1
       !
       ! c^2 of the state alpha
       ! e of the state alpha
       ca2 = cgs
       kinksa  = kinksgs
       !
       DO j=1,L/2
          IF ( btest(i,j-1) ) THEN
             kinksa = kinksa + 2.0_dp*kinksvect(j)
             ca2 = ca2*coeff(j)
          END IF
       END DO
       !
       mean_kinks = mean_kinks + ca2*kinksa
       !
       ibin = int((kinksa-kinksmin)/delta)+1
       if ( (ibin <= nbin).and.( ibin > 0 ) ) then
          dist(ibin,2) = dist(ibin,2) + ca2
       end if
       !
    END DO
    !
    print*, mean_kinks
    !
    dist(:,2) = dist(:,2)/delta
    !
    !
    DEALLOCATE( kinksvect, coeff)
    !
  END SUBROUTINE kinks_dist_calc
  !
  !
  !
  FUNCTION long_time_mag( )
    !
    IMPLICIT NONE
    !
    REAL( kind = dp ) :: k, long_time_mag, e0, e
    INTEGER :: i
    !
    long_time_mag = 0.0d0
    !
    DO i=1,L/2
       !
       k = (pi*dble(2*i-1))/dble(L)
       !
       e0 = energy( k, h0, gamma0 )
       !
       e = energy( k, h, gamma )
       !
       long_time_mag = long_time_mag - 4.0d0* ( dcos(k) - h0 )/e0
       long_time_mag = long_time_mag - 16.0d0*(h0-h)*dsin(k)*dsin(k)/(e*e0*e)
       !
    ENDDO
    !
    long_time_mag = long_time_mag / dble(L)
    !
  END FUNCTION long_time_mag
  !
  !
  FUNCTION long_time_kinks( )
    !
    IMPLICIT NONE
    !
    REAL ( kind = dp ) :: k, long_time_kinks, e0, e
    INTEGER :: i
    !
    long_time_kinks = 0.0d0
    !
    DO i=1,L/2
       !
       k  = (pi*dble(2*i-1))/dble(L)
       !
       e0 = energy( k, h0, gamma0 )
       !
       e = energy( k, h, gamma )
       !
       long_time_kinks = long_time_kinks - 2.0d0*( 1.0d0 - h0*dcos(k) )/e0
       !
       long_time_kinks = long_time_kinks - 8.0d0*h*(h0-h)*dsin(k)*dsin(k)/(e*e0*e)
       !
    ENDDO
    !
    long_time_kinks = long_time_kinks / dble(L)
    long_time_kinks = 0.5d0 + long_time_kinks
    !
  END FUNCTION long_time_kinks
  !
  !
  FUNCTION sigma_mag_quench( ener, deltae, nbin )
    !
    IMPLICIT NONE
    !
    REAL(kind = dp), INTENT(in) :: ener, deltae
    INTEGER, INTENT(in) :: nbin
    !
    REAL(kind = 8), DIMENSION (:), ALLOCATABLE :: evect, mvect
    !
    REAL(kind = 8) :: egs, mgs_loc, k, sigma_mag_quench,delta,mgs_max
    INTEGER :: i, ndata, ibin
    !
    REAL(kind = 8) :: malpha, ealpha
    REAL(kind = 8) :: mean_m, mean2_m
    !
    INTEGER :: j
    !
    allocate( evect(L/2), mvect(L/2) )
    !
    IF ( allocated(dist) ) deallocate( dist )
    allocate( dist(nbin,2) )
    dist = 0.0
    !
    !
    !k points for an even number of fermions
    egs = 0.0_dp
    mgs_loc = 0.0_dp
    mgs_max = 0.0_dp
    !
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
       !
       mgs_max = mgs_max + abs( mvect(i) )
       !
    END DO
    !
    !
    delta = 2.0*mgs_max/dble(nbin)
    !
    !
    ndata = 0
    mean2_m = 0.0d0
    mean_m = 0.0d0
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
       IF ( dabs(ener-ealpha) < deltae )  THEN
          mean_m = mean_m + malpha
          mean2_m = mean2_m + malpha*malpha
          ndata   = ndata + 1
          ibin = int( (malpha + mgs_max)/delta ) + 1
          if (ibin > nbin ) ibin = nbin
          dist(ibin,2) = dist(ibin,2) + 1.0
       END IF
       !
    END DO
    !
    DO i=1,nbin
       dist(i,1) = - mgs_max + ( dble(i) - 0.5d0 ) * delta
    END DO
    !
    mean_m = mean_m / dble(ndata)
    dist(:,2) = dist(:,2)/(delta*dble(ndata))
    !
    print*,"number of points:",ndata
    !
    sigma_mag_quench = dsqrt(mean2_m/dble(ndata) - mean_m*mean_m)
    !
    DEALLOCATE(mvect, evect)
    !
    !
    RETURN
    !
  END FUNCTION sigma_mag_quench
  !
  !
  !
  FUNCTION sigma_mag_quench2( ener, deltae, nbin )
    !
    IMPLICIT NONE
    !
    REAL(kind = dp), INTENT(in) :: ener, deltae
    INTEGER, INTENT(in) :: nbin
    !
    REAL(kind = 8), DIMENSION (:), ALLOCATABLE :: evect, mvect
    !
    REAL(kind = 8) :: egs, mgs_loc, k, sigma_mag_quench2,delta,mgs_max
    INTEGER :: i, ndata, ibin
    !
    REAL(kind = 8) :: malpha, ealpha
    REAL(kind = 8) :: mean_m, mean2_m
    !
    INTEGER :: j
    !
    allocate( evect(L), mvect(L) )
    !
    IF ( allocated(dist) ) deallocate( dist )
    allocate( dist(nbin,2) )
    dist = 0.0
    !
    !
    !k points for an even number of fermions
    egs = 0.0_dp
    mgs_loc = 0.0_dp
    mgs_max = 0.0_dp
    !
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
       !
       mgs_max = mgs_max + 0.5d0*abs( mvect(i) )
       !
    END DO
    !
    !
    delta = 2.0*mgs_max/dble(nbin)
    !
    !
    ndata = 0
    mean2_m = 0.0d0
    mean_m = 0.0d0
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
       !
       IF ( dabs(ener-ealpha) < deltae )  THEN
          mean_m = mean_m + malpha
          mean2_m = mean2_m + malpha*malpha
          ndata   = ndata + 1
          ibin = int( (malpha + mgs_max)/delta ) + 1
          if (ibin > nbin ) ibin = nbin
          dist(ibin,2) = dist(ibin,2) + 1.0
       END IF
       !
    END DO
    !
    DO i=1,nbin
       dist(i,1) = - mgs_max + ( dble(i) - 0.5d0 ) * delta
    END DO
    !
    mean_m = mean_m / dble(ndata)
    dist(:,2) = dist(:,2)/(delta*dble(ndata))
    !
    print*,"number of points:",ndata
    !
    sigma_mag_quench2 = dsqrt(mean2_m/dble(ndata) - mean_m*mean_m)
    !
    DEALLOCATE(mvect, evect)
    !
    !
    RETURN
    !
    !
  END FUNCTION sigma_mag_quench2
  !
  !
  !
  FUNCTION sigma_kink_quench( ener, deltae, nbin )
    !
    IMPLICIT NONE
    !
    REAL(kind = dp), INTENT(in) :: ener, deltae
    INTEGER, INTENT(in) :: nbin
    !
    REAL(kind = 8), DIMENSION (:), ALLOCATABLE :: evect, kinkvect
    !
    REAL(kind = 8) :: egs, kinkgs_loc, k, sigma_kink_quench,delta,kinkgs_max
    INTEGER :: i, ndata, ibin
    !
    REAL(kind = 8) :: kinkalpha, ealpha
    REAL(kind = 8) :: mean_kink, mean2_kink
    !
    INTEGER :: j
    !
    allocate( evect(L/2), kinkvect(L/2) )
    !
    IF ( allocated(dist) ) deallocate( dist )
    allocate( dist(nbin,2) )
    dist = 0.0
    !
    !
    !k points for an even number of fermions
    egs = 0.0_dp
    kinkgs_loc = 0.0_dp
    kinkgs_max = 0.0_dp
    mean_kink = 0.0_dp
    mean2_kink = 0.0_dp
    !
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
       !
       kinkgs_max = kinkgs_max + abs( kinkvect(i) )
       !
    END DO
    !
    !
    delta = 2.0*kinkgs_max/dble(nbin)
    !
    !
    ndata = 0
    !
    DO i=0,2**(L/2)-1
       !
       ealpha = egs
       kinkalpha  = kinkgs_loc
       !
       DO j=1,L/2
          IF ( btest(i,j-1) ) THEN
             !
             kinkalpha = kinkalpha + 2.0_dp*kinkvect(j)
             ealpha = ealpha + 2.0_dp*evect(j)
             !
          END IF
       END DO
       !
       IF ( dabs(ener-ealpha) < deltae )  THEN
          mean_kink = mean_kink + kinkalpha
          mean2_kink = mean2_kink + kinkalpha*kinkalpha
          ndata   = ndata + 1
          ibin = int( (kinkalpha + kinkgs_max)/delta ) + 1
          if (ibin > nbin ) ibin = nbin
          dist(ibin,2) = dist(ibin,2) + 1.0
       END IF
       !
    END DO
    !
    DO i=1,nbin
       dist(i,1) = - kinkgs_max + ( dble(i) - 0.5d0 ) * delta + 0.5d0
    END DO
    !
    mean_kink = mean_kink / dble(ndata)
    dist(:,2) = dist(:,2)/(delta*dble(ndata))
    !
    print*,"number of points:",ndata
    !
    sigma_kink_quench = dsqrt(mean2_kink/dble(ndata) - mean_kink*mean_kink)
    !
    DEALLOCATE(kinkvect, evect)
    !
    !
    RETURN
    !
  END FUNCTION sigma_kink_quench
  !
  !
  !
END MODULE quench
