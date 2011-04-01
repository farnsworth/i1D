MODULE microcanonical
  !
  USE system
  !
  REAL(kind = 8), DIMENSION (:,:), ALLOCATABLE :: dist
  INTEGER :: ndata
  REAL(kind = 8) :: obsmax,obsmin
  !
  ! ... exact microcanonical distribution around given energy ...
  !
CONTAINS
  !
  FUNCTION xxcorrelation_sigma_mc( d, ener, deltae, xxmax, nbin )
    !
    IMPLICIT NONE
    !
    REAL(kind = dp), INTENT(in) :: ener, deltae,xxmax
    INTEGER, INTENT(in) :: nbin,d
    !
    REAL(kind = 8), DIMENSION (:), ALLOCATABLE :: evect
    !
    REAL(kind = 8) :: egs, k ,delta,xxcorrelation_sigma_mc
    INTEGER :: i, j, ibin
    !
    REAL(kind = 8) :: xxalpha, ealpha
    REAL(kind = 8) :: mean_xx, mean2_xx
    !
    LOGICAL, DIMENSION(L/2) :: state
    LOGICAL :: first
    !
    allocate( evect(L/2) )
    !
    IF ( allocated(dist) ) deallocate( dist )
    allocate( dist(nbin,2) )
    dist = 0.0
    !
    obsmax = 0.0d0
    obsmin = 0.0d0
    !
    !k points for an even number of fermions
    egs = 0.0_dp
    !
    DO i=1,L/2
       !
       k = (pi*dble(2*i-1))/dble(L)
       !
       evect(i) = energy( k, h, gamma )/dble(L)
       egs = egs - evect(i)
       !
    END DO
    !
    !
    delta = 2.0*xxmax/dble(nbin)
    !
    !
    ndata = 0
    mean2_xx = 0.0d0
    mean_xx = 0.0d0
    !
    first = .true.
    !
    DO i=0,2**(L/2)-1
       ealpha = egs
       state = .false.
       !
       DO j=1,L/2
          IF ( btest(i,j-1) ) THEN
             state(j) = .true.
             ealpha = ealpha + 2.0_dp*evect(j)
          END IF
       END DO
       !
       IF ( dabs(ener-ealpha) < deltae )  THEN
          !
          IF (first) THEN
             obsmax = xxalpha
             obsmin = xxalpha
             first  = .false.
          ENDIF
          !
          xxalpha = xxcorrelation_red(d,state,L/2)
          IF ( obsmax < xxalpha ) obsmax = xxalpha
          IF ( obsmin > xxalpha ) obsmin = xxalpha
          !
          mean_xx = mean_xx + xxalpha
          mean2_xx = mean2_xx + xxalpha*xxalpha
          ndata   = ndata + 1
          ibin = int( (xxalpha + xxmax)/delta ) + 1
          if (ibin > nbin ) ibin = nbin
          dist(ibin,2) = dist(ibin,2) + 1.0
       END IF
       !
    END DO
    !
    DO i=1,nbin
       dist(i,1) = - xxmax + ( dble(i) - 0.5d0 ) * delta
    END DO
    !
    mean_xx = mean_xx / dble(ndata)
    dist(:,2) = dist(:,2)/(delta*dble(ndata))
    !
    print*,"number of points:",ndata
    !
    xxcorrelation_sigma_mc = dsqrt(mean2_xx/dble(ndata) - mean_xx*mean_xx)
    !
    DEALLOCATE(evect)
    !
    !
    RETURN
    !
  END FUNCTION xxcorrelation_sigma_mc
  !
  !
END MODULE microcanonical
