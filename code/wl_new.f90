!
MODULE wl_new
  !
  USE system
  USE wl, ONLY : read_conf, save_conf, check_obs 
  !
  ! numerical parameters
  REAL (kind=dp) :: f = 2.718281828459_dp, threshold = 0.2_dp, accuracy = 1.0d-6,delta
  !
  REAL (kind=dp), DIMENSION(:,:),ALLOCATABLE :: logdosf
  REAL (kind=dp), DIMENSION(:,:),ALLOCATABLE :: dist
  !
  INTEGER :: ierr
  !
  !
CONTAINS
  !
  ! ... here we use a gaussian instead of a step function
  ! ... probably doesn't have a big gain
  !
  FUNCTION wl_mdos( ener, deltae, mnbin, seed, readconf )
    !
    IMPLICIT NONE
    !
    !
    REAL (kind=8), INTENT(IN) :: ener,deltae
    INTEGER, DIMENSION(8),INTENT(IN) :: seed
    LOGICAL, OPTIONAL,INTENT(IN) :: readconf
    INTEGER,INTENT(IN) :: mnbin
    INTEGER :: wl_mdos
    !
    !
    LOGICAL, DIMENSION(:), ALLOCATABLE :: conf
    INTEGER :: new_conf
    !
    REAL(kind=8), DIMENSION(:), ALLOCATABLE :: logdos
    INTEGER, DIMENSION(:), ALLOCATABLE :: histo
    INTEGER, DIMENSION(:), ALLOCATABLE :: filled
    !
    REAL (kind = 8), DIMENSION(:), ALLOCATABLE :: egrid, mgrid
    REAL (kind = 8) :: econf, new_econf
    REAL (kind = 8) :: mconf, new_mconf
    REAL (kind = 8) :: mmax,mmin,logftemp
    REAL (kind = 8) :: logf,deltam
    !
    INTEGER :: binconf, new_binconf
    INTEGER :: i
    INTEGER :: i_swe, i_swe_tot, step
    INTEGER :: npoints_inside_histo
    !
    !
    wl_mdos = 0
    !
    CALL initialize
    IF (ierr /= 0) THEN
       wl_mdos = ierr
       RETURN
    ENDIF
    !
!    OPEN(unit=2,file=trim(datadir)//'intconf.out',action='write',status='replace')
    !
    !
    i_swe_tot = 0
    step = 0
    !
    DO WHILE ( logf > accuracy )
       !
       step = step + 1
       i_swe = 0
       histo = 0
       npoints_inside_histo = 0
       !
       DO
          i_swe = i_swe + 1
          !
          DO i=1,L
             CALL get_new_conf
             CALL mc_step_1
             CALL update
!             CALL save_intermediate_conf(conf,L)
          END DO
          !
          IF ( mod(i_swe,1000)==0 ) THEN
             !
             CALL check_obs( conf, econf, egrid , L )
             !
             CALL check_obs( conf, mconf, mgrid , L )
             !
             !
             IF (check(histo,mnbin,npoints_inside_histo)) THEN
!             IF (check(histo,mnbin,i_swe*L)) THEN
                EXIT
             END IF
             !
             IF (i_swe>1000000) THEN
                PRINT*,"Some problem in the convergence"
                ierr = 5
                EXIT
             END IF
             !
          END IF
          !
       END DO
       !
       IF (ierr /= 0) GOTO 100
       !
       i_swe_tot = i_swe_tot + i_swe
       !
       IF (mod(step,5) == 0) THEN
          PRINT*,"Step:",step,"logf :",logf
          PRINT*,"Flat histo in",i_swe,"sweps"
          PRINT*
       END IF
       !
       logf = logf/2.0_dp
       !
    END DO
    
    PRINT*,"----------------- RUN IS ENDED -------------------"
    PRINT*
    print*,"Obtained convergence in ",i_swe_tot,"sweps"
    CALL save

100 CALL finalize
    wl_mdos = ierr

    RETURN
    
  CONTAINS
    
    FUNCTION check(histo_in, mnbin_in,npoint_in)
      
      IMPLICIT NONE
      
      LOGICAL :: check
      INTEGER, INTENT(IN) :: mnbin_in , npoint_in
      INTEGER, DIMENSION(mnbin_in), INTENT(IN) :: histo_in
      REAL(kind=8) :: meanhisto
      INTEGER :: i_tmp
      INTEGER :: nbin_tmp

      nbin_tmp = 0
      DO i_tmp=1,mnbin_in
         IF (histo_in(i_tmp) > 0) THEN
            nbin_tmp = nbin_tmp + 1
            IF (filled(i_tmp) == 0) THEN
               filled(i_tmp) = 1
               IF (step/=1) THEN
                  PRINT*,"ERROR : in the first step you did not scan all bins"
                  PRINT*, "--- restart from zero ---"
                  logf = dlog(f)
                  check = .true.
                  GO TO 11
               END IF
            END IF
         ELSE
            IF (filled(i_tmp) == 1) THEN
               check = .false.
               GO TO 11
            END IF
         END IF
      END DO
      
      meanhisto = dble(npoint_in)/dble(nbin_tmp)
      
      check = .true.
      DO i_tmp=1,mnbin_in
         IF (histo_in(i_tmp) /= 0) THEN
            IF (abs(meanhisto-dble(histo_in(i_tmp)))/meanhisto > threshold) THEN
               check = .false.
               EXIT
            END IF
         END IF
      END DO
      
11    CONTINUE

      DO i_tmp=1,mnbin_in
         WRITE(unit=1,fmt=*) histo_in(i_tmp)
      ENDDO
      WRITE(unit=1,fmt=*)
      WRITE(unit=1,fmt=*)
      
    END FUNCTION check
    !
    ! ... non funziona
    !
    SUBROUTINE mc_step_2
      !
      IMPLICIT NONE
      !
      REAL(kind=8) :: rnd_tmp
      REAL(kind=8) :: rate_tmp
      !
      rate_tmp = logdos(binconf) - logdos(new_binconf)
      rate_tmp = dexp(rate_tmp)

      CALL random_number(rnd_tmp)
      IF (rnd_tmp > rate_tmp) GO TO 10
         
      binconf = new_binconf
      econf = new_econf
      mconf = new_mconf
      logftemp = logf*gaussian(new_econf,ener,deltae)
      
      conf( new_conf ) = .not.conf( new_conf )

      
10    RETURN

    END SUBROUTINE mc_step_2
    !
    !
    !
    SUBROUTINE mc_step_1
      !
      IMPLICIT NONE
      !
      REAL(kind=8) :: rnd_tmp
      REAL(kind=8) :: rate_tmp,rate_tmp2
      !
      !
      rate_tmp2 = gaussian(new_econf, ener ,deltae) / gaussian(econf, ener, deltae)
      !CALL random_number(rnd_tmp)
      !
      !IF (rnd_tmp<rate_tmp) THEN
         rate_tmp = logdos(binconf) - logdos(new_binconf)
         rate_tmp = dexp(rate_tmp)
         !
         rate_tmp = rate_tmp*rate_tmp2
         !
         CALL random_number(rnd_tmp)
         !
         ! probabilita' di non accettare
         IF (rnd_tmp > rate_tmp) GO TO 10
         !
         binconf = new_binconf
         econf = new_econf
         mconf = new_mconf
         logftemp = logf
         !
         IF (new_conf /= 0) THEN
            conf( new_conf ) = .not.conf( new_conf )
         ENDIF
      !ENDIF
      !
10    RETURN
      !
    END SUBROUTINE mc_step_1
    !
    !
    !
    SUBROUTINE mc_step_old
      
      IMPLICIT NONE
      
      REAL(kind=8) :: rnd_tmp
      REAL(kind=8) :: rate_tmp
      
      
      IF ( dabs(new_econf-ener) <= deltae ) THEN
         rate_tmp = logdos(binconf) - logdos(new_binconf)
         IF ( rate_tmp < 0.0_dp ) THEN
            rate_tmp = dexp(rate_tmp)
            CALL random_number(rnd_tmp)
            IF (rnd_tmp > rate_tmp) GO TO 10
         END IF
         
         binconf = new_binconf
         econf = new_econf
         mconf = new_mconf
         
!         conf(new_conf(1)) = .not.conf(new_conf(1))
!         conf(new_conf(2)) = .not.conf(new_conf(2))
         conf( new_conf ) = .not.conf( new_conf )
         logftemp = logf
         
      END IF

!      WRITE(unit=2,fmt=*) (dabs(new_econf-ener) <= deltae)
      
10    RETURN
      
    END SUBROUTINE mc_step_old
    
    
    
    SUBROUTINE update
      IMPLICIT NONE
      
      logdos( binconf ) = logdos( binconf ) + logftemp
      ! tentativo: faccio check solo su stati che stanno
      ! entro 3 sigma di energia
      ! togli questo if se non funziona
!      IF ( dabs(econf - ener) < 3.0*deltae ) THEN
         histo ( binconf ) =  histo( binconf ) + 1
         npoints_inside_histo = npoints_inside_histo + 1
!      ENDIF
      
    END SUBROUTINE update

    SUBROUTINE update_1
      IMPLICIT NONE
      
      logdos( binconf ) = logdos( binconf ) + logftemp

      IF ( dabs(econf - ener) < 3.0*deltae ) THEN
         histo ( binconf ) =  histo( binconf ) + 1
         npoints_inside_histo = npoints_inside_histo + 1
      ENDIF
      
    END SUBROUTINE update_1
    
    
    
    SUBROUTINE get_new_conf
      IMPLICIT NONE
      
      REAL(kind=8) :: rnd_tmp
      
      new_econf = econf
      new_mconf = mconf


      CALL random_number( rnd_tmp )
      new_conf = int(dble(L)*rnd_tmp) + 1
      
      IF (conf( new_conf )) THEN
         new_econf = new_econf - egrid( new_conf )
         new_mconf = new_mconf - mgrid( new_conf )
      ELSE
         new_econf = new_econf + egrid( new_conf )
         new_mconf = new_mconf + mgrid( new_conf )
      ENDIF
      
      new_binconf = int( (new_mconf - mmin)/deltam ) + 1
      if (new_binconf > mnbin) new_binconf = mnbin 
      
    END SUBROUTINE get_new_conf
    !
    !
    SUBROUTINE get_new_conf_1
      IMPLICIT NONE
      
      REAL(kind=8) :: rnd_tmp,rate_tmp
      
      new_econf = econf
      new_mconf = mconf


      CALL random_number( rnd_tmp )
      new_conf = int(dble(L)*rnd_tmp) + 1
      
      IF (conf( new_conf )) THEN
         new_econf = new_econf - egrid( new_conf )
         new_mconf = new_mconf - mgrid( new_conf )
      ELSE
         new_econf = new_econf + egrid( new_conf )
         new_mconf = new_mconf + mgrid( new_conf )
      ENDIF

      CALL random_number( rnd_tmp )
      rate_tmp = gaussian(new_econf, ener ,deltae) / gaussian(econf, ener, deltae)
      !
      IF (rnd_tmp<rate_tmp) THEN
         new_econf = econf
         new_mconf = mconf
         new_conf = 0
      ENDIF
      
      new_binconf = int( (new_mconf - mmin)/deltam ) + 1
      if (new_binconf > mnbin) new_binconf = mnbin 
      
    END SUBROUTINE get_new_conf_1
    !
    !
    !--------------------------------------------------------------------!
    SUBROUTINE initialize
      !
      IMPLICIT NONE
      !
      REAL(kind=8) :: delta_tmp,rnd_tmp,emax_tmp
      REAL(kind=8) :: k_tmp
      INTEGER :: i_tmp,j_tmp,l1_tmp,l2_tmp, l0_tmp
      LOGICAL :: readconf_tmp
      !
      ierr = 0
      !
      IF (mod(L,2) /= 0) THEN
         CALL error('initialize','only even L are admitted')
         ierr = 1
         RETURN
      END IF
      !
      logf = dlog(f)
      !
      IF ( deltae <= 0.0d0) THEN
         CALL error('initialize','negative deltae')
         ierr = 1
         RETURN
      END IF
      !
      CALL random_seed(put=seed)
      !
      ALLOCATE(conf(L))
      ALLOCATE(egrid(L))
      ALLOCATE(mgrid(L))
      !
      delta_tmp = 2.0_dp*pi/dble(L)
      
      k_tmp = pi*(-1.0_dp-1.0_dp/dble(L))
      !
      !
      emax_tmp = egrid(1)
      !
      mmax = 0.0_dp
      mmin = 0.0_dp
      mconf = 0.0_dp
      econf = 0.0_dp
      !
      DO i_tmp=1,L
         !
         k_tmp = k_tmp+delta_tmp
         !energy associated to point k_tmp
         egrid(i_tmp) = energy(k_tmp, h, gamma )
         mgrid(i_tmp) = sigmaz(k_tmp)
         !
         !calculation of max sigmaz in order to calculare range of mean value
         IF ( mgrid(i_tmp) > 0.0_dp) mmax = mmax + mgrid(i_tmp)
         IF ( mgrid(i_tmp) < 0.0_dp) mmin = mmin + mgrid(i_tmp)
         !
         !ground state energy
         econf = econf - 0.5_dp*egrid(i_tmp)
         ! ground state magnetization
         mconf = mconf - 0.5_dp*mgrid(i_tmp)
         !
         ! i compute the maximum energy
         IF (i==1) THEN
            emax_tmp = egrid(1)
         ELSE
            IF ( emax_tmp < egrid(i_tmp) )  emax_tmp = egrid(i_tmp)
         ENDIF
         !
      END DO
      !
      conf = .false.
      !
      mmin = (mconf+mmin) / dble(L)
      mmax = (mconf+mmax)/dble(L)
      !
      !
      deltam = (mmax-mmin)/dble(mnbin)
      !
      !
      !check if emax and emin are compatible with gs energy
      IF ((ener+deltae < econf/dble(L) ).or.(ener-deltae > -econf/dble(L) )) THEN
         PRINT*,"gs energy",econf/dble(L)
         CALL error("initialize", "energy interval outside &
              &energy spectra" )
         ierr = 1
         RETURN
      END IF
      !
      !
      IF (  dabs( econf/dble(L)-ener) <= deltae ) THEN
         GO TO 10
      END IF
      !
      readconf_tmp = .false.
      !
      IF (present(readconf)) THEN
         readconf_tmp = readconf
      ENDIF
      !
      IF (readconf_tmp) THEN
         !
         CALL read_conf(conf,L)
         DO l1_tmp = 1, L
            IF ( conf(l1_tmp) ) THEN
               econf = econf + egrid(l1_tmp)
               mconf = mconf + mgrid(l1_tmp)
            ENDIF
         END DO
         IF (  dabs( econf/dble(L)-ener) > deltae ) THEN
            ierr = 2
            RETURN
         END IF
         !
         GOTO 10
         !
      END IF
      !
      !
      !
      !check if i can always find a good initial state
      IF ( (2.0d0*deltae)<(emax_tmp/dble(L)) ) THEN
         PRINT*,"Warning:"," i am not sure to get an initial state"
         PRINT*,"Window :",2.0d0*deltae
         PRINT*,"Maximum gap:",(2.0*emax_tmp/dble(L))
         ierr = -1
         RETURN
      END IF
      !
      !
      l0_tmp = 0
      !
      !
      DO
         CALL random_number(rnd_tmp)
         j_tmp = int(rnd_tmp*(L-l0_tmp)) + 1
         IF (j_tmp>L-l0_tmp) j_tmp = L
         !
         l2_tmp = 0
         !
         DO l1_tmp = 1, L
            !
            IF ( .not.conf(l1_tmp) ) THEN
               l2_tmp = l2_tmp + 1
               IF (l2_tmp == j_tmp) THEN
                  conf(l1_tmp) = .true.
                  econf = econf + egrid(l1_tmp)
                  mconf = mconf + mgrid(l1_tmp)
                  l0_tmp = l0_tmp + 1
                  print*,"econf :",econf/dble(L)
                  IF ( dabs( econf/dble(L) - ener ) <= deltae ) THEN
                     CALL save_conf(conf,L)
                     GO TO 10
                  ENDIF
                  GO TO 11
               END IF
            END IF
            !
11          CONTINUE
            !
         END DO
         !
         IF (l0_tmp == L) EXIT
         !
      END DO
      !
      CALL error("intialize","I have not been able to generate the initial configuration")
      ierr = 1
      RETURN
      !
10    CONTINUE
      !
      PRINT*,"finded configuration with energy",econf/dble(L)
      !
      ! energy per site
      egrid = egrid/dble(L)
      econf = econf/dble(L)
      mgrid = mgrid/dble(L)
      mconf = mconf/dble(L)
      !
      !
      !
      !
      ALLOCATE(histo(mnbin))
      ALLOCATE(filled(mnbin))
      ALLOCATE(logdos(mnbin))
      !
      histo = 0
      filled = 0
      logdos = 0.0_dp
      !
      binconf = int((mconf-mmin) / deltam )+1
      IF (binconf > mnbin) binconf = mnbin 
      !
      !
      IF ((binconf<1).or.(binconf>mnbin)) THEN
         CALL error("initialize", "magnetization outside its interval" )
         ierr = 1
         RETURN
      END IF
      !
      PRINT*
      PRINT*,"----------------- RUN IS STARTED -----------------"
      !
!      CALL check_obs( conf, econf, egrid , L )
!      CALL check_obs( conf, mconf, mgrid , L )
      !
      RETURN
      !
      !
    END SUBROUTINE initialize
    !
    !
    SUBROUTINE finalize
      IMPLICIT NONE
      
      !    CLOSE(unit=2)
      delta = deltam
      DEALLOCATE(conf,egrid,logdos,histo,filled)
      
    END SUBROUTINE finalize
    !
    !
    SUBROUTINE save
      !
      IMPLICIT NONE
      !
      INTEGER :: i_tmp,neffbin,ieff
      REAL(kind=8) :: alpha_tmp
      !
      neffbin = 0
      alpha_tmp = -1.0
      DO i_tmp=1,mnbin
         IF (filled(i_tmp) == 1) THEN
            neffbin = neffbin + 1
            IF (alpha_tmp < 0) alpha_tmp = logdos( i_tmp )
         END IF
      END DO
      !
      IF (allocated(logdosf)) deallocate(logdosf)
      !
      ALLOCATE( logdosf( neffbin, 2 ) )
      !
      ieff = 0
      DO i_tmp=1,mnbin
         IF (filled(i_tmp) == 1) THEN
            ieff = ieff + 1
            logdosf( ieff, 1) = mmin+(dble(i_tmp)-0.5)*deltam
            logdosf( ieff, 2)  = logdos(i_tmp)-alpha_tmp
         END IF
      END DO
      !
    END SUBROUTINE save
    !
  END FUNCTION wl_mdos
  !
  !
  ! ... exact calculation ...
  !
  SUBROUTINE exact( ener, deltae, nbin )
    !
    IMPLICIT NONE
    !
    REAL(kind = dp), INTENT(in) :: ener, deltae
    INTEGER, INTENT(in) :: nbin
    !
    REAL(kind = 8), DIMENSION (:), ALLOCATABLE :: evect, mvect
    !
    REAL(kind = 8) :: egs, mgs_loc, k,delta,mgs_max,ndata,weight
    INTEGER :: i, ibin
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
    ndata = 0.0d0
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
       weight = gaussian(ealpha,ener,deltae)
       !
       mean_m = mean_m + weight * malpha
       mean2_m = mean2_m + weight * malpha * malpha
       ndata   = ndata + weight
       ibin = int( (malpha + mgs_max)/delta ) + 1
       if (ibin > nbin ) ibin = nbin
       dist(ibin,2) = dist(ibin,2) + weight
       !
    END DO
    !
    DO i=1,nbin
       dist(i,1) = - mgs_max + ( dble(i) - 0.5d0 ) * delta
    END DO
    !
    mean_m = mean_m / ndata
    mean2_m = mean2_m / ndata
    dist(:,2) = dist(:,2)/(delta*ndata)
    !
    !
    DEALLOCATE(mvect, evect)
    !
    !
    RETURN
    !
    !
  END SUBROUTINE exact
  !
  !
  ! ... P_D distributions
  !
  FUNCTION wl_ediag( enbins, seed )
    !
    IMPLICIT NONE
    !
    !
    INTEGER, DIMENSION(8),INTENT(IN) :: seed
    INTEGER,INTENT(IN) :: enbins
    INTEGER :: wl_ediag
    !
    !
    LOGICAL, DIMENSION(:), ALLOCATABLE :: conf
    INTEGER :: new_conf
    !
    REAL(kind=8), DIMENSION(:), ALLOCATABLE :: logdos
!    INTEGER, DIMENSION(:), ALLOCATABLE :: histo
    REAL(kind=8), DIMENSION(:), ALLOCATABLE :: histo
    INTEGER, DIMENSION(:), ALLOCATABLE :: filled
    !
    REAL (kind = 8), DIMENSION(:), ALLOCATABLE :: egrid, coeff
    REAL (kind = 8) :: econf, new_econf, emax,emin
    REAL (kind = 8) :: logftemp
    REAL (kind = 8) :: logf,deltae,normalization
    !
    INTEGER :: binconf, new_binconf
    INTEGER :: i
    INTEGER :: i_swe, i_swe_tot, step
    !
    REAL (kind = 8) :: mean_value
    REAL (kind = 8) :: calpha,norm
    !
    !
    wl_ediag = 0
    !
    CALL initialize
    IF (ierr /= 0) THEN
       wl_ediag = ierr
       RETURN
    ENDIF
    !
    !
    i_swe_tot = 0
    step = 0
    !
    DO WHILE ( logf > accuracy )
       !
       step = step + 1
       i_swe = 0
       histo = 0.0d0
       normalization = 0.0d0
       mean_value = 0.0
       norm = 0.0d0
       !
       DO
          i_swe = i_swe + 1
          !
          DO i=1,L
             CALL get_new_conf
             CALL mc_step
             CALL update
          END DO
          !
          !
          IF (ierr /= 0) THEN
             wl_ediag = ierr
             RETURN
          ENDIF
          !
          IF ( mod(i_swe,1000)==0 ) THEN
             !
             CALL check_obs( conf, econf, egrid , L/2 )
             !
             !
             IF (check(histo,enbins,i_swe*L)) THEN
                print*,"mean value:",mean_value/norm
                EXIT
             END IF
             !
             IF (i_swe>1000000) THEN
                PRINT*,"Some problem in the convergence"
                ierr = 5
                EXIT
             END IF
             !
          END IF
          !
       END DO
       !
       IF (ierr /= 0) GOTO 100
       !
       i_swe_tot = i_swe_tot + i_swe
       !
       IF (mod(step,5) == 0) THEN
          PRINT*,"Step:",step,"logf :",logf
          PRINT*,"Flat histo in",i_swe,"sweps"
          PRINT*
       END IF
       !
       logf = logf/2.0_dp
       !
    END DO
    
    PRINT*,"----------------- RUN IS ENDED -------------------"
    PRINT*
    print*,"Obtained convergence in ",i_swe_tot,"sweps"
    CALL save

100 CALL finalize
    wl_ediag = ierr

    RETURN
    
  CONTAINS
    
    FUNCTION check(histo_in, nbins_in, npoints_in)
      
      IMPLICIT NONE
      
      LOGICAL :: check
      INTEGER, INTENT(IN) :: nbins_in , npoints_in
      REAL(kind=8), DIMENSION(nbins_in), INTENT(IN) :: histo_in
      REAL(kind=8) :: meanhisto
      INTEGER :: i_tmp
      INTEGER :: nbins_tmp
      !
      !
      nbins_tmp = 0
      DO i_tmp=1,nbins_in
         IF (histo_in(i_tmp) > 0) THEN
            nbins_tmp = nbins_tmp + 1
            IF (filled(i_tmp) == 0) THEN
               filled(i_tmp) = 1
               IF (step/=1) THEN
                  PRINT*,"ERROR : in the first step you did not scan all bins"
                  PRINT*, "--- restart from zero ---"
                  logf = dlog(f)
                  check = .true.
                  GO TO 11
               END IF
            END IF
         ELSE
            IF (filled(i_tmp) == 1) THEN
               check = .false.
               GO TO 11
            END IF
         END IF
      END DO
      !
!      meanhisto = dble(npoints_in)/dble(nbins_tmp)
      meanhisto = dble(normalization)/dble(nbins_tmp)
      !
      check = .true.
      DO i_tmp=1,nbins_in
         IF (histo_in(i_tmp) /= 0) THEN
            IF (abs(meanhisto-dble(histo_in(i_tmp)))/meanhisto > threshold) THEN
               check = .false.
               EXIT
            END IF
         END IF
      END DO
      
11    CONTINUE
      
    END FUNCTION check
    !
    !
    !
    SUBROUTINE mc_step
      !
      IMPLICIT NONE
      !
      REAL(kind=8) :: rnd_tmp
      REAL(kind=8) :: rate_tmp,rate_tmp2
      !
      rate_tmp2 = 1.0d0
      !
      IF ( conf( new_conf ) ) THEN
         rate_tmp2 = rate_tmp2/coeff( new_conf )
      ELSE
         rate_tmp2 = rate_tmp2*coeff( new_conf )
      ENDIF
      !
      rate_tmp = logdos(binconf) - logdos(new_binconf)
      rate_tmp = dexp(rate_tmp)
      !
      rate_tmp = rate_tmp*rate_tmp2
      !
      CALL random_number(rnd_tmp)
      !
      ! probabilita' di non accettare
      IF (rnd_tmp > rate_tmp) GO TO 10
      !
      binconf = new_binconf
      econf = new_econf
      logftemp = logf
      !
      IF ( conf( new_conf ) ) THEN
         calpha = calpha/coeff( new_conf )
      ELSE
         calpha = calpha*coeff( new_conf )
      ENDIF

      !
      conf( new_conf ) = .not.conf( new_conf )
      !
10    RETURN
      !
    END SUBROUTINE mc_step
    !
    !
    SUBROUTINE update
      IMPLICIT NONE
      !
!      INTEGER :: i_tmp
      REAL (kind=8) :: c2alpha
      !
      c2alpha = 1.0d0
      !
!      DO i_tmp = 1,L/2
!         IF (conf(i_tmp)) THEN
!            c2alpha = c2alpha*coeff(i_tmp)
!         ENDIF
!      ENDDO
      !
      normalization = normalization + 1.0d0/c2alpha
      !
      logdos( binconf ) = logdos( binconf ) + logftemp
      histo ( binconf ) =  histo( binconf ) + 1.0d0/c2alpha
      !
      mean_value = mean_value + calpha * econf
      norm = norm + calpha
      !
    END SUBROUTINE update
    !
    !
    SUBROUTINE get_new_conf
      IMPLICIT NONE
      !
      REAL(kind=8) :: rnd_tmp
      !
      new_econf = econf
      !
      CALL random_number( rnd_tmp )
      new_conf = int(dble(L/2)*rnd_tmp) + 1
      !
      IF ( new_conf == L/2+1 ) new_conf = L/2 
      !
      IF (conf( new_conf )) THEN
         new_econf = new_econf - egrid( new_conf )
      ELSE
         new_econf = new_econf + egrid( new_conf )
      ENDIF
      !
      !
      new_binconf = int( (new_econf - emin)/deltae ) + 1
      !
      IF (new_binconf > enbins) THEN
         new_binconf = enbins
      ENDIF
      !
      RETURN
      !
    END SUBROUTINE get_new_conf
    !
    !
    !
    !--------------------------------------------------------------------!
    SUBROUTINE initialize
      !
      IMPLICIT NONE
      !
      REAL(kind=8) :: delta_tmp,rnd_tmp
      REAL(kind=8) :: k_tmp, e0_tmp
      INTEGER :: i_tmp,j_tmp,l_tmp,l1_tmp
      !
      ierr = 0
      !
      IF (mod(L,2) /= 0) THEN
         CALL error('initialize','only even L are admitted')
         ierr = 1
         RETURN
      END IF
      !
      logf = dlog(f)
      !
      CALL random_seed(put=seed)
      !
      !
      ALLOCATE(conf(L/2))
      ALLOCATE(egrid(L/2))
      ALLOCATE(coeff(L/2))
      !
      !
      delta_tmp = 2.0_dp*pi/dble(L)
      k_tmp = -pi/dble(L)
      !
      !
      econf = 0.0_dp
      coeff = 0.0_dp
      !
      DO i_tmp=1,L/2
         !
         k_tmp = k_tmp+delta_tmp
         egrid(i_tmp) = energy(k_tmp, h, gamma )
         !
         !ground state energy
         econf = econf - egrid(i_tmp)
         !
         e0_tmp = energy(k_tmp, h0, gamma )
         coeff(i_tmp) = ( 4.0_dp* ( h0-h )**2 - ( e0_tmp - egrid(i_tmp) )**2 ) /&
              ( -4.0_dp*(h0-h)**2  +  ( e0_tmp + egrid(i_tmp) )**2 )
         !
      END DO
      !
      !
      conf = .false.
      !
      emin = econf/dble(L)
      emax = -econf/dble(L)
      deltae = (emax-emin)/dble(enbins)
      ! because each element corresponds to two quasiparticles
      egrid = 2.0*egrid 
      !
      !
      calpha = 1.0d0
      !
      DO i_tmp = 1,L/4
         !
         CALL random_number(rnd_tmp)
         j_tmp = int(rnd_tmp*( L/2 - i_tmp + 1 )) + 1
         IF (j_tmp > (L/2-i_tmp + 1) ) j_tmp = L/2 - i_tmp + 1
         !
         l_tmp = 0
         !
         DO l1_tmp=1,L/2
            !
            IF ( .not.conf(l1_tmp) ) THEN
               l_tmp = l_tmp + 1
               IF (l_tmp == j_tmp) THEN
                  conf(l1_tmp) = .true.
                  econf = econf + egrid(l1_tmp) 
                  calpha = calpha*coeff(l1_tmp)
                  GO TO 11
               END IF
            END IF
            !
11          CONTINUE
            !
         END DO
         !
         !
      END DO
      !
      !
      PRINT*," energy",econf/dble(L)
      !
      ! energy per site
      egrid = egrid/dble(L)
      econf = econf/dble(L)
      !
      ALLOCATE(histo(enbins))
      ALLOCATE(filled(enbins))
      ALLOCATE(logdos(enbins))
      !
      histo = 0
      filled = 0
      logdos = 0.0_dp
      !
      binconf = int( (econf-emin) / deltae )+1
      IF (binconf > enbins) binconf = enbins 
      !
      !
      IF ((binconf<1).or.(binconf>enbins)) THEN
         CALL error("initialize", "energy outside its interval" )
         ierr = 1
         RETURN
      END IF
      !
      PRINT*
      PRINT*,"----------------- RUN IS STARTED -----------------"
      !
      !
      CALL check_obs( conf, econf, egrid , L/2 )
      !
      !
      RETURN
      !
      !
    END SUBROUTINE initialize
    !
    !
    SUBROUTINE finalize
      IMPLICIT NONE
      
      delta = deltae
      DEALLOCATE(conf,egrid,logdos,histo,filled)
      
    END SUBROUTINE finalize
    !
    !
    SUBROUTINE save
      !
      IMPLICIT NONE
      !
      INTEGER :: i_tmp,neffbin,ieff
      REAL(kind=8) :: alpha_tmp
      !
      ! ... decidere come salvare ... !
      !
      neffbin = 0
      alpha_tmp = -1.0
      DO i_tmp=1,enbins
         IF (filled(i_tmp) == 1) THEN
            neffbin = neffbin + 1
            IF (alpha_tmp < 0) alpha_tmp = logdos( i_tmp )
         END IF
      END DO
      !
      IF (allocated(logdosf)) deallocate(logdosf)
      !
      ALLOCATE( logdosf( neffbin, 2 ) )
      !
      ieff = 0
      DO i_tmp=1,enbins
         IF (filled(i_tmp) == 1) THEN
            ieff = ieff + 1
            logdosf( ieff, 1) = emin+(dble(i_tmp)-0.5)*deltae
            logdosf( ieff, 2)  = logdos(i_tmp)-alpha_tmp
         END IF
      END DO
      !
    END SUBROUTINE save
    !
  END FUNCTION wl_ediag
  !
  !
  !
  FUNCTION wl_mdiag( mnbins, seed )
    !
    IMPLICIT NONE
    !
    INTEGER, DIMENSION(8),INTENT(IN) :: seed
    INTEGER,INTENT(IN) :: mnbins
    INTEGER :: wl_mdiag
    !
    LOGICAL, DIMENSION(:), ALLOCATABLE :: conf
    INTEGER :: new_conf
    !
    REAL(kind=8), DIMENSION(:), ALLOCATABLE :: logdos
    INTEGER, DIMENSION(:), ALLOCATABLE :: histo
    INTEGER, DIMENSION(:), ALLOCATABLE :: filled
    !
    REAL (kind = 8), DIMENSION(:), ALLOCATABLE :: mgrid, coeff
    REAL (kind = 8) :: mconf, new_mconf, mmax,mmin
    REAL (kind = 8) :: logftemp
    REAL (kind = 8) :: logf,deltam
    !
    INTEGER :: binconf, new_binconf
    INTEGER :: i
    INTEGER :: i_swe, i_swe_tot, step
    !
    !
    wl_mdiag = 0
    !
    CALL initialize
    IF (ierr /= 0) THEN
       wl_mdiag = ierr
       RETURN
    ENDIF
    !
    !
    i_swe_tot = 0
    step = 0
    !
    DO WHILE ( logf > accuracy )
       !
       step = step + 1
       i_swe = 0
       histo = 0
       !
       DO
          i_swe = i_swe + 1
          !
          DO i=1,L
             CALL get_new_conf
             CALL mc_step
             CALL update
          END DO
          !
          !
          IF (ierr /= 0) THEN
             wl_mdiag = ierr
             RETURN
          ENDIF
          !
          IF ( mod(i_swe,1000)==0 ) THEN
             !
             CALL check_obs( conf, mconf, mgrid , L/2 )
             !
             !
             IF (check(histo,mnbins,i_swe*L)) THEN
                EXIT
             END IF
             !
             IF (i_swe>1000000) THEN
                PRINT*,"Some problem in the convergence"
                ierr = 5
                EXIT
             END IF
             !
          END IF
          !
       END DO
       !
       IF (ierr /= 0) GOTO 100
       !
       i_swe_tot = i_swe_tot + i_swe
       !
       IF (mod(step,5) == 0) THEN
          PRINT*,"Step:",step,"logf :",logf
          PRINT*,"Flat histo in",i_swe,"sweps"
          PRINT*
       END IF
       !
       logf = logf/2.0_dp
       !
    END DO
    
    PRINT*,"----------------- RUN IS ENDED -------------------"
    PRINT*
    print*,"Obtained convergence in ",i_swe_tot,"sweps"
    CALL save

100 CALL finalize
    wl_mdiag = ierr
    !
    RETURN
    !
  CONTAINS
    !
    FUNCTION check(histo_in, nbins_in, npoints_in)
      !
      IMPLICIT NONE
      !
      LOGICAL :: check
      INTEGER, INTENT(IN) :: nbins_in , npoints_in
      INTEGER, DIMENSION(nbins_in), INTENT(IN) :: histo_in
      REAL(kind=8) :: meanhisto
      INTEGER :: i_tmp
      INTEGER :: nbins_tmp
      !
      !
      nbins_tmp = 0
      DO i_tmp=1,nbins_in
         IF (histo_in(i_tmp) > 0) THEN
            nbins_tmp = nbins_tmp + 1
            IF (filled(i_tmp) == 0) THEN
               filled(i_tmp) = 1
               IF (step/=1) THEN
                  PRINT*,"ERROR : in the first step you did not scan all bins"
                  PRINT*, "--- restart from zero ---"
                  logf = dlog(f)
                  check = .true.
                  GO TO 11
               END IF
            END IF
         ELSE
            IF (filled(i_tmp) == 1) THEN
               check = .false.
               GO TO 11
            END IF
         END IF
      END DO
      !
      meanhisto = dble(npoints_in)/dble(nbins_tmp)
      !
      check = .true.
      DO i_tmp=1,nbins_in
         IF (histo_in(i_tmp) /= 0) THEN
            IF (abs(meanhisto-dble(histo_in(i_tmp)))/meanhisto > threshold) THEN
               check = .false.
               EXIT
            END IF
         END IF
      END DO
      
11    CONTINUE
      
    END FUNCTION check
    !
    !
    !
    SUBROUTINE mc_step
      !
      IMPLICIT NONE
      !
      REAL(kind=8) :: rnd_tmp
      REAL(kind=8) :: rate_tmp,rate_tmp2
      !
      rate_tmp2 = 1.0d0
      !
      IF ( conf( new_conf ) ) THEN
         rate_tmp2 = rate_tmp2/coeff( new_conf )
      ELSE
         rate_tmp2 = rate_tmp2*coeff( new_conf )
      ENDIF
      !
      rate_tmp = logdos(binconf) - logdos(new_binconf)
      rate_tmp = dexp(rate_tmp)
      !
      rate_tmp = rate_tmp*rate_tmp2
      !
      CALL random_number(rnd_tmp)
      !
      ! probabilita' di non accettare
      IF (rnd_tmp > rate_tmp) GO TO 10
      !
      binconf = new_binconf
      mconf = new_mconf
      logftemp = logf
      !
      conf( new_conf ) = .not.conf( new_conf )
      !
10    RETURN
      !
    END SUBROUTINE mc_step
    !
    !
    SUBROUTINE update
      IMPLICIT NONE
      !
      !
      logdos( binconf ) = logdos( binconf ) + logftemp
      histo ( binconf ) =  histo( binconf ) + 1
      !
    END SUBROUTINE update
    !
    !
    SUBROUTINE get_new_conf
      IMPLICIT NONE
      !
      REAL(kind=8) :: rnd_tmp
      !
      new_mconf = mconf
      !
      CALL random_number( rnd_tmp )
      new_conf = int(dble(L/2)*rnd_tmp) + 1
      !
      IF ( new_conf == L/2+1 ) new_conf = L/2 
      !
      IF (conf( new_conf )) THEN
         new_mconf = new_mconf - mgrid( new_conf )
      ELSE
         new_mconf = new_mconf + mgrid( new_conf )
      ENDIF
      !
      !
      new_binconf = int( (new_mconf - mmin)/deltam ) + 1
      !
      IF (new_binconf > mnbins) THEN
         new_binconf = mnbins
      ENDIF
      !
      RETURN
      !
    END SUBROUTINE get_new_conf
    !
    !
    !
    !--------------------------------------------------------------------!
    SUBROUTINE initialize
      !
      IMPLICIT NONE
      !
      REAL(kind=8) :: delta_tmp,rnd_tmp
      REAL(kind=8) :: k_tmp, e0_tmp, e_tmp
      INTEGER :: i_tmp,j_tmp,l_tmp,l1_tmp
      !
      ierr = 0
      !
      IF (mod(L,2) /= 0) THEN
         CALL error('initialize','only even L are admitted')
         ierr = 1
         RETURN
      END IF
      !
      logf = dlog(f)
      !
      CALL random_seed(put=seed)
      !
      !
      ALLOCATE(conf(L/2))
      ALLOCATE(mgrid(L/2))
      ALLOCATE(coeff(L/2))
      !
      !
      delta_tmp = 2.0_dp*pi/dble(L)
      k_tmp = -pi/dble(L)
      !
      !
      mconf = 0.0_dp
      coeff = 0.0_dp
      mmin = 0.0_dp
      !
      DO i_tmp=1,L/2
         !
         k_tmp = k_tmp+delta_tmp
         mgrid(i_tmp) = sigmaz( k_tmp )
         !
         !ground state energy
         mconf = mconf - mgrid(i_tmp)
         !
         e_tmp = energy(k_tmp, h , gamma )
         e0_tmp = energy(k_tmp, h0, gamma )
         !
         coeff(i_tmp) = ( 4.0_dp * ( h0-h )**2 - ( e0_tmp - e_tmp )**2 ) /&
              ( -4.0_dp * ( h0 - h )**2  +  ( e0_tmp + e_tmp )**2 )
         !
         mmin = mmin - dabs( mgrid(i_tmp) )
         !
      END DO
      !
      !
      conf = .false.
      !
      mmin = mmin/dble(L)
      mmax = -mmin
      !
      !
      deltam = (mmax-mmin)/dble(mnbins)
      ! because each element corresponds to two quasiparticles
      mgrid = 2.0*mgrid
      !
      !
      DO i_tmp = 1,L/4
         !
         CALL random_number(rnd_tmp)
         j_tmp = int(rnd_tmp*( L/2 - i_tmp + 1 )) + 1
         IF (j_tmp > (L/2-i_tmp + 1) ) j_tmp = L/2 - i_tmp + 1
         !
         l_tmp = 0
         !
         DO l1_tmp=1,L/2
            !
            IF ( .not.conf(l1_tmp) ) THEN
               l_tmp = l_tmp + 1
               IF (l_tmp == j_tmp) THEN
                  conf(l1_tmp) = .true.
                  mconf = mconf + mgrid(l1_tmp) 
                  GO TO 11
               END IF
            END IF
            !
11          CONTINUE
            !
         END DO
         !
         !
      END DO
      !
      !
      PRINT*," magnetization",mconf/dble(L)
      !
      ! energy per site
      mgrid = mgrid/dble(L)
      mconf = mconf/dble(L)
      !
      ALLOCATE(histo(mnbins))
      ALLOCATE(filled(mnbins))
      ALLOCATE(logdos(mnbins))
      !
      histo = 0
      filled = 0
      logdos = 0.0_dp
      !
      binconf = int( (mconf-mmin) / deltam )+1
      IF (binconf > mnbins) binconf = mnbins 
      !
      !
      IF ((binconf<1).or.(binconf>mnbins)) THEN
         CALL error("initialize", "energy outside its interval" )
         ierr = 1
         RETURN
      END IF
      !
      PRINT*
      PRINT*,"----------------- RUN IS STARTED -----------------"
      !
      !
      CALL check_obs( conf, mconf, mgrid , L/2 )
      !
      !
      RETURN
      !
      !
    END SUBROUTINE initialize
    !
    !
    SUBROUTINE finalize
      IMPLICIT NONE
      
      delta = deltam
      DEALLOCATE(conf,mgrid,logdos,histo,filled)
      
    END SUBROUTINE finalize
    !
    !
    SUBROUTINE save
      !
      IMPLICIT NONE
      !
      INTEGER :: i_tmp,neffbin,ieff
      REAL(kind=8) :: alpha_tmp
      !
      ! ... decidere come salvare ... !
      !
      neffbin = 0
      alpha_tmp = -1.0
      DO i_tmp=1,mnbins
         IF (filled(i_tmp) == 1) THEN
            neffbin = neffbin + 1
            IF (alpha_tmp < 0) alpha_tmp = logdos( i_tmp )
         END IF
      END DO
      !
      IF (allocated(logdosf)) deallocate(logdosf)
      !
      ALLOCATE( logdosf( neffbin, 2 ) )
      !
      ieff = 0
      DO i_tmp=1,mnbins
         IF (filled(i_tmp) == 1) THEN
            ieff = ieff + 1
            logdosf( ieff, 1) = mmin+(dble(i_tmp)-0.5)*deltam
            logdosf( ieff, 2)  = logdos(i_tmp)-alpha_tmp
         END IF
      END DO
      !
    END SUBROUTINE save
    !
  END FUNCTION wl_mdiag
  !
  !
  !
  FUNCTION gaussian(x,mu,sigma)
    !
    IMPLICIT NONE
    !
    REAL(kind = dp) :: x,mu,sigma,gaussian
    !
    gaussian = (x-mu)*(x-mu)/(2.0_dp*sigma*sigma)
    !
    gaussian = dexp(-gaussian)/(sigma*dsqrt(2.0_dp*pi))
    !
    RETURN
    !
  END FUNCTION gaussian
  !
  !
  !
  !
END MODULE wl_new
