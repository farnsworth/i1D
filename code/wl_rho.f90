!it is written to work with even number of fermions
!
! ... bugs: for now it is conserving the parity of the number of fermions
!           no more necessary
!
!
MODULE wl_rho
  !
  USE system
  !
  ! numerical parameters
  REAL (kind=dp) :: f = 2.718281828459_dp, threshold = 0.2_dp, accuracy = 1.0d-6,delta
  !
  REAL (kind=dp), DIMENSION(:,:),ALLOCATABLE :: logdosf
  INTEGER,DIMENSION(:,:), ALLOCATABLE :: domain
  !
  REAL (kind=dp) :: rhomax_out,rhomin_out
  !
CONTAINS
  !
  ! it works in the reduced space
  !
  FUNCTION wl_rhojdos( d, emin, emax, rhomin,rhomax, rhonbin, seed )
    !
    IMPLICIT NONE
    !
    !
    INTEGER,INTENT(in) :: rhonbin,d
    REAL(kind=8),INTENT(in) :: emax, emin
    REAL(kind=8),INTENT(in) :: rhomax,rhomin
    INTEGER, DIMENSION(8),INTENT(in) :: seed
    !
    INTEGER :: wl_rhojdos
    !
    LOGICAL, DIMENSION(:), ALLOCATABLE :: conf,new_conf
    INTEGER :: inew_conf
    !
    REAL(kind=8), DIMENSION(:,:), ALLOCATABLE :: logdos
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: histo
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: filled
    !
    !
    REAL(kind=8), DIMENSION(:), ALLOCATABLE :: egrid
    REAL(kind=8) :: econf, new_econf, deltae
    REAL(kind=8) :: rhoconf, new_rhoconf, deltarho
    REAL(kind=8) :: logf
    !
    INTEGER, DIMENSION(2) :: binconf, new_binconf
    INTEGER :: i,ierr
    INTEGER :: enbin, i_swe, i_swe_tot, step
    !
    !
    wl_rhojdos = 0
    !
    CALL initialize
    IF (ierr /= 0) THEN
       wl_rhojdos = ierr
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
          !
          i_swe = i_swe + 1
          ! now a sweep is composed by L/2 flips
          DO i=1,L/2
             CALL get_new_conf
             CALL mc_step
             CALL update
          END DO
          !
          IF ( mod(i_swe,100)==0 ) THEN
             !
             IF (check(histo,enbin,rhonbin,i_swe*L/2)) THEN
                EXIT
             END IF
             !
             IF (i_swe>10000000) THEN
                CALL error('wl_rhojdos','some problem in the convergence')
                CALL save
                wl_rhojdos = 2
                RETURN
             END IF
             !
          END IF
          !
       END DO
       !
       i_swe_tot = i_swe_tot + i_swe
       !
       IF (mod(step,5) == 0) THEN
          PRINT*,"Step:",step,"logf :",logf
          PRINT*,"Flat histo in",i_swe,"sweps"
          PRINT*
       END IF
       
       logf = logf/2.0_dp
       
    END DO
    
    PRINT*,"----------------- RUN IS ENDED -------------------"
    PRINT*
    print*,"Obtained convergence in ",i_swe_tot,"sweps"
    CALL save
    CALL finalize
    RETURN
    
  CONTAINS
    
    FUNCTION check(histo_in,enbin_in, rhonbin_in,npoint_in)
      
      IMPLICIT NONE
      
      LOGICAL :: check
      INTEGER, INTENT(IN) :: enbin_in, rhonbin_in , npoint_in
      INTEGER, DIMENSION(enbin_in,rhonbin_in), INTENT(IN) :: histo_in
      REAL(kind=8) :: meanhisto
      INTEGER :: i_tmp,j_tmp
      INTEGER :: nbin_tmp
      
      nbin_tmp = 0
      DO i_tmp=1,enbin_in
         DO j_tmp =1,rhonbin_in
            IF (histo_in(i_tmp, j_tmp) > 0) THEN
               nbin_tmp = nbin_tmp + 1
               IF (filled(i_tmp, j_tmp) == 0) THEN
                  filled(i_tmp,j_tmp) = 1
                  IF (step/=1) THEN
                     PRINT*,"ERROR : in the first step you did not scan all energies"
                     PRINT*, "--- restart from zero ---"
                     logf = dlog(f)
                     check = .true.
                     GO TO 11
                  END IF
               END IF
            ELSE
               IF (filled(i_tmp,j_tmp) == 1) THEN
                  check = .false.
                  GO TO 11
               END IF
            END IF
         END DO
      END DO
      
      meanhisto = dble(npoint_in)/dble(nbin_tmp)
      
      check = .true.
      DO i_tmp=1,enbin_in
         DO j_tmp=1,rhonbin_in
            IF (histo_in(i_tmp,j_tmp) /= 0) THEN
               IF (abs(meanhisto-dble(histo_in(i_tmp, j_tmp)))/meanhisto > threshold) THEN
                  !print*,"threshold",abs(meanhisto-dble(histo_in(i_tmp, j_tmp)))/meanhisto
                  check = .false.
                  EXIT
               END IF
            END IF
         END DO
      END DO
      
11    CONTINUE
      
    END FUNCTION check
    
    
    
    SUBROUTINE mc_step
      
      IMPLICIT NONE
      
      REAL(kind=8) :: rnd_tmp
      REAL(kind=8) :: rate_tmp
      
      
      IF ((new_binconf(1)>0).and.(new_binconf(1)<enbin+1)) THEN
         rate_tmp = logdos(binconf(1),binconf(2)) - logdos(new_binconf(1),new_binconf(2))
         IF ( rate_tmp < 0.0_dp ) THEN
            rate_tmp = dexp(rate_tmp)
            CALL random_number(rnd_tmp)
            IF (rnd_tmp > rate_tmp) GO TO 10
         END IF
         
         binconf = new_binconf
         econf = new_econf
         rhoconf = new_rhoconf
         conf(inew_conf) = new_conf(inew_conf)
         
      END IF
      
10    RETURN
      
    END SUBROUTINE mc_step
    !
    !
    SUBROUTINE update
      IMPLICIT NONE
      
      logdos( binconf(1), binconf(2) ) = logdos( binconf(1), binconf(2) ) + logf
      histo ( binconf(1), binconf(2) ) =  histo( binconf(1), binconf(2) ) + 1
      
    END SUBROUTINE update
    !
    !
    SUBROUTINE get_new_conf
      !
      IMPLICIT NONE
      !
      REAL(kind=8) :: rnd_tmp
      !
      new_econf = econf
      !
      !i do twice because i want to conser the parity of fermions
      CALL random_number( rnd_tmp )
      inew_conf = int(dble(L/2)*rnd_tmp) + 1
      !
      new_conf = conf
      new_conf( inew_conf ) = .not. conf( inew_conf )
      new_rhoconf = xxcorrelation_red(d, new_conf , L/2 )
      !
      IF ( conf( inew_conf ) ) THEN
         new_econf = new_econf - 2.0*egrid( inew_conf )
      ELSE
         new_econf = new_econf + 2.0*egrid( inew_conf )
      ENDIF
      
      new_binconf(1) = int( (new_econf - emin)/deltae ) + 1
      new_binconf(2) = int( (new_rhoconf - rhomin)/deltarho ) + 1
      
    END SUBROUTINE get_new_conf
    
    
!-----------------------------------------------------------------------!
    
    SUBROUTINE initialize
      IMPLICIT NONE
      
      REAL(kind=8) :: delta_tmp,rnd_tmp,emax_tmp
      REAL(kind=8) :: k_tmp,gapmin_tmp,gap_tmp
      INTEGER :: i_tmp,j_tmp,l1_tmp,l2_tmp, l0_tmp
      
      
      ierr = 0
      IF (mod(L,2) /= 0) THEN
         PRINT*,"Error: only even L are admitted"
         STOP
      END IF
      
      logf = dlog(f)
      
      IF ( emax < emin ) THEN
         CALL error('initialize','error emax is smaller than emin')
         ierr = 1
         RETURN
      ENDIF
      
      CALL random_seed(put=seed)
      
      ALLOCATE(conf(L/2))
      ALLOCATE(new_conf(L/2))
      ALLOCATE(egrid(L/2))
      !
      delta_tmp = 2.0_dp*pi/dble(L)
      !
      k_tmp = pi/dble(L)
      !
      egrid(1) = energy(k_tmp, h, gamma)
      econf = - egrid(1)
      !
      gapmin_tmp = egrid(1)
      emax_tmp = egrid(1)
      !
      DO i_tmp=2,L/2
         !
         k_tmp = k_tmp+delta_tmp
         !energy associated to point k_tmp
         egrid(i_tmp) = energy(k_tmp, h, gamma)
         !
         !minimum gap using single particle state
         !this method is good for gamma = 1 and h \= 0 where the energy is monotonic for k>0 or k<0
         !in this way at every mcstep the random wolker change bin
         gap_tmp = abs( egrid(i_tmp) - egrid(i_tmp-1) )
         IF ( (gap_tmp > 1e-12) .and. (gap_tmp < gapmin_tmp) ) THEN
            gapmin_tmp = gap_tmp
         END IF
         !
         !ground state energy
         econf = econf - egrid(i_tmp)
         !
         ! i compute the maximum energy
         IF ( emax_tmp < egrid(i_tmp) )  emax_tmp = egrid(i_tmp)
      END DO
      !
      !check if emax and emin are compatible with gs energy
      IF ((emax < econf/dble(L) ).or.(emin > -econf/dble(L) )) THEN
         CALL error("initialize", "energy interval outside energy spectra" )
         ierr = 1
         RETURN
      END IF
      !
      !
      IF ( ( econf/dble(L) < emax ) .and. ( econf/dble(L) > emin )) THEN
         GO TO 10
      END IF
      !
      !check if i can always find a good initial state
      IF ( (emax-emin)<(2.0*emax_tmp/dble(L)) ) THEN
         PRINT*,"Warning:"," i am not sure to get an initial state"
      END IF
      !
      !
      conf = .false.
      l0_tmp = 0
      !
      DO
         CALL random_number(rnd_tmp)
         j_tmp = int(rnd_tmp*(L/2-l0_tmp)) + 1
         
         l2_tmp = 0
         
         DO l1_tmp = 1, L/2
            
            IF ( .not.conf(l1_tmp) ) THEN
               l2_tmp = l2_tmp + 1
               IF (l2_tmp == j_tmp) THEN
                  conf(l1_tmp) = .true.
                  econf = econf + 2.0_dp*egrid(l1_tmp)
                  l0_tmp = l0_tmp + 1
                  IF ( ( econf/dble(L) < emax ) .and. ( econf/dble(L) > emin )) THEN
                        GO TO 10
                  END IF
                  GO TO 11
               END IF
            END IF
11          CONTINUE
            
         END DO
         
         IF (l0_tmp == L) EXIT
         
      END DO
      !
      CALL error("intialize","I was not able to generate the initial configuration")
      !
10    CONTINUE
      !
      PRINT*,"finded configuration with energy",econf/dble(L)
      !
      ! compute necessary number of bins
      gapmin_tmp = gapmin_tmp/dble(L)
      gapmin_tmp = 0.5_dp*gapmin_tmp
      !
      enbin = int((emax-emin)/gapmin_tmp)+1
      if (enbin < 2) CALL error("initialize","too small interval")
      !
      ! energy per unit of volume
      egrid = egrid/dble(L)
      econf = econf/dble(L)
      rhoconf = xxcorrelation_red(d, conf , L/2 )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! i put one more bin because otherwise there are problem at the edge bin
      deltae = (emax-emin)/dble(enbin)
      !emax = emax + deltae*0.5_dp
      !emin = emin - deltae*0.5_dp
      !enbin = enbin + 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      deltarho = (rhomax-rhomin)/dble(rhonbin)
      
      
      ALLOCATE(histo(enbin,rhonbin))
      ALLOCATE(filled(enbin,rhonbin))
      ALLOCATE(logdos(enbin,rhonbin))
      
      histo = 0
      filled = 0
      logdos = 0.0_dp
      
      
      binconf(1) = int((econf-emin) / deltae )+1
      IF ((binconf(1)<1).or.(binconf(1)>enbin)) THEN
         ierr = 1
         CALL error("initialize", "energy interval outside energy spectra" )
         RETURN
      END IF
      
      binconf(2) = int((rhoconf-rhomin) / deltarho )+1
      IF ((binconf(2)<1).or.(binconf(2)>rhonbin)) THEN
         CALL error("initialize", "magnetization outside its interval" )
      END IF
      
      PRINT*,enbin,"X",rhonbin," bins"
      PRINT*
      PRINT*,"----------------- RUN IS STARTED -----------------"
      
    END SUBROUTINE initialize
    !
    SUBROUTINE finalize
      IMPLICIT NONE
      !    CLOSE(unit=2)
      DEALLOCATE(conf,egrid,logdos,histo,filled,new_conf)      
    END SUBROUTINE finalize
    !
    SUBROUTINE save
      
      IMPLICIT NONE
      
      INTEGER :: i_tmp,neffbin,ieff,j_tmp
      REAL(kind=8) :: alpha_tmp
      
      neffbin = 0
      alpha_tmp = -1.0
      DO i_tmp=1,enbin
         DO j_tmp=1,rhonbin
            IF (filled(i_tmp,j_tmp) == 1) THEN
               neffbin = neffbin + 1
               IF (alpha_tmp < 0) alpha_tmp = logdos(i_tmp,j_tmp)
            END IF
         END DO
      END DO
      
      IF (allocated(logdosf)) deallocate(logdosf)
      IF (allocated(domain)) deallocate(domain)
      
      ALLOCATE( logdosf( neffbin,3 ) )
      ALLOCATE( domain(enbin,rhonbin) )
      domain = 0
      ieff = 0
      DO i_tmp=1,enbin
         DO j_tmp=1,rhonbin
            
            IF (filled(i_tmp,j_tmp) == 1) THEN
               domain(i_tmp,j_tmp) = 1
               ieff = ieff + 1
               logdosf( ieff, 1) = emin+(dble(i_tmp)-0.5)*deltae
               logdosf( ieff, 2) = rhomin+(dble(j_tmp)-0.5)*deltarho
               logdosf( ieff, 3)  = logdos(i_tmp,j_tmp)-alpha_tmp
            END IF
         END DO
      END DO
      
    END SUBROUTINE save
    
  END FUNCTION wl_rhojdos
  !
  !
  !
  !
  !
  FUNCTION wl_rhodos( d, ener, deltae, rhomin, rhomax, rhonbin, seed, readconf )
    !
    IMPLICIT NONE
    !
    !
    REAL (kind=8), INTENT(IN) :: ener,deltae
    INTEGER, INTENT(IN) :: d
    INTEGER, DIMENSION(8),INTENT(IN) :: seed
    LOGICAL, OPTIONAL,INTENT(IN) :: readconf
    INTEGER,INTENT(IN) :: rhonbin
    REAL (kind=8),INTENT(IN):: rhomax, rhomin
    INTEGER :: wl_rhodos
    !
    !
    LOGICAL, DIMENSION(:), ALLOCATABLE :: conf
    LOGICAL, DIMENSION(:), ALLOCATABLE :: new_conf
    !
    REAL(kind=8), DIMENSION(:), ALLOCATABLE :: logdos
    INTEGER, DIMENSION(:), ALLOCATABLE :: histo
    INTEGER, DIMENSION(:), ALLOCATABLE :: filled
    !
    REAL (kind = 8), DIMENSION(:), ALLOCATABLE :: egrid
    REAL (kind = 8) :: econf, new_econf
    REAL (kind = 8) :: rhoconf, new_rhoconf

    REAL (kind = 8) :: logf,deltarho
    !
    INTEGER :: binconf, new_binconf
    INTEGER :: i,inew_conf,ierr
    INTEGER :: i_swe, i_swe_tot, step
    !
    !
    wl_rhodos = 0
    !
    CALL initialize
    IF (ierr /= 0) THEN
       wl_rhodos = ierr
       RETURN
    ENDIF
    !
  !!    OPEN(unit=2,file=trim(datadir)//'intconf.out',action='write',status='replace')
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
          DO i=1,L/2
             !print*,"state",conf
             !print*,"enconf",econf
             !print*
             CALL get_new_conf
             CALL mc_step
             CALL update
             !CALL save_intermediate_conf(conf,L)
          END DO
          !
          IF ( mod(i_swe,1000)==0 ) THEN
             !
             CALL check_obs( conf, econf, 2.0d0*egrid , L/2 )
             !
             IF (check(histo,rhonbin,i_swe*L/2)) THEN
                EXIT
             END IF
             !
             IF (i_swe>1000000) THEN
                CALL error('wl_rhodos','Some problem in the convergence')
                wl_rhodos = 2
                RETURN
             END IF
             !
          END IF
          !
       END DO
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
    !  
    PRINT*,"----------------- RUN IS ENDED -------------------"
    PRINT*
    print*,"Obtained convergence in ",i_swe_tot,"sweps"
    CALL save
    CALL finalize
    
! !   close(unit=2)
    !
    RETURN
    !  
  CONTAINS
    !
    FUNCTION check(histo_in, rhonbin_in,npoint_in)
      !
      IMPLICIT NONE
      !
      LOGICAL :: check
      INTEGER, INTENT(IN) :: rhonbin_in , npoint_in
      INTEGER, DIMENSION(rhonbin_in), INTENT(IN) :: histo_in
      REAL(kind=8) :: meanhisto
      INTEGER :: i_tmp
      INTEGER :: nbin_tmp
      !  
      nbin_tmp = 0
      DO i_tmp=1,rhonbin_in
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
      !
      meanhisto = dble(npoint_in)/dble(nbin_tmp)
      !   
      check = .true.
      DO i_tmp=1,rhonbin_in
         IF (histo_in(i_tmp) /= 0) THEN
            IF (abs(meanhisto-dble(histo_in(i_tmp)))/meanhisto > threshold) THEN
               check = .false.
               EXIT
            END IF
         END IF
      END DO
      !
11    CONTINUE
      !
    END FUNCTION check
    !
    SUBROUTINE mc_step
      !   
      IMPLICIT NONE
      !
      REAL(kind=8) :: rnd_tmp
      REAL(kind=8) :: rate_tmp
      !
      IF ( dabs(new_econf-ener) <= deltae ) THEN
         ! update max and min value of rho_xx
         IF ( new_rhoconf > rhomax_out ) rhomax_out = new_rhoconf
         IF ( new_rhoconf < rhomin_out ) rhomin_out = new_rhoconf
         !
         rate_tmp = logdos(binconf) - logdos(new_binconf)
         IF ( rate_tmp < 0.0_dp ) THEN
            rate_tmp = dexp(rate_tmp)
            CALL random_number(rnd_tmp)
            IF (rnd_tmp > rate_tmp) GO TO 10
         END IF
         
         binconf = new_binconf
         econf = new_econf
         rhoconf = new_rhoconf         
         conf( inew_conf ) = new_conf( inew_conf )
         
      END IF

 !     WRITE(unit=2,fmt=*) (dabs(new_econf-ener) <= deltae)
      
10    RETURN
      
    END SUBROUTINE mc_step
    !
    SUBROUTINE update
      !
      IMPLICIT NONE
      !
      logdos( binconf ) = logdos( binconf ) + logf
      histo ( binconf ) =  histo( binconf ) + 1
      !
    END SUBROUTINE update
    !
    SUBROUTINE get_new_conf
      !
      IMPLICIT NONE
      !    
      REAL(kind=8) :: rnd_tmp
      !
      new_econf = econf
      CALL random_number( rnd_tmp )
      inew_conf = int(dble(L/2)*rnd_tmp) + 1
      !print*,"trial",inew_conf
      IF (inew_conf>L/2) inew_conf = L/2
      !
      new_conf = conf
      new_conf( inew_conf ) = .not. conf( inew_conf )
      new_rhoconf = xxcorrelation_red(d, new_conf , L/2 )
      !
      IF ( conf( inew_conf ) ) THEN
         new_econf = new_econf - 2.0_dp*egrid( inew_conf )
      ELSE
         new_econf = new_econf + 2.0_dp*egrid( inew_conf )
      ENDIF
      !
      new_binconf = int( (new_rhoconf - rhomin)/deltarho ) + 1
      if (new_binconf > rhonbin) new_binconf = rhonbin 
      ! 
    END SUBROUTINE get_new_conf
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
      ALLOCATE(conf(L/2))
      ALLOCATE(new_conf(L/2))
      ALLOCATE(egrid(L/2))
      !
      delta_tmp = 2.0_dp*pi/dble(L)
      !
      k_tmp = -pi/dble(L)
      !
      !
      DO i_tmp=1,L/2
         !
         k_tmp = k_tmp+delta_tmp
         !energy associated to point k_tmp
         egrid(i_tmp) = energy(k_tmp, h, gamma )
         !
         !ground state energy
         econf = econf - egrid(i_tmp)
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
      !
      deltarho = (rhomax-rhomin)/dble(rhonbin)
      !
      !
      !check if emax and emin are compatible with gs energy
      IF ((ener+deltae < econf/dble(L) ).or.(ener-deltae > -econf/dble(L) )) THEN
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
               econf = econf + 2.0_dp*egrid(l1_tmp)
            ENDIF
         END DO
         IF (  dabs( econf/dble(L)-ener) > deltae ) THEN
            ierr = 1
            RETURN
         END IF
         !
         GOTO 10
         !
      END IF
      !
      !
      !check if i can always find a good initial state
      IF ( (2.0d0*deltae)<(2.0d0*emax_tmp/dble(L)) ) THEN
         PRINT*,"Warning:"," i am not sure to get an initial state"
         PRINT*,"Window :",2.0d0*deltae
         PRINT*,"Maximum gap:",(2.0*emax_tmp/dble(L))
      END IF
      !
      !
      l0_tmp = 0
      !
      !
      DO
         CALL random_number(rnd_tmp)
         j_tmp = int(rnd_tmp*(L/2-l0_tmp)) + 1
         IF (j_tmp>L/2) j_tmp = L/2
         !
         l2_tmp = 0
         !
         DO l1_tmp = 1, L/2
            !
            IF ( .not.conf(l1_tmp) ) THEN
               l2_tmp = l2_tmp + 1
               IF (l2_tmp == j_tmp) THEN
                  conf(l1_tmp) = .true.
                  econf = econf + 2.0_dp*egrid(l1_tmp)
                  l0_tmp = l0_tmp + 1
                  IF (mod(l0_tmp,2) == 0) THEN
                     print*,"econf :",econf/dble(L)
                     IF ( dabs( econf/dble(L) - ener ) <= deltae ) THEN
                        CALL save_conf(conf,L/2)
                        GO TO 10
                     ENDIF
                  END IF
                  GO TO 11
               END IF
            END IF
            !
11          CONTINUE
            !
         END DO
         !
         IF (l0_tmp == L/2) EXIT
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
      rhoconf = xxcorrelation_red(d, conf , L/2 )
      !
      rhomax_out = rhoconf
      rhomin_out = rhoconf
      !
      !
      !
      !
      ALLOCATE(histo(rhonbin))
      ALLOCATE(filled(rhonbin))
      ALLOCATE(logdos(rhonbin))
      !
      histo = 0
      filled = 0
      logdos = 0.0_dp
      !
      binconf = int((rhoconf-rhomin) / deltarho )+1
      IF (binconf > rhonbin) binconf = rhonbin 
      !
      !
      IF ((binconf<1).or.(binconf>rhonbin)) THEN
         CALL error("initialize", "magnetization outside its interval" )
         ierr = 1
         RETURN
      END IF
      !
      PRINT*
      PRINT*,"----------------- RUN IS STARTED -----------------"
      !
      CALL check_obs( conf, econf, 2.0_dp*egrid , L/2 )
      !
      RETURN
      !
      !
    END SUBROUTINE initialize
    !
    !
    SUBROUTINE finalize
      IMPLICIT NONE
      !
      !CLOSE(unit=2)
      !delta = deltarho
      DEALLOCATE(conf,egrid,logdos,histo,filled,new_conf)
      !
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
      !
      neffbin = 0
      alpha_tmp = -1.0
      DO i_tmp=1,rhonbin
         IF (filled(i_tmp) == 1) THEN
            neffbin = neffbin + 1
            IF (alpha_tmp < 0) alpha_tmp = logdos( i_tmp )
         END IF
      END DO
      !
      !
      IF (allocated(logdosf)) deallocate(logdosf)
      allocate( logdosf( neffbin, 2 ) )
      !
      ieff = 0
      DO i_tmp=1,rhonbin
         IF (filled(i_tmp) == 1) THEN
            ieff = ieff + 1
            logdosf( ieff, 1) = rhomin+(dble(i_tmp)-0.5)*deltarho
            logdosf( ieff, 2)  = logdos(i_tmp)-alpha_tmp
         END IF
      END DO
      !
    END SUBROUTINE save
    !
  END FUNCTION wl_rhodos
  !
  !
  ! ... it does a check of the observable
  SUBROUTINE check_obs( conf, obsConf, obs , dim )
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: dim
    LOGICAL, DIMENSION(dim), INTENT(IN) :: conf
    REAL(kind=8), INTENT(INOUT) :: obsConf
    REAL(kind=8),DIMENSION(dim), INTENT(IN) :: obs
    INTEGER :: i
    REAL(kind=8) :: obsConf_check
    !
    obsConf_check = 0.0d0
    DO i=1,dim
       !
       IF (conf(i)) THEN
          obsConf_check = obsConf_check + 0.5d0*obs(i)
       ELSE
          obsConf_check = obsConf_check - 0.5d0*obs(i)
       ENDIF
       !
    ENDDO
    !
    IF ( abs(obsConf - obsConf_check)/abs(obsConf_check) > 1.0e-5 ) THEN
       IF (abs(obsConf)>1e-10) PRINT '( a, e10.3, a )',"Warning: correction of &
            &observable ( relative error",abs(obsConf - obsConf_check)/abs(obsConf_check)*100.0,"%)"
       obsConf = obsConf_check
       
    END IF

  END SUBROUTINE check_obs
  !
  !
  !
  SUBROUTINE save_conf(conf, dim)
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: dim
    LOGICAL, DIMENSION(dim), INTENT(IN) :: conf
    !
    OPEN(unit=1,file=trim(datadir)//'initialconf.out',action='write',status='replace')
    !
    WRITE(unit=1,fmt=*) conf
    !
    CLOSE(unit=1)
    !
    RETURN
    !
  END SUBROUTINE save_conf
  !
  !
  !
  SUBROUTINE save_intermediate_conf(conf, dim)
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: dim
    LOGICAL, DIMENSION(dim), INTENT(IN) :: conf
    !
    !
    WRITE(unit=2,fmt=*) conf
    !
    !
    RETURN
    !
  END SUBROUTINE save_intermediate_conf
  !
  !
  !
  SUBROUTINE read_conf(conf, dim)
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: dim
    LOGICAL, DIMENSION(dim), INTENT(OUT) :: conf
    !
    OPEN(unit=1,file=trim(datadir)//'initialconf.out',action='read',status='old')
    !
    READ(unit=1,fmt=*) conf
    !
    CLOSE(unit=1)
    !
    RETURN
    !
  END SUBROUTINE read_conf
  !
  !
  !
END MODULE wl_rho
