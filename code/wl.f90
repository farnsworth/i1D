!it is written to work with even number of fermions -> add for odd 
!    number of fermions
!
! ... bugs: for now it is conserving the parity of the number of fermions
!           no more necessary
!
!
! Functions:
!    wl_edos : calculates the denisty of states in a fixed 
!              energy interval
!    wl_jdos : calculates the joined density of state using an
!              interval of energy
!    wl_mdos : does a direct calculation of mdos in a fixed energy
!              interval
!    check_obs : does a check of the calculated observable
!
!
MODULE wl
  !
  USE system
  !
  ! numerical parameters
  REAL (kind=dp) :: f = 2.718281828459_dp, threshold = 0.2_dp, accuracy = 1.0d-6,delta
  !
  REAL (kind=dp), DIMENSION(:,:),ALLOCATABLE :: logdosf
  INTEGER :: ierr
  !
  PUBLIC :: read_conf, save_intermediate_conf, save_conf, check_obs
  !
  !
CONTAINS
  !
  ! it calculates the denisty of states in the interval emin emax
  !
  FUNCTION wl_edos( emin, emax , nbin , seed , readconf )
    !
    IMPLICIT NONE
    !
    REAL(kind=8),INTENT(IN) :: emax, emin
    INTEGER, INTENT(IN) :: nbin
    INTEGER, DIMENSION(8),INTENT(IN),OPTIONAL :: seed
    LOGICAL, OPTIONAL,INTENT(IN) :: readconf
    INTEGER :: wl_edos
    !
    !
    LOGICAL, DIMENSION(:), ALLOCATABLE :: conf
    INTEGER :: new_conf
    !
    REAL(kind=8), DIMENSION(:), ALLOCATABLE :: engrid,logdos
    INTEGER, DIMENSION(:), ALLOCATABLE :: histo
    INTEGER, DIMENSION(:), ALLOCATABLE :: filled
    !
    REAL(kind=8) :: enconf, new_enconf, logf, deltaen,norm
    INTEGER :: i,binconf, new_binconf, index_min
    INTEGER :: i_swe,i_swe_tot,step,npoints,npoints_inside
    !
    CHARACTER(len=10) :: string1
    CHARACTER(len=10) :: string2
    !
    !
    CALL initialize
    IF (ierr /= 0) THEN
       wl_edos = ierr
       RETURN
    ENDIF
    !
    !
    i_swe_tot = 0
    step = 0
    !
    !
    DO WHILE ( logf > accuracy )
       !
       step = step + 1
       i_swe = 0
       !
       norm = 0.0d0
       DO i=1,nbin
          !
          IF ( histo(i) > 0 ) THEN
             norm = logdos(i)
             EXIT
          ENDIF
          !
       END DO
       !
       DO i=1,nbin
          logdos(i) = logdos(i) - norm
       ENDDO
       !
       !       IF ( ( mod( i_swe,10 )==0 ) .and. ( mod(step,1)==0 ) ) THEN
          WRITE(string1,'(i2)') step
          WRITE(string2,'(i6)') i_swe
          string1 = adjustl(string1)
          string2 = adjustl(string2)
          OPEN(unit=1,file=trim(string1)//"-"//trim(string2)//".out",action="write",status="replace")
          CALL snapshot()
          CLOSE(unit=1)
!       ENDIF
       !
       histo = 0
       npoints = 0
       npoints_inside = 0
       !
       !
       !
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
          IF ( ( mod( i_swe,10 )==0 ) .and. ( mod(step,1)==0 ) ) THEN
             WRITE(string1,'(i2)') step
             WRITE(string2,'(i6)') i_swe
             string1 = adjustl(string1)
             string2 = adjustl(string2)
             OPEN(unit=1,file=trim(string1)//"-"//trim(string2)//".out",action="write",status="replace")
             CALL snapshot()
             CLOSE(unit=1)
          ENDIF
          !
          IF ( mod(i_swe,10)==0 ) THEN
             !
             IF (check(histo,nbin,npoints)) THEN
                EXIT
             END IF
             !
             IF (i_swe>100000000) THEN
                PRINT*,"Some problem in the convergence"
                STOP 1
             END IF
             !
          END IF
          !
       END DO
       !
       !
       i_swe_tot = i_swe_tot + i_swe
       !
       IF (mod(step,5) == 0) THEN
          PRINT*,"Step:",step,"logf :",logf
          PRINT*,"Flat histo in",i_swe,"sweps"
          print*,"ratio inside = ", dble(npoints_inside)/dble(npoints)
          PRINT*
       END IF
       
       logf = logf/2.0_dp
       
    END DO
    !
    PRINT*,"----------------- RUN IS ENDED -------------------"
    PRINT*
    PRINT*,"Obtained convergence in ",i_swe_tot,"sweps"
    CALL save
    CALL finalize
    !
    wl_edos = ierr
    !
  CONTAINS
    !
    !
    SUBROUTINE snapshot()
      !
      INTEGER :: i_tmp
      LOGICAL :: first
      !
      first = .true.
      DO i_tmp=1,nbin
         !
         IF ( histo(i_tmp) > 0 ) THEN
         !
!            IF (first) THEN
!               first = .false.
!               norm = logdos(i_tmp)
!            ENDIF
            WRITE(unit=1,fmt=*) emin+(dble(i_tmp)-0.5)*deltaen, logdos(i_tmp), histo(i_tmp)
                  ENDIF
         !
      END DO
      !
      WRITE(unit=1,fmt=*)
      WRITE(unit=1,fmt=*)
      !
      RETURN
      !
    END SUBROUTINE snapshot
    !
    !
    SUBROUTINE initialize
      
      REAL(kind=8) :: delta_tmp,rnd_tmp
      REAL(kind=8) :: k_tmp,data_tmp,sp_emax,sp_emin
      INTEGER :: i_tmp,j_tmp,l1_tmp,l2_tmp, l0_tmp
      LOGICAL :: readconf_tmp
      
      ierr = 0

      IF (mod(L,2) /= 0) THEN
         PRINT*,"Error: only even L are admitted"
         STOP
      END IF
      
      logf = dlog(f)
      
      print*

      IF ( emax < emin ) THEN
         CALL error('initialize','error emax is smaller than emin')
         ierr = 1
         RETURN
      ENDIF
      
      IF (present(seed)) THEN
         CALL random_seed(put=seed)
      ELSE
         CALL random_seed(put=(/ 123,4345,25423,522,524,234,9653,9473 /) )
      END IF
      !
      ALLOCATE(conf(L))
      ALLOCATE(engrid(L))
      !
      delta_tmp = 2.0_dp*pi/dble(L)
      !
      conf = .false.
      !
      k_tmp = pi*(-1.0_dp+1.0_dp/dble(L))
      !
      engrid(1) = energy(k_tmp, h, gamma)
      enconf = -engrid(1)/2.0_dp
      sp_emax = engrid(1)
      sp_emin = engrid(1)
      !
      !
      DO i_tmp=2,L
         !
         k_tmp = k_tmp+delta_tmp
         !energy associated to point k_tmp
         engrid(i_tmp) = energy(k_tmp, h, gamma)
         !
         ! minimum gap using single particle state
         ! this method is good for gamma = 1 and h \= 0 where the energy 
         ! is monotonic for k>0 or k<0
         ! in this way at every mcstep the random wolker change bin
         !gap_tmp = abs( engrid(i_tmp) - engrid(i_tmp-1) )
         !IF ( (gap_tmp > 1e-12) .and. (gap_tmp < gapmin_tmp) ) THEN
         !   gapmin_tmp = abs( engrid(i_tmp) - engrid(i_tmp-1) )
         !END IF
         !
         ! ground state energy
         enconf = enconf - engrid(i_tmp)/2.0_dp
         !
         ! I compute the maximum energy
         IF ( sp_emax < engrid(i_tmp) ) sp_emax = engrid(i_tmp)
         IF ( sp_emin > engrid(i_tmp) ) sp_emin = engrid(i_tmp)
         !
      END DO
      !
      !check if emax and emin are compatible with gs energy
      IF ((emax < enconf/dble(L) ).or.(emin > -enconf/dble(L) )) THEN
         PRINT*,"gs energy",enconf/dble(L)
         CALL error("initialize", "energy interval outside energy spectra" )
         ierr = 7
         RETURN
      END IF
      !
      PRINT*,"Min single part en:",sp_emin/dble(L)
      PRINT*,"Max single part en:",sp_emax/dble(L)
      !
      data_tmp = engrid(1)
      index_min = 1
      !
      DO i=2,L
         IF (engrid(i)< data_tmp) THEN
            data_tmp = engrid(i)
            index_min = i
         ENDIF
      ENDDO
      !
      print*,"minimum energy",data_tmp/dble(L)
      !
      !
      IF ( ( enconf/dble(L) < emax ) .and. ( enconf/dble(L) > emin )) THEN
         GO TO 10
      END IF
      !
      readconf_tmp = .false.
      IF (present(readconf)) THEN
         readconf_tmp = readconf
      ENDIF
      !
      IF (readconf_tmp) THEN
         !
         CALL read_conf(conf,L)
         DO l1_tmp = 1, L
            IF ( conf(l1_tmp) ) THEN
               enconf = enconf + engrid(l1_tmp)
            ENDIF
         END DO
         !
         IF ( ( enconf/dble(L) > emax ) .or. ( enconf/dble(L) < emin )) THEN
            ierr = 2
            RETURN
         END IF
         !
         GOTO 10
         !
      END IF
      !
      !check if i can always find a good initial state
      IF ( (emax-emin) < (sp_emax/dble(L)) ) THEN
         PRINT*,"Warning:"," i am not sure to get an initial state"
      END IF
      !
      !
      !
      !
      l0_tmp = 0
      
      DO
         CALL random_number(rnd_tmp)
         j_tmp = int(rnd_tmp*(L-l0_tmp)) + 1
         IF (j_tmp>(L-l0_tmp)) j_tmp = L - l0_tmp
         
         l2_tmp = 0
         
         DO l1_tmp = 1, L
            
            IF ( .not.conf(l1_tmp) ) THEN
               l2_tmp = l2_tmp + 1
               IF (l2_tmp == j_tmp) THEN
                  conf(l1_tmp) = .true.
                  enconf = enconf + engrid(l1_tmp)
                  l0_tmp = l0_tmp + 1
                  !print*,"enconf :",enconf/dble(L)
                  IF ( ( enconf/dble(L) < emax ) .and. ( enconf/dble(L) > emin )) THEN
                     CALL save_conf(conf,L)
                     GO TO 10
                  END IF
                  GO TO 11
               END IF
            END IF
11          CONTINUE
            
         END DO
         
         IF (l0_tmp == L) EXIT
         
      END DO
      
      CALL error("intialize","I was not able to generate the initial configuration")
      ierr = 6
      RETURN

10    CONTINUE
      
      PRINT*,"finded configuration with energy",enconf/dble(L)
      
      ! compute necessary number of bins
      !gapmin_tmp = gapmin_tmp/dble(L)
      !gapmin_tmp = 0.5_dp*gapmin_tmp
      !
      !nbin = int((emax-emin)/gapmin_tmp)+1
      !if (nbin < 2) CALL error("initialize","too small interval")
      
      ! energy per unit of volume
      engrid = engrid/dble(L)
      enconf = enconf/dble(L)
      
      !if the interval is symmetric respect to 0 i want to put 0 in the center
      !of a bin, therefore i need odd number of bin
      !pay attention : after i add one more bin
      !IF (((emax+emin)<1e-8).and.(mod(nbin,2) == 0)) THEN
      !   nbin = nbin +1
      !   PRINT*,"unlacky case"
      !END IF
      !
      deltaen = (emax-emin)/dble(nbin)
      !
      ALLOCATE(histo(nbin))
      ALLOCATE(filled(nbin))
      ALLOCATE(logdos(nbin))
      
      histo = 0
      filled = 0
      logdos = 0.0_dp
      
      
      binconf = int((enconf-emin)/deltaen)+1
      IF (binconf == nbin + 1) binconf = nbin
      !
      IF ((binconf<0).or.(binconf>nbin)) THEN
         CALL error("initialize", "energy interval outside energy spectra" )
         ierr = 5
         RETURN
      END IF
      !
      PRINT*,nbin,"bins"
      PRINT*
      PRINT*,"----------------- RUN IS STARTED -----------------"
      
      
    END SUBROUTINE initialize
    
    
    
    SUBROUTINE finalize
      
      !    CLOSE(unit=2)
      delta = deltaen
      DEALLOCATE(conf,engrid,logdos,histo,filled)
      
    END SUBROUTINE finalize
    
    
    
    FUNCTION check(histo_in,nbin_in,npoint_in)
      
      LOGICAL :: check
      INTEGER, INTENT(IN) :: nbin_in,npoint_in
      INTEGER, DIMENSION(nbin_in), INTENT(IN) :: histo_in
      REAL(kind=8) :: meanhisto
      INTEGER :: i_tmp
      INTEGER :: nbin_tmp
      
      nbin_tmp = 0
      DO i_tmp=1,nbin_in
         IF (histo_in(i_tmp) > 0) THEN
            nbin_tmp = nbin_tmp + 1
            IF (filled(i_tmp) == 0) THEN
               filled(i_tmp) = 1
               IF (step/=1) THEN
                  PRINT*,"ERROR : in the first step you did not scan all energies"
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
      
      check = .true.
      DO i_tmp=1,nbin_in
         IF (histo_in(i_tmp) /= 0) THEN
            IF (abs(meanhisto-dble(histo_in(i_tmp)))/meanhisto > threshold) THEN
               check = .false.
!               print*,meanhisto,histo_in(i_tmp),(meanhisto-histo_in(i_tmp))/meanhisto
               EXIT
            END IF
         END IF
      END DO
      
11    CONTINUE
      
    END FUNCTION check
    
    
    
    SUBROUTINE save
      
      IMPLICIT NONE
      
      INTEGER :: i_tmp,neffbin,ieff
      REAL(kind=8) :: alpha_tmp
      
      neffbin = 0
      alpha_tmp = -1.0
      DO i_tmp=1,nbin
         IF (filled(i_tmp) == 1) THEN
            neffbin = neffbin + 1
            IF (alpha_tmp < 0) alpha_tmp = logdos(i_tmp)
!            IF (i_tmp == 20) alpha_tmp = logdos(i_tmp)
         END IF
      END DO
      
      IF (allocated(logdosf)) deallocate(logdosf)
      
      ALLOCATE( logdosf(neffbin, 2) )
      !
      ieff = 0
      DO i_tmp=1,nbin
         IF (filled(i_tmp) == 1) THEN
            ieff = ieff + 1
            logdosf(ieff,1) = emin+(dble(i_tmp)-0.5)*deltaen
            logdosf(ieff,2) = logdos(i_tmp)-alpha_tmp
         END IF
      END DO
      !
    END SUBROUTINE save
    !
    !
    !
    SUBROUTINE mc_step
      !
      REAL(kind=8) :: rnd_tmp
      REAL(kind=8) :: rate_tmp
      !
      !
! it doesn't work because if  emin-delta < new_enconf < emin
! we get new_binconf = 0!!!! due to how int() works
!      IF ((new_binconf>0).and.(new_binconf<nbin+1)) THEN
      IF ((new_enconf >= emin ).and.(new_enconf <= emax )) THEN
         !
         npoints_inside = npoints_inside + 1
         rate_tmp = logdos(binconf) - logdos(new_binconf)
         IF ( rate_tmp < 0.0_dp ) THEN
            rate_tmp = dexp(rate_tmp)
            CALL random_number(rnd_tmp)
            IF (rnd_tmp > rate_tmp) GO TO 10
         END IF
         !
         binconf = new_binconf
         enconf = new_enconf
         !
         conf(new_conf ) = .not.conf( new_conf )
         !
      END IF
      !
10    RETURN
      !
    END SUBROUTINE mc_step
    !
    !
    SUBROUTINE update
      !
      logdos(binconf) = logdos(binconf)+logf
      histo(binconf) = histo(binconf)+1
      npoints = npoints+1
      !
    END SUBROUTINE update
    !
    !
    SUBROUTINE get_new_conf
      !
      REAL(kind=8) :: rnd_tmp
      !
      new_enconf = enconf
      !
      CALL random_number(rnd_tmp)
      new_conf = int(dble(L)*rnd_tmp) + 1
      IF (new_conf == L+1 ) new_conf = L
      !
      ! 
      !
      !
      IF ( conf( new_conf ) ) THEN
         new_enconf = new_enconf - engrid( new_conf )
      ELSE
         new_enconf = new_enconf + engrid( new_conf )
      ENDIF
      !
!      IF (new_conf == index_min) THEN
!         PRINT*,"min trans:",enconf, new_enconf
!      ENDIF
!      !
!      IF ( (new_enconf>=emin) .and. (new_enconf<=emax) ) THEN
!         PRINT*,"c'e' almeno uno stato"
!         STOP
!      ENDIF
      !
      new_binconf = int( (new_enconf -emin)/deltaen )+1
      IF (new_enconf == emax) new_binconf = nbin
      !
    END SUBROUTINE get_new_conf
    
    
  END FUNCTION wl_edos
  !
  ! It calculate the density of state in the interval emin e max
  !
  SUBROUTINE wl_jdos( emin, emax , mnbin, seed )
    !
    IMPLICIT NONE
    !
    !
    INTEGER :: mnbin
    REAL(kind=8) :: emax, emin
    INTEGER, DIMENSION(8) :: seed
    !
    !f2py INTEGER, INTENT(IN) :: mnbin
    !f2py INTEGER, DIMENSION(8), INTENT(IN) :: seed
    !f2py REAL (kind=8), INTENT(IN) :: emax, emin
    !
    LOGICAL, DIMENSION(:), ALLOCATABLE :: conf
    INTEGER, DIMENSION(2) :: new_conf
    !
    REAL(kind=8), DIMENSION(:,:), ALLOCATABLE :: logdos
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: histo
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: filled
    !
    !
    REAL(kind=8), DIMENSION(:), ALLOCATABLE :: egrid, mgrid
    REAL(kind=8) :: econf, new_econf, deltae
    REAL(kind=8) :: mconf, new_mconf, deltam,mmin
    REAL(kind=8) :: logf
    !
    INTEGER, DIMENSION(2) :: binconf, new_binconf
    INTEGER :: i
    INTEGER :: enbin, i_swe, i_swe_tot, step
    !
    !
    CALL initialize
    
    
    i_swe_tot = 0
    step = 0
    
    DO WHILE ( logf > accuracy )
       
       step = step + 1
       i_swe = 0
       histo = 0
       
       DO
          i_swe = i_swe + 1
          
          DO i=1,L
             CALL get_new_conf
             CALL mc_step
             CALL update
          END DO
          
          IF ( mod(i_swe,100)==0 ) THEN
             
             IF (check(histo,enbin,mnbin,i_swe*L)) THEN
                EXIT
             END IF
             
             IF (i_swe>100000000) THEN
                PRINT*,"Some problem in the convergence"
                STOP 1
             END IF
             
          END IF
          
       END DO
       
       i_swe_tot = i_swe_tot + i_swe
       
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
    
  CONTAINS
    
    FUNCTION check(histo_in,enbin_in, mnbin_in,npoint_in)
      
      IMPLICIT NONE
      
      LOGICAL :: check
      INTEGER, INTENT(IN) :: enbin_in, mnbin_in , npoint_in
      INTEGER, DIMENSION(enbin_in,mnbin_in), INTENT(IN) :: histo_in
      REAL(kind=8) :: meanhisto
      INTEGER :: i_tmp,j_tmp
      INTEGER :: nbin_tmp
      
      nbin_tmp = 0
      DO i_tmp=1,enbin_in
         DO j_tmp =1,mnbin_in
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
         DO j_tmp=1,mnbin_in
            IF (histo_in(i_tmp,j_tmp) /= 0) THEN
               IF (abs(meanhisto-dble(histo_in(i_tmp, j_tmp)))/meanhisto > threshold) THEN
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
         mconf = new_mconf
         
         conf(new_conf(1)) = .not.conf(new_conf(1))
         conf(new_conf(2)) = .not.conf(new_conf(2))
         
      END IF
      
10    RETURN
      
    END SUBROUTINE mc_step
    
    
    
    SUBROUTINE update
      IMPLICIT NONE
      
      logdos( binconf(1), binconf(2) ) = logdos( binconf(1), binconf(2) ) + logf
      histo ( binconf(1), binconf(2) ) =  histo( binconf(1), binconf(2) ) + 1
      
    END SUBROUTINE update
    
    
    
    SUBROUTINE get_new_conf
      IMPLICIT NONE
      
      REAL(kind=8), DIMENSION(2) :: rnd_tmp
      INTEGER :: j_tmp
      
      new_econf = econf
      new_mconf = mconf
      
      
      
      !i do twice because i want to conser the parity of fermions
      CALL random_number( rnd_tmp )
      new_conf = int(dble(L)*rnd_tmp) + 1
      
      IF (new_conf(1) /= new_conf(2)) THEN
         DO j_tmp =1,2
            IF (conf( new_conf(j_tmp) )) THEN
               new_econf = new_econf - egrid( new_conf(j_tmp) )
               new_mconf = new_mconf - mgrid( new_conf(j_tmp) )
            ELSE
               new_econf = new_econf + egrid( new_conf(j_tmp) )
               new_mconf = new_mconf + mgrid( new_conf(j_tmp) )
            ENDIF
         END DO
      END IF
      
      new_binconf(1) = int( (new_econf - emin)/deltae ) + 1
      new_binconf(2) = int( (new_mconf - mmin)/deltam ) + 1
      
    END SUBROUTINE get_new_conf
    
    
    !----------------------------------------------------------------------------------------------!
    
    SUBROUTINE initialize
      IMPLICIT NONE
      
      REAL(kind=8) :: delta_tmp,rnd_tmp,emax_tmp
      REAL(kind=8) :: k_tmp,gapmin_tmp,gap_tmp,mmax
      INTEGER :: i_tmp,j_tmp,l1_tmp,l2_tmp, l0_tmp
      
      
      IF (mod(L,2) /= 0) THEN
         PRINT*,"Error: only even L are admitted"
         STOP
      END IF
      
      logf = dlog(f)
      
      IF ( emax < emin ) CALL error('initialize','error emax is smaller than emin')
      
      CALL random_seed(put=seed)
      
      ALLOCATE(conf(L))
      ALLOCATE(egrid(L))
      ALLOCATE(mgrid(L))
      
      delta_tmp = 2.0_dp*pi/dble(L)
      
      k_tmp = pi*(-1.0_dp+1.0_dp/dble(L))
      !
      egrid(1) = energy(k_tmp, h, gamma)
      econf = - egrid(1)/2.0_dp
      !
      mgrid(1) = sigmaz(k_tmp)
      mconf = - mgrid(1)/2.0_dp
      
      gapmin_tmp = egrid(1)
      emax_tmp = egrid(1)
      
      mmax = 0.0_dp
      mmin = 0.0_dp
      IF (mgrid(1)>0.0_dp) mmax = mgrid(1)
      IF (mgrid(1)<0.0_dp) mmin = mgrid(1)
      !
      !
      DO i_tmp=2,L
         !
         k_tmp = k_tmp+delta_tmp
         !energy associated to point k_tmp
         egrid(i_tmp) = energy(k_tmp, h, gamma)
         mgrid(i_tmp) = sigmaz(k_tmp)
         
         !minimum gap using single particle state
         !this method is good for gamma = 1 and h \= 0 where the energy is monotonic for k>0 or k<0
         !in this way at every mcstep the random wolker change bin
         gap_tmp = abs( egrid(i_tmp) - egrid(i_tmp-1) )
         IF ( (gap_tmp > 1e-12) .and. (gap_tmp < gapmin_tmp) ) THEN
            gapmin_tmp = gap_tmp
         END IF
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
         IF ( emax_tmp < egrid(i_tmp) )  emax_tmp = egrid(i_tmp)
      END DO
      
      mmin = mmin + mconf
      mmax = mmax + mconf
      
      !check if emax and emin are compatible with gs energy
      IF ((emax < econf/dble(L) ).or.(emin > -econf/dble(L) )) THEN
         PRINT*,"gs energy",econf/dble(L)
         CALL error("initialize", "energy interval outside energy spectra" )
      END IF
      
      
      IF ( ( econf/dble(L) < emax ) .and. ( econf/dble(L) > emin )) THEN
         GO TO 10
      END IF
      
      !check if i can always find a good initial state
      IF ( (emax-emin)<(2.0*emax_tmp/dble(L)) ) THEN
         PRINT*,"Warning:"," i am not sure to get an initial state"
      END IF
      
      
      conf = .false.
      l0_tmp = 0
      
      DO
         CALL random_number(rnd_tmp)
         j_tmp = int(rnd_tmp*(L-l0_tmp)) + 1
         
         l2_tmp = 0
         
         DO l1_tmp = 1, L
            
            IF ( .not.conf(l1_tmp) ) THEN
               l2_tmp = l2_tmp + 1
               IF (l2_tmp == j_tmp) THEN
                  conf(l1_tmp) = .true.
                  econf = econf + egrid(l1_tmp)
                  mconf = mconf + mgrid(l1_tmp)
                  l0_tmp = l0_tmp + 1
                  IF (mod(l0_tmp,2) == 0) THEN
                     IF ( ( econf/dble(L) < emax ) .and. ( econf/dble(L) > emin )) THEN
                        GO TO 10
                     END IF
                  END IF
                  GO TO 11
               END IF
            END IF
11          CONTINUE
            
         END DO
         
         IF (l0_tmp == L) EXIT
         
      END DO
      
      CALL error("intialize","I was not able to generate the initial configuration")
      
10    CONTINUE
      
      PRINT*,"finded configuration with energy",econf/dble(L)
      
      ! compute necessary number of bins
      gapmin_tmp = gapmin_tmp/dble(L)
      gapmin_tmp = 0.5_dp*gapmin_tmp
      !
      enbin = int((emax-emin)/gapmin_tmp)+1
      if (enbin < 2) CALL error("initialize","too small interval")
      
      ! energy per unit of volume
      egrid = egrid/dble(L)
      econf = econf/dble(L)
      mgrid = mgrid/dble(L)
      mconf = mconf/dble(L)
      mmax = mmax/dble(L)
      mmin = mmin/dble(L)
      
      !if the interval is symmetric respect to 0 i want to put 0 in the center
      !of a bin, therefore i need odd number of bin
      !pay attention : after i add one more bin
      !IF (((emax+emin)<1e-8).and.(mod(enbin,2) == 1)) THEN
      !   enbin = enbin +1
      !   PRINT*,"unlacky case"
      !END IF
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !i put one more bin because otherwise there are problem at the edge bin
      deltae = (emax-emin)/dble(enbin)
      emax = emax + deltae*0.5_dp
      emin = emin - deltae*0.5_dp
      enbin = enbin + 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      deltam = (emax-emin)/dble(mnbin)
      mmax = mmax + deltam*0.5_dp
      mmin = mmin - deltam*0.5_dp
      
      
      
      ALLOCATE(histo(enbin,mnbin))
      ALLOCATE(filled(enbin,mnbin))
      ALLOCATE(logdos(enbin,mnbin))
      
      histo = 0
      filled = 0
      logdos = 0.0_dp
      
      
      binconf(1) = int((econf-emin) / deltae )+1
      IF ((binconf(1)<1).or.(binconf(1)>enbin)) THEN
         CALL error("initialize", "energy interval outside energy spectra" )
      END IF
      
      binconf(2) = int((mconf-mmin) / deltam )+1
      IF ((binconf(2)<1).or.(binconf(2)>mnbin)) THEN
         CALL error("initialize", "magnetization outside its interval" )
      END IF
      
      PRINT*,enbin,"X",mnbin," bins"
      PRINT*
      PRINT*,"----------------- RUN IS STARTED -----------------"
      
      
    END SUBROUTINE initialize
    
    
    SUBROUTINE finalize
      IMPLICIT NONE
      !    CLOSE(unit=2)
      DEALLOCATE(conf,egrid,logdos,histo,filled)
      
    END SUBROUTINE finalize
    
    
    SUBROUTINE save
      
      IMPLICIT NONE
      
      INTEGER :: i_tmp,neffbin,ieff,j_tmp
      REAL(kind=8) :: alpha_tmp
      
      neffbin = 0
      alpha_tmp = -1.0
      DO i_tmp=1,enbin
         DO j_tmp=1,mnbin
            IF (filled(i_tmp,j_tmp) == 1) THEN
               neffbin = neffbin + 1
               IF (alpha_tmp < 0) alpha_tmp = logdos(i_tmp,j_tmp)
            END IF
         END DO
      END DO
      
      IF (allocated(logdosf)) deallocate(logdosf)
      
      ALLOCATE( logdosf( neffbin,3 ) )
      
      ieff = 0
      DO i_tmp=1,enbin
         DO j_tmp=1,mnbin
            
            IF (filled(i_tmp,j_tmp) == 1) THEN
               ieff = ieff + 1
               logdosf( ieff, 1) = emin+(dble(i_tmp)-0.5)*deltae
               logdosf( ieff, 2) = mmin+(dble(j_tmp)-0.5)*deltam
               logdosf( ieff, 3)  = logdos(i_tmp,j_tmp)-alpha_tmp
            END IF
         END DO
      END DO
      
    END SUBROUTINE save
    
  END SUBROUTINE wl_jdos
  !
  !
  !
  !
  !
  FUNCTION wl_mdos( ener, deltae , mnbin, seed, readconf )
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
    REAL (kind = 8) :: mmax,mmin
    REAL (kind = 8) :: logf,deltam
    !
    INTEGER :: binconf, new_binconf
    INTEGER :: i
    INTEGER :: i_swe, i_swe_tot, step
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
       !
       DO
          i_swe = i_swe + 1
          !
          DO i=1,L
             CALL get_new_conf
             CALL mc_step
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
             IF (check(histo,mnbin,i_swe*L)) THEN
                EXIT
             END IF
             !
             IF (i_swe>100000000) THEN
                PRINT*,"Some problem in the convergence"
                STOP 1
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
    
    PRINT*,"----------------- RUN IS ENDED -------------------"
    PRINT*
    print*,"Obtained convergence in ",i_swe_tot,"sweps"
    CALL save
    CALL finalize

!    close(unit=2)

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
      
    END FUNCTION check
    
    
    
    SUBROUTINE mc_step
      
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
         
         conf( new_conf ) = .not.conf( new_conf )
         
      END IF

!      WRITE(unit=2,fmt=*) (dabs(new_econf-ener) <= deltae)
      
10    RETURN
      
    END SUBROUTINE mc_step
    
    
    
    SUBROUTINE update
      IMPLICIT NONE
      
      logdos( binconf ) = logdos( binconf ) + logf
      histo ( binconf ) =  histo( binconf ) + 1
      
    END SUBROUTINE update
    
    
    
    SUBROUTINE get_new_conf
      IMPLICIT NONE
      
!      REAL(kind=8), DIMENSION(2) :: rnd_tmp
      REAL(kind=8) :: rnd_tmp
!      INTEGER :: j_tmp
      
      new_econf = econf
      new_mconf = mconf


      CALL random_number( rnd_tmp )
      new_conf = int(dble(L)*rnd_tmp) + 1
      
!      IF (new_conf(1) /= new_conf(2)) THEN
!         DO j_tmp =1,2
!            IF (conf( new_conf(j_tmp) )) THEN
!               new_econf = new_econf - egrid( new_conf(j_tmp) )
!               new_mconf = new_mconf - mgrid( new_conf(j_tmp) )
!            ELSE
!               new_econf = new_econf + egrid( new_conf(j_tmp) )
!               new_mconf = new_mconf + mgrid( new_conf(j_tmp) )
!            ENDIF
!         END DO
!      END IF
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
         j_tmp = int(rnd_tmp*(L-l0_tmp)) + 1
         IF (j_tmp>L) j_tmp = L
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
                  IF (mod(l0_tmp,2) == 0) THEN
                     print*,"econf :",econf/dble(L)
                     IF ( dabs( econf/dble(L) - ener ) <= deltae ) THEN
                        CALL save_conf(conf,L)
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
  !
  ! it does a check of the observable
  SUBROUTINE check_obs( conf, obsConf, obs , dim )
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: dim
    LOGICAL, DIMENSION(dim), INTENT(IN) :: conf
    REAL(kind=8), INTENT(INOUT) :: obsConf
    REAL(kind=8),DIMENSION(dim), INTENT(IN) :: obs
    INTEGER :: i
    REAL(kind=8) :: obsConf_check

    obsConf_check = 0.0d0
    DO i=1,dim

       IF (conf(i)) THEN
          obsConf_check = obsConf_check + 0.5d0*obs(i)
       ELSE
          obsConf_check = obsConf_check - 0.5d0*obs(i)
       ENDIF

    ENDDO

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
END MODULE wl
