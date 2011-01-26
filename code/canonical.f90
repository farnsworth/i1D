MODULE canonical
  !
  ! temp = k_b * T / J
  !
  USE system
  !
  !
CONTAINS
  !
  !
  FUNCTION thermal_ex_energy( temperature )
    !
    IMPLICIT NONE
    !
    REAL( kind = 8), INTENT(IN) :: temperature
    REAL( kind = 8) :: k,e, thermal_ex_energy
    INTEGER :: i
    !
    thermal_ex_energy = 0.0d0
    !
    ! here i use antiperiodic boundary condition
    !
    DO i=1,L/2
       !
       !
       k = ( pi*dble(2*i-1) ) / dble(L)
       e = energy( k, h, gamma )
       !
       thermal_ex_energy = thermal_ex_energy + e *2.0d0*fermi_dist(e, temperature)
       !
    END DO
    !
    thermal_ex_energy = thermal_ex_energy/dble(L)
    !
  END FUNCTION thermal_ex_energy
  !
  !
  ! ... compute the canonical magnetization at temperature T
  ! ... using antiperiodic boundary condition on c^+ and c
  ! ... and without any restriction on K tot
  !
  FUNCTION thermal_magnetization( temperature )
    !
    IMPLICIT NONE
    !
    REAL (kind = dp), INTENT(IN) :: temperature
    REAL (kind = dp) :: k, thermal_magnetization,mzk,ek
    INTEGER :: i
    !
    thermal_magnetization = 0.0d0
    !
    ! here i use antiperiodic boundary condition
    ! without any condition on total momentum
    !
    DO i=1,L/2
       !
       !
       k = ( pi*dble(2*i-1) ) / dble(L)
       mzk = sigmaz( k )
       ek = energy( k, h, gamma )
       !
       thermal_magnetization = thermal_magnetization + mzk * ( 2.0d0*fermi_dist( ek, temperature ) - 1.0d0 )
       !
    END DO
    !
    thermal_magnetization = thermal_magnetization/dble(L)
    !
  END FUNCTION thermal_magnetization
  !
  !
  !
  ! ... compute the canonical magnetization at temperature T
  ! ... using antiperiodic boundary condition on c^+ and c
  ! ... and without any restriction on K tot
  !
  FUNCTION thermal_kinks( temperature )
    !
    IMPLICIT NONE
    !
    REAL (kind = dp), INTENT(IN) :: temperature
    REAL (kind = dp) :: k, thermal_kinks,kinksk,ek
    INTEGER :: i
    !
    thermal_kinks = 0.0d0
    !
    DO i=1,L/2
       !
       k = ( pi*dble(2*i-1) ) / dble(L)
       kinksk = kinks( k )
       ek = energy( k, h, gamma )
       !
       thermal_kinks = thermal_kinks + kinksk * ( 2.0d0*fermi_dist( ek, temperature ) - 1.0d0 )
       !
    END DO
    !
    thermal_kinks = thermal_kinks/dble(L) + 0.5d0
    !
  END FUNCTION thermal_kinks
  !
  !
  FUNCTION fermi_dist ( en, temp )
    !
    IMPLICIT NONE
    !
    REAL(kind = dp), INTENT(IN) :: en,temp
    REAL(kind = dp) :: fermi_dist
    !
    IF (temp == 0.0d0) THEN
       fermi_dist = 0.0d0
    ELSE
       fermi_dist = 1.0d0/( dexp(en/temp) + 1.0d0 )
    END IF
    !
  END FUNCTION fermi_dist
  !
  !
  ! la temperatura che sputa fuori è meta o doppia di quella
  ! che vogliamo noi, da quella canonica classica
  !
  !
  FUNCTION find_temperature( en, eps )
    !
    IMPLICIT NONE
    !
    REAL(kind = dp), INTENT(IN) :: en, eps
    REAL(kind = dp) :: find_temperature, emin, etemp, Tmin, Tmax, Ttemp,en_ex
    !
    !
    ! ... work with thermal exitation energy: more precise
    ! ... at low temperature
    !
    en_ex = en - Egs()
    emin = 0.0d0
    !
    !
    IF ( (en_ex < 0.0d0) .or. (en_ex >= -Egs() ) ) THEN
       !
       CALL error ('find_temperature','energy out of range')
       find_temperature = -1.0
       RETURN
       !
    END IF
    !
    Tmin = 0.0d0
    !
    !
    ! ... a way to find Tmax
    !
    Tmax = 1.0d0
    !
    ! ... find temperature interval for which en is inside 
    !
    DO WHILE ( Tmax < 1e10 )
       !
       etemp = thermal_ex_energy( Tmax )
       IF (etemp > en_ex) THEN
          EXIT
       ELSE
          Tmin = Tmax
          Tmax = 2.0d0*Tmax
       END IF
       !
    END DO
    !
    IF ( Tmax > 1e10 ) THEN
       CALL error('find_temperature','too high energy')
       find_temperature = -1.0
       RETURN
    END IF
    !
    ! ... bisection method
    !
    DO WHILE ( ( (Tmax-Tmin)/Tmax ) > eps )
       !
       Ttemp = (Tmax + Tmin)/2.0d0
       etemp = thermal_ex_energy( Ttemp )
       !
       IF ( abs(etemp - en_ex) < 1e-12 ) THEN
          CALL warning('find_temperature', "energies are too similar, I&
          & can not go beyond")
          EXIT
       ENDIF
       !
       IF (etemp > en_ex ) THEN
          Tmax = Ttemp
       ELSE
          Tmin = Ttemp
       ENDIF
       !
    END DO
    !
    print*,"error on temperature",(Tmax - Tmin)/2.0
    !
    find_temperature = (Tmax + Tmin) / 2.0d0
    !
    !
  END FUNCTION find_temperature
  !
  !
  FUNCTION thermal_BiAj (d,temperature)
    !
    INTEGER, INTENT(IN) :: d
    REAL (kind = dp), INTENT(IN) :: temperature
    REAL (kind = 8) :: thermal_BiAj, term, k,ek
    !
    thermal_BiAj = 0.0_dp
    !
    DO i_tmp=1,L/2
       !
       k = (pi*dble(2*i_tmp-1))/dble(L)
       ek = energy(k,h,gamma)
       !
       term = (dcos(k*dble(d-1))-h*dcos(k*dble(d))) / ek
       !
       thermal_BiAj = thermal_BiAj + term * ( 1.0d0 - 2.0d0*fermi_dist( ek, temperature ) )
       !
    ENDDO
    !
    thermal_BiAj = thermal_BiAj*4.0d0/dble(L)
    !
    !
  END FUNCTION thermal_BiAj
  !
  !
  FUNCTION thermal_xxcorrelation (d, temperature)
    !
    INTEGER, INTENT(IN) :: d
    REAL(kind = dp),INTENT(IN) :: temperature
    REAL(kind = dp), DIMENSION ( d, d ) :: matrix
    REAL(kind=dp), DIMENSION (-d+2:d ) :: vectorBiAj
    INTEGER :: j,k,info,sign
    INTEGER,DIMENSION(d) :: pivot
    REAL(kind = dp) :: det
    !
    REAL(kind= 8) :: thermal_xxcorrelation
    !
    !
    matrix = 0.0_dp
    det = 1.0_dp
    !
    DO k=-d+2,d
       !
       vectorBiAj(k) = thermal_BiAj(k, temperature)
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
    thermal_xxcorrelation = det
    !
  END FUNCTION THERMAL_XXCORRELATION
  !
  !
END MODULE canonical
