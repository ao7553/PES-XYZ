program potential_xyz
  implicit none

  integer, parameter :: ark = selected_real_kind(25,32)
  real(ark) :: deg, pi
  integer ipar, ieq, parmax, info
  double precision  f, local(3), par(9), force(84)
  double precision  r1, r2, alpha
  character(9) buf30

  interface
  function pes_xyz(local,parmax,par,force)
    integer parmax
    double precision local(3), par(9), force(84), pes_xyz
  end function pes_xyz
  end interface

  pi = 4.0d0 * datan2(1.0d0,1.0d0) ! Value of pi
  deg = pi/180.0d0 ! Conversion factor degrees to radians

  ! Read number of potential energy surface expansion parameters
  read*, parmax

  ! Read equilibrium geometry parameters
  do ieq =1, 9
    read*, buf30, par(ieq)
  enddo

  ! Read potential energy surface expansion parameters
  do ipar=1, parmax
    read*, buf30, force(ipar)
  enddo

  ! Read grid of user-defined geometries and print potential energy of molecule
  do
    read(*,*,end=4,err=4) local(1:3)
    r1  = local(1) ! r1 internal stretching coordinate in Angstrom
    r2  = local(2) ! r2 internal stretching coordinate in Angstrom
    alpha  = local(3) ! Bending internal coordinate in degrees
    f = pes_xyz(local,parmax,par,force)
    write(*,'(3(1x,f10.4),2x,f14.6)') local(1:3), f
    cycle
4  exit
  enddo

end program potential_xyz


function pes_xyz(local,parmax,par,force) result (f)
  implicit none

  integer, parameter :: ark = selected_real_kind(25,32)
  integer parmax, ipar
  double precision local(3), force(84), par(9)
  double precision f

  real(ark) :: y1, y2, y3, r1e, r2e, alphae, deg, pi, a01, a02, rho, rhoe
  real(ark) :: r1,r2,alpha,v0,v1,v2,v3,v4,v5,v6,b1,b2,g1,g2,roo,voo

  real(ark), parameter :: tocm = 219474.63067_ark

  pi = 4.0d0 * datan2(1.0d0,1.0d0) ! Value of pi
  deg = pi/180.0d0 ! Conversion factor degrees to radians

  r1e = par(1)
  r2e = par(2)
  alphae = par(3) * deg
  a01 = par(4)
  a02 = par(5)
  b1 = par(6)
  b2 = par(7)
  g1 = par(8)
  g2 = par(9)

  r1 = local(1)
  r2 = local(2)
  alpha = local(3) * deg
  rhoe = pi - alphae
  rho = pi - alpha

  y1 = 1.0_ark - exp(-a01 * (r1 - r1e)) ! 1st vibrational stretching coordinate
  y2 = 1.0_ark - exp(-a02 * (r2 - r2e)) ! 2nd vibrational stretching coordinate
  y3 = sin(rho) - sin(rhoe) ! Bending coordinate

  roo = sqrt(r1**2 + r2**2 - 2.0_ark * r1 * r2 * cos(alpha)) ! Distance between atoms
  voo = b1 * exp(-g1 * roo) + b2 * exp(-g2 * roo**2) ! Damping term in potential when roo small

  v0 = force(1) * y1**0 * y2**0 * y3**0
  v1 = force(2) * y1**0 * y2**0 * y3**1 &
     + force(3) * y1**1 * y2**0 * y3**0 &
     + force(4) * y1**0 * y2**1 * y3**0
  v2 = force(5) * y1**0 * y2**0 * y3**2 &
     + force(6) * y1**1 * y2**0 * y3**1 &
     + force(7) * y1**0 * y2**1 * y3**1 &
     + force(8) * y1**1 * y2**1 * y3**0 &
     + force(9) * y1**2 * y2**0 * y3**0 &
     + force(10) * y1**0 * y2**2 * y3**0
  v3 = force(11) * y1**0 * y2**0 * y3**3 &
     + force(12) * y1**1 * y2**0 * y3**2 &
     + force(13) * y1**0 * y2**1 * y3**2 &
     + force(14) * y1**1 * y2**1 * y3**1 &
     + force(15) * y1**2 * y2**0 * y3**1 &
     + force(16) * y1**0 * y2**2 * y3**1 &
     + force(17) * y1**2 * y2**1 * y3**0 &
     + force(18) * y1**1 * y2**2 * y3**0 &
     + force(19) * y1**3 * y2**0 * y3**0 &
     + force(20) * y1**0 * y2**3 * y3**0
  v4 = force(21) * y1**0 * y2**0 * y3**4 &
     + force(22) * y1**1 * y2**0 * y3**3 &
     + force(23) * y1**0 * y2**1 * y3**3 &
     + force(24) * y1**1 * y2**1 * y3**2 &
     + force(25) * y1**2 * y2**0 * y3**2 &
     + force(26) * y1**0 * y2**2 * y3**2 &
     + force(27) * y1**2 * y2**1 * y3**1 &
     + force(28) * y1**1 * y2**2 * y3**1 &
     + force(29) * y1**2 * y2**2 * y3**0 &
     + force(30) * y1**3 * y2**0 * y3**1 &
     + force(31) * y1**0 * y2**3 * y3**1 &
     + force(32) * y1**3 * y2**1 * y3**0 &
     + force(33) * y1**1 * y2**3 * y3**0 &
     + force(34) * y1**4 * y2**0 * y3**0 &
     + force(35) * y1**0 * y2**4 * y3**0
  v5 = force(36) * y1**0 * y2**0 * y3**5 &
     + force(37) * y1**1 * y2**0 * y3**4 &
     + force(38) * y1**0 * y2**1 * y3**4 &
     + force(39) * y1**1 * y2**1 * y3**3 &
     + force(40) * y1**2 * y2**0 * y3**3 &
     + force(41) * y1**0 * y2**2 * y3**3 &
     + force(42) * y1**2 * y2**1 * y3**2 &
     + force(43) * y1**1 * y2**2 * y3**2 &
     + force(44) * y1**2 * y2**2 * y3**1 &
     + force(45) * y1**3 * y2**0 * y3**2 &
     + force(46) * y1**0 * y2**3 * y3**2 &
     + force(47) * y1**3 * y2**1 * y3**1 &
     + force(48) * y1**1 * y2**3 * y3**1 &
     + force(49) * y1**3 * y2**2 * y3**0 &
     + force(50) * y1**2 * y2**3 * y3**0 &
     + force(51) * y1**4 * y2**0 * y3**1 &
     + force(52) * y1**0 * y2**4 * y3**1 &
     + force(53) * y1**4 * y2**1 * y3**0 &
     + force(54) * y1**1 * y2**4 * y3**0 &
     + force(55) * y1**5 * y2**0 * y3**0 &
     + force(56) * y1**0 * y2**5 * y3**0
  v6 = force(57) * y1**0 * y2**0 * y3**6 &
     + force(58) * y1**1 * y2**0 * y3**5 &
     + force(59) * y1**0 * y2**1 * y3**5 &
     + force(60) * y1**1 * y2**1 * y3**4 &
     + force(61) * y1**2 * y2**0 * y3**4 &
     + force(62) * y1**0 * y2**2 * y3**4 &
     + force(63) * y1**2 * y2**1 * y3**3 &
     + force(64) * y1**1 * y2**2 * y3**3 &
     + force(65) * y1**2 * y2**2 * y3**2 &
     + force(66) * y1**3 * y2**0 * y3**3 &
     + force(67) * y1**0 * y2**3 * y3**3 &
     + force(68) * y1**3 * y2**1 * y3**2 &
     + force(69) * y1**1 * y2**3 * y3**2 &
     + force(70) * y1**3 * y2**2 * y3**1 &
     + force(71) * y1**2 * y2**3 * y3**1 &
     + force(72) * y1**3 * y2**3 * y3**0 &
     + force(73) * y1**4 * y2**0 * y3**2 &
     + force(74) * y1**0 * y2**4 * y3**2 &
     + force(75) * y1**4 * y2**1 * y3**1 &
     + force(76) * y1**1 * y2**4 * y3**1 &
     + force(77) * y1**4 * y2**2 * y3**0 &
     + force(78) * y1**2 * y2**4 * y3**0 &
     + force(79) * y1**5 * y2**0 * y3**1 &
     + force(80) * y1**0 * y2**5 * y3**1 &
     + force(81) * y1**5 * y2**1 * y3**0 &
     + force(82) * y1**1 * y2**5 * y3**0 &
     + force(83) * y1**6 * y2**0 * y3**0 &
     + force(84) * y1**0 * y2**6 * y3**0

  f = v0 + v1 + v2 + v3 + v4 + v5 + v6 + voo ! Total potential energy

end function pes_xyz
