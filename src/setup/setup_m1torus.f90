!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Self-gravitating collisionless torus around a live central mass, built from
! Keplerian elements as in Bannikova et al. (2026, Sect. 2): semi-major axes
! uniform in [1 - R0, 1 + R0], eccentricities in [0, e_max], inclinations in
! [-i_max, i_max], node and true anomaly random, argument of periapsis zero.
! Units G = M_c = R_tor = 1; the orbit at R_tor lasts 2 pi.
!
! :References:
!   - Bannikova et al. (2026), A&A, arXiv:2604.11528
!
! :Owner: etrevol
!
! :Runtime parameters:
!   - M_tor            : *torus mass [M_c]*
!   - R0               : *half-width of the semi-major axis range [R_tor]*
!   - e_max            : *maximum initial eccentricity*
!   - i_max            : *maximum initial inclination [deg]*
!   - omega_random     : *random argument of periapsis (paper's control run)*
!   - np               : *number of torus particles*
!   - eps_c            : *softening length of the central mass (0 = Kepler)*
!   - eps_star         : *fixed Plummer-equivalent softening of the torus (0 = adaptive)*
!   - r_acc            : *accretion radius of the central mass*
!   - iseed            : *random seed*
!   - n_orbits         : *run length [orbits at R_tor]*
!   - n_analysis_orbit : *analysis calls per orbit*
!   - n_dumps_orbit    : *full dumps per orbit*
!
! :Dependencies: infile_utils, io, io_control, options, part, physcon,
!   random, timestep, units
!
 implicit none
 public :: setpart

 private

 real, parameter :: pi = 4.*atan(1.)

 real    :: M_tor = 0.1, R0 = 0.3, e_max = 0.5, i_max = 60.
 logical :: omega_random = .false.
 integer :: np = 131072, iseed = 1
 ! eps_c > 0 makes phantom soften the particle-sink force over max(eps_c, h_i)
 ! (ptmass.F90) with the particle's adaptive h_i, and that potential changes with
 ! h_i without a correction term: in the empty centre h_i ~ 0.2, energy drifted
 ! 1.4e-4 per orbit. eps_c = 0 gives the Kepler force, as the paper's eps = 0.01
 ! does for pericentres >= 0.05, and a drift of 1.7e-5.
 real    :: eps_c = 0., r_acc = 1.e-3
 ! eps_star > 0 fixes the torus softening (paper: Plummer eps = 0.01) instead of phantom's
 ! adaptive h; the cubic spline matches a Plummer potential at r = 0 for h = 1.4 eps
 real    :: eps_star = 0.
 real    :: n_orbits = 300.
 integer :: n_analysis_orbit = 4, n_dumps_orbit = 1

contains

subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,&
                   gamma,hfact,time,fileprefix)
 use units,        only:set_units
 use physcon,      only:pc,solarm
 use part,         only:istar,set_particle_type,nptmass,xyzmh_ptmass,vxyz_ptmass, &
                        ihacc,ihsoft,hsoft_nbody
 use options,      only:iexternalforce
 use timestep,     only:tmax,dtmax
 use io_control,   only:nout,nfulldump
 use random,       only:ran2
 use infile_utils, only:get_options
 use io,           only:master,fatal
 integer,           intent(in)  :: id
 integer,           intent(out) :: npart
 integer,           intent(out) :: npartoftype(:)
 real,              intent(out) :: xyzh(:,:)
 real,              intent(out) :: polyk,gamma,hfact
 real,              intent(out) :: vxyzu(:,:)
 real,              intent(out) :: massoftype(:)
 real,              intent(in)  :: time
 character(len=20), intent(in)  :: fileprefix
 integer :: ierr,i,is
 real :: pmass,a,e,inc,node,omega,nu,p,r,xo,yo,vxo,vyo,h0,vol
 real :: c1,s1,ci,si,co,so,xcm(3),vcm(3),mtot,p_orb

 call get_options(trim(fileprefix)//'.setup',id==master,ierr,&
                  read_setupfile,write_setupfile)
 if (ierr /= 0) stop 'rerun phantomsetup after editing .setup file'
 if (np > size(xyzh,2)) call fatal('setup','np > maxp: rerun with --maxp')

 ! one length unit = 1 pc and M_c = 1e8 M_sun as in the paper's M31 scaling; the
 ! dynamics are scale free
 call set_units(dist=pc,mass=1.d8*solarm,G=1.d0)

 gamma = 1.
 polyk = 0.
 hfact = 1.2
 iexternalforce = 0

 npart = np
 pmass = M_tor/np
 npartoftype(:) = 0
 npartoftype(istar) = np
 massoftype(:) = 0.
 massoftype(istar) = pmass

 ! initial guess of h from the mean density of a ring of radius 1, half-width R0
 ! and height sin(i_max); phantom iterates it to its own value
 vol = 2.*pi*(2.*R0)*(2.*max(sin(i_max*pi/180.),R0))
 h0  = hfact*(pmass*vol/M_tor)**(1./3.)
 if (eps_star > 0.) then
    hsoft_nbody = 1.4*eps_star
    h0 = hsoft_nbody
 endif

 is = -abs(iseed)
 do i=1,np
    a     = (1. - R0) + 2.*R0*ran2(is)
    e     = e_max*ran2(is)
    inc   = (2.*ran2(is) - 1.)*i_max*pi/180.
    node  = 2.*pi*ran2(is)
    nu    = 2.*pi*ran2(is)
    omega = 0.
    if (omega_random) omega = 2.*pi*ran2(is)

    ! perifocal frame about M_c = 1
    p   = a*(1. - e*e)
    r   = p/(1. + e*cos(nu))
    xo  = r*cos(nu)
    yo  = r*sin(nu)
    vxo = -sqrt(1./p)*sin(nu)
    vyo =  sqrt(1./p)*(e + cos(nu))

    ! rotate by omega, inclination, node: R = Rz(node) Rx(inc) Rz(omega)
    c1 = cos(node)
    s1 = sin(node)
    ci = cos(inc)
    si = sin(inc)
    co = cos(omega)
    so = sin(omega)
    xyzh(1,i) = (c1*co - s1*so*ci)*xo + (-c1*so - s1*co*ci)*yo
    xyzh(2,i) = (s1*co + c1*so*ci)*xo + (-s1*so + c1*co*ci)*yo
    xyzh(3,i) = (so*si)*xo + (co*si)*yo
    xyzh(4,i) = h0
    vxyzu(1,i) = (c1*co - s1*so*ci)*vxo + (-c1*so - s1*co*ci)*vyo
    vxyzu(2,i) = (s1*co + c1*so*ci)*vxo + (-s1*so + c1*co*ci)*vyo
    vxyzu(3,i) = (so*si)*vxo + (co*si)*vyo
    call set_particle_type(i,istar)
 enddo

 nptmass = 1
 xyzmh_ptmass(:,1) = 0.
 xyzmh_ptmass(4,1) = 1.
 xyzmh_ptmass(ihacc,1)  = r_acc
 xyzmh_ptmass(ihsoft,1) = eps_c
 vxyz_ptmass(:,1) = 0.

 ! the system barycentre at rest at the origin, so that r_c = -M_tor r_tb holds
 mtot = 1. + M_tor
 do i=1,3
    xcm(i) = pmass*sum(xyzh(i,1:np))/mtot
    vcm(i) = pmass*sum(vxyzu(i,1:np))/mtot
    xyzh(i,1:np)  = xyzh(i,1:np) - xcm(i)
    vxyzu(i,1:np) = vxyzu(i,1:np) - vcm(i)
    xyzmh_ptmass(i,1) = xyzmh_ptmass(i,1) - xcm(i)
    vxyz_ptmass(i,1)  = vxyz_ptmass(i,1) - vcm(i)
 enddo

 p_orb     = 2.*pi
 tmax      = n_orbits*p_orb
 dtmax     = p_orb/n_analysis_orbit
 nout      = n_analysis_orbit/n_dumps_orbit
 nfulldump = 1

 if (id==master) then
    print "(/,a)",' Collisionless torus (Bannikova et al. 2026)'
    print "(a,i10,a,es10.3,a,f6.3,a,f6.3,a,f6.1)",'  N = ',np,'  M_tor = ',M_tor, &
          '  R0 = ',R0,'  e_max = ',e_max,'  i_max = ',i_max
    print "(a,es10.3,a,es10.3,a,es10.3,a,es10.3)",'  eps_c = ',eps_c,'  r_acc = ',r_acc, &
          '  h0 = ',h0,'  eps_star = ',eps_star
    print "(a,3es11.3)",'  barycentre shift removed: ',xcm
    print "(a,f8.1,a,i4,a,i4,a)",'  tmax = ',n_orbits,' orbits, ',n_analysis_orbit, &
          ' analysis calls and ',n_dumps_orbit,' dumps per orbit'
 endif

end subroutine setpart

subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for the collisionless torus setup (m1torus)'
 write(iunit,"(/,a)") '# torus (Keplerian elements about the central mass)'
 call write_inopt(M_tor,'M_tor','torus mass [M_c]',iunit)
 call write_inopt(R0,'R0','half-width of the semi-major axis range [R_tor]',iunit)
 call write_inopt(e_max,'e_max','maximum initial eccentricity',iunit)
 call write_inopt(i_max,'i_max','maximum initial inclination [deg]',iunit)
 call write_inopt(omega_random,'omega_random','random argument of periapsis',iunit)
 call write_inopt(np,'np','number of torus particles',iunit)
 call write_inopt(iseed,'iseed','random seed',iunit)
 call write_inopt(eps_star,'eps_star','fixed Plummer softening of the torus (0 = adaptive)',iunit)
 write(iunit,"(/,a)") '# central mass'
 call write_inopt(eps_c,'eps_c','softening length of the central mass (0 = Kepler)',iunit)
 call write_inopt(r_acc,'r_acc','accretion radius of the central mass',iunit)
 write(iunit,"(/,a)") '# run'
 call write_inopt(n_orbits,'n_orbits','run length [orbits at R_tor]',iunit)
 call write_inopt(n_analysis_orbit,'n_analysis_orbit','analysis calls per orbit',iunit)
 call write_inopt(n_dumps_orbit,'n_dumps_orbit','full dumps per orbit',iunit)
 close(iunit)

end subroutine write_setupfile

subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use io,           only:error
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)

 nerr = 0
 ierr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(M_tor,'M_tor',db,min=0.,errcount=nerr)
 call read_inopt(R0,'R0',db,min=0.,max=1.,errcount=nerr)
 call read_inopt(e_max,'e_max',db,min=0.,max=1.,errcount=nerr)
 call read_inopt(i_max,'i_max',db,min=0.,max=90.,errcount=nerr)
 call read_inopt(omega_random,'omega_random',db,errcount=nerr)
 call read_inopt(np,'np',db,min=1,errcount=nerr)
 call read_inopt(iseed,'iseed',db,errcount=nerr)
 call read_inopt(eps_c,'eps_c',db,min=0.,errcount=nerr)
 call read_inopt(eps_star,'eps_star',db,min=0.)
 call read_inopt(r_acc,'r_acc',db,min=0.,errcount=nerr)
 call read_inopt(n_orbits,'n_orbits',db,min=0.,errcount=nerr)
 call read_inopt(n_analysis_orbit,'n_analysis_orbit',db,min=1,errcount=nerr)
 call read_inopt(n_dumps_orbit,'n_dumps_orbit',db,min=1,errcount=nerr)
 call close_db(db)
 if (mod(n_analysis_orbit,n_dumps_orbit) /= 0) then
    call error('setup','n_dumps_orbit must divide n_analysis_orbit')
    nerr = nerr + 1
 endif
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

end module setup
