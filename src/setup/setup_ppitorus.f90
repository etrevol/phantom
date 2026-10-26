!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Papaloizou-Pringle torus (constant specific angular momentum, polytropic)
! around a fixed point mass, unseeded or with a controlled seed for the
! non-axisymmetric instability.
!
! Units follow the Athena++ torus of athena-collab: L_0 = chi r_g, mass M_bh,
! G = 1, so the density maximum sits at r_center and the orbit there is
! 2 pi r_center^(3/2). The torus is not self-gravitating; rho_0 sets the
! mass scale only.
!
! :References:
!   - Papaloizou & Pringle (1984), MNRAS 208, 721
!
! :Owner: etrevol
!
! :Runtime parameters:
!   - C_prime          : *surface constant (0 < C_prime < 0.5)*
!   - M_bh             : *central mass [M_sun]*
!   - chi              : *length unit L_0 = chi r_g*
!   - gamma            : *adiabatic index*
!   - lattice          : *rings | random*
!   - n_analysis_orbit : *analysis calls per orbit*
!   - n_dumps_orbit    : *full dumps per orbit*
!   - n_orbits         : *run length [orbits at r_center]*
!   - np               : *target number of particles*
!   - nsym             : *ring counts are multiples of nsym (1 = no symmetry)*
!   - pert_amp         : *delta v_R in units of the local c_s*
!   - pert_kind        : *none | single | multi | noise*
!   - pert_m           : *azimuthal mode (single) or highest mode (multi)*
!   - pert_seed        : *random seed (multi, noise, random)*
!   - relax            : *relax before spinning up (needs r_center = 1)*
!   - n_relax          : *relaxation time [orbits at r_center]*
!   - damp_relax       : *velocity damping rate [code units]*
!   - r_acc            : *accretion radius of the central mass [L_0]*
!   - r_center         : *radius of the density maximum [L_0]*
!   - rho_0            : *density at the maximum [g/cm^3]*
!
! :Dependencies: damping, eos, externalforces, infile_utils, io, io_control,
!   options, part, random, shock_capturing, timestep, units
!
 implicit none
 public :: setpart,spin_up

 private

 ! physcon's pi has 11 digits: rings of 2 pi/N would not close, and the gap seeds
 ! every m at 3e-12, far above round-off
 real, parameter :: pi = 4.*atan(1.)

 ! CGS constants of acc_disk_visc.cpp, so L_0 matches the grid runs bit for bit
 real(kind=8), parameter :: c_cgs = 2.99792458d10, g_cgs = 6.67430d-8
 real(kind=8), parameter :: msun_cgs = 1.98841d33

 real    :: M_bh = 4.5e7, chi = 500., rho_0 = 1.e-13
 real    :: r_center = 1., C_prime = 0.2, gamma_torus = 1.3, r_acc = 0.35
 integer :: np = 100000, nsym = 32
 character(len=20) :: lattice = 'rings', pert_kind = 'none'
 real    :: pert_amp = 0.
 integer :: pert_m = 2, pert_seed = 1
 logical :: relax = .true.
 real    :: n_relax = 2., damp_relax = 0.3
 real    :: n_orbits = 40.
 integer :: n_analysis_orbit = 20, n_dumps_orbit = 2

 real :: polyn, f_c, rho_c, r_in, r_out, z_top, pmass, p_orb
 integer :: iseed

contains

! With relax, stage 1 of two: the torus at rest in the effective potential
! -G M/r + l^2/(2 R^2) of iexternalforce = 5, velocities damped, so the particles
! settle into the SPH equilibrium rather than the analytic one; spin_up (called by
! moddump_ppitorus) then gives them v_phi = l/R and the seed. Without relax the
! torus starts spinning on the analytic profile.
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,&
                   gamma,hfact,time,fileprefix)
 use units,          only:set_units,umass,udist,utime,unit_density
 use part,           only:igas
 use options,        only:iexternalforce
 use externalforces, only:iext_torus
 use damping,        only:idamp,damp
 use timestep,       only:tmax
 use eos,            only:ieos,gmw
 use infile_utils,   only:get_options
 use io,             only:master,fatal
 integer,           intent(in)  :: id
 integer,           intent(out) :: npart
 integer,           intent(out) :: npartoftype(:)
 real,              intent(out) :: xyzh(:,:)
 real,              intent(out) :: polyk,gamma,hfact
 real,              intent(out) :: vxyzu(:,:)
 real,              intent(out) :: massoftype(:)
 real,              intent(in)  :: time
 character(len=20), intent(in)  :: fileprefix
 integer :: ierr,i
 real(kind=8) :: r_g
 real :: m_torus,m_placed

 call get_options(trim(fileprefix)//'.setup',id==master,ierr,&
                  read_setupfile,write_setupfile)
 if (ierr /= 0) stop 'rerun phantomsetup after editing .setup file'

 r_g = 2.d0*g_cgs*M_bh*msun_cgs/c_cgs**2
 call set_units(dist=chi*r_g,mass=M_bh*msun_cgs,G=1.d0)

 gamma = gamma_torus
 hfact = 1.2
 polyk = 0.
 ieos  = 2
 gmw   = 0.6     ! as athena-collab; enters only phantom's temperature output
 call derive_shape()
 z_top = surface_height_max()
 rho_c = real(rho_0/unit_density)
 m_torus = rho_c*torus_volume()
 pmass   = m_torus/real(np)

 npart = 0
 select case(trim(lattice))
 case('rings')
    call set_rings(npart,xyzh,hfact)
 case('random')
    call set_random(npart,xyzh,hfact)
 case default
    call fatal('setup','unknown lattice '//trim(lattice))
 end select
 if (npart > size(xyzh,2)) call fatal('setup','npart > maxp: rerun with --maxp')

 ! Ring masses are matched to the analytic column, so the mass is pmass*npart
 ! up to the surface rings; keep pmass rather than rescaling the density.
 m_placed = pmass*npart
 npartoftype(:) = 0
 npartoftype(igas) = npart
 massoftype(igas)  = pmass

 iseed = -abs(pert_seed) - 1
 do i=1,npart
    call set_velocity_and_energy(xyzh(1:3,i),vxyzu(:,i),gamma,spin=.not.relax, &
                                 set_u=.true.)
 enddo

 call set_run_options()
 if (relax) then
    ! the potential of iext_torus is written for G M = 1, R_torus = 1
    if (abs(r_center - 1.) > epsilon(1.)) call fatal('setup','relax needs r_center = 1')
    iexternalforce = iext_torus
    idamp = 1
    damp  = damp_relax
    tmax  = n_relax*p_orb
 endif

 if (id==master) then
    print "(/,a)",' Papaloizou-Pringle torus'
    print "(a,f8.5,a,f8.5,a,f8.5,a,f9.6)",'  C'' = ',C_prime,'  r_in = ',r_in, &
          '  r_out = ',r_out,'  z_top = ',z_top
    print "(a,es12.5,a,es12.5,a,f7.4)",'  M_torus = ',m_torus,' placed = ',m_placed, &
          '  error = ',m_placed/m_torus - 1.
    print "(a,i10,a,i10,a,a,a,i4)",'  npart = ',npart,'  target = ',np, &
          '  lattice = ',trim(lattice),'  nsym = ',nsym
    print "(a,es12.5,a,es12.5,a)",'  L_0 = ',udist,' cm   t_unit = ',utime,' s'
    print "(a,es12.5,a,es12.5,a)",'  M_unit = ',umass,' g   rho_c = ',rho_c,' code'
    print "(a,f10.6,a,f8.2,a)",'  P_orb = ',p_orb,'   tmax = ',tmax/p_orb,' orbits'
    if (relax) then
       print "(a,f6.2,a,f6.3)",'  stage 1: relaxation at rest for',n_relax, &
             ' orbits, damp = ',damp_relax
    elseif (trim(pert_kind) /= 'none') then
       print "(a,a,a,es10.3,a,i3)",'  seed: ',trim(pert_kind),'  amp = ',pert_amp,'  m = ',pert_m
    else
       print "(a)",'  seed: none'
    endif
 endif

end subroutine setpart

! Stage 2: rotation l/R and the seed on the relaxed particles, point mass on,
! damping off. u is kept: it is the relaxed equilibrium, not the analytic one.
subroutine spin_up(npart,xyzh,vxyzu,setupfile)
 use eos,     only:gamma
 use damping, only:idamp,damp
 use io,      only:fatal
 integer,          intent(in)    :: npart
 real,             intent(in)    :: xyzh(:,:)
 real,             intent(inout) :: vxyzu(:,:)
 character(len=*), intent(in)    :: setupfile
 integer :: i,ierr

 call read_setupfile(setupfile,ierr)
 if (ierr /= 0) call fatal('spin_up','cannot read '//trim(setupfile))
 call derive_shape()
 iseed = -abs(pert_seed) - 1
 do i=1,npart
    if (xyzh(4,i) <= 0.) cycle
    call set_velocity_and_energy(xyzh(1:3,i),vxyzu(:,i),gamma,spin=.true.,set_u=.false.)
 enddo
 call set_run_options()
 idamp = 0
 damp  = 0.
 print "(a,i10,a,a,a,es10.3,a,i3)",' spin_up: ',npart,' particles, seed ', &
       trim(pert_kind),' amp = ',pert_amp,' m = ',pert_m

end subroutine spin_up

subroutine derive_shape()

 polyn = 1./(gamma_torus - 1.)
 f_c   = 0.5 - C_prime
 r_in  = r_center/(1. + sqrt(1. - 2.*C_prime))
 r_out = r_center/(1. - sqrt(1. - 2.*C_prime))
 p_orb = 2.*pi*r_center**1.5

end subroutine derive_shape

! production run: fixed point mass, cadence of analysis and dumps
subroutine set_run_options()
 use options,         only:iexternalforce
 use externalforces,  only:mass1,accradius1,iext_star
 use timestep,        only:tmax,dtmax
 use io_control,      only:nout,nfulldump
 use shock_capturing, only:alphau

 iexternalforce = iext_star
 mass1      = 1.
 accradius1 = r_acc
 tmax       = n_orbits*p_orb
 dtmax      = p_orb/n_analysis_orbit
 nout       = n_analysis_orbit/n_dumps_orbit
 nfulldump  = 1
 ! Artificial conductivity sees the pressure jump between neighbours in a stratified
 ! equilibrium and carries heat from the core to the envelope: at alphau = 1 the
 ! relaxing torus contracted 32% in rho_max in 6 orbits while its envelope blew out.
 ! The torus is homentropic, so there is no contact for it to smooth.
 alphau     = 0.

end subroutine set_run_options

! f(R,z) of PP84 eq. (2.9) in units G M = 1: (n+1) p/rho = f/r_center
real function fshape(r,z)
 real, intent(in) :: r,z

 fshape = r_center/sqrt(r*r + z*z) - 0.5*(r_center/r)**2 - C_prime

end function fshape

real function rhohat(r,z)
 real, intent(in) :: r,z
 real :: f

 f = fshape(r,z)
 if (f > 0.) then
    rhohat = (f/f_c)**polyn
 else
    rhohat = 0.
 endif

end function rhohat

real function surface_height(r)
 real, intent(in) :: r
 real :: s

 s = 1./(0.5*r_center/r**2 + C_prime/r_center)
 surface_height = sqrt(max((r_center*s)**2 - r*r,0.))

end function surface_height

real function surface_height_max()
 integer, parameter :: n = 4000
 integer :: i
 real :: r

 surface_height_max = 0.
 do i=0,n
    r = r_in + (r_out - r_in)*i/real(n)
    surface_height_max = max(surface_height_max,surface_height(r))
 enddo

end function surface_height_max

real function torus_volume()
 integer, parameter :: nr = 800, nz = 400
 integer :: i,j
 real :: dr,dz,r,z

 dr = (r_out - r_in)/nr
 dz = z_top/nz
 torus_volume = 0.
 do i=1,nr
    r = r_in + (i - 0.5)*dr
    do j=1,nz
       z = (j - 0.5)*dz
       torus_volume = torus_volume + rhohat(r,z)*r
    enddo
 enddo
 torus_volume = 4.*pi*torus_volume*dr*dz

end function torus_volume

! x at which the monotone table y(x) reaches yt
real function invert(n,x,y,yt)
 integer, intent(in) :: n
 real,    intent(in) :: x(n),y(n),yt
 integer :: lo,hi,mid
 real :: w

 if (yt <= y(1)) then
    invert = x(1)
    return
 elseif (yt >= y(n)) then
    invert = x(n)
    return
 endif
 lo = 1
 hi = n
 do while (hi - lo > 1)
    mid = (lo + hi)/2
    if (y(mid) > yt) then
       hi = mid
    else
       lo = mid
    endif
 enddo
 w = (yt - y(lo))/max(y(hi) - y(lo),tiny(w))
 invert = x(lo) + w*(x(hi) - x(lo))

end function invert

! Rings of particles: columns in R by stretch mapping the mid-plane spacing,
! rings in z by stretch mapping the column mass. Every ring holds a multiple of
! nsym particles, so the whole set is nsym-fold symmetric and modes with m not a
! multiple of nsym start at round-off. Ring masses match the analytic column, so
! a ring that cannot hold an isotropic count is stretched in z, not misweighted.
subroutine set_rings(npart,xyzh,hfact)
 integer, intent(inout) :: npart
 real,    intent(inout) :: xyzh(:,:)
 real,    intent(in)    :: hfact
 integer, parameter :: nt = 20000, nzt = 4000
 real :: rt(nt),lam(nt),zt(nzt),mz(nzt)
 real :: dr,lam_tot,r_lo,r_hi,rc,drc,zs,dzt,m_half,mlo,mring,zlo,zk
 integer :: i,k,ncol,icol,nphi,ilayer

 dr = (r_out - r_in)/(nt - 1)
 rt(1)  = r_in
 lam(1) = 0.
 do i=2,nt
    rt(i)  = r_in + (i - 1)*dr
    lam(i) = lam(i-1) + 0.5*dr*(spacing_inv(rt(i-1),0.) + spacing_inv(rt(i),0.))
 enddo
 lam_tot = lam(nt)
 ncol = max(1,nint(lam_tot))

 do icol=1,ncol
    r_lo = invert(nt,rt,lam,lam_tot*(icol - 1)/ncol)
    r_hi = invert(nt,rt,lam,lam_tot*icol/ncol)
    rc   = 0.5*(r_lo + r_hi)
    drc  = r_hi - r_lo
    zs   = surface_height(rc)
    if (zs <= 0.) cycle

    ! half-column mass above z = 0, per ring of radius rc and width drc
    dzt = zs/(nzt - 1)
    zt(1) = 0.
    mz(1) = 0.
    do k=2,nzt
       zt(k) = (k - 1)*dzt
       mz(k) = mz(k-1) + 0.5*dzt*(rhohat(rc,zt(k-1)) + rhohat(rc,zt(k)))
    enddo
    mz = 2.*pi*rc*drc*rho_c*mz
    m_half = mz(nzt)
    if (2.*m_half < 0.5*nsym*pmass) cycle

    ! mid-plane ring, cell [-zlo, zlo]
    nphi = ring_count(rc,0.)
    if (nphi*pmass >= 2.*m_half) nphi = nsym*max(1,nint(2.*m_half/(nsym*pmass)))
    call add_ring(npart,xyzh,rc,0.,nphi,mod(icol,2),hfact)
    mlo = min(0.5*nphi*pmass,m_half)
    zlo = invert(nzt,zt,mz,mlo)

    ilayer = 0
    do while (m_half - mlo >= 0.5*nsym*pmass)
       ilayer = ilayer + 1
       zk = invert(nzt,zt,mz,mlo + 0.5*ring_count(rc,zlo)*pmass)
       nphi = ring_count(rc,zk)
       mring = nphi*pmass
       if (mlo + mring > m_half) then
          nphi  = nsym*max(1,nint((m_half - mlo)/(nsym*pmass)))
          mring = m_half - mlo
       endif
       zk = invert(nzt,zt,mz,mlo + 0.5*mring)
       call add_ring(npart,xyzh,rc, zk,nphi,mod(icol+ilayer,2),hfact)
       call add_ring(npart,xyzh,rc,-zk,nphi,mod(icol+ilayer,2),hfact)
       mlo = mlo + mring
       zlo = invert(nzt,zt,mz,mlo)
    enddo
 enddo

 if (npart == 0) stop 'setup: no rings placed, increase np'

end subroutine set_rings

! 1/d for the isotropic spacing d = (m/rho)^(1/3); zero outside the torus
real function spacing_inv(r,z)
 real, intent(in) :: r,z
 real :: rh

 rh = rhohat(r,z)
 spacing_inv = 0.
 if (rh > 0.) spacing_inv = (rho_c*rh/pmass)**(1./3.)

end function spacing_inv

! isotropic particle count of a ring, rounded to a multiple of nsym
integer function ring_count(r,z)
 real, intent(in) :: r,z
 real :: rh,d

 rh = max(rhohat(r,z),1.e-30)
 d  = (pmass/(rho_c*rh))**(1./3.)
 ring_count = nsym*max(1,nint(2.*pi*r/(nsym*d)))

end function ring_count

subroutine add_ring(npart,xyzh,r,z,nphi,ishift,hfact)
 integer, intent(inout) :: npart
 real,    intent(inout) :: xyzh(:,:)
 real,    intent(in)    :: r,z,hfact
 integer, intent(in)    :: nphi,ishift
 integer :: i
 real :: phi,dphi,h

 dphi = 2.*pi/nphi
 h = hfact*(pmass/(rho_c*max(rhohat(r,z),1.e-6)))**(1./3.)
 do i=1,nphi
    npart = npart + 1
    if (npart > size(xyzh,2)) return
    phi = (i - 1 + 0.5*ishift)*dphi
    xyzh(1,npart) = r*cos(phi)
    xyzh(2,npart) = r*sin(phi)
    xyzh(3,npart) = z
    xyzh(4,npart) = h
 enddo

end subroutine add_ring

! Monte Carlo positions: noise in every mode at the N^(-1/2) level
subroutine set_random(npart,xyzh,hfact)
 use random, only:ran2
 integer, intent(inout) :: npart
 real,    intent(inout) :: xyzh(:,:)
 real,    intent(in)    :: hfact
 real :: r,z,phi

 iseed = -abs(pert_seed)
 do while (npart < np)
    r = r_in + (r_out - r_in)*ran2(iseed)
    z = z_top*(2.*ran2(iseed) - 1.)
    if (ran2(iseed)*r_out > rhohat(r,z)*r) cycle
    phi = 2.*pi*ran2(iseed)
    npart = npart + 1
    if (npart > size(xyzh,2)) return
    xyzh(1,npart) = r*cos(phi)
    xyzh(2,npart) = r*sin(phi)
    xyzh(3,npart) = z
    xyzh(4,npart) = hfact*(pmass/(rho_c*max(rhohat(r,z),1.e-6)))**(1./3.)
 enddo

end subroutine set_random

! l = sqrt(G M r_center) = const, u = n p/rho, and the v_R seed of athena-vv:
! pert_amp * c_s * (f/f_c) * shape(phi). shape is what acc_disk_visc_vv.cpp
! computes, peak 1 and rms 1/sqrt(2), not the unit rms its comment claims.
subroutine set_velocity_and_energy(xyz,vxyzu,gamma,spin,set_u)
 use random, only:ran2
 real,    intent(in)    :: xyz(3),gamma
 real,    intent(inout) :: vxyzu(:)
 logical, intent(in)    :: spin,set_u
 real, save :: phase(64)
 logical, save :: first = .true.
 real :: r,phi,f,vphi,dvr,cs,shape
 integer :: k

 if (first) then
    do k=1,size(phase)
       phase(k) = pi*(2.*ran2(iseed) - 1.)
    enddo
    first = .false.
 endif

 r   = sqrt(xyz(1)**2 + xyz(2)**2)
 phi = atan2(xyz(2),xyz(1))
 f   = max(fshape(r,xyz(3)),0.)
 vphi = sqrt(r_center)/r
 cs   = sqrt(gamma*f/((polyn + 1.)*r_center))

 select case(trim(pert_kind))
 case('single')
    shape = sin(pert_m*phi)
 case('multi')
    shape = 0.
    do k=1,min(pert_m,size(phase))
       shape = shape + sin(k*phi + phase(k))
    enddo
    shape = shape/sqrt(real(min(pert_m,size(phase))))
 case('noise')
    shape = 2.*ran2(iseed) - 1.     ! uniform in [-1, 1), as HashDeviate
 case default
    shape = 0.
 end select
 dvr = pert_amp*cs*(f/f_c)*shape

 if (spin) then
    vxyzu(1) = -vphi*sin(phi) + dvr*cos(phi)
    vxyzu(2) =  vphi*cos(phi) + dvr*sin(phi)
 else
    vxyzu(1:2) = 0.
 endif
 vxyzu(3) = 0.
 if (set_u .and. size(vxyzu) >= 4) vxyzu(4) = polyn*f/((polyn + 1.)*r_center)

end subroutine set_velocity_and_energy

subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for the Papaloizou-Pringle torus setup (ppitorus)'
 write(iunit,"(/,a)") '# units'
 call write_inopt(M_bh,'M_bh','central mass [M_sun]',iunit)
 call write_inopt(chi,'chi','length unit L_0 = chi r_g',iunit)
 call write_inopt(rho_0,'rho_0','density at the maximum [g/cm^3]',iunit)
 write(iunit,"(/,a)") '# torus'
 call write_inopt(r_center,'r_center','radius of the density maximum [L_0]',iunit)
 call write_inopt(C_prime,'C_prime','surface constant (0 < C_prime < 0.5)',iunit)
 call write_inopt(gamma_torus,'gamma','adiabatic index',iunit)
 call write_inopt(r_acc,'r_acc','accretion radius of the central mass [L_0]',iunit)
 write(iunit,"(/,a)") '# particles'
 call write_inopt(np,'np','target number of particles',iunit)
 call write_inopt(lattice,'lattice','rings | random',iunit)
 call write_inopt(nsym,'nsym','ring counts are multiples of nsym (1 = no symmetry)',iunit)
 write(iunit,"(/,a)") '# seed perturbation'
 call write_inopt(pert_kind,'pert_kind','none | single | multi | noise',iunit)
 call write_inopt(pert_amp,'pert_amp','delta v_R in units of the local c_s',iunit)
 call write_inopt(pert_m,'pert_m','azimuthal mode (single) or highest mode (multi)',iunit)
 call write_inopt(pert_seed,'pert_seed','random seed (multi, noise, random)',iunit)
 write(iunit,"(/,a)") '# relaxation (stage 1, at rest in the effective potential)'
 call write_inopt(relax,'relax','relax before spinning up (needs r_center = 1)',iunit)
 call write_inopt(n_relax,'n_relax','relaxation time [orbits at r_center]',iunit)
 call write_inopt(damp_relax,'damp_relax','velocity damping rate [code units]',iunit)
 write(iunit,"(/,a)") '# run'
 call write_inopt(n_orbits,'n_orbits','run length [orbits at r_center]',iunit)
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
 call read_inopt(M_bh,'M_bh',db,min=0.,errcount=nerr)
 call read_inopt(chi,'chi',db,min=0.,errcount=nerr)
 call read_inopt(rho_0,'rho_0',db,min=0.,errcount=nerr)
 call read_inopt(r_center,'r_center',db,min=0.,errcount=nerr)
 call read_inopt(C_prime,'C_prime',db,min=0.,max=0.5,errcount=nerr)
 call read_inopt(gamma_torus,'gamma',db,min=1.,errcount=nerr)
 call read_inopt(r_acc,'r_acc',db,min=0.,errcount=nerr)
 call read_inopt(np,'np',db,min=1,errcount=nerr)
 call read_inopt(lattice,'lattice',db,errcount=nerr)
 call read_inopt(nsym,'nsym',db,min=1,errcount=nerr)
 call read_inopt(pert_kind,'pert_kind',db,errcount=nerr)
 call read_inopt(pert_amp,'pert_amp',db,min=0.,errcount=nerr)
 call read_inopt(pert_m,'pert_m',db,min=0,errcount=nerr)
 call read_inopt(pert_seed,'pert_seed',db,errcount=nerr)
 call read_inopt(relax,'relax',db,errcount=nerr)
 call read_inopt(n_relax,'n_relax',db,min=0.,errcount=nerr)
 call read_inopt(damp_relax,'damp_relax',db,min=0.,errcount=nerr)
 call read_inopt(n_orbits,'n_orbits',db,min=0.,errcount=nerr)
 call read_inopt(n_analysis_orbit,'n_analysis_orbit',db,min=1,errcount=nerr)
 call read_inopt(n_dumps_orbit,'n_dumps_orbit',db,min=1,errcount=nerr)
 call close_db(db)
 if (mod(n_analysis_orbit,n_dumps_orbit) /= 0) then
    call error('setup','n_dumps_orbit must divide n_analysis_orbit')
    nerr = nerr + 1
 endif
 if (r_acc >= r_center/(1. + sqrt(1. - 2.*C_prime))) then
    call error('setup','r_acc must lie inside the torus inner edge')
    nerr = nerr + 1
 endif
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

end module setup
