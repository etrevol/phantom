!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Numerical tests on an evolved m1torus dump. Reads moddump_m1torus.opts in the
! working directory (key = value):
!   mode  = fixed   replaces the central sink by a fixed point mass
!                   (iexternalforce = 1) at the origin, particles taken to the
!                   sink's frame, so the sink integration drops out of the problem
!   mode  = soften  sets the sink softening length to eps_c
!   eps_c = value   softening for mode = soften
!   mode  = gas     turns the torus particles into SPH gas, locally isothermal about
!                   the sink (ieos = 6) with c_s = hor v_K (q = 1/2, so H/R = hor),
!                   and sets the sink accretion radius to r_acc: the fluid limit of
!                   colliding clouds
!   hor   = value   H/R of the gas for mode = gas
!   r_acc = value   accretion radius for mode = gas or clouds
!   mode  = clouds  makes the torus particles colliding clouds of radius r_cloud
!                   (softened over it) with restitution e_rest, and sets the sink
!                   accretion radius to r_acc: collisions switched on in an evolved torus
!   r_cloud, e_rest = values for mode = clouds
!   heat_s, heat_tau = cloud thermostat for mode = clouds (heat_s = 0: off)
!   sf_n = overlapping clouds that make a bound clump a sink (0: off)
!
! :References: None
!
! :Owner: etrevol
!
! :Runtime parameters: None
!
! :Dependencies: cloudcollisions, eos, externalforces, infile_utils,
!   options, part
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,           only:nptmass,xyzmh_ptmass,vxyz_ptmass,ihsoft,ihacc,igas,istar,&
                           set_particle_type,iphase,iamtype
 use eos,            only:ieos,polyk,qfacdisc,isink
 use part,           only:hsoft_nbody,accrete_stars
 use cloudcollisions, only:r_cloud,e_rest,heat_s,heat_tau,sf_n
 use options,        only:iexternalforce
 use externalforces, only:mass1,accradius1,iext_star
 use infile_utils,   only:open_db_from_file,inopts,read_inopt,close_db
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 type(inopts), allocatable :: db(:)
 character(len=20) :: mode
 real :: eps_c,xc(3),vc(3),hor,r_acc
 integer :: ierr,i

 mode  = 'none'
 eps_c = 1.e-2
 hor   = 0.15
 r_acc = 0.05
 call open_db_from_file(db,'moddump_m1torus.opts',22,ierr)
 if (ierr /= 0) stop 'moddump_m1torus: no moddump_m1torus.opts'
 call read_inopt(mode,'mode',db,errcount=ierr)
 call read_inopt(eps_c,'eps_c',db,errcount=ierr)
 call read_inopt(hor,'hor',db,min=0.,default=0.15)
 call read_inopt(r_acc,'r_acc',db,min=0.,default=0.05)
 ! collisions only for mode = clouds; the other modes keep the dump's collisionless torus
 if (trim(mode)=='clouds') then
    call read_inopt(r_cloud,'r_cloud',db,min=0.,default=0.01)
    call read_inopt(e_rest,'e_rest',db,min=0.,max=1.,default=0.5)
    call read_inopt(heat_s,'heat_s',db,min=0.,default=0.)
    call read_inopt(heat_tau,'heat_tau',db,min=0.,default=1.)
    call read_inopt(sf_n,'sf_n',db,min=0,default=0)
 else
    r_cloud = 0.
 endif
 call close_db(db)

 select case(trim(mode))
 case('fixed')
    xc = xyzmh_ptmass(1:3,1)
    vc = vxyz_ptmass(1:3,1)
    do i=1,npart
       xyzh(1:3,i)  = xyzh(1:3,i) - xc
       vxyzu(1:3,i) = vxyzu(1:3,i) - vc
    enddo
    iexternalforce = iext_star
    mass1      = xyzmh_ptmass(4,1)
    accradius1 = xyzmh_ptmass(5,1)
    nptmass    = 0
    print "(a,es10.3)",' moddump_m1torus: sink -> fixed point mass, M = ',mass1
 case('soften')
    xyzmh_ptmass(ihsoft,1:nptmass) = eps_c
    print "(a,es10.3)",' moddump_m1torus: sink softening = ',eps_c
 case('gas')
    do i=1,npart
       if (iamtype(iphase(i))==istar) call set_particle_type(i,igas)
    enddo
    npartoftype(igas)  = npartoftype(igas) + npartoftype(istar)
    massoftype(igas)   = massoftype(istar)
    npartoftype(istar) = 0
    massoftype(istar)  = 0.
    ! c_s^2 = polyk r^(-2 q) about the sink; G = M_sink = 1, so c_s = hor v_K for q = 1/2
    ieos     = 6
    isink    = 1
    qfacdisc = 0.5
    polyk    = hor**2
    xyzmh_ptmass(ihacc,1) = r_acc
    print "(a,i8,a,f6.3,a,f6.3)",' moddump_m1torus: ',npartoftype(igas),' particles to gas, H/R = ', &
          hor,', r_acc = ',r_acc
 case('clouds')
    hsoft_nbody   = 1.4*r_cloud
    accrete_stars = .true.
    xyzmh_ptmass(ihacc,1) = r_acc
    print "(a,f6.3,a,f5.2,a,f6.3)",' moddump_m1torus: clouds of radius ',r_cloud,', e_rest = ', &
          e_rest,', r_acc = ',r_acc
 case default
    stop 'moddump_m1torus: mode must be fixed, soften, gas or clouds'
 end select

end subroutine modify_dump

end module moddump
