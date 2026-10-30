!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module inject
!
! Inflow of gas clouds (star particles, as in cloudcollisions) onto a torus
! around sink 1: clouds enter on a ring of radius feed_r (Gaussian width
! feed_dr) at a rate feed_mdot, on orbits tilted by up to feed_inc about a
! random line of nodes, with feed_frot times the circular speed and an
! isotropic dispersion feed_s v_K. Clouds come in point-symmetric pairs
! (x, v) and (-x, -v): both carry the same angular momentum and the pair adds
! no m = 1 harmonic and no centre-of-mass shift, so a lopsided mode cannot be
! seeded by the feeding itself. Stands for the bar-driven inflow from the
! circumnuclear disc that replaces gas lost to the black hole and spreading.
!
! :References: Hopkins & Quataert (2010), MNRAS 407, 1529
!
! :Owner: etrevol
!
! :Runtime parameters:
!   - feed_dr   : *Gaussian width of the injection ring [R_tor]*
!   - feed_frot : *rotation speed of injected clouds / v_K*
!   - feed_inc  : *maximum tilt of the injected orbits [deg]*
!   - feed_mdot : *injected mass per orbit at R_tor (2 pi code time) [M_BH]*
!   - feed_r    : *radius of the injection ring [R_tor]*
!   - feed_s    : *isotropic velocity dispersion of injected clouds / v_K*
!
! :Dependencies: infile_utils, io, part, partinject, random
!
 implicit none
 character(len=*), parameter, public :: inject_type = 'cloudfeed'

 public :: init_inject,inject_particles,write_options_inject,read_options_inject,&
      set_default_options_inject,update_injected_par

 real, public :: feed_mdot = 0.     ! 0 = no inflow
 real, public :: feed_r    = 2.
 real, public :: feed_dr   = 0.
 real, public :: feed_frot = 1.
 real, public :: feed_s    = 0.
 real, public :: feed_inc  = 0.
 integer(kind=8), public :: nfeed = 0   ! clouds injected since the start of this run
 integer, private :: iseed = -4271

 private

contains

subroutine init_inject(ierr)
 use part, only:nptmass
 integer, intent(out) :: ierr

 ierr = 0
 if (feed_mdot > 0. .and. nptmass < 1) ierr = 1

end subroutine init_inject

subroutine set_default_options_inject(flag)
 integer, optional, intent(in) :: flag

end subroutine set_default_options_inject

!-----------------------------------------------------------------------
!+
!  inject the clouds due over the last step, dtlast
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                            npart,npart_old,npartoftype,dtinject)
 use part,       only:massoftype,istar,iphase,iamtype,isdead_or_accreted
 use partinject, only:add_or_update_particle
 use random,     only:ran2,gauss_random
 real,    intent(in)    :: time,dtlast
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: npart,npart_old
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject
 real, parameter :: pi = 4.*atan(1.)
 real :: x0(3),v0(3),mbh,dmin,d,hguess,npairs,r,phi,vk,inc,node,xi(3),vi(3),u0(3),u1(3)
 integer :: i,k,ni

 dtinject = huge(dtinject)
 if (feed_mdot <= 0. .or. dtlast <= 0. .or. massoftype(istar) <= 0.) return
 x0  = xyzmh_ptmass(1:3,1)
 v0  = vxyz_ptmass(1:3,1)
 mbh = xyzmh_ptmass(4,1)

 ! pairs due this step; the fraction of a pair is drawn
 npairs = feed_mdot*dtlast/(2.*pi)/(2.*massoftype(istar))
 ni = int(npairs)
 if (ran2(iseed) < npairs - ni) ni = ni + 1
 if (ni == 0) return

 ! h of the cloud nearest the ring, so the density iteration starts close
 hguess = 0.1
 dmin = huge(dmin)
 do i=1,npart
    if (iamtype(iphase(i)) /= istar .or. isdead_or_accreted(xyzh(4,i))) cycle
    d = abs(norm2(xyzh(1:3,i) - x0) - feed_r)
    if (d < dmin) then
       dmin = d
       hguess = xyzh(4,i)
    endif
 enddo

 do k=1,ni
    r   = abs(feed_r + feed_dr*gauss_random(iseed))
    phi = 2.*pi*ran2(iseed)
    vk  = sqrt(mbh/r)
    ! circular orbit in the plane z = 0, then tilted about a random line of nodes
    xi = r*(/cos(phi),sin(phi),0./)
    vi = feed_frot*vk*(/-sin(phi),cos(phi),0./)
    inc  = feed_inc*pi/180.*(2.*ran2(iseed) - 1.)
    node = 2.*pi*ran2(iseed)
    u0 = (/cos(node),sin(node),0./)
    u1 = (/-sin(node),cos(node),0./)
    xi = rotate(xi,u0,u1,inc)
    vi = rotate(vi,u0,u1,inc)
    vi = vi + feed_s*vk*(/gauss_random(iseed),gauss_random(iseed),gauss_random(iseed)/)
    call add_or_update_particle(istar,x0 + xi,v0 + vi,hguess,0.,npart+1,npart,npartoftype,xyzh,vxyzu)
    call add_or_update_particle(istar,x0 - xi,v0 - vi,hguess,0.,npart+1,npart,npartoftype,xyzh,vxyzu)
    nfeed = nfeed + 2
 enddo

end subroutine inject_particles

! rotation by angle a about the unit vector u0 (in the plane z = 0, u1 = z x u0)
pure function rotate(x,u0,u1,a) result(y)
 real, intent(in) :: x(3),u0(3),u1(3),a
 real :: y(3),p,q

 p = dot_product(x,u0)
 q = dot_product(x,u1)
 y = p*u0 + (q*cos(a) - x(3)*sin(a))*u1 + (q*sin(a) + x(3)*cos(a))*(/0.,0.,1./)

end function rotate

subroutine update_injected_par
end subroutine update_injected_par

subroutine write_options_inject(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(feed_mdot,'feed_mdot','injected mass per orbit at R_tor (2 pi code time) [M_BH]',iunit)
 call write_inopt(feed_r,'feed_r','radius of the injection ring [R_tor]',iunit)
 call write_inopt(feed_dr,'feed_dr','Gaussian width of the injection ring [R_tor]',iunit)
 call write_inopt(feed_frot,'feed_frot','rotation speed of injected clouds / v_K',iunit)
 call write_inopt(feed_s,'feed_s','isotropic velocity dispersion of injected clouds / v_K',iunit)
 call write_inopt(feed_inc,'feed_inc','maximum tilt of the injected orbits [deg]',iunit)

end subroutine write_options_inject

subroutine read_options_inject(db,nerr)
 use infile_utils, only:inopts,read_inopt
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_inopt(feed_mdot,'feed_mdot',db,errcount=nerr,min=0.)
 call read_inopt(feed_r,'feed_r',db,errcount=nerr,min=0.)
 call read_inopt(feed_dr,'feed_dr',db,errcount=nerr,min=0.)
 call read_inopt(feed_frot,'feed_frot',db,errcount=nerr,min=0.)
 call read_inopt(feed_s,'feed_s',db,errcount=nerr,min=0.)
 call read_inopt(feed_inc,'feed_inc',db,errcount=nerr,min=0.,max=180.)

end subroutine read_options_inject

end module inject
