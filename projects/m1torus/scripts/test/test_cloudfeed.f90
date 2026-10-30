! Tests of the cloud inflow (inject_cloudfeed): the number of clouds per step,
! ring radius, circular speed, orbit tilts, point-symmetric pairs (centre of mass
! on the sink), the mean rate with fractional pairs, and the injected dispersion.
! Build and run: projects/m1torus/scripts/test/run_test.sh test_cloudfeed
! (needs a build with SRCINJECT=inject_cloudfeed.f90)
program test_cloudfeed
 use io,     only:nprocs
 use memory, only:allocate_memory
 use part,   only:istar,massoftype,npart,npartoftype,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass
 use inject, only:inject_particles,feed_mdot,feed_r,feed_dr,feed_frot,feed_s,feed_inc,inject_type
 implicit none
 real, parameter :: pi = 4.*atan(1.)
 integer :: nfail

 nfail = 0
 nprocs = 1   ! set by init_mpi in phantom
 call allocate_memory(int(4000,kind=8))
 call check('module is cloudfeed',merge(1.,0.,inject_type=='cloudfeed'),1.,0.,nfail)
 call test_ring(nfail)
 call test_rate(nfail)
 if (nfail == 0) then
    print "(a)",'ALL PASSED'
 else
    print "(a,i3,a)",'FAILED ',nfail,' checks'
    stop 1
 endif

contains

subroutine check(label,val,ref,tol,nfail)
 character(len=*), intent(in) :: label
 real,    intent(in)    :: val,ref,tol
 integer, intent(inout) :: nfail
 logical :: ok

 ok = abs(val - ref) <= tol*max(1.,abs(ref))
 if (.not.ok) nfail = nfail + 1
 print "(a,a40,2es16.8)",merge('  ok   ','  FAIL ',ok),label,val,ref

end subroutine check

subroutine reset(x0,v0)
 real, intent(in) :: x0(3),v0(3)

 npart = 0
 npartoftype = 0
 massoftype = 0.
 massoftype(istar) = 1.e-5
 nptmass = 1
 xyzmh_ptmass(:,1) = 0.
 xyzmh_ptmass(1:3,1) = x0
 xyzmh_ptmass(4,1) = 1.
 vxyz_ptmass(:,1) = 0.
 vxyz_ptmass(1:3,1) = v0

end subroutine reset

! 100 whole pairs on a thin tilted ring, no dispersion
subroutine test_ring(nfail)
 integer, intent(inout) :: nfail
 real :: x0(3),v0(3),dt,dtinj,r(3),v(3),l(3),errr,errv,errrad,cmin,tiltmax
 integer :: i,nold

 x0 = (/0.3,-0.2,0.1/)
 v0 = (/0.05,0.,0./)
 call reset(x0,v0)
 feed_mdot = 0.01
 feed_r    = 2.
 feed_dr   = 0.
 feed_frot = 0.8
 feed_s    = 0.
 feed_inc  = 30.
 ! npairs = feed_mdot dt / (2 pi) / (2 m) = 100
 dt = 2.*pi*200.*massoftype(istar)/feed_mdot
 nold = 0
 call inject_particles(0.,dt,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,nold,npartoftype,dtinj)
 call check('clouds injected',real(npart),200.,0.,nfail)
 call check('npartoftype(istar)',real(npartoftype(istar)),200.,0.,nfail)
 errr = 0.; errv = 0.; errrad = 0.; cmin = 1.; tiltmax = 0.
 do i=1,npart
    r = xyzh(1:3,i) - x0
    v = vxyzu(1:3,i) - v0
    errr   = max(errr,abs(norm2(r) - feed_r))
    errv   = max(errv,abs(norm2(v) - feed_frot*sqrt(1./feed_r)))
    errrad = max(errrad,abs(dot_product(r,v)))
    l = (/r(2)*v(3) - r(3)*v(2), r(3)*v(1) - r(1)*v(3), r(1)*v(2) - r(2)*v(1)/)
    cmin = min(cmin,l(3)/norm2(l))
    tiltmax = max(tiltmax,acos(l(3)/norm2(l))*180./pi)
 enddo
 call check('max | |r| - feed_r |',errr,0.,1.e-12,nfail)
 call check('max | |v| - frot v_K |',errv,0.,1.e-12,nfail)
 call check('max |r.v| (circular)',errrad,0.,1.e-12,nfail)
 call check('all tilts <= feed_inc',merge(1.,0.,cmin >= cos(feed_inc*pi/180.) - 1.e-12),1.,0.,nfail)
 call check('largest tilt near feed_inc (deg)',tiltmax,feed_inc,0.1,nfail)
 do i=1,3
    call check('centre of mass on the sink',sum(xyzh(i,1:npart))/npart,x0(i),1.e-12,nfail)
    call check('mean velocity of the sink',sum(vxyzu(i,1:npart))/npart,v0(i),1.e-12,nfail)
 enddo
 call check('no timestep limit',merge(1.,0.,dtinj > 1.e30),1.,0.,nfail)

end subroutine test_ring

! 0.37 pairs per call over 400 calls, with a wide ring and dispersion 0.2 v_K
subroutine test_rate(nfail)
 integer, intent(inout) :: nfail
 real :: x0(3),v0(3),dt,dtinj,r(3),v(3),vc(3),rr,s2,expect
 integer :: i,k,nold

 x0 = 0.
 v0 = 0.
 call reset(x0,v0)
 feed_mdot = 0.01
 feed_r    = 2.
 feed_dr   = 0.3
 feed_frot = 1.
 feed_s    = 0.2
 feed_inc  = 0.
 dt = 2.*pi*2.*0.37*massoftype(istar)/feed_mdot
 do k=1,400
    nold = npart
    call inject_particles(0.,dt,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,nold,npartoftype,dtinj)
 enddo
 expect = 2.*400*0.37
 ! pairs drawn with p = 0.37: sigma of the clouds = 2 sqrt(400 p (1 - p))
 call check('clouds over 400 calls, 3 sigma',real(npart),expect,3.*2.*sqrt(400*0.37*0.63)/expect,nfail)
 s2 = 0.
 do i=1,npart
    r = xyzh(1:3,i)
    v = vxyzu(1:3,i)
    rr = norm2(r)
    vc = sqrt(1./rr)*(/-r(2),r(1),0./)/rr
    s2 = s2 + sum((v - vc)**2)*rr/3.
 enddo
 call check('dispersion / v_K',sqrt(s2/npart),feed_s,0.1,nfail)
 call check('ring radius (mean |r|)',sum(norm2(xyzh(1:3,1:npart),dim=1))/npart,feed_r,0.05,nfail)

end subroutine test_rate

end program test_cloudfeed
