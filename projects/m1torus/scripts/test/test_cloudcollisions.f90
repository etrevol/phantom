! Tests of the cloud-collision module against analytic results:
! a head-on and an oblique pair (velocities, momentum, angular momentum, energy),
! pairs that must not collide, and the collision rate of elastic clouds in a box
! against kinetic theory, n sigma <v_rel> N / 2.
! Build and run: projects/m1torus/scripts/test/run_test_cloudcollisions.sh
program test_cloudcollisions
 use part,            only:istar,isetphase,maxtypes,ediss_coll
 use cloudcollisions, only:r_cloud,e_rest,C_coll,ncoll,collide_clouds
 implicit none
 real, parameter :: pi = 4.*atan(1.)
 integer :: nfail
 nfail = 0
 call test_pairs(nfail)
 call test_box(nfail)
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
 print "(a,a40,2es16.8,a)",merge('  ok   ','  FAIL ',ok),label,val,ref,''

end subroutine check

subroutine two(x1,v1,x2,v2,xyzh,vxyzu,iphase,mtype)
 real, intent(in)  :: x1(3),v1(3),x2(3),v2(3)
 real, intent(out) :: xyzh(4,2),vxyzu(4,2),mtype(maxtypes)
 integer(kind=1), intent(out) :: iphase(2)

 xyzh(1:3,1) = x1
 xyzh(1:3,2) = x2
 xyzh(4,:)   = 0.014
 vxyzu(1:3,1) = v1
 vxyzu(1:3,2) = v2
 vxyzu(4,:)   = 0.
 iphase = isetphase(istar,iactive=.true.)
 mtype  = 0.
 mtype(istar) = 1.

end subroutine two

subroutine test_pairs(nfail)
 integer, intent(inout) :: nfail
 real :: xyzh(4,2),vxyzu(4,2),mtype(maxtypes),dt,p0(3),l0(3),e0,p1(3),l1(3),e1
 integer(kind=1) :: iphase(2)

 r_cloud = 0.01
 print "(a)",'head-on, e = 0.5: v -> -e v, E lost = (1-e^2)/2 mu v_n^2 = 0.75'
 e_rest = 0.5
 ediss_coll = 0.
 call two((/-0.009,0.,0./),(/1.,0.,0./),(/0.009,0.,0./),(/-1.,0.,0./),xyzh,vxyzu,iphase,mtype)
 call collide_clouds(2,xyzh,vxyzu,iphase,mtype,dt)
 call check('v1x',vxyzu(1,1),-0.5,1.e-14,nfail)
 call check('v2x',vxyzu(1,2),0.5,1.e-14,nfail)
 call check('E lost',ediss_coll,0.75,1.e-14,nfail)

 print "(a)",'oblique, e = 0.3: p, L and E + E_lost conserved'
 e_rest = 0.3
 ediss_coll = 0.
 call two((/0.3,0.2,0.1/),(/0.7,-1.1,0.4/),(/0.31,0.205,0.112/),(/-0.5,0.9,-0.2/),xyzh,vxyzu,iphase,mtype)
 call budget(xyzh,vxyzu,p0,l0,e0)
 call collide_clouds(2,xyzh,vxyzu,iphase,mtype,dt)
 call budget(xyzh,vxyzu,p1,l1,e1)
 call check('ncoll this call',real(ncoll),2.,0.,nfail)
 call check('px',p1(1),p0(1),1.e-14,nfail)
 call check('py',p1(2),p0(2),1.e-14,nfail)
 call check('Lx',l1(1),l0(1),1.e-14,nfail)
 call check('Ly',l1(2),l0(2),1.e-14,nfail)
 call check('Lz',l1(3),l0(3),1.e-14,nfail)
 call check('E + E_lost',e1 + ediss_coll,e0,1.e-14,nfail)
 call check('E lost > 0',merge(1.,0.,ediss_coll > 0.),1.,0.,nfail)

 print "(a)",'receding and separated pairs: no change'
 ediss_coll = 0.
 call two((/-0.009,0.,0./),(/-1.,0.,0./),(/0.009,0.,0./),(/1.,0.,0./),xyzh,vxyzu,iphase,mtype)
 call collide_clouds(2,xyzh,vxyzu,iphase,mtype,dt)
 call check('receding v1x',vxyzu(1,1),-1.,0.,nfail)
 call two((/-0.011,0.,0./),(/1.,0.,0./),(/0.011,0.,0./),(/-1.,0.,0./),xyzh,vxyzu,iphase,mtype)
 call collide_clouds(2,xyzh,vxyzu,iphase,mtype,dt)
 call check('separated v1x',vxyzu(1,1),1.,0.,nfail)
 call check('E lost',ediss_coll,0.,0.,nfail)

end subroutine test_pairs

subroutine budget(xyzh,vxyzu,p,l,e)
 real, intent(in)  :: xyzh(:,:),vxyzu(:,:)
 real, intent(out) :: p(3),l(3),e
 integer :: i

 p = 0.
 l = 0.
 e = 0.
 do i=1,size(xyzh,2)
    p = p + vxyzu(1:3,i)
    l(1) = l(1) + xyzh(2,i)*vxyzu(3,i) - xyzh(3,i)*vxyzu(2,i)
    l(2) = l(2) + xyzh(3,i)*vxyzu(1,i) - xyzh(1,i)*vxyzu(3,i)
    l(3) = l(3) + xyzh(1,i)*vxyzu(2,i) - xyzh(2,i)*vxyzu(1,i)
    e = e + 0.5*dot_product(vxyzu(1:3,i),vxyzu(1:3,i))
 enddo

end subroutine budget

! N elastic clouds, uniform in a unit cube with specular walls, Maxwellian with
! 1D dispersion 1, moving ballistically (partners lost within 2r of a wall: <1%): collisions per unit time = N^2/(2V) pi (2r)^2 <v_rel>,
! <v_rel> = 4/sqrt(pi) for sigma_1D = 1
subroutine test_box(nfail)
 integer, intent(inout) :: nfail
 integer, parameter :: n = 20000
 real, allocatable :: xyzh(:,:),vxyzu(:,:)
 integer(kind=1), allocatable :: iphase(:)
 real :: mtype(maxtypes),t,tend,dt,dt0,expected,e0,e1,u(6)
 integer :: i,nstep
 integer(kind=8) :: nfirst

 allocate(xyzh(4,n),vxyzu(4,n),iphase(n))
 call random_seed()
 do i=1,n
    call random_number(u)
    xyzh(1:3,i) = u(1:3)
    ! Box-Muller for three Gaussian components
    vxyzu(1,i) = sqrt(-2.*log(1.-u(4)))*cos(2.*pi*u(5))
    vxyzu(2,i) = sqrt(-2.*log(1.-u(4)))*sin(2.*pi*u(5))
    call random_number(u(1:2))
    vxyzu(3,i) = sqrt(-2.*log(1.-u(1)))*cos(2.*pi*u(2))
 enddo
 xyzh(4,:)  = 0.003
 vxyzu(4,:) = 0.
 iphase = isetphase(istar,iactive=.true.)
 mtype = 0.
 mtype(istar) = 1.
 r_cloud = 0.002
 e_rest  = 1.
 C_coll  = 0.25
 ediss_coll = 0.
 e0 = 0.5*sum(vxyzu(1:3,:)**2)

 ! the first call resolves the pairs that overlap at t = 0; count from there
 ncoll = 0
 call collide_clouds(n,xyzh,vxyzu,iphase,mtype,dt0)
 nfirst = ncoll
 t = 0.
 tend = 0.05
 nstep = 0
 dt = dt0
 do while (t < tend)
    dt = min(dt,tend - t)
    xyzh(1:3,:) = xyzh(1:3,:) + dt*vxyzu(1:3,:)
    ! specular walls keep the density uniform (free streaming thins it at the faces)
    where (xyzh(1:3,:) < 0.)
       xyzh(1:3,:)  = -xyzh(1:3,:)
       vxyzu(1:3,:) = -vxyzu(1:3,:)
    endwhere
    where (xyzh(1:3,:) > 1.)
       xyzh(1:3,:)  = 2. - xyzh(1:3,:)
       vxyzu(1:3,:) = -vxyzu(1:3,:)
    endwhere
    call collide_clouds(n,xyzh,vxyzu,iphase,mtype,dt)
    t = t + dt
    nstep = nstep + 1
 enddo
 e1 = 0.5*sum(vxyzu(1:3,:)**2)
 expected = 0.5*real(n)**2*pi*(2.*r_cloud)**2*4./sqrt(pi)*tend
 print "(a,i6,a,i8,a,f9.1,a,i5)",'box: steps ',nstep,'  collisions ',ncoll-nfirst, &
       '  expected ',expected,'  first call ',nfirst
 call check('collisions / expected (3 sigma)',real(ncoll-nfirst)/expected,1.,3./sqrt(expected),nfail)
 call check('elastic: kinetic energy',e1,e0,1.e-12,nfail)

end subroutine test_box

end program test_cloudcollisions
