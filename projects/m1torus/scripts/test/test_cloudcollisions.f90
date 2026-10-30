! Tests of the cloud-collision module against analytic results:
! a head-on and an oblique pair (velocities, momentum, angular momentum, energy),
! pairs that must not collide, and the collision rate of elastic clouds in a box
! against kinetic theory, n sigma <v_rel> N / 2.
! Build and run: projects/m1torus/scripts/test/run_test.sh
program test_cloudcollisions
 use part,            only:istar,isetphase,maxtypes,ediss_coll,eheat_coll,esf_coll,nptmass,&
                            xyzmh_ptmass,vxyz_ptmass,ispinx,ispinz,nsinkproperties
 use cloudcollisions, only:r_cloud,e_rest,C_coll,ncoll,collide_clouds,heat_s,heat_tau,heat_phi,heat_clouds,v_alfven,&
                            sf_n,form_clump_sinks
 implicit none
 real, parameter :: pi = 4.*atan(1.)
 integer :: nfail
 nfail = 0
 call test_pairs(nfail)
 call test_box(nfail)
 call test_heat(nfail,.false.)
 call test_heat(nfail,.true.)
 call test_sinks(nfail)
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

 print "(a)",'magnetised head-on, v_n = 2: M_A = 10 keeps 1/(sqrt(2) 10 - 1/2), M_A = 0.5 is elastic'
 v_alfven = 0.2
 call two((/-0.009,0.,0./),(/1.,0.,0./),(/0.009,0.,0./),(/-1.,0.,0./),xyzh,vxyzu,iphase,mtype)
 call collide_clouds(2,xyzh,vxyzu,iphase,mtype,dt)
 call check('v1x, M_A = 10',vxyzu(1,1),-sqrt(1./(sqrt(2.)*10. - 0.5)),1.e-12,nfail)
 v_alfven = 4.
 call two((/-0.009,0.,0./),(/1.,0.,0./),(/0.009,0.,0./),(/-1.,0.,0./),xyzh,vxyzu,iphase,mtype)
 call collide_clouds(2,xyzh,vxyzu,iphase,mtype,dt)
 call check('v1x, M_A = 0.5 (elastic)',vxyzu(1,1),-1.,1.e-12,nfail)
 v_alfven = 0.

 print "(a)",'oblique, e = 0.3: p, L and E + E_lost conserved'
 e_rest = 0.3
 ediss_coll = 0.
 call two((/0.3,0.2,0.1/),(/0.7,-1.1,0.4/),(/0.31,0.205,0.112/),(/-0.5,0.9,-0.2/),xyzh,vxyzu,iphase,mtype)
 call budget(xyzh,vxyzu,p0,l0,e0)
 ncoll = 0
 call collide_clouds(2,xyzh,vxyzu,iphase,mtype,dt)
 call budget(xyzh,vxyzu,p1,l1,e1)
 call check('ncoll this call',real(ncoll),1.,0.,nfail)
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

! thermostat: a cold Keplerian ring about a unit mass; one call with a vanishing
! relaxation time must give sigma_R = sigma_z = heat_s v_K in each annulus and leave
! sigma_phi cold (or heat it too with heat_phi), keep momentum and L_z, and book the
! kinetic energy it adds
subroutine test_heat(nfail,phi)
 integer, intent(inout) :: nfail
 logical, intent(in)    :: phi
 integer, parameter :: n = 40000
 real, allocatable :: xyzh(:,:),vxyzu(:,:)
 integer(kind=1), allocatable :: iphase(:)
 real :: mtype(maxtypes),u(2),rr,ph,p0(3),p1(3),lz0,lz1,e0,e1,bh(4),vbh(3)
 real :: sel_r(3),s2(3),m1(3),cnt,cp,sp,dv(3),cyl(3),vk
 integer :: i

 allocate(xyzh(4,n),vxyzu(4,n),iphase(n))
 do i=1,n
    call random_number(u)
    rr = 0.9 + 0.2*u(1)
    ph = 2.*pi*u(2)
    xyzh(:,i)  = (/rr*cos(ph), rr*sin(ph), 0., 0.01/)
    vxyzu(:,i) = (/-sin(ph)/sqrt(rr), cos(ph)/sqrt(rr), 0., 0./)
 enddo
 iphase = isetphase(istar,iactive=.true.)
 mtype = 0.
 mtype(istar) = 1./n
 bh  = (/0., 0., 0., 1./)
 vbh = 0.
 heat_s = 0.1
 heat_tau = 1.e-8
 heat_phi = phi
 eheat_coll = 0.
 p0 = sum(vxyzu(1:3,:),dim=2)
 lz0 = sum(xyzh(1,:)*vxyzu(2,:) - xyzh(2,:)*vxyzu(1,:))
 e0 = 0.5*sum(vxyzu(1:3,:)**2)/n
 call heat_clouds(n,xyzh,vxyzu,iphase,mtype,bh,vbh,1.)
 p1 = sum(vxyzu(1:3,:),dim=2)
 lz1 = sum(xyzh(1,:)*vxyzu(2,:) - xyzh(2,:)*vxyzu(1,:))
 e1 = 0.5*sum(vxyzu(1:3,:)**2)/n
 ! dispersions in the annulus 0.97 < R < 1.03, in units of v_K
 s2 = 0.
 m1 = 0.
 cnt = 0.
 do i=1,n
    rr = sqrt(xyzh(1,i)**2 + xyzh(2,i)**2)
    if (abs(rr - 1.) > 0.03) cycle
    cp = xyzh(1,i)/rr
    sp = xyzh(2,i)/rr
    vk = 1./sqrt(rr)
    dv = vxyzu(1:3,i)
    cyl = (/ (dv(1)*cp + dv(2)*sp)/vk, (-dv(1)*sp + dv(2)*cp)/vk - 1., dv(3)/vk /)
    m1 = m1 + cyl
    s2 = s2 + cyl**2
    cnt = cnt + 1.
 enddo
 sel_r = sqrt(s2/cnt - (m1/cnt)**2)
 print "(a,l2,i7,a,3f8.4)",'thermostat heat_phi',phi,int(cnt),' clouds near R = 1, sigma/v_K (R, phi, z) =',sel_r
 call check('sigma_R / v_K (5%)',sel_r(1),0.1,0.05,nfail)
 if (phi) then
    call check('sigma_phi / v_K (5%)',sel_r(2),0.1,0.05,nfail)
 else
    call check('sigma_phi / v_K stays cold',sel_r(2),0.,0.01,nfail)
 endif
 call check('sigma_z / v_K (5%)',sel_r(3),0.1,0.05,nfail)
 call check('momentum x',p1(1),p0(1),1.e-10,nfail)
 call check('momentum y',p1(2),p0(2),1.e-10,nfail)
 call check('L_z',lz1,lz0,1.e-10,nfail)
 call check('E added = eheat',e1 - e0,eheat_coll,1.e-10,nfail)
 heat_s = 0.
 heat_phi = .false.

end subroutine test_heat

! clump sinks: a cold, bound clump of 40 overlapping clouds becomes one sink with its
! mass, momentum and angular momentum (orbital plus spin); a hot clump does not
subroutine test_sinks(nfail)
 integer, intent(inout) :: nfail
 integer, parameter :: n = 80
 real :: xyzh(4,n),vxyzu(4,n),mtype(maxtypes),u(6),dt,p0(3),l0(3),p1(3),l1(3),m0,m1s
 integer(kind=1) :: iphase(n)
 integer :: i

 allocate(xyzmh_ptmass(nsinkproperties,10),vxyz_ptmass(3,10))
 xyzmh_ptmass = 0.
 vxyz_ptmass = 0.
 nptmass = 0
 r_cloud = 0.01
 e_rest = 0.5
 sf_n = 16
 mtype = 0.
 mtype(istar) = 1.e-3
 iphase = isetphase(istar,iactive=.true.)
 do i=1,n
    call random_number(u)
    if (i <= 40) then       ! cold clump near (1,0,0), slow internal motion
       xyzh(1:3,i) = (/1.,0.,0./) + 0.006*(u(1:3) - 0.5)
       vxyzu(1:3,i) = (/0.,1.,0./) + 1.e-3*(u(4:6) - 0.5)
    else                    ! hot clump near (-1,0,0): internal speeds far above escape
       xyzh(1:3,i) = (/-1.,0.,0./) + 0.006*(u(1:3) - 0.5)
       vxyzu(1:3,i) = (/0.,-1.,0./) + 10.*(u(4:6) - 0.5)
    endif
 enddo
 xyzh(4,:) = 0.014
 vxyzu(4,:) = 0.
 call budget_m(xyzh,vxyzu,mtype(istar),p0,l0,m0)
 ! collisions fill the overlap counts the sink finder uses (they also dissipate a bit)
 e_rest = 1.
 call collide_clouds(n,xyzh,vxyzu,iphase,mtype,dt)
 call budget_m(xyzh,vxyzu,mtype(istar),p0,l0,m0)
 call form_clump_sinks(n,xyzh,vxyzu,iphase,mtype,0.)
 call check('one sink formed',real(nptmass),1.,0.,nfail)
 call check('sink mass = 40 clouds',xyzmh_ptmass(4,1),40.e-3,1.e-12,nfail)
 call budget_m(xyzh,vxyzu,mtype(istar),p1,l1,m1s)
 p1 = p1 + xyzmh_ptmass(4,1)*vxyz_ptmass(1:3,1)
 l1 = l1 + xyzmh_ptmass(4,1)*(/xyzmh_ptmass(2,1)*vxyz_ptmass(3,1) - xyzmh_ptmass(3,1)*vxyz_ptmass(2,1), &
                              xyzmh_ptmass(3,1)*vxyz_ptmass(1,1) - xyzmh_ptmass(1,1)*vxyz_ptmass(3,1), &
                              xyzmh_ptmass(1,1)*vxyz_ptmass(2,1) - xyzmh_ptmass(2,1)*vxyz_ptmass(1,1)/) &
         + xyzmh_ptmass(ispinx:ispinz,1)
 m1s = m1s + xyzmh_ptmass(4,1)
 call check('mass',m1s,m0,1.e-14,nfail)
 call check('momentum y',p1(2),p0(2),1.e-12,nfail)
 call check('L_z (orbital + spin)',l1(3),l0(3),1.e-12,nfail)
 call check('E booked > 0 (bound group removed)',merge(1.,0.,esf_coll > 0.),1.,0.,nfail)
 sf_n = 0
 nptmass = 0

end subroutine test_sinks

subroutine budget_m(xyzh,vxyzu,m,p,l,mt)
 real, intent(in)  :: xyzh(:,:),vxyzu(:,:),m
 real, intent(out) :: p(3),l(3),mt
 integer :: i

 p = 0.
 l = 0.
 mt = 0.
 do i=1,size(xyzh,2)
    if (xyzh(4,i) <= 0.) cycle
    p = p + m*vxyzu(1:3,i)
    l(1) = l(1) + m*(xyzh(2,i)*vxyzu(3,i) - xyzh(3,i)*vxyzu(2,i))
    l(2) = l(2) + m*(xyzh(3,i)*vxyzu(1,i) - xyzh(1,i)*vxyzu(3,i))
    l(3) = l(3) + m*(xyzh(1,i)*vxyzu(2,i) - xyzh(2,i)*vxyzu(1,i))
    mt = mt + m
 enddo

end subroutine budget_m

end program test_cloudcollisions
