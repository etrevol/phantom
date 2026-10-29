!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module cloudcollisions
!
! Partially inelastic collisions between star particles treated as
! clouds of radius r_cloud ("sticky particles"). Two clouds collide when
! they overlap and approach; the impulse acts along the line of centres
! with restitution e_rest, so momentum and angular momentum are exact and
! the kinetic energy lost is added to ediss_coll (radiated).
!
! Optional heating (heat_s > 0): a thermostat that keeps the radial and vertical cloud
! velocity dispersions in each annulus at heat_s times the local Keplerian speed,
! standing for the turbulence that feedback and inflow drive in real tori; sigma_phi
! follows from the epicycles (sigma_phi/sigma_R = kappa/2 Omega = 0.5 for Keplerian), so
! it is not forced: holding it at heat_s too drives sigma_R to ~2 heat_s and unbinds
! the torus at heat_s = 0.3. heat_phi = .true. restores the old thermostat that also
! held sigma_phi (runs of 2026-10-29 before 21:00). It relaxes towards the target over
! heat_tau local orbits with zero-mean random kicks and books the energy in eheat_coll.
!
! :References: Vollmer et al. (2022), A&A 665, A102; Ostriker & Shetty (2011), ApJ 731, 41;
!              Binney & Tremaine (2008), Galactic Dynamics, Sect. 4.4.3
!
! :Owner: etrevol
!
! :Runtime parameters: None
!
! :Dependencies: dim, io, part, sortutils
!
 implicit none

 real,            public :: r_cloud = 0.    ! cloud radius, 0 = no collisions
 real,            public :: e_rest  = 0.5   ! normal restitution coefficient
 real,            public :: C_coll  = 0.25  ! dt < C_coll r_cloud / v_max
 integer(kind=8), public :: ncoll   = 0     ! collisions since the start of this run
 real,            public :: heat_s  = 0.    ! target sigma / v_K of the thermostat, 0 = off
 real,            public :: heat_tau = 1.   ! relaxation time of the thermostat [local orbits]
 logical,         public :: heat_phi = .false. ! also hold sigma_phi (old thermostat)
 integer,         public :: sf_n    = 0     ! clouds overlapping one cloud to form a sink, 0 = off
 integer, allocatable    :: nover(:)        ! clouds overlapping each particle, last step

 public :: collide_clouds,heat_clouds,form_clump_sinks

 private

contains

!-----------------------------------------------------------------------
!+
!  find overlapping, approaching pairs of clouds and apply the impulse;
!  returns the timestep that keeps clouds from passing through each other
!+
!-----------------------------------------------------------------------
subroutine collide_clouds(npart,xyzh,vxyzu,iphase,massoftype,dtcoll)
 use part,      only:istar,iamtype,isdead_or_accreted,ediss_coll
 use sortutils, only:indexx
 integer,         intent(in)    :: npart
 real,            intent(in)    :: xyzh(:,:),massoftype(:)
 real,            intent(inout) :: vxyzu(:,:)
 integer(kind=1), intent(in)    :: iphase(:)
 real,            intent(out)   :: dtcoll
 integer,         allocatable   :: list(:),idx(:)
 integer(kind=8), allocatable   :: key(:),skey(:)
 integer         :: n,i,j,k,kk,ix,iy,iz,dx,dy,dz,lo,hi
 integer(kind=8) :: kn
 real            :: dc,rsum2,xij(3),vij(3),r2,rij,vn,mcl,mu,dvi,vmax2

 dtcoll = huge(dtcoll)
 if (r_cloud <= 0.) return

 n = 0
 do i=1,npart
    if (iamtype(iphase(i))==istar .and. .not.isdead_or_accreted(xyzh(4,i))) n = n + 1
 enddo
 if (n < 2) return
 allocate(list(n),idx(n),key(n),skey(n))
 if (allocated(nover)) then
    if (size(nover) < npart) deallocate(nover)
 endif
 if (.not.allocated(nover)) allocate(nover(npart))
 nover = 0

 ! cells of size 2 r_cloud: a colliding partner lies in one of the 27 around
 dc = 2.*r_cloud
 n  = 0
 do i=1,npart
    if (iamtype(iphase(i))==istar .and. .not.isdead_or_accreted(xyzh(4,i))) then
       n = n + 1
       list(n) = i
       call cell_of(xyzh(1:3,i),dc,ix,iy,iz)
       key(n) = cell_key(ix,iy,iz)
    endif
 enddo
 call indexx(n,key,idx)
 skey = key(idx)

 ! all clouds have the same mass; m_i m_j / (m_i + m_j) = m / 2
 mcl   = massoftype(istar)
 mu    = 0.5*mcl
 rsum2 = (2.*r_cloud)**2

 ! sequential pass in cell order, each pair once (j > i); a cloud in several
 ! pairs feels them one after another within the step
 do k=1,n
    i = list(idx(k))
    call cell_of(xyzh(1:3,i),dc,ix,iy,iz)
    do dz=-1,1
       do dy=-1,1
          do dx=-1,1
             kn = cell_key(ix+dx,iy+dy,iz+dz)
             call key_range(n,skey,kn,lo,hi)
             do kk=lo,hi
                j = list(idx(kk))
                if (j <= i) cycle
                xij = xyzh(1:3,i) - xyzh(1:3,j)
                r2  = dot_product(xij,xij)
                if (r2 >= rsum2 .or. r2 <= 0.) cycle
                nover(i) = nover(i) + 1
                nover(j) = nover(j) + 1
                rij = sqrt(r2)
                vij = vxyzu(1:3,i) - vxyzu(1:3,j)
                vn  = dot_product(vij,xij)/rij
                if (vn >= 0.) cycle
                ! impulse J = (1 + e) mu |v_n| along the line of centres
                dvi = (1. + e_rest)*mu*(-vn)/mcl
                vxyzu(1:3,i) = vxyzu(1:3,i) + dvi*xij/rij
                vxyzu(1:3,j) = vxyzu(1:3,j) - dvi*xij/rij
                ediss_coll = ediss_coll + 0.5*(1. - e_rest**2)*mu*vn**2
                ncoll = ncoll + 1
             enddo
          enddo
       enddo
    enddo
 enddo

 ! relative speeds are at most 2 v_max; a pair closes 2 r_cloud in r_cloud/v_max
 vmax2 = 0.
 do k=1,n
    i = list(k)
    vmax2 = max(vmax2,dot_product(vxyzu(1:3,i),vxyzu(1:3,i)))
 enddo
 if (vmax2 > 0.) dtcoll = C_coll*r_cloud/sqrt(vmax2)

 deallocate(list,idx,key,skey)

end subroutine collide_clouds

!-----------------------------------------------------------------------
!+
!  integer cell coordinates of a position
!+
!-----------------------------------------------------------------------
pure subroutine cell_of(x,dc,ix,iy,iz)
 real,    intent(in)  :: x(3),dc
 integer, intent(out) :: ix,iy,iz

 ix = floor(x(1)/dc)
 iy = floor(x(2)/dc)
 iz = floor(x(3)/dc)

end subroutine cell_of

!-----------------------------------------------------------------------
!+
!  one integer key per cell, 21 bits per axis (2^21 cells of 2 r_cloud)
!+
!-----------------------------------------------------------------------
pure integer(kind=8) function cell_key(ix,iy,iz)
 integer, intent(in) :: ix,iy,iz
 integer(kind=8), parameter :: nb = 2097152_8, off = 1048576_8
 integer(kind=8) :: jx,jy,jz

 jx = min(max(int(ix,kind=8) + off,0_8),nb-1)
 jy = min(max(int(iy,kind=8) + off,0_8),nb-1)
 jz = min(max(int(iz,kind=8) + off,0_8),nb-1)
 cell_key = jx + nb*(jy + nb*jz)

end function cell_key

!-----------------------------------------------------------------------
!+
!  positions lo..hi of key kn in the sorted keys (empty if hi < lo)
!+
!-----------------------------------------------------------------------
pure subroutine key_range(n,skey,kn,lo,hi)
 integer,         intent(in)  :: n
 integer(kind=8), intent(in)  :: skey(n),kn
 integer,         intent(out) :: lo,hi
 integer :: a,b,m

 ! first position with skey >= kn
 a = 1
 b = n + 1
 do while (a < b)
    m = (a + b)/2
    if (skey(m) < kn) then
       a = m + 1
    else
       b = m
    endif
 enddo
 lo = a
 hi = lo - 1
 do while (hi < n)
    if (skey(hi+1) /= kn) exit
    hi = hi + 1
 enddo

end subroutine key_range

!-----------------------------------------------------------------------
!+
!  thermostat: in each annulus about the black hole, raise sigma_R and sigma_z (and
!  sigma_phi if heat_phi) that lie below heat_s v_K towards it over heat_tau local orbits, with
!  random kicks that change neither the momentum nor L_z of the annulus; the kinetic
!  energy added goes to eheat_coll
!+
!-----------------------------------------------------------------------
subroutine heat_clouds(npart,xyzh,vxyzu,iphase,massoftype,xyzmh_bh,vxyz_bh,dt)
 use part, only:istar,iamtype,isdead_or_accreted,eheat_coll
 integer,         intent(in)    :: npart
 real,            intent(in)    :: xyzh(:,:),massoftype(:),xyzmh_bh(:),vxyz_bh(:),dt
 real,            intent(inout) :: vxyzu(:,:)
 integer(kind=1), intent(in)    :: iphase(:)
 integer, parameter :: nb = 40, nmin = 16
 real,    parameter :: rlo = 0.02, rhi = 8., twopi = 2.*3.14159265358979323846
 logical, save :: seeded = .false.
 integer, allocatable :: ib(:)
 real,    allocatable :: kick(:,:),cyl(:,:)
 integer :: i,b,c,cnt(nb),nseed,it
 integer, allocatable :: seed(:)
 real :: dx(3),dv(3),rcyl,cph,sph,rb,vk2,frac,add(nb,3),sum1(nb,3),sum2(nb,3)
 real :: g1,g2,u1,u2,mcl,ek0,ek1,v(3),lz(nb),rsum(nb),kmean(nb,3)

 if (heat_s <= 0.) return
 if (.not.seeded) then
    call random_seed(size=nseed)
    allocate(seed(nseed))
    seed = 20261029 + 37*[(i,i=1,nseed)]
    call random_seed(put=seed)
    seeded = .true.
 endif
 allocate(ib(npart),kick(3,npart),cyl(3,npart))
 mcl = massoftype(istar)
 cnt = 0
 sum1 = 0.
 sum2 = 0.
 ib = 0
 do i=1,npart
    if (iamtype(iphase(i))/=istar .or. isdead_or_accreted(xyzh(4,i))) cycle
    dx = xyzh(1:3,i) - xyzmh_bh(1:3)
    dv = vxyzu(1:3,i) - vxyz_bh(1:3)
    rcyl = sqrt(dx(1)**2 + dx(2)**2)
    if (rcyl < rlo .or. rcyl >= rhi) cycle
    b = 1 + int(nb*log(rcyl/rlo)/log(rhi/rlo))
    ib(i) = b
    cph = dx(1)/rcyl
    sph = dx(2)/rcyl
    cyl(:,i) = (/ dv(1)*cph + dv(2)*sph, -dv(1)*sph + dv(2)*cph, dv(3) /)
    cnt(b) = cnt(b) + 1
    sum1(b,:) = sum1(b,:) + cyl(:,i)
    sum2(b,:) = sum2(b,:) + cyl(:,i)**2
 enddo
 ! variance to add per component in each annulus
 add = 0.
 do b=1,nb
    if (cnt(b) < nmin) cycle
    rb  = rlo*(rhi/rlo)**((b - 0.5)/nb)
    vk2 = xyzmh_bh(4)/rb
    frac = min(1., dt/(heat_tau*twopi*sqrt(rb**3/xyzmh_bh(4))))
    do c=1,3
       if (c == 2 .and. .not.heat_phi) cycle      ! sigma_phi is left to the epicycles
       add(b,c) = max(0., heat_s**2*vk2 - (sum2(b,c)/cnt(b) - (sum1(b,c)/cnt(b))**2))*frac
    enddo
 enddo
 ! Gaussian kicks in R, phi, z
 kick = 0.
 do i=1,npart
    b = ib(i)
    if (b == 0) cycle
    if (cnt(b) < nmin) cycle
    do c=1,3
       call random_number(u1)
       call random_number(u2)
       g1 = sqrt(-2.*log(max(1.-u1,1.e-300)))*cos(twopi*u2)
       kick(c,i) = sqrt(add(b,c))*g1
    enddo
 enddo
 ! to Cartesian kicks; then per annulus remove their mean (momentum) and their L_z (a
 ! uniform shift along e_phi), alternately until both vanish to round-off
 do i=1,npart
    b = ib(i)
    if (b == 0) cycle
    if (cnt(b) < nmin) cycle
    dx = xyzh(1:3,i) - xyzmh_bh(1:3)
    rcyl = sqrt(dx(1)**2 + dx(2)**2)
    cyl(:,i) = kick(:,i)
    kick(1,i) = (cyl(1,i)*dx(1) - cyl(2,i)*dx(2))/rcyl
    kick(2,i) = (cyl(1,i)*dx(2) + cyl(2,i)*dx(1))/rcyl
 enddo
 do it=1,4
    kmean = 0.
    do i=1,npart
       b = ib(i)
       if (b == 0) cycle
       if (cnt(b) < nmin) cycle
       kmean(b,:) = kmean(b,:) + kick(:,i)
    enddo
    lz = 0.
    rsum = 0.
    do i=1,npart
       b = ib(i)
       if (b == 0) cycle
       if (cnt(b) < nmin) cycle
       kick(:,i) = kick(:,i) - kmean(b,:)/cnt(b)
       dx = xyzh(1:3,i) - xyzmh_bh(1:3)
       lz(b) = lz(b) + dx(1)*kick(2,i) - dx(2)*kick(1,i)
       rsum(b) = rsum(b) + sqrt(dx(1)**2 + dx(2)**2)
    enddo
    do i=1,npart
       b = ib(i)
       if (b == 0) cycle
       if (cnt(b) < nmin) cycle
       dx = xyzh(1:3,i) - xyzmh_bh(1:3)
       rcyl = sqrt(dx(1)**2 + dx(2)**2)
       kick(1,i) = kick(1,i) + lz(b)/rsum(b)*dx(2)/rcyl
       kick(2,i) = kick(2,i) - lz(b)/rsum(b)*dx(1)/rcyl
    enddo
 enddo
 g2 = 0.
 do i=1,npart
    b = ib(i)
    if (b == 0) cycle
    if (cnt(b) < nmin) cycle
    v = vxyzu(1:3,i)
    ek0 = dot_product(v,v)
    v = v + kick(:,i)
    ek1 = dot_product(v,v)
    vxyzu(1:3,i) = v
    g2 = g2 + 0.5*mcl*(ek1 - ek0)
 enddo
 eheat_coll = eheat_coll + g2
 deallocate(ib,kick,cyl)

end subroutine heat_clouds

!-----------------------------------------------------------------------
!+
!  star formation in clumps: a cloud that sf_n or more clouds overlap, with the clouds
!  within 2 r_cloud of it gravitationally bound (kinetic energy about their centre of
!  mass plus their softened mutual potential < 0), becomes a sink at the centre of mass;
!  mass, momentum and angular momentum (orbital plus the sink's spin) are kept, the
!  group's internal energy goes to esf_coll; the sink accretes clouds within 2 r_cloud
!+
!-----------------------------------------------------------------------
subroutine form_clump_sinks(npart,xyzh,vxyzu,iphase,massoftype,time)
 use part, only:istar,iamtype,isdead_or_accreted,nptmass,xyzmh_ptmass,vxyz_ptmass, &
                ihacc,ihsoft,ispinx,ispiny,ispinz,itbirth,esf_coll
 use dim,  only:maxptmass
 use io,   only:iprint
 integer,         intent(in)    :: npart
 real,            intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real,            intent(in)    :: massoftype(:),time
 integer(kind=1), intent(in)    :: iphase(:)
 integer, allocatable :: cand(:),grp(:)
 integer :: i,j,k,l,nc,ng,nsink0
 real :: mcl,r2max,xcm(3),vcm(3),ek,ep,dx(3),spin(3),d2

 if (sf_n <= 0 .or. .not.allocated(nover)) return
 nc = count(nover(1:npart) >= sf_n)
 if (nc == 0) return
 mcl = massoftype(istar)
 r2max = (2.*r_cloud)**2
 allocate(cand(nc),grp(npart))
 nc = 0
 do i=1,npart
    if (nover(i) >= sf_n .and. .not.isdead_or_accreted(xyzh(4,i))) then
       nc = nc + 1
       cand(nc) = i
    endif
 enddo
 nsink0 = nptmass
 do k=1,nc
    i = cand(k)
    if (isdead_or_accreted(xyzh(4,i))) cycle      ! taken by an earlier sink this step
    ng = 0
    do j=1,npart
       if (iamtype(iphase(j))/=istar .or. isdead_or_accreted(xyzh(4,j))) cycle
       dx = xyzh(1:3,j) - xyzh(1:3,i)
       if (dot_product(dx,dx) < r2max) then
          ng = ng + 1
          grp(ng) = j
       endif
    enddo
    if (ng < sf_n + 1) cycle
    xcm = 0.
    vcm = 0.
    do l=1,ng
       xcm = xcm + xyzh(1:3,grp(l))
       vcm = vcm + vxyzu(1:3,grp(l))
    enddo
    xcm = xcm/ng
    vcm = vcm/ng
    ek = 0.
    ep = 0.
    spin = 0.
    do l=1,ng
       j = grp(l)
       dx = xyzh(1:3,j) - xcm
       ek = ek + 0.5*mcl*dot_product(vxyzu(1:3,j) - vcm,vxyzu(1:3,j) - vcm)
       spin = spin + mcl*cross(dx,vxyzu(1:3,j) - vcm)
       do i=1,l-1
          d2 = sum((xyzh(1:3,j) - xyzh(1:3,grp(i)))**2)
          ep = ep - mcl*mcl/sqrt(d2 + r_cloud**2)      ! Plummer-softened over r_cloud
       enddo
    enddo
    if (ek + ep >= 0.) cycle
    if (nptmass + 1 > maxptmass) then
       write(iprint,"(a)") ' form_clump_sinks: maxptmass reached, no more sinks'
       exit
    endif
    nptmass = nptmass + 1
    xyzmh_ptmass(:,nptmass)   = 0.
    xyzmh_ptmass(1:3,nptmass) = xcm
    xyzmh_ptmass(4,nptmass)   = ng*mcl
    xyzmh_ptmass(ihacc,nptmass)  = 2.*r_cloud
    xyzmh_ptmass(ihsoft,nptmass) = r_cloud
    xyzmh_ptmass(ispinx,nptmass) = spin(1)
    xyzmh_ptmass(ispiny,nptmass) = spin(2)
    xyzmh_ptmass(ispinz,nptmass) = spin(3)
    xyzmh_ptmass(itbirth,nptmass) = time
    vxyz_ptmass(1:3,nptmass)  = vcm
    do l=1,ng
       xyzh(4,grp(l)) = -abs(xyzh(4,grp(l)))
    enddo
    esf_coll = esf_coll - (ek + ep)
 enddo
 if (nptmass > nsink0) write(iprint,"(a,i5,a,i6)") ' form_clump_sinks: ',nptmass - nsink0, &
                                                   ' new sinks, total ',nptmass
 deallocate(cand,grp)

contains

 pure function cross(a,b)
  real, intent(in) :: a(3),b(3)
  real :: cross(3)
  cross = (/a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)/)
 end function cross

end subroutine form_clump_sinks

end module cloudcollisions
