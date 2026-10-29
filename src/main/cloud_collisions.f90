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
! :References: Vollmer et al. (2022), A&A 665, A102
!
! :Owner: etrevol
!
! :Runtime parameters: None
!
! :Dependencies: part, sortutils
!
 implicit none

 real,            public :: r_cloud = 0.    ! cloud radius, 0 = no collisions
 real,            public :: e_rest  = 0.5   ! normal restitution coefficient
 real,            public :: C_coll  = 0.25  ! dt < C_coll r_cloud / v_max
 integer(kind=8), public :: ncoll   = 0     ! collisions since the start of this run

 public :: collide_clouds

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

end module cloudcollisions
