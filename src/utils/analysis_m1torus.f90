!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module analysis
!
! Diagnostics of the collisionless torus (setup m1torus), in the units of
! Bannikova et al. (2026, Sect. 4): azimuthal harmonics of the particle count
! in N_sec = 1200 sectors, A_m in particles per sector, so that A_1 = 40 at
! N = 128k reads as in their Table 2. Built into phantom with LIVE_ANALYSIS it
! runs every dtmax; phantomanalysis runs the same code on dumps.
!
! Appends one row per call to <prefix>_m1.ev: central mass and torus barycentre,
! A_m and phases for m = 1..5, and osculating elements about the central mass.
!
! :References:
!   - Bannikova et al. (2026), A&A, arXiv:2604.11528
!
! :Owner: etrevol
!
! :Runtime parameters: None
!
! :Dependencies: part
!
 implicit none
 character(len=20), parameter, public :: analysistype = 'm1torus'
 public :: do_analysis

 private

 integer, parameter :: mmax = 5
 real,    parameter :: nsec = 1200.
 real,    parameter :: pi = 4.*atan(1.)

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use part, only:nptmass,xyzmh_ptmass,vxyz_ptmass
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 character(len=120) :: prefix
 real(kind=8) :: ca(mmax),sa(mmax),amp(mmax),phase(mmax),xtb(3),xc(3),vc(3)
 real(kind=8) :: r(3),v(3),hv(3),ev(3),rr,v2,rv,a,e,inc,se,se2,si,sa_,nel,mu
 real(kind=8) :: phi,c1,s1,cm,sm,tmp
 integer :: i,m,n,iu
 logical :: exists

 prefix = dumpfile(1:index(dumpfile,'_',back=.true.)-1)
 if (len_trim(prefix)==0) prefix = 'm1torus'

 xc = 0.d0
 vc = 0.d0
 if (nptmass > 0) then
    xc = xyzmh_ptmass(1:3,1)
    vc = vxyz_ptmass(1:3,1)
 endif
 mu = 1.d0
 if (nptmass > 0) mu = xyzmh_ptmass(4,1)

 ca = 0.d0
 sa = 0.d0
 xtb = 0.d0
 se = 0.d0
 se2 = 0.d0
 si = 0.d0
 sa_ = 0.d0
 nel = 0.d0
 n = 0
 do i=1,npart
    if (xyzh(4,i) <= 0.) cycle
    n = n + 1
    xtb = xtb + xyzh(1:3,i)

    ! harmonics about the origin, the system barycentre (setup puts it there)
    phi = atan2(xyzh(2,i),xyzh(1,i))
    c1 = cos(phi)
    s1 = sin(phi)
    cm = c1
    sm = s1
    do m=1,mmax
       ca(m) = ca(m) + cm
       sa(m) = sa(m) + sm
       tmp = cm*c1 - sm*s1
       sm  = sm*c1 + cm*s1
       cm  = tmp
    enddo

    ! osculating elements about the central mass; unbound orbits are skipped
    r  = xyzh(1:3,i) - xc
    v  = vxyzu(1:3,i) - vc
    rr = sqrt(sum(r*r))
    v2 = sum(v*v)
    rv = sum(r*v)
    if (2.d0/rr - v2/mu <= 0.d0) cycle
    a  = 1.d0/(2.d0/rr - v2/mu)
    ev = ((v2 - mu/rr)*r - rv*v)/mu
    e  = sqrt(sum(ev*ev))
    hv = (/r(2)*v(3) - r(3)*v(2), r(3)*v(1) - r(1)*v(3), r(1)*v(2) - r(2)*v(1)/)
    inc = acos(max(-1.d0,min(1.d0,hv(3)/sqrt(sum(hv*hv)))))
    se  = se + e
    se2 = se2 + e*e
    si  = si + inc
    sa_ = sa_ + a
    nel = nel + 1.d0
 enddo
 if (n > 0) xtb = xtb/n

 ! least squares over equal sectors = discrete Fourier sums: a_m = (2/N_sec) sum cos
 do m=1,mmax
    amp(m)   = 2.d0/nsec*sqrt(ca(m)**2 + sa(m)**2)
    phase(m) = atan2(sa(m),ca(m))
 enddo
 if (nel > 0.d0) then
    se = se/nel
    se2 = sqrt(se2/(2.d0*nel))
    si = si/nel*180.d0/pi
    sa_ = sa_/nel
 endif

 inquire(file=trim(prefix)//'_m1.ev',exist=exists)
 open(newunit=iu,file=trim(prefix)//'_m1.ev',status='unknown',position='append', &
      form='formatted')
 if (.not.exists) then
    write(iu,"(a)",advance='no') '# time n xc yc zc xtb ytb ztb sigma0'
    do m=1,mmax
       write(iu,"(a,i0)",advance='no') ' A',m
    enddo
    do m=1,mmax
       write(iu,"(a,i0)",advance='no') ' Phi',m
    enddo
    write(iu,"(a)") ' e_mean e_sigma i_mean a_mean'
 endif
 write(iu,"(es24.16,i10,*(es24.16))") time,n,xc,xtb,n/nsec,amp,phase,se,se2,si,sa_
 close(iu)

end subroutine do_analysis

end module analysis
