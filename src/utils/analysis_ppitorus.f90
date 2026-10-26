!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module analysis
!
! Diagnostics of the Papaloizou-Pringle torus (setup ppitorus). Built into
! phantom with LIVE_ANALYSIS, it runs every dtmax whether or not a dump is
! written, on the double precision particle arrays; phantomanalysis runs the
! same code on dumps.
!
! Appends to <prefix>_ppi.ev one row per call: torus mass, L_z, centre of
! mass, azimuthal mode amplitudes A_m of the surface density and of the radial
! mass flux, m = 0..mmax, and the phase of m = 1..4 at r_center. Appends to
! <prefix>_ppimaps.bin one frame per call: SPH-interpolated mid-plane maps
! (rho, p, v_R, v_phi) on an (R, phi) grid, rho in the y = 0 plane, and the
! azimuthally averaged rho(R, z). Layout: scripts/ppi/vis/ppi_data.py.
!
! :References:
!   - Papaloizou & Pringle (1984), MNRAS 208, 721
!
! :Owner: etrevol
!
! :Runtime parameters: None
!
! :Dependencies: eos, infile_utils, kernel, part
!
 implicit none
 character(len=20), parameter, public :: analysistype = 'ppitorus'
 public :: do_analysis

 private

 integer, parameter :: mmax = 16, nbin = 64
 integer, parameter :: nr = 192, nphi = 256, nx = 384
 integer, parameter :: maps_magic = 1984208, maps_version = 1
 real, parameter :: pi = 4.*atan(1.)
 logical, save :: initialised = .false.
 integer, save :: nz
 real, save :: r_center = 1., C_prime = 0.2
 real, save :: r_in, r_out, z_top, r_lo, r_hi, x_half, z_half, cnorm2d

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 character(len=120) :: prefix
 integer :: iu

 prefix = dumpfile(1:index(dumpfile,'_',back=.true.)-1)
 if (len_trim(prefix)==0) prefix = 'ppitorus'
 if (.not.initialised) call init_geometry(prefix)

 call write_modes(trim(prefix)//'_ppi.ev',xyzh,vxyzu,particlemass,npart,time)

 open(newunit=iu,file=trim(prefix)//'_ppimaps.bin',access='stream',form='unformatted',&
      status='unknown',position='append')
 call write_maps(iu,xyzh,vxyzu,particlemass,npart,time)
 close(iu)

end subroutine do_analysis

! Torus geometry from <prefix>.setup; without one, r_center = 1 and C' = 0.2.
subroutine init_geometry(prefix)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in) :: prefix
 type(inopts), allocatable :: db(:)
 integer :: ierr,nerr,i
 real :: r,s

 nerr = 0
 call open_db_from_file(db,trim(prefix)//'.setup',27,ierr)
 if (ierr == 0) then
    call read_inopt(r_center,'r_center',db,errcount=nerr)
    call read_inopt(C_prime,'C_prime',db,errcount=nerr)
    call close_db(db)
 else
    print "(a)",' ppitorus analysis: no '//trim(prefix)//'.setup, assuming r_center = 1, C'' = 0.2'
 endif
 r_in  = r_center/(1. + sqrt(1. - 2.*C_prime))
 r_out = r_center/(1. - sqrt(1. - 2.*C_prime))
 z_top = 0.
 do i=0,4000
    r = r_in + (r_out - r_in)*i/4000.
    s = 1./(0.5*r_center/r**2 + C_prime/r_center)
    z_top = max(z_top,sqrt(max((r_center*s)**2 - r*r,0.)))
 enddo
 r_lo   = 0.5*r_in
 r_hi   = r_out + 0.25*(r_out - r_in)
 x_half = r_hi
 z_half = 1.25*z_top
 nz = max(16,nint(nx*z_half/x_half))
 cnorm2d = kernel_norm_2d()
 initialised = .true.

end subroutine init_geometry

! 1/int w(q) 2 pi q dq for the compiled kernel, so the ring projection below
! does not depend on which kernel phantom was built with
real function kernel_norm_2d()
 use kernel, only:wkern,radkern
 integer, parameter :: n = 4000
 integer :: i
 real :: q,s

 s = 0.
 do i=1,n
    q = (i - 0.5)*radkern/n
    s = s + wkern(q*q,q)*2.*pi*q
 enddo
 kernel_norm_2d = 1./(s*radkern/n)

end function kernel_norm_2d

! A_m = [sum_k |Sigma_m(R_k)|^2 R_k dR]^(1/2) with Sigma_m the m-th azimuthal
! coefficient of the surface density in annulus k; Fourier first, radius second,
! as vvlib/ppi.py does, so radial lobes of opposite sign cannot cancel. Particles
! are shared between the two nearest annuli (cloud in cell): with hard bins the
! copies of a ring cross an edge on different steps, and the incomplete ring
! puts a flat spectrum into every m.
subroutine write_modes(filename,xyzh,vxyzu,pmass,npart,time)
 character(len=*), intent(in) :: filename
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:),pmass,time
 integer,          intent(in) :: npart
 complex(kind=8) :: cs(0:mmax,nbin),cf(0:mmax,nbin),e1,em
 real(kind=8) :: amp(0:mmax),ampf(0:mmax),phase(4),w,rb,dr
 real(kind=8) :: mtot,lz,xcm,ycm,r,vr
 integer :: i,k,m,nalive,iu,kc,kb
 logical :: exists
 real(kind=8) :: s,wk(2)

 cs = (0.d0,0.d0)
 cf = (0.d0,0.d0)
 mtot = 0.d0
 lz = 0.d0
 xcm = 0.d0
 ycm = 0.d0
 nalive = 0
 dr = (r_hi - r_lo)/nbin
 do i=1,npart
    if (xyzh(4,i) <= 0.) cycle
    nalive = nalive + 1
    lz  = lz + xyzh(1,i)*vxyzu(2,i) - xyzh(2,i)*vxyzu(1,i)
    xcm = xcm + xyzh(1,i)
    ycm = ycm + xyzh(2,i)
    r = sqrt(xyzh(1,i)**2 + xyzh(2,i)**2)
    s = (r - r_lo)/dr + 0.5d0
    k = floor(s)
    if (k < 0 .or. k > nbin) cycle
    wk(2) = s - k
    wk(1) = 1.d0 - wk(2)
    vr = (xyzh(1,i)*vxyzu(1,i) + xyzh(2,i)*vxyzu(2,i))/r
    e1 = cmplx(xyzh(1,i)/r,-xyzh(2,i)/r,kind=8)
    em = (1.d0,0.d0)
    do m=0,mmax
       do kb=1,2
          if (k + kb - 1 < 1 .or. k + kb - 1 > nbin) cycle
          cs(m,k+kb-1) = cs(m,k+kb-1) + wk(kb)*em
          cf(m,k+kb-1) = cf(m,k+kb-1) + wk(kb)*vr*em
       enddo
       em = em*e1
    enddo
 enddo
 mtot = pmass*nalive
 lz   = pmass*lz
 if (nalive > 0) then
    xcm = xcm/nalive
    ycm = ycm/nalive
 endif

 amp  = 0.d0
 ampf = 0.d0
 do k=1,nbin
    rb = r_lo + (k - 0.5d0)*dr
    w  = pmass/(2.d0*pi*rb*dr)
    amp(:)  = amp(:)  + abs(w*cs(:,k))**2*rb*dr
    ampf(:) = ampf(:) + abs(w*cf(:,k))**2*rb*dr
 enddo
 amp  = sqrt(amp)
 ampf = sqrt(ampf)
 kc = min(nbin,max(1,int((r_center - r_lo)/dr) + 1))
 do m=1,4
    phase(m) = atan2(aimag(cs(m,kc)),real(cs(m,kc)))
 enddo

 inquire(file=filename,exist=exists)
 open(newunit=iu,file=filename,status='unknown',position='append',form='formatted')
 if (.not.exists) then
    write(iu,"(a)",advance='no') '# time mass nalive Lz xcm ycm'
    do m=0,mmax
       write(iu,"(a,i0)",advance='no') ' A',m
    enddo
    do m=0,mmax
       write(iu,"(a,i0)",advance='no') ' F',m
    enddo
    write(iu,"(a)") ' phi1 phi2 phi3 phi4'
 endif
 write(iu,"(es24.16,es24.16,i10,3es24.16,*(es24.16))") time,mtot,nalive,lz,xcm,ycm, &
       amp,ampf,phase
 close(iu)

end subroutine write_modes

! Scatter interpolation: rho = sum m W, other fields Shepard-normalised by
! sum (m/rho) W, so they stay finite at the surface where the sum is incomplete.
subroutine write_maps(iu,xyzh,vxyzu,pmass,npart,time)
 use kernel, only:wkern,radkern,radkern2,cnormk
 use part,   only:rhoh
 use eos,    only:gamma
 integer, intent(in) :: iu,npart
 real,    intent(in) :: xyzh(:,:),vxyzu(:,:),pmass,time
 real(kind=8) :: pol(nr,nphi,5),sl(nx,nz),av(nr,nz)
 real(kind=8) :: rg(nr),cphi(nphi),sphi(nphi),xg(nx),zg(nz)
 real :: dr,dphi,dx,dz,x,y,z,h,hi1,r,phi,rk,s,q2,q,wk,rhoi,vr,vphi,pri,dphimax
 integer :: i,ir,jp,j0,j1,jj,ix,iz,ir0,ir1,ix0,ix1,iz0,iz1
 real(kind=4) :: buf_pol(nr,nphi,4),buf_sl(nx,nz),buf_av(nr,nz)

 dr   = (r_hi - r_lo)/nr
 dphi = 2.*pi/nphi
 dx   = 2.*x_half/nx
 dz   = 2.*z_half/nz
 do ir=1,nr
    rg(ir) = r_lo + (ir - 0.5)*dr
 enddo
 do jp=1,nphi
    cphi(jp) = cos((jp - 0.5)*dphi)
    sphi(jp) = sin((jp - 0.5)*dphi)
 enddo
 do ix=1,nx
    xg(ix) = -x_half + (ix - 0.5)*dx
 enddo
 do iz=1,nz
    zg(iz) = -z_half + (iz - 0.5)*dz
 enddo

 pol = 0.d0
 sl  = 0.d0
 av  = 0.d0
!$omp parallel do default(none) schedule(guided) &
!$omp shared(npart,xyzh,vxyzu,pmass,gamma,rg,cphi,sphi,xg,zg,dr,dphi,dx,dz,r_lo,r_hi,x_half,z_half,nz) &
!$omp shared(cnorm2d) &
!$omp private(i,x,y,z,h,hi1,r,phi,rk,s,q2,q,wk,rhoi,vr,vphi,pri,dphimax) &
!$omp private(ir,jp,j0,j1,jj,ix,iz,ir0,ir1,ix0,ix1,iz0,iz1) &
!$omp reduction(+:pol,sl,av)
 do i=1,npart
    h = xyzh(4,i)
    if (h <= 0.) cycle
    x = xyzh(1,i)
    y = xyzh(2,i)
    z = xyzh(3,i)
    hi1  = 1./h
    rk   = radkern*h
    r    = sqrt(x*x + y*y)
    rhoi = rhoh(h,pmass)

    ! azimuthal average: the particle smeared round its ring, i.e. the 3D kernel
    ! integrated along phi, which is the 2D kernel in (R, z) while h << R
    ir0 = max(1, int((r - rk - r_lo)/dr) + 1)
    ir1 = min(nr,int((r + rk - r_lo)/dr) + 1)
    iz0 = max(1, int((z - rk + z_half)/dz) + 1)
    iz1 = min(nz,int((z + rk + z_half)/dz) + 1)
    do iz=iz0,iz1
       do ir=ir0,ir1
          q2 = ((rg(ir) - r)**2 + (zg(iz) - z)**2)*hi1*hi1
          if (q2 >= radkern2) cycle
          q  = sqrt(q2)
          av(ir,iz) = av(ir,iz) + pmass/(2.*pi*rg(ir))*cnorm2d*wkern(q2,q)*hi1*hi1
       enddo
    enddo

    ! mid-plane (R, phi) map
    if (abs(z) < rk) then
       phi  = atan2(y,x)
       vr   = (x*vxyzu(1,i) + y*vxyzu(2,i))/r
       vphi = (x*vxyzu(2,i) - y*vxyzu(1,i))/r
       pri  = (gamma - 1.)*rhoi*vxyzu(4,i)
       s    = sqrt(rk*rk - z*z)
       ir0  = max(1, int((r - s - r_lo)/dr) + 1)
       ir1  = min(nr,int((r + s - r_lo)/dr) + 1)
       do ir=ir0,ir1
          if (s >= rg(ir)) then
             j0 = 1
             j1 = nphi
          else
             dphimax = asin(s/rg(ir))
             j0 = floor((phi - dphimax)/dphi) + 1
             j1 = floor((phi + dphimax)/dphi) + 1
          endif
          do jj=j0,j1
             jp = modulo(jj - 1,nphi) + 1
             q2 = ((rg(ir)*cphi(jp) - x)**2 + (rg(ir)*sphi(jp) - y)**2 + z*z)*hi1*hi1
             if (q2 >= radkern2) cycle
             q  = sqrt(q2)
             wk = cnormk*wkern(q2,q)*hi1**3*pmass
             pol(ir,jp,1) = pol(ir,jp,1) + wk
             pol(ir,jp,2) = pol(ir,jp,2) + wk/rhoi*pri
             pol(ir,jp,3) = pol(ir,jp,3) + wk/rhoi*vr
             pol(ir,jp,4) = pol(ir,jp,4) + wk/rhoi*vphi
             pol(ir,jp,5) = pol(ir,jp,5) + wk/rhoi
          enddo
       enddo
    endif

    ! y = 0 slice (x, z)
    if (abs(y) < rk) then
       s   = sqrt(rk*rk - y*y)
       ix0 = max(1, int((x - s + x_half)/dx) + 1)
       ix1 = min(nx,int((x + s + x_half)/dx) + 1)
       iz0 = max(1, int((z - s + z_half)/dz) + 1)
       iz1 = min(nz,int((z + s + z_half)/dz) + 1)
       do iz=iz0,iz1
          do ix=ix0,ix1
             q2 = ((xg(ix) - x)**2 + y*y + (zg(iz) - z)**2)*hi1*hi1
             if (q2 >= radkern2) cycle
             q  = sqrt(q2)
             sl(ix,iz) = sl(ix,iz) + cnormk*wkern(q2,q)*hi1**3*pmass
          enddo
       enddo
    endif
 enddo
!$omp end parallel do

 buf_pol(:,:,1) = real(pol(:,:,1),kind=4)
 do jp=2,4
    where (pol(:,:,5) > 1.d-3)
       buf_pol(:,:,jp) = real(pol(:,:,jp)/pol(:,:,5),kind=4)
    elsewhere
       buf_pol(:,:,jp) = 0.
    end where
 enddo
 buf_sl = real(sl,kind=4)
 buf_av = real(av,kind=4)

 write(iu) maps_magic,maps_version,real(time,kind=8),nr,nphi,nx,nz, &
           real(r_lo,kind=8),real(r_hi,kind=8),real(x_half,kind=8),real(z_half,kind=8)
 write(iu) buf_pol,buf_sl,buf_av

end subroutine write_maps

end module analysis
