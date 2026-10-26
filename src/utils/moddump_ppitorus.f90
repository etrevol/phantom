!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Stage 2 of the relaxed Papaloizou-Pringle torus: spins up a torus relaxed at
! rest in the effective potential (setup ppitorus, relax = T), adds the seed and
! switches to the point mass. Reads the parameters from ppitorus.setup in the
! working directory.
!
! :References:
!   - Papaloizou & Pringle (1984), MNRAS 208, 721
!
! :Owner: etrevol
!
! :Runtime parameters: None
!
! :Dependencies: setup
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use setup, only:spin_up
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)

 call spin_up(npart,xyzh,vxyzu,'ppitorus.setup')

end subroutine modify_dump

end module moddump
