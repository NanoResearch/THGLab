c      
c
c     Program to run a LJ fluid simulation
c
c

      program LJFluid
      integer a, b, c
      real*8 blahblah
      
      real*8, allocatable :: derivs(:,:)

      integer N_PARTICLES = 216
      real*8 mp = 18/1000
      real*8 dt = 1e-6
      integer nstep = 10000

c     constants
      real*8 kb = 1.3806488e-23
      real*8 av = 6.0221413e23

c     For lennard-jones potential

      real*8 e_lj = .16275*4184
      real*8 sigma_lj = 3.16435/10;

c     Box dimensions
      real*8 xbox = 2.2
      real*8 ybox = 2.2
      real*8 zbox = 2.2

c     arrays to hold position, velocity, and energy derivatives

      real*8, allocatable :: x(:)
      real*8, allocatable :: y(:)
      real*8, allocatable :: z(:)

      real*8, allocatable :: vx(:)
      real*8, allocatable :: vy(:)
      real*8, allocatable :: vz(:)

      real*8, allocatable :: derivs(:,:)
      real*8, allocatable :: derivs_old(:,:)

      integer T
      T = nint(N_PARTICLES**(1.0/3.0))

      do i = 1, T
         print*, i
      end do
      end
