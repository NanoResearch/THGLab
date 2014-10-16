c      
c
c     Program to run a LJ fluid simulation
c
c

      program LJFluid
      implicit none
      integer :: N_PARTICLES = 216
      real*8 :: mp = 18.0/1000
      real*8 :: dt = 1e-6
      integer :: nstep = 10000

c     constants
      real*8 :: kb = 1.3806488e-23
      real*8 :: av = 6.0221413e23

c     For lennard-jones potential

      real*8 :: e_lj = .16275*4184
      real*8 :: sigma_lj = 3.16435/10;

c     Box dimensions
      real*8 :: xbox = 2.2
      real*8 :: ybox = 2.2
      real*8 :: zbox = 2.2

c     arrays to hold position, velocity, and energy derivatives

      real*8, allocatable :: x(:)
      real*8, allocatable :: y(:)
      real*8, allocatable :: z(:)

      real*8, allocatable :: vx(:)
      real*8, allocatable :: vy(:)
      real*8, allocatable :: vz(:)

      real*8 dx, dy, dz, R, R2, R6

      real*8, allocatable :: derivs(:,:)
      real*8, allocatable :: derivs_old(:,:)

      integer T, i, j, k, ii
      T = nint(N_PARTICLES**(1.0/3.0))

      allocate(x(N_PARTICLES))
      allocate(y(N_PARTICLES))
      allocate(z(N_PARTICLES))

      allocate(vx(N_PARTICLES))
      allocate(vy(N_PARTICLES))
      allocate(vz(N_PARTICLES))

      allocate(derivs(N_PARTICLES,3))
      allocate(derivs_old(N_PARTICLES,3))

      
      do i = 1, T
         do j = 1,T
            do k = 1,T
               ii = (i-1)*T*T + (j-1)*T + k
               x(ii) = xbox/T*i
               y(ii) = ybox/T*j
               z(ii) = zbox/T*k
            end do
         end do
      end do

      do i = 1, N_PARTICLES
         vx(i) = 0
         vy(i) = 0
         vz(i) = 0
      end do

      do i = 1,nstep
         derivs_old = derivs

c        Update Positions
         x = x + vx*dt + derivs(:,1)/(2*mp)*dt*dt
         y = y + vy*dt + derivs(:,2)/(2*mp)*dt*dt
         z = z + vz*dt + derivs(:,3)/(2*mp)*dt*dt

         do j = 1, N_PARTICLES
            do k = j+1, N_PARTICLES
               dx = x(k) - x(j)
               dy = y(k) - y(j)
               dz = z(k) - z(j)

               R = (dx**2 + dy**2 + dz**2)**0.5

               R2 = R**2
               R6 = R**6


               R = 4*e_lj*(-12*sigma_lj**12/(R6*R2)
     &            + 6*sigma_lj**6/R2)/R6



c           computer energy derivates and update

               derivs(j,1) = derivs(j,1) + R*dx
               derivs(k,1) = derivs(k,1) + derivs(j,1)*-1
               derivs(j,2) = derivs(j,2) + R*dy
               derivs(k,2) = derivs(k,2) + derivs(j,2)*-1
               derivs(j,3) = derivs(j,3) + R*dz
               derivs(k,3) = derivs(k,3) + derivs(j,3)*-1
            end do
         end do

c        update velocities

         derivs_old = derivs + derivs_old

         vx = vx + derivs_old(:,1)/(2*mp)*dt
         vy = vy + derivs_old(:,2)/(2*mp)*dt
         vz = vz + derivs_old(:,3)/(2*mp)*dt


         if(mod(i,100) == 0) then

            print*, i

         end if
      end do

      end
