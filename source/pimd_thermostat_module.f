      module pimd_thermo_module

c**********************************************************************
c     
c     dl_poly_quantum module for PIGLET thermostat for path integral
c     molecular dynamics
c     
c     reference: Ceriotti, Bussi and Parrinello 
c                Phys. Rev. Lett. 102, 020601 (2009),                         
c                Ceriotti, Bussi and Parrinello
c                Phys. Rev. Lett. 103, 030603 (2009),  
c                Ceriotti, G. Bussi and M. Parrinello  
c                J. Chem. Theory Comput. 6, 1170 (2010)     
c     
c     copyright - 
c     authors   - Dil Limbu and Nathan London 2023
c     
c**********************************************************************
      
      implicit none
      use setup_module,   only : pi,boltz,hbar,mspimd,nrite
      use pimd_module,    only : zmass,rzmass,pxx,pyy,pzz,nbeads
     x                           uxx,uyy,uzz,
     x                           wxx,wyy,wzz
      use error_module,   only : error
      use utility_module, only : puni
      
      public 
      
      contains
      
c     Initialize the GLE thermostat by allocating and populating several
c     temporary arrays.
      subroutine gle_initialize(dt, Natoms, Nbeads, A, C, Ns)

      integer, intent(in) :: Natoms, Nbeads, Ns
      double precision, intent(in) :: dt, A(Ns+1,Ns+1), C(Ns+1,Ns+1)

      double precision :: gr(Ns+1), C1(Ns+1,Ns+1)
      integer :: i, j, k, s

c     Allocate arrays
      allocate(gle_S(Ns+1,Ns+1))
      allocate(gle_T(Ns+1,Ns+1))
      allocate(gle_p(3,Natoms,Nbeads,Ns+1))
      allocate(gle_np(3,Natoms,Nbeads,Ns+1))

c      Determine the deterministic part of the propagator
      call matrix_exp(-dt*A, Ns+1, 15, 15, gle_T)

c        ! Determine the stochastic part of the propagator
      call cholesky(C - matmul(gle_T, matmul(C, transpose(gle_T))), gle_S, Ns+1)

c        ! Initialize the auxiliary noise vectors
c        ! To stay general, we use the Cholesky decomposition of C; this allows
c        ! for use of non-diagonal C to break detailed balance
c        ! We also use an extra slot for the physical momentum, as we could then
c        ! use it to initialize the momentum in the calling code
      call cholesky(C, C1, Ns+1)
      do i = 1, 3
        do j = 1, Natoms
          do k = 1, Nbeads
            do s = 1, Ns+1
              call randomn(gr(s))
            end do
            gle_p(i,j,k,:) = matmul(C1, gr)
          end do
        end do
      end do

      end subroutine gle_initialize

c    ! Apply the GLE thermostat to the momentum.
     subroutine gle_thermostat(p,mass,beta,dxi,Natoms,Nbeads,Ns,constrain,result)

      implicit none
      integer, intent(in) :: Natoms, Nbeads, Ns
      double precision, intent(in) :: mass(Natoms), beta, dxi(3,Natoms)
      double precision, intent(inout) :: p(3,Natoms,Nbeads)
      integer, intent(in) :: constrain
      integer, intent(inout) :: result

      double precision :: p0(3,Natoms,Nbeads), p00(3,Natoms,Nbeads), cgj
      integer :: i, j, k, s
      integer :: N, nmrep

      nmrep = 0

      N = 3 * Natoms * Nbeads

c      ! Switch to mass-scaled coordinates when storing momenta in gle_p
      do j = 1, Natoms
         gle_p(:,j,:,1) = p(:,j,:) / dsqrt(mass(j))
      end do

c      ! We pretend that gp is a (3*Natoms*Nbeads)x(Ns+1) matrix, which should be fine...
      call dgemm('N','T', N, Ns+1, Ns+1, 1.0d0, gle_p, N, gle_T, Ns+1, 0.0d0, gle_np, N)

c      ! Compute the random part
      do s = 1, Ns+1
        do i = 1, 3
          do k = 1, Nbeads
            cgj = 1.0d0    
c          ! This is to make it work when applied in NM representation
            if (nmrep .gt. 0 .and. k .ne. 1 .and. (mod(Nbeads,2) .ne. 0 .or. k .ne. (Nbeads/2+1))) then
              cgj = dsqrt(0.5d0)
            end if
            do j = 1, Natoms
              call randomn(gle_p(i,j,k,s))
              gle_p(i,j,k,s) = gle_p(i,j,k,s) * cgj
            end do
          end do
        end do
      end do

c      !Again we pretend that gp is a (3*Natoms*Nbeads)x(Ns+1) matrix, which should be fine...
      call dgemm('N','T', N, Ns+1, Ns+1, 1.0d0, gle_p, N, gle_S, Ns+1, 1.0d0, gle_np, N)
      gle_p = gle_np

c      !Switch back from mass-scaled coordinates when recovering momenta from gle_p
      do j = 1, Natoms
        p(:,j,:) = gle_p(:,j,:,1) * dsqrt(mass(j))
      end do

      end subroutine gle_thermostat

c    ! Clean up the GLE thermostat by deallocating temporary arrays.
      subroutine gle_cleanup()

      deallocate(gle_S, gle_T, gle_p, gle_np)

      end subroutine gle_cleanup


      subroutine matrix_exp(M,n,j,k,EM)

c**********************************************************************
c     
c     routine to compute the expoential of a (square) matrix using 
c     scale and square algorithm.
c     Parameters:
c     M - The square matrix to exponentiate
c     n - The size of the matrix
c     j - The number of terms in the exponential to evaluate
c     k - The exponent of the power of two used to scale the matrix
c     Returns:
c     EM - The expoential of the given matrix
c     
c     copyright - 
c     authors   - 
c     
c**********************************************************************

      implicit none
      integer, intent(in)  :: n, j, k
      real(8), intent(in)  :: M(n,n)
      real(8), intent(out) :: EM(n,n)

      integer :: i,p
      real(8) :: tc(j+1),SM(n,n)

      tc(1)=1.d0
      do i=1,j
        tc(i+1)=tc(i)/dble(i)
      enddo

c    !scale
      SM=M*(1.d0/2.d0**k)
      EM=0.d0
      do i=1,n
        EM(i,i)=tc(j+1)
      enddo

c    !taylor exp of scaled matrix
      do p=j,1,-1
        EM=matmul(SM,EM);
        do i=1,n
          EM(i,i)=EM(i,i)+tc(p)
        enddo
      enddo

c    !square
      do p=1,k
        EM=matmul(EM,EM)
      enddo

      end subroutine matrix_exp


      subroutine cholesky(SST,S,n)

c**********************************************************************
c     
c     routine to compute the brute-force stabilized Cholesky 
c     decomposition of a square matrix.
c     Parameters:
c     SST - The matrix to determine the Cholesky decomposition of
c     n - The size of the matrix
c     Returns:
c     S - The computed Cholesky decomposition
c     
c     copyright - 
c     authors   - 
c     
c**********************************************************************

      implicit none
      integer, intent(in)  :: n
      real(8), intent(in)  :: SST(n,n)
      real(8), intent(out) :: S(n,n)
    
      integer  :: i,j,k
      real(8)  :: D(n), L(n,n)
   
      D=0.d0
      L=0.d0
      S=0.d0

      do i=1,n
        L(i,i)=1.d0
        D(i)=SST(i,i)
        do j=1,i-1
          L(i,j)=SST(i,j)
          do k=1,j-1
            L(i,j)=L(i,j)-L(i,k)*L(j,k)*D(k)
          end do
          if (D(j).ne.0.) L(i,j)=L(i,j)/D(j)
        end do
        do k=1,i-1
          D(i)=D(i)-L(i,k)*L(i,k)*D(k)
        end do
      end do

      do i=1,n
        do j=1,i
          if (D(j)>0.d0) S(i,j)=S(i,j)+L(i,j)*sqrt(D(j))
        end do
      end do

      end subroutine cholesky




      end module pimd_thermo_module
