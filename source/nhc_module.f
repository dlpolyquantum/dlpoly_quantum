      module nhc_module
      
c**********************************************************************
c     
c     dl_poly_quantum module for defining Nose-Hoover Chain arrays for
c     different ensembles
c     
c     authors    - M.R. Momeni & F.A. Shakib     
c     copyright  - M.R. Momeni & F.A. Shakib
c
c     Method Development and Materials Simulation Laboratory
c     New Jersey Institute of Technology
c     
c**********************************************************************

      use error_module,   only : error
      
      implicit none

      real(8) v_epsilon
      real(8), allocatable, save :: eta_nhc(:)
      real(8), allocatable, save :: peta(:)
      real(8), allocatable, save :: ksi(:)
      real(8), allocatable, save :: pksi(:)
      
      public alloc_nhc_arrays,dealloc_nhc_arrays
      
      contains
      
      subroutine alloc_nhc_arrays(idnode,mxnode,nchain)
      
c**********************************************************************
c     
c     dl_poly quantum routine to allocate arrays for Nose-Hoover Chain
c     
c     authors    - M.R. Momeni & F.A. Shakib     
c     copyright  - M.R. Momeni & F.A. Shakib 2021
c
c     Method Development and Materials Simulation Laboratory
c     
c**********************************************************************

      implicit none

      logical safe
      integer, intent(in) :: idnode,mxnode,nchain
      integer, dimension(1:4) :: fail

      safe=.true.

c     allocate arrays
      fail(:)=0
      
      allocate (eta_nhc(1:nchain),stat=fail(1))
      allocate (peta(1:nchain),stat=fail(2))
      allocate (ksi(1:nchain),stat=fail(3))
      allocate (pksi(1:nchain),stat=fail(4))
      
      if(any(fail.gt.0))safe=.false.
      if(mxnode.gt.1)call gstate(safe)
      if(.not.safe)call error(idnode,3001)
      
      end subroutine alloc_nhc_arrays
      
      subroutine dealloc_nhc_arrays(idnode,mxnode)
      
c**********************************************************************
c     
c     dl_poly quantum routine to deallocate arrays for Nose-Hoover Chain
c     
c     authors    - M.R. Momeni & F.A. Shakib     
c     copyright  - M.R. Momeni & F.A. Shakib 2021
c
c     Method Development and Materials Simulation Laboratory
c     
c**********************************************************************
      
      implicit none

      logical safe
      integer, intent(in) :: idnode,mxnode
      integer, dimension(2) :: fail
      
      fail(:)=0
      safe=.true.
      
      deallocate(eta_nhc,peta,stat=fail(1))
      deallocate(ksi,pksi,stat=fail(2))
      
      if(any(fail.gt.0))safe=.false.
      if(mxnode.gt.1)call gstate(safe)
      if(.not.safe)call error(idnode,3002)
        
      end subroutine dealloc_nhc_arrays
    
      subroutine nhc_init(idnode,mxnode,nchain)
      
c**********************************************************************
c     
c     dl_poly quantum routine to initialise NVT-NHC thermostat
c     
c     authors    - M.R. Momeni & F.A. Shakib     
c     copyright  - M.R. Momeni & F.A. Shakib 2021
c
c     Method Development and Materials Simulation Laboratory
c     
c**********************************************************************
      
      implicit none
      
      integer i
      integer, intent(in) :: idnode,mxnode,nchain
      

        do i=1,nchain
          eta_nhc(i)=0.d0
          peta(i)=0.d0
          ksi(i)=0.d0
          pksi(i)=0.d0
        enddo

        v_epsilon=0.d0

      end subroutine nhc_init

      end module nhc_module 
