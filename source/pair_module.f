      module pair_module

c***********************************************************************
c     
c     dl_poly module for defining atom pair data
c     copyright - daresbury laboratory
c     author    - w. smith    mar 2004
c     
c***********************************************************************

      use setup_module
      use error_module

      implicit none

      integer, allocatable :: ilist(:),jlist(:)
      real(8), allocatable :: xdf(:),ydf(:),zdf(:)
      real(8), allocatable :: rsqdf(:)

      save ilist,jlist,xdf,ydf,zdf,rsqdf

      contains

      subroutine alloc_pair_arrays(idnode,mxnode)
      
      implicit none
      
      integer, parameter :: nnn=6
      
      logical safe
      integer i,fail,idnode,mxnode
      dimension fail(nnn)
      
      safe=.true.
      
c     allocate arrays
      
      fail(:)=0
      
      allocate (ilist(mxxdf),stat=fail(1))
      allocate (jlist(mxxdf),stat=fail(2))
      allocate (xdf(mxxdf),stat=fail(3))
      allocate (ydf(mxxdf),stat=fail(4))
      allocate (zdf(mxxdf),stat=fail(5))
      allocate (rsqdf(mxxdf),stat=fail(6))
      
      if(any(fail.gt.0))safe=.false.      
      if(mxnode.gt.1)call gstate(safe)    
      if(.not.safe)call error(idnode,1940)
      
      end subroutine alloc_pair_arrays

      end module pair_module
