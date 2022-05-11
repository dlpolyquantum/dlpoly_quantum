      module config_module

c***********************************************************************
c     
c     dl_poly module for defining simulation configuration data
c     copyright - daresbury laboratory
c     author    - w. smith    sep 2003
c     
c***********************************************************************

      use setup_module
      use error_module
      
      implicit none

      character*1 cfgname(80)
      character*1 sysname(80)
      real(8) cell(9),rcell(9),celprp(10)
      real(8) eta(9),stress(9),stresl(9),strcns(9),strbod(9)
      
      character*8, allocatable :: atmnam(:)
      real(8), allocatable :: xxx(:),yyy(:),zzz(:)
      real(8), allocatable :: vxx(:),vyy(:),vzz(:)
      real(8), allocatable :: fxx(:),fyy(:),fzz(:)
      real(8), allocatable :: flx(:),fly(:),flz(:)
      real(8), allocatable :: chge(:),weight(:),rmass(:)
      integer, allocatable :: ltype(:),lstfrz(:)
      integer, allocatable :: neulst(:),lstneu(:)
      integer, allocatable :: lentry(:),list(:,:)
      integer, allocatable :: lstout(:),link(:)
      integer, allocatable :: lct(:),lst(:)

      real(8), allocatable :: buffer(:)

      save atmnam,neulst,lstneu,cfgname,sysname
      save cell,xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz
      save buffer,weight,chge,ltype,lstfrz,flx,fly,flz
      save lentry,list,lstout,link,lct,lst,celprp,rmass
      save eta,stress,stresl,strcns,rcell
      
      contains
      
      subroutine alloc_config_arrays(idnode,mxnode)

c***********************************************************************
c     
c     dl_poly subroutine for defining simulation configuration arrays
c     copyright - daresbury laboratory
c     author    - w. smith    sep 2003
c     
c***********************************************************************
      
      integer, parameter :: nnn=27

      logical safe
      integer i,fail,idnode,mxnode,numatm
      dimension fail(nnn)

      safe=.true.
      numatm=mxatms*nbeads
      
c     allocate arrays
      
      fail(:)=0
      
      allocate (xxx(numatm),stat=fail(1))
      allocate (yyy(numatm),stat=fail(2))
      allocate (zzz(numatm),stat=fail(3))
      allocate (vxx(numatm),stat=fail(4))
      allocate (vyy(numatm),stat=fail(5))
      allocate (vzz(numatm),stat=fail(6))
      allocate (fxx(numatm),stat=fail(7))
      allocate (fyy(numatm),stat=fail(8))
      allocate (fzz(numatm),stat=fail(9))
      allocate (weight(numatm),stat=fail(11))
      allocate (chge(numatm),stat=fail(12))
      allocate (ltype(numatm),stat=fail(13))
      allocate (lstfrz(numatm),stat=fail(14))
      allocate (flx(numatm),stat=fail(15))
      allocate (fly(numatm),stat=fail(16))
      allocate (flz(numatm),stat=fail(17))
      allocate (atmnam(numatm),stat=fail(18))
      allocate (neulst(mxneut),stat=fail(19))
      allocate (lstneu(numatm),stat=fail(20))
      allocate (lstout(numatm),stat=fail(21))
      allocate (lentry(mslist),stat=fail(22))
      allocate (list(mslist,mxlist),stat=fail(23))
      allocate (link(numatm),stat=fail(24))
      allocate (lct(mxcell),stat=fail(25))
      allocate (lst(mxcell),stat=fail(26))
      allocate (rmass(numatm),stat=fail(27))
      allocate (buffer(mxbuff),stat=fail(10))

      if(any(fail.gt.0))safe=.false.      
      if(mxnode.gt.1)call gstate(safe)    
      if(.not.safe)then
        
        if(idnode.eq.0)write(nrite,'(10i5)')fail
        call error(idnode,1000)

      endif
      
      end subroutine alloc_config_arrays

      end module config_module
