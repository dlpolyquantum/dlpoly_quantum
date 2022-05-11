      module rigid_body_module
      
c***********************************************************************
c     
c     dl_poly module for defining rigid body arrays
c     copyright - daresbury laboratory
c     author    - w. smith    oct 2003
c     
c***********************************************************************
      
      use config_module
      use error_module
      use parse_module
      use setup_module
      use site_module
      use solvation_module
      
      implicit none
      
      real(8), allocatable :: omx(:),omy(:),omz(:)
      real(8), allocatable :: gcmx(:),gcmy(:),gcmz(:)
      real(8), allocatable :: gvxx(:),gvyy(:),gvzz(:),gmass(:)
      real(8), allocatable :: q0(:),q1(:),q2(:),q3(:)
      real(8), allocatable :: gxx(:,:),gyy(:,:),gzz(:,:)
      real(8), allocatable :: rotinx(:,:),rotiny(:,:),rotinz(:,:)
      integer, allocatable :: lstrgd(:),numgsit(:),lstgtp(:)
      integer, allocatable :: listyp(:),lstgst(:,:),lstfre(:)
      integer, allocatable :: lstme(:),lstbod(:),lstcsit(:)
      
      save omx,omy,omz,gcmx,gcmy,gcmz,gvxx,gvyy,gvzz,gmass
      save q0,q1,q2,q3,gxx,gyy,gzz,rotinx,rotiny,rotinz
      save lstrgd,numgsit,lstgtp,listyp,lstgst,lstfre,lstme
      save lstbod,lstcsit
      
      contains
      
      subroutine alloc_rgbdy_arrays(idnode,mxnode)
      
      implicit none
      
      integer, parameter :: nnn=12
      
      logical safe
      integer i,fail,idnode,mxnode
      dimension fail(nnn)

      safe=.true.

c     allocate arrays

      fail(:)=0
      
      allocate (omx(mxgrp),omy(mxgrp),omz(mxgrp),stat=fail(1))
      allocate (gcmx(mxgrp),gcmy(mxgrp),gcmz(mxgrp),stat=fail(2))
      allocate (gvxx(mxgrp),gvyy(mxgrp),gvzz(mxgrp),stat=fail(3))
      allocate (q0(mxgrp),q1(mxgrp),q2(mxgrp),q3(mxgrp),stat=fail(4))
      allocate (gxx(mxungp,mxngp),gyy(mxungp,mxngp),stat=fail(5))
      allocate (gzz(mxungp,mxngp),gmass(mxungp),stat=fail(6))
      allocate (rotinx(mxungp,2),rotiny(mxungp,2),stat=fail(7))
      allocate (rotinz(mxungp,2),lstgtp(mxgrp),stat=fail(8))
      allocate (lstrgd(mxgatm),numgsit(mxungp),stat=fail(9))
      allocate (listyp(mxungp),lstgst(mxungp,mxngp),stat=fail(10))
      allocate (lstfre(mxatms),lstme(mxatms),stat=fail(11))
      allocate (lstbod(mxatms),lstcsit(2*mxcons),stat=fail(12))

      if(any(fail.gt.0))safe=.false.
      if(mxnode.gt.1)call gstate(safe)
      if(.not.safe)call error(idnode,1013)
      
      end subroutine alloc_rgbdy_arrays
      
      subroutine define_rigid_body
     x  (safe,lghost,idnode,itmols,ngrp,natmsr)
      
c***********************************************************************
c     
c     dl_poly subroutine for defining rigid bodies
c     
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2003
c     adapted   - p-a cazade  oct 2007,  solvation, excitation etc.
c     
c***********************************************************************
      
      implicit none
      
      logical safe,lghost,site_test
      integer idnode,itmols,ngrp,ntmp,idum,igrp
      integer j,site,natmsr
      
      ngrp_ghost=0
      ntmp=intstr(record,lenrec,idum)
      numgrp(itmols)=numgrp(itmols)+ntmp
      if(idnode.eq.0) then
        write(nrite,"(/,1x,'number of rigid units    ',
     x    6x,i10)") ntmp
        write(nrite,"(/,' rigid body details:',/,/,21x,
     x    6x,'unit',3x,'indices',/) ")
      endif
      
      do igrp=1,numgrp(itmols)
        
        ngrp=ngrp+1
        site_test=.true.
        
        if(ngrp.gt.mxungp) call error(idnode,301)
        
        call getrec(safe,idnode,nfield)
        if(.not.safe)return
        numgsit(ngrp)=intstr(record,lenrec,idum)
        if(numgsit(ngrp).gt.mxngp) 
     x    call error (idnode,302)
        
        listyp(ngrp)=ngrp
        
        do j=1,numgsit(ngrp)
          
          site=intstr(record,lenrec,idum)
          
          if(site.eq.0)then
            
            call getrec(safe,idnode,nfield)
            if(.not.safe)return
            site=intstr(record,lenrec,idum)
            
          endif
          
          lstgst(ngrp,j)=site
          
          if(lghost)then
            
            if(site_test)then
              
              if(site+natmsr.ge.ind_fre(3))then
                
                site_test=.false.
                ngrp_ghost=ngrp_ghost+1
                
              endif
              
            endif
            
          endif
          
        enddo
        
        if(idnode.eq.0) 
     x    write(nrite,"(21x,10i10,100(/,21x,10i10))")
     x    listyp(ngrp),(lstgst(ngrp,j),j=1,
     x    numgsit(ngrp))
        
      enddo
      
      numgrp(itmols)=numgrp(itmols)-ngrp_ghost
      ngrp=ngrp-ngrp_ghost
        
      return
      end subroutine define_rigid_body
      
      subroutine bodystress
     x  (idnode,mxnode,ngrp,vircom,strbod,dtx,dty,dtz)
      
c***********************************************************************
c     
c     dlpoly routine to calculate rigid body contributions to the 
c     stress tensor
c     
c     copyright daresbury laboratory
c     author      w.smith   aug 2005
c     
c**********************************************************************
      
      implicit none
      
      integer i,j,ig,id,jr,igrp1,igrp2,idnode,mxnode,ngrp
      real(8) vircom,strbod,dtx,dty,dtz
      
      dimension dtx(mxatms),dty(mxatms),dtz(mxatms),strbod(9)
      
c     group block indices
      
      igrp1=(idnode*ngrp)/mxnode+1
      igrp2=((idnode+1)*ngrp)/mxnode
      
c     zero stress tensor accumulators
      
      vircom=0.d0
      do i=1,9
        strbod(i)=0.d0
      enddo
      
c     convert atomic virial to molecular
c     note convention: virial(atom-atom)=-sum(Ri.Fi)
c     : virial(com-com)=-sum(Rcom.Fcom) so
c     virial(com-com)=virial(atom-atom)+sum((Ri-Rcom).Fi)
      
      jr=0
      do ig=igrp1,igrp2
        
        id=lstgtp(ig)
        
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          
          vircom=vircom+
     x      (dtx(jr)*fxx(i)+dty(jr)*fyy(i)+dtz(jr)*fzz(i))
          
c     stress tensor : rigid body contributions
          
          strbod(1)=strbod(1)-dtx(jr)*fxx(i)
          strbod(2)=strbod(2)-dtx(jr)*fyy(i)
          strbod(3)=strbod(3)-dtx(jr)*fzz(i)
          strbod(4)=strbod(4)-dty(jr)*fxx(i)
          strbod(5)=strbod(5)-dty(jr)*fyy(i)
          strbod(6)=strbod(6)-dty(jr)*fzz(i)
          strbod(7)=strbod(7)-dtz(jr)*fxx(i)
          strbod(8)=strbod(8)-dtz(jr)*fyy(i)
          strbod(9)=strbod(9)-dtz(jr)*fzz(i)
          
        enddo
        
      enddo
      
      if(mxnode.gt.1)then
        
        call gdsum(strbod,9,buffer)
        buffer(1)=vircom
        call gdsum(buffer(1),1,buffer(2))
        vircom=buffer(1)
        
      endif
      
c     symmetrise stress tensor
      
      strbod(2)=0.5d0*(strbod(2)+strbod(4))
      strbod(4)=strbod(2)
      strbod(3)=0.5d0*(strbod(3)+strbod(7))
      strbod(7)=strbod(3)
      strbod(6)=0.5d0*(strbod(6)+strbod(8))
      strbod(8)=strbod(6)
      
      return
      end subroutine bodystress
      
      end module rigid_body_module
