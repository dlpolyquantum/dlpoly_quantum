      module shake_module
      
c***********************************************************************
c     
c     dl_poly module for defining bond shake arrays
c     
c     copyright - daresbury laboratory
c     author    - w. smith    sep 2003
c     
c***********************************************************************
      
      use error_module
      use parse_module
      use setup_module
      use site_module
      use solvation_module
      
      implicit none
      
      real(8), allocatable :: prmcon(:)
      integer, allocatable :: listcon(:,:),listot(:)
      integer, allocatable :: numcon(:),lstcon(:,:)
      integer, allocatable :: listme(:),lishap(:),lashap(:)
      
      save prmcon,listcon,listot,numcon,lstcon,listme,lishap,lashap
      
      contains
      
      subroutine alloc_shake_arrays(idnode,mxnode)
      
      implicit none
      
      integer, parameter :: nnn=8
      
      logical safe
      integer i,fail,idnode,mxnode
      dimension fail(nnn)

      safe=.true.

c     allocate arrays
      
      fail(:)=0
      
      allocate (prmcon(mxtcon),stat=fail(1))
      allocate (numcon(mxtmls),stat=fail(2))
      allocate (lstcon(mxtcon,2),stat=fail(3))
      allocate (listcon(mxcons,3),stat=fail(4))
      allocate (listme(mxatms),stat=fail(5))
      allocate (lishap(mxlshp),stat=fail(6))
      allocate (lashap(mxnode),stat=fail(7))
      allocate (listot(mxatms),stat=fail(8))
      
      if(any(fail.gt.0))safe=.false.      
      if(mxnode.gt.1)call gstate(safe)    
      if(.not.safe)call error(idnode,1070)

c     initialise numcon array
      
      do i=1,mxtmls
        numcon(i)=0
      enddo
      
      end subroutine alloc_shake_arrays
      
      subroutine define_constraints
     x  (safe,lghost,idnode,itmols,nconst,nsite,natmsr)
      
c***********************************************************************
c     
c     dl_poly subroutine for defining constraints
c     
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2003
c     adapted   - p-a cazade  oct 2007,  solvation, excitation etc.
c     
c***********************************************************************
      
      implicit none
      
      logical safe,lghost
      character*1 message(80)
      integer idnode,itmols,nconst,nsite,ntmp,icnst
      integer icnst1,iatm1,iatm2,isite1,isite2,idum,i
      integer isol1,isol2,natmsr
      
      ntmp=intstr(record,lenrec,idum)
      numcon(itmols)=numcon(itmols)+ntmp
      if(idnode.eq.0) then
        write(nrite,"(/,1x,'number of bond constraints',
     x    5x,i10)") ntmp
        write(nrite,"(/,/,1x,'constraint bond details:',
     x    /,/,21x,5x,'index',5x,'index',2x,'bondlength',/)
     x    ")
      endif
      
      icnst1 = numcon(itmols)
      do icnst=1,icnst1
        
        call getrec(safe,idnode,nfield)
        if(.not.safe)return
        
        call copystring(record,message,80)
        iatm1 = intstr(record,lenrec,idum)
        iatm2 = intstr(record,lenrec,idum)
        
c     test for frozen atom pairs(:)
        
        isite1 = nsite - numsit(itmols) + iatm1
        isite2 = nsite - numsit(itmols) + iatm2
        if(lghost)then
          
          isol1=natmsr+iatm1
          isol2=natmsr+iatm2
          
        endif
        
        if(lfzsit(isite1)*lfzsit(isite2).ne.0) then
          
          numcon(itmols) = numcon(itmols) -1
          if(idnode.eq.0) write(nrite,'(14x,a16,40a1)')
     x      '*** frozen *** ',(message(i),i=1,40)
          
        else
          
          nconst=nconst+1
          
          if(nconst.gt.mxtcon) call error(idnode,40)
          
          lstcon(nconst,1)= iatm1
          lstcon(nconst,2)= iatm2
          prmcon(nconst)=dblstr(record,lenrec,idum)
          
          if(lghost)then
            
            if((isol1.ge.ind_fre(3)).or.(isol2.ge.ind_fre(3)))then
              
              numcon(itmols)=numcon(itmols)-1
              ntcons_ghost=ntcons_ghost+1
              
            endif
            
          endif
          
          if(idnode.eq.0) 
     x      write(nrite,"(21x,2i10,f12.6)")
     x      lstcon(nconst,1),lstcon(nconst,2),
     x      prmcon(nconst)
          
        endif
        
      enddo
      
      return
      end subroutine define_constraints
      
      end module shake_module
