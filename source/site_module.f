      module site_module

c***********************************************************************
c     
c     dl_poly module for defining atomic/site arrays
c     
c     copyright - daresbury laboratory
c     author    - w. smith    sep 2003
c     
c***********************************************************************

      use error_module
      use parse_module
      use setup_module

      implicit none

      character*1, allocatable :: molnam(:,:)
      character*8, allocatable :: sitnam(:),unqatm(:)
      real(8), allocatable :: dens(:),chgsit(:),wgtsit(:)
      integer, allocatable :: nexsit(:),lfzsit(:),numsit(:),ltpsit(:)
      integer, allocatable :: nugrp(:),lexsit(:,:),numgrp(:)
      integer, allocatable :: numtyp(:),numfrz(:),nummols(:)

      save numtyp,numfrz,dens,chgsit,wgtsit,sitnam,unqatm,nexsit
      save lfzsit,numsit,ltpsit,nugrp,lexsit,numgrp,molnam,nummols

      contains

      subroutine alloc_site_arrays(idnode,mxnode)
      
      implicit none
      
      integer, parameter :: nnn=16
      
      logical safe
      integer i,fail,idnode,mxnode
      dimension fail(nnn)
      
      safe=.true.

c     allocate arrays

      fail(:)=0
      
      allocate (chgsit(mxsite),stat=fail(1))
      allocate (wgtsit(mxsite),stat=fail(2))
      allocate (nexsit(mxsite),stat=fail(3))
      allocate (lfzsit(mxsite),stat=fail(4))
      allocate (nugrp(mxsite) ,stat=fail(5))
      allocate (ltpsit(mxsite),stat=fail(6))
      allocate (numsit(mxtmls),stat=fail(7))
      allocate (lexsit(mxsite,mxexcl),stat=fail(8))
      allocate (sitnam(mxsite),stat=fail(9))
      allocate (unqatm(mxsite),stat=fail(10))
      allocate (numgrp(mxtmls),stat=fail(11))
      allocate (numtyp(mxatyp),stat=fail(12))
      allocate (numfrz(mxatyp),stat=fail(13))
      allocate (dens(mxatyp),stat=fail(14))
      allocate (nummols(mxtmls),stat=fail(15))
      allocate (molnam(40,mxtmls),stat=fail(16))
      
      if(any(fail.gt.0))safe=.false.      
      if(mxnode.gt.1)call gstate(safe)    
      if(.not.safe)call error(idnode,1090)

c     initialise numsit array
      
      do i=1,mxtmls
         numsit(i)=0
      enddo
      
      end subroutine alloc_site_arrays

      subroutine define_atoms
     x  (safe,lneut,lfreeze,idnode,itmols,nsite,ksite,ntpatm)

c***********************************************************************
c     
c     dl_poly subroutine for  defining atom types in system
c     
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2003
c     
c***********************************************************************

      implicit none

      character*8 atom1
      character*1 message(80)
      logical lneut,safe,atmchk,lfreeze
      integer idnode,itmols,nsite,ksite,ntpatm,isite,nrept
      integer ifrz,neugp,irept,jsite,idum,ntmp
      real(8) weight,charge

      ntmp=intstr(record,lenrec,idum)
      numsit(itmols)=ntmp

      if(idnode.eq.0) then
        write(nrite,"(/,1x,'number of atoms/sites',
     x    10x,i10)") numsit(itmols)
        if(.not.lneut)
     x    write(nrite,"(/,/,1x,'atomic characteristics:',
     x    /,/,21x,' site',5x,'name',10x,'mass',8x,
     x    'charge',4x,'repeat',4x,'freeze'/)")
        if(lneut)
     x    write(nrite,"(/,/,1x,'atomic characteristics:',/
     x    /,21x,' site',5x,'name',10x,'mass',8x,'charge',
     x    4x,'repeat',4x,'freeze',3x,'chg grp')")
      endif
      
      do isite=1,ntmp

        if(ksite.lt.numsit(itmols))then

c     read atom name, site number, mass, charge, freeze option
          
          call getrec(safe,idnode,nfield)
          if(.not.safe)return

          call copystring(record,message,80)
          call getword(atom1,record,8,lenrec)
          weight=dblstr(record,lenrec,idum)
          charge=dblstr(record,lenrec,idum)
          nrept=intstr(record,lenrec,idum)
          ifrz =intstr(record,lenrec,idum)
          if(ifrz.gt.0)lfreeze=.true.
          neugp=intstr(record,lenrec,idum)
          if(nrept.eq.0)nrept=1
          ksite=ksite+nrept
          
          if(idnode.eq.0) then
            
            if(.not.lneut) then

              write(nrite,
     x          "(21x,i5,5x,a8,2f12.5,2i10)")
     x          nsite+1,atom1,weight,charge,nrept,
     x          ifrz

            else

              write(nrite,
     x          "(21x,i5,5x,a8,2f12.5,3i10)")
     x          nsite+1,atom1,weight,charge,nrept,
     x          ifrz,neugp

            endif

          endif
          
          do irept=1,nrept
            
            nsite=nsite+1
            if(nsite.gt.mxsite) call error(idnode,20)
            
            sitnam(nsite)=atom1
            wgtsit(nsite)=weight
            chgsit(nsite)=charge
            lfzsit(nsite)=ifrz
            nugrp(nsite)=neugp
            
          enddo
          
c     establish list of unique atom types
          
          atmchk=.true.
          
          do jsite=1,ntpatm
            
            if(atom1.eq.unqatm(jsite)) then
              
              atmchk=.false.
              do irept=nsite,nsite-nrept+1,-1
                
                ltpsit(irept)=jsite
                
              enddo
              
            endif
            
          enddo
          
          if(atmchk)then
            
            ntpatm=ntpatm+1
            if(ntpatm.gt.mxatyp)call error(idnode,14)
            unqatm(ntpatm)=atom1
            
            do irept=nsite,nsite-nrept+1,-1
              
              ltpsit(irept)=ntpatm
              
            enddo
            
          endif
          
        endif
        
      enddo
      
      return
      end subroutine define_atoms

      subroutine check_syschg(idnode,ntpmls,sumchg)

c***********************************************************************
c     
c     dl_poly subroutine for checking the system charge
c     
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2003
c     
c***********************************************************************

      implicit none

      integer idnode,ntpmls,jsite,itmols,lsite
      real(8) sumchg

      jsite=0
      do itmols=1,ntpmls
        
        do lsite=1,numsit(itmols)
          
          jsite=jsite+1
          sumchg=sumchg+dble(nummols(itmols))*chgsit(jsite)
          
        enddo
        
      enddo
      
      if(abs(sumchg).gt.1.0d-6) then
        
        call warning(idnode,60,sumchg,0.d0,0.d0)
        
      endif
      
      return
      end subroutine check_syschg
      
      end module site_module
