      module exclude_module

c***********************************************************************
c     
c     dl_poly module for defining excluded atom arrays
c     copyright - daresbury laboratory
c     author    - w. smith    sep 2003
c     
c***********************************************************************

      use angles_module
      use bonds_module
      use config_module
      use core_shell_module
      use dihedral_module
      use error_module
      use inversion_module
      use rigid_body_module
      use setup_module
      use shake_module
      use site_module
      use water_module

      implicit none

      integer, allocatable :: lexatm(:,:)
      integer, allocatable :: nexatm(:),noxatm(:)

      save lexatm,nexatm,noxatm

      contains
      
      subroutine alloc_exc_arrays(idnode,mxnode)

c***********************************************************************
c     
c     dl_poly subroutine for allocating excluded atom arrays
c     copyright - daresbury laboratory
c     author    - w. smith    sep 2003
c     
c***********************************************************************
      
      implicit none
      
      integer, parameter :: nnn=3
      
      logical safe
      integer i,fail,idnode,mxnode
      dimension fail(nnn)
      
      safe=.true.
      
c     allocate arrays
      
      fail(:)=0
      allocate (lexatm(msatms,mxexcl),stat=fail(1))
      allocate (nexatm(msatms),stat=fail(2))
      allocate (noxatm(msatms),stat=fail(3))
      
      if(any(fail.gt.0))safe=.false.      
      if(mxnode.gt.1)call gstate(safe)    
      if(.not.safe)call error(idnode,1012)
      
      end subroutine alloc_exc_arrays

      subroutine exclude(idnode,mxnode,natms,ntpmls,lmsite)

c***********************************************************************
c     
c     dl_poly subroutine for constructing the excluded pair
c     interaction list of the system to be simulated
c     
c     keybnd < 0 distance restraint so not excluded
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith        june 1992
c     
c     rigid body exclusions added : t.forester nov 1993
c     check on 1..4 scale factors : t.forester feb 1994
c     inversion terms added       : w.smith    jul 1996
c     
c***********************************************************************
      
      implicit none

      logical newjob,check,safe,lchk,lmsite
      integer idnode,mxnode,natms,ntpmls,nsatms,ibig,i,ibonds,ivoids
      integer iangle,iconst,idihdr,invers,igrp,isite,ishels,ia,ib,k
      integer ic,id,jz,j,jj,ia1,ib1,kk,jk,ntpsit,nlast,lexsav,itmols
      real(8) a1,a2,a3
      
      save newjob
      
      data newjob/.true./
      
      if(newjob)then

c     check on array allocations
        
        nsatms=(natms-1)/mxnode+1
        if(nsatms.gt.msatms)then
          
          if(idnode.eq.0)write(nrite,*)'make msatms >= ',nsatms
          call error(idnode,100)
          
        endif
        
        newjob=.false.
        
      endif
      
c     variables for array bound checking

      ibig=0
      safe=.true.

c     initialise excluded atom arrays
      
      do i=1,mxsite
        
        nexsit(i)=0
        
      enddo
      
      do i=1,msatms
        
        nexatm(i)=0
        
      enddo

      do j=1,mxexcl
        
        do i=1,mxsite
          
          lexsit(i,j)=0
          
        enddo
        
        do i=1,msatms
          
          lexatm(i,j)=0
          
        enddo
        
      enddo
      
c     loop over molecules in system
      
      ibonds=0
      ivoids=0
      iangle=0
      iconst=0
      idihdr=0
      invers=0
      igrp  =0
      isite =0
      ishels=0
      
      do itmols=1,ntpmls
        
c     exclude sites on basis of chemical bonds
        
        do i=1,numbonds(itmols)
          
          ibonds=ibonds+1
          
          if(keybnd(ibonds).gt.0)then
            
            ia=lstbnd(ibonds,1)+isite
            ib=lstbnd(ibonds,2)+isite

            call check_exclude(lchk,safe,ia,ib,isite,ibig,nexsit,lexsit)
            
          endif
          
        enddo
        
c *******************************************************************      
c     M.R.Momeni & F.A.Shakib
c     Method Development and Materials Simulation Laboratory
c
c     exclude intramolecular interactions for M-site in
c     water molecule

        if(lmsite.and.itmols.eq.ntpmls)then

          do i=1,numvoids(ntpmls)

            ivoids=ivoids+1

            if(keyvod(ivoids).gt.0)then

            ia=lstvod(ivoids,1)+isite
            ib=lstvod(ivoids,2)+isite

            call check_exclude(lchk,safe,ia,ib,isite,ibig,nexsit,lexsit)

            endif  

          enddo

        endif

c *******************************************************************      

c     exclude sites on basis of bond constraints
        
        do i=1,numcon(itmols)
          
          iconst=iconst+1
          ia=lstcon(iconst,1)+isite
          ib=lstcon(iconst,2)+isite

          call check_exclude(lchk,safe,ia,ib,isite,ibig,nexsit,lexsit)
          
        enddo
        
c     exclude sites on basis of bond angles
        
        do i=1,numang(itmols)
          
          iangle=iangle+1

          if(keyang(iangle).gt.0)then
            
            ia=lstang(iangle,1)+isite
            ib=lstang(iangle,2)+isite
            ic=lstang(iangle,3)+isite

            call check_exclude(lchk,safe,ia,ib,isite,ibig,nexsit,lexsit)
            call check_exclude(lchk,safe,ib,ic,isite,ibig,nexsit,lexsit)
            call check_exclude(lchk,safe,ia,ic,isite,ibig,nexsit,lexsit)
            
          endif          
          
        enddo
        
c     exclude on basis of rigid groups
        
        do i=1,numgrp(itmols)
          
          igrp=igrp+1
          id=listyp(igrp)
          
          do jj=1,numgsit(id)-1
            
            ia=lstgst(igrp,jj)+isite
            
            do jk=jj+1,numgsit(id)
              
              ib=lstgst(igrp,jk)+isite
              
              call check_exclude(lchk,safe,ia,ib,isite,ibig,nexsit,
     x          lexsit)
              
            enddo
            
          enddo
          
        enddo
        
c     exclude sites on basis of 1-4 dihedral angles
        
        do i=1,numdih(itmols)
          
          idihdr=idihdr+1
          ia=lstdih(idihdr,1)+isite
          ib=lstdih(idihdr,2)+isite
          ic=lstdih(idihdr,3)+isite
          id=lstdih(idihdr,4)+isite

          call check_exclude(lchk,safe,ia,ib,isite,ibig,nexsit,lexsit)
          call check_exclude(lchk,safe,ib,ic,isite,ibig,nexsit,lexsit)
          call check_exclude(lchk,safe,ia,ic,isite,ibig,nexsit,lexsit)
          call check_exclude(lchk,safe,id,ib,isite,ibig,nexsit,lexsit)
          call check_exclude(lchk,safe,id,ic,isite,ibig,nexsit,lexsit)
          call check_exclude(lchk,safe,ia,id,isite,ibig,nexsit,lexsit)
          
c     ia-id interaction: may need to reset vdw and elec scale factors
          
          if(.not.lchk.and.keydih(idihdr).ne.7)then
            
c     if already excluded reset 1..4 vdw and coulombic scale factors
            
            check=((abs(prmdih(idihdr,4)).gt.1.d-10).or.
     x        (abs(prmdih(idihdr,5)).gt.1.d-10))
            
            if(check)then
              
              a1=dble(itmols)
              a2=dble(ia)
              a3=dble(id)
              call warning(idnode,20,a1,a2,a3)
              
              prmdih(idihdr,4)=0.d0
              prmdih(idihdr,5)=0.d0
              
            endif
            
          endif
          
        enddo
        
c     exclude sites on basis of inversion potentials
        
        do i=1,numinv(itmols)
          
          invers=invers+1
          ia=lstinv(invers,1)+isite
          ib=lstinv(invers,2)+isite
          ic=lstinv(invers,3)+isite
          id=lstinv(invers,4)+isite
          
          call check_exclude(lchk,safe,ia,ib,isite,ibig,nexsit,lexsit)
          call check_exclude(lchk,safe,ib,ic,isite,ibig,nexsit,lexsit)
          call check_exclude(lchk,safe,ia,ic,isite,ibig,nexsit,lexsit)
          call check_exclude(lchk,safe,id,ib,isite,ibig,nexsit,lexsit)
          call check_exclude(lchk,safe,id,ic,isite,ibig,nexsit,lexsit)
          call check_exclude(lchk,safe,ia,id,isite,ibig,nexsit,lexsit)
          
        enddo
        
c     exclude sites on basis of core-shell units
        
        do i=1,numshl(itmols)
          
          ishels=ishels+1
          
          ia=lstshl(ishels,1)+isite
          ib=lstshl(ishels,2)+isite

          call check_exclude(lchk,safe,ia,ib,isite,ibig,nexsit,lexsit)
          
c     exclude sites on basis of bonds to core-shell units

          ibonds=ibonds-numbonds(itmols)
          do kk=1,numbonds(itmols)
            
            ibonds=ibonds+1
            
            if(keybnd(ibonds).gt.0)then
              
              ia1=lstbnd(ibonds,1)+isite
              ib1=lstbnd(ibonds,2)+isite

              if(ia.eq.ia1)then

                call check_exclude
     x            (lchk,safe,ib1,ib,isite,ibig,nexsit,lexsit)

              endif

              if(ia.eq.ib1)then

                call check_exclude
     x            (lchk,safe,ia1,ib,isite,ibig,nexsit,lexsit)

              endif

              if(ib.eq.ia1)then

                call check_exclude
     x            (lchk,safe,ia,ib1,isite,ibig,nexsit,lexsit)

              endif
              
              if(ib.eq.ib1)then

                call check_exclude
     x            (lchk,safe,ia,ia1,isite,ibig,nexsit,lexsit)

              endif

            endif

          enddo
          
c     exclude sites on basis of constraint bonds to core-shell units
          
          iconst=iconst-numcon(itmols)
          do kk=1,numcon(itmols)
            
            iconst=iconst+1
            
            ia1=lstcon(iconst,1)+isite
            ib1=lstcon(iconst,2)+isite

            if(ia.eq.ia1)then

              call check_exclude
     x          (lchk,safe,ib1,ib,isite,ibig,nexsit,lexsit)

            endif

            if(ia.eq.ib1)then

              call check_exclude
     x          (lchk,safe,ia1,ib,isite,ibig,nexsit,lexsit)

            endif

            if(ib.eq.ia1)then

              call check_exclude
     x          (lchk,safe,ia,ib1,isite,ibig,nexsit,lexsit)

            endif
            
            if(ib.eq.ib1)then

              call check_exclude
     x          (lchk,safe,ia,ia1,isite,ibig,nexsit,lexsit)

            endif

          enddo
          
c     exclude sites on basis of rigid units involving  core or shell
          
          igrp=igrp-numgrp(itmols)
          do kk=1,numgrp(itmols)
            
            igrp=igrp+1
            
            id=listyp(igrp)
            
            do jj=1,numgsit(id)
              
              ia1=lstgst(igrp,jj)+isite
              if(ia1.eq.ia)then

                do jk=1,numgsit(id)
                  
                  if(jk.ne.jj)then
                    ib1=lstgst(igrp,jk)+isite

                    call check_exclude
     x                (lchk,safe,ib1,ib,isite,ibig,nexsit,lexsit)

                  endif

                enddo
                
              endif

              if(ia1.eq.ib)then

                do jk=1,numgsit(id)
                  
                  if(jk.ne.jj)then
                    ib1=lstgst(igrp,jk)+isite

                    call check_exclude
     x                (lchk,safe,ia,ib1,isite,ibig,nexsit,lexsit)

                  endif

                enddo
                
              endif

            enddo

          enddo

        enddo
        
        isite=isite+numsit(itmols)
        
      enddo
      
      ntpsit=isite

c     check for exceeded array bounds
      
      if(mxnode.gt.1)call gstate(safe)
      if(.not.safe)then
        if(mxnode.gt.1)call gimax(ibig,1,jj)
        if(idnode.eq.0)then
          write(nrite,*)'mxexcl must be at least ',ibig
          write(nrite,*)'mxexcl is currently     ',mxexcl
        endif
        call error(idnode,65)
      endif
      
c     remove redundant entries from exclusion list (none expected!)
      
      do i=1,ntpsit
        
        nlast=nexsit(i)
        do j=1,nexsit(i)-1
          
          if(j.lt.nlast)then
            
            kk=j
            do k=j+1,nexsit(i)
              
              if(lexsit(i,j).eq.lexsit(i,k))then
                
                nlast=nlast-1
                lexsit(i,k)=0
                
              else if(lexsit(i,k).gt.0)then
                
                kk=kk+1
                lexsav=lexsit(i,k)
                lexsit(i,k)=0
                lexsit(i,kk)=lexsav
                
              endif
              
            enddo
            
          endif
          
        enddo
        
        nexsit(i)=nlast
        
      enddo
      
      return
      end subroutine exclude

      subroutine check_exclude(lchk,safe,ia,ib,isite,ibig,nexsit,lexsit)
c***********************************************************************
c     
c     dl_poly subroutine for checking that atom pair ia and ib have
c     a bonded interaction and must be incorporated in the excluded
c     atom list unless identified previously
c     
c     copyright - daresbury laboratory     2016
c     author    - w. smith 2016
c     
c***********************************************************************
      
      implicit none

      logical lchk,safe
      integer ia,ib,isite,ibig,jz
      integer nexsit(mxsite),lexsit(mxsite,mxexcl)
      
      lchk=.true.
      do jz=1,min(nexsit(ia),mxexcl)
        if(lexsit(ia,jz).eq.ib-isite)lchk=.false.
      enddo
      if(lchk)then
        nexsit(ia)=nexsit(ia)+1
        nexsit(ib)=nexsit(ib)+1
        if(max(nexsit(ia),nexsit(ib)).gt.mxexcl)then
          ibig=max(ibig,nexsit(ia),nexsit(ib))
          safe=.false.
        else
          lexsit(ia,nexsit(ia))=ib-isite
          lexsit(ib,nexsit(ib))=ia-isite
        endif
      endif
      
      return
      end subroutine check_exclude

      subroutine excludeneu(idnode,mxnode,nneut)

c***********************************************************************
c     
c     dl_poly subroutine for constructing the excluded pair
c     interaction list of the system to be simulated
c     part 2 - neutral group implementation
c     
c     copyright - daresbury laboratory     1994
c     author    - t. forester        march 1994
c     
c***********************************************************************
      implicit none

      logical lchk
      integer idnode,mxnode,nneut,ibig,iatom,jatom,last,mpm2
      integer npm2,m,ii,im,itmols,inoff,isoff,isit,iolsit,jm,jtmols
      integer jnoff,jsoff,jsit,jolsit,jn1,jno1,jsite,jsite0,in1,ino1
      integer jj0,isite,ij,idum,it
      
c     construct excluded pair list for verlet neighbour correction

      ibig=0
      iatom=0
      jatom=0

c     generate all atomic pairs and check for exclusions 
c     with Brode Ahlrichs ordering of groups
      
      last=nneut
      lchk=.true.
      mpm2=nneut/2+1
      npm2=(nneut-1)/2+1
      
c     outer loop over groups
      
      do m=1,mpm2
        
        if(m.gt.npm2)last=mpm2-1
        
c     inner loop over groups - include intragroup interactions
        
        ii=0
        
        do im=idnode+1,last,mxnode
          
          ii=ii+1

c     first site in neutral group
          
          itmols=1
          inoff=0
          isoff=0
          isit=numsit(itmols)*nummols(itmols)
          iolsit=numsit(itmols)
          
c     calculate j group indices
          
          jm=im+m-1
          if(jm.gt.nneut)jm=jm-nneut
          
c     inner loop over neutral groups
          
          jtmols=1
          jnoff=0
          jsoff=0
          jsit=numsit(jtmols)*nummols(jtmols)
          jolsit=numsit(jtmols)
          
c     test first sites in neutral group
          
          jatom=neulst(jm)        

c     establish pointer to sets
          
          do while(jatom.gt.jsit)
            
            jtmols=jtmols+1
            jnoff=jsit
            jsoff=jsoff+jolsit
            jsit=jsit+nummols(jtmols)*numsit(jtmols)
            jolsit=numsit(jtmols)
            
          enddo
          
          jn1=jatom-jnoff
          jno1=(jn1/jolsit)*jolsit
          jsite=jn1-jno1
          if(jsite.eq.0)then 
            jsite=jolsit
            jno1=jno1-jolsit
          endif
          jsite=jsite+jsoff
          jsite0=jsite-1
          
          do iatom=neulst(im),neulst(im+1)-1
            
c     establish pointer to sets
            
            do while(iatom.gt.isit)
              
              itmols=itmols+1
              inoff=isit
              isoff=isoff+iolsit
              isit=isit+nummols(itmols)*numsit(itmols)
              iolsit=numsit(itmols)
              
            enddo
            
            in1=iatom-inoff
            ino1=(in1/iolsit)*iolsit
            isite=in1-ino1
            if(isite.eq.0)then 
              isite=iolsit
              ino1=ino1-iolsit
            endif
            isite=isite+isoff
            
c     test im and jm are neutral groups on same molecule
            
            if((jnoff.eq.inoff).and.(ino1.eq.jno1))then
            if(abs(im-jm).lt.iolsit)then
              
              jj0=neulst(jm)
              jsite=jsite0
              
c     special case for im=jm (ie. same group)
              
              if(im.eq.jm)then 
                
                jj0=iatom+1
                jsite=isite
                
              endif

c     test for excluded interaction
              
              do jatom=jj0,neulst(jm+1)-1
                
                jsite=jsite+1
                
                do ij=1,nexsit(isite)
                  
                  if(lexsit(isite,ij).eq.jsite-jsoff)then
                    
                    it=nexatm(ii)
                    
                    if(it+2.gt.mxexcl)then
                      
                      ibig=max(it+2,ibig)
                      nexatm(ii)=it+2
                      lchk=.false.
                      
                    else
                      
                      lexatm(ii,it+1)=iatom
                      lexatm(ii,it+2)=jatom
                      nexatm(ii)=nexatm(ii)+2
                      
                    endif

                  endif
                  
                enddo
                
              enddo
              
            endif
            endif
            
          enddo
          
        enddo
        
      enddo

c     global check
      
      call gstate(lchk)
      if(.not.lchk)then
        
        if(mxnode.gt.1)call gimax(ibig,1,idum)
        if(idnode.eq.0)write(nrite,*)'mxexcl must be at least ',ibig
        if(idnode.eq.0)write(nrite,*)'mxexcl is currently     ',mxexcl
        call error(idnode,260)

      endif
      
      return
      end subroutine excludeneu

      subroutine exclude_link(idnode,mxnode,ntpmls)

c***********************************************************************
c     
c     dl_poly subroutine for constructing the excluded pair
c     interaction list of the system to be simulated
c     
c     part 2 - link cell implementation
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith        june 1992
c     
c***********************************************************************
      
      implicit none

      integer idnode,mxnode,ntpmls,iatom,jatom,lsite,ksite
      integer itmols,imols,isite,kk,newatm,k

c     construct excluded pair list for verlet neighbour correction
      
      iatom=0
      jatom=0
      lsite=0
      ksite=0
      
      do itmols=1,ntpmls
        
        do imols=1,nummols(itmols)
          
          do isite=1,numsit(itmols)
            
            iatom=iatom+1
            
            if(mod(iatom-1,mxnode).eq.idnode)then
              
              kk=0
              jatom=jatom+1
              
              do k=1,nexsit(ksite+isite)
                
                newatm=lexsit(ksite+isite,k)+lsite
                
                kk=kk+1
                lexatm(jatom,kk)=newatm
                
              enddo
              
              nexatm(jatom)=kk
              
            endif
            
          enddo
          
          lsite=lsite+numsit(itmols)
          
        enddo
        
        ksite=ksite+numsit(itmols)
        
      enddo
      
      return
      end subroutine exclude_link

      subroutine exclude_atom(idnode,mxnode,natms,ntpmls)

c***********************************************************************
c     
c     dl_poly subroutine for constructing the excluded pair
c     interaction list of the system to be simulated
c     part 2 
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith        june 1992
c     
c***********************************************************************

      implicit none

      integer idnode,mxnode,natms,ntpmls,iatom,jatom,lsite
      integer ksite,itmols,isite,imols,k,kk,newatm,j,latom,i,ii

c     construct excluded pair list for verlet neighbour correction
      
      iatom=0
      jatom=0
      lsite=0
      ksite=0
      
      do itmols=1,ntpmls
        
        do imols=1,nummols(itmols)
          
          do isite=1,numsit(itmols)
            
            iatom=iatom+1
            
            if(mod(iatom-1,mxnode).eq.idnode)then
              
              kk=0
              jatom=jatom+1
              
              do k=1,nexsit(ksite+isite)
                
                newatm=lexsit(ksite+isite,k)+lsite

c     keep only brode-ahlrichs combinations of indices
                
                if(((newatm.gt.iatom).and.
     x            (newatm-iatom.le.natms/2)).or.
     x            ((newatm.lt.iatom).and.
     x            (newatm+natms-iatom.le.(natms-1)/2)))then
                  
                  kk=kk+1
                  lexatm(jatom,kk)=newatm
                  
                  if(kk.gt.1)then
                    
c     sort the excluded atom list in ascending indices

                    do j=kk,2,-1
                      
                      if(lexatm(jatom,j).lt.lexatm(jatom,j-1))
     x                  then
                        latom=lexatm(jatom,j)
                        lexatm(jatom,j)=lexatm(jatom,j-1)
                        lexatm(jatom,j-1)=latom
                      endif
                      
                    enddo
                    
                  endif
                  
                endif
                
              enddo
              
              nexatm(jatom)=kk
              
            endif
            
          enddo
          
          lsite=lsite+numsit(itmols)
          
        enddo
        
        ksite=ksite+numsit(itmols)
        
      enddo
      
c     final sort into brode-ahlrichs ordering
      
      ii=0
      do i=1+idnode,natms,mxnode

        ii=ii+1
        do j=1,nexatm(ii)
          
          if(lexatm(ii,1).lt.i)then
            
            latom=lexatm(ii,1)
            
            do k=1,nexatm(ii)-1
              
              lexatm(ii,k)=lexatm(ii,k+1)
              
            enddo
            
            lexatm(ii,nexatm(ii))=latom
            
          endif
          
        enddo
        
      enddo
      
      return
      end subroutine exclude_atom

      subroutine exclude_copy_mtd(idnode)

c***********************************************************************
c     
c     dl_poly subroutine for copying excluded atom arrays into 
c     the metadynamics module for use in computing order parameters
c
c     author    - d. quigley    April 2012
c     
c***********************************************************************
      use metafreeze_module, only : mtd_lexatm,mtd_nexatm
      implicit none
      integer, parameter :: nnn=2

      integer i,fail,idnode
      dimension fail(nnn)

      do i=1,nnn
        fail(i)=0
      enddo

c     data needed by metadynamics module
      allocate (mtd_lexatm(msatms,mxexcl),stat=fail(1))
      allocate (mtd_nexatm(msatms)       ,stat=fail(2))

      do i=1,nnn
        if(fail(i).gt.0) call error(idnode,1012)
      enddo


c     copy exclude list into metafreeze module
      mtd_nexatm = nexatm
      mtd_lexatm = lexatm

      return
      end subroutine exclude_copy_mtd
      
      end module exclude_module
