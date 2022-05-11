      module core_shell_module

c***********************************************************************
c     
c     dl_poly module for defining core_shell arrays
c     copyright - daresbury laboratory
c     author    - w. smith    sep 2003
c     
c***********************************************************************

      use config_module
      use error_module
      use pair_module
      use parse_module
      use property_module
      use rigid_body_module
      use setup_module
      use site_module
      use solvation_module

      implicit none

      real(8), allocatable :: prmshl(:,:)
      integer, allocatable :: listshl(:,:)
      integer, allocatable :: numshl(:),lstshl(:,:)

      save prmshl,listshl,numshl,lstshl

      contains
      
      subroutine alloc_csh_arrays(idnode,mxnode)

c***********************************************************************
c     
c     dl_poly subroutine for defining core_shell arrays
c     copyright - daresbury laboratory
c     author    - w. smith    sep 2003
c     
c***********************************************************************
      
      implicit none

      logical safe
      integer i,fail,idnode,mxnode
      dimension fail(4)

      safe=.true.
      
c     allocate arrays
      
      fail(:)=0

      allocate (prmshl(mxtshl,2),stat=fail(1))
      allocate (numshl(mxtmls),stat=fail(2))
      allocate (lstshl(mxtshl,2),stat=fail(3))
      allocate (listshl(mxshl,3),stat=fail(4))
      
      if(any(fail.gt.0))safe=.false.      
      if(mxnode.gt.1)call gstate(safe)    
      if(.not.safe)call error(idnode,1100)

      do i=1,mxtmls
         numshl(i)=0
      enddo

      end subroutine alloc_csh_arrays

      subroutine define_core_shell
     x  (safe,idnode,itmols,nshels,nsite,keyshl,engunit)

c***********************************************************************
c     
c     dl_poly subroutine for defining core-shell units
c     
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2003
c     
c***********************************************************************

      logical safe
      integer idnode,itmols,nshels,nsite,ntmp,ishls
      integer iatm1,iatm2,isite1,isite2,keyshl,kshl,idum
      real(8) engunit

      ntmp=intstr(record,lenrec,idum)
      numshl(itmols)=numshl(itmols)+ntmp
      kshl=intstr(record,lenrec,idum)
      if(keyshl.eq.0)then
        keyshl=kshl
      elseif(kshl.ne.keyshl)then
        call error(idnode,1960)
      endif
      if(idnode.eq.0) then
        
        write(nrite,
     x    "(/,1x,'number of core-shell units',5x,i10)")
     x    ntmp
        if(keyshl.eq.1)then

          write(nrite,
     x       "(/,/,1x,'core-shell details:',/,/,21x,
     x       5x,'index',5x,'index',6x,'parameter')")

        else

          write(nrite,
     x       "(/,/,1x,'core-shell details:',/,/,21x,
     x       6x,'core',5x,'shell',6x,'parameter')")
        
        endif

      endif
      
      do ishls=1,numshl(itmols)
        
        nshels=nshels+1
        if(nshels.gt.mxtshl) call error(idnode,57)
        
        call getrec(safe,idnode,nfield)
        if(.not.safe)return

        iatm1=intstr(record,lenrec,idum)
        iatm2=intstr(record,lenrec,idum)
        lstshl(nshels,1)=iatm1
        lstshl(nshels,2)=iatm2
        prmshl(nshels,1)=dblstr(record,lenrec,idum)
        prmshl(nshels,2)=dblstr(record,lenrec,idum)
        if(idnode.eq.0) write(nrite,
     x    "(21x,2i10,2f15.4)")
     x    lstshl(nshels,1),lstshl(nshels,2),
     x    prmshl(nshels,1),prmshl(nshels,2)

c     test for frozen cores or shells
        
        isite1=nsite-numsit(itmols)+iatm1
        isite2=nsite-numsit(itmols)+iatm2
        if(lfzsit(isite1)*lfzsit(isite2).ne.0)
     x    call error(idnode,49)
        
c     convert energy units to internal units
        
        prmshl(nshels,1)=prmshl(nshels,1)*engunit
        prmshl(nshels,2)=prmshl(nshels,2)*engunit
        
      enddo

      return
      end subroutine define_core_shell

      subroutine corshl(idnode,mxnode,ntshl,shlke)
      
c***********************************************************************
c     
c     dl_poly subroutine for calculating the internal kinetic
c     energy of core-shell units in the shell polarisation model
c     
c     copyright - daresbury laboratory 1994
c     author    - w. smith july 1994
c
c***********************************************************************
      
      implicit none

      integer  idnode,mxnode,ntshl,ishl1,ishl2,i,j,k,m
      real(8) shlke,ppp,ccc,sss
      
      shlke=0.d0

c     block indices

      ishl1=(idnode*ntshl)/mxnode+1
      ishl2=((idnode+1)*ntshl)/mxnode

c     loop over all specified core-shell pairs
      
      m=0

      do k=ishl1,ishl2
        
        m=m+1
        
c     indices of atoms involved
        
        i=listshl(m,2)
        j=listshl(m,3)

c     calculate atom translational kinetic energy
        
        ppp=((weight(i)*vxx(i)+weight(j)*vxx(j))**2
     x      +(weight(i)*vyy(i)+weight(j)*vyy(j))**2
     x      +(weight(i)*vzz(i)+weight(j)*vzz(j))**2)
     x      /(weight(i)+weight(j))

c     calculate individual core and shell kinetic energies
        
        ccc=weight(i)*(vxx(i)**2+vyy(i)**2+vzz(i)**2)
        sss=weight(j)*(vxx(j)**2+vyy(j)**2+vzz(j)**2)
        
c     calculate core-shell internal kinetic energy
        
        shlke=shlke+0.5d0*(ccc+sss-ppp)
        
      enddo

c     global average of core-shell internal kinetic energy
        
      if(mxnode.gt.1)then
        buffer(1)=shlke
        call gdsum(buffer(1),1,buffer(2))
        shlke=buffer(1)
      endif
      
      return
      end subroutine corshl

      subroutine put_shells_on_cores(idnode,mxnode,ntshl)
      
c***********************************************************************
c     
c     dl_poly subroutine for placing shells on top of cores in the
c     shell model at the start of a simulation
c     
c     copyright - daresbury laboratory
c     author    - w. smith feb 2006
c     
c***********************************************************************
            
      implicit none
      
      integer, allocatable :: ltop(:)
      integer  idnode,mxnode,ntshl,ishl1,ishl2,i,j,k,m,fail,mtshl
      
c     block indices
      
      ishl1=(idnode*ntshl)/mxnode+1
      ishl2=((idnode+1)*ntshl)/mxnode
      mtshl=ishl2-ishl1+1
      
c     allocate ltop array

      allocate(ltop(mtshl),stat=fail)
      
c     zero ltop array
      
      do i=1,mtshl
        ltop(i)=0
      enddo

c     loop over all specified core-shell pairs
      
      do m=1,mtshl
        
c     indices of atoms involved
        
        i=listshl(m,2)
        j=listshl(m,3)
        
c     set shell and core positions equal
        
        ltop(m)=j
        xxx(j)=xxx(i)
        yyy(j)=yyy(i)
        zzz(j)=zzz(i)
        
      enddo
      
c     merge data on different processors
      
      if(mxnode.gt.1)call merge1
     x  (idnode,mxnode,ntshl,ltop,xxx,yyy,zzz,buffer)
      
c     deallocate ltop array
      
      deallocate(ltop,stat=fail)
      
      return
      end subroutine put_shells_on_cores

      subroutine shlfrc
     x  (lsolva,lfree,lexcite,idnode,imcon,mxnode,ntshl,engshl,virshl)

c***********************************************************************
c     
c     dl_poly subroutine for calculating shell model spring energy and 
c     force terms in molecular dynamics.
c     
c     copyright - daresbury laboratory 1994
c     author    - w. smith        july 1994
c
c***********************************************************************
      
      implicit none

      logical lsolva,lfree,lexcite,lselect
      integer idnode,imcon,mxnode,ntshl,ishl1,ishl2,i,j,k,m,kk
      real(8) engshl,virshl,rij2,omega,gamma,ffx,ffy,ffz,strs

      dimension strs(6)
      
c     check adequate workspace is available

      if(mxxdf.lt.mxshl)call error(idnode,423)

c     block indices

      ishl1=(idnode*ntshl)/mxnode+1
      ishl2=((idnode+1)*ntshl)/mxnode

c     initialise accumulators
      
      engshl=0.d0
      virshl=0.d0
      shl_fre=0.d0
      shl_vir=0.d0
      
      do i=1,6
        strs(i)=0.d0
      enddo

      if(lsolva)then
        
        lcomp(5)=.true.
        shl_sol(:)=0.d0
        if(lexcite)shl_exc(:)=0.d0
        
      endif
      
c     calculate core-shell separation vectors
      
      m=0
      do k=ishl1,ishl2
        
        m=m+1

c     indices of core and shell
        
        i=listshl(m,2)
        j=listshl(m,3)
        
c     components of bond vector
        
        xdf(m)=xxx(i)-xxx(j)
        ydf(m)=yyy(i)-yyy(j)
        zdf(m)=zzz(i)-zzz(j)
        
      enddo
      
c     periodic boundary condition
      
      call images(imcon,0,1,m,cell,xdf,ydf,zdf)

c     loop over all specified core-shell units
      
      m=0
      do k=ishl1,ishl2
        
        m=m+1

c     index of potential parameters

        kk=listshl(m,1)

c     core-shell separation
        
        rij2=xdf(m)**2+ydf(m)**2+zdf(m)**2
        
c     calculate scalar constant terms
        
        omega=(0.5d0*prmshl(kk,1)+0.25d0*prmshl(kk,2)*rij2)*rij2
        gamma=prmshl(kk,1)+prmshl(kk,2)*rij2

c     indices of core and shell
        
        i=listshl(m,2)
        j=listshl(m,3)
        
c     set selection control
        
        lselect=.true.
        
        if(lexcite)then
          
c     selected excitation option
        
          if((atm_fre(i).ne.1).and.(atm_fre(j).ne.1))then
            
c     reset selection control
            
            lselect=(atm_fre(i)+atm_fre(i).eq.0)
            
            if(lsolva)then
              shl_exc(atmolt(i))=shl_exc(atmolt(i))+omega
            endif
            
          endif
          
        elseif(lfree)then
          
c     selected free energy option
          
          if((atm_fre(i).eq.1).or.(atm_fre(j).eq.1))then
            
c     set hamiltonian mixing parameter

            shl_fre=shl_fre-omega
            shl_vir=shl_vir-gamma*rij2
            omega=lambda1*omega
            gamma=lambda1*gamma
            
          elseif((atm_fre(i).eq.2).or.(atm_fre(j).eq.2))then
            
c     set hamiltonian mixing parameter

            shl_fre=shl_fre+omega
            shl_vir=shl_vir+gamma*rij2
            omega=lambda2*omega
            gamma=lambda2*gamma
                        
          endif
          
        endif
        
        if(lselect)then
          
c     calculate spring energy and virial
          
          engshl=engshl+omega
          virshl=virshl+gamma*rij2
          
          if(lsolva)then
            shl_sol(atmolt(i))=shl_sol(atmolt(i))+omega
          endif
          
c     calculate spring forces
          
          ffx=-gamma*xdf(m)
          ffy=-gamma*ydf(m) 
          ffz=-gamma*zdf(m)
          
          fxx(i)=fxx(i)+ffx
          fyy(i)=fyy(i)+ffy
          fzz(i)=fzz(i)+ffz
          
          fxx(j)=fxx(j)-ffx
          fyy(j)=fyy(j)-ffy
          fzz(j)=fzz(j)-ffz
          
c     calculate stress tensor
          
          strs(1)=strs(1)+xdf(m)*ffx
          strs(2)=strs(2)+xdf(m)*ffy
          strs(3)=strs(3)+xdf(m)*ffz
          strs(4)=strs(4)+ydf(m)*ffy
          strs(5)=strs(5)+ydf(m)*ffz
          strs(6)=strs(6)+zdf(m)*ffz
          
        endif
        
      enddo
        
c     complete stress tensor

      stress(1)=stress(1)+strs(1)
      stress(2)=stress(2)+strs(2)
      stress(3)=stress(3)+strs(3)
      stress(4)=stress(4)+strs(2)
      stress(5)=stress(5)+strs(4)
      stress(6)=stress(6)+strs(5)
      stress(7)=stress(7)+strs(3)
      stress(8)=stress(8)+strs(5)
      stress(9)=stress(9)+strs(6)

c     sum contributions to potential and virial
      
      if(mxnode.gt.1) then
        
        buffer(1)=engshl
        buffer(2)=virshl
        buffer(3)=shl_fre
        buffer(4)=shl_vir
        call gdsum(buffer(1),4,buffer(5))
        engshl=buffer(1)
        virshl=buffer(2)
        shl_fre=buffer(3)
        shl_vir=buffer(4)
        
c     sum up solvation energies
        
        if(lsolva)then

          call gdsum(shl_sol(1),mxtmls,buffer(1))
          if(lexcite)call gdsum(shl_exc(1),mxtmls,buffer(1))
          
        endif
        
      endif
      
      return
      end subroutine shlfrc

      subroutine check_shells(idnode,itmols,nshels,ngrp)

c***********************************************************************
c     
c     dl_poly subroutine to check no core-shell units are in 
c     rigid bodies
c     
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2003
c     
c***********************************************************************

      implicit none

      integer idnode,itmols,nshels,ngrp,k1,ia,ib,kk
      integer id,jj,ia1,jk,ib1

      nshels=nshels-numshl(itmols)
      do k1=1,numshl(itmols)

        nshels=nshels+1
        ia=lstshl(nshels,1)
        ib=lstshl(nshels,2)

        ngrp=ngrp-numgrp(itmols)

        do kk=1,numgrp(itmols)
          
          ngrp=ngrp+1
          id=listyp(ngrp)
          
          do jj=1,numgsit(id)-1
            
            ia1=lstgst(ngrp,jj)
            if(ia1.eq.ia) then

              do jk=jj,numgsit(id)
                
                ib1=lstgst(ngrp,jk)
                if(ib1.eq.ib) then 
                  
                  if(idnode.eq.0)write(nrite,'(/,13x,a,2i10)')
     x              'error: sites ',ia,ib
                  call error(idnode,456)

                endif

              enddo

            elseif(ia1.eq.ib) then

              do jk=jj,numgsit(id)
                
                ib1=lstgst(ngrp,jk)
                if(ib1.eq.ia) then 
                  
                  if(idnode.eq.0)write(nrite,'(/,13x,a,2i10)')
     x              'error: sites ',ia,ib
                  call error(idnode,456)

                endif

              enddo

            endif

          enddo
        enddo
      enddo

      return
      end subroutine check_shells

      subroutine relax_shells
     x   (relaxed,keyrlx,idnode,mxnode,natms,ntpmls,tstep,rlxtol)

c***********************************************************************
c     
c     dl_poly subroutine for relaxing shells to zero force
c     
c     copyright - daresbury laboratory
c     author    - w. smith    sep 2004
c     
c***********************************************************************

      implicit none

      logical relaxed,newjob
      integer keyrlx,idnode,mxnode,natms,i,j,itmols,imols
      integer iatm1,iatm2,fail,numopt,isite,ishls,jshls,lshls,ntpmls
      integer nbuff
      real(8) hnorm,grad0,grad1,grad2,stride,tstep,step
      real(8) ggg,gam2,sgn,rlxtol

      dimension fail(4)

      integer, allocatable :: lstopt(:)
      real(8), allocatable :: ggx(:),ggy(:),ggz(:)
      real(8), allocatable :: hhx(:),hhy(:),hhz(:)
      real(8), allocatable :: oxx(:),oyy(:),ozz(:)

      save hnorm,grad0,grad1,grad2,stride,lstopt
      save ggx,ggy,ggz,hhx,hhy,hhz,oxx,oyy,ozz,numopt,nbuff,sgn

      data newjob/.true./,fail/0,0,0,0/

c     define initial data

      if(newjob)then

        newjob=.false.
        allocate(lstopt(mxatms),stat=fail(1))
        allocate(ggx(mxatms),ggy(mxatms),ggz(mxatms),stat=fail(2))
        allocate(hhx(mxatms),hhy(mxatms),hhz(mxatms),stat=fail(3))
        allocate(oxx(mxatms),oyy(mxatms),ozz(mxatms),stat=fail(4))
        do i=1,4
          if(fail(i).ne.0)call error(idnode,1970)
        enddo

c     identify the shells

        isite=0
        ishls=0
        jshls=0
        do i=1,natms

          lstopt(i)=0

        enddo
        do itmols=1,ntpmls

          do imols=1,nummols(itmols)

            do lshls=1,numshl(itmols)
              
              ishls=ishls+1
              lstopt(lstshl(lshls+jshls,2)+isite)=1
              
            enddo
            
            isite=isite+numsit(itmols)

          enddo

          jshls=jshls+numshl(itmols)

        enddo

        numopt=ishls

      endif

c     load coordinates of shells

      j=0
      do i=1,natms

        if(lstopt(i).gt.0)then

          j=j+1
          oxx(j)=xxx(i)
          oyy(j)=yyy(i)
          ozz(j)=zzz(i)
          ggx(j)=fxx(i)
          ggy(j)=fyy(i)
          ggz(j)=fzz(i)

        endif

      enddo

c     step length for relaxation

      step=tstep**2

c     define atoms for this nodes

      iatm1=(idnode*numopt)/mxnode+1
      iatm2=((idnode+1)*numopt)/mxnode

      ggg=0.d0
      do i=iatm1,iatm2
        ggg=ggg+ggx(i)**2+ggy(i)**2+ggz(i)**2
      enddo
      if(mxnode.gt.1)then
        buffer(1)=ggg
        call gdsum(buffer(1),1,buffer(2))
        ggg=buffer(1)
      endif
      ggg=sqrt(ggg)

c     check convergence

      if(abs(ggg)/dble(numopt).lt.rlxtol)then

        relaxed=.true.
        return
        
      endif

      if(keyrlx.eq.0) then

c     set original search direction

        hnorm=ggg
        grad0=ggg
        grad2=ggg
        do i=iatm1,iatm2

          hhx(i)=ggx(i)
          hhy(i)=ggy(i)
          hhz(i)=ggz(i)
          oxx(i)=oxx(i)+step*hhx(i)
          oyy(i)=oyy(i)+step*hhy(i)
          ozz(i)=ozz(i)+step*hhz(i)

        enddo
        keyrlx=1
        sgn=1.d0

      elseif(keyrlx.eq.1)then

c     line search along chosen direction

        grad1=grad2
        grad2=0.d0
        do i=iatm1,iatm2
          grad2=grad2+hhx(i)*ggx(i)+hhy(i)*ggy(i)+hhz(i)*ggz(i)
        enddo
        if(mxnode.gt.1)then
          buffer(1)=grad2
          call gdsum(buffer(1),1,buffer(2))
          grad2=buffer(1)
        endif
        grad2=sgn*grad2/hnorm

c     linear extrapolation to minimum

        stride=sgn*step
        if(grad2.lt.0.d0)then

          keyrlx=2
          stride=sgn*step*grad2/(grad1-grad2)

        endif
        
        do i=iatm1,iatm2
          
          oxx(i)=oxx(i)+stride*hhx(i)
          oyy(i)=oyy(i)+stride*hhy(i)
          ozz(i)=ozz(i)+stride*hhz(i)
          
        enddo

      elseif(keyrlx.eq.2)then

c     construct conjugate search vector

        gam2=(ggg/grad0)**2
        hnorm=0.d0
        grad0=ggg
        grad2=0.d0
        do i=iatm1,iatm2
          
          hhx(i)=ggx(i)+gam2*hhx(i)
          hhy(i)=ggy(i)+gam2*hhy(i)
          hhz(i)=ggz(i)+gam2*hhz(i)
          hnorm=hnorm+hhx(i)**2+hhy(i)**2+hhz(i)**2
          grad2=grad2+hhx(i)*ggx(i)+hhy(i)*ggy(i)+hhz(i)*ggz(i)
          
        enddo
        if(mxnode.gt.1)then
          
          buffer(1)=hnorm
          buffer(2)=grad2
          call gdsum(buffer(1),2,buffer(3))
          hnorm=buffer(1)
          grad2=buffer(2)
          
        endif
        hnorm=sqrt(hnorm)
        grad2=grad2/hnorm
        sgn=sign(1.d0,grad2)
        grad2=sgn*grad2

        do i=iatm1,iatm2
          
          oxx(i)=oxx(i)+sgn*step*hhx(i)
          oyy(i)=oyy(i)+sgn*step*hhy(i)
          ozz(i)=ozz(i)+sgn*step*hhz(i)
          
        enddo
        
        keyrlx=1

      endif 
      
c     merge coordinate arrays
      
      if(mxnode.gt.1)then
        
        nbuff=6*(numopt+mxnode-1)/mxnode
        call merge(idnode,mxnode,numopt,nbuff,oxx,oyy,ozz,buffer)
        
      endif

c     unload coordinates of shells

      j=0
      do i=1,natms

        if(lstopt(i).gt.0)then

          j=j+1
          xxx(i)=oxx(j)
          yyy(i)=oyy(j)
          zzz(i)=ozz(j)

        endif

      enddo
      
      return
      end subroutine relax_shells
      
      end module core_shell_module
