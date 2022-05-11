      module metal_module

c***********************************************************************
c     
c     dl_poly module for defining metal potential arrays
c     copyright - daresbury laboratory
c     author    - w. smith    sep 2003
c     
c***********************************************************************

      use config_module
      use error_module
      use pair_module
      use parse_module
      use property_module
      use setup_module
      use site_module
      use utility_module
      
      implicit none

      logical lmetab
      integer, allocatable :: ltpmet(:),lstmet(:)
      real(8), allocatable :: prmmet(:,:),vmet(:,:,:),dmet(:,:,:)
      real(8), allocatable :: rho(:),elrcm(:),vlrcm(:),fmet(:,:,:)

      save lmetab,ltpmet,lstmet,prmmet,vmet,dmet,fmet,rho,elrcm,vlrcm

      contains
      
      subroutine alloc_met_arrays(idnode,mxnode)
      
      implicit none
      
      logical safe
      integer, parameter :: nnn=8
      integer i,fail,idnode,mxnode
      
      dimension fail(nnn)
      
      safe=.true.
      
c     allocate arrays
      
      fail(:)=0
      
      allocate (ltpmet(mxmet),stat=fail(1))
      allocate (lstmet(mxmet),stat=fail(2))
      allocate (prmmet(mxmet,mxpmet),stat=fail(3))
      allocate (vmet(mxgrid,mxmet,2),stat=fail(4))
      allocate (dmet(mxgrid,mxmet,2),stat=fail(5))
      allocate (rho(mxatms),stat=fail(6))
      allocate (elrcm(0:mxsmet),stat=fail(7))
      allocate (vlrcm(0:mxsmet),stat=fail(8))
      
      if(any(fail.gt.0))safe=.false.      
      if(mxnode.gt.1)call gstate(safe)    
      if(.not.safe)call error(idnode,1280)
      
      end subroutine alloc_met_arrays

      subroutine define_metals
     x   (safe,lunits,lmols,idnode,ntpmet,ntpatm,rmet,engunit)

c***********************************************************************
c     
c     dl_poly subroutine for defining metal potentials
c     
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2003
c     amended   - w. smith  march 2006
c     
c***********************************************************************

      implicit none

      logical safe,lunits,lmols
      character*8 keyword
      character*8 atom1,atom2
      character*1 message(80)
      integer idnode,ntpmet,ntpatm,idum,imet,j
      integer keypot,numpar,katm1,katm2,keymet,ntab,i,fail,itpmet
      integer jtpatm
      real(8) rmet,engunit

      real(8), allocatable :: parpot(:)
      allocate (parpot(mxpmet),stat=fail)

      ntpmet=intstr(record,lenrec,idum)
      
      lmetab=.false.
      
      if(idnode.eq.0) then
        
        write(nrite,"(/,/,1x,'number of specified metal ',
     x    'potentials',i10)") ntpmet
        write(nrite,"(/,/,16x,'atom 1  ','atom 2  ',3x,
     x    ' key',30x,'parameters'/,/)")
        
      endif      

      if(ntpmet.ge.mxmet) call error(idnode,71)
      if(.not.lunits) call error(idnode,6)
      if(.not.lmols) call error(idnode,13)
      
      do imet=1,mxmet
        lstmet(imet)=0
        ltpmet(imet)=0
      enddo
      
      do itpmet=1,ntpmet
        
        do i=1,mxpmet
          parpot(i)=0.d0
        enddo
        
        call getrec(safe,idnode,nfield)
        if(.not.safe)return

        call copystring(record,message,80)
        call getword(atom1,record,8,lenrec)
        call getword(atom2,record,8,lenrec)
        call lowcase(record,4)
        call getword(keyword,record,4,lenrec)
        
        if(keyword(1:4).eq.'stch') then
          keypot=1
          numpar=5
        else if(keyword(1:4).eq.'fnsc') then
          keypot=2
          numpar=7
        else if(keyword(1:4).eq.'gupt') then
          keypot=3
          numpar=5
        elseif(keyword(1:4).eq.'eam ') then
          keypot=999
          lmetab=.true.
          numpar=0
        else
          if(idnode.eq.0) write(nrite,*) message
          call error(idnode,461)
        endif
        
        if(.not.lmetab)then
          
          parpot(1)=dblstr(record,lenrec,idum)
          parpot(2)=dblstr(record,lenrec,idum)
          parpot(3)=dblstr(record,lenrec,idum)
          parpot(4)=dblstr(record,lenrec,idum)
          parpot(5)=dblstr(record,lenrec,idum)
          parpot(6)=dblstr(record,lenrec,idum)
          parpot(7)=dblstr(record,lenrec,idum)
          
          if(idnode.eq.0)
     x      write(nrite,"(16x,2a8,2x,a4,3x,1p,9e13.5)") 
     x      atom1,atom2,keyword(1:4),(parpot(j),j=1,numpar)
        
        endif
        
        katm1=0
        katm2=0
        
        do jtpatm=1,ntpatm

          if(atom1.eq.unqatm(jtpatm))katm1=jtpatm
          if(atom2.eq.unqatm(jtpatm))katm2=jtpatm
          
        enddo
        
        if(katm1.eq.0.or.katm2.eq.0) then
          call  error(idnode,463)
        endif
        
        keymet=loc2(katm1,katm2)
        
c     convert energies to internal unit

        if(keymet.ge.mxmet) call error(idnode,465)
        
        parpot(1)=parpot(1)*engunit

        if(keypot.eq.2)then

          parpot(2)=parpot(2)*engunit
          parpot(3)=parpot(3)*engunit
          parpot(5)=parpot(5)*engunit

        endif

        if(keypot.eq.3)then
          parpot(4)=parpot(4)*engunit
        endif

        if(lstmet(keymet).ne.0) call error(idnode,141)
        lstmet(keymet)=itpmet
        ltpmet(itpmet)=keypot
        if(itpmet.gt.1)then
          if(keypot.ne.ltpmet(itpmet-1))call error(idnode,72)
        endif        
        
        if(.not.lmetab)then
          
          do i=1,numpar
            prmmet(itpmet,i)=parpot(i)
          enddo
          
        endif
        
      enddo

c     check for unspecified atom-atom potentials
      
      ntab=(ntpatm*(ntpatm+1))/2
      
      if(ntpmet.lt.ntab) then
        
        call warning(idnode,110,0.d0,0.d0,0.d0)

        do i=1,ntab
          if(lstmet(i).eq.0) lstmet(i)=ntpmet+1
        enddo

c     set zero potential for undefined interactions
        
        do i=1,mxmet
          
          vmet(1,i,1)=0.d0
          vmet(1,i,2)=0.d0
          dmet(1,i,1)=0.d0
          dmet(1,i,2)=0.d0
          
        enddo
        
      endif

c     generate metal force arrays
      
      call metgen(idnode,ntpatm,rmet)
      
      if(lmetab)
     x  call mettab(ntpatm,idnode,rmet,engunit)
      
      deallocate (parpot,stat=fail)

      return
      end subroutine define_metals

      subroutine metdens
     x  (idnode,imcon,mxnode,natms,engmet,virden)
      
c***********************************************************************
c     
c     dl_poly subroutine for calculating local density in metals
c     using the verlet neighbour list and sutton-chen potentials
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1995
c     author    - w. smith june 1995
c     amended   - w. smith  march 2006
c     
c***********************************************************************
      
      implicit none
      
      logical safe
      integer idnode,imcon,mxnode,natms,i,j,k,ii,k0,l
      real(8) engmet,engtmp,virden,rhosqr,rrr,ppp,fk0,fk1,fk2,t1,t2
      
      safe=.true.
      
c     initialise energy accumulator
      
      engmet=0.d0
      virden=0.d0

c     initialise density array
      
      do i=1,natms
        rho(i)=0.d0
      enddo

c     calculate local atomic density
      
      ii=0
      
c     outer loop over atoms

      do i=idnode+1,natms,mxnode
        
        ii=ii+1
        
c     calculate interatomic distances
        
        do k=1,lentry(ii)
          
          j=list(ii,k)
          ilist(k)=j
          xdf(k)=xxx(i)-xxx(j)
          ydf(k)=yyy(i)-yyy(j)
          zdf(k)=zzz(i)-zzz(j)
          
        enddo
        
c     periodic boundary conditions
        
        call images(imcon,0,1,lentry(ii),cell,xdf,ydf,zdf)
        
c     square of distances
        
        do k=1,lentry(ii)
          rsqdf(k)=xdf(k)**2+ydf(k)**2+zdf(k)**2
        enddo
        
        if(lmetab)then
          
c     calculate EAM density contributions
          
          call eamden(safe,i,lentry(ii))
          
        else
          
c     calculate FS density contributions
          
          call fsden(safe,i,lentry(ii))
          
        endif
        
      enddo

c     global sum of local atomic densities
      
      if(mxnode.gt.1)call gdsum(rho,natms,buffer)
      
c     calculate embedding energy

      if(lmetab)then
        
c     density terms for eam potentials
        
        do i=1,natms
          
          if(rho(i).gt.0.d0)then
            
            k0=ltype(i)
            rrr=rho(i)-fmet(2,k0,1)
            l=min(nint(rrr/fmet(4,k0,1)),int(fmet(1,k0,1))-1)
            if(l.lt.1)then
              
              safe=.false.
              
            else
              
              ppp=(rrr/fmet(4,k0,1))-dble(l)
              
c     calculate embedding energy using 3-point interpolation
              
              fk0=fmet(l-1,k0,1)
              fk1=fmet(l,k0,1)
              fk2=fmet(l+1,k0,1)
              
              t1=fk1+(fk1-fk0)*ppp
              t2=fk1+(fk2-fk1)*ppp
              if(ppp.lt.0.d0)then
                engtmp=-(t1+0.5d0*(t2-t1)*(ppp+1.d0))
              else
                engtmp=-(t2+0.5d0*(t2-t1)*(ppp-1.d0))
              endif

              engmet=engmet+engtmp

c     calculate derivative of embedding function wrt density using 3-point
c     interpolation - store result in rho array
              
              fk0=fmet(l-1,k0,2)
              fk1=fmet(l,k0,2)
              fk2=fmet(l+1,k0,2)
              
              t1=fk1+(fk1-fk0)*ppp
              t2=fk1+(fk2-fk1)*ppp
              if(ppp.lt.0.d0)then
                rho(i)=(t1+0.5d0*(t2-t1)*(ppp+1.d0))
              else
                rho(i)=(t2+0.5d0*(t2-t1)*(ppp-1.d0))
              endif
            
            endif
            
          endif
          
        enddo
        
      else
        
c     analytical square root of density dependence
        
        do i=1,natms
          
          if(rho(i).gt.0.d0)then
            
            rhosqr=sqrt(rho(i)+elrcm(ltype(i)))
            engmet=engmet+rhosqr
            rho(i)=0.5d0/rhosqr
            virden=virden+vlrcm(ltype(i))/rhosqr

          endif
          
        enddo
        
      endif
      
      engmet=-engmet/dble(mxnode)
      virden=virden/dble(mxnode)
      
c     check interpolation is safe
      
      if(mxnode.gt.1)call gstate(safe)
      if(.not.safe)call error(idnode,142)
      
      return
      end subroutine metdens

      subroutine fsden(safe,iatm,ik)

c***********************************************************************
c     
c     dl_poly subroutine for calculating local atomic density
c     for FS type metal potentials
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory
c     author    - w. smith  june  1995
c     amended   - w. smith  march 2006
c     
c***********************************************************************
      
      implicit none
      
      logical safe
      integer iatm,jatm,ik,m,k0,l
      real(8) rdr,ai,aj,ab,rsq,rrr,ppp,t1,t2
      real(8) vk0,vk1,vk2,density
      
c     start of primary loop for density

      ai=dble(ltype(iatm))

      do m=1,ik

c     atomic and potential function indices
        
        jatm=ilist(m)

        aj=dble(ltype(jatm))
        if(ai.gt.aj) then
          ab=ai*(ai-1.d0)*0.5d0+aj+0.5d0
        else
          ab=aj*(aj-1.d0)*0.5d0+ai+0.5d0
        endif

        k0=lstmet(int(ab))
        
        if((ltpmet(k0).ge.1).and.(abs(dmet(1,k0,1)).gt.0.d0))then

c     apply truncation of potential
          
          rsq=rsqdf(m)

c     apply cutoff condition
      
          if(rsq.le.dmet(3,k0,1)**2)then
            
c     interpolation parameters

            rdr=1.d0/dmet(4,k0,1)
            rrr=sqrt(rsq)-dmet(2,k0,1)
            l=min(nint(rrr*rdr),int(dmet(1,k0,1))-1)
            if(l.lt.1)then
              
              safe=.false.
              
            else
              
              ppp=rrr*rdr-dble(l)
              
c     calculate density using 3-point interpolation
              
              vk0=dmet(l-1,k0,1)
              vk1=dmet(l,k0,1)
              vk2=dmet(l+1,k0,1)
              
              t1=vk1+ppp*(vk1-vk0)
              t2=vk1+ppp*(vk2-vk1)
              if(ppp.lt.0.d0)then
                density=t1+0.5d0*(t2-t1)*(ppp+1.d0)
              else
                density=t2+0.5d0*(t2-t1)*(ppp-1.d0)
              endif
              
              if(ai.gt.aj)then
                
                rho(iatm)=rho(iatm)+density*dmet(1,k0,2)
                rho(jatm)=rho(jatm)+density*dmet(2,k0,2)
                
              else
                
                rho(iatm)=rho(iatm)+density*dmet(2,k0,2)
                rho(jatm)=rho(jatm)+density*dmet(1,k0,2)
                
              endif
              
            endif
            
          endif
        
        endif
        
      enddo
      
      return
      end subroutine fsden

      subroutine eamden(safe,iatm,ik)

c***********************************************************************
c     
c     dl_poly subroutine for calculating local atomic density
c     for EAM type metal potentials
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory
c     author    - w. smith  june  1995
c     amended   - w. smith  march 2006
c     
c***********************************************************************
            
      implicit none
      logical safe
      integer iatm,jatm,ik,m,l,ktyp1,ktyp2
      real(8) rdr,rsq,rrr,ppp,t1,t2
      real(8) vk0,vk1,vk2,density
      
c     start of primary loop for density

      ktyp1=ltype(iatm)

      do m=1,ik

c     first metal atom density
        
        jatm=ilist(m)
        ktyp2=ltype(jatm)
        
        if(abs(dmet(1,ktyp2,1)).gt.0.d0)then
          
c     apply truncation of potential
          
          rsq=rsqdf(m)

          if(rsq.le.dmet(3,ktyp2,1)**2)then
            
c     interpolation parameters
            
            rdr=1.d0/dmet(4,ktyp2,1)
            rrr=sqrt(rsq)-dmet(2,ktyp2,1)
            l=min(nint(rrr*rdr),int(dmet(1,ktyp2,1))-1)
            if(l.lt.1)then
              
              safe=.false.
              
            else
              
              ppp=rrr*rdr-dble(l)
              
c     calculate density using 3-point interpolation
              
              vk0=dmet(l-1,ktyp2,1)
              vk1=dmet(l,ktyp2,1)
              vk2=dmet(l+1,ktyp2,1)
              
              t1=vk1+ppp*(vk1-vk0)
              t2=vk1+ppp*(vk2-vk1)
              if(ppp.lt.0.d0)then
                density=t1+0.5d0*(t2-t1)*(ppp+1.d0)
              else
                density=t2+0.5d0*(t2-t1)*(ppp-1.d0)
              endif
              
              rho(iatm)=rho(iatm)+density
              if(ktyp1.eq.ktyp2)rho(jatm)=rho(jatm)+density
              
            endif
            
          endif

        endif
        
c     second  metal atom density
            
        if(ktyp1.ne.ktyp2)then
          
          if(abs(dmet(1,ktyp1,1)).gt.0.d0)then
            
c     apply truncation of potential
            
            if(rsq.le.dmet(3,ktyp1,1)**2)then
              
c     interpolation parameters
              
              rdr=1.d0/dmet(4,ktyp1,1)
              rrr=sqrt(rsq)-dmet(2,ktyp1,1)
              l=min(nint(rrr*rdr),int(dmet(1,ktyp1,1))-1)
              if(l.lt.1)then
                
                safe=.false.
                
              else
                
                ppp=rrr*rdr-dble(l)
                
c     calculate density using 3-point interpolation
                
                vk0=dmet(l-1,ktyp1,1)
                vk1=dmet(l,ktyp1,1)
                vk2=dmet(l+1,ktyp1,1)
                
                t1=vk1+(vk1-vk0)*ppp
                t2=vk1+(vk2-vk1)*ppp
                if(ppp.lt.0.d0)then
                  density=t1+0.5d0*(t2-t1)*(ppp+1.d0)
                else
                  density=t2+0.5d0*(t2-t1)*(ppp-1.d0)
                endif
                
                rho(jatm)=rho(jatm)+density
                
              endif
              
            endif
            
          endif
        
        endif
        
      enddo
      
      return
      end subroutine eamden

      subroutine metfrc(safe,iatm,ik,engmet,virmet)

c***********************************************************************
c     
c     dl_poly subroutine for calculating metal forces
c     for EAM and FS potentials using a verlet neighbour list
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory
c     author    - w. smith   june 1995
c     amended   - w. smith  march 2006
c     
c***********************************************************************
      
      implicit none
      
      logical safe
      integer iatm,jatm,ik,m,k0,l,ld,ktyp1,ktyp2
      real(8) engmet,virmet,strs
      real(8) rdr,rsq,rrr,ppp,vk0,vk1,vk2,t1,t2,gk0,gk1,gk2
      real(8) gamma,gamma1,gamma2,gamma3,fx,fy,fz,fi
      dimension fi(3),strs(6)

CDIR$ CACHE_ALIGN fi

c     initialise potential energy and virial
      
      engmet=0.d0
      virmet=0.d0

c     initialise stress tensor accumulators
      
      strs(:)=0.d0

c     store forces for iatm 
      
      fi(1)=fxx(iatm)
      fi(2)=fyy(iatm)
      fi(3)=fzz(iatm)
      ktyp1=ltype(iatm)
      
c     start of primary loop for forces evaluation
      
      do m=1,ik
        
c     atomic and potential function indices
        
        jatm=ilist(m)
        ktyp2=ltype(jatm)
        k0=lstmet(loc2(ktyp1,ktyp2))
        
        if((ltpmet(k0).gt.0).and.(abs(vmet(1,k0,1)).gt.0.d0))then

c     apply truncation of potential
          
          rsq=rsqdf(m)

          if(rsq.le.vmet(3,k0,1)**2)then
            
c     interpolation parameters
      
            rdr=1.d0/vmet(4,k0,1)
            rrr=sqrt(rsq)-vmet(2,k0,1)
            l=min(nint(rrr*rdr),int(vmet(1,k0,1))-1)
            if(l.lt.1)then
              
              safe=.false.
              gamma1=0.d0
              
            else
              
              ppp=rrr*rdr-dble(l)
              
c     calculate interaction energy using 3-point interpolation
              
              vk0=vmet(l-1,k0,1)
              vk1=vmet(l,k0,1)
              vk2=vmet(l+1,k0,1)
              
              t1=vk1+(vk1-vk0)*ppp
              t2=vk1+(vk2-vk1)*ppp
              if(ppp.lt.0.d0)then
                engmet=engmet+t1+0.5d0*(t2-t1)*(ppp+1.d0)
              else
                engmet=engmet+t2+0.5d0*(t2-t1)*(ppp-1.d0)
              endif
              
c     calculate pair forces using 3-point interpolation
              
              gk0=vmet(l-1,k0,2)
              gk1=vmet(l,k0,2)
              gk2=vmet(l+1,k0,2)
              
              t1=gk1+(gk1-gk0)*ppp
              t2=gk1+(gk2-gk1)*ppp
              if(ppp.lt.0.d0)then
                gamma1=t1+0.5d0*(t2-t1)*(ppp+1.d0)
              else
                gamma1=t2+0.5d0*(t2-t1)*(ppp-1.d0)
              endif
              
            endif
            
c     calculate embedding forces using 3-point interpolation

            if(lmetab)then
              
              if(rsq.le.dmet(3,ktyp2,1)**2)then
                
                rdr=1.d0/dmet(4,ktyp2,1)
                rrr=sqrt(rsq)-dmet(2,ktyp2,1)
                ld=min(nint(rrr*rdr),int(dmet(1,ktyp2,1))-1)
                if(ld.lt.1)then
                  
                  safe=.false.
                  gamma2=0.d0
                  
                else
                  
                  ppp=rrr*rdr-dble(ld)
                  
                  gk0=dmet(ld-1,ktyp2,2)
                  gk1=dmet(ld,ktyp2,2)
                  gk2=dmet(ld+1,ktyp2,2)
                  
                  t1=gk1+(gk1-gk0)*ppp
                  t2=gk1+(gk2-gk1)*ppp
                  if(ppp.lt.0.d0)then
                    gamma2=t1+0.5d0*(t2-t1)*(ppp+1.d0)
                  else
                    gamma2=t2+0.5d0*(t2-t1)*(ppp-1.d0)
                  endif
                
                endif
                
              else
                
                gamma2=0.d0
                
              endif
              
              if(ktyp1.eq.ktyp2)then
                
                gamma3=gamma2
                
              elseif(rsq.le.dmet(3,ktyp1,1)**2)then
                
                rdr=1.d0/dmet(4,ktyp1,1)
                rrr=sqrt(rsq)-dmet(2,ktyp1,1)
                ld=min(nint(rrr*rdr),int(dmet(1,ktyp1,1))-1)
                if(ld.lt.1)then
                  
                  safe=.false.
                  gamma3=0.d0
                  
                else
                  
                  ppp=rrr*rdr-dble(ld)
                  gk0=dmet(ld-1,ktyp1,2)
                  gk1=dmet(ld,ktyp1,2)
                  gk2=dmet(ld+1,ktyp1,2)
                  
                  t1=gk1+(gk1-gk0)*ppp
                  t2=gk1+(gk2-gk1)*ppp
                  if(ppp.lt.0.d0)then
                    gamma3=t1+0.5d0*(t2-t1)*(ppp+1.d0)
                  else
                    gamma3=t2+0.5d0*(t2-t1)*(ppp-1.d0)
                  endif
                  
                endif
                
              else
                
                gamma3=0.d0
                
              endif
                
              gamma=(gamma1+(gamma2*rho(iatm)+gamma3*rho(jatm)))/rsq
              
            else
              
              if(safe.and.rsq.le.dmet(3,k0,1)**2)then
                
                gk0=dmet(l-1,k0,2)
                gk1=dmet(l,k0,2)
                gk2=dmet(l+1,k0,2)
                
                t1=gk1+(gk1-gk0)*ppp
                t2=gk1+(gk2-gk1)*ppp
                if(ppp.lt.0.d0)then
                  gamma2=t1+0.5d0*(t2-t1)*(ppp+1.d0)
                else
                  gamma2=t2+0.5d0*(t2-t1)*(ppp-1.d0)
                endif
                
              else
                
                gamma2=0.d0
                
              endif
              
              if(ktyp1.gt.ktyp2)then
                
                gamma=(gamma1-gamma2*(rho(iatm)*dmet(1,k0,2)+
     x            rho(jatm)*dmet(2,k0,2)))/rsq
                
              else
                
                gamma=(gamma1-gamma2*(rho(iatm)*dmet(2,k0,2)+
     x            rho(jatm)*dmet(1,k0,2)))/rsq
                
              endif
              
            endif
            
c     calculate forces
            
            fx=gamma*xdf(m)
            fy=gamma*ydf(m)
            fz=gamma*zdf(m)
            
            fi(1)=fi(1)+fx
            fi(2)=fi(2)+fy
            fi(3)=fi(3)+fz
            
            fxx(jatm)=fxx(jatm)-fx
            fyy(jatm)=fyy(jatm)-fy
            fzz(jatm)=fzz(jatm)-fz
            
c     calculate stress tensor
            
            strs(1)=strs(1)+xdf(m)*fx
            strs(2)=strs(2)+xdf(m)*fy
            strs(3)=strs(3)+xdf(m)*fz
            strs(4)=strs(4)+ydf(m)*fy
            strs(5)=strs(5)+ydf(m)*fz
            strs(6)=strs(6)+zdf(m)*fz
            
c     calculate virial
            
            virmet=virmet-gamma*rsq
            
          endif
          
        endif
        
      enddo

c     load temps back to fxx(iatm) etc
      
      fxx(iatm)=fi(1)
      fyy(iatm)=fi(2)
      fzz(iatm)=fi(3)

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
      
      return
      end subroutine metfrc

      subroutine metgen(idnode,ntpatm,rmet)

c***********************************************************************
c     
c     dl_poly subroutine for generating potential energy and 
c     force arrays for metal potentials
c     
c     copyright - daresbury laboratory
c     author    - w. smith june 1995
c     amended   - w. smith  march 2006
c     
c***********************************************************************
      
      implicit none

      integer i,kmet,mmet,katm1,katm2,ntpatm,imet,nmet
      integer idnode,nnn,mmm
      real(8) dlrpot,rmet,rrr,eps,sig,cc0,cc1,cc2,ccc
      real(8) ddd,bet,cut1,cut2,aaa,rr0,ppp,qqq
      
c     define grid resolution for potential arrays
      
      dlrpot=rmet/dble(mxgrid-1)

c     construct arrays for metal potentials
      
      kmet=0
      do katm1=1,ntpatm

        do katm2=1,katm1

          kmet=kmet+1
          imet=lstmet(kmet)

          if(ltpmet(imet).gt.0)then

c     store array specification parameters
            
            vmet(1,imet,1)=dble(mxgrid)
            vmet(2,imet,1)=0.d0
            vmet(3,imet,1)=rmet
            vmet(4,imet,1)=dlrpot
            
            do i=1,4
              
              vmet(i,imet,2)=vmet(i,imet,1)
              dmet(i,imet,1)=vmet(i,imet,1)
              dmet(i,imet,2)=0.d0
              
            enddo

            if(ltpmet(imet).eq.1)then
              
c     sutton-chen potentials

              eps=prmmet(imet,1)
              sig=prmmet(imet,2)
              nnn=nint(prmmet(imet,3))
              mmm=nint(prmmet(imet,4))

              do i=5,mxgrid
                
                rrr=dble(i)*dlrpot
                vmet(i,imet,1)=eps*(sig/rrr)**nnn
                vmet(i,imet,2)=dble(nnn)*eps*(sig/rrr)**nnn
                dmet(i,imet,1)=(sig/rrr)**mmm
                dmet(i,imet,2)=dble(mmm)*(sig/rrr)**mmm

              enddo
              
              if(katm1.eq.katm2)then

                dmet(1,imet,2)=(prmmet(imet,1)*prmmet(imet,5))**2
                dmet(2,imet,2)=(prmmet(imet,1)*prmmet(imet,5))**2

              else

                nmet=lstmet((katm1*(katm1+1))/2)
                mmet=lstmet((katm2*(katm2+1))/2)
                dmet(1,imet,2)=(prmmet(nmet,1)*prmmet(nmet,5))**2
                dmet(2,imet,2)=(prmmet(mmet,1)*prmmet(mmet,5))**2

              endif

            else if(ltpmet(imet).eq.2)then
              
c     finnis sinclair potentials

              cc0=prmmet(imet,1)
              cc1=prmmet(imet,2)
              cc2=prmmet(imet,3)
              ccc=prmmet(imet,4)
              ddd=prmmet(imet,6)
              bet=prmmet(imet,7)
              cut1=ccc+4.d0*dlrpot
              cut2=ddd+4.d0*dlrpot

              do i=5,mxgrid
                
                rrr=dble(i)*dlrpot
                vmet(i,imet,1)=0.d0
                vmet(i,imet,2)=0.d0
                dmet(i,imet,1)=0.d0
                dmet(i,imet,2)=0.d0
                
                if(rrr.le.cut1)then
                  
                  vmet(i,imet,1)=(cc0+cc1*rrr+cc2*rrr*rrr)*(rrr-ccc)**2
                  vmet(i,imet,2)=-rrr*(2.d0*(cc0+cc1*rrr+cc2*rrr*rrr)*
     x              (rrr-ccc)+(cc1+2.d0*cc2*rrr)*(rrr-ccc)**2)
                  
                endif

                if(rrr.le.cut2)then
                  
                  dmet(i,imet,1)=(rrr-ddd)**2+bet*(rrr-ddd)**3/ddd
                  dmet(i,imet,2)=-rrr*(2.d0*(rrr-ddd)+
     x              3.d0*bet*(rrr-ddd)**2/ddd)
                  
                endif
                
              enddo
              
              if(katm1.eq.katm2)then

                dmet(1,imet,2)=prmmet(imet,5)**2
                dmet(2,imet,2)=prmmet(imet,5)**2

              else

                nmet=lstmet((katm1*(katm1+1))/2)
                mmet=lstmet((katm2*(katm2+1))/2)
                dmet(1,imet,2)=prmmet(nmet,5)**2
                dmet(2,imet,2)=prmmet(mmet,5)**2

              endif
              
            else if(ltpmet(imet).eq.3)then
              
c     gupta potentials

              aaa=prmmet(imet,1)
              rr0=prmmet(imet,2)
              ppp=prmmet(imet,3)
              qqq=prmmet(imet,5)

              do i=5,mxgrid
                
                rrr=dble(i)*dlrpot
                vmet(i,imet,1)=aaa*exp(-ppp*(rrr-rr0)/rr0)
                vmet(i,imet,2)=vmet(i,imet,1)*rrr*ppp/rr0
                dmet(i,imet,1)=exp(-2.d0*qqq*(rrr-rr0)/rr0)
                dmet(i,imet,2)=2.d0*dmet(i,imet,1)*rrr*qqq/rr0
                
              enddo
              
              dmet(1,imet,2)=prmmet(imet,4)**2
              dmet(2,imet,2)=prmmet(imet,4)**2

            else if(.not.lmetab)then
              
              call error(idnode,151)
              
            endif

          endif
          
        enddo

      enddo
      
      return
      end subroutine metgen

      subroutine lrcmetal
     x  (idnode,imcon,natms,ntpatm,engunit,rmet,volm)
      
c*************************************************************************
c     
c     DL_POLY subroutine to evaluate long-range corrections to
c     pressure and energy in a periodic metal system.
c     
c     copyright daresbury laboratory
c     author -  w. smith   june 1995
c     amended - w. smith  march 2006
c     
c***************************************************************************
      
      implicit none
      
      logical newjob
      integer idnode,imcon,natms,ntpatm,i,ka,j
      integer kmet,k0,k1,k2
      real(8) engunit,rmet,volm,twopi,forpi,eps,sig,nnn,mmm,ccc
      real(8) elrcm0,elrcm1,elrcm2,vlrcm0,vlrcm1,vlrcm2,aaa,rr0,ppp
      real(8) zet,qqq,eee

      save newjob
      data newjob/.true./
      
      twopi=2.0d0*pi
      forpi=4.0d0*pi

c     initalise counter arrays
      
      do i=1,mxsmet
        numtyp(i)=0
      enddo

c     evaluate species populations in system
      
      do i=1,natms
        
        ka=ltype(i)
        numtyp(ka)=numtyp(ka)+1
        
      enddo
      
c     number densities
      
      do i=1,ntpatm
        dens(i)=dble(numtyp(i))/volm
      enddo
      
c     long range corrections to density, energy and pressure
      
      do i=0,mxsmet

        elrcm(i)=0.d0
        vlrcm(i)=0.d0

      enddo
      
      if(imcon.ne.0.and.imcon.ne.6) then
        
        kmet=0
        do i=1,ntpatm
          
          do j=1,i
            
            elrcm0=0.d0
            elrcm1=0.d0
            elrcm2=0.d0
            vlrcm0=0.d0
            vlrcm1=0.d0
            vlrcm2=0.d0
            
            kmet=kmet+1
            k0=lstmet(kmet)
            
            if(ltpmet(k0).eq.1) then
              
c     sutton-chen potentials
              
              eps=prmmet(k0,1)
              sig=prmmet(k0,2)
              nnn=prmmet(k0,3)
              mmm=prmmet(k0,4)
              ccc=prmmet(k0,5)

              elrcm0=eps*sig**3*(sig/rmet)**(nnn-3.d0)/(nnn-3.d0)
              vlrcm0=eps*nnn*sig**3*(sig/rmet)**(nnn-3.d0)/(nnn-3.d0)
              if(i.ne.j) then
                elrcm0=elrcm0*2.d0
                vlrcm0=vlrcm0*2.d0
              endif
              elrcm(0)=elrcm(0)+twopi*volm*dens(i)*dens(j)*elrcm0
              vlrcm(0)=vlrcm(0)-twopi*volm*dens(i)*dens(j)*vlrcm0
              
              if(i.eq.j) then
                
                elrcm1=sig**3*(sig/rmet)**(mmm-3.d0)/(mmm-3.d0)*
     x            (eps*ccc)**2
                elrcm(i)=elrcm(i)+forpi*dens(i)*elrcm1
                
                vlrcm1=mmm*sig**3*(sig/rmet)**(mmm-3.d0)/(mmm-3.d0)*
     x            (eps*ccc)**2
                vlrcm(i)=vlrcm(i)+twopi*dens(i)*vlrcm1

              else

                k1=lstmet((i*(i+1))/2)
                k2=lstmet((j*(j+1))/2)
                elrcm1=sig**3*(sig/rmet)**(mmm-3.d0)/(mmm-3.d0)*
     x            (prmmet(k1,1)*prmmet(k1,5))**2
                elrcm2=sig**3*(sig/rmet)**(mmm-3.d0)/(mmm-3.d0)*
     x            (prmmet(k2,1)*prmmet(k2,5))**2
                elrcm(i)=elrcm(i)+forpi*dens(j)*elrcm1
                elrcm(j)=elrcm(j)+forpi*dens(i)*elrcm2

                vlrcm1=mmm*sig**3*(sig/rmet)**(mmm-3.d0)/(mmm-3.d0)*
     x            (prmmet(k1,1)*prmmet(k1,5))**2
                vlrcm2=mmm*sig**3*(sig/rmet)**(mmm-3.d0)/(mmm-3.d0)*
     x            (prmmet(k2,1)*prmmet(k2,5))**2
                vlrcm(i)=vlrcm(i)+twopi*dens(j)*vlrcm1
                vlrcm(j)=vlrcm(j)+twopi*dens(i)*vlrcm2

              endif
              
            else if(ltpmet(k0).eq.3) then
              
c     gupta potentials
              
              aaa=prmmet(k0,1)
              rr0=prmmet(k0,2)
              ppp=prmmet(k0,3)
              zet=prmmet(k0,4)
              qqq=prmmet(k0,5)
              eee=exp(-ppp*(rmet-rr0)/rr0)

              elrcm0=aaa*(rr0/ppp)*(rmet**2+2.d0*rmet*(rr0/ppp)+
     x          2.d0*(rr0/ppp)**2)*eee
              vlrcm0=aaa*(rmet**3+3.d0*rmet**2*(rr0/ppp)+
     x          6.d0*rmet*(rr0/ppp)**2+6.d0*(rr0/rmet)**3)*eee
              if(i.ne.j) then
                elrcm0=elrcm0*2.d0
                vlrcm0=vlrcm0*2.d0
              endif
              elrcm(0)=elrcm(0)+twopi*volm*dens(i)*dens(j)*elrcm0
              vlrcm(0)=vlrcm(0)-twopi*volm*dens(i)*dens(j)*vlrcm0
              
              eee=exp(-2.d0*qqq*(rmet-rr0)/rr0)

              if(i.eq.j) then
                
                elrcm1=(rmet**2+2.d0*rmet*(0.5d0*rr0/qqq)+
     x            2.d0*(0.5d0*rr0/qqq)**2)*(0.5d0*rr0/qqq)*eee*zet**2
                elrcm(i)=elrcm(i)+forpi*dens(i)*elrcm1

                vlrcm1=(rmet**3+3.d0*rmet**2*(0.5d0*rr0/qqq)+
     x            6.d0*rmet*(0.5d0*rr0/qqq)**2+(0.5d0*rr0/qqq)**3)*
     x            eee*zet**2
                vlrcm(i)=vlrcm(i)+twopi*dens(i)*vlrcm1

              else

                elrcm1=(rmet**2+2.d0*rmet*(0.5d0*rr0/qqq)+
     x            2.d0*(0.5d0*rr0/qqq)**2)*(0.5d0*rr0/qqq)*eee*
     x            zet**2
                elrcm2=(rmet**2+2.d0*rmet*(0.5d0*rr0/qqq)+
     x            2.d0*(0.5d0*rr0/qqq)**2)*(0.5d0*rr0/qqq)*eee*
     x            zet**2
                elrcm(i)=elrcm(i)+forpi*dens(j)*elrcm1
                elrcm(j)=elrcm(j)+forpi*dens(i)*elrcm2

                vlrcm1=(rmet**3+3.d0*rmet**2*(0.5d0*rr0/qqq)+
     x            6.d0*rmet*(0.5d0*rr0/qqq)**2+(0.5d0*rr0/qqq)**3)*
     x            eee*zet**2
                vlrcm2=(rmet**3+3.d0*rmet**2*(0.5d0*rr0/qqq)+
     x            6.d0*rmet*(0.5d0*rr0/qqq)**2+(0.5d0*rr0/qqq)**3)*
     x            eee*zet**2
                vlrcm(i)=vlrcm(i)+twopi*dens(j)*vlrcm1
                vlrcm(j)=vlrcm(j)+twopi*dens(i)*vlrcm2

              endif
              
            endif
            
          enddo
          
        enddo
        
      endif
      
      if(newjob)then
        
        newjob=.false.
        
        if(idnode.eq.0)then
          
          write(nrite,"(/,/,
     x      'long range corrections for metal potentials',/)")
          write(nrite,
     x      "('short range energy and virial corrections:',
     x      1p,2e15.6,/)")
     x    elrcm(0)/engunit,vlrcm(0)/engunit
          write(nrite,
     x      "('density dependent energy and virial corrections',/)")
          
          do i=1,ntpatm
            
            kmet=lstmet((i*(i+1))/2)
            if(lstmet(kmet).gt.0)then
              
              write(nrite,"(25x,a8,1p,2e15.6)")unqatm(i),
     x          elrcm(i)/engunit,vlrcm(i)/engunit
              
            endif
            
          enddo
          
        endif
        
      endif
      
      return
      end subroutine lrcmetal

      subroutine mettab(ntpatm,idnode,rmet,engunit)

c***********************************************************************
c     
c     dl_poly subroutine for reading potential energy and 
c     force arrays for EAM metal forces only
c     
c     copyright - daresbury laboratory
c     author    - w. smith march 2006
c     
c***********************************************************************
      
      implicit none

      logical safe
      character*8 atom1,atom2,type
      integer idnode,ntpatm,idum,fail
      integer imet,katm1,katm2,jtpatm,i,j,k,ktype
      integer numpot,numpts,ipot
      real(8) rmet,engunit,start,finish
      data fail/0/
      
c     allocate embedding array
      
      allocate(fmet(mxgrid,mxmet,2),stat=fail)
      if(fail.ne.0)call error(idnode,36)
      
c     define zero function for undefined interactions
      
      do i=1,mxmet
        
        fmet(1,i,1)=0.d0
        fmet(1,i,2)=0.d0
        
      enddo
      
      if(idnode.eq.0)open (ntable,file='TABEAM')
      
c     skip header record
      
      call getrec(safe,idnode,ntable)
      if(.not.safe)call abort_eamtable_read(idnode,ntable)
      
c     read number of potential functions in file
      
      call getrec(safe,idnode,ntable)
      if(.not.safe)call abort_eamtable_read(idnode,ntable)
      numpot=intstr(record,lenrec,idum)
      
      do ipot=1,numpot
        
c     read data type, atom labels, number of points, start and end
        
        call getrec(safe,idnode,ntable)
        if(.not.safe)call abort_eamtable_read(idnode,ntable)
        
c     indentify data type
        
        ktype=1
        if(findstring('dens',record,idum).or.
     x    findstring('DENS',record,idum))ktype=2
        if(findstring('embe',record,idum).or.
     x    findstring('EMBE',record,idum))ktype=3
        call getword(type,record,8,lenrec)
        
c     identify atom types
        
        call getword(atom1,record,8,lenrec)
        if(ktype.eq.1)then
          call getword(atom2,record,8,lenrec)
        else
          atom2=atom1
        endif
        
c     data specifiers
        
        numpts=intstr(record,lenrec,idum)
        start=dblstr(record,lenrec,idum)
        finish=dblstr(record,lenrec,idum)
        
c     check atom indentities
        
        katm1=0
        katm2=0
        
        do jtpatm=1,ntpatm
          
          if(atom1.eq.unqatm(jtpatm))katm1=jtpatm
          if(atom2.eq.unqatm(jtpatm))katm2=jtpatm
          
        enddo
        
        if(katm1.eq.0.or.katm2.eq.0) then
          if(idnode.eq.0) 
     x      write(nrite,'(a)') ' **** '//atom1//' *** '//atom2//' ****'
          call  error(idnode,81)
        endif
        
c     check array dimensions
        
        if(mxbuff.lt.numpts+4)then
          
          if(idnode.eq.0)
     x      write(nrite,*) 'mxbuff must be >=',numpts+4,' in mettab'
          call error(idnode,28)
          
        endif
        
c     store working parameters (start shifted for DL_POLY interpolation)
        
        buffer(1)=dble(numpts+4)
        buffer(4)=(finish-start)/dble(numpts-1)
        buffer(2)=start-5.d0*buffer(4)
        buffer(3)=finish
        if(idnode.eq.0)
     x    write(nrite,"(16x,2a8,2x,a4,3x,1p,4e13.5)") 
     x    atom1,atom2,type,dble(numpts),start,finish,buffer(4)

c     read potential arrays
        
        k=4
        do j=1,(numpts+3)/4
          
          call getrec(safe,idnode,ntable)
          if(.not.safe)call abort_eamtable_read(idnode,ntable)
          buffer(k+1)=dblstr(record,lenrec,idum)
          buffer(k+2)=dblstr(record,lenrec,idum)
          buffer(k+3)=dblstr(record,lenrec,idum)
          buffer(k+4)=dblstr(record,lenrec,idum)
          k=k+4
          
        enddo
        
c     copy data to internal arrays
        
        if(ktype.eq.1)then
          
c     check range against specified cutoff
        
          if(rmet.lt.finish)call error(idnode,26)
        
c     identify potential

          imet=lstmet(loc2(katm1,katm2))
        
c     pair potential terms
          
          vmet(1,imet,1)=buffer(1)
          vmet(2,imet,1)=buffer(2)
          vmet(3,imet,1)=buffer(3)
          vmet(4,imet,1)=buffer(4)
          
          do i=5,mxgrid
            
            if(i-4.gt.numpts)then
              vmet(i,imet,1)=0.d0
            else
              vmet(i,imet,1)=buffer(i)*engunit
              buffer(i)=buffer(i)*engunit
            endif
            
          enddo
          
c     calculate derivative of pair potential function
          
          call metal_deriv(imet,vmet,buffer)
          
c     adapt derivatives for use in interpolation
      
          do i=5,numpts+4
            vmet(i,imet,2)=-(dble(i)*buffer(4)+buffer(2))*
     x        vmet(i,imet,2)
          enddo          
          
        else if(ktype.eq.2)then
          
c     check range against specified cutoff
        
          if(rmet.lt.finish)call error(idnode,26)
        
c     density  terms
          
          dmet(1,katm1,1)=buffer(1)
          dmet(2,katm1,1)=buffer(2)
          dmet(3,katm1,1)=buffer(3)
          dmet(4,katm1,1)=buffer(4)
          
          do i=5,mxgrid
            
            if(i-4.gt.numpts)then
              dmet(i,katm1,1)=0.d0
            else
              dmet(i,katm1,1)=buffer(i)
            endif
            
          enddo
          
c     calculate derivative of density function
          
          call metal_deriv(katm1,dmet,buffer)
          
c     adapt derivatives for use in interpolation
      
          dmet(1,katm1,2)=0.d0
          dmet(2,katm1,2)=0.d0
          dmet(3,katm1,2)=0.d0
          dmet(4,katm1,2)=0.d0
          do i=5,numpts+4
            dmet(i,katm1,2)=-(dble(i)*buffer(4)+buffer(2))*
     x        dmet(i,katm1,2)
          enddo          
          
        else if(ktype.eq.3)then
          
c     embedding function arrays
          
          fmet(1,katm1,1)=buffer(1)
          fmet(2,katm1,1)=buffer(2)
          fmet(3,katm1,1)=buffer(3)
          fmet(4,katm1,1)=buffer(4)
          
          do i=5,mxgrid
            
            if(i-4.gt.numpts)then
              fmet(i,katm1,1)=0.d0
            else
              fmet(i,katm1,1)=buffer(i)*engunit
              buffer(i)=buffer(i)*engunit
            endif
            
          enddo
          
c     calculate derivative of embedding function
          
          call metal_deriv(katm1,fmet,buffer)
          
        endif
        
      enddo

      if(idnode.eq.0)close (ntable)
      
      if(idnode.eq.0)write(nrite,'(/,/,1x,a)')
     x  'potential tables read from TABEAM file'
      
      return
      end subroutine mettab

      subroutine metal_deriv(ityp,vvv,buffer)

c**********************************************************************
c
c     calculate numerical derivatives of tabulated EAM metal potentials
c     
c     copyright - daresbury laboratory
c     author    - w.smith march 2006
c
c**********************************************************************
      
      implicit none
      
      integer ityp,i,npt
      real(8) vvv,buffer,delmet,aa0,aa1,aa2,aa3,aa4,d1y,d2y,d3y,d4y
      real(8) f0,f1,f2,f3,f4

      dimension vvv(mxgrid,mxmet,2),buffer(mxbuff)
      
c     interpolation parameters
      
      vvv(1,ityp,2)=buffer(1)
      vvv(2,ityp,2)=buffer(2)
      vvv(3,ityp,2)=buffer(3)
      vvv(4,ityp,2)=buffer(4)
      
c     construct interpolation table

      delmet=buffer(4)
      npt=nint(buffer(1))-2
      do i=7,npt

        aa0=buffer(i)
        f0=buffer(i-2)/aa0
        f1=buffer(i-1)/aa0
        f2=1.d0
        f3=buffer(i+1)/aa0
        f4=buffer(i+2)/aa0
        
c     calculate numerical differences for 5-point interpolation
        
        d1y=(f1-f0)
        d2y=(f2-f1)-(f1-f0)
        d3y=(f3-f0)+3.d0*(f1-f2)
        d4y=(f4-f3)+3.d0*(f2-f3)+3.d0*(f2-f1)+(f0-f1)
        
c     calculate polynomial coefficients
        
        aa0=aa0/delmet
        aa4=d4y/24.d0
        aa3=(d3y+12.d0*aa4)/6.d0
        aa2=(d2y+6.d0*aa3-14.d0*aa4)/2.d0
        aa1=d1y+3.d0*aa2-7.d0*aa3+15.d0*aa4
        
c     calculate derivatives
        
        vvv(i,ityp,2)=aa1*aa0
        
c     derivatives at extremes of range
        
        if(i.eq.7)then
          
          vvv(5,ityp,2)=(aa1-4.d0*aa2+12.d0*aa3-32.d0*aa4)*aa0
          vvv(6,ityp,2)=(aa1-2.d0*aa2+3.d0*aa3-4.d0*aa4)*aa0
          
        else if(i.eq.npt)then
          
          vvv(npt+1,ityp,2)=(aa1+2.d0*aa2+3.d0*aa3+4.d0*aa4)*aa0
          vvv(npt+2,ityp,2)=(aa1+4.d0*aa2+12.d0*aa3+32.d0*aa4)*aa0
          
        endif
          
      enddo
      
c     set derivatives to zero beyond end point of function
      
      do i=npt+3,mxgrid
        vvv(i,ityp,2)=0.d0
      enddo
      
      return
      end subroutine metal_deriv

      subroutine abort_eamtable_read(idnode,ntable)

c***********************************************************************
c     
c     dl_poly error exit subroutine for reading TABEAM file
c     
c     copyright - daresbury laboratory 
c     author    - w. smith   mar 2006
c     
c***********************************************************************

      implicit none
      integer idnode,ntable

      if(idnode.eq.0)close (ntable)
      
      call error(idnode,29)
      
      end subroutine abort_eamtable_read
      
      end module metal_module
