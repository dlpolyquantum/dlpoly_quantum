      module vdw_module

c***********************************************************************
c     
c     dl_poly module for defining van der waals potential arrays
c     copyright - daresbury laboratory
c     author    - w. smith    sep 2003
c     adapted   - p.-a. cazade oct 2007: solvation, free energy etc
c     adapted   - d. quigley : metadynamics
c     
c***********************************************************************

      use config_module
      use error_module
      use metafreeze_module
      use pair_module
      use parse_module
      use setup_module
      use site_module
      use solvation_module
      use utility_module

      implicit none

      integer, allocatable :: ltpvdw(:),lstvdw(:)
      real(8), allocatable :: vvv(:,:),ggg(:,:),prmvdw(:,:)

      save ltpvdw,lstvdw,prmvdw,vvv,ggg

      contains
      
      subroutine alloc_vdw_arrays(idnode,mxnode)
      
      implicit none
      
      integer, parameter :: nnn=5
      
      logical safe
      integer i,fail,idnode,mxnode
      dimension fail(nnn)
      
      safe=.true.
      
c     allocate arrays
      
      fail(:)=0
      
      allocate (ltpvdw(mxvdw),stat=fail(1))
      allocate (lstvdw(mxvdw),stat=fail(2))
      allocate (prmvdw(mxvdw,mxpvdw),stat=fail(3))
      allocate (vvv(mxgrid,mxvdw),stat=fail(4))
      allocate (ggg(mxgrid,mxvdw),stat=fail(5))
      
      if(any(fail.gt.0))safe=.false.      
      if(mxnode.gt.1)call gstate(safe)    
      if(.not.safe)call error(idnode,1014)
      
      end subroutine alloc_vdw_arrays

      subroutine define_van_der_waals
     x  (safe,ltable,lunits,lmols,idnode,ntpvdw,
     x  ntpatm,keyfce,dlrpot,rvdw,engunit)

c***********************************************************************
c     
c     dl_poly subroutine for defining van der Waals potentials
c     
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2003
c     
c***********************************************************************

      implicit none

      logical safe,ltable,lunits,lmols
      character*1 message(80)
      character*8 atom1,atom2,keyword

      integer ntpvdw,ntpatm,keyfce,fail,idum,ivdw
      integer itpvdw,keypot,numpar,katom1,katom2,jtpatm,keyvdw,i
      integer ntab,idnode,j
      real(8) dlrpot,rvdw,engunit
      real(8), allocatable :: parpot(:)

      allocate (parpot(mxpvdw),stat=fail)

      ntpvdw=intstr(record,lenrec,idum)

      ltable=findstring('table',record,idum)
      
      if(idnode.eq.0) then
        
        write(nrite,"(/,/,1x,'number of specified pair ',
     x    'potentials',i10)") ntpvdw
        write(nrite,"(/,/,16x,'atom 1  ','atom 2  ',3x,
     x    ' key',30x,'parameters'/,/)")
        
      endif      

      if(ntpvdw.gt.mxvdw) call error(idnode,80)
      if(.not.lunits) call error(idnode,6)
      if(.not.lmols) call error(idnode,13)
      
      do ivdw=1,mxvdw
        
        lstvdw(ivdw)=0
        ltpvdw(ivdw)=-1
        
      enddo
      
      do itpvdw=1,ntpvdw
        
        do i=1,mxpvdw
          parpot(i)=0.d0
        enddo
        
        call getrec(safe,idnode,nfield)
        if(.not.safe)return

        call copystring(record,message,80)
        call getword(atom1,record,8,lenrec)
        call getword(atom2,record,8,lenrec)
        call lowcase(record,lenrec-16)
        call getword(keyword,record,4,lenrec)

        if(keyword(1:4).eq.'12-6') then
          keypot=1
          numpar=2
        elseif(keyword(1:4).eq.'lj  ') then
          keypot=2
          numpar=2
        elseif(keyword(1:4).eq.'nm  ') then
          keypot=3
          numpar=4
        elseif(keyword(1:4).eq.'buck') then
          keypot=4
          numpar=3
        elseif(keyword(1:4).eq.'bhm ') then
          keypot=5
          numpar=5
        elseif(keyword(1:4).eq.'hbnd') then
          keypot=6
          numpar=2
        elseif(keyword(1:4).eq.'snm ') then
          keypot=7
          numpar=5
        elseif(keyword(1:4).eq.'mors') then
          keypot=8
          numpar=3
        elseif(keyword(1:4).eq.'wca ') then
          keypot=9
          numpar=3
        elseif(keyword(1:4).eq.'gaus') then
          keypot=10
          numpar=6
        elseif(keyword(1:4).eq.'tab ') then
          keypot=0
          numpar=0
        else
          if(idnode.eq.0) write(nrite,*) message
          call error(idnode,452)
        endif

        parpot(1)=dblstr(record,lenrec,idum)
        parpot(2)=dblstr(record,lenrec,idum)
        parpot(3)=dblstr(record,lenrec,idum)
        parpot(4)=dblstr(record,lenrec,idum)
        parpot(5)=dblstr(record,lenrec,idum)
        parpot(6)=dblstr(record,lenrec,idum)
        
        if(idnode.eq.0) 
     x    write(nrite,"(16x,2a8,2x,a4,3x,1p,9e13.5)") 
     x    atom1,atom2,keyword(1:4),(parpot(j),j=1,numpar)
        
        katom1=0
        katom2=0
        
        do jtpatm=1,ntpatm

          if(atom1.eq.unqatm(jtpatm))katom1=jtpatm
          if(atom2.eq.unqatm(jtpatm))katom2=jtpatm
          
        enddo
        
        if(katom1.eq.0.or.katom2.eq.0) then
          call  error(idnode,81)
        endif
        
        keyvdw=loc2(katom1,katom2)

c     convert energies to internal unit

        if(keyvdw.gt.mxvdw) call error(idnode,82)
        
        parpot(1)=parpot(1)*engunit
        
        if(keypot.eq.1) then
          
          parpot(2)=parpot(2)*engunit
          
        else if(keypot.eq.4) then
          
          parpot(3)=parpot(3)*engunit
          
        else if(keypot.eq.5) then
          
          parpot(4)=parpot(4)*engunit
          parpot(5)=parpot(5)*engunit
          
        else if(keypot.eq.6) then
          
          parpot(2)=parpot(2)*engunit
          
        else if(keypot.eq.10) then
          
          parpot(3)=parpot(3)*engunit
          parpot(5)=parpot(5)*engunit
          
        endif

        ltable=(ltable.or.(keypot.eq.0))

        if(lstvdw(keyvdw).ne.0) call error(idnode,15)
        lstvdw(keyvdw)=itpvdw
        ltpvdw(itpvdw)=keypot
        
        do i=1,mxpvdw
          
          prmvdw(itpvdw,i)=parpot(i)
          
        enddo
        
      enddo

c     generate nonbonded force arrays

      if((ntpvdw.gt.0.and.mod(keyfce,2).eq.1).or.(keyfce.eq.2))
     x  then
        
        call forgen(ltable,idnode,ntpvdw,dlrpot,rvdw)
        
        if(ltable)then
          
          call fortab
     x      (idnode,ntpvdw,ntpatm,dlrpot,rvdw,engunit)
          
        endif
        
      endif

c     check for unspecified atom-atom potentials
      
      ntab=(ntpatm*(ntpatm+1))/2
      
      if(ntpvdw.lt.ntab) then
        
        call warning(idnode,110,0.d0,0.d0,0.d0)

        if(mxvdw.le.ntpvdw) call error(idnode,82)

        do i=1,ntab
          
          if(lstvdw(i).eq.0)then
            
            lstvdw(i)=ntpvdw+1
            
          endif
          
        enddo

c     define zero potential for undefined interactions
        
        do i=1,mxgrid
          
          ggg(i,ntpvdw+1)=0.d0
          vvv(i,ntpvdw+1)=0.d0
          
        enddo
        
      endif

      deallocate (parpot,stat=fail)

      return
      end subroutine define_van_der_waals

      subroutine forgen(ltable,idnode,ntpvdw,dlrpot,rcut)

c***********************************************************************
c     
c     dl_poly subroutine for generating potential energy and 
c     force arrays for van der waals forces only
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith may 1992.
c     
c***********************************************************************
      
      implicit none

      logical ltable
      integer i,ivdw,ntpvdw,idnode
      real(8) dlrpot,rcut,rrr,ann,amm,gam,bet,eps,rr0,aaa,bbb
      real(8) ccc,ddd,eee,sig,rho,rrc,aa1,aa2,aa3,ee1,ee2,ee3
      real(8) rsq,ex1,ex2,ex3

c     define grid resolution for potential arrays
      
      dlrpot=rcut/dble(mxgrid-4)

c     construct arrays for all types of short ranged  potential
      
      do ivdw=1,ntpvdw
        
        if(ltpvdw(ivdw).eq.1)then
          
c       12 - 6 potential
      
          aaa=prmvdw(ivdw,1)
          bbb=prmvdw(ivdw,2)
          
          do i=1,mxgrid
            
            rrr=dble(i)*dlrpot
            vvv(i,ivdw)=(aaa/rrr**6-bbb)/rrr**6
            ggg(i,ivdw)=6.d0*(2.d0*aaa/rrr**6-bbb)/rrr**6
            
          enddo
          
        else if(ltpvdw(ivdw).eq.2)then
          
c       lennard-jones potential
      
          eps=prmvdw(ivdw,1)
          sig=prmvdw(ivdw,2)

          do i=1,mxgrid
            
            rrr=dble(i)*dlrpot
            vvv(i,ivdw)=4.d0*eps*(sig/rrr)**6*((sig/rrr)**6-1.d0)
            ggg(i,ivdw)=24.d0*eps*(sig/rrr)**6*(2.d0*(sig/rrr)**6-1.d0)
            
          enddo
          
        else if(ltpvdw(ivdw).eq.3)then

c       n - m potential
      
          eps=prmvdw(ivdw,1)
          ann=max(prmvdw(ivdw,2),prmvdw(ivdw,3))
          amm=min(prmvdw(ivdw,2),prmvdw(ivdw,3))
          rr0=prmvdw(ivdw,4)

          do i=1,mxgrid
            
            rrr=dble(i)*dlrpot
            vvv(i,ivdw)=eps/(ann-amm)*(amm*(rr0/rrr)**ann-
     x        ann*(rr0/rrr)**amm)
            ggg(i,ivdw)=eps*amm*ann/(ann-amm)*((rr0/rrr)**ann-
     x        (rr0/rrr)**amm)
            
          enddo
          
        else if(ltpvdw(ivdw).eq.4)then
          
c       buckingham exp - 6 potential
      
          aaa=prmvdw(ivdw,1)
          rho=prmvdw(ivdw,2)
          ccc=prmvdw(ivdw,3)

          do i=1,mxgrid
            
            rrr=dble(i)*dlrpot
            vvv(i,ivdw)=aaa*exp(-rrr/rho)-ccc/rrr**6
            ggg(i,ivdw)=rrr*aaa*exp(-rrr/rho)/rho-6.d0*ccc/rrr**6
            
          enddo
          
        else if(ltpvdw(ivdw).eq.5)then
          
c       born-huggins-meyer exp - 6 - 8 potential
      
          aaa=prmvdw(ivdw,1)
          bbb=prmvdw(ivdw,2)
          ccc=prmvdw(ivdw,3)
          ddd=prmvdw(ivdw,4)
          eee=prmvdw(ivdw,5)

          do i=1,mxgrid
            
            rrr=dble(i)*dlrpot
            vvv(i,ivdw)=aaa*exp(bbb*(ccc-rrr))-ddd/rrr**6-eee/rrr**8
            ggg(i,ivdw)=rrr*aaa*bbb*exp(bbb*(ccc-rrr))-6.d0*ddd/rrr**6
     x        -8.d0*eee/rrr**8
            
          enddo
          
        else if(ltpvdw(ivdw).eq.6) then
          
c       Hydrogen-bond 12 - 10 potential
      
          aaa=prmvdw(ivdw,1)
          bbb=prmvdw(ivdw,2)

          do i=1,mxgrid
            
            rrr=dble(i)*dlrpot
            vvv(i,ivdw)=aaa/rrr**12-bbb/rrr**10
            ggg(i,ivdw)=12.0d0*aaa/rrr**12-10.d0*bbb/rrr**10
            
          enddo
          
        else if(ltpvdw(ivdw).eq.7) then
          
c       shifted and force corrected n - m potential (w. smith)
      
          eps=prmvdw(ivdw,1)
          ann=prmvdw(ivdw,2)
          amm=prmvdw(ivdw,3)
          rr0=prmvdw(ivdw,4)
          rrc=prmvdw(ivdw,5)
          if(rrc.lt.1.d-6)rrc=rcut
 
          if(ann.le.amm) call error(idnode,470)

          gam=rrc/rr0
          if(gam.lt.1.d0) call error(idnode,468)
          bet=gam*((gam**(amm+1.d0)-1.d0)/(gam**(ann+1.d0)-1.d0))
     x      **(1.d0/(ann-amm))
          eps=-eps*(ann-amm)/(amm*(bet**ann)*(1.d0+(ann/gam-ann-1.d0)
     x      /gam**ann)-ann*(bet**amm)*(1.d0+(amm/gam-amm-1.d0)
     x      /gam**amm))

          do i=1,mxgrid
            
            rrr=dble(i)*dlrpot
            if(rrr.gt.rrc)then

              vvv(i,ivdw)=0.d0
              ggg(i,ivdw)=0.d0

            else

              vvv(i,ivdw)=eps/(ann-amm)*(amm*(bet**ann)*((rr0/rrr)**ann-
     x          (1.d0/gam)**ann)-ann*(bet**amm)*((rr0/rrr)**amm-
     x          (1.d0/gam)**amm)+ann*amm*((rrr/(gam*rr0)-1.d0)*
     x          ((bet/gam)**ann-(bet/gam)**amm)))
              ggg(i,ivdw)=eps*amm*ann/(ann-amm)*((bet**ann)*
     x          (rr0/rrr)**ann-(bet**amm)*(rr0/rrr)**amm-rrr/
     x          (gam*rr0)*((bet/gam)**ann-(bet/gam)**amm))

            endif

          enddo
          
        else if(ltpvdw(ivdw).eq.8) then
          
c       morse potential
          
          eps=prmvdw(ivdw,1)
          rr0=prmvdw(ivdw,2)
          sig=prmvdw(ivdw,3)
          
          do i=1,mxgrid
            
            rrr=dble(i)*dlrpot
            vvv(i,ivdw)=eps*((1.d0-exp(-sig*(rrr-rr0)))**2-1.d0)
            ggg(i,ivdw)=-2.d0*rrr*eps*sig*(1.d0-exp(-sig*(rrr-rr0)))*
     x        exp(-sig*(rrr-rr0))
            
          enddo
          
        else if(ltpvdw(ivdw).eq.9) then
          
c       weeks-chandler-anderson potential
          
          eps=prmvdw(ivdw,1)
          sig=prmvdw(ivdw,2)
          rr0=prmvdw(ivdw,3)
          ddd=sig*2.d0**(1.d0/6.d0)
          
          do i=1,mxgrid
            
            rrr=dble(i)*dlrpot-rr0
            if(rrr.gt.ddd)then
              
              vvv(i,ivdw)=0.d0
              ggg(i,ivdw)=0.d0

            else if(rrr.gt.0.d0)then
              
              vvv(i,ivdw)=4.d0*eps*(sig/rrr)**6*
     x          ((sig/rrr)**6-1.d0)+eps
              ggg(i,ivdw)=24.d0*eps*(1.d0+rr0/rrr)*(sig/rrr)**6*
     x          (2.d0*(sig/rrr)**6-1.d0)
            
            endif
              
          enddo
          
        else if(ltpvdw(ivdw).eq.10)then
          
c       gaussian potential
      
          aa1=prmvdw(ivdw,1)
          ee1=prmvdw(ivdw,2)
          aa2=prmvdw(ivdw,3)
          ee2=prmvdw(ivdw,4)
          aa3=prmvdw(ivdw,5)
          ee3=prmvdw(ivdw,6)

          do i=1,mxgrid
            
            rsq=(dble(i)*dlrpot)**2
            ex1=aa1*exp(-rsq*ee1)
            ex2=aa2*exp(-rsq*ee2)
            ex3=aa3*exp(-rsq*ee3)
            vvv(i,ivdw)=ex1+ex2+ex3
            ggg(i,ivdw)=2.d0*rsq*(ee1*ex1+ee2*ex2+ee3*ex3)
            
          enddo
          
        else if(ltpvdw(ivdw).lt.100) then
          
          if(.not.ltable)call error(idnode,150)
          
        endif
        
      enddo
      
      return
      end subroutine forgen

      subroutine fortab
     x  (idnode,ntpvdw,ntpatm,dlrpot,rcut,engunit)

c***********************************************************************
c     
c     dl_poly subroutine for reading potential energy and 
c     force arrays for van der waals forces only
c     
c     copyright - daresbury laboratory 1994
c     author    - w. smith march 1994
c     
c***********************************************************************
      
      implicit none

      logical safe
      character*8 atom1,atom2
      integer idnode,ntpvdw,ntpatm,idum,ngrid
      integer ivdw,katom1,katom2,jtpatm,l,i,j,k,keyvdw
      real(8) dlrpot,rcut,engunit,delpot,cutpot,rdr,rrr,ppp
      real(8) vk0,vk1,vk2,t1,t2

      if(idnode.eq.0)open (ntable,file='TABLE')

c     skip header record
      
      call getrec(safe,idnode,ntable)
      if(.not.safe)call abort_table_read(idnode,ntable)

c     read mesh resolution
      
      call getrec(safe,idnode,ntable)
      if(.not.safe)call abort_table_read(idnode,ntable)
      delpot=dblstr(record,lenrec,idum)
      cutpot=dblstr(record,lenrec,idum)
      ngrid=intstr(record,lenrec,idum)

      dlrpot=rcut/dble(mxgrid-4)

      if (abs(delpot-dlrpot) <= 1.0d-8) delpot=dlrpot
      if ((delpot>dlrpot) .or. (ngrid-4 /= nint(cutpot/delpot))) then
        
         if (idnode == 0) write(nrite,"(                 
     x    'expected radial increment : ',1p,e15.7,/,     
     x    'TABLE    radial increment : ',1p,e15.7,/,/,   
     x    'expected number of grid points : ',0p,i10,/,  
     x    'grid points in TABLE           : ',i10)")     
     x    dlrpot, delpot, mxgrid, ngrid
         
         call error(idnode,22)
         
      endif

      if(cutpot.lt.rcut) call error(idnode,504)
      if(abs(1.d0-(delpot/dlrpot)).gt.1.0d-8) then
        if(idnode.eq.0) write(nrite,
     x    "(/,' TABLE arrays resized for mxgrid=',i10)") mxgrid
      endif

c     read potential arrays for all pairs
      
      do ivdw=1,ntpvdw

c     read potential arrays if potential not already defined
        
        if(ltpvdw(ivdw).eq.0)then
          
c     read pair potential labels and long range corrections
          
          call getrec(safe,idnode,ntable)
          if(.not.safe)call abort_table_read(idnode,ntable)

          call getword(atom1,record,8,lenrec)
          call getword(atom2,record,8,lenrec)
          prmvdw(ivdw,1)=dblstr(record,lenrec,idum)
          prmvdw(ivdw,2)=dblstr(record,lenrec,idum)
          
          katom1=0
          katom2=0
          
          do jtpatm=1,ntpatm
            
            if(atom1.eq.unqatm(jtpatm))katom1=jtpatm
            if(atom2.eq.unqatm(jtpatm))katom2=jtpatm
            
          enddo
          
          if(katom1.eq.0.or.katom2.eq.0)then
            if(idnode.eq.0) 
     x        write(nrite,'(a)') '****',atom1,'***',atom2,'****'
            call  error(idnode,81)
          endif
          
          keyvdw=loc2(katom1,katom2)
          
          if(lstvdw(keyvdw).ne.ivdw) call error(idnode,23)
          
c     read potential arrays
          
          if(mxbuff.lt.ngrid)  then
              
            if(idnode.eq.0)
     x         write(nrite,*) 'mxbuff must be >=',ngrid,' in fortab'
            call error(idnode,48)
              
          endif

c     read in potential arrays

          do i=1,(ngrid+3)/4
            
             l=min(4,ngrid-(i-1)*4)
             if (idnode == 0) then
                read(unit=ntable, fmt=*, end=100)
     x              (buffer((i-1)*4+j),j=1,l)
             else
                buffer((i-1)*4+1:(i-1)*4+l)=0.0d0
             endif
             
          enddo
          call gdsum(buffer(1:ngrid),ngrid,buffer(ngrid+1:2*ngrid))

c     reconstruct arrays using 3pt interpolation

          rdr=1.d0/delpot
          vvv(1,ivdw)=1.d0
          ggg(1,ivdw)=0.d0
          do i=2,mxgrid
            
            rrr=dble(i)*dlrpot
            l=int(rrr*rdr)
            ppp=rrr*rdr-dble(l)
            vk0=buffer(l)
            vk1=buffer(l+1)
            vk2=buffer(l+2)
            
            t1=vk0+(vk1-vk0)*ppp
            t2=vk1+(vk2-vk1)*(ppp-1.0d0)
            vvv(i,ivdw)=t1+(t2-t1)*ppp*0.5d0

          enddo

c     read in force arrays

          do i=1,(ngrid+3)/4
            
             l=min(4,ngrid-(i-1)*4)
             if (idnode == 0) then
                read(unit=ntable, fmt=*, end=100)
     x              (buffer((i-1)*4+j),j=1,l)
             else
                buffer((i-1)*4+1:(i-1)*4+l)=0.0d0
             endif
             
          enddo
          call gdsum(buffer(1:ngrid),ngrid,buffer(ngrid+1:2*ngrid))

c     reconstruct ggg arrays using 3pt interpolation

          do i=2,mxgrid

            rrr=dble(i)*dlrpot
            l=int(rrr*rdr)
            ppp=rrr*rdr-dble(l)
            vk0=buffer(l)
            vk1=buffer(l+1)
            vk2=buffer(l+2)
            
            t1=vk0+(vk1-vk0)*ppp
            t2=vk1+(vk2-vk1)*(ppp-1.0d0)
            
            ggg(i,ivdw)=t1+(t2-t1)*ppp*0.5d0

          enddo

        endif
        
      enddo

c     convert to internal units
      
      do k=1,ntpvdw
        
        if(ltpvdw(k).eq.0)then

          do i=1,mxgrid
            
            vvv(i,k)=vvv(i,k)*engunit
            ggg(i,k)=ggg(i,k)*engunit
            
          enddo
          
        endif
        
      enddo
      
      if(idnode.eq.0)close (ntable)
      
      if(idnode.eq.0)write(nrite,'(/,/,1x,a)')
     x  'potential tables read from TABLE file'
      
      return
      
c     end of file error exit
      
  100 call abort_table_read(idnode,ntable)

      end subroutine fortab

      subroutine abort_table_read(idnode,ntable)

c***********************************************************************
c     
c     dl_poly error exit subroutine for reading TABLE file
c     
c     copyright - daresbury laboratory 
c     author    - w. smith   sept 2005
c     
c***********************************************************************

      implicit none
      integer idnode,ntable

      if(idnode.eq.0)close (ntable)
      
      call error(idnode,24)
      
      end subroutine abort_table_read

      subroutine srfrce
     x  (lsolva,lfree,lghost,iatm,ik,engsrp,virsrp,rcut,dlrpot)

c***********************************************************************
c     
c     dl_poly subroutine for calculating short range force and
c     potential energy terms using verlet neighbour list
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith       march 1992
c     
c     version 3
c     author    - t. forester    june  1993
c     stress tensor added t.forester may 1994
c     adapted   - p.-a. cazade oct 2007: solvation, free energy etc
c     adapted   - d. quigley - metadynamics
c     
c***********************************************************************
      
      implicit none

      logical lsolva,lfree,lghost,lselect,lskip,idrive,jdrive
      integer iatm,ik,m,jatm,k,l,kkk
      real(8) engsrp,virsrp,rcut,dlrpot
      real(8) ab,rrr,rsq,ppp,t1,t2,vk0,vk1,vk2,gk0,gk1,gk2,gamma
      real(8) fi,rcsq,rdr,ai,aj,fx,fy,fz,omega
      real(8) strs(6),strs_loc(6)

      dimension fi(3)

CDIR$ CACHE_ALIGN fi
      
      lskip=(lfree.or.lghost)
      if(lmetadyn)idrive=driven(ltype(iatm))
      
c     set cutoff condition for pair forces

      rcsq=rcut**2

c     interpolation spacing
      
      rdr=1.d0/dlrpot

c     initialise stress tensor accumulators

      strs(:)=0.d0
      strs_loc(:)=0.d0

c     initialise potential energy and virial
      
      engsrp=0.d0
      virsrp=0.d0

c     store forces for iatm 
      
      ai=dble(ltype(iatm))
      fi(1)=fxx(iatm)
      fi(2)=fyy(iatm)
      fi(3)=fzz(iatm)

c     start of primary loop for forces evaluation
      
      do m=1,ik
        
c     atomic and potential function indices
        
        jatm=ilist(m)
        if(lmetadyn)jdrive=driven(ltype(jatm))
        
        if(lskip)then
          if(atm_fre(iatm)*atm_fre(jatm).eq.2)cycle
        endif
        
        aj=dble(ltype(jatm))
        
        if(ai.gt.aj) then
          ab=ai*(ai-1.d0)*0.5d0+aj+0.5d0
        else
          ab=aj*(aj-1.d0)*0.5d0+ai+0.5d0
        endif
        
        k=lstvdw(int(ab))
        
        if((ltpvdw(k).lt.100).and.(abs(vvv(1,k)).gt.1.d-10))then
          
c     apply truncation of potential
          
          rsq=rsqdf(m)
          
          if(rcsq.gt.rsq)then
            
            rrr=sqrt(rsq)               
            l=int(rrr*rdr)
            ppp=rrr*rdr-dble(l)
            
            if(l.eq.0)then
              
              omega=vvv(1,k)
              gamma=ggg(1,k)
              
            else
              
c     calculate interaction energy using 3-point interpolation
              
              vk0=vvv(l,k)
              vk1=vvv(l+1,k)
              vk2=vvv(l+2,k)
              t1=vk0+(vk1-vk0)*ppp
              t2=vk1+(vk2-vk1)*(ppp-1.0d0)
              omega=t1+(t2-t1)*ppp*0.5d0
              
c     calculate forces using 3-point interpolation
              
              gk0=ggg(l,k)
              gk1=ggg(l+1,k)
              gk2=ggg(l+2,k)
              t1=gk0+(gk1-gk0)*ppp
              t2=gk1+(gk2-gk1)*(ppp-1.0d0)
              gamma=(t1+(t2-t1)*ppp*0.5d0)/rsq
              
            endif
            
c     set selection control
            
            lselect=.true.
            
c     set double index
            
            if(lsolva)kkk=loc2(atmolt(iatm),atmolt(jatm))
            
            if(lghost)then
              
c     selected excitation option
              
              if((atm_fre(iatm).ne.1).and.(atm_fre(jatm).ne.1))then
                
c     reset selection control
                
                lselect=(atm_fre(iatm)+atm_fre(jatm).eq.0)
                
                if(lsolva)vdw_exc(kkk)=vdw_exc(kkk)+omega
                
              endif
              
            elseif(lfree)then
              
c     selected free energy option
              
              if((atm_fre(iatm).eq.1).or.(atm_fre(jatm).eq.1))then
                
c     set hamiltonian mixing parameter
                
                vdw_fre=vdw_fre-omega
                vdw_vir=vdw_vir+gamma*rsq
                omega=lambda1*omega
                gamma=lambda1*gamma
                
              elseif((atm_fre(iatm).eq.2).or.(atm_fre(jatm).eq.2))then
                
c     set hamiltonian mixing parameter
                
                vdw_fre=vdw_fre+omega
                vdw_vir=vdw_vir-gamma*rsq
                omega=lambda2*omega
                gamma=lambda2*gamma
                
              endif
              
            endif
            
            if(lselect)then
              
c     calculate potential and virial
              
              engsrp=engsrp+omega
              virsrp=virsrp-gamma*rsq
              
              if(lsolva)vdw_sol(kkk)=vdw_sol(kkk)+omega
              
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
              
            endif
            
            if(lmetadyn.and.(idrive.or.jdrive))then
              
              eng_loc=eng_loc+omega
              vir_loc=vir_loc-gamma*rsq
              
              fxx_loc(iatm)=fxx_loc(iatm)+fx
              fyy_loc(iatm)=fyy_loc(iatm)+fy
              fzz_loc(iatm)=fzz_loc(iatm)+fz
              
              fxx_loc(jatm)=fxx_loc(jatm)-fx
              fyy_loc(jatm)=fyy_loc(jatm)-fy
              fzz_loc(jatm)=fzz_loc(jatm)-fz
              
              strs_loc(1)=strs_loc(1)+xdf(m)*fx
              strs_loc(2)=strs_loc(2)+xdf(m)*fy
              strs_loc(3)=strs_loc(3)+xdf(m)*fz
              strs_loc(4)=strs_loc(4)+ydf(m)*fy
              strs_loc(5)=strs_loc(5)+ydf(m)*fz
              strs_loc(6)=strs_loc(6)+zdf(m)*fz
              
            endif
            
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

      if(lmetadyn)then
        
        stress_loc(1)=stress_loc(1)+strs_loc(1)
        stress_loc(2)=stress_loc(2)+strs_loc(2)
        stress_loc(3)=stress_loc(3)+strs_loc(3)
        stress_loc(4)=stress_loc(4)+strs_loc(2)
        stress_loc(5)=stress_loc(5)+strs_loc(4)
        stress_loc(6)=stress_loc(6)+strs_loc(5)
        stress_loc(7)=stress_loc(7)+strs_loc(3)
        stress_loc(8)=stress_loc(8)+strs_loc(5)
        stress_loc(9)=stress_loc(9)+strs_loc(6)
        
      endif

      return
      end subroutine srfrce
      
      subroutine lrcorrect
     x  (lsolva,lfree,lghost,idnode,imcon,keyfce,natms,
     x  ntpatm,ntpvdw,elrc,engunit,virlrc,rcut,volm)
      
c*************************************************************************
c     
c     DL_POLY subroutine to evaluate long-range corrections to
c     pressure and energy in a periodic system.
c     
c     copyright daresbury laboratory 1993
c     author    - t. forester may 1993
c     adapted   - p.-a. cazade oct 2007: solvation, free energy etc
c     
c***************************************************************************
      
      implicit none

      integer, parameter :: nnn=10
      logical lsolva,lfree,lghost
      integer idnode,imcon,keyfce,natms,ntpatm,i,ka,ntpvdw
      integer ivdw,j,k,it,jt,kt,fail
      real(8) natyp,nbtyp,nctyp,ndtyp,nafrz,nbfrz,ncfrz,ndfrz
      real(8) elrc,engunit,virlrc,rcut,volm,twopi,eadd,padd
      real(8) denprd,aaa,bbb,ccc,ddd,eee,eps,sig,rr0,ann,amm
      real(8) denprd1,denprd2,denprd3,denprdf
      integer, allocatable :: numtyp_sol0(:,:),numfrz_sol0(:,:)
      integer, allocatable :: numtyp_sol1(:,:),numfrz_sol1(:,:)
      integer, allocatable :: numtyp_sol2(:,:),numfrz_sol2(:,:)
      integer, allocatable :: numtyp_fre(:,:),numfrz_fre(:,:)
      real(8), allocatable :: elrc_sol0(:),elrc_exc0(:)
      
      dimension fail(nnn)
      
      twopi=2.0d0*pi
      
c     allocate working arrays
      
      do i=1,nnn
        fail(i)=0
      enddo
      
      if(lfree.or.lghost)then
        
        allocate (numtyp_fre(mxatyp,0:2),stat=fail(1))
        allocate (numfrz_fre(mxatyp,0:2),stat=fail(2))
        allocate (elrc_exc0(mxtmls_exc2),stat=fail(3))
        
      endif
      
      if(lsolva)then
        
        allocate (elrc_sol0(mxtmls_sol2),stat=fail(4))
        allocate (numtyp_sol0(mxatyp,mxtmls),stat=fail(5))
        allocate (numfrz_sol0(mxatyp,mxtmls),stat=fail(6))
        
        if(lghost)then
          
          allocate (numtyp_sol1(mxatyp,mxtmls),stat=fail(7))
          allocate (numfrz_sol1(mxatyp,mxtmls),stat=fail(8))
          allocate (numtyp_sol2(mxatyp,mxtmls),stat=fail(9))
          allocate (numfrz_sol2(mxatyp,mxtmls),stat=fail(10))
          
        endif
        
      endif
      
      do i=1,nnn
        if(fail(i).ne.0)call error(idnode,1015)
      enddo
      
c     initalise counter arrays
      
      do i=1,ntpatm
        
        numtyp(i)=0
        numfrz(i)=0
        
      enddo
      
      if(lfree.or.lghost)then
        
        numtyp_fre(:,:)=0
        numfrz_fre(:,:)=0
        
      endif
      
      if(lsolva)then
        
        numtyp_sol0(:,:)=0
        numfrz_sol0(:,:)=0
        
        if(lghost)then
          
          numtyp_sol1(:,:)=0
          numfrz_sol1(:,:)=0
          numtyp_sol2(:,:)=0
          numfrz_sol2(:,:)=0
          
        endif
        
      endif
      
c     evaluate number density in system
      
      do i=1,natms
        
        ka=ltype(i)
        numtyp(ka)=numtyp(ka)+1
        if(lstfrz(i).ne.0)numfrz(ka)=numfrz(ka)+1
        
      enddo
      
      if(lfree.or.lghost)then
         
        do i=1,natms
          
          ka=ltype(i)
          numtyp_fre(ka,atm_fre(i))=numtyp_fre(ka,atm_fre(i))+1
          if(lstfrz(i).ne.0)
     x      numfrz_fre(ka,atm_fre(i))=numfrz_fre(ka,atm_fre(i))+1          
          
        enddo
        
      endif
      
      if(lsolva)then
        
        if(lghost)then
          
          do i=1,natms
            
            ka=ltype(i)
            
            if(atm_fre(i).eq.0)then
              
              numtyp_sol0(ka,atmolt(i))=numtyp_sol0(ka,atmolt(i))+1
              if(lstfrz(i).ne.0)
     x          numfrz_sol0(ka,atmolt(i))=numfrz_sol0(ka,atmolt(i))+1
              
            elseif(atm_fre(i).eq.1)then
              
              numtyp_sol1(ka,atmolt(i))=numtyp_sol1(ka,atmolt(i))+1
              if(lstfrz(i).ne.0)
     x          numfrz_sol1(ka,atmolt(i))=numfrz_sol1(ka,atmolt(i))+1
              
            elseif(atm_fre(i).eq.2)then
              
              numtyp_sol2(ka,atmolt(i))=numtyp_sol2(ka,atmolt(i))+1
              if(lstfrz(i).ne.0)
     x          numfrz_sol2(ka,atmolt(i))=numfrz_sol2(ka,atmolt(i))+1
              
            endif
            
          enddo
          
        else
          
          do i=1,natms
            
            ka=ltype(i)
            numtyp_sol0(ka,atmolt(i))=numtyp_sol0(ka,atmolt(i))+1
            if(lstfrz(i).ne.0)
     x        numfrz_sol0(ka,atmolt(i))=numfrz_sol0(ka,atmolt(i))+1
            
          enddo
          
        endif
        
      endif
      
c     number densities
      
      do i=1,ntpatm
        dens(i)=dble(numtyp(i))/volm
      enddo
      
c     long range corrections to energy and pressure
      
      elrc=0.d0
      elrc2=0.d0
      virlrc=0.d0
      virlrc2=0.d0
      denprdf=0.d0
      volm_sav=0.d0
      elrc_fre=0.d0
      vlrc_fre=0.d0
      
      if(imcon.ne.0.and.imcon.ne.6.and.ntpvdw.gt.0) then 
         
        if(mod(keyfce,2).eq.1) then
          
          ivdw=0
          
          do i=1,ntpatm
            
            do j=1,i
               
              eadd=0.d0
              padd=0.d0
              
              ivdw=ivdw+1
              k=lstvdw(ivdw)
              
              if(ltpvdw(k).eq.0) then
                
c     tabulated potential
                
                eadd=prmvdw(k,1)
                padd=-prmvdw(k,2)
                
              else if(ltpvdw(k).eq.1) then
                
c     12-6 potential
                
                aaa=prmvdw(k,1)
                bbb=prmvdw(k,2)
                
                eadd=aaa/(9.d0*rcut**9)-bbb/(3.d0*rcut**3)
                padd=12.d0*aaa/(9.d0*rcut**9)-6.d0*bbb/(3.d0*rcut**3)
                
              else if(ltpvdw(k).eq.2) then
                
c     Lennard Jones potential
      
                eps=prmvdw(k,1)
                sig=prmvdw(k,2)
                
                eadd=4.d0*eps*(sig**12/(9.d0*rcut**9)-
     x            sig**6/(3.d0*rcut**3))
                padd=4.d0*eps*(12.d0*sig**12/(9.d0*rcut**9)-
     x            2.d0*sig**6/(rcut**3))
                
              else if(ltpvdw(k).eq.3) then
                
c     n - m potential
                
                eps=prmvdw(k,1)
                ann=prmvdw(k,2)
                amm=prmvdw(k,3)
                rr0=prmvdw(k,4)
                
                eadd=eps/(ann-amm)*(amm*rr0**ann/((ann-3.d0)*
     x            rcut**(ann-3.d0))-ann*rr0**amm/((amm-3.0d0)*
     x            rcut**(amm-3.d0)))
                padd=eps/(ann-amm)*ann*amm*(rr0**ann/((ann-3.d0)*
     x            rcut**(ann-3.d0))-rr0**amm/((amm-3.0d0)*
     x            rcut**(amm-3.d0)))
                
              else if(ltpvdw(k).eq.4) then
                
c     buckingham exp - 6 potential
                
                ccc=prmvdw(k,3)
                
                eadd=-ccc/(3.d0*rcut**3)
                padd=-2.d0*ccc/(rcut**3)
                
              else if(ltpvdw(k).eq.5) then
                
c     born huggins meyer exp -6 - 8  potential
                
                ddd=prmvdw(k,4)
                eee=prmvdw(k,5)
                
                eadd=-ddd/(3.d0*rcut**3)-eee/(5.d0*rcut**5)
                padd=-2.d0*ddd/(rcut**3)-8.d0*eee/(5.d0*rcut**5)
                
              else if(ltpvdw(k).eq.6) then
                
c     hydrogen bond  12 - 10 potential
                
                aaa=prmvdw(k,1)
                bbb=prmvdw(k,2)
                
                eadd=aaa/(9.d0*rcut**9)-bbb/(7.d0*rcut**7)
                padd=12.d0*aaa/(9.d0*rcut**9)-1.d1*bbb/(7.d0*rcut**7)
                
              endif
              
              if(i.ne.j) then
                
                eadd=eadd*2.d0
                padd=padd*2.d0
                
              endif
              
              if(.not.(lfree.or.lghost))then
                
                denprd=twopi*(dble(numtyp(i))*dble(numtyp(j))-
     x            dble(numfrz(i))*dble(numfrz(j)))/volm**2
                
              else
                 
                nafrz=dble(numfrz_fre(i,0)+numfrz_fre(i,1))
                natyp=dble(numtyp_fre(i,0)+numtyp_fre(i,1))
                nbfrz=dble(numfrz_fre(j,0)+numfrz_fre(j,1))
                nbtyp=dble(numtyp_fre(j,0)+numtyp_fre(j,1))
                ncfrz=dble(numfrz_fre(i,0)+numfrz_fre(i,2))
                nctyp=dble(numtyp_fre(i,0)+numtyp_fre(i,2))
                ndfrz=dble(numfrz_fre(j,0)+numfrz_fre(j,2))
                ndtyp=dble(numtyp_fre(j,0)+numtyp_fre(j,2))
                
                if(lghost)then
                  
                  denprd=twopi*(natyp*nbtyp-nafrz*nbfrz)/volm**2
                  denprd3=twopi*(nctyp*ndtyp-ncfrz*ndfrz)/volm**2
                  
                elseif(lfree)then
                  
                  denprd1=twopi*(natyp*nbtyp-nafrz*nbfrz)/volm**2
                  denprd2=twopi*(nctyp*ndtyp-ncfrz*ndfrz)/volm**2
                  denprd=lambda1*denprd1+lambda2*denprd2
                  denprd3=lambda2*denprd1+lambda1*denprd2
                  denprdf=denprd2-denprd1
                  
                endif
                
              endif
              
              elrc=elrc+volm*denprd*eadd
              virlrc=virlrc-denprd*padd*volm
              
              if(lfree.or.lghost)then
                
                elrc2=elrc2+volm*denprd3*eadd
                virlrc2=virlrc2-denprd3*padd*volm
                if(lfree)then
                  elrc_fre=elrc_fre+volm*denprdf*eadd
                  vlrc_fre=vlrc_fre-denprdf*padd*volm
                endif
                
              endif
              
              if(lsolva)then
                
                elrc_sol0(:)=0.d0
                if(lghost)elrc_exc0(:)=0.d0
                
                do it=1,mxtmls
                  
                  do jt=1,mxtmls
                    
                    kt=loc2(it,jt)
                    
                    if(lghost)then
                       
                      natyp=dble(numtyp_sol0(i,it)+numtyp_sol1(i,it))
                      nbtyp=dble(numtyp_sol0(j,jt)+numtyp_sol1(j,jt))
                      nafrz=dble(numfrz_sol0(i,it)+numfrz_sol1(i,it))
                      nbfrz=dble(numfrz_sol0(j,jt)+numfrz_sol1(j,jt))
                      
                      elrc_sol0(kt)=elrc_sol0(kt)+twopi*(natyp*
     x                nbtyp-nafrz*nbfrz)/volm**2
                      
                      nctyp=dble(numtyp_sol0(i,it)+numtyp_sol2(i,it))
                      ndtyp=dble(numtyp_sol0(j,jt)+numtyp_sol2(j,jt))
                      ncfrz=dble(numfrz_sol0(i,it)+numfrz_sol2(i,it))
                      ndfrz=dble(numfrz_sol0(j,jt)+numfrz_sol2(j,jt))
                      
                      elrc_exc0(kt)=elrc_exc0(kt)+twopi*(nctyp*
     x                ndtyp-ncfrz*ndfrz)/volm**2
                      
                    else
                      
                      natyp=dble(numtyp_sol0(i,it))
                      nbtyp=dble(numtyp_sol0(j,jt))
                      nafrz=dble(numfrz_sol0(i,it))
                      nbfrz=dble(numfrz_sol0(j,jt))
                      
                      elrc_sol0(kt)=elrc_sol0(kt)+twopi*(natyp*
     x                nbtyp-nafrz*nbfrz)/volm**2             
                      
                    endif
                    
                  enddo
                  
                enddo
                
                if(lghost)then
                   
                  elrc_sol(:)=elrc_sol(:)+volm*eadd*elrc_sol0(:)
                  elrc_exc(:)=elrc_exc(:)+volm*eadd*elrc_exc0(:)
                  
                else
                  
                  elrc_sol(:)=elrc_sol(:)+volm*eadd*elrc_sol0(:)
                  
                endif
                
              endif
              
            enddo
            
          enddo
          
          if(lfree.or.lghost)then
             
            elrc_sav=elrc
            elrc2_sav=elrc2
            virlrc_sav=virlrc
            virlrc2_sav=virlrc2
            elrc_fre_sav=elrc_fre
            vlrc_fre_sav=vlrc_fre
            
          endif
          
          volm_sav=volm
          
          if(lghost)then
             
            elrc_sol_sav(:)=elrc_sol(:)
            elrc_exc_sav(:)=elrc_exc(:)
            
          elseif(lsolva)then
            
            elrc_sol_sav(:)=elrc_sol(:)
            
          endif
          
        endif
        
      endif
      
      if(idnode.eq.0)then
        
        write(nrite,
     x    "(/,/,'long range correction for: vdw energy  ',e15.6,/,
     x    25x,': vdw pressure',e15.6)")elrc/engunit,
     x    prsunt*virlrc/(-3.d0*volm)
      
        if(lghost)
     x    write(nrite,
     x    "(/,/,'long range correction for: vdw energy  ',e15.6,/,
     x    25x,': vdw pressure',e15.6)")elrc2/engunit,
     x    prsunt*virlrc2/(-3.d0*volm)
        
      endif
      
c     deallocate work arrays
      
      if(lfree.or.lghost)
     x  deallocate (elrc_exc0,numtyp_fre,numfrz_fre,stat=fail(1))
      
      if(lsolva)then
        
        deallocate (elrc_sol0,numtyp_sol0,numfrz_sol0,stat=fail(2))
        
        if(lghost)then
          
          deallocate (numtyp_sol1,numfrz_sol1,stat=fail(3))
          deallocate (numtyp_sol2,numfrz_sol2,stat=fail(4))
          
        endif
        
      endif
      
      return
      end subroutine lrcorrect


      subroutine srfrceneu
     x  (lsolva,lfree,lghost,ik,engsrp,virsrp,dlrpot,rcut)

c***********************************************************************
c     
c     dl_poly subroutine for calculating short range force and
c     potential energy terms using verlet neighbour list
c     neutral group implementation
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith       march 1992
c     
c     neutral groups
c     author    - t. forester    march  1994
c     
c     adapted   - p.-a. cazade oct 2007: solvation, free energy etc
c     
c***********************************************************************
      
      implicit none

      logical lsolva,lfree,lghost,lselect,lskip,idrive,jdrive
      integer ik,m,iatm,jatm,l,k,kkk
      real(8) engsrp,virsrp,dlrpot,rcut,rcsq,fx,fy,fz,omega,omega_exc
      real(8) rrr,ppp,vk0,vk1,vk2,t1,t2,gk0,gk1,gk2,rdlpot,gamma
      real(8) ai,aj,ak,rsq,strs(6),strs_loc(6)
      
      lskip=(lfree.or.lghost)

c     set cutoff condition for pair forces
      
      rcsq=rcut**2

c     reciprocal of interpolation spacing

      rdlpot=1.d0/dlrpot

c     initialise stress tensor accumulators

      strs(:)=0.d0
      strs_loc(:)=0.d0
      
c     initialise potential energy and virial
      
      engsrp=0.d0
      virsrp=0.d0

c     start of primary loop for forces evaluation
      
      do m=1,ik

c     atomic and potential function indices
        
        iatm=ilist(m)
        jatm=jlist(m)
        
c     metadynamics local definitions
        
        if(lmetadyn)then
          
          idrive=driven(ltype(iatm))
          jdrive=driven(ltype(jatm))
          
        endif
        
        if(lskip)then
          if(atm_fre(iatm)*atm_fre(jatm).eq.2)cycle
        endif
        
        aj=ltype(jatm)
        ai=ltype(iatm)

        if(ai.gt.aj) then
          ak=(ai*(ai-1.d0)*0.5d0+aj+0.5d0)
        else
          ak=(aj*(aj-1.d0)*0.5d0+ai+0.5d0)
        endif
        k=lstvdw(int(ak))

        if(abs(vvv(1,k)).gt.1.d-10)then

          rsq=rsqdf(m)

          if(rsq.lt.rcsq) then
              
            rrr=sqrt(rsq)

c     determine interpolation panel for force arrays
            
            l=int(rrr*rdlpot)
            ppp=rrr*rdlpot-dble(l)

c     calculate interaction energy using 3-point interpolation
            
            vk0=vvv(l,k)
            vk1=vvv(l+1,k)
            vk2=vvv(l+2,k)
            t1=vk0+(vk1-vk0)*ppp
            t2=vk1+(vk2-vk1)*(ppp-1.0d0)
            omega=t1+(t2-t1)*ppp*0.5d0

c     calculate forces using 3-point interpolation
            
            gk0=ggg(l,k)
            gk1=ggg(l+1,k)
            gk2=ggg(l+2,k)
            t1=gk0+(gk1-gk0)*ppp
            t2=gk1+(gk2-gk1)*(ppp-1.0d0)
            gamma=(t1+(t2-t1)*ppp*0.5d0)/rsq

c     set selection control
              
            lselect=.true.
            
c     set double index
            
            if(lsolva)kkk=loc2(atmolt(iatm),atmolt(jatm))
            
            if(lghost)then
              
c     selected excitation option
              
              if((atm_fre(iatm).ne.1).and.(atm_fre(jatm).ne.1))then
                
c     reset selection control
                
                lselect=(atm_fre(iatm)+atm_fre(jatm).eq.0)
                
                if(lsolva)vdw_exc(kkk)=vdw_exc(kkk)+omega
                
              endif
              
            elseif(lfree)then
              
c     selected free energy option
              
              if((atm_fre(iatm).eq.1).or.(atm_fre(jatm).eq.1))then
                
c     set hamiltonian mixing parameter
                
                omega=lambda1*omega
                gamma=lambda1*gamma
                
              elseif((atm_fre(iatm).eq.2).or.(atm_fre(jatm).eq.2))then
                
c     set hamiltonian mixing parameter
                
                omega=lambda2*omega
                gamma=lambda2*gamma
                
              endif
              
            endif
            
            if(lselect)then
              
c     calculate potential energy and virial
            
              engsrp=omega+engsrp
              virsrp=virsrp-gamma*rsq
              
              if(lsolva)vdw_sol(kkk)=vdw_sol(kkk)+omega
              
              fx=gamma*xdf(m)
              fy=gamma*ydf(m)
              fz=gamma*zdf(m)
              
              fxx(iatm)=fxx(iatm)+fx
              fyy(iatm)=fyy(iatm)+fy
              fzz(iatm)=fzz(iatm)+fz
              
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
              
            endif

c     metadynamics local parameters
        
            if(lmetadyn.and.(idrive.or.jdrive))then
              
c     local energy and virial
          
              eng_loc=eng_loc+omega
              vir_loc=vir_loc-gamma*rsq
              
c     local forces
          
              fxx_loc(iatm)=fxx_loc(iatm)+fx
              fyy_loc(iatm)=fyy_loc(iatm)+fy
              fzz_loc(iatm)=fzz_loc(iatm)+fz
              
              fxx_loc(jatm)=fxx_loc(jatm)-fx
              fyy_loc(jatm)=fyy_loc(jatm)-fy
              fzz_loc(jatm)=fzz_loc(jatm)-fz
              
c     local stress tensor
          
              strs_loc(1)=strs_loc(1)+xdf(m)*fx
              strs_loc(2)=strs_loc(2)+xdf(m)*fy
              strs_loc(3)=strs_loc(3)+xdf(m)*fz
              strs_loc(4)=strs_loc(4)+ydf(m)*fy
              strs_loc(5)=strs_loc(5)+ydf(m)*fz
              strs_loc(6)=strs_loc(6)+zdf(m)*fz
              
            endif
            
          endif
          
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

      if(lmetadyn)then
        
        stress_loc(1)=stress_loc(1)+strs_loc(1)
        stress_loc(2)=stress_loc(2)+strs_loc(2)
        stress_loc(3)=stress_loc(3)+strs_loc(3)
        stress_loc(4)=stress_loc(4)+strs_loc(2)
        stress_loc(5)=stress_loc(5)+strs_loc(4)
        stress_loc(6)=stress_loc(6)+strs_loc(5)
        stress_loc(7)=stress_loc(7)+strs_loc(3)
        stress_loc(8)=stress_loc(8)+strs_loc(5)
        stress_loc(9)=stress_loc(9)+strs_loc(6)
        
      endif
      
      return
      end subroutine srfrceneu

      end module vdw_module
