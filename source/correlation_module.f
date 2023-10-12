      module correlation_module

c***********************************************************************
c     
c     dl_poly_quantum module for calculation real time correlation
c     functions
c
c     authors - Nathan London and Dil Limbu 2023
c
c***********************************************************************

      use setup_module, only: mspimd,nbeads,corr,mxbuff
      use config_module, only: xxx,yyy,zzz,vxx,vyy,vzz,buffer,cell,
     x                         rcell
      use site_module, only: chgsit,wgtsit,nummols,numsit
      use utility_module, only: sdot0,invert
      use error_module, only: error
c      use merge_tools, only: merge

      implicit none

      real(8), allocatable, save :: xvalinit(:),yvalinit(:),zvalinit(:)
      real(8), allocatable, save :: corrt,xoper(:),yoper(:),zoper(:)

      contains

      subroutine corr_init
     x  (idnode,mxnode,natms,ntpmls,itmols,keyens,keycorr,nummols,
     x  numsit)
      character*6 name
      integer, intent(in) :: idnode,mxnode,natms,ntpmls,itmols,keycorr
      integer, intent(in) :: keyens
      integer, intent(in) :: nummols(ntpmls),numsit(ntpmls)
      integer :: i

c      if(keycorr.eq.1)then
c        allocate(xvalinit(nummols(itmols)))
c        allocate(yvalinit(nummols(itmols)))
c        allocate(zvalinit(nummols(itmols)))
c      else
        allocate(xoper(nummols(itmols)))
        allocate(yoper(nummols(itmols)))
        allocate(zoper(nummols(itmols)))
c      endif

      if(idnode.eq.0)then
         
        if(keycorr.eq.1)then
          name='CORVEL'
          open(unit=corr,file=name,status='replace')
          write(corr,'(A8,I5)') "velocity",nummols(itmols)
        elseif(keycorr.eq.2)then
          name='CORDIP'
          open(unit=corr,file=name,status='replace')
          write(corr,'(A8,I5)') "dipole",nummols(itmols)
        endif
      endif

c      if(keyens.eq.45)then
c        call  calc_val_cent
c     x    (idnode,mxnode,natms,ntpmls,itmols,keycorr,0,nummols,numsit)
c      else 
        call calc_val_bead
     x    (idnode,mxnode,natms,ntpmls,itmols,keycorr,0,nummols,numsit)
c      endif

      if(idnode.eq.0)then
        
c        if(keycorr.eq.1)then   
c          write(corr,'(2e14.6)') 0.d0,corrt
c        else
          write(corr,'(1e14.6)') 0.d0
          do i=1,nummols(itmols)
            write(corr,'(3e14.6)') xoper(i),yoper(i),zoper(i)
          enddo
c        endif
         
      endif
      
      end subroutine corr_init
      
      subroutine correlation
     x  (idnode,mxnode,natms,ntpmls,itmols,keyens,keycorr,nstep,nummols,
     x  numsit,tstep)
      integer, intent(in) :: idnode,mxnode,natms,ntpmls,itmols,nstep
      integer, intent(in) :: keyens,keycorr
      integer, intent(in) :: nummols(ntpmls),numsit(ntpmls)
      real(8), intent(in) :: tstep
      integer :: i

c      if(keyens.eq.45)then
c        call  calc_val_cent
c     x    (idnode,mxnode,natms,ntpmls,itmols,keycorr,nstep,nummols,
c     x    numsit)
c      else 
        call calc_val_bead
     x    (idnode,mxnode,natms,ntpmls,itmols,keycorr,nstep,nummols,
     x    numsit)
c      endif
      
      if(idnode.eq.0)then
         
c        if(keycorr.eq.1)then
c          write(corr,'(2e14.6)') nstep*tstep,corrt
c        else
          write(corr,'(1e14.6)') nstep*tstep
          do i=1,nummols(itmols)
            write(corr,'(3e14.6)') xoper(i),yoper(i),zoper(i)
          enddo
c        endif
         
      endif
      
      end subroutine correlation

      subroutine calc_val_cent
     x  (idnode,mxnode,natms,ntpmls,itmols,keycorr,nstep,nummols,numsit)

      integer, intent(in) :: idnode,mxnode,natms,ntpmls,itmols,nstep
      integer, intent(in) :: keycorr
      integer, intent(in) :: nummols(ntpmls),numsit(ntpmls)
      integer i,j,k,imol0,imol1,jsite,ksite
      real(8) datx(nbeads),daty(nbeads),datz(nbeads)
      real(8) mass(numsit(itmols)),charge(numsit(itmols))
      real(8) atmx(numsit(itmols)),atmy(numsit(itmols))
      real(8) atmz(numsit(itmols))
      real(8) molmass,avgx,avgy,avgz,xval,yval,zval
      real(8) comx,comy,comz

      imol0=(idnode*nummols(itmols))/mxnode+1
      imol1=((idnode+1)*nummols(itmols))/mxnode
      
      jsite=0
      ksite=0
c      write(6,*)"jsite",jsite
      
      if(itmols.gt.1)then
        do i=1,itmols-1
          jsite=jsite+nummols(i)*numsit(i)
          ksite=ksite+numsit(i)
c        write(6,*)"jsite",jsite
        enddo
      endif
c      write(6,*) "nummols",nummols(itmols)
c      write(6,*) "numsit",numsit(itmols)
c      imol0=imol0+jsite
c      imol1=imol1+jsite
c      write(6,*)"imol0",imol0
c      write(6,*)"imol1",imol1
      
      molmass=0.d0
      do i=1,numsit(itmols)
        mass(i)=wgtsit(ksite+i)
        molmass=molmass+mass(i)
      enddo
      
      if(keycorr.eq.2)then
        do i=1,numsit(itmols)
          charge(i)=chgsit(ksite+i)
        enddo
      endif
c      write(6,*)"molmass",molmass
c      write(6,*)"mass",mass(:)

      corrt=0.d0
      do i=imol0,imol1
        xval=0.d0
        yval=0.d0
        zval=0.d0
        do j=1,numsit(itmols)
          do k=1,nbeads
            if(keycorr.eq.1)then 
              datx(k)=vxx((k-1)*natms+(i-1)*numsit(itmols)+jsite+j)
              daty(k)=vyy((k-1)*natms+(i-1)*numsit(itmols)+jsite+j)
              datz(k)=vzz((k-1)*natms+(i-1)*numsit(itmols)+jsite+j)
            elseif(keycorr.eq.2)then
              datx(k)=xxx((k-1)*natms+(i-1)*numsit(itmols)+jsite+j)
              daty(k)=yyy((k-1)*natms+(i-1)*numsit(itmols)+jsite+j)
              datz(k)=zzz((k-1)*natms+(i-1)*numsit(itmols)+jsite+j)
            endif
          enddo

          call bead_average(datx,daty,datz,avgx,avgy,avgz)
          atmx(j)=avgx
          atmy(j)=avgy
          atmz(j)=avgz
        enddo
        if(keycorr.eq.2) call mol_gather(numsit(itmols),atmx,atmy,atmz)
        do j=1,numsit(itmols)
          xval=xval+mass(j)*atmx(j)/molmass
          yval=yval+mass(j)*atmy(j)/molmass
          zval=zval+mass(j)*atmz(j)/molmass
        enddo

        if(keycorr.eq.2)then
          comx=xval
          comy=yval
          comz=zval

          xval=0.d0
          yval=0.d0
          zval=0.d0

          do j=1,numsit(itmols)
            xval=xval+charge(j)*(atmx(j)-comx)
            yval=yval+charge(j)*(atmy(j)-comy)
            zval=zval+charge(j)*(atmz(j)-comz)
          enddo
        
        endif

c        if(nstep.eq.0.and.keycorr.eq.1)then
c          xvalinit(i)=xval
c          yvalinit(i)=yval
c          zvalinit(i)=zval
c        endif
        
c        if(keycorr.eq.1)then
c          corrt=corrt+xval*xvalinit(i)+yval*yvalinit(i)+zval*zvalinit(i)
c        else
          xoper(i)=xval
          yoper(i)=yval
          zoper(i)=zval
c        endif
      enddo
      if(mxnode.gt.1)then
c        if(keycorr.eq.1)then
c          buffer(1)=corrt

c          call gsync()
c          call gdsum(buffer(1),1,buffer(2))

c          corrt=buffer(1)
c        elseif(keycorr.eq.2)then
          call merge
     x      (idnode,mxnode,nummols,mxbuff,xoper,yoper,zoper,buffer)
c        endif
      endif
      
      corrt=corrt/dble(nummols(itmols))
c      write(6,*) "corrt",corrt

      end subroutine calc_val_cent

      subroutine calc_val_bead
     x  (idnode,mxnode,natms,ntpmls,itmols,keycorr,nstep,nummols,numsit)

      integer, intent(in) :: idnode,mxnode,natms,ntpmls,itmols,nstep
      integer, intent(in) :: keycorr
      integer, intent(in) :: nummols(ntpmls),numsit(ntpmls)
      integer i,j,k,imol0,imol1,jsite,ksite
      real(8) datx(nbeads),daty(nbeads),datz(nbeads)
      real(8) mass(numsit(itmols)),charge(numsit(itmols))
      real(8) atmx(numsit(itmols)),atmy(numsit(itmols))
      real(8) atmz(numsit(itmols))
      real(8) molmass,avgx,avgy,avgz,xval,yval,zval
      real(8) comx,comy,comz

      imol0=(idnode*nummols(itmols))/mxnode+1
      imol1=((idnode+1)*nummols(itmols))/mxnode
      
      jsite=0
      ksite=0
c      write(6,*)"jsite",jsite
      if(itmols.gt.1)then
        do i=1,itmols-1
          jsite=jsite+nummols(i)*numsit(i)
          ksite=ksite+numsit(i)
        enddo
      endif
c      write(6,*) "nummols",nummols(itmols)
c      write(6,*) "numsit",numsit(itmols)
c      imol0=imol0+shift
c      imol1=imol1+shift
c      write(6,*)"imol0",idnode,imol0
c      write(6,*)"imol1",idnode,imol1
      
      molmass=0.d0
      do i=1,numsit(itmols)
        mass(i)=wgtsit(ksite+i)
        molmass=molmass+mass(i)
      enddo
      
      if(keycorr.eq.2)then
        do i=1,numsit(itmols)
          charge(i)=chgsit(ksite+i)
        enddo
      endif
c      write(6,*)"molmass",molmass
c      write(6,*)"mass",mass(:)
c      write(6,*)"charge",charge(:)

      corrt=0.d0
c      xoper=0.d0
c      yoper=0.d0
c      zoper=0.d0
      do i=imol0,imol1
      xval=0.d0
        yval=0.d0
        zval=0.d0
        do k=1,nbeads
          datx(k)=0.d0
          daty(k)=0.d0
          datz(k)=0.d0
          do j=1,numsit(itmols)
            if(keycorr.eq.1)then 
              atmx(j)=vxx((k-1)*natms+(i-1)*numsit(itmols)+jsite+j)
              atmy(j)=vyy((k-1)*natms+(i-1)*numsit(itmols)+jsite+j)
              atmz(j)=vzz((k-1)*natms+(i-1)*numsit(itmols)+jsite+j)
            elseif(keycorr.eq.2)then
              atmx(j)=xxx((k-1)*natms+(i-1)*numsit(itmols)+jsite+j)
              atmy(j)=yyy((k-1)*natms+(i-1)*numsit(itmols)+jsite+j)
              atmz(j)=zzz((k-1)*natms+(i-1)*numsit(itmols)+jsite+j)
c              atmx(j)=xxx((k-1)*natms+(i-1)*numsit(itmols)+j)
c              atmy(j)=yyy((k-1)*natms+(i-1)*numsit(itmols)+j)
c              atmz(j)=zzz((k-1)*natms+(i-1)*numsit(itmols)+j)
            endif
          enddo
c          write(6,*) "datx",i,atmx(1)
c          write(6,*) "daty",i,atmx(1)
c          write(6,*) "datz",i,atmx(1)

c         make sure full molecule is in same periodic image
          if(keycorr.eq.2)then
            call mol_gather(numsit(itmols),atmx,atmy,atmz)
          endif
          
c         get center of mass value          
          if(keycorr.eq.1)then
          do j=1,numsit(itmols)
            datx(k)=datx(k)+mass(j)*atmx(j)/molmass
            daty(k)=daty(k)+mass(j)*atmy(j)/molmass
            datz(k)=datz(k)+mass(j)*atmz(j)/molmass
          enddo
          endif

          if(keycorr.eq.2)then
            comx=datx(k)
            comy=daty(k)
            comz=daty(k)

            datx(k)=0.d0
            daty(k)=0.d0
            datz(k)=0.d0

            do j=1,numsit(itmols)
c              datx(k)=datx(k)+charge(j)*(atmx(j))
c              daty(k)=daty(k)+charge(j)*(atmy(j))
c              datz(k)=datz(k)+charge(j)*(atmz(j))
              datx(k)=datx(k)+charge(j)*(atmx(j)-comx)
              daty(k)=daty(k)+charge(j)*(atmy(j)-comy)
              datz(k)=datz(k)+charge(j)*(atmz(j)-comz)
            enddo
          
          endif

        enddo
        
        call bead_average(datx,daty,datz,avgx,avgy,avgz)
        xval=avgx
        yval=avgy
        zval=avgz
c          write(6,*) "avgx",i,j,avgx
c          write(6,*) "avgy",i,j,avgy
c          write(6,*) "avgz",i,j,avgz

c        if(nstep.eq.0.and.keycorr.eq.1)then
c          xvalinit(i)=xval
c          yvalinit(i)=yval
c          zvalinit(i)=zval
c        endif
c        call min_image
c     x    (xval,yval,zval,xvalinit(i),yvalinit(i),zvalinit(i))
c        write(6,*) "mol",i
c        write(6,*) "xval",xval
c        write(6,*) "yval",yval
c        write(6,*) "zval",zval
c        if(keycorr.eq.1)then
c          corrt=corrt+xval*xvalinit(i)+yval*yvalinit(i)+zval*zvalinit(i)
c        else
          xoper(i)=xval
          yoper(i)=yval
          zoper(i)=zval
c          write(6,*) "xoper",xoper(i)
c          xoper=xoper+xval
c          yoper=yoper+yval
c          zoper=zoper+zval
c        endif
      enddo
      if(mxnode.gt.1)then
c        if(keycorr.eq.1)then
c          buffer(1)=corrt

c          call gsync()
c          call gdsum(buffer(1),1,buffer(2))

c          corrt=buffer(1)
c        elseif(keycorr.eq.2)then
          call merge
     x      (idnode,mxnode,nummols(itmols),mxbuff,xoper,yoper,zoper,
     x      buffer)
c       endif

      endif
     
c      write(6,*) xoper(1) 
      corrt=corrt/dble(nummols(itmols))
c      write(6,*) "corrt",corrt

      end subroutine calc_val_bead

      subroutine bead_average(datx,daty,datz,avgx,avgy,avgz)

      implicit none
      integer i
      real(8), intent(in) :: datx(nbeads),daty(nbeads),datz(nbeads)
      real(8) :: avgx,avgy,avgz

      avgx=0.d0
      avgy=0.d0
      avgz=0.d0

      do i=1,nbeads
         avgx=avgx+datx(i)
         avgy=avgy+daty(i)
         avgz=avgz+datz(i)
      enddo
  
      avgx=avgx/dble(nbeads)
      avgy=avgy/dble(nbeads)
      avgz=avgz/dble(nbeads)

      end subroutine bead_average

      subroutine mol_gather(numsit,atmx,atmy,atmz)
        
      integer, intent(in) :: numsit
      integer i
      real(8) :: atmx(numsit),atmy(numsit),atmz(numsit)
      real(8) :: dxx,dyy,dzz,sxx,syy,szz,det

      call invert(cell,rcell,det)

      sxx=atmx(1)*rcell(1)+atmy(1)*rcell(4)+atmz(1)*rcell(7)
      syy=atmx(1)*rcell(2)+atmy(1)*rcell(5)+atmz(1)*rcell(8)
      szz=atmx(1)*rcell(3)+atmy(1)*rcell(6)+atmz(1)*rcell(9)
      sxx=sxx-anint(sxx)
      syy=syy-anint(syy)
      szz=szz-anint(szz)
      atmx(1)=sxx*cell(1)+syy*cell(4)+szz*cell(7)
      atmy(1)=sxx*cell(2)+syy*cell(5)+szz*cell(8)
      atmz(1)=sxx*cell(3)+syy*cell(6)+szz*cell(9)

      do i=2,numsit
        
        dxx=atmx(i)-atmx(1)
        dyy=atmy(i)-atmy(1)
        dzz=atmz(i)-atmz(1)
        sxx=dxx*rcell(1)+dyy*rcell(4)+dzz*rcell(7)
        syy=dxx*rcell(2)+dyy*rcell(5)+dzz*rcell(8)
        szz=dxx*rcell(3)+dyy*rcell(6)+dzz*rcell(9)
        sxx=sxx-anint(sxx)
        syy=syy-anint(syy)
        szz=szz-anint(szz)
        dxx=sxx*cell(1)+syy*cell(4)+szz*cell(7)
        dyy=sxx*cell(2)+syy*cell(5)+szz*cell(8)
        dzz=sxx*cell(3)+syy*cell(6)+szz*cell(9)
        atmx(i)=atmx(1)+dxx
        atmy(i)=atmy(1)+dyy
        atmz(i)=atmz(1)+dzz

      enddo
       
      end subroutine mol_gather
      
      subroutine min_image(xval,yval,zval,xvalinit,yvalinit,zvalinit)
        
      real(8) :: xval,yval,zval,xvalinit,yvalinit,zvalinit
      real(8) :: dxx,dyy,dzz,sxx,syy,szz,det

      call invert(cell,rcell,det)

      dxx=xval-xvalinit
      dyy=yval-yvalinit
      dzz=zval-zvalinit
      sxx=dxx*rcell(1)+dyy*rcell(4)+dzz*rcell(7)
      syy=dxx*rcell(2)+dyy*rcell(5)+dzz*rcell(8)
      szz=dxx*rcell(3)+dyy*rcell(6)+dzz*rcell(9)
      sxx=sxx-anint(sxx)
      syy=syy-anint(syy)
      szz=szz-anint(szz)
      dxx=sxx*cell(1)+syy*cell(4)+szz*cell(7)
      dyy=sxx*cell(2)+syy*cell(5)+szz*cell(8)
      dzz=sxx*cell(3)+syy*cell(6)+szz*cell(9)
      xval=xvalinit+dxx
      yval=yvalinit+dxx
      zval=zvalinit+dxx
       
      end subroutine min_image
      
      end module correlation_module
