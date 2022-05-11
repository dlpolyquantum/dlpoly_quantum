      module vv_pimd_module

c***********************************************************************
c     
c     dl_poly module for path integral velocity verlet integration
c     copyright - daresbury laboratory
c     author    - w. smith    jun 2016
c     
c***********************************************************************
      
      use setup_module,   only : pi,boltz,hbar,mspimd,nrite
      use config_module,  only : stress,buffer,xxx,yyy,zzz,vxx,vyy,vzz
      use pimd_module,    only : zmass,rzmass,etx,ety,etz,pcx,pcy,pcz,
     x                           pxx,pyy,pzz,uxx,uyy,uzz,
     x                           wxx,wyy,wzz,nbeads,nchain,
     x                           unstage_coords,unstage_momenta
      use error_module,   only : error
      use utility_module, only : puni
      
      public pimd_nvt,pimd_nvt_nhc,pimd_nvt_gth1,pimd_nvt_gth2
      
      contains
      
      subroutine pimd_nvt
     x  (lmsite,isw,idnode,mxnode,natms,imcon,ntpmls,tstep,taut,temp,
     x  g_qt4f,engke,engthe)
      
c**********************************************************************
c     
c     dl_poly_classic subroutine for integrating the equations of 
c     motion for path intergral molecular dynamics - using staging 
c     variables thermostated with the nose-hoover thermostat and 
c     integrated with the velocity verlet algorithm
c     
c     copyright - daresbury laboratory
c     author    - w.smith sep 2016
c     
c**********************************************************************
      
      implicit none
      
      logical lmsite
      integer, intent(in) :: imcon,ntpmls
      real(8), intent(in) :: g_qt4f
      integer, intent(in ) :: isw,idnode,mxnode,natms
      real(8), intent(in)  :: tstep,taut,temp
      real(8), intent(out) :: engke,engthe
      
      integer i,k,iatm0,iatm1
      real(8) qqq,qq1,qqk,ppp,hstep,qstep,strkin(9)
      
      iatm0=nbeads*((idnode*natms)/mxnode)
      iatm1=nbeads*(((idnode+1)*natms)/mxnode)
      
c     time step parameters
      
      hstep=tstep/2.d0
      qstep=tstep/4.d0
      
c     mass parameters of thermostats
      
      qq1=boltz*temp*taut**2
      qqk=hbar**2/(boltz*temp*dble(nbeads))
      
c     verlet first part
      
      if(isw.eq.1)then
        
c     update thermostat momenta - 1/4 step
        
        do i=1,iatm1-iatm0
          
          pcx(1,i)=pcx(1,i)+qstep*(pxx(i)**2*rzmass(i)-boltz*temp)
          pcy(1,i)=pcy(1,i)+qstep*(pyy(i)**2*rzmass(i)-boltz*temp)
          pcz(1,i)=pcz(1,i)+qstep*(pzz(i)**2*rzmass(i)-boltz*temp)
          
        enddo
        
        qqq=qq1
        
        do k=1,nbeads
          
c     update thermostats - 1/2 step
          
          do i=k,iatm1-iatm0,nbeads
            
            etx(1,i)=etx(1,i)+(hstep/qqq)*pcx(1,i)
            ety(1,i)=ety(1,i)+(hstep/qqq)*pcy(1,i)
            etz(1,i)=etz(1,i)+(hstep/qqq)*pcz(1,i)
            
          enddo
          
          qqq=qqk
          
        enddo
        
c     update thermostat momenta - 1/4 step
        
        do i=1,iatm1-iatm0
          
          pcx(1,i)=pcx(1,i)+qstep*(pxx(i)**2*rzmass(i)-boltz*temp)
          pcy(1,i)=pcy(1,i)+qstep*(pyy(i)**2*rzmass(i)-boltz*temp)
          pcz(1,i)=pcz(1,i)+qstep*(pzz(i)**2*rzmass(i)-boltz*temp)
          
        enddo
        
c     apply thermostat to bead momenta - 1/2 step
        
        qqq=qq1
        
        do k=1,nbeads
          
          do i=k,iatm1-iatm0,nbeads
            
            pxx(i)=pxx(i)*exp(-(hstep/qqq)*pcx(1,i))
            pyy(i)=pyy(i)*exp(-(hstep/qqq)*pcy(1,i))
            pzz(i)=pzz(i)*exp(-(hstep/qqq)*pcz(1,i))
            
          enddo
          
          qqq=qqk
          
        enddo
        
c     update bead momenta - 1/2 step
        
        do i=1,iatm1-iatm0
          
          pxx(i)=pxx(i)+hstep*wxx(i)
          pyy(i)=pyy(i)+hstep*wyy(i)
          pzz(i)=pzz(i)+hstep*wzz(i)
          
        enddo
        
c     integrate bead positions - 1 step
        
        do i=1,iatm1-iatm0
          
          uxx(i)=uxx(i)+(tstep*rzmass(i))*pxx(i)
          uyy(i)=uyy(i)+(tstep*rzmass(i))*pyy(i)
          uzz(i)=uzz(i)+(tstep*rzmass(i))*pzz(i)
          
        enddo
        
c     unstage coordinates
        
        call unstage_coords(lmsite,idnode,mxnode,natms,imcon,nbeads,
     x       ntpmls,g_qt4f)
        
c     verlet second part
        
      else
        
c     update bead momenta - 1/2 step
        
        do i=1,iatm1-iatm0
          
          pxx(i)=pxx(i)+hstep*wxx(i)
          pyy(i)=pyy(i)+hstep*wyy(i)
          pzz(i)=pzz(i)+hstep*wzz(i)
          
        enddo
        
c     apply thermostat to bead momenta - 1/2 step
        
        qqq=qq1
        
        do k=1,nbeads
          
          do i=k,iatm1-iatm0,nbeads
            
            pxx(i)=pxx(i)*exp(-(hstep/qqq)*pcx(1,i))
            pyy(i)=pyy(i)*exp(-(hstep/qqq)*pcy(1,i))
            pzz(i)=pzz(i)*exp(-(hstep/qqq)*pcz(1,i))
            
          enddo
          
          qqq=qqk
          
        enddo
        
c     update thermostat momenta - 1/4 step

        do i=1,iatm1-iatm0
          
          pcx(1,i)=pcx(1,i)+qstep*(pxx(i)**2*rzmass(i)-boltz*temp)
          pcy(1,i)=pcy(1,i)+qstep*(pyy(i)**2*rzmass(i)-boltz*temp)
          pcz(1,i)=pcz(1,i)+qstep*(pzz(i)**2*rzmass(i)-boltz*temp)
          
        enddo
        
c     update thermostats - 1/2 step
        
        qqq=qq1
        
        do k=1,nbeads
          
          do i=k,iatm1-iatm0,nbeads
            
            etx(1,i)=etx(1,i)+(hstep/qqq)*pcx(1,i)
            ety(1,i)=ety(1,i)+(hstep/qqq)*pcy(1,i)
            etz(1,i)=etz(1,i)+(hstep/qqq)*pcz(1,i)
            
          enddo
          
          qqq=qqk
          
        enddo
        
c     update thermostat momenta - 1/4 step
      
        do i=1,iatm1-iatm0
          
          pcx(1,i)=pcx(1,i)+qstep*(pxx(i)**2*rzmass(i)-boltz*temp)
          pcy(1,i)=pcy(1,i)+qstep*(pyy(i)**2*rzmass(i)-boltz*temp)
          pcz(1,i)=pcz(1,i)+qstep*(pzz(i)**2*rzmass(i)-boltz*temp)
          
        enddo
        
c     calculate thermostat energy
        
        call thermostat_energy
     x    (idnode,mxnode,natms,nbeads,nchain,qq1,qqk,temp,engthe)
        
c     unstage momenta (not strictly necessary)
        
        call unstage_momenta(idnode,mxnode,natms)
        
c     calculate  kinetic tensor and kinetic energy
        
        strkin(:)=0.d0
        
        do i=1,iatm1-iatm0
          
          strkin(1)=strkin(1)+pxx(i)*pxx(i)*rzmass(i)
          strkin(2)=strkin(2)+pxx(i)*pyy(i)*rzmass(i)
          strkin(3)=strkin(3)+pxx(i)*pzz(i)*rzmass(i)
          strkin(5)=strkin(5)+pyy(i)*pyy(i)*rzmass(i)
          strkin(6)=strkin(6)+pyy(i)*pzz(i)*rzmass(i)
          strkin(9)=strkin(9)+pzz(i)*pzz(i)*rzmass(i)
          
        enddo
        
        strkin(4)=strkin(2)
        strkin(7)=strkin(3)
        strkin(8)=strkin(6)
        
        if(mxnode.gt.1)then
          
          buffer(1:9)=strkin(1:9)
          buffer(10)=engthe
          call gdsum(buffer(1),10,buffer(11))
          strkin(1:9)=buffer(1:9)
          engthe=buffer(10)
          
        endif
        
        engke=0.5d0*(strkin(1)+strkin(5)+strkin(9))
        stress(:)=stress(:)+strkin(:)
        
      endif
      
      end subroutine pimd_nvt
      
      subroutine pimd_nvt_gth1
     x  (lmsite,isw,idnode,mxnode,natms,imcon,ntpmls,tstep,taut,temp,
     x  g_qt4f,engke,engthe,chi,uuu)
      
c**********************************************************************
c     
c     dl_poly_classic subroutine for integrating the equations of 
c     motion for path intergral molecular dynamics - using staging 
c     variables thermostated with the gentle thermostat and integrated 
c     with the velocity verlet algorithm
c     
c     references:
c     tuckerman, berne, martyna, klein, j. chem. phys., 99, 2796 (1993)
c     leimkuhler, noorizadeh, thiel, j. stat. phys., 139, 261 (2009)
c     
c     copyright - daresbury laboratory
c     author    - w.smith sep 2016
c     
c**********************************************************************
      
      implicit none
      
      logical noskip,lmsite
      integer, intent(in) :: imcon,ntpmls
      real(8), intent(in) :: g_qt4f
      integer, intent(in ) :: isw,idnode,mxnode,natms
      real(8), intent(in)  :: tstep,taut,temp,chi,uuu(102)
      real(8), intent(out) :: engke,engthe
      
      integer i,k,iatm0,iatm1
      real(8) qqq,qq1,qqk,ppp,hstep,qstep
      real(8) elam,elam1,elamk
      real(8) strkin(9)
      
c     block indices
      
      iatm0=nbeads*((idnode*natms)/mxnode)
      iatm1=nbeads*(((idnode+1)*natms)/mxnode)
      
c     time step parameters
      
      hstep=tstep/2.d0
      qstep=tstep/4.d0
      
c     parameters for thermostats
      
      qq1=boltz*temp*taut**2
      qqk=hbar**2/(boltz*temp*dble(nbeads))
      elam1=exp(-0.25d0*qq1*chi**2/(boltz*temp))
      elamk=exp(-0.25d0*qqk*chi**2/(boltz*temp))
      
c     verlet first part
      
      if(isw.eq.1)then
        
c     update thermostat momenta - 1/2 step
        
        qqq=qq1
        elam=elam1
        noskip=.true.
        
        do k=1,nbeads
          
          do i=k,iatm1-iatm0,nbeads
            
            pcx(1,i)=elam*pcx(1,i)+qqq*chi*gssrnd(noskip,uuu)
            pcy(1,i)=elam*pcy(1,i)+qqq*chi*gssrnd(noskip,uuu)
            pcz(1,i)=elam*pcz(1,i)+qqq*chi*gssrnd(noskip,uuu)
            
          enddo
          
          do i=k,iatm1-iatm0,nbeads
            
           pcx(1,i)=elam*pcx(1,i)+qstep*(pxx(i)**2*rzmass(i)-boltz*temp)
           pcy(1,i)=elam*pcy(1,i)+qstep*(pyy(i)**2*rzmass(i)-boltz*temp)
           pcz(1,i)=elam*pcz(1,i)+qstep*(pzz(i)**2*rzmass(i)-boltz*temp)
            
          enddo
          
c     thermostat the bead momentum
          
          do i=k,iatm1-iatm0,nbeads
            
            pxx(i)=pxx(i)*(1.d0-(hstep/qqq)*pcx(1,i))
            pyy(i)=pyy(i)*(1.d0-(hstep/qqq)*pcy(1,i))
            pzz(i)=pzz(i)*(1.d0-(hstep/qqq)*pcz(1,i))
            
          enddo

          do i=k,iatm1-iatm0,nbeads
            
            etx(1,i)=etx(1,i)+(hstep/qqq)*pcx(1,i)
            ety(1,i)=ety(1,i)+(hstep/qqq)*pcy(1,i)
            etz(1,i)=etz(1,i)+(hstep/qqq)*pcz(1,i)
            
          enddo
          
          do i=k,iatm1-iatm0,nbeads
            
            pcx(1,i)=elam*(pcx(1,i)+qstep*(pxx(i)**2*rzmass(i)-
     x        boltz*temp))
            pcy(1,i)=elam*(pcy(1,i)+qstep*(pyy(i)**2*rzmass(i)-
     x        boltz*temp))
            pcz(1,i)=elam*(pcz(1,i)+qstep*(pzz(i)**2*rzmass(i)-
     x        boltz*temp))
            
          enddo
          
          do i=k,iatm1-iatm0,nbeads
            
            pcx(1,i)=elam*(pcx(1,i)+qqq*chi*gssrnd(noskip,uuu))
            pcy(1,i)=elam*(pcy(1,i)+qqq*chi*gssrnd(noskip,uuu))
            pcz(1,i)=elam*(pcz(1,i)+qqq*chi*gssrnd(noskip,uuu))
            
          enddo
          
          qqq=qqk
          elam=elamk
          
        enddo
        
c     update bead momenta - 1/2 step
        
        do i=1,iatm1-iatm0
          
          pxx(i)=pxx(i)+hstep*wxx(i)
          pyy(i)=pyy(i)+hstep*wyy(i)
          pzz(i)=pzz(i)+hstep*wzz(i)
          
        enddo
        
c     integrate bead positions - full step

        do i=1,iatm1-iatm0
          
          uxx(i)=uxx(i)+(tstep*rzmass(i))*pxx(i)
          uyy(i)=uyy(i)+(tstep*rzmass(i))*pyy(i)
          uzz(i)=uzz(i)+(tstep*rzmass(i))*pzz(i)
          
        enddo
        
c     unstage coordinates
        
        call unstage_coords(lmsite,idnode,mxnode,natms,imcon,nbeads,
     x       ntpmls,g_qt4f)
        
c     verlet second part
        
      else
        
c     update bead momenta - 1/2 step
        
        do i=1,iatm1-iatm0
          
          pxx(i)=pxx(i)+hstep*wxx(i)
          pyy(i)=pyy(i)+hstep*wyy(i)
          pzz(i)=pzz(i)+hstep*wzz(i)
          
        enddo
        
c     update thermostat momenta - 1/2 step
        
        qqq=qq1
        elam=elam1
        
        do k=1,nbeads
          
          do i=k,iatm1-iatm0,nbeads
            
            pcx(1,i)=elam*pcx(1,i)+qqq*chi*gssrnd(noskip,uuu)
            pcy(1,i)=elam*pcy(1,i)+qqq*chi*gssrnd(noskip,uuu)
            pcz(1,i)=elam*pcz(1,i)+qqq*chi*gssrnd(noskip,uuu)
            
          enddo
          
          do i=k,iatm1-iatm0,nbeads
            
           pcx(1,i)=elam*pcx(1,i)+qstep*(pxx(i)**2*rzmass(i)-boltz*temp)
           pcy(1,i)=elam*pcy(1,i)+qstep*(pyy(i)**2*rzmass(i)-boltz*temp)
           pcz(1,i)=elam*pcz(1,i)+qstep*(pzz(i)**2*rzmass(i)-boltz*temp)
            
          enddo
          
          do i=k,iatm1-iatm0,nbeads
            
            etx(1,i)=etx(1,i)+(hstep/qqq)*pcx(1,i)
            ety(1,i)=ety(1,i)+(hstep/qqq)*pcy(1,i)
            etz(1,i)=etz(1,i)+(hstep/qqq)*pcz(1,i)
            
          enddo
          
c     thermostat the bead momentum
          
          do i=k,iatm1-iatm0,nbeads
            
            pxx(i)=(1.d0-(hstep/qqq)*pcx(1,i))*pxx(i)
            pyy(i)=(1.d0-(hstep/qqq)*pcy(1,i))*pyy(i)
            pzz(i)=(1.d0-(hstep/qqq)*pcz(1,i))*pzz(i)
            
          enddo
          
          do i=k,iatm1-iatm0,nbeads
            
            pcx(1,i)=elam*(pcx(1,i)+qstep*(pxx(i)**2*rzmass(i)-
     x        boltz*temp))
            pcy(1,i)=elam*(pcy(1,i)+qstep*(pyy(i)**2*rzmass(i)-
     x       boltz*temp))
            pcz(1,i)=elam*(pcz(1,i)+qstep*(pzz(i)**2*rzmass(i)-
     x        boltz*temp))
            
          enddo
          
          do i=k,iatm1-iatm0,nbeads
            
            pcx(1,i)=elam*(pcx(1,i)+qqq*chi*gssrnd(noskip,uuu))
            pcy(1,i)=elam*(pcy(1,i)+qqq*chi*gssrnd(noskip,uuu))
            pcz(1,i)=elam*(pcz(1,i)+qqq*chi*gssrnd(noskip,uuu))
            
          enddo
          
          qqq=qqk
          elam=elamk
          
        enddo
        
c     calculate thermostat energy
        
        call thermostat_energy
     x    (idnode,mxnode,natms,nbeads,nchain,qq1,qqk,temp,engthe)
        
c     unstage momenta (needed for REVCON file)
        
        call unstage_momenta(idnode,mxnode,natms)
        
c     calculate  kinetic tensor and kinetic energy
        
        strkin(:)=0.d0
        
        do i=1,iatm1-iatm0
          
          strkin(1)=strkin(1)+pxx(i)*pxx(i)*rzmass(i)
          strkin(2)=strkin(2)+pxx(i)*pyy(i)*rzmass(i)
          strkin(3)=strkin(3)+pxx(i)*pzz(i)*rzmass(i)
          strkin(5)=strkin(5)+pyy(i)*pyy(i)*rzmass(i)
          strkin(6)=strkin(6)+pyy(i)*pzz(i)*rzmass(i)
          strkin(9)=strkin(9)+pzz(i)*pzz(i)*rzmass(i)
          
        enddo
        
        strkin(4)=strkin(2)
        strkin(7)=strkin(3)
        strkin(8)=strkin(6)
        
        if(mxnode.gt.1)then
          
          buffer(1:9)=strkin(1:9)
          buffer(10)=engthe
          call gdsum(buffer(1),10,buffer(11))
          strkin(1:9)=buffer(1:9)
          engthe=buffer(10)
          
        endif
        
        engke=0.5d0*(strkin(1)+strkin(5)+strkin(9))
        stress(:)=stress(:)+strkin(:)
        
      endif
      
      end subroutine pimd_nvt_gth1

      subroutine pimd_nvt_gth2
     x  (lmsite,isw,idnode,mxnode,natms,imcon,ntpmls,tstep,taut,temp,
     x  g_qt4f,engke,engthe,chi,uuu)
      
c**********************************************************************
c     
c     dl_poly_classic subroutine for integrating the equations of 
c     motion for path intergral molecular dynamics - using staging 
c     variables thermostated with the gentle thermostat and integrated 
c     with the velocity verlet algorithm
c     
c     references:
c     tuckerman, berne, martyna, klein, j. chem. phys., 99, 2796 (1993)
c     leimkuhler, noorizadeh, thiel, j. stat. phys., 139, 261 (2009)
c     
c     copyright - daresbury laboratory
c     author    - w.smith sep 2016
c     
c**********************************************************************
      
      implicit none
      
      logical noskip,lmsite
      integer, intent(in) :: imcon,ntpmls
      real(8), intent(in) :: g_qt4f
      integer, intent(in ) :: isw,idnode,mxnode,natms
      real(8), intent(in)  :: tstep,taut,temp,chi,uuu(102)
      real(8), intent(out) :: engke,engthe
      
      integer i,k,iatm0,iatm1
      real(8) qqq,qq1,qqk,ppp,hstep,rstep,facm,facp
      real(8) gaus1,gaus2
      real(8) strkin(9)
      
      iatm0=nbeads*((idnode*natms)/mxnode)
      iatm1=nbeads*(((idnode+1)*natms)/mxnode)
      
c     mass parameters of thermostats
      
      qq1=boltz*temp*taut**2
      qqk=hbar**2/(boltz*temp*dble(nbeads))
      
c     time step parameters
      
      hstep=tstep/2.d0
      rstep=chi*sqrt(tstep)
      
c     verlet first part
      
      if(isw.eq.1)then
        
c     integrate bead momenta - 1/2 step
        
        do i=1,iatm1-iatm0
          
          pxx(i)=pxx(i)+hstep*wxx(i)
          pyy(i)=pyy(i)+hstep*wyy(i)
          pzz(i)=pzz(i)+hstep*wzz(i)
          
        enddo
        
c     integrate bead positions - 1/2 step
        
        do i=1,iatm1-iatm0
          
          uxx(i)=uxx(i)+(hstep*rzmass(i))*pxx(i)
          uyy(i)=uyy(i)+(hstep*rzmass(i))*pyy(i)
          uzz(i)=uzz(i)+(hstep*rzmass(i))*pzz(i)
          
        enddo
        
c     thermalise bead momenta - 1/2 step
        
        qqq=qq1
        
        do k=1,nbeads
          
          do i=k,iatm1-iatm0,nbeads
            
            pxx(i)=pxx(i)*exp(-(hstep/qqq)*pcx(1,i))
            pyy(i)=pyy(i)*exp(-(hstep/qqq)*pcy(1,i))
            pzz(i)=pzz(i)*exp(-(hstep/qqq)*pcz(1,i))
            
          enddo
          
          qqq=qqk
          
        enddo
        
c     integrate thermostats - 1/2 step
        
        qqq=qq1
        
        do k=1,nbeads
          
          do i=k,iatm1-iatm0,nbeads
            
            etx(1,i)=etx(1,i)+(hstep/qqq)*pcx(1,i)
            ety(1,i)=ety(1,i)+(hstep/qqq)*pcy(1,i)
            etz(1,i)=etz(1,i)+(hstep/qqq)*pcz(1,i)
            
          enddo

          qqq=qqk

        enddo
        
c     update thermostat momenta - full step

        qqq=qq1
        noskip=.true.
        
        do k=1,nbeads
          
          facm=1.d0-0.25d0*tstep*chi**2/qqq
          facp=1.d0/(1.d0+0.25d0*tstep*chi**2/qqq)
          
          do i=k,iatm1-iatm0,nbeads
            
            pcx(1,i)=facp*(pcx(1,i)*facm+
     x        tstep*(pxx(i)**2*rzmass(i)-boltz*temp)+
     x        rstep*qqq*gssrnd(noskip,uuu))
            pcy(1,i)=facp*(pcy(1,i)*facm+
     x        tstep*(pyy(i)**2*rzmass(i)-boltz*temp)+
     x        rstep*qqq*gssrnd(noskip,uuu))
            pcz(1,i)=facp*(pcz(1,i)*facm+
     x        tstep*(pzz(i)**2*rzmass(i)-boltz*temp)+
     x        rstep*qqq*gssrnd(noskip,uuu))
            
          enddo
          
          qqq=qqk
          
        enddo
        
c     integrate thermostats - 1/2 step
        
        qqq=qq1
        
        do k=1,nbeads
          
          do i=k,iatm1-iatm0,nbeads
            
            etx(1,i)=etx(1,i)+(hstep/qqq)*pcx(1,i)
            ety(1,i)=ety(1,i)+(hstep/qqq)*pcy(1,i)
            etz(1,i)=etz(1,i)+(hstep/qqq)*pcz(1,i)
            
          enddo
          
          qqq=qqk
          
        enddo
        
c     thermalise bead momenta - 1/2 step
        
        qqq=qq1
        
        do k=1,nbeads
          
          do i=k,iatm1-iatm0,nbeads
            
            pxx(i)=pxx(i)*exp(-(hstep/qqq)*pcx(1,i))
            pyy(i)=pyy(i)*exp(-(hstep/qqq)*pcy(1,i))
            pzz(i)=pzz(i)*exp(-(hstep/qqq)*pcz(1,i))
            
          enddo
          
          qqq=qqk
          
        enddo
        
c     integrate bead positions - 1/2 step
          
        do i=1,iatm1-iatm0
          
          uxx(i)=uxx(i)+(hstep*rzmass(i))*pxx(i)
          uyy(i)=uyy(i)+(hstep*rzmass(i))*pyy(i)
          uzz(i)=uzz(i)+(hstep*rzmass(i))*pzz(i)
          
        enddo
        
c     unstage coordinates

        call unstage_coords(lmsite,idnode,mxnode,natms,imcon,nbeads,
     x       ntpmls,g_qt4f)
        
c     verlet second part
        
      else
        
c     integrate momenta - 1/2 step
        
        do i=1,iatm1-iatm0
          
          pxx(i)=pxx(i)+hstep*wxx(i)
          pyy(i)=pyy(i)+hstep*wyy(i)
          pzz(i)=pzz(i)+hstep*wzz(i)
          
        enddo

c     calculate thermostat energy
        
        call thermostat_energy
     x    (idnode,mxnode,natms,nbeads,nchain,qq1,qqk,temp,engthe)
        
c     unstage momenta (needed for REVCON file)
        
        call unstage_momenta(idnode,mxnode,natms)
        
c     calculate  kinetic tensor and kinetic energy
        
        strkin(:)=0.d0
        
        do i=1,iatm1-iatm0
          
          strkin(1)=strkin(1)+pxx(i)*pxx(i)*rzmass(i)
          strkin(2)=strkin(2)+pxx(i)*pyy(i)*rzmass(i)
          strkin(3)=strkin(3)+pxx(i)*pzz(i)*rzmass(i)
          strkin(5)=strkin(5)+pyy(i)*pyy(i)*rzmass(i)
          strkin(6)=strkin(6)+pyy(i)*pzz(i)*rzmass(i)
          strkin(9)=strkin(9)+pzz(i)*pzz(i)*rzmass(i)
          
        enddo
        
        strkin(4)=strkin(2)
        strkin(7)=strkin(3)
        strkin(8)=strkin(6)
        
        if(mxnode.gt.1)then
          
          buffer(1:9)=strkin(1:9)
          buffer(10)=engthe
          call gdsum(buffer(1),10,buffer(11))
          strkin(1:9)=buffer(1:9)
          engthe=buffer(10)
          
        endif
        
        engke=0.5d0*(strkin(1)+strkin(5)+strkin(9))
        stress(:)=stress(:)+strkin(:)
        
      endif
      
      end subroutine pimd_nvt_gth2

      subroutine pimd_nvt_nhc(lmsite,isw,idnode,mxnode,natms,imcon,
     x  ntpmls,tstep,taut,g_qt4f,temp,engke,engthe)
      
c**********************************************************************
c     
c     dl_poly_classic subroutine for integrating the equations of 
c     motion for path intergral molecular dynamics - using staging 
c     variables thermostated with nose-hoover chains and integrated 
c     with the velocity verlet algorithm
c     
c     reference: tuckerman, berne, martyna, klein
c     j. chem. phys. vol. 99 (4) p. 2796
c     
c     copyright - daresbury laboratory
c     author    - w.smith july 2016
c     
c**********************************************************************
      
      implicit none

      logical lmsite
      integer, intent(in) :: imcon,ntpmls
      real(8), intent(in) :: g_qt4f
      integer, intent(in)  :: isw,idnode,mxnode,natms
      real(8), intent(in)  :: tstep,taut,temp
      real(8), intent(out) :: engke,engthe
      
      integer i,j,k,iatm0,iatm1
      real(8) ppp,qqq,qq1,qqk,hstep,qstep,strkin(9)
      
c     block indices
      
      iatm0=nbeads*((idnode*natms)/mxnode)
      iatm1=nbeads*(((idnode+1)*natms)/mxnode)
      
c     thermostat parameters
      
      hstep=tstep/2.d0
      qstep=tstep/4.d0
      qq1=boltz*temp*taut**2
      qqk=hbar**2/(boltz*temp*dble(nbeads))
      
c     verlet first part
      
      if(isw.eq.1)then
        
c     integrate thermostats - 1/4 step
        
        if(nchain.gt.0)call thermo_chain_sy
     x    (idnode,mxnode,natms,qstep,taut,temp)
        
c     integrate bead momenta - 1/2 step
        
        do i=1,iatm1-iatm0
          
          pxx(i)=pxx(i)+hstep*wxx(i)
          pyy(i)=pyy(i)+hstep*wyy(i)
          pzz(i)=pzz(i)+hstep*wzz(i)
          
        enddo
        
c     integrate thermostats - 1/4 step
        
        if(nchain.gt.0)call thermo_chain_sy
     x    (idnode,mxnode,natms,qstep,taut,temp)
        
c     integrate positions (full step)
        
        do i=1,iatm1-iatm0
         
          uxx(i)=uxx(i)+(tstep*rzmass(i))*pxx(i)
          uyy(i)=uyy(i)+(tstep*rzmass(i))*pyy(i)
          uzz(i)=uzz(i)+(tstep*rzmass(i))*pzz(i)
          
        enddo
        
c     unstage coordinates

        call unstage_coords(lmsite,idnode,mxnode,natms,imcon,nbeads,
     x       ntpmls,g_qt4f)
        
c     verlet second part
        
      else
        
c     integrate thermostats - 1/4 step
        
        if(nchain.gt.0)call thermo_chain_sy
     x    (idnode,mxnode,natms,qstep,taut,temp)
        
c     integrate bead momenta - 1/2 step
        
        do i=1,iatm1-iatm0
          
          pxx(i)=pxx(i)+hstep*wxx(i)
          pyy(i)=pyy(i)+hstep*wyy(i)
          pzz(i)=pzz(i)+hstep*wzz(i)
          
        enddo
        
c     integrate thermostats - 1/4 step
        
        if(nchain.gt.0)call thermo_chain_sy
     x    (idnode,mxnode,natms,qstep,taut,temp)
        
c     calculate thermostat energy
        
        call thermostat_energy
     x    (idnode,mxnode,natms,nbeads,nchain,qq1,qqk,temp,engthe)
        
c     unstage momenta (needed for REVCON file)
        
        call unstage_momenta(idnode,mxnode,natms)
        
c     calculate  kinetic tensor and kinetic energy
        
        strkin(:)=0.d0
        
        do i=1,iatm1-iatm0
          
          strkin(1)=strkin(1)+pxx(i)*pxx(i)*rzmass(i)
          strkin(2)=strkin(2)+pxx(i)*pyy(i)*rzmass(i)
          strkin(3)=strkin(3)+pxx(i)*pzz(i)*rzmass(i)
          strkin(5)=strkin(5)+pyy(i)*pyy(i)*rzmass(i)
          strkin(6)=strkin(6)+pyy(i)*pzz(i)*rzmass(i)
          strkin(9)=strkin(9)+pzz(i)*pzz(i)*rzmass(i)
          
        enddo
        
        strkin(4)=strkin(2)
        strkin(7)=strkin(3)
        strkin(8)=strkin(6)
        
        if(mxnode.gt.1)then
          
          buffer(1:9)=strkin(1:9)
          buffer(10)=engthe
          call gdsum(buffer(1),10,buffer(11))
          strkin(1:9)=buffer(1:9)
          engthe=buffer(10)
          
        endif
        
        engke=0.5d0*(strkin(1)+strkin(5)+strkin(9))
        stress(:)=stress(:)+strkin(:)
        
      endif
      
      end subroutine pimd_nvt_nhc
      
      function gssrnd(noskip,uuu)
      
c**********************************************************************
c     
c     subroutine returns a gaussian random number selected from a
c     distribution with zero mean and unit variance
c     
c     box muller method - parallel version
c     
c     copyright - daresbury laboratory
c     author    - w.smith feb 2017
c     
c**********************************************************************

      implicit none

      logical noskip
      real(8) rr0,rr1,rr2,gssrnd,uuu(102)
      save rr1,rr2
      
      if(noskip)then
        
        rr0=puni(0,uuu)
        if(rr0.lt.1.d-15)rr0=puni(0,uuu)
        rr1=sqrt(-2.d0*log(rr0))
        rr2=2.d0*pi*puni(0,uuu)
        gssrnd=rr1*cos(rr2)
        noskip=.false.
        
      else
        
        gssrnd=rr1*sin(rr2)
        noskip=.true.
        
      endif

      return
      end function gssrnd
      
      subroutine thermo_chain_vv
     x  (idnode,mxnode,natms,tstep,taut,temp)
      
c**********************************************************************
c     
c     dl_poly classic routine for updating the thermostat chain momenta in
c     path integral molecular dynamics (velocity verlet version)
c     
c     copyright - daresbury laboratory
c     author    - w.smith aug 2016
c     
c**********************************************************************
      
      implicit none
      
      integer, intent(in) :: idnode,mxnode,natms
      real(8), intent(in) :: tstep,taut,temp
      integer i,j,k,m,n,fail,iatm0,iatm1
      real(8) hstep,qstep
      real(8), allocatable :: qqq(:)

      logical safe
      iatm0=nbeads*((idnode*natms)/mxnode)
      iatm1=nbeads*(((idnode+1)*natms)/mxnode)

c     allocate working arrays
      
      fail=0
      safe=.true.
      
      allocate (qqq(nbeads),stat=fail)
      
      safe=(fail.eq.0)
      if(mxnode.gt.1)call gstate(safe)
      if(.not.safe)call error(idnode,536)
      
c     define scaled time steps
      
      hstep=0.5d0*tstep
      qstep=0.25d0*tstep
      
c     define thermostat masses
      
      qqq(1)=boltz*temp*taut**2
      do k=2,nbeads
        qqq(k)=hbar**2/(boltz*temp*dble(nbeads))
      enddo
      
c     integrate thermostat momenta - 1/2 step
      
      do k=1,nbeads
        
        do i=k,iatm1-iatm0,nbeads
          
          if(nchain.gt.1)then
            
            pcx(nchain,i)=pcx(nchain,i)+
     x        hstep*(pcx(nchain-1,i)**2/qqq(k)-boltz*temp)
            pcy(nchain,i)=pcy(nchain,i)+
     x        hstep*(pcy(nchain-1,i)**2/qqq(k)-boltz*temp)
            pcz(nchain,i)=pcz(nchain,i)+
     x        hstep*(pcz(nchain-1,i)**2/qqq(k)-boltz*temp)

          endif
          
          if(nchain.gt.2)then
            
            do j=nchain-1,2
              
              pcx(j,i)=pcx(j,i)*(1.d0-qstep*pcx(j+1,i)/qqq(k))
              pcy(j,i)=pcy(j,i)*(1.d0-qstep*pcy(j+1,i)/qqq(k))
              pcz(j,i)=pcz(j,i)*(1.d0-qstep*pcz(j+1,i)/qqq(k))

              pcx(j,i)=pcx(j,i)+hstep*(pcx(j-1,i)**2/qqq(k)-boltz*temp)
              pcy(j,i)=pcy(j,i)+hstep*(pcy(j-1,i)**2/qqq(k)-boltz*temp)
              pcz(j,i)=pcz(j,i)+hstep*(pcz(j-1,i)**2/qqq(k)-boltz*temp)
              
              pcx(j,i)=pcx(j,i)*(1.d0-qstep*pcx(j+1,i)/qqq(k))
              pcy(j,i)=pcy(j,i)*(1.d0-qstep*pcy(j+1,i)/qqq(k))
              pcz(j,i)=pcz(j,i)*(1.d0-qstep*pcz(j+1,i)/qqq(k))
              
            enddo

          endif
          
          if(nchain.gt.1)then
            
            pcx(1,i)=pcx(1,i)*(1.d0-qstep*pcx(2,i)/qqq(k))
            pcy(1,i)=pcy(1,i)*(1.d0-qstep*pcy(2,i)/qqq(k))
            pcz(1,i)=pcz(1,i)*(1.d0-qstep*pcz(2,i)/qqq(k))
            
          endif
          
          pcx(1,i)=pcx(1,i)+hstep*(pxx(i)**2*rzmass(i)-boltz*temp)
          pcy(1,i)=pcy(1,i)+hstep*(pyy(i)**2*rzmass(i)-boltz*temp)
          pcz(1,i)=pcz(1,i)+hstep*(pzz(i)**2*rzmass(i)-boltz*temp)
          
          if(nchain.gt.1)then
            
            pcx(1,i)=pcx(1,i)*(1.d0-qstep*pcx(2,i)/qqq(k))
            pcy(1,i)=pcy(1,i)*(1.d0-qstep*pcy(2,i)/qqq(k))
            pcz(1,i)=pcz(1,i)*(1.d0-qstep*pcz(2,i)/qqq(k))
            
          endif
         
        enddo
        
      enddo
      
c     apply thermostats to bead momenta
      
      do k=1,nbeads
        
        do i=k,iatm1-iatm0,nbeads
          
          pxx(i)=pxx(i)*(1.d0-(tstep/qqq(k))*pcx(1,i))
          pyy(i)=pyy(i)*(1.d0-(tstep/qqq(k))*pcy(1,i))
          pzz(i)=pzz(i)*(1.d0-(tstep/qqq(k))*pcz(1,i))
          
        enddo
        
      enddo
      
c     integrate thermostats - full step
      
      do k=1,nbeads
        
        do i=k,iatm1-iatm0,nbeads
          
          do j=1,nchain
            
            etx(j,i)=etx(j,i)+tstep*pcx(j,i)/qqq(k)
            ety(j,i)=ety(j,i)+tstep*pcy(j,i)/qqq(k)
            etz(j,i)=etz(j,i)+tstep*pcz(j,i)/qqq(k)
            
          enddo
          
        enddo
        
      enddo
      
c     integrate thermostat momenta - 1/2 step
      
      do k=1,nbeads
        
        do i=k,iatm1-iatm0,nbeads

          if(nchain.gt.1)then
            
            pcx(1,i)=pcx(1,i)*(1.d0-qstep*pcx(2,i)/qqq(k))
            pcy(1,i)=pcy(1,i)*(1.d0-qstep*pcy(2,i)/qqq(k))
            pcz(1,i)=pcz(1,i)*(1.d0-qstep*pcz(2,i)/qqq(k))
            
          endif
          
          pcx(1,i)=pcx(1,i)+hstep*(pxx(i)**2*rzmass(i)-boltz*temp)
          pcy(1,i)=pcy(1,i)+hstep*(pyy(i)**2*rzmass(i)-boltz*temp)
          pcz(1,i)=pcz(1,i)+hstep*(pzz(i)**2*rzmass(i)-boltz*temp)

          if(nchain.gt.1)then
            
            pcx(1,i)=pcx(1,i)*(1.d0-qstep*pcx(2,i)/qqq(k))
            pcy(1,i)=pcy(1,i)*(1.d0-qstep*pcy(2,i)/qqq(k))
            pcz(1,i)=pcz(1,i)*(1.d0-qstep*pcz(2,i)/qqq(k))

          endif
          
          if(nchain.gt.2)then
            
            do j=2,nchain-1
              
              pcx(j,i)=pcx(j,i)*(1.d0-qstep*pcx(j+1,i)/qqq(k))
              pcy(j,i)=pcy(j,i)*(1.d0-qstep*pcy(j+1,i)/qqq(k))
              pcz(j,i)=pcz(j,i)*(1.d0-qstep*pcz(j+1,i)/qqq(k))
              
              pcx(j,i)=pcx(j,i)+hstep*(pcx(j-1,i)**2/qqq(k)-boltz*temp)
              pcy(j,i)=pcy(j,i)+hstep*(pcy(j-1,i)**2/qqq(k)-boltz*temp)
              pcz(j,i)=pcz(j,i)+hstep*(pcz(j-1,i)**2/qqq(k)-boltz*temp)
              
              pcx(j,i)=pcx(j,i)*(1.d0-qstep*pcx(j+1,i)/qqq(k))
              pcy(j,i)=pcy(j,i)*(1.d0-qstep*pcy(j+1,i)/qqq(k))
              pcz(j,i)=pcz(j,i)*(1.d0-qstep*pcz(j+1,i)/qqq(k))
              
            enddo

          endif
          
          if(nchain.gt.1)then
            
            pcx(nchain,i)=pcx(nchain,i)+
     x        hstep*(pcx(nchain-1,i)**2/qqq(k)-boltz*temp)
            pcy(nchain,i)=pcy(nchain,i)+
     x        hstep*(pcy(nchain-1,i)**2/qqq(k)-boltz*temp)
            pcz(nchain,i)=pcz(nchain,i)+
     x        hstep*(pcz(nchain-1,i)**2/qqq(k)-boltz*temp)
            
          endif
         
        enddo
        
      enddo
      
      deallocate (qqq,stat=fail)
      
      safe=(fail.eq.0)
      if(mxnode.gt.1)call gstate(safe)
      if(.not.safe)call error(idnode,537)
      
      end subroutine thermo_chain_vv

      subroutine thermo_chain_sy
     x  (idnode,mxnode,natms,tstep,taut,temp)
      
c**********************************************************************
c     
c     dl_poly classic routine for updating the thermostat chain momenta in
c     path integral molecular dynamics (suzuki-yoshida version)
c     
c     suzuki-yoshida sixth order scheme
c     
c     copyright - daresbury laboratory
c     author    - w.smith aug 2016
c     
c**********************************************************************
      
      implicit none

      logical safe
      integer, intent(in) :: idnode,mxnode,natms
      real(8), intent(in) :: tstep,taut,temp
      integer i,j,k,m,n,fail,iatm0,iatm1,nsy,ncyc
      real(8) fstep,hstep,qstep,w1,w2,w3,w4,w5,w6,w7
      real(8) wsy(7)
      real(8), allocatable :: qqq(:)

      safe=.true.
      
c     suzuki-yoshida parameters

      nsy=7
      ncyc=10
      wsy(1)=0.784513610477560d0
      wsy(2)=0.235573213359357d0
      wsy(3)=-1.17767998417887d0
      wsy(4)=1.d0-2.d0*(wsy(1)+wsy(2)+wsy(3))
      wsy(5)=-1.17767998417887d0
      wsy(6)=0.235573213359357d0
      wsy(7)=0.784513610477560d0
      
      iatm0=nbeads*((idnode*natms)/mxnode)
      iatm1=nbeads*(((idnode+1)*natms)/mxnode)

c     allocate working arrays
      
      fail=0
      
      allocate (qqq(nbeads),stat=fail)
      
      safe=(fail.eq.0)
      if(mxnode.gt.1)call gstate(safe)
      if(.not.safe)call error(idnode,538)
      
c     define thermostat masses
      
      qqq(1)=boltz*temp*taut**2
      do k=2,nbeads
        qqq(k)=hbar**2/(boltz*temp*dble(nbeads))
      enddo
      
      do m=1,nsy
        
        do n=1,ncyc
          
c     define scaled time steps
          
          fstep=tstep*wsy(m)/dble(ncyc)
          hstep=0.5d0*tstep*wsy(m)/dble(ncyc)
          qstep=0.25d0*tstep*wsy(m)/dble(ncyc)
          
c     integrate thermostat momenta
          
          do k=1,nbeads
            
            do i=k,iatm1-iatm0,nbeads
              
              if(nchain.gt.1)then
                
                pcx(nchain,i)=pcx(nchain,i)+
     x            hstep*(pcx(nchain-1,i)**2/qqq(k)-boltz*temp)
                pcy(nchain,i)=pcy(nchain,i)+
     x            hstep*(pcy(nchain-1,i)**2/qqq(k)-boltz*temp)
                pcz(nchain,i)=pcz(nchain,i)+
     x            hstep*(pcz(nchain-1,i)**2/qqq(k)-boltz*temp)

              endif
              
              if(nchain.gt.2)then
                
                do j=nchain-1,2
                  
                  pcx(j,i)=pcx(j,i)*(1.d0-qstep*pcx(j+1,i)/qqq(k))
                  pcy(j,i)=pcy(j,i)*(1.d0-qstep*pcy(j+1,i)/qqq(k))
                  pcz(j,i)=pcz(j,i)*(1.d0-qstep*pcz(j+1,i)/qqq(k))

                  pcx(j,i)=pcx(j,i)+hstep*(pcx(j-1,i)**2/qqq(k)-
     x              boltz*temp)
                  pcy(j,i)=pcy(j,i)+hstep*(pcy(j-1,i)**2/qqq(k)-
     x              boltz*temp)
                  pcz(j,i)=pcz(j,i)+hstep*(pcz(j-1,i)**2/qqq(k)-
     x              boltz*temp)
                  
                  pcx(j,i)=pcx(j,i)*(1.d0-qstep*pcx(j+1,i)/qqq(k))
                  pcy(j,i)=pcy(j,i)*(1.d0-qstep*pcy(j+1,i)/qqq(k))
                  pcz(j,i)=pcz(j,i)*(1.d0-qstep*pcz(j+1,i)/qqq(k))
                  
                enddo

              endif
              
              if(nchain.gt.1)then
                
                pcx(1,i)=pcx(1,i)*(1.d0-qstep*pcx(2,i)/qqq(k))
                pcy(1,i)=pcy(1,i)*(1.d0-qstep*pcy(2,i)/qqq(k))
                pcz(1,i)=pcz(1,i)*(1.d0-qstep*pcz(2,i)/qqq(k))
                
              endif
              
              pcx(1,i)=pcx(1,i)+hstep*(pxx(i)**2*rzmass(i)-boltz*temp)
              pcy(1,i)=pcy(1,i)+hstep*(pyy(i)**2*rzmass(i)-boltz*temp)
              pcz(1,i)=pcz(1,i)+hstep*(pzz(i)**2*rzmass(i)-boltz*temp)
              
              if(nchain.gt.1)then
                
                pcx(1,i)=pcx(1,i)*(1.d0-qstep*pcx(2,i)/qqq(k))
                pcy(1,i)=pcy(1,i)*(1.d0-qstep*pcy(2,i)/qqq(k))
                pcz(1,i)=pcz(1,i)*(1.d0-qstep*pcz(2,i)/qqq(k))
                
              endif
              
            enddo
            
          enddo
          
c     integrate thermostat variables
          
          do k=1,nbeads
            
            do i=k,iatm1-iatm0,nbeads
              
              do j=1,nchain
                
                etx(j,i)=etx(j,i)+fstep*pcx(j,i)/qqq(k)
                ety(j,i)=ety(j,i)+fstep*pcy(j,i)/qqq(k)
                etz(j,i)=etz(j,i)+fstep*pcz(j,i)/qqq(k)
                
              enddo
              
            enddo
            
          enddo
          
c     apply thermostats to bead momenta
          
          do k=1,nbeads
            
            do i=k,iatm1-iatm0,nbeads
              
              pxx(i)=pxx(i)*(1.d0-(fstep/qqq(k))*pcx(1,i))
              pyy(i)=pyy(i)*(1.d0-(fstep/qqq(k))*pcy(1,i))
              pzz(i)=pzz(i)*(1.d0-(fstep/qqq(k))*pcz(1,i))
              
            enddo
            
          enddo
        
c     integrate thermostat momenta
          
          do k=1,nbeads
            
            do i=k,iatm1-iatm0,nbeads
              
              if(nchain.gt.1)then
                
                pcx(1,i)=pcx(1,i)*(1.d0-qstep*pcx(2,i)/qqq(k))
                pcy(1,i)=pcy(1,i)*(1.d0-qstep*pcy(2,i)/qqq(k))
                pcz(1,i)=pcz(1,i)*(1.d0-qstep*pcz(2,i)/qqq(k))
                
              endif
              
              pcx(1,i)=pcx(1,i)+hstep*(pxx(i)**2*rzmass(i)-boltz*temp)
              pcy(1,i)=pcy(1,i)+hstep*(pyy(i)**2*rzmass(i)-boltz*temp)
              pcz(1,i)=pcz(1,i)+hstep*(pzz(i)**2*rzmass(i)-boltz*temp)
              
              if(nchain.gt.1)then
                
                pcx(1,i)=pcx(1,i)*(1.d0-qstep*pcx(2,i)/qqq(k))
                pcy(1,i)=pcy(1,i)*(1.d0-qstep*pcy(2,i)/qqq(k))
                pcz(1,i)=pcz(1,i)*(1.d0-qstep*pcz(2,i)/qqq(k))
                
              endif
              
              if(nchain.gt.2)then
                
                do j=2,nchain-1
                  
                  pcx(j,i)=pcx(j,i)*(1.d0-qstep*pcx(j+1,i)/qqq(k))
                  pcy(j,i)=pcy(j,i)*(1.d0-qstep*pcy(j+1,i)/qqq(k))
                  pcz(j,i)=pcz(j,i)*(1.d0-qstep*pcz(j+1,i)/qqq(k))
                  
                  pcx(j,i)=pcx(j,i)+hstep*(pcx(j-1,i)**2/qqq(k)-
     x              boltz*temp)
                  pcy(j,i)=pcy(j,i)+hstep*(pcy(j-1,i)**2/qqq(k)-
     x              boltz*temp)
                  pcz(j,i)=pcz(j,i)+hstep*(pcz(j-1,i)**2/qqq(k)-
     x              boltz*temp)
                  
                  pcx(j,i)=pcx(j,i)*(1.d0-qstep*pcx(j+1,i)/qqq(k))
                  pcy(j,i)=pcy(j,i)*(1.d0-qstep*pcy(j+1,i)/qqq(k))
                  pcz(j,i)=pcz(j,i)*(1.d0-qstep*pcz(j+1,i)/qqq(k))
                  
                enddo
                
              endif
              
              if(nchain.gt.1)then
                
                pcx(nchain,i)=pcx(nchain,i)+
     x            hstep*(pcx(nchain-1,i)**2/qqq(k)-boltz*temp)
                pcy(nchain,i)=pcy(nchain,i)+
     x            hstep*(pcy(nchain-1,i)**2/qqq(k)-boltz*temp)
                pcz(nchain,i)=pcz(nchain,i)+
     x            hstep*(pcz(nchain-1,i)**2/qqq(k)-boltz*temp)
                
              endif
              
            enddo
            
          enddo
          
        enddo
        
      enddo
      
      deallocate (qqq,stat=fail)
      
      safe=(fail.eq.0)
      if(mxnode.gt.1)call gstate(safe)
      if(.not.safe)call error(idnode,539)
      
      end subroutine thermo_chain_sy

      subroutine thermostat_energy
     x  (idnode,mxnode,natms,nbeads,nchain,qq1,qqk,temp,engthe)

c**********************************************************************
c     
c     dl_poly_classic subroutine for calculating the total thermostat
c     energy in a path intergral molecular dynamics simulation
c     
c     copyright - daresbury laboratory
c     author    - w.smith dec 2016
c     
c**********************************************************************
      
      implicit none
      
      integer, intent(in ) :: idnode,mxnode,natms,nbeads,nchain
      real(8), intent(in)  :: qq1,qqk,temp
      real(8), intent(out) :: engthe

      integer i,j,k,iatm0,iatm1
      real(8) qqq

      iatm0=nbeads*((idnode*natms)/mxnode)
      iatm1=nbeads*(((idnode+1)*natms)/mxnode)
      
      qqq=qq1
      engthe=0.d0
      
      do k=1,nbeads
        
        do i=k,iatm1-iatm0,nbeads
          
          do j=1,nchain
            
            engthe=engthe+boltz*temp*(etx(j,i)+ety(j,i)+etz(j,i))+
     x        (pcx(j,i)**2+pcy(j,i)**2+pcz(j,i)**2)/(2.d0*qqq)
            
          enddo
          
        enddo
        
        qqq=qqk
        
      enddo

      end subroutine thermostat_energy
            
      end module vv_pimd_module
