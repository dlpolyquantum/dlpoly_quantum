      module vv_pimd_module

c***********************************************************************
c     
c     dl_poly module for path integral velocity verlet integration
c     copyright - daresbury laboratory
c     author    - w. smith    jun 2016
c     
c***********************************************************************
      
      use setup_module,   only : pi,boltz,hbar,mspimd,nrite
      use config_module,  only : stress,buffer,xxx,yyy,zzz,vxx,vyy,vzz,
     x                           weight,fxx,fyy,fzz
      use pimd_module,    only : zmass,rzmass,etx,ety,etz,pcx,pcy,pcz,
     x                           pxx,pyy,pzz,uxx,uyy,uzz,
     x                           wxx,wyy,wzz,nbeads,nchain,
     x                           unstage_coords,unstage_momenta,
     x                           norm2coord,norm2momenta,
     x                           coord2norm,momenta2norm,
     x                           freerp,pileC1,pileC2,nmfreq,pmerge,
     x                           freerp_noc,cvec1,cvec2
      use error_module,   only : error
      use utility_module, only : puni
      use pimd_piglet_module, only : piglet_thermo_step

ccc npt variable       
      use config_module, only:cell,strcns
      use ensemble_tools_module,only:nhc_baro
c      use error_module
c      use metafreeze_module,only : lmetadyn
c      use property_module
      use setup_module, only:mxatyp,mxbuff
c      use shake_module
      use site_module, only:dens
c      use utility_module, only:images
      use nhc_module, only:ksi,pksi,eta_nhc,peta,v_epsilon
c      use water_module, only:qtip4pf,qt4_force_redist
     
      real(8) :: press2 
      
      public pimd_nvt,pimd_nvt_nhc,pimd_nvt_gth1,pimd_nvt_gth2,press2
      
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
        
c        if(nchain.gt.0)call thermo_chain_sy
c     x    (idnode,mxnode,natms,qstep,taut,temp)
        
c     integrate bead momenta - 1/2 step
        
        do i=1,iatm1-iatm0
          
          pxx(i)=pxx(i)+hstep*wxx(i)
          pyy(i)=pyy(i)+hstep*wyy(i)
          pzz(i)=pzz(i)+hstep*wzz(i)
          
        enddo
        
c     integrate thermostats - 1/4 step
        
c        if(nchain.gt.0)call thermo_chain_sy
c     x    (idnode,mxnode,natms,qstep,taut,temp)
        
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
        
c        if(nchain.gt.0)call thermo_chain_sy
c     x    (idnode,mxnode,natms,qstep,taut,temp)
        
c     integrate bead momenta - 1/2 step
        
        do i=1,iatm1-iatm0
          
          pxx(i)=pxx(i)+hstep*wxx(i)
          pyy(i)=pyy(i)+hstep*wyy(i)
          pzz(i)=pzz(i)+hstep*wzz(i)
          
        enddo
        
c     integrate thermostats - 1/4 step
        
c        if(nchain.gt.0)call thermo_chain_sy
c     x    (idnode,mxnode,natms,qstep,taut,temp)
        
c     calculate thermostat energy
        
c        call thermostat_energy
c     x    (idnode,mxnode,natms,nbeads,nchain,qq1,qqk,temp,engthe)
        
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
      
      subroutine pimd_nvt_nhc_nm(lmsite,isw,idnode,mxnode,natms,imcon,
     x  ntpmls,tstep,taut,g_qt4f,temp,engke,engthe)
      
c**********************************************************************
c     
c     dl_poly_quantum subroutine for integrating the equations of 
c     motion for path intergral molecular dynamics - using normal mode 
c     variables thermostated with nose-hoover chains and integrated 
c     with the velocity verlet algorithm
c     
c     reference: tuckerman, berne, martyna, klein
c     j. chem. phys. vol. 99 (4) p. 2796
c     
c     copyright - 
c     authors   - Dil Limbu and Nathan London 2023
c     
c**********************************************************************
      
      implicit none

      logical lmsite
      integer, intent(in) :: imcon,ntpmls
      real(8), intent(in) :: g_qt4f
      integer, intent(in)  :: isw,idnode,mxnode,natms
      real(8), intent(in)  :: tstep,taut,temp
      real(8), intent(out) :: engke,engthe
      
      integer i,j,k,iatm0,iatm1,nstart
      real(8) ppp,qqq,qq1,qqk,hstep,qstep,strkin(9)
      real(8) :: mass(nbeads),rmass(nbeads)

c     block indices
      
      iatm0=nbeads*((idnode*natms)/mxnode)
      iatm1=nbeads*(((idnode+1)*natms)/mxnode)
      
c     thermostat parameters
      
      hstep=tstep/2.d0
      qstep=tstep/4.d0
      qq1=boltz*temp*taut**2
      qqk=hbar**2/(boltz*temp*dble(nbeads))
      nstart=1

c     verlet first part

      if(isw.eq.1)then
       
c      write(6,*) "vv x pos"
c      write(6,*) uxx(:) 
c      write(6,*) "vv y pos"
c      write(6,*) uyy(:) 
c      write(6,*) "vv z pos"
c      write(6,*) uzz(:) 
c      write(6,*) "vv mom"
c      write(6,*) pxx(:) 
c     integrate thermostats - 1/4 step
        
        if(nchain.gt.0)call thermo_chain_nm
     x    (idnode,mxnode,natms,nstart,qstep,taut,temp)
        
c      write(6,*) "vv mom after thermo"
c      write(6,*) pxx(:) 
c      write(6,*) "vv force"
c      write(6,*) wxx(:) 
c     integrate bead momenta - 1/2 step
        
        do i=1,iatm1-iatm0
          
          pxx(i)=pxx(i)+hstep*wxx(i)
          pyy(i)=pyy(i)+hstep*wyy(i)
          pzz(i)=pzz(i)+hstep*wzz(i)
          
        enddo
c      write(6,*) "vv mom after force"
c      write(6,*) pxx(:) 
        
c     integrate thermostats - 1/4 step
        
        if(nchain.gt.0)call thermo_chain_nm
     x    (idnode,mxnode,natms,nstart,qstep,taut,temp)
        
c      write(6,*) "vv mom after thermo 2"
c      write(6,*) pxx(:) 
c     integrate free ring polymer positions (full step) 
        
        do i=1,iatm1-iatm0,nbeads
          
          j=i-1+nbeads

          mass(1:nbeads)=zmass(i:j)
          rmass(1:nbeads)=rzmass(i:j)

          call freerp(pxx(i:j),uxx(i:j),mass,rmass,tstep,temp)
          call freerp(pyy(i:j),uyy(i:j),mass,rmass,tstep,temp)
          call freerp(pzz(i:j),uzz(i:j),mass,rmass,tstep,temp)
          
        enddo
c      write(6,*) "vv x pos after free"
c      write(6,*) uxx(:) 
c      write(6,*) "vv y pos after free"
c      write(6,*) uyy(:) 
c      write(6,*) "vv z pos after free"
c      write(6,*) uzz(:) 
c      write(6,*) "vv mom after free"
c      write(6,*) pxx(:) 
        
c        do i=1,(iatm1-iatm0)/nbeads

c          mass(1:nbeads)=zmass((i-1)*nbeads+1:i*nbeads)
c          rmass(1:nbeads)=rzmass((i-1)*nbeads+1:i*nbeads)

c          call freerp (pxx((i-1)*nbeads+1:i*nbeads),
c     x            uxx((i-1)*nbeads+1:i*nbeads),mass,rmass,tstep,temp)

c          call freerp (pyy((i-1)*nbeads+1:i*nbeads),
c     x            uyy((i-1)*nbeads+1:i*nbeads),mass,rmass,tstep,temp)

c          call freerp (pzz((i-1)*nbeads+1:i*nbeads),
c     x            uzz((i-1)*nbeads+1:i*nbeads),mass,rmass,tstep,temp)

c        enddo
        
c     unstage coordinates

        call norm2coord(lmsite,idnode,mxnode,natms,imcon,nbeads,
     x       ntpmls,g_qt4f)

c        call unstage_coords(lmsite,idnode,mxnode,natms,imcon,nbeads,
c     x       ntpmls,g_qt4f)
        
c     verlet second part
        
      else
        
c     integrate thermostats - 1/4 step
        
        if(nchain.gt.0)call thermo_chain_nm
     x    (idnode,mxnode,natms,nstart,qstep,taut,temp)
        
c     integrate bead momenta - 1/2 step
        
        do i=1,iatm1-iatm0
          
          pxx(i)=pxx(i)+hstep*wxx(i)
          pyy(i)=pyy(i)+hstep*wyy(i)
          pzz(i)=pzz(i)+hstep*wzz(i)
          
        enddo
        
c     integrate thermostats - 1/4 step
        
        if(nchain.gt.0)call thermo_chain_nm
     x    (idnode,mxnode,natms,nstart,qstep,taut,temp)
        
c     calculate thermostat energy
        
        call thermostat_energy_nm
     x    (idnode,mxnode,natms,nbeads,nchain,nstart,qq1,qqk,taut,temp,
     x    engthe)
        
c     unstage momenta (needed for REVCON file)
        
        call norm2momenta(idnode,mxnode,natms)
        
c        call unstage_momenta(idnode,mxnode,natms)
        
c     calculate  kinetic tensor and kinetic energy
        
        strkin(:)=0.d0
        
        do i=1,iatm1-iatm0
          
          strkin(1)=strkin(1)+pxx(i)*pxx(i)*rzmass(i)
          strkin(2)=strkin(2)+pxx(i)*pyy(i)*rzmass(i)
          strkin(3)=strkin(3)+pxx(i)*pzz(i)*rzmass(i)
          strkin(5)=strkin(5)+pyy(i)*pyy(i)*rzmass(i)
          strkin(6)=strkin(6)+pyy(i)*pzz(i)*rzmass(i)
          strkin(9)=strkin(9)+pzz(i)*pzz(i)*rzmass(i)
          
c          strkin(1)=strkin(1)+vxx(i)*vxx(i)*weight(i)
c          strkin(2)=strkin(2)+vxx(i)*vyy(i)*weight(i)
c          strkin(3)=strkin(3)+vxx(i)*vzz(i)*weight(i)
c          strkin(5)=strkin(5)+vyy(i)*vyy(i)*weight(i)
c          strkin(6)=strkin(6)+vyy(i)*vzz(i)*weight(i)
c          strkin(9)=strkin(9)+vzz(i)*vzz(i)*weight(i)
          
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
      
      end subroutine pimd_nvt_nhc_nm
      
      subroutine pimd_nve(lmsite,isw,idnode,mxnode,natms,imcon,
     x  ntpmls,tstep,g_qt4f,temp,engke)
      
c**********************************************************************
c     
c     dl_poly_quantum subroutine for integrating the equations of 
c     motion for ring polymer molecular dynamics - using normal mode 
c     variables 
c     
c     reference: tuckerman, berne, martyna, klein
c     j. chem. phys. vol. 99 (4) p. 2796
c     
c     copyright - 
c     authors   - Dil Limbu and Nathan London 2023
c     
c**********************************************************************
      
      implicit none

      logical lmsite
      integer, intent(in) :: imcon,ntpmls
      real(8), intent(in) :: g_qt4f
      integer, intent(in)  :: isw,idnode,mxnode,natms
      real(8), intent(in)  :: tstep,temp
      real(8), intent(out) :: engke
      
      integer i,j,k,iatm0,iatm1
      real(8) hstep,strkin(9)

      real(8) :: mass(nbeads),rmass(nbeads)

c     block indices
      
      iatm0=nbeads*((idnode*natms)/mxnode)
      iatm1=nbeads*(((idnode+1)*natms)/mxnode)
      
c     thermostat parameters
      
      hstep=tstep/2.d0
      
c     verlet first part
      
      if(isw.eq.1)then
       
c      write(6,*) "vv x pos"
c      write(6,*) uxx(:) 
c      write(6,*) "vv y pos"
c      write(6,*) uyy(:) 
c      write(6,*) "vv z pos"
c      write(6,*) uzz(:) 
c      write(6,*) "vv mom"
c      write(6,*) pxx(:) 
        
c      write(6,*) "vv mom after thermo"
c      write(6,*) pxx(:) 
c      write(6,*) "vv force"
c      write(6,*) wxx(:) 
c     integrate bead momenta - 1/2 step
        
        do i=1,iatm1-iatm0
          
          pxx(i)=pxx(i)+hstep*wxx(i)
          pyy(i)=pyy(i)+hstep*wyy(i)
          pzz(i)=pzz(i)+hstep*wzz(i)
          
        enddo
c      write(6,*) "vv mom after force"
c      write(6,*) pxx(:) 
        
c      write(6,*) "vv mom after thermo 2"
c      write(6,*) pxx(:) 
c     integrate free ring polymer positions (full step) 
        
        do i=1,(iatm1-iatm0)/nbeads
c          write(6,*) "free RP atom",i
          mass(1:nbeads)=zmass((i-1)*nbeads+1:i*nbeads)
          rmass(1:nbeads)=rzmass((i-1)*nbeads+1:i*nbeads)

          call freerp (pxx((i-1)*nbeads+1:i*nbeads),
     x            uxx((i-1)*nbeads+1:i*nbeads),mass,rmass,tstep,temp)
          call freerp (pyy((i-1)*nbeads+1:i*nbeads),
     x            uyy((i-1)*nbeads+1:i*nbeads),mass,rmass,tstep,temp)
          call freerp (pzz((i-1)*nbeads+1:i*nbeads),
     x            uzz((i-1)*nbeads+1:i*nbeads),mass,rmass,tstep,temp)

c          uxx(i)=uxx(i)+(tstep*rzmass(i))*pxx(i)
c          uyy(i)=uyy(i)+(tstep*rzmass(i))*pyy(i)
c          uzz(i)=uzz(i)+(tstep*rzmass(i))*pzz(i)
          
        enddo
c      write(6,*) "vv x pos after free"
c      write(6,*) uxx(:) 
c      write(6,*) "vv y pos after free"
c      write(6,*) uyy(:) 
c      write(6,*) "vv z pos after free"
c      write(6,*) uzz(:) 
c      write(6,*) "vv mom after free"
c      write(6,*) pxx(:) 
        
c     unstage coordinates

      call norm2coord(lmsite,idnode,mxnode,natms,imcon,nbeads,
     x       ntpmls,g_qt4f)

c        call unstage_coords(lmsite,idnode,mxnode,natms,imcon,nbeads,
c     x       ntpmls,g_qt4f)
        
c     verlet second part
        
      else
        
c     integrate bead momenta - 1/2 step
        
        do i=1,iatm1-iatm0
          
          pxx(i)=pxx(i)+hstep*wxx(i)
          pyy(i)=pyy(i)+hstep*wyy(i)
          pzz(i)=pzz(i)+hstep*wzz(i)
          
        enddo
        
c     unstage momenta (needed for REVCON file)
        
        call norm2momenta(idnode,mxnode,natms)
        
c        call unstage_momenta(idnode,mxnode,natms)
        
c     calculate  kinetic tensor and kinetic energy
        
        strkin(:)=0.d0
        
        do i=1,iatm1-iatm0
          
          strkin(1)=strkin(1)+pxx(i)*pxx(i)*rzmass(i)
          strkin(2)=strkin(2)+pxx(i)*pyy(i)*rzmass(i)
          strkin(3)=strkin(3)+pxx(i)*pzz(i)*rzmass(i)
          strkin(5)=strkin(5)+pyy(i)*pyy(i)*rzmass(i)
          strkin(6)=strkin(6)+pyy(i)*pzz(i)*rzmass(i)
          strkin(9)=strkin(9)+pzz(i)*pzz(i)*rzmass(i)
          
c          strkin(1)=strkin(1)+vxx(i)*vxx(i)*weight(i)
c          strkin(2)=strkin(2)+vxx(i)*vyy(i)*weight(i)
c          strkin(3)=strkin(3)+vxx(i)*vzz(i)*weight(i)
c          strkin(5)=strkin(5)+vyy(i)*vyy(i)*weight(i)
c          strkin(6)=strkin(6)+vyy(i)*vzz(i)*weight(i)
c          strkin(9)=strkin(9)+vzz(i)*vzz(i)*weight(i)
          
        enddo
        
        strkin(4)=strkin(2)
        strkin(7)=strkin(3)
        strkin(8)=strkin(6)
        
        if(mxnode.gt.1)then
          
          buffer(1:9)=strkin(1:9)
          call gdsum(buffer(1),10,buffer(11))
          strkin(1:9)=buffer(1:9)
          
        endif
        
        engke=0.5d0*(strkin(1)+strkin(5)+strkin(9))
        stress(:)=stress(:)+strkin(:)
        
      endif
      
      end subroutine pimd_nve
      
      subroutine pacmd(lmsite,isw,idnode,mxnode,natms,imcon,
     x  ntpmls,tstep,taut,g_qt4f,temp,engke,engthe)
      
c**********************************************************************
c     
c     dl_poly_quantum subroutine for integrating the equations of 
c     motion for path intergral molecular dynamics - using normal mode 
c     variables thermostated with nose-hoover chains and integrated 
c     with the velocity verlet algorithm
c     
c     reference: tuckerman, berne, martyna, klein
c     j. chem. phys. vol. 99 (4) p. 2796
c     
c     copyright - 
c     authors   - Dil Limbu and Nathan London 2023
c     
c**********************************************************************
      
      implicit none

      logical lmsite
      integer, intent(in) :: imcon,ntpmls
      real(8), intent(in) :: g_qt4f
      integer, intent(in)  :: isw,idnode,mxnode,natms
      real(8), intent(in)  :: tstep,taut,temp
      real(8), intent(out) :: engke,engthe
      
      integer i,j,k,iatm0,iatm1,nstart
      real(8) ppp,qqq,qq1,qqk,hstep,qstep,strkin(9)

      real(8) :: mass(nbeads),rmass(nbeads)

c     block indices
      
      iatm0=nbeads*((idnode*natms)/mxnode)
      iatm1=nbeads*(((idnode+1)*natms)/mxnode)
      
c     thermostat parameters
      
      hstep=tstep/2.d0
      qstep=tstep/4.d0
      qq1=boltz*temp*taut**2
      qqk=hbar**2/(boltz*temp*dble(nbeads))
      nstart=2

c     verlet first part
      
      if(isw.eq.1)then
       
c      write(6,*) "vv x pos"
c      write(6,*) uxx(:) 
c      write(6,*) "vv y pos"
c      write(6,*) uyy(:) 
c      write(6,*) "vv z pos"
c      write(6,*) uzz(:)
c      write(6,*) "mass"
c      write(6,*) zmass(:) 
c      write(6,*) "vv mom"
c      write(6,*) pxx(:) 
c     integrate thermostats - 1/4 step
        
        if(nchain.gt.0)call thermo_chain_nm
     x    (idnode,mxnode,natms,nstart,qstep,taut,temp)
        
c      write(6,*) "vv mom after thermo"
c      write(6,*) pxx(:) 
c      write(6,*) "vv force"
c      write(6,*) wxx(:) 
c     integrate bead momenta - 1/2 step
        
        do i=1,iatm1-iatm0
          
          pxx(i)=pxx(i)+hstep*wxx(i)
          pyy(i)=pyy(i)+hstep*wyy(i)
          pzz(i)=pzz(i)+hstep*wzz(i)
          
        enddo
c      write(6,*) "vv mom after force"
c      write(6,*) pxx(:) 
        
c     integrate thermostats - 1/4 step
        
        if(nchain.gt.0)call thermo_chain_nm
     x    (idnode,mxnode,natms,nstart,qstep,taut,temp)
        
c      write(6,*) "vv mom after thermo 2"
c      write(6,*) pxx(:) 
c     integrate free ring polymer positions (full step) 
        
        do i=1,(iatm1-iatm0)/nbeads
c          write(6,*) "free RP atom",i
          mass(1:nbeads)=zmass((i-1)*nbeads+1:i*nbeads)
          rmass(1:nbeads)=rzmass((i-1)*nbeads+1:i*nbeads)

          call freerp (pxx((i-1)*nbeads+1:i*nbeads),
     x            uxx((i-1)*nbeads+1:i*nbeads),mass,rmass,tstep,temp)
          call freerp (pyy((i-1)*nbeads+1:i*nbeads),
     x            uyy((i-1)*nbeads+1:i*nbeads),mass,rmass,tstep,temp)
          call freerp (pzz((i-1)*nbeads+1:i*nbeads),
     x            uzz((i-1)*nbeads+1:i*nbeads),mass,rmass,tstep,temp)

c          uxx(i)=uxx(i)+(tstep*rzmass(i))*pxx(i)
c          uyy(i)=uyy(i)+(tstep*rzmass(i))*pyy(i)
c          uzz(i)=uzz(i)+(tstep*rzmass(i))*pzz(i)
          
        enddo
c      write(6,*) "vv x pos after free"
c      write(6,*) uxx(:) 
c      write(6,*) "vv y pos after free"
c      write(6,*) uyy(:) 
c      write(6,*) "vv z pos after free"
c      write(6,*) uzz(:) 
c      write(6,*) "vv mom after free"
c      write(6,*) pxx(:) 
        
c     unstage coordinates

      call norm2coord(lmsite,idnode,mxnode,natms,imcon,nbeads,
     x       ntpmls,g_qt4f)

c        call unstage_coords(lmsite,idnode,mxnode,natms,imcon,nbeads,
c     x       ntpmls,g_qt4f)
        
c     verlet second part
        
      else
        
c     integrate thermostats - 1/4 step
        
        if(nchain.gt.0)call thermo_chain_nm
     x    (idnode,mxnode,natms,nstart,qstep,taut,temp)
        
c      write(6,*) "vv mom after thermo 3"
c       write(6,*) pxx(:) 
c     integrate bead momenta - 1/2 step
        
        do i=1,iatm1-iatm0
          
          pxx(i)=pxx(i)+hstep*wxx(i)
          pyy(i)=pyy(i)+hstep*wyy(i)
          pzz(i)=pzz(i)+hstep*wzz(i)
          
        enddo
c     write(6,*) "vv mom after force 2"
c      write(6,*) pxx(:) 
        
c     integrate thermostats - 1/4 step
        
        if(nchain.gt.0)call thermo_chain_nm
     x    (idnode,mxnode,natms,nstart,qstep,taut,temp)
        
c      write(6,*) "vv mom after thermo 4"
c      write(6,*) pxx(:) 
c      write(6,*) ""
c     calculate thermostat energy
        
        call thermostat_energy_nm
     x    (idnode,mxnode,natms,nbeads,nchain,nstart,qq1,qqk,taut,temp,
     x    engthe)
        
c     unstage momenta (needed for REVCON file)
        
        call norm2momenta(idnode,mxnode,natms)
        
c        call unstage_momenta(idnode,mxnode,natms)
        
c     calculate  kinetic tensor and kinetic energy
        
        strkin(:)=0.d0
        
        do i=1,iatm1-iatm0
          
          strkin(1)=strkin(1)+pxx(i)*pxx(i)*rzmass(i)
          strkin(2)=strkin(2)+pxx(i)*pyy(i)*rzmass(i)
          strkin(3)=strkin(3)+pxx(i)*pzz(i)*rzmass(i)
          strkin(5)=strkin(5)+pyy(i)*pyy(i)*rzmass(i)
          strkin(6)=strkin(6)+pyy(i)*pzz(i)*rzmass(i)
          strkin(9)=strkin(9)+pzz(i)*pzz(i)*rzmass(i)
          
c          strkin(1)=strkin(1)+vxx(i)*vxx(i)*weight(i)
c          strkin(2)=strkin(2)+vxx(i)*vyy(i)*weight(i)
c          strkin(3)=strkin(3)+vxx(i)*vzz(i)*weight(i)
c          strkin(5)=strkin(5)+vyy(i)*vyy(i)*weight(i)
c          strkin(6)=strkin(6)+vyy(i)*vzz(i)*weight(i)
c          strkin(9)=strkin(9)+vzz(i)*vzz(i)*weight(i)
          
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
      
      end subroutine pacmd
      
      subroutine trpmd(lmsite,isw,idnode,mxnode,natms,imcon,
     x  ntpmls,tstep,g_qt4f,temp,engke)
      
c**********************************************************************
c     
c     dl_poly_quantum subroutine for integrating the equations of 
c     motion for ring polymer molecular dynamics - using normal mode 
c     variables 
c     
c     reference: tuckerman, berne, martyna, klein
c     j. chem. phys. vol. 99 (4) p. 2796
c     
c     copyright - 
c     authors   - Dil Limbu and Nathan London 2023
c     
c**********************************************************************
      
      implicit none

      logical lmsite,noskip
      integer, intent(in) :: imcon,ntpmls
      real(8), intent(in) :: g_qt4f
      integer, intent(in)  :: isw,idnode,mxnode,natms
      real(8), intent(in)  :: tstep,temp
      real(8), intent(out) :: engke
      
      integer i,j,k,iatm0,iatm1
      real(8) hstep,strkin(9),uuu(102),beta

      real(8) :: mass(nbeads),rmass(nbeads)

c     block indices
      
      iatm0=nbeads*((idnode*natms)/mxnode)
      iatm1=nbeads*(((idnode+1)*natms)/mxnode)
      
c     thermostat parameters
      
      hstep=tstep/2.d0
      beta=1.d0/(temp*boltz)
      noskip=.true.
      
c     verlet first part
      
      if(isw.eq.1)then
       
c      write(6,*) "vv x pos"
c      write(6,*) uxx(:) 
c      write(6,*) "vv y pos"
c      write(6,*) uyy(:) 
c      write(6,*) "vv z pos"
c      write(6,*) uzz(:) 
c      write(6,*) "vv mom"
c      write(6,*) pxx(:) 
        

c     thermostat momentum - 1/2 step

        do i=1,(iatm1-iatm0)/nbeads
          do j=1,nbeads
c            if(idnode.eq.0)then
c              write(6,*)"index",(i-1)*nbeads+j
c          write(6,*)"j",j
c          write(6,*)"mass",zmass((i-1)*nbeads+j)
c              write(6,*)"cvec1",cvec1(j)
c              write(6,*)"cvec2",cvec2(j)
c            endif
c            write(6,*)"factor",sqrt(zmass((i-1)*nbeads+j)/beta)
c            write(6,*)"facC2",sqrt(zmass((i-1)*nbeads+j)/beta)*cvec2(j)
            pxx((i-1)*nbeads+j)=pxx((i-1)*nbeads+j)*cvec1(j)+
     x        sqrt(zmass((i-1)*nbeads+j)/beta)*cvec2(j)*
     x        gssrnd(noskip,uuu)
            pyy((i-1)*nbeads+j)=pyy((i-1)*nbeads+j)*cvec1(j)+
     x        sqrt(zmass((i-1)*nbeads+j)/beta)*cvec2(j)*
     x        gssrnd(noskip,uuu)
             pzz((i-1)*nbeads+j)=pzz((i-1)*nbeads+j)*cvec1(j)+
     x        sqrt(zmass((i-1)*nbeads+j)/beta)*cvec2(j)*
     x        gssrnd(noskip,uuu)
          enddo
        enddo
c      write(6,*) "vv mom after thermo"
c      write(6,*) pxx(:) 
c      write(6,*) "vv force"
c      write(6,*) wxx(:) 

c     integrate bead momenta - 1/2 step
        
        do i=1,iatm1-iatm0
          
          pxx(i)=pxx(i)+hstep*wxx(i)
          pyy(i)=pyy(i)+hstep*wyy(i)
          pzz(i)=pzz(i)+hstep*wzz(i)
          
        enddo
c      write(6,*) "vv mom after force"
c      write(6,*) pxx(:) 
        
c      write(6,*) "vv mom after thermo 2"
c      write(6,*) pxx(:) 
c     integrate free ring polymer positions (full step) 
        
        do i=1,(iatm1-iatm0)/nbeads
c          write(6,*) "free RP atom",i
          mass(1:nbeads)=zmass((i-1)*nbeads+1:i*nbeads)
          rmass(1:nbeads)=rzmass((i-1)*nbeads+1:i*nbeads)

          call freerp (pxx((i-1)*nbeads+1:i*nbeads),
     x            uxx((i-1)*nbeads+1:i*nbeads),mass,rmass,tstep,temp)
          call freerp (pyy((i-1)*nbeads+1:i*nbeads),
     x            uyy((i-1)*nbeads+1:i*nbeads),mass,rmass,tstep,temp)
          call freerp (pzz((i-1)*nbeads+1:i*nbeads),
     x            uzz((i-1)*nbeads+1:i*nbeads),mass,rmass,tstep,temp)

c          uxx(i)=uxx(i)+(tstep*rzmass(i))*pxx(i)
c          uyy(i)=uyy(i)+(tstep*rzmass(i))*pyy(i)
c          uzz(i)=uzz(i)+(tstep*rzmass(i))*pzz(i)
          
        enddo
c      write(6,*) "vv x pos after free"
c      write(6,*) uxx(:) 
c      write(6,*) "vv y pos after free"
c      write(6,*) uyy(:) 
c      write(6,*) "vv z pos after free"
c      write(6,*) uzz(:) 
c      write(6,*) "vv mom after free"
c      write(6,*) pxx(:) 
        
c     unstage coordinates

      call norm2coord(lmsite,idnode,mxnode,natms,imcon,nbeads,
     x       ntpmls,g_qt4f)

c        call unstage_coords(lmsite,idnode,mxnode,natms,imcon,nbeads,
c     x       ntpmls,g_qt4f)
        
c     verlet second part
        
      else
        
c     integrate bead momenta - 1/2 step
        
        do i=1,iatm1-iatm0
          
          pxx(i)=pxx(i)+hstep*wxx(i)
          pyy(i)=pyy(i)+hstep*wyy(i)
          pzz(i)=pzz(i)+hstep*wzz(i)
          
        enddo
        
c     thermostat momentum - 1/2 step

        do i=1,(iatm1-iatm0)/nbeads
          do j=1,nbeads
            pxx((i-1)*nbeads+j)=pxx((i-1)*nbeads+j)*cvec1(j)+
     x        sqrt(zmass((i-1)*nbeads+j)/beta)*cvec2(j)*
     x        gssrnd(noskip,uuu)
            pyy((i-1)*nbeads+j)=pyy((i-1)*nbeads+j)*cvec1(j)+
     x        sqrt(zmass((i-1)*nbeads+j)/beta)*cvec2(j)*
     x        gssrnd(noskip,uuu)
            pzz((i-1)*nbeads+j)=pzz((i-1)*nbeads+j)*cvec1(j)+
     x        sqrt(zmass((i-1)*nbeads+j)/beta)*cvec2(j)*
     x        gssrnd(noskip,uuu)
          enddo
        enddo

c     unstage momenta (needed for REVCON file)
        
        call norm2momenta(idnode,mxnode,natms)
        
c        call unstage_momenta(idnode,mxnode,natms)
        
c     calculate  kinetic tensor and kinetic energy
        
        strkin(:)=0.d0
        
        do i=1,iatm1-iatm0
          
          strkin(1)=strkin(1)+pxx(i)*pxx(i)*rzmass(i)
          strkin(2)=strkin(2)+pxx(i)*pyy(i)*rzmass(i)
          strkin(3)=strkin(3)+pxx(i)*pzz(i)*rzmass(i)
          strkin(5)=strkin(5)+pyy(i)*pyy(i)*rzmass(i)
          strkin(6)=strkin(6)+pyy(i)*pzz(i)*rzmass(i)
          strkin(9)=strkin(9)+pzz(i)*pzz(i)*rzmass(i)
          
c          strkin(1)=strkin(1)+vxx(i)*vxx(i)*weight(i)
c          strkin(2)=strkin(2)+vxx(i)*vyy(i)*weight(i)
c          strkin(3)=strkin(3)+vxx(i)*vzz(i)*weight(i)
c          strkin(5)=strkin(5)+vyy(i)*vyy(i)*weight(i)
c          strkin(6)=strkin(6)+vyy(i)*vzz(i)*weight(i)
c          strkin(9)=strkin(9)+vzz(i)*vzz(i)*weight(i)
          
        enddo
        
        strkin(4)=strkin(2)
        strkin(7)=strkin(3)
        strkin(8)=strkin(6)
        
        if(mxnode.gt.1)then
          
          buffer(1:9)=strkin(1:9)
          call gdsum(buffer(1),10,buffer(11))
          strkin(1:9)=buffer(1:9)
          
        endif
        
        engke=0.5d0*(strkin(1)+strkin(5)+strkin(9))
        stress(:)=stress(:)+strkin(:)
        
      endif
      
      end subroutine trpmd

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
            
            do j=nchain-1,2,-1
              
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
      
c      qqq(1)=boltz*temp*taut**2
c      do k=2,nbeads
c        qqq(k)=hbar**2/(boltz*temp*dble(nbeads))
c      enddo
      
      qqq(1)=4.d0*nbeads*boltz*temp*taut**2
      do k=2,nbeads
        qqq(k)=(boltz*temp)/nmfreq(k)**2
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
                
                do j=nchain-1,2,-1
                  
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


      subroutine pimd_nvt_pile_nm(lmsite,isw,idnode,mxnode,natms,imcon,
     x  ntpmls,tstep,taut,g_qt4f,temp,engke,engthe,uuu)
      
c**********************************************************************
c     
c     dl_poly_quantum subroutine for integrating the equations of 
c     motion for path intergral molecular dynamics - using normal mode 
c     variables thermostated with PILE  and integrated 
c     with the velocity verlet algorithm
c     
c     reference: Ceriotti,Parrinello, Markland and Manolopoulos
c     j. chem. phys. 133, 124104 (2010)
c     
c     copyright - 
c     authors   - Dil Limbu and Nathan London 2023
c     
c**********************************************************************
      
      implicit none

      logical :: lmsite,noskip
      integer, intent(in)  :: imcon,ntpmls
      real(8), intent(in)  :: g_qt4f
      integer, intent(in)  :: isw,idnode,mxnode,natms
      real(8), intent(in)  :: tstep,taut,temp,uuu(102)
      real(8), intent(out) :: engke,engthe
      
      integer :: init,nnn
      integer :: i,j,k,iatm0,iatm1
      real(8) :: pC2,qq1,qqk,hstep,strkin(9)
      real(8) :: mass(nbeads),rmass(nbeads),gama(nbeads)

      data init /0/
      save init
c     block indices
      
      iatm0=nbeads*((idnode*natms)/mxnode)
      iatm1=nbeads*(((idnode+1)*natms)/mxnode)
      
c     mass parameters of thermostats
      
      qq1=boltz*temp*taut**2
      qqk=hbar**2/(boltz*temp*dble(nbeads))
      
c     thermostat parameters
c     initialize C1 and C2 vector for nbeads

      if (init.eq.0) then

        hstep=tstep/2.d0

        do k=1,nbeads

          gama(k)=2.d0*nmfreq(k)
          if((nbeads.gt.1).and.(k.eq.1).and.(taut.ne.1.d0))
     x           gama(1)=1.d0/taut
          pileC1(k)=dexp(-hstep*gama(k))
          pileC2(k)=dsqrt(1.d0-pileC1(k)*pileC1(k))

c          write(6,*) 'c1,c2:', pileC1(k),pileC2(k)

        enddo

        init=1
      endif

c     verlet first part
      
      if(isw.eq.1)then
        
c     apply PILE thermostat to bead momenta - 1/2 step
        
        noskip=.true.

        do i=1,iatm1-iatm0,nbeads

          do k=1,nbeads

            j=(i-1)+k
            pC2=sqrt(zmass(j)*boltz*temp)*pileC2(k)

            pxx(j)=pileC1(k)*pxx(j)+pC2*gssrnd(noskip,uuu)
            pyy(j)=pileC1(k)*pyy(j)+pC2*gssrnd(noskip,uuu)
            pzz(j)=pileC1(k)*pzz(j)+pC2*gssrnd(noskip,uuu)

          enddo

        enddo
        
c     integrate bead momenta - 1/2 step
        
        do i=1,iatm1-iatm0
          
          pxx(i)=pxx(i)+hstep*wxx(i)
          pyy(i)=pyy(i)+hstep*wyy(i)
          pzz(i)=pzz(i)+hstep*wzz(i)
          
        enddo
        
c     integrate free ring polymer positions (full step) 
        
        do i=1,iatm1-iatm0,nbeads
          
          j=i-1+nbeads

          mass(1:nbeads)=zmass(i:j)
          rmass(1:nbeads)=rzmass(i:j)

          call freerp(pxx(i:j),uxx(i:j),mass,rmass,tstep,temp)
          call freerp(pyy(i:j),uyy(i:j),mass,rmass,tstep,temp)
          call freerp(pzz(i:j),uzz(i:j),mass,rmass,tstep,temp)
          
        enddo
        
c     unstage coordinates

        call norm2coord(lmsite,idnode,mxnode,natms,imcon,nbeads,
     x       ntpmls,g_qt4f)
        
c     verlet second part
        
      else
        
c     integrate bead momenta - 1/2 step
        
        do i=1,iatm1-iatm0
          
          pxx(i)=pxx(i)+hstep*wxx(i)
          pyy(i)=pyy(i)+hstep*wyy(i)
          pzz(i)=pzz(i)+hstep*wzz(i)
          
        enddo
        
c     apply PILE thermostat to bead momenta - 1/2 step
        
        noskip=.true.

        do i=1,iatm1-iatm0,nbeads
          
          do k=1,nbeads

            j=(i-1)+k
            pC2=sqrt(zmass(j)*boltz*temp)*pileC2(k)

            pxx(j)=pileC1(k)*pxx(j)+pC2*gssrnd(noskip,uuu)
            pyy(j)=pileC1(k)*pyy(j)+pC2*gssrnd(noskip,uuu)
            pzz(j)=pileC1(k)*pzz(j)+pC2*gssrnd(noskip,uuu)

          enddo

        enddo
        
c     calculate thermostat energy
        
c        call thermostat_energy
c     x    (idnode,mxnode,natms,nbeads,nchain,qq1,qqk,temp,engthe)
        
c     unstage momenta (needed for REVCON file)
        
        call norm2momenta(idnode,mxnode,natms)
        
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
      
      end subroutine pimd_nvt_pile_nm
      
      
      subroutine pimd_nvt_piglet(lmsite,isw,idnode,mxnode,natms,imcon,
     x  ntpmls,tstep,taut,g_qt4f,temp,engke,engthe,uuu,
     x  virtot,vircon,virrng)
      
c**********************************************************************
c     
c     dl_poly_quantum subroutine for integrating the equations of 
c     motion for path intergral molecular dynamics - using normal mode 
c     variables thermostated with PIGLET THERMOSTAT and integrated 
c     with the velocity verlet algorithm
c     
c     reference:
c     
c     copyright - 
c     authors   - Dil Limbu and Nathan London 2023
c     
c**********************************************************************

      implicit none

      logical              :: lmsite
      integer, intent(in)  :: imcon,ntpmls
      real(8), intent(in)  :: g_qt4f
      integer, intent(in)  :: isw,idnode,mxnode,natms
      real(8), intent(in)  :: tstep,taut,temp,uuu(102)
      real(8), intent(out) :: engke,engthe

      integer :: i,j,k,iatm0,iatm1
      real(8) :: pC2,qq1,qqk,hstep,strkin(9)
      real(8) :: mass(nbeads),rmass(nbeads),gama(nbeads)

      real(8) :: virtot,vircon,virrng
      real(8) :: engkec,vir2
c     block indices
      
      iatm0=nbeads*((idnode*natms)/mxnode)
      iatm1=nbeads*(((idnode+1)*natms)/mxnode)
      
      hstep=0.5d0*tstep

c     mass parameters of thermostats
      
c      qq1=boltz*temp*taut**2
c      qqk=hbar**2/(boltz*temp*dble(nbeads))
      
c     verlet first part
      
      if(isw.eq.1)then
        
c     apply piglet thermostat to bead momenta - 1/2 step
        
        call piglet_thermo_step(idnode,mxnode,natms,tstep,temp,uuu)
        
c     integrate bead momenta - 1/2 step
        
        do i=1,iatm1-iatm0
          
          pxx(i)=pxx(i)+hstep*wxx(i)
          pyy(i)=pyy(i)+hstep*wyy(i)
          pzz(i)=pzz(i)+hstep*wzz(i)
          
        enddo
        
c     integrate free ring polymer positions (full step) 
        
        do i=1,iatm1-iatm0,nbeads
          
          j=i-1+nbeads

          mass(1:nbeads)=zmass(i:j)
          rmass(1:nbeads)=rzmass(i:j)

          call freerp(pxx(i:j),uxx(i:j),mass,rmass,tstep,temp)
          call freerp(pyy(i:j),uyy(i:j),mass,rmass,tstep,temp)
          call freerp(pzz(i:j),uzz(i:j),mass,rmass,tstep,temp)
          
        enddo
        
c     unstage coordinates

        call norm2coord(lmsite,idnode,mxnode,natms,imcon,nbeads,
     x       ntpmls,g_qt4f)
        
c     verlet second part
        
      else
        
c     integrate bead momenta - 1/2 step
        
        do i=1,iatm1-iatm0
          
          pxx(i)=pxx(i)+hstep*wxx(i)
          pyy(i)=pyy(i)+hstep*wyy(i)
          pzz(i)=pzz(i)+hstep*wzz(i)
          
        enddo
        
c     apply piglet thermostat to bead momenta - 1/2 step
        
        call piglet_thermo_step(idnode,mxnode,natms,tstep,temp,uuu)
        
c     calculate thermostat energy
        
c        call thermostat_energy
c     x    (idnode,mxnode,natms,nbeads,nchain,qq1,qqk,temp,engthe)
        
c     unstage momenta (needed for REVCON file)
        
        call norm2momenta(idnode,mxnode,natms)
        
c     calculate  kinetic tensor and kinetic energy
        
        strkin(:)=0.d0

c     for centroid mode only

        do i=1,iatm1-iatm0,nbeads
          
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
        
        engke=0.5d0*(strkin(1)+strkin(5)+strkin(9))*nbeads
        stress(:)=stress(:)+strkin(:)*nbeads
        
c      endif
c  testing virial as in i-Pi

c        engkec=getkec(natms,idnode,mxnode)
        
        call getvir2(natms,idnode,mxnode,vir2)

c        if(idnode.eq.0)write(2001,'(9e15.6)')engke,virtot,vircon,virrng,
c     x     vir2
        
      endif

c 
      end subroutine pimd_nvt_piglet
      
      
      subroutine pimd_npt_nhc_nm
     x    (safe,lmsite,isw,idnode,mxnode,natms,imcon,
     x    nrespa,ntpmls,        
     x    ntpatm,ntshl,keyshl,tstep,taut,taup,sigma,
     x    sigma_nhc,sigma_volm,alpha_volm,virtot,vircon,virlrc,
     x    g_qt4f,press,volm,chit,consv,
     x    conint,engke,elrc,chit_shl,sigma_shl,temp,engthe)
      
c**********************************************************************
c     
c     dl_poly_quantum subroutine for integrating the equations of 
c     motion for path intergral molecular dynamics - using normal mode 
c     variables thermostated with nose-hoover chains barostat and 
c     integrated with the velocity verlet algorithm
c     
c     reference: tuckerman, berne, martyna, klein
c     j. chem. phys. vol. 99 (4) p. 2796
c     
c     copyright - 
c     authors   - Dil Limbu and Nathan London 2023
c     
c**********************************************************************
      
      implicit none

      logical :: safe,newjob,lmsite

      integer, parameter :: nnn=7
      integer :: fail(nnn)
      integer :: ntshl,keyshl
      integer :: isw,idnode,mxnode,natms,imcon,ntpatm
      integer :: i,j,k,iatm0,iatm1,kk
      integer :: nrespa,ntpmls

      real(8) :: temp,engthe
      real(8) :: qq1,qqk
      real(8) :: qstep
      real(8) :: mass(nbeads),rmass(nbeads)
      real(8) :: chit_shl,sigma_shl
      real(8) :: tstep,taut,taup,chit,consv,conint,engke,elrc
      real(8) :: hstep,sigma,qmass_t
      real(8) :: sigma_nhc,qmass_baro,qmass_part
      real(8) :: volm_mass,sigma_volm,alpha_volm,g_qt4f
      real(8) :: press,volm,virtot,vircon,virlrc
      real(8) :: heta_nhc,hpeta_nhc,heta_1,heta_rest,hpeta_1,hpeta_rest
      real(8) :: hepsilon,hksi,hpksi
      real(8) :: volm0,elrc0,virlrc0,scale
      real(8) :: a_2n(6),sinh_v,sinh_r
      real(8) :: cell0(9),com(3),vom(3),strkin(9),uni(9)
      real(8) :: Pint,part1,part2,part3,engkec,exp_a
      real(8), allocatable :: dens0(:)

      real(8) :: vir2

      save newjob,volm0,elrc0,virlrc0,dens0
      save cell0,iatm0,iatm1,hstep,qstep
      save qmass_t,qmass_baro,qmass_part,volm_mass

      data uni/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/
      data newjob/.true./

c     thermostat parameters
      
c      qq1=boltz*temp*taut**2
c      qqk=hbar**2/(boltz*temp*dble(nbeads))

      if(newjob)then

c     block indices

        iatm0=nbeads*((idnode*natms)/mxnode)
        iatm1=nbeads*(((idnode+1)*natms)/mxnode)

c     allocate density storage array

        fail(1)=0
        allocate(dens0(mxatyp),stat=fail(1))
        if(fail(1).gt.0)call error(idnode,2040)

c     store intitial parameters

        volm0=volm
        elrc0=elrc
        virlrc0=virlrc

        hstep=tstep/2.d0
        qstep=tstep/4.d0

        do i=1,ntpatm
          dens0(i)=dens(i)
        enddo

        do i=1,9
          cell0(i)=cell(i)
        enddo
        
c     inertia parameter for Nose-Hoover thermostat and barostat

        qmass_t=2.d0*sigma*taut**2
        qmass_part=2.d0*sigma_nhc*taut**2
        qmass_baro=2.d0*sigma_nhc*taup**2
        volm_mass=sigma_volm*taup**2

        newjob=.false.

      endif
        
c     a_2n parameters for NHC barostat

        a_2n(1)=1.d0
        a_2n(2)=1.d0/6.d0
        a_2n(3)=1.d0/120.d0
        a_2n(4)=1.d0/5040.d0
        a_2n(5)=1.d0/362880.d0
        a_2n(6)=1.d0/39916800.d0

c      if(lmetadyn.and.lfirst.and.(ntshl>0))then
c        if(idnode.eq.0)then
c          write(*,*)"Warning - Metadynamics Modification"
c          write(*,*)"========================="
c          write(*,*)"Coupling core-shell motion thermostat at 1 K"
c        endif
c        lfirst=.false.
c     use same relaxation time for global and core-shell?
c        qmass_shl=2.d0*sigma_shl*taut**2
c      endif

c     verlet first part
      
      if(isw.eq.1)then

c     integrate and apply nhc barostat - 1/2 step

        call nhc_baro
     x     (idnode,mxnode,natms,nchain,nrespa,sigma_nhc,
     x     hstep,volm_mass,qmass_baro,taup,v_epsilon)
 
c     integrate thermostats - 1/2 step
        
        if(nchain.gt.0)call thermo_chain_sy
     x    (idnode,mxnode,natms,hstep,taut,temp)
        
        engke=getke(natms,idnode,mxnode)
        
c        call getvir2(natms,idnode,mxnode,vir2)

c        Pint=(2.d0*engke/nbeads+vir2-vircon)/(3.d0*volm)

        Pint=(2.d0*engke-virtot-vircon)/(3.d0*volm)

        part1=3.d0*(volm*(Pint-press)+boltz*temp)

        part2=0.d0
        part3=0.d0

        do i=1,iatm1-iatm0,nbeads
          
          part2=part2+rzmass(i)*(wxx(i)*pxx(i)+wyy(i)*pyy(i)+
     x          wzz(i)*pzz(i))

          part3=part3+rzmass(i)/3.d0*(wxx(i)*wxx(i)+wyy(i)*wyy(i)+
     x          wzz(i)*wzz(i))
          
        enddo
        
        if(mxnode.gt.1)then
          
          buffer(1)=part2
          buffer(2)=part3
          call gdsum(buffer(1),2,buffer(3))
          part2=buffer(1)
          part3=buffer(2)
          
        endif
        
c     update the momentum variable of the volume

        v_epsilon=v_epsilon+(1.d0/volm_mass)*
     x     (hstep*part1+(hstep**2)*part2+(hstep**3)*part3)

c        v_epsilon=v_epsilon+hstep*(1.d0/volm_mass)*
c     x     ((alpha_volm*2.0*engke)-virtot-vircon-(press*volm)) 

c     integrate momenta - 1/2 step
        
        do i=1,iatm1-iatm0
          
          pxx(i)=pxx(i)+hstep*wxx(i)
          pyy(i)=pyy(i)+hstep*wyy(i)
          pzz(i)=pzz(i)+hstep*wzz(i)
          
        enddo

c     update the volume based on v_epsilon

        volm=volm*exp(3.d0*tstep*v_epsilon)

        if(mxnode.gt.1)then
          
          buffer(1)=volm
          call gdsum(buffer(1),1,buffer(2))
          volm=buffer(1)
          volm=volm/(mxnode)
          
        endif

c     scale cell vectors - isotropic

          scale=(volm/volm0)**(1.d0/3.d0)

          do i=1,9
            cell(i)=cell0(i)*scale
          enddo

c     first integrate sinh(alpha*v_epsilon*tstep/4)

c        sinh_v=0.d0
c        do i=1,6
c          sinh_v=sinh_v+a_2n(i)
c     x          *(alpha_volm*v_epsilon*hstep*0.5d0)**(2*(i-1))
c        enddo

c     first integrate sinh(v_epsilon*tstep)

        sinh_r=0.d0 
        do i=1,6
          sinh_r=sinh_r+a_2n(i)*(v_epsilon*tstep)**(2*(i-1))
        enddo   

c     integrate free ring polymer centroid mode (full step) 

        exp_a=exp(-v_epsilon*tstep)

        do i=1,iatm1-iatm0,nbeads

          pxx(i)=pxx(i)*exp_a
          pyy(i)=pyy(i)*exp_a
          pzz(i)=pzz(i)*exp_a

          uxx(i)=uxx(i)/exp_a+tstep*pxx(i)*rzmass(i)*sinh_r
          uyy(i)=uyy(i)/exp_a+tstep*pyy(i)*rzmass(i)*sinh_r
          uzz(i)=uzz(i)/exp_a+tstep*pzz(i)*rzmass(i)*sinh_r

        enddo

c     integrate free ring polymer (non-centroid) positions (full step) 
        
        do i=1,iatm1-iatm0,nbeads
          
          j=i-1+nbeads

          mass(1:nbeads)=zmass(i:j)
          rmass(1:nbeads)=rzmass(i:j)

          call freerp_noc(pxx(i:j),uxx(i:j),mass,rmass,tstep,temp)
          call freerp_noc(pyy(i:j),uyy(i:j),mass,rmass,tstep,temp)
          call freerp_noc(pzz(i:j),uzz(i:j),mass,rmass,tstep,temp)
          
        enddo
       
c     unstage coordinates

        call norm2coord(lmsite,idnode,mxnode,natms,imcon,nbeads,
     x       ntpmls,g_qt4f)

c     verlet second part
        
      else

c     integrate momenta - 1/2 step
        
        do i=1,iatm1-iatm0
          
          pxx(i)=pxx(i)+hstep*wxx(i)
          pyy(i)=pyy(i)+hstep*wyy(i)
          pzz(i)=pzz(i)+hstep*wzz(i)
          
        enddo
        
c     update the momentum variable of the volume

        engke=getke(natms,idnode,mxnode)

c        call getvir2(natms,idnode,mxnode,vir2)

c        Pint=(2.d0*engke/nbeads+vir2-vircon)/(3.d0*volm)

        Pint=(2.d0*engke-virtot-vircon)/(3.d0*volm)

        part1=3.d0*(volm*(Pint-press)+boltz*temp)

        part2=0.d0
        part3=0.d0

        do i=1,iatm1-iatm0,nbeads
          
          part2=part2+rzmass(i)*(wxx(i)*pxx(i)+wyy(i)*pyy(i)+
     x          wzz(i)*pzz(i))

          part3=part3+rzmass(i)/3.d0*(wxx(i)*wxx(i)+wyy(i)*wyy(i)+
     x          wzz(i)*wzz(i))
          
        enddo

        if(mxnode.gt.1)then
          
          buffer(1)=part2
          buffer(2)=part3
          call gdsum(buffer(1),2,buffer(3))
          part2=buffer(1)
          part3=buffer(2)
          
        endif
        
c     update the momentum variable of the volume

        v_epsilon=v_epsilon+(1.d0/volm_mass)*
     x     (hstep*part1+(hstep**2)*part2+(hstep**3)*part3)

c        v_epsilon=v_epsilon+hstep*(1.d0/volm_mass)*
c     x     ((alpha_volm*2.0*engke)-virtot-vircon-(press*volm)) 

c     integrate thermostats - 1/2 step
        
        if(nchain.gt.0)call thermo_chain_sy
     x    (idnode,mxnode,natms,hstep,taut,temp)

c     integrate and apply nhc barostat - 1/2 step

        call nhc_baro
     x     (idnode,mxnode,natms,nchain,nrespa,sigma_nhc,
     x     hstep,volm_mass,qmass_baro,taup,v_epsilon)

c     conserved quantity less kinetic and potential energy terms

       hepsilon=0.5d0*volm_mass*(v_epsilon**2)

       hksi=0.d0
         do j=1,nchain
            hksi=hksi+ksi(j)
         enddo

       hksi=2.d0*sigma_nhc*hksi

       hpksi=0.d0
         do j=1,nchain
            hpksi=hpksi+0.5d0*(pksi(j)**2/qmass_baro)
         enddo

       heta_1=2.d0*sigma*eta_nhc(1)

       heta_rest=0.d0
         do j=2,nchain
            heta_rest=heta_rest+eta_nhc(j)
         enddo

       heta_nhc=(2.d0*sigma_nhc*heta_rest)+heta_1

       hpeta_1=0.5d0*(peta(1)**2/qmass_t)

       hpeta_rest=0.d0
         do j=2,nchain
            hpeta_rest=hpeta_rest+0.5d0*(peta(j)**2/qmass_part)
         enddo

       hpeta_nhc=hpeta_1+hpeta_rest

       consv=hepsilon+hksi+hpksi+heta_nhc+hpeta_nhc+(press*volm)

c     kinetic contribution to stress tensor
        
c        call kinstress(natms,idnode,mxnode,strkin)
c        engke=0.5d0*(strkin(1)+strkin(5)+strkin(9))

c     metadynamics shell thermostat

c        if(lmetadyn.and.keyshl.eq.1)then
c           consv=consv+0.5d0*qmass_shl*chit_shl**2
c        endif

c     add contributions to stress tensor
        
c        do i=1,9
c          stress(i)=stress(i)+strkin(i)+strcns(i)
c        enddo

c     calculate thermostat energy
        
c        call thermostat_energy
c     x    (idnode,mxnode,natms,nbeads,nchain,qq1,qqk,temp,engthe)
        
c     unstage momenta (needed for REVCON file)
        
        call norm2momenta(idnode,mxnode,natms)
        
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
        
c        engke=0.5d0*(strkin(1)+strkin(5)+strkin(9))
        stress(:)=stress(:)+strkin(:)
 
       press2=Pint       
c        if(idnode.eq.0)write(2001,'(6e15.6)')engke,virtot,vircon,
c     x      vir2,Pint,(2.d0*engke+vir2)/(3.d0*volm)

      endif

c     adjust long range corrections and number density

      elrc=elrc0*(volm0/volm)
      virlrc=virlrc0*(volm0/volm)
      do kk=1,ntpatm
        dens(kk)=dens0(kk)*(volm0/volm)
      enddo

      return
      
      end subroutine pimd_npt_nhc_nm
     
 
      subroutine pimd_npt_piglet_nm
     x    (safe,lmsite,isw,idnode,mxnode,natms,imcon,
     x    nrespa,ntpmls,        
     x    ntpatm,ntshl,keyshl,tstep,taut,taup,sigma,
     x    sigma_nhc,sigma_volm,alpha_volm,virtot,vircon,virlrc,
     x    g_qt4f,press,volm,chit,consv,
     x    conint,engke,elrc,chit_shl,sigma_shl,temp,engthe,uuu)
      
c**********************************************************************
c     
c     dl_poly_quantum subroutine for integrating the equations of 
c     motion for path intergral molecular dynamics - using normal mode 
c     variables thermostated with 
c
c     nose-hoover chains barostat and 
c     piglet thermostat 
c     integrated with the velocity verlet algorithm
c     
c     reference: tuckerman, berne, martyna, klein
c     j. chem. phys. vol. 99 (4) p. 2796
c     
c     copyright - 
c     authors   - Dil Limbu and Nathan London 2023
c     
c**********************************************************************
      
      implicit none

      logical :: safe,newjob,lmsite

      integer, parameter :: nnn=7
      integer :: fail(nnn)
      integer :: ntshl,keyshl
      integer :: isw,idnode,mxnode,natms,imcon,ntpatm
      integer :: i,j,k,iatm0,iatm1,kk
      integer :: nrespa,ntpmls

      real(8) :: temp,engthe
      real(8) :: qq1,qqk
      real(8) :: qstep
      real(8) :: mass(nbeads),rmass(nbeads)
      real(8) :: chit_shl,sigma_shl
      real(8) :: tstep,taut,taup,chit,consv,conint,engke,elrc
      real(8) :: hstep,sigma,qmass_t
      real(8) :: sigma_nhc,qmass_baro,qmass_part
      real(8) :: volm_mass,sigma_volm,alpha_volm,g_qt4f
      real(8) :: press,volm,virtot,vircon,virlrc
      real(8) :: heta_nhc,hpeta_nhc,heta_1,heta_rest,hpeta_1,hpeta_rest
      real(8) :: hepsilon,hksi,hpksi
      real(8) :: volm0,elrc0,virlrc0,scale
      real(8) :: a_2n(6),sinh_v,sinh_r
      real(8) :: cell0(9),com(3),vom(3),strkin(9),uni(9)
      real(8) :: Pint,part1,part2,part3,engkec,exp_a
      real(8), allocatable :: dens0(:)

      real(8) :: vir2
      real(8) :: uuu(102)

      save newjob,volm0,elrc0,virlrc0,dens0
      save cell0,iatm0,iatm1,hstep,qstep
      save qmass_t,qmass_baro,qmass_part,volm_mass

      data uni/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/
      data newjob/.true./

c     thermostat parameters
      
c      qq1=boltz*temp*taut**2
c      qqk=hbar**2/(boltz*temp*dble(nbeads))

      if(newjob)then

c     block indices

        iatm0=nbeads*((idnode*natms)/mxnode)
        iatm1=nbeads*(((idnode+1)*natms)/mxnode)

c     allocate density storage array

        fail(1)=0
        allocate(dens0(mxatyp),stat=fail(1))
        if(fail(1).gt.0)call error(idnode,2040)

c     store intitial parameters

        volm0=volm
        elrc0=elrc
        virlrc0=virlrc

        hstep=tstep/2.d0
        qstep=tstep/4.d0

        do i=1,ntpatm
          dens0(i)=dens(i)
        enddo

        do i=1,9
          cell0(i)=cell(i)
        enddo
        
c     inertia parameter for Nose-Hoover thermostat and barostat

        qmass_t=2.d0*sigma*taut**2
        qmass_part=2.d0*sigma_nhc*taut**2
        qmass_baro=2.d0*sigma_nhc*taup**2
        volm_mass=sigma_volm*taup**2

        newjob=.false.

      endif
        
c     a_2n parameters for NHC barostat

        a_2n(1)=1.d0
        a_2n(2)=1.d0/6.d0
        a_2n(3)=1.d0/120.d0
        a_2n(4)=1.d0/5040.d0
        a_2n(5)=1.d0/362880.d0
        a_2n(6)=1.d0/39916800.d0

c      if(lmetadyn.and.lfirst.and.(ntshl>0))then
c        if(idnode.eq.0)then
c          write(*,*)"Warning - Metadynamics Modification"
c          write(*,*)"========================="
c          write(*,*)"Coupling core-shell motion thermostat at 1 K"
c        endif
c        lfirst=.false.
c     use same relaxation time for global and core-shell?
c        qmass_shl=2.d0*sigma_shl*taut**2
c      endif

c     verlet first part
      
      if(isw.eq.1)then

        call nhc_baro
     x     (idnode,mxnode,natms,nchain,nrespa,sigma_nhc,
     x     hstep,volm_mass,qmass_baro,taup,v_epsilon)
 
c     apply piglet thermostat to bead momenta - 1/2 step
        
        call piglet_thermo_step(idnode,mxnode,natms,tstep,temp,uuu)
        
        engke=getke(natms,idnode,mxnode)
        
c        call getvir2(natms,idnode,mxnode,vir2)

c        Pint=(2.d0*engke/nbeads+vir2-vircon)/(3.d0*volm)

        Pint=(2.d0*engke-virtot-vircon)/(3.d0*volm)

        part1=3.d0*(volm*(Pint-press)+boltz*temp)

        part2=0.d0
        part3=0.d0

        do i=1,iatm1-iatm0,nbeads
          
          part2=part2+rzmass(i)*(wxx(i)*pxx(i)+wyy(i)*pyy(i)+
     x          wzz(i)*pzz(i))

          part3=part3+rzmass(i)/3.d0*(wxx(i)*wxx(i)+wyy(i)*wyy(i)+
     x          wzz(i)*wzz(i))
          
        enddo
        
        if(mxnode.gt.1)then
          
          buffer(1)=part2
          buffer(2)=part3
          call gdsum(buffer(1),2,buffer(3))
          part2=buffer(1)
          part3=buffer(2)
          
        endif
        
c     update the momentum variable of the volume

        v_epsilon=v_epsilon+(1.d0/volm_mass)*
     x     (hstep*part1+(hstep**2)*part2+(hstep**3)*part3)

c        v_epsilon=v_epsilon+hstep*(1.d0/volm_mass)*
c     x     ((alpha_volm*2.0*engke)-virtot-vircon-(press*volm)) 

c     integrate bead momenta - 1/2 step
        
        do i=1,iatm1-iatm0
          
          pxx(i)=pxx(i)+hstep*wxx(i)
          pyy(i)=pyy(i)+hstep*wyy(i)
          pzz(i)=pzz(i)+hstep*wzz(i)
          
        enddo
        
c     update the volume based on v_epsilon

        volm=volm*exp(3.d0*tstep*v_epsilon)

        if(mxnode.gt.1)then
          
          buffer(1)=volm
          call gdsum(buffer(1),1,buffer(2))
          volm=buffer(1)
          volm=volm/(mxnode)
          
        endif

c     scale cell vectors - isotropic

          scale=(volm/volm0)**(1.d0/3.d0)

          do i=1,9
            cell(i)=cell0(i)*scale
          enddo

c     first integrate sinh(alpha*v_epsilon*tstep/4)

c        sinh_v=0.d0
c        do i=1,6
c          sinh_v=sinh_v+a_2n(i)
c     x          *(alpha_volm*v_epsilon*hstep*0.5d0)**(2*(i-1))
c        enddo

c     first integrate sinh(v_epsilon*tstep)

        sinh_r=0.d0 
        do i=1,6
          sinh_r=sinh_r+a_2n(i)*(v_epsilon*tstep)**(2*(i-1))
        enddo   

c     integrate free ring polymer centroid mode (full step) 

        exp_a=exp(-v_epsilon*tstep)

        do i=1,iatm1-iatm0,nbeads

          pxx(i)=pxx(i)*exp_a
          pyy(i)=pyy(i)*exp_a
          pzz(i)=pzz(i)*exp_a

          uxx(i)=uxx(i)/exp_a+tstep*pxx(i)*rzmass(i)*sinh_r
          uyy(i)=uyy(i)/exp_a+tstep*pyy(i)*rzmass(i)*sinh_r
          uzz(i)=uzz(i)/exp_a+tstep*pzz(i)*rzmass(i)*sinh_r

        enddo

c     integrate free ring polymer (non-centroid) positions (full step) 
        
        do i=1,iatm1-iatm0,nbeads
          
          j=i-1+nbeads

          mass(1:nbeads)=zmass(i:j)
          rmass(1:nbeads)=rzmass(i:j)

          call freerp_noc(pxx(i:j),uxx(i:j),mass,rmass,tstep,temp)
          call freerp_noc(pyy(i:j),uyy(i:j),mass,rmass,tstep,temp)
          call freerp_noc(pzz(i:j),uzz(i:j),mass,rmass,tstep,temp)
          
        enddo
       
c     unstage coordinates

        call norm2coord(lmsite,idnode,mxnode,natms,imcon,nbeads,
     x       ntpmls,g_qt4f)

c     verlet second part
        
      else

c     integrate bead momenta - 1/2 step
        
        do i=1,iatm1-iatm0
          
          pxx(i)=pxx(i)+hstep*wxx(i)
          pyy(i)=pyy(i)+hstep*wyy(i)
          pzz(i)=pzz(i)+hstep*wzz(i)
          
        enddo
        
c     update the momentum variable of the volume

        engke=getke(natms,idnode,mxnode)

c        call getvir2(natms,idnode,mxnode,vir2)

c        Pint=(2.d0*engke/nbeads+vir2-vircon)/(3.d0*volm)

        Pint=(2.d0*engke-virtot-vircon)/(3.d0*volm)

        part1=3.d0*(volm*(Pint-press)+boltz*temp)

        part2=0.d0
        part3=0.d0

        do i=1,iatm1-iatm0,nbeads
          
          part2=part2+rzmass(i)*(wxx(i)*pxx(i)+wyy(i)*pyy(i)+
     x          wzz(i)*pzz(i))

          part3=part3+rzmass(i)/3.d0*(wxx(i)*wxx(i)+wyy(i)*wyy(i)+
     x          wzz(i)*wzz(i))
          
        enddo

        if(mxnode.gt.1)then
          
          buffer(1)=part2
          buffer(2)=part3
          call gdsum(buffer(1),2,buffer(3))
          part2=buffer(1)
          part3=buffer(2)
          
        endif
        

c     update the momentum variable of the volume

        v_epsilon=v_epsilon+(1.d0/volm_mass)*
     x     (hstep*part1+(hstep**2)*part2+(hstep**3)*part3)

c        v_epsilon=v_epsilon+hstep*(1.d0/volm_mass)*
c     x     ((alpha_volm*2.0*engke)-virtot-vircon-(press*volm)) 

c     apply piglet thermostat to bead momenta - 1/2 step
        
        call piglet_thermo_step(idnode,mxnode,natms,tstep,temp,uuu)
        
c     integrate and apply nhc barostat

        call nhc_baro
     x     (idnode,mxnode,natms,nchain,nrespa,sigma_nhc,
     x     hstep,volm_mass,qmass_baro,taup,v_epsilon)

c     conserved quantity less kinetic and potential energy terms

       hepsilon=0.5d0*volm_mass*(v_epsilon**2)

       hksi=0.d0
         do j=1,nchain
            hksi=hksi+ksi(j)
         enddo

       hksi=2.d0*sigma_nhc*hksi

       hpksi=0.d0
         do j=1,nchain
            hpksi=hpksi+0.5d0*(pksi(j)**2/qmass_baro)
         enddo

       heta_1=2.d0*sigma*eta_nhc(1)

       heta_rest=0.d0
         do j=2,nchain
            heta_rest=heta_rest+eta_nhc(j)
         enddo

       heta_nhc=(2.d0*sigma_nhc*heta_rest)+heta_1

       hpeta_1=0.5d0*(peta(1)**2/qmass_t)

       hpeta_rest=0.d0
         do j=2,nchain
            hpeta_rest=hpeta_rest+0.5d0*(peta(j)**2/qmass_part)
         enddo

       hpeta_nhc=hpeta_1+hpeta_rest

       consv=hepsilon+hksi+hpksi+heta_nhc+hpeta_nhc+(press*volm)

c     kinetic contribution to stress tensor
        
c        call kinstress(natms,idnode,mxnode,strkin)
c        engke=0.5d0*(strkin(1)+strkin(5)+strkin(9))

c     metadynamics shell thermostat

c        if(lmetadyn.and.keyshl.eq.1)then
c           consv=consv+0.5d0*qmass_shl*chit_shl**2
c        endif

c     add contributions to stress tensor
        
c        do i=1,9
c          stress(i)=stress(i)+strkin(i)+strcns(i)
c        enddo

c     calculate thermostat energy
        
c        call thermostat_energy
c     x    (idnode,mxnode,natms,nbeads,nchain,qq1,qqk,temp,engthe)
        
c     unstage momenta (needed for REVCON file)
        
        call norm2momenta(idnode,mxnode,natms)
        
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
        
c        engke=0.5d0*(strkin(1)+strkin(5)+strkin(9))
        stress(:)=stress(:)+strkin(:)
 
       press2=Pint       
c        if(idnode.eq.0)write(2001,'(6e15.6)')engke,virtot,vircon,
c     x      vir2,Pint,(2.d0*engke+vir2)/(3.d0*volm)

      endif

c     adjust long range corrections and number density

      elrc=elrc0*(volm0/volm)
      virlrc=virlrc0*(volm0/volm)
      do kk=1,ntpatm
        dens(kk)=dens0(kk)*(volm0/volm)
      enddo

      return
      
      end subroutine pimd_npt_piglet_nm

      
      function getke(natms,idnode,mxnode)

c*********************************************************************
c
c     dl_poly routine to calculate system kinetic energy
c
c     copyright daresbury laboratory
c     author - w.smith october 2002
c     amended - w. smith january 2005 : f90 conversion
c
c*********************************************************************

      implicit none

      integer natms,idnode,mxnode,i,iatm0,iatm1
      real(8) getke,engke
      
c     block indices
      
      iatm0=nbeads*((idnode*natms)/mxnode)
      iatm1=nbeads*(((idnode+1)*natms)/mxnode)

      engke=0.d0
      
      do i=1,iatm1-iatm0
        engke=engke+rzmass(i)*(pxx(i)**2+pyy(i)**2+pzz(i)**2)
      enddo

      if(mxnode.gt.1)then
        
        buffer(1)=engke
        call gdsum(buffer(1),1,buffer(2))
        engke=buffer(1)
        
      endif

      getke=0.5d0*engke

      return
      end function getke

      function getkec(natms,idnode,mxnode)

c*********************************************************************
c
c     dl_poly routine to calculate system centroid kinetic energy
c
c     copyright daresbury laboratory
c     author - w.smith october 2002
c     amended - w. smith january 2005 : f90 conversion
c
c*********************************************************************

      implicit none

      integer natms,idnode,mxnode,i,iatm0,iatm1
      real(8) getkec,engkec
      
c     block indices
      
      iatm0=nbeads*((idnode*natms)/mxnode)
      iatm1=nbeads*(((idnode+1)*natms)/mxnode)

      engkec=0.d0
      
      do i=1,iatm1-iatm0,nbeads
        engkec=engkec+rzmass(i)*(pxx(i)**2+pyy(i)**2+pzz(i)**2)
      enddo

      if(mxnode.gt.1)then
        
        buffer(1)=engkec
        call gdsum(buffer(1),1,buffer(2))
        engkec=buffer(1)
        
      endif

      getkec=0.5d0*engkec

      return
      end function getkec


c  testing virial as in i-Pi
      subroutine getvir2(natms,idnode,mxnode,vir2)

c*********************************************************************
c
c     dl_poly routine to calculate system kinetic energy
c
c     copyright daresbury laboratory
c     author - w.smith october 2002
c     amended - w. smith january 2005 : f90 conversion
c
c*********************************************************************

      implicit none

      integer :: i,j
      integer :: natms,idnode,mxnode,iatm0,iatm1
      real(8) :: xxc,yyc,zzc
      real(8) :: vir2

c     block indices
      
      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode

      vir2=0.d0

      do i=iatm0,iatm1

        xxc=0.d0
        yyc=0.d0
        zzc=0.d0

        do j=1,nbeads

          xxc=xxc+xxx((j-1)*natms+i)
          yyc=yyc+yyy((j-1)*natms+i)
          zzc=zzc+zzz((j-1)*natms+i)

        enddo
        xxc=xxc/nbeads
        yyc=yyc/nbeads
        zzc=zzc/nbeads

        do j=1,nbeads

          vir2=vir2+(xxx((j-1)*natms+i)-xxc)*fxx((j-1)*natms+i)
          vir2=vir2+(yyy((j-1)*natms+i)-yyc)*fyy((j-1)*natms+i)
          vir2=vir2+(zzz((j-1)*natms+i)-zzc)*fzz((j-1)*natms+i)

        enddo

      enddo

      if(mxnode.gt.1)then

        buffer(1)=vir2
        call gdsum(buffer(1),1,buffer(2))
        vir2=buffer(1)

      endif

      return

      end subroutine getvir2

c**********************************************************************
            
      subroutine thermo_chain_nm
     x  (idnode,mxnode,natms,nstart,tstep,taut,temp)
      
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
      integer, intent(in) :: idnode,mxnode,natms,nstart
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
      
      qqq(1)=4.d0*boltz*temp*taut**2
      do k=2,nbeads
        qqq(k)=boltz*temp/nmfreq(k)**2
      enddo
      
      do m=1,nsy
        
        do n=1,ncyc
          
c     define scaled time steps
          
          fstep=tstep*wsy(m)/dble(ncyc)
          hstep=0.5d0*tstep*wsy(m)/dble(ncyc)
          qstep=0.25d0*tstep*wsy(m)/dble(ncyc)
          
c     integrate thermostat momenta
          
          do k=nstart,nbeads
            
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
                
                do j=nchain-1,2,-1
                  
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
          
          do k=nstart,nbeads
            
            do i=k,iatm1-iatm0,nbeads
              
              do j=1,nchain
                
                etx(j,i)=etx(j,i)+fstep*pcx(j,i)/qqq(k)
                ety(j,i)=ety(j,i)+fstep*pcy(j,i)/qqq(k)
                etz(j,i)=etz(j,i)+fstep*pcz(j,i)/qqq(k)
                
              enddo
              
            enddo
            
          enddo
          
c     apply thermostats to bead momenta
          
          do k=nstart,nbeads
            
            do i=k,iatm1-iatm0,nbeads
              
              pxx(i)=pxx(i)*(1.d0-(fstep/qqq(k))*pcx(1,i))
              pyy(i)=pyy(i)*(1.d0-(fstep/qqq(k))*pcy(1,i))
              pzz(i)=pzz(i)*(1.d0-(fstep/qqq(k))*pcz(1,i))
              
            enddo
            
          enddo
        
c     integrate thermostat momenta
          
          do k=nstart,nbeads
            
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
      
      end subroutine thermo_chain_nm

      subroutine thermostat_energy_nm
     x  (idnode,mxnode,natms,nbeads,nchain,nstart,qq1,qqk,taut,temp,
     x  engthe)

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
      
      integer, intent(in ) :: idnode,mxnode,natms,nbeads,nchain,nstart
      real(8), intent(in)  :: qq1,qqk,taut,temp
      real(8), intent(out) :: engthe

      integer i,j,k,iatm0,iatm1
      real(8), allocatable :: qqq(:)

      iatm0=nbeads*((idnode*natms)/mxnode)
      iatm1=nbeads*(((idnode+1)*natms)/mxnode)
      
      allocate (qqq(nbeads))
       
c     define thermostat masses
      
      qqq(1)=4.d0*boltz*temp*taut**2
      do k=2,nbeads
        qqq(k)=boltz*temp/nmfreq(k)**2
      enddo
      engthe=0.d0
      
      do k=nstart,nbeads
        
        do i=k,iatm1-iatm0,nbeads
          
          do j=1,nchain
            
            engthe=engthe+boltz*temp*(etx(j,i)+ety(j,i)+etz(j,i))+
     x        (pcx(j,i)**2+pcy(j,i)**2+pcz(j,i)**2)/(2.d0*qqq(k))
            
          enddo
          
        enddo
        
c        qqq=qqk
        
      enddo

      end subroutine thermostat_energy_nm
      end module vv_pimd_module
      
