      module ensemble_tools_module

c***********************************************************************
c     
c     dl_poly module defining tools for ensemble simulations
c     copyright - daresbury laboratory
c     author    - w. smith    sep 2003
c     
c     nhc_part and nhc_baro subroutines are added for NVT and
c     NPT ensembles together with Nose-Hoover Chain thermostat/barostat
c
c     copyright - M.R.Momeni and F.A. Shakib
c     authors   - M.R.Momeni and F.A.Shakib 2021
c
c     Method Development and Materials Simulation Laboratory
c
c***********************************************************************

      use setup_module, only: nrite
      use config_module
      use core_shell_module
      use property_module
      use rigid_body_module
      use utility_module
      use nhc_module

      contains
      
      function getmass(natms,idnode,mxnode)

c*********************************************************************
c
c     dl_poly routine to calculate total system mass
c
c     copyright daresbury laboratory
c     author - w.smith october 2002
c     amended - w.smith january 2005 : f90 conversion
c
c*********************************************************************
      
      implicit none

      integer natms,idnode,mxnode,i,iatm0,iatm1
      real(8) getmass
      
      iatm0 = (idnode*natms)/mxnode+1
      iatm1 = ((idnode+1)*natms)/mxnode

      getmass=0.d0

      do i=iatm0,iatm1

        getmass=getmass+weight(i)

      enddo

      if(mxnode.gt.1)then
        
        buffer(1)=getmass
        call gdsum(buffer(1),1,buffer(2))
        getmass=buffer(1)
        
      endif
      
      return
      end function getmass

      subroutine getcom(natms,idnode,mxnode,totmas,com)

c*********************************************************************
c
c     dl_poly routine to calculate system centre of mass
c
c     copyright daresbury laboratory
c     author - w.smith october 2002
c     amended - w.smith january 2005 : f90 conversion
c
c*********************************************************************

      implicit none

      integer natms,idnode,mxnode,i,iatm0,iatm1
      real(8) totmas
      
      real(8) com(3)
      
      iatm0 = (idnode*natms)/mxnode+1
      iatm1 = ((idnode+1)*natms)/mxnode

      com(1)=0.d0
      com(2)=0.d0
      com(3)=0.d0
      
      do i=iatm0,iatm1
        
        com(1)=com(1)+weight(i)*xxx(i)
        com(2)=com(2)+weight(i)*yyy(i)
        com(3)=com(3)+weight(i)*zzz(i)

      enddo

      if(mxnode.gt.1) call gdsum(com,3,buffer)

      com(1)=com(1)/totmas
      com(2)=com(2)/totmas
      com(3)=com(3)/totmas

      return
      end subroutine getcom

      subroutine getvom(natms,idnode,mxnode,totmas,vom)

c*********************************************************************
c
c     dl_poly routine to calculate system centre of mass
c
c     copyright daresbury laboratory
c     author - w.smith october 2002
c     amended - w.smith january 2005 : f90 conversion
c
c*********************************************************************

      implicit none

      integer natms,idnode,mxnode,i,iatm0,iatm1
      real(8) totmas
      
      real(8) vom(3)
      
      iatm0 = (idnode*natms)/mxnode+1
      iatm1 = ((idnode+1)*natms)/mxnode

      vom(1)=0.d0
      vom(2)=0.d0
      vom(3)=0.d0
      
      do i=iatm0,iatm1
        
        vom(1)=vom(1)+weight(i)*vxx(i)
        vom(2)=vom(2)+weight(i)*vyy(i)
        vom(3)=vom(3)+weight(i)*vzz(i)

      enddo

      if(mxnode.gt.1) call gdsum(vom,3,buffer)

      vom(1)=vom(1)/totmas
      vom(2)=vom(2)/totmas
      vom(3)=vom(3)/totmas

      return
      end subroutine getvom

      subroutine nvtscale
     x  (idnode,mxnode,natms,engke,sigma,tstep,qmass,taut,chit,conint)

c*********************************************************************
c
c     dl_poly routine to integrate and apply NVT thermostat
c
c     copyright daresbury laboratory
c     author - w.smith october 2002
c     amended - w.smith january 2005 : f90 conversion
c
c*********************************************************************

      implicit none

      integer idnode,mxnode,natms,i,iatm0,iatm1
      real(8) engke,sigma,tstep,qmass,chit,conint,scale,taut

      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode

c     calculate kinetic energy
      
      engke=getkin(natms,idnode,mxnode)

c     update chit to 1/2 step
      
      chit=chit+tstep*(engke-sigma)/qmass

c     thermostat the velocities
      
      scale=exp(-tstep*chit)

      do i=iatm0,iatm1
        
        vxx(i)=scale*vxx(i)
        vyy(i)=scale*vyy(i)
        vzz(i)=scale*vzz(i)
        
      enddo
      engke=engke*scale**2

c     update chi to full step
      
      conint=conint+tstep*chit*qmass/taut**2

c     update chit to full step
      
      chit=chit+tstep*(engke-sigma)/qmass
      
      return
      end subroutine nvtscale


      subroutine nhc_part
     x     (idnode,mxnode,natms,nchain,nrespa,engke,sigma,sigma_nhc,
     x     tstep,qmass_t,qmass_part,taut)

c*********************************************************************
c
c     dl_poly quantum routine to integrate and apply NHC thermostat
c     together with NVT ensemble 
c
c     copyright - M.R.Momeni and F.A.Shakib
c     authors   - M.R.Momeni and F.A.Shakib 2021
c
c     Method Development and Materials Simulation Laboratory
c
c*********************************************************************

      implicit none

      integer idnode,mxnode,natms,i,j,iatm0,iatm1
      integer nw,nsy,respa,nrespa,nchain
      parameter (nsy=7)
      real(8) engke,sigma,tstep,qmass_t,scale,taut
      real(8) sigma_nhc,qmass_part
      real(8) dt2,dt4,dt8,weight(nsy),kT,dNkT

c     define block indices

      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode

c     define NHC time variables

      dt2=tstep
      dt4=0.5d0*tstep
      dt8=0.25d0*tstep
      kT=2.d0*sigma_nhc
      dNkT=2.d0*sigma

c     Assign Suzuki-Yoshida weights

      weight(1)=0.784513610477560d0
      weight(2)=0.235573213359357d0
      weight(3)=-1.17767998417887d0
      weight(4)=1.315186321d0
      weight(5)=-1.17767998417887d0
      weight(6)=0.235573213359357d0
      weight(7)=0.784513610477560d0

c     Start Suzuki-Yoshida scheme      

       DO nw=nsy,1,-1

c     Start RESPA loop

       do respa=1,nrespa

c     Calculate kinetic energy
      
      engke=getkin(natms,idnode,mxnode)

c     Start 1st Suzuki-Yoshida scheme

        peta(nchain) = peta(nchain) +
     x  (weight(nw)*dt4/nrespa)*((peta(nchain-1)**2/qmass_part)-kT)
         
        do j=nchain-1,1,-1
           peta(j) = peta(j)*
     x     exp((-weight(nw)*dt8/nrespa)*peta(j+1)/qmass_part)

           if (j.eq.1) then
           peta(j)=peta(j) + (weight(nw)*dt4/nrespa)*((2.d0*engke)-dNkT)
           else
           peta(j) = peta(j) + 
     x     (weight(nw)*dt4/nrespa)*((peta(j-1)**2/qmass_part)-kT)
           endif
           peta(j) = peta(j)*
     x     exp((-weight(nw)*dt8/nrespa)*peta(j+1)/qmass_part)
        enddo      
  
c     Thermostat the velocities
      
      scale=exp((-weight(nw)*dt2/nrespa)*peta(1)/qmass_t)

      do i=iatm0,iatm1
        
        vxx(i)=scale*vxx(i)
        vyy(i)=scale*vyy(i)
        vzz(i)=scale*vzz(i)
        
      enddo

      do j=nchain,1,-1
         if (j.eq.1) then
           eta_nhc(j) = eta_nhc(j) +
     x     (weight(nw)*dt2/nrespa)*peta(j)/qmass_t
         else
           eta_nhc(j) = eta_nhc(j) +
     x     (weight(nw)*dt2/nrespa)*peta(j)/qmass_part
         endif
      enddo

c     Start 2nd Suzuki-Yoshida scheme

      do j=1,nchain-1
           peta(j) = peta(j)*
     x     exp((-weight(nw)*dt8/nrespa)*peta(j+1)/qmass_part)
           if (j.eq.1) then
           peta(j)=peta(j) + (weight(nw)*dt4/nrespa)*((2.d0*engke)-dNkT)
           else
           peta(j) = peta(j) + 
     x     (weight(nw)*dt4/nrespa)*((peta(j-1)**2/qmass_part)-kT)
           endif
           peta(j) = peta(j)*
     x     exp((-weight(nw)*dt8/nrespa)*peta(j+1)/qmass_part)
       enddo      

       peta(nchain) = peta(nchain) + 
     x (weight(nw)*dt4/nrespa)*((peta(nchain-1)**2/qmass_part)-kT)
      
       ENDDO

       ENDDO
      
      return

      end subroutine nhc_part

      subroutine nhc_baro
     x     (idnode,mxnode,natms,nchain,nrespa,sigma_nhc,
     x     tstep,volm_mass,qmass_baro,taup,v_epsilon)

c*********************************************************************
c
c     dl_poly quantum routine to integrate and apply NHC barostat
c     on volume dynamics variables in NPT ensemble 
c
c     copyright - M.R.Momeni and F.A.Shakib
c     authors   - M.R.Momeni and F.A.Shakib 2021
c
c     Method Development and Materials Simulation Laboratory
c
c*********************************************************************

      implicit none

      integer idnode,mxnode,natms,i,j,iatm0,iatm1
      integer nw,nsy,respa,nrespa,nchain
      parameter (nsy=7)
      real(8) tstep,scale,taup
      real(8) sigma_nhc,volm_mass,qmass_baro,v_epsilon
      real(8) dt2,dt4,dt8,weight(nsy),kT


c     define NHC time variables

      dt2=tstep
      dt4=0.5d0*tstep
      dt8=0.25d0*tstep
      kT=2.0*sigma_nhc

c     Assign Suzuki-Yoshida weights

      weight(1)=0.784513610477560d0
      weight(2)=0.235573213359357d0
      weight(3)=-1.17767998417887d0
      weight(4)=1.315186321d0
      weight(5)=-1.17767998417887d0
      weight(6)=0.235573213359357d0
      weight(7)=0.784513610477560d0

c     Start Suzuki-Yoshida scheme      

      DO nw=nsy,1,-1

c     Start RESPA loop

       do respa=1,nrespa

c     Start 1st Suzuki-Yoshida scheme

        pksi(nchain) = pksi(nchain) +
     x  (weight(nw)*dt4/nrespa)*((pksi(nchain-1)**2/qmass_baro)-kT)
         
        do j=nchain-1,1,-1
           pksi(j) = pksi(j)*
     x     exp((-weight(nw)*dt8/nrespa)*pksi(j+1)/qmass_baro)

           if (j.eq.1) then
           pksi(j) = pksi(j) + (weight(nw)*dt4/nrespa)
     x     *((volm_mass*v_epsilon**2)-kT)
           else
           pksi(j) = pksi(j) + 
     x     (weight(nw)*dt4/nrespa)*((pksi(j-1)**2/qmass_baro)-kT)
           endif
           pksi(j) = pksi(j)*
     x     exp((-weight(nw)*dt8/nrespa)*pksi(j+1)/qmass_baro)
        enddo      
  

c     Thermostat the velocities
      
        scale=exp((-weight(nw)*dt2/nrespa)*pksi(1)/qmass_baro)

        v_epsilon=scale*v_epsilon

        do j=nchain,1,-1
           ksi(j) = ksi(j) + 
     x     (weight(nw)*dt2/nrespa)*pksi(j)/qmass_baro
        enddo

c     Start 2nd Suzuki-Yoshida scheme

      do j=1,nchain-1
           pksi(j) = pksi(j)*
     x     exp((-weight(nw)*dt8/nrespa)*pksi(j+1)/qmass_baro)
           if (j.eq.1) then
           pksi(j) = pksi(j) + (weight(nw)*dt4/nrespa)
     x     *((volm_mass*v_epsilon**2)-kT)
           else
           pksi(j) = pksi(j) + 
     x     (weight(nw)*dt4/nrespa)*((pksi(j-1)**2/qmass_baro)-kT)
           endif
           pksi(j) = pksi(j)*
     x     exp((-weight(nw)*dt8/nrespa)*pksi(j+1)/qmass_baro)
       enddo      

       pksi(nchain) = pksi(nchain) + 
     x (weight(nw)*dt4/nrespa)*((pksi(nchain-1)**2/qmass_baro)-kT)
      
       ENDDO

       ENDDO
       
      return

      end subroutine nhc_baro

      subroutine nptscale_t
     x  (idnode,mxnode,natms,engke,temp,sigma,tstep,pmass,qmass,taut,
     x  chip,chit,conint)

c*********************************************************************
c
c     dl_poly routine to integrate and apply NPT thermostat
c
c     copyright daresbury laboratory
c     author - w.smith july 2005
c
c*********************************************************************

      implicit none

      integer idnode,mxnode,natms,i,iatm0,iatm1
      real(8) engke,temp,sigma,tstep,pmass,qmass,chip,chit,conint,scale
      real(8) taut

      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode

c     calculate kinetic energy
      
      engke=getkin(natms,idnode,mxnode)

c     update chit to 1/2 step
      
      chit=chit+0.5d0*tstep*(2.d0*(engke-sigma)+
     x  pmass*chip**2-boltz*temp)/qmass

c     thermostat the velocities
      
      scale=exp(-tstep*chit)

      do i=iatm0,iatm1
        
        vxx(i)=scale*vxx(i)
        vyy(i)=scale*vyy(i)
        vzz(i)=scale*vzz(i)
        
      enddo
      engke=engke*scale**2

c     update chi to full step
      
      conint=conint+tstep*chit*(qmass/taut**2+boltz*temp)

c     update chit to full step
      
      chit=chit+0.5d0*tstep*(2.d0*(engke-sigma)+
     x  pmass*chip**2-boltz*temp)/qmass

      return
      end subroutine nptscale_t

      subroutine nptscale_p
     x  (idnode,mxnode,natms,engke,tstep,pmass,chip,chit,
     x  volm,press,vircon,virtot)
      
c*********************************************************************
c     
c     dl_poly routine to integrate and apply NPT barostat
c     
c     copyright daresbury laboratory
c     author - w.smith october 2002
c     amended - w.smith january 2005 : f90 conversion
c     
c*********************************************************************
      
      implicit none

      integer idnode,mxnode,natms,i,iatm0,iatm1
      real(8) engke,tstep,pmass,chip,press,vircon,virtot
      real(8) volm,scale,chit

      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode

c     propagate chip to 1/2 step
      
      chip=chip+0.5d0*tstep*(((2.d0*engke-virtot-vircon)-
     x  3.d0*press*volm)/pmass-chip*chit)

c     barostat the velocities
      
      scale=exp(-tstep*chip)

      do i=iatm0,iatm1
        
        vxx(i)=scale*vxx(i)
        vyy(i)=scale*vyy(i)
        vzz(i)=scale*vzz(i)
        
      enddo
      engke=engke*scale**2

c     update volume parameter to full step
      
      volm=volm*exp(3.d0*tstep*chip)

c     update chip to full step
      
      chip=chip+0.5d0*tstep*(((2.d0*engke-virtot-vircon)-
     x  3.d0*press*volm)/pmass-chip*chit)
      
      return
      end subroutine nptscale_p

      subroutine nvtqscl
     x  (idnode,mxnode,ntfree,ngrp,engfke,engtrn,engrot,sigma,
     x   tstep,qmass,taut,chit,conint)

c*********************************************************************
c
c     dl_poly routine to integrate and apply NVT thermostat
c     to atomic, group and quaternion momenta
c
c     copyright daresbury laboratory
c     author - w.smith april 2005
c
c*********************************************************************

      implicit none

      integer idnode,mxnode,ntfree,ngrp,i,j,igrp1,igrp2,ifre1,ifre2
      integer ig
      real(8) engke,engtrn,engrot,engfke,sigma,tstep,qmass,chit,taut
      real(8) conint,scale

c     group block indices
      
      igrp1=(idnode*ngrp)/mxnode+1
      igrp2=((idnode+1)*ngrp)/mxnode

c     free atom block indices
      
      ifre1=(idnode*ntfree)/mxnode+1
      ifre2=((idnode+1)*ntfree)/mxnode
      
c     calculate kinetic energy
      
      engfke=getkinf(ntfree,idnode,mxnode)
      call getking(ngrp,idnode,mxnode,engtrn,engrot)
      engke=engfke+engtrn+engrot

c     update chit to 1/2 step
      
      chit=chit+tstep*(engke-sigma)/qmass

c     thermostat scale parameter
      
      scale=exp(-tstep*chit)

c     thermostat free atoms

      do j=ifre1,ifre2

         i=lstfre(j)
         vxx(i)=scale*vxx(i)
         vyy(i)=scale*vyy(i)
         vzz(i)=scale*vzz(i)
        
      enddo

c     thermostat rigid body velocities

      do ig=igrp1,igrp2
         
         omx(ig)=scale*omx(ig)
         omy(ig)=scale*omy(ig)
         omz(ig)=scale*omz(ig)
         gvxx(ig)=scale*gvxx(ig)
         gvyy(ig)=scale*gvyy(ig)
         gvzz(ig)=scale*gvzz(ig)
         
      enddo

c     scale kinetic energy

      engfke=engfke*scale**2
      engtrn=engtrn*scale**2
      engrot=engrot*scale**2

c     update chi to full step
      
      conint=conint+tstep*chit*qmass/taut**2

c     update chit to full step
      
      engke=engfke+engtrn+engrot
      chit=chit+tstep*(engke-sigma)/qmass

      return
      end subroutine nvtqscl

      subroutine nptqscl_t
     x  (idnode,mxnode,ntfree,ngrp,engfke,engtrn,engrot,temp,sigma,
     x  tstep,pmass,qmass,taut,chip,chit,conint)

c*********************************************************************
c     
c     dl_poly routine to integrate and apply NPT thermostat
c     to atomic, group and quaternion momenta
c     
c     copyright daresbury laboratory
c     author - w.smith april 2005
c     
c*********************************************************************

      implicit none

      integer idnode,mxnode,ntfree,ngrp,i,j,igrp1,igrp2,ifre1,ifre2
      integer ig
      real(8) engke,engtrn,engrot,engfke,sigma,tstep,qmass,chit,taut
      real(8) conint,scale,chip,pmass,temp

c     group block indices
      
      igrp1=(idnode*ngrp)/mxnode+1
      igrp2=((idnode+1)*ngrp)/mxnode

c     free atom block indices
      
      ifre1=(idnode*ntfree)/mxnode+1
      ifre2=((idnode+1)*ntfree)/mxnode
      
c     calculate kinetic energy
      
      engfke=getkinf(ntfree,idnode,mxnode)
      call getking(ngrp,idnode,mxnode,engtrn,engrot)
      engke=engfke+engtrn+engrot

c     update chit to 1/2 tstep
      
      chit=chit+0.5d0*tstep*(2.d0*(engke-sigma)+
     x  pmass*chip**2-boltz*temp)/qmass
      
c     thermostat scale parameter
      
      scale=exp(-tstep*chit)
      
c     thermostat free atoms
      
      do j=ifre1,ifre2
        
        i=lstfre(j)
        vxx(i)=scale*vxx(i)
        vyy(i)=scale*vyy(i)
        vzz(i)=scale*vzz(i)
        
      enddo
      
c     thermostat rigid body velocities
      
      do ig=igrp1,igrp2
        
        omx(ig)=scale*omx(ig)
        omy(ig)=scale*omy(ig)
        omz(ig)=scale*omz(ig)
        gvxx(ig)=scale*gvxx(ig)
        gvyy(ig)=scale*gvyy(ig)
        gvzz(ig)=scale*gvzz(ig)
        
      enddo
      
c     scale kinetic energy
      
      engfke=engfke*scale**2
      engtrn=engtrn*scale**2
      engrot=engrot*scale**2
      
c     update chi to full tstep
      
      conint=conint+tstep*chit*(qmass/taut**2+boltz*temp)
      
c     update chit to full tstep
      
      engke=engfke+engtrn+engrot
      chit=chit+0.5d0*tstep*(2.d0*(engke-sigma)+
     x  pmass*chip**2-boltz*temp)/qmass
      
      return
      end subroutine nptqscl_t

      subroutine nptqscl_p
     x  (idnode,mxnode,ntfree,ngrp,engfke,engtrn,tstep,pmass,
     x  chip,chit,volm,press,vircon,virtot,vircom)
      
c*********************************************************************
c     
c     dl_poly routine to integrate and apply NPT barostat
c     for system with atomic sites and rigid bodies
c     
c     copyright daresbury laboratory
c     author - w.smith october 2002
c     amended - w.smith january 2005 : f90 conversion
c     
c*********************************************************************
      
      implicit none

      integer idnode,mxnode,i,ngrp,ntfree,igrp1,igrp2,ifre1,ifre2
      integer j,ig
      real(8) engke,tstep,pmass,chip,press,vircon,virtot
      real(8) vircom,volm,scale,engtrn,engfke,chit

c     group block indices
      
      igrp1=(idnode*ngrp)/mxnode+1
      igrp2=((idnode+1)*ngrp)/mxnode

c     free atom block indices
      
      ifre1=(idnode*ntfree)/mxnode+1
      ifre2=((idnode+1)*ntfree)/mxnode

c     propagate chip to 1/2 tstep
      
      engke=engfke+engtrn
      chip=chip+0.5d0*tstep*(((2.d0*engke-virtot-vircon-vircom)-
     x  3.d0*press*volm)/pmass-chip*chit)
      
c     barostat the free atom velocities
      
      scale=exp(-tstep*chip)
      
      do j=ifre1,ifre2
        
        i=lstfre(j)
        vxx(i)=scale*vxx(i)
        vyy(i)=scale*vyy(i)
        vzz(i)=scale*vzz(i)
        
      enddo
      
c     barostat the group translational velocities
      
      do ig=igrp1,igrp2
        
        gvxx(ig)=scale*gvxx(ig)
        gvyy(ig)=scale*gvyy(ig)
        gvzz(ig)=scale*gvzz(ig)
        
      enddo
      
c     scale kinetic energy
      
      engfke=engfke*scale**2
      engtrn=engtrn*scale**2
      
c     update volume parameter to full tstep
      
      volm=volm*exp(3.d0*tstep*chip)
      
c     update chip to full tstep
      
      engke=engfke+engtrn
      chip=chip+0.5d0*tstep*(((2.d0*engke-virtot-vircon-vircom)-
     x  3.d0*press*volm)/pmass-chip*chit)
      
      return
      end subroutine nptqscl_p

      subroutine nstscale_t
     x  (idnode,mxnode,natms,mode,engke,temp,sigma,tstep,
     x  pmass,qmass,taut,chit,conint)

c*********************************************************************
c
c     dl_poly routine to integrate and apply NST thermostat
c
c     copyright daresbury laboratory
c     author - w.smith july 2005
c
c*********************************************************************

      implicit none

      integer idnode,mxnode,natms,mode,i,iatm0,iatm1
      real(8) engke,temp,sigma,tstep,pmass,qmass,chip2,chit,conint,scale
      real(8) taut,fac(0:3)
      data fac/9.d0,3.d0,2.d0,5.d0/

      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode

c     calculate kinetic energy
      
      chip2=sdot0(9,eta,eta)
      if(mode.eq.2)chip2=chip2-eta(1)**2
      engke=getkin(natms,idnode,mxnode)

c     update chit to 1/2 step
      
      chit=chit+0.5d0*tstep*(2.d0*(engke-sigma)+
     x  pmass*chip2-boltz*temp*fac(mode))/qmass

c     thermostat the velocities
      
      scale=exp(-tstep*chit)
      do i=iatm0,iatm1
        
        vxx(i)=scale*vxx(i)
        vyy(i)=scale*vyy(i)
        vzz(i)=scale*vzz(i)
        
      enddo
      engke=engke*scale**2

c     update chi to full step
      
      conint=conint+tstep*chit*(qmass/taut**2+boltz*temp*fac(mode))

c     update chit to full step
      
      chit=chit+0.5d0*tstep*(2.d0*(engke-sigma)+
     x  pmass*chip2-boltz*temp*fac(mode))/qmass

      return
      end subroutine nstscale_t

      subroutine nstscale_p
     x  (idnode,mxnode,natms,mode,tstep,pmass,chit,press,volm)
      
c*********************************************************************
c     
c     dl_poly routine to integrate and apply NPT anisotropic barostat
c     
c     copyright daresbury laboratory
c     author - w.smith october 2002
c     amended - w.smith january 2005 : f90 conversion
c
c*********************************************************************
      
      implicit none

      integer idnode,mxnode,natms,mode,i,iatm0,iatm1
      real(8) tstep,pmass,press,volm,txx,tyy,tzz,chit
      real(8) strkin(9),uni(9),celp(10)

      data uni/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/
      
      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode

c     calculate kinetic contribution to stress tensor

      call kinstress(natms,idnode,mxnode,strkin)

c     propagate barostat momentum to 1/2 step

      do i=1,9
        eta(i)=eta(i)+0.5d0*tstep*((stress(i)+strkin(i)-
     x    press*volm*uni(i))/pmass-chit*eta(i))
      enddo
      if(mode.gt.0)then
        eta(3)=0.d0
        eta(6)=0.d0
        eta(7)=0.d0
        eta(8)=0.d0
        if(mode.lt.3)then
          eta(2)=0.d0
          eta(4)=0.d0
          if(mode.eq.2)then
            eta(1)=0.5d0*(eta(1)+eta(5))
            eta(5)=eta(1)
          endif
        endif
      endif

c     barostat the velocities
      
      do i=iatm0,iatm1
        
        txx=vxx(i)
        tyy=vyy(i)
        tzz=vzz(i)
        vxx(i)=txx-tstep*(eta(1)*txx+eta(4)*tyy+eta(7)*tzz)
        vyy(i)=tyy-tstep*(eta(2)*txx+eta(5)*tyy+eta(8)*tzz)
        vzz(i)=tzz-tstep*(eta(3)*txx+eta(6)*tyy+eta(9)*tzz)

      enddo

c     new cell vectors
          
      call cell_update(tstep,cell,eta)

c     update volume to full time step

      call dcell(cell,celp)
      volm=celp(10)

c     calculate kinetic energy and contribution to stress tensor

      call kinstress(natms,idnode,mxnode,strkin)

c     propagate barostat momentum to full step

      do i=1,9
        eta(i)=eta(i)+0.5d0*tstep*((stress(i)+strkin(i)-
     x    press*volm*uni(i))/pmass-chit*eta(i))
      enddo
      if(mode.gt.0)then
        eta(3)=0.d0
        eta(6)=0.d0
        eta(7)=0.d0
        eta(8)=0.d0
        if(mode.lt.3)then
          eta(2)=0.d0
          eta(4)=0.d0
          if(mode.eq.2)then
            eta(1)=0.5d0*(eta(1)+eta(5))
            eta(5)=eta(1)
          endif
        endif
      endif
      
      return
      end subroutine nstscale_p

      subroutine nstqscl_t
     x  (idnode,mxnode,ntfree,ngrp,mode,engfke,engtrn,engrot,temp,
     x  sigma,tstep,pmass,qmass,taut,chit,conint)

c*********************************************************************
c
c     dl_poly routine to integrate and apply NPT thermostat
c     to atomic, group and quaternion momenta
c
c     copyright daresbury laboratory
c     author - w.smith april 2005
c
c*********************************************************************

      implicit none

      integer idnode,mxnode,ntfree,ngrp,i,j,igrp1,igrp2,ifre1,ifre2
      integer ig,mode
      real(8) engke,engtrn,engrot,engfke,sigma,tstep,qmass,chit,taut
      real(8) conint,scale,chip2,pmass,temp,fac(0:3)
      data fac/9.d0,3.d0,2.d0,5.d0/

c     group block indices
      
      igrp1=(idnode*ngrp)/mxnode+1
      igrp2=((idnode+1)*ngrp)/mxnode

c     free atom block indices
      
      ifre1=(idnode*ntfree)/mxnode+1
      ifre2=((idnode+1)*ntfree)/mxnode
      
c     calculate kinetic energy
      
      chip2=sdot0(9,eta,eta)
      if(mode.eq.2)chip2=chip2-eta(1)**2
      engfke=getkinf(ntfree,idnode,mxnode)
      call getking(ngrp,idnode,mxnode,engtrn,engrot)
      engke=engfke+engtrn+engrot

c     update chit to 1/2 step
      
      chit=chit+0.5d0*tstep*(2.d0*(engke-sigma)+
     x  pmass*chip2-boltz*temp*fac(mode))/qmass

c     thermostat scale parameter
      
      scale=exp(-tstep*chit)

c     thermostat free atoms

      do j=ifre1,ifre2

         i=lstfre(j)
         vxx(i)=scale*vxx(i)
         vyy(i)=scale*vyy(i)
         vzz(i)=scale*vzz(i)
        
      enddo

c     thermostat rigid body velocities

      do ig=igrp1,igrp2
         
         omx(ig)=scale*omx(ig)
         omy(ig)=scale*omy(ig)
         omz(ig)=scale*omz(ig)
         gvxx(ig)=scale*gvxx(ig)
         gvyy(ig)=scale*gvyy(ig)
         gvzz(ig)=scale*gvzz(ig)
         
      enddo

c     scale kinetic energy

      engfke=engfke*scale**2
      engtrn=engtrn*scale**2
      engrot=engrot*scale**2

c     update chi to full step
      
      conint=conint+tstep*chit*(qmass/taut**2+boltz*temp*fac(mode))

c     update chit to full step
      
      engke=engfke+engtrn+engrot
      chit=chit+0.5d0*tstep*(2.d0*(engke-sigma)+
     x  pmass*chip2-boltz*temp*fac(mode))/qmass

      return
      end subroutine nstqscl_t

      subroutine nstqscl_p
     x  (idnode,mxnode,ntfree,ngrp,mode,tstep,pmass,chit,press,volm)
      
c*********************************************************************
c     
c     dl_poly routine to integrate and apply NPT anisotropic barostat
c     for system with atomic sites and rigid bodies
c     
c     copyright daresbury laboratory
c     author - w.smith may 2005
c
c*********************************************************************
      
      implicit none

      integer idnode,mxnode,ntfree,ngrp,i,igrp1,igrp2,ifre1,ifre2,ig,j
      integer mode
      real(8) tstep,pmass,press,volm,txx,tyy,tzz,chit
      real(8) strkin(9),strgrp(9),uni(9),celp(10)
      
      data uni/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/
      
c     group block indices
      
      igrp1=(idnode*ngrp)/mxnode+1
      igrp2=((idnode+1)*ngrp)/mxnode

c     free atom block indices
      
      ifre1=(idnode*ntfree)/mxnode+1
      ifre2=((idnode+1)*ntfree)/mxnode

c     propagate barostat momentum to 1/2 step
      
      do i=1,9
        eta(i)=eta(i)+0.5d0*tstep*((stress(i)-
     x    press*volm*uni(i))/pmass-chit*eta(i))
      enddo
      if(mode.gt.0)then
        eta(3)=0.d0
        eta(6)=0.d0
        eta(7)=0.d0
        eta(8)=0.d0
        if(mode.lt.3)then
          eta(2)=0.d0
          eta(4)=0.d0
          if(mode.eq.2)then
            eta(1)=0.5d0*(eta(1)+eta(5))
            eta(5)=eta(1)
          endif
        endif
      endif

c     subtract kinetic contribution from stress tensor

      call kinstressf(ntfree,idnode,mxnode,strkin)        
      call kinstressg(ngrp,idnode,mxnode,strgrp)

      do i=1,9
        stress(i)=stress(i)-strkin(i)-strgrp(i)
      enddo
      
c     barostat the free atom velocities
      
      do j=ifre1,ifre2
        
        i=lstfre(j)
        txx=vxx(i)
        tyy=vyy(i)
        tzz=vzz(i)
        vxx(i)=txx-tstep*(eta(1)*txx+eta(4)*tyy+eta(7)*tzz)
        vyy(i)=tyy-tstep*(eta(2)*txx+eta(5)*tyy+eta(8)*tzz)
        vzz(i)=tzz-tstep*(eta(3)*txx+eta(6)*tyy+eta(9)*tzz)

      enddo

c     barostat the group translational velocities

      do ig=igrp1,igrp2
         
         txx=gvxx(ig)
         tyy=gvyy(ig)
         tzz=gvzz(ig)
         gvxx(ig)=txx-tstep*(eta(1)*txx+eta(4)*tyy+eta(7)*tzz)
         gvyy(ig)=tyy-tstep*(eta(2)*txx+eta(5)*tyy+eta(8)*tzz)
         gvzz(ig)=tzz-tstep*(eta(3)*txx+eta(6)*tyy+eta(9)*tzz)
         
      enddo

c     new cell vectors

      call cell_update(tstep,cell,eta)
      
c     new system volume
      
      call dcell(cell,celp)
      volm=celp(10)
      
c     add new kinetic contribution to stress tensor

      call kinstressf(ntfree,idnode,mxnode,strkin)        
      call kinstressg(ngrp,idnode,mxnode,strgrp)

      do i=1,9
        stress(i)=stress(i)+strkin(i)+strgrp(i)
      enddo
      
c     propagate barostat momentum to full step

      do i=1,9
        eta(i)=eta(i)+0.5d0*tstep*((stress(i)-
     x    press*volm*uni(i))/pmass-chit*eta(i))
      enddo
      if(mode.gt.0)then
        eta(3)=0.d0
        eta(6)=0.d0
        eta(7)=0.d0
        eta(8)=0.d0
        if(mode.lt.3)then
          eta(2)=0.d0
          eta(4)=0.d0
          if(mode.eq.2)then
            eta(1)=0.5d0*(eta(1)+eta(5))
            eta(5)=eta(1)
          endif
        endif
      endif

      return
      end subroutine nstqscl_p

      subroutine nstqscl_t2
     x  (idnode,mxnode,ntfree,ngrp,mode,engfke,engtrn,engrot,temp,
     x  sigma,tstep,pmass,qmass,taut,chit,conint,strkin,strgrp)

c*********************************************************************
c
c     dl_poly routine to integrate and apply NPT thermostat
c     to atomic, group and quaternion momenta
c
c     copyright daresbury laboratory
c     author - w.smith april 2005
c
c*********************************************************************

      implicit none

      integer idnode,mxnode,ntfree,ngrp,i,j,igrp1,igrp2,ifre1,ifre2
      integer ig,mode
      real(8) engke,engtrn,engrot,engfke,sigma,tstep,qmass,chit,taut
      real(8) conint,scale,chip2,pmass,temp,fac(0:3)
      real(8) strkin(9),strgrp(9)
      data fac/9.d0,3.d0,2.d0,5.d0/

c     group block indices
      
      igrp1=(idnode*ngrp)/mxnode+1
      igrp2=((idnode+1)*ngrp)/mxnode

c     free atom block indices
      
      ifre1=(idnode*ntfree)/mxnode+1
      ifre2=((idnode+1)*ntfree)/mxnode
      
c     remove old kinetic term from stress tensor

      do i=1,9
        stress(i)=stress(i)-strkin(i)-strgrp(i)
      enddo

c     calculate kinetic energy
      
      chip2=sdot0(9,eta,eta)
      if(mode.eq.2)chip2=chip2-eta(1)**2
      engfke=getkinf(ntfree,idnode,mxnode)
      call getking(ngrp,idnode,mxnode,engtrn,engrot)
      engke=engfke+engtrn+engrot

c     update chit to 1/2 step
      
      chit=chit+0.5d0*tstep*(2.d0*(engke-sigma)+
     x  pmass*chip2-boltz*temp*fac(mode))/qmass

c     thermostat scale parameter
      
      scale=exp(-tstep*chit)

c     thermostat free atoms

      do j=ifre1,ifre2

         i=lstfre(j)
         vxx(i)=scale*vxx(i)
         vyy(i)=scale*vyy(i)
         vzz(i)=scale*vzz(i)
        
      enddo

c     thermostat rigid body velocities

      do ig=igrp1,igrp2
         
         omx(ig)=scale*omx(ig)
         omy(ig)=scale*omy(ig)
         omz(ig)=scale*omz(ig)
         gvxx(ig)=scale*gvxx(ig)
         gvyy(ig)=scale*gvyy(ig)
         gvzz(ig)=scale*gvzz(ig)
         
      enddo

c     scale kinetic energy

      engfke=engfke*scale**2
      engtrn=engtrn*scale**2
      engrot=engrot*scale**2

c     scale kinetic energy tensors

      do i=1,9

        strkin(i)=strkin(i)*scale**2
        strgrp(i)=strgrp(i)*scale**2

      enddo

c     update chi to full step
      
      conint=conint+tstep*chit*(qmass/taut**2+boltz*temp*fac(mode))

c     update chit to full step
      
      engke=engfke+engtrn+engrot
      chit=chit+0.5d0*tstep*(2.d0*(engke-sigma)+
     x  pmass*chip2-boltz*temp*fac(mode))/qmass

c     add new kinetic terms to stress tensor

      do i=1,9
        stress(i)=stress(i)+strkin(i)+strgrp(i)
      enddo

      return
      end subroutine nstqscl_t2

      subroutine nstqscl_p2
     x  (idnode,mxnode,ntfree,ngrp,mode,tstep,pmass,chit,press,volm,
     x  strkin,strgrp)
      
c*********************************************************************
c     
c     dl_poly routine to integrate and apply NPT anisotropic barostat
c     for system with atomic sites and rigid bodies
c     
c     copyright daresbury laboratory
c     author - w.smith may 2005
c
c*********************************************************************
      
      implicit none

      integer idnode,mxnode,ntfree,ngrp,i,igrp1,igrp2,ifre1,ifre2,ig,j
      integer mode
      real(8) tstep,pmass,press,volm,txx,tyy,tzz,chit
      real(8) strkin(9),strgrp(9),uni(9),celp(10)
      
      data uni/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/
      
c     group block indices
      
      igrp1=(idnode*ngrp)/mxnode+1
      igrp2=((idnode+1)*ngrp)/mxnode

c     free atom block indices
      
      ifre1=(idnode*ntfree)/mxnode+1
      ifre2=((idnode+1)*ntfree)/mxnode

c     propagate barostat momentum to 1/2 step
      
      do i=1,9
        eta(i)=eta(i)+0.5d0*tstep*((stress(i)-
     x    press*volm*uni(i))/pmass-chit*eta(i))
      enddo
      if(mode.gt.0)then
        eta(3)=0.d0
        eta(6)=0.d0
        eta(7)=0.d0
        eta(8)=0.d0
        if(mode.lt.3)then
          eta(2)=0.d0
          eta(4)=0.d0
          if(mode.eq.2)then
            eta(1)=0.5d0*(eta(1)+eta(5))
            eta(5)=eta(1)
          endif
        endif
      endif

c     subtract kinetic contribution from stress tensor

      do i=1,9
        stress(i)=stress(i)-strkin(i)-strgrp(i)
      enddo
      
c     barostat the free atom velocities
      
      do j=ifre1,ifre2
        
        i=lstfre(j)
        txx=vxx(i)
        tyy=vyy(i)
        tzz=vzz(i)
        vxx(i)=txx-tstep*(eta(1)*txx+eta(4)*tyy+eta(7)*tzz)
        vyy(i)=tyy-tstep*(eta(2)*txx+eta(5)*tyy+eta(8)*tzz)
        vzz(i)=tzz-tstep*(eta(3)*txx+eta(6)*tyy+eta(9)*tzz)

      enddo

c     barostat the group translational velocities

      do ig=igrp1,igrp2
         
         txx=gvxx(ig)
         tyy=gvyy(ig)
         tzz=gvzz(ig)
         gvxx(ig)=txx-tstep*(eta(1)*txx+eta(4)*tyy+eta(7)*tzz)
         gvyy(ig)=tyy-tstep*(eta(2)*txx+eta(5)*tyy+eta(8)*tzz)
         gvzz(ig)=tzz-tstep*(eta(3)*txx+eta(6)*tyy+eta(9)*tzz)
         
      enddo

c     new cell vectors

      call cell_update(tstep,cell,eta)
      
c     new system volume
      
      call dcell(cell,celp)
      volm=celp(10)
      
c     add new kinetic contribution to stress tensor

      call kinstressf(ntfree,idnode,mxnode,strkin)        
      call kinstressg(ngrp,idnode,mxnode,strgrp)

      do i=1,9
        stress(i)=stress(i)+strkin(i)+strgrp(i)
      enddo
      
c     propagate barostat momentum to full step

      do i=1,9
        eta(i)=eta(i)+0.5d0*tstep*((stress(i)-
     x    press*volm*uni(i))/pmass-chit*eta(i))
      enddo
      if(mode.gt.0)then
        eta(3)=0.d0
        eta(6)=0.d0
        eta(7)=0.d0
        eta(8)=0.d0
        if(mode.lt.3)then
          eta(2)=0.d0
          eta(4)=0.d0
          if(mode.eq.2)then
            eta(1)=0.5d0*(eta(1)+eta(5))
            eta(5)=eta(1)
          endif
        endif
      endif

      return
      end subroutine nstqscl_p2

      subroutine cell_update(tstep,cell,eta)

c***********************************************************************
c     
c     dlpoly utility to update the cell vectors in the hoover 
c     nst algorithms (velocity verlet version)
c
c     copyright daresbury laboratory
c     author      w.smith  july 2005
c     
c**********************************************************************

      implicit none

      integer i
      real(8) tstep,cell(9),eta(9),ctmp(9),uni(9)
      data uni/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/
      
      do i=1,9
        ctmp(i)=uni(i)+tstep*eta(i)
      enddo

      call mat_mul(ctmp,cell,cell)
      
      return
      end subroutine cell_update

      subroutine cell_propagate(tstep,cell,eta)

c***********************************************************************
c     
c     dlpoly utility to update the cell vectors in the hoover 
c     nst algorithms (leapfrog version)
c
c     copyright daresbury laboratory
c     author      w.smith  july 2005
c     
c**********************************************************************

      implicit none

      integer i
      real(8) tstep
      real(8) cell(9),eta(9),aaa(9),bbb(9),uni(9)
      data uni/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/

      do i=1,9
        aaa(i)=tstep*eta(i)
      enddo
      
      call mat_mul(aaa,aaa,bbb)
      
      do i=1,9
        bbb(i)=uni(i)+aaa(i)+0.5d0*bbb(i)
      enddo
      
      call mat_mul(bbb,cell,cell)

      return
      end subroutine cell_propagate

      subroutine nstqmtk_p
     x  (idnode,mxnode,ntfree,ngrp,mode,tstep,pmass,chit,press,volm,
     x  engfke,engtrn,engrot,temp,sigma)
      
c*********************************************************************
c     
c     dl_poly routine to integrate and apply NPT anisotropic barostat
c     of martyna tobias and klein to atomic, group and quaternion 
c     system with atomic sites and rigid bodies
c     
c     copyright daresbury laboratory
c     author - w.smith may 2005
c
c*********************************************************************
      
      implicit none

      integer idnode,mxnode,ntfree,ngrp,i,igrp1,igrp2,ifre1,ifre2,ig,j
      integer mode
      real(8) tstep,pmass,press,volm,txx,tyy,tzz,chit,temp,sigma,degfre
      real(8) engtke,engfke,engtrn,engrot,trace
      real(8) strkin(9),strgrp(9),uni(9),ctmp(9)

      data uni/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/
      
      degfre=2.d0*sigma/(temp*boltz)

c     group block indices
      
      igrp1=(idnode*ngrp)/mxnode+1
      igrp2=((idnode+1)*ngrp)/mxnode

c     free atom block indices
      
      ifre1=(idnode*ntfree)/mxnode+1
      ifre2=((idnode+1)*ntfree)/mxnode

c     calculate kinetic contribution to stress tensor

      call kinstressf(ntfree,idnode,mxnode,strkin)        
      call kinstressg(ngrp,idnode,mxnode,strgrp)

c     propagate barostat momentum to 1/2 step

      call invert(cell,ctmp,volm)
      volm=abs(volm)
      engtke=2.d0*(engfke+engtrn+engrot)/degfre
      do i=1,9
        eta(i)=eta(i)+0.5d0*tstep*((stress(i)+strkin(i)+strgrp(i)+
     x    (engtke-press*volm)*uni(i))/pmass-chit*eta(i))
      enddo
      if(mode.gt.0)then
        eta(3)=0.d0
        eta(6)=0.d0
        eta(7)=0.d0
        eta(8)=0.d0
        if(mode.lt.3)then
          eta(2)=0.d0
          eta(4)=0.d0
          if(mode.eq.2)then
            eta(1)=0.5d0*(eta(1)+eta(5))
            eta(5)=eta(1)
          endif
        endif
      endif

c     barostat the free atom velocities
      
      trace=(eta(1)+eta(5)+eta(9))/degfre

      do j=ifre1,ifre2
        
        i=lstfre(j)
        txx=vxx(i)
        tyy=vyy(i)
        tzz=vzz(i)
        vxx(i)=txx-tstep*(eta(1)*txx+eta(4)*tyy+eta(7)*tzz+trace)
        vyy(i)=tyy-tstep*(eta(2)*txx+eta(5)*tyy+eta(8)*tzz+trace)
        vzz(i)=tzz-tstep*(eta(3)*txx+eta(6)*tyy+eta(9)*tzz+trace)

      enddo

c     barostat the group translational velocities

      do ig=igrp1,igrp2
         
         txx=gvxx(ig)
         tyy=gvyy(ig)
         tzz=gvzz(ig)
         gvxx(ig)=txx-tstep*(eta(1)*txx+eta(4)*tyy+eta(7)*tzz+trace)
         gvyy(ig)=tyy-tstep*(eta(2)*txx+eta(5)*tyy+eta(8)*tzz+trace)
         gvzz(ig)=tzz-tstep*(eta(3)*txx+eta(6)*tyy+eta(9)*tzz+trace)
         
      enddo

c     update volume to full time step

      volm=volm*exp(tstep*(eta(1)+eta(5)+eta(9)))

c     calculate kinetic contribution to stress tensor

      call kinstressf(ntfree,idnode,mxnode,strkin)        
      call kinstressg(ngrp,idnode,mxnode,strgrp)

c     calculate new kinetic energy

      engfke=getkinf(ntfree,idnode,mxnode)
      call getking(ngrp,idnode,mxnode,engtrn,engrot)
      engtke=2.d0*(engfke+engtrn+engrot)/degfre

c     propagate barostat momentum to full step

      do i=1,9
        eta(i)=eta(i)+0.5d0*tstep*((stress(i)+strkin(i)+strgrp(i)+
     x    (engtke-press*volm)*uni(i))/pmass-chit*eta(i))
      enddo
      if(mode.gt.0)then
        eta(3)=0.d0
        eta(6)=0.d0
        eta(7)=0.d0
        eta(8)=0.d0
        if(mode.lt.3)then
          eta(2)=0.d0
          eta(4)=0.d0
          if(mode.eq.2)then
            eta(1)=0.5d0*(eta(1)+eta(5))
            eta(5)=eta(1)
          endif
        endif
      endif

      return
      end subroutine nstqmtk_p

      subroutine kinstr(idnode,mxnode,natms,tstep)
      
c***********************************************************************
c     
c     dlpoly routine to calculate the kinetic energy contribution to
c     the stress tensor
c     
c     assumes velocities are half-timestep behind forces
c     
c     replicated data version / block data
c     
c     copyright daresbury laboratory 1994
c     author t.forester may 1994
c     amended t.forester dec 1994 : block data
c     
c***********************************************************************

      implicit none

      integer idnode,mxnode,natms,i,iatm1,iatm2
      real(8) tstep,vxt,vyt,vzt
      
c     block indices

      iatm1 = (idnode*natms)/mxnode + 1
      iatm2 = ((idnode+1)*natms)/mxnode

      do i = iatm1,iatm2

        if(rmass(i).gt.0.d0) then

          vxt = vxx(i)+fxx(i)*rmass(i)*tstep*0.5d0
          vyt = vyy(i)+fyy(i)*rmass(i)*tstep*0.5d0
          vzt = vzz(i)+fzz(i)*rmass(i)*tstep*0.5d0

          stress(1)=stress(1)-weight(i)*vxt*vxt
          stress(2)=stress(2)-weight(i)*vxt*vyt
          stress(3)=stress(3)-weight(i)*vxt*vzt
          stress(4)=stress(4)-weight(i)*vyt*vxt
          stress(5)=stress(5)-weight(i)*vyt*vyt
          stress(6)=stress(6)-weight(i)*vyt*vzt
          stress(7)=stress(7)-weight(i)*vzt*vxt
          stress(8)=stress(8)-weight(i)*vzt*vyt
          stress(9)=stress(9)-weight(i)*vzt*vzt

        endif

      enddo
      
      return
      end subroutine kinstr
      
      function getkin(natms,idnode,mxnode)

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
      real(8) getkin,engke
      
      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode

      engke=0.d0
      
      do i=iatm0,iatm1
        engke=engke+weight(i)*(vxx(i)**2+vyy(i)**2+vzz(i)**2)
      enddo

      if(mxnode.gt.1)then
        
        buffer(1)=engke
        call gdsum(buffer(1),1,buffer(2))
        engke=buffer(1)
        
      endif

      getkin=0.5d0*engke

      return
      end function getkin

      function getkinf(ntfree,idnode,mxnode)

c*********************************************************************
c
c     dl_poly routine to calculate kinetic energy of atoms not in
c     rigid bodies
c
c     copyright daresbury laboratory
c     author - w.smith october 2002
c     amended - w. smith january 2005 : f90 conversion
c
c*********************************************************************

      implicit none

      integer ntfree,idnode,mxnode,i,j,ifre0,ifre1
      real(8) getkinf,engke
      
      ifre0=(idnode*ntfree)/mxnode+1
      ifre1=((idnode+1)*ntfree)/mxnode

      engke=0.d0
      
      do j=ifre0,ifre1
        
        i=lstfre(j)
        engke=engke+weight(i)*(vxx(i)**2+vyy(i)**2+vzz(i)**2)

      enddo

      if(mxnode.gt.1)then
        
        buffer(1)=engke
        call gdsum(buffer(1),1,buffer(2))
        engke=buffer(1)
        
      endif

      getkinf=0.5d0*engke

      return
      end function getkinf

      subroutine getking(ngrp,idnode,mxnode,engtrn,engrot)

c*********************************************************************
c
c     dl_poly routine to calculate system kinetic energy
c     for rigid groups only
c
c     copyright daresbury laboratory
c     author - m.leslie february 2003
c     amended - w.smith january 2005 : f90 conversion
c
c*********************************************************************

      implicit none

      integer ngrp,idnode,mxnode,igrp1,igrp2,ig,id
      real(8) engtrn,engrot

      igrp1=(idnode*ngrp)/mxnode+1
      igrp2=((idnode+1)*ngrp)/mxnode

      engtrn=0.d0
      engrot=0.d0
      do ig=igrp1,igrp2

        id=lstgtp(ig)

c     group kinetic energy

        engtrn=engtrn+
     x    gmass(id)*(gvxx(ig)**2+gvyy(ig)**2+gvzz(ig)**2)

c     rotational kinetic energy
        
        engrot=engrot+(rotinx(id,1)*omx(ig)**2
     x    +rotiny(id,1)*omy(ig)**2
     x    +rotinz(id,1)*omz(ig)**2)

      enddo

      if(mxnode.gt.1) then
        
        buffer(5)=engtrn
        buffer(6)=engrot
        call  gdsum(buffer(5),2,buffer(1))
        engtrn=buffer(5)
        engrot=buffer(6)
        
      endif

      engtrn=0.5d0*engtrn
      engrot=0.5d0*engrot

      return
      end subroutine getking

      function getkint(ngrp,idnode,mxnode)

c*********************************************************************
c
c     dl_poly routine to calculate translational kinetic energy
c     for rigid groups only
c
c     copyright daresbury laboratory
c     author  - w.smith october 2005
c
c*********************************************************************

      implicit none

      integer ngrp,idnode,mxnode,igrp1,igrp2,ig,id
      real(8) engtrn,getkint

      igrp1=(idnode*ngrp)/mxnode+1
      igrp2=((idnode+1)*ngrp)/mxnode

      engtrn=0.d0
      do ig=igrp1,igrp2

        id=lstgtp(ig)

c     group kinetic energy

        engtrn=engtrn+
     x    gmass(id)*(gvxx(ig)**2+gvyy(ig)**2+gvzz(ig)**2)

      enddo

      if(mxnode.gt.1)then
        
        buffer(1)=engtrn
        call gdsum(buffer(1),1,buffer(2))
        engtrn=buffer(1)
        
      endif
      
      getkint=0.5d0*engtrn

      return
      end function getkint

      function getkinr(ngrp,idnode,mxnode)

c*********************************************************************
c
c     dl_poly routine to calculate rotational kinetic energy
c     for rigid groups only
c
c     copyright daresbury laboratory
c     author  - w.smith october 2005
c
c*********************************************************************

      implicit none

      integer ngrp,idnode,mxnode,igrp1,igrp2,ig,id
      real(8) engrot,getkinr

      igrp1=(idnode*ngrp)/mxnode+1
      igrp2=((idnode+1)*ngrp)/mxnode

      engrot=0.d0
      do ig=igrp1,igrp2

        id=lstgtp(ig)

c     rotational kinetic energy
        
        engrot=engrot+(rotinx(id,1)*omx(ig)**2
     x    +rotiny(id,1)*omy(ig)**2
     x    +rotinz(id,1)*omz(ig)**2)

      enddo

      if(mxnode.gt.1)then
        
        buffer(1)=engrot
        call gdsum(buffer(1),1,buffer(2))
        engrot=buffer(1)
        
      endif
      
      getkinr=0.5d0*engrot
      
      return
      end function getkinr

      subroutine kinstress(natms,idnode,mxnode,stresh)

c*********************************************************************
c     
c     dl_poly routine to calculate kinetic contribution to the 
c     stress tensor
c     
c     copyright daresbury laboratory
c     author - w.smith november 2002
c
c*********************************************************************
      
      implicit none

      integer natms,idnode,mxnode,iatm0,iatm1,i
      real(8) stresh(9)
      
      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode

c     initialise stress tensor

      do i=1,9
        stresh(i)=0.d0
      enddo

c     kinetic contribution to stress tensor
        
      do i=iatm0,iatm1
        
        stresh(1)=stresh(1)+weight(i)*vxx(i)*vxx(i)
        stresh(2)=stresh(2)+weight(i)*vxx(i)*vyy(i)
        stresh(3)=stresh(3)+weight(i)*vxx(i)*vzz(i)
        stresh(5)=stresh(5)+weight(i)*vyy(i)*vyy(i)
        stresh(6)=stresh(6)+weight(i)*vyy(i)*vzz(i)
        stresh(9)=stresh(9)+weight(i)*vzz(i)*vzz(i)
        
      enddo

      stresh(4)=stresh(2)
      stresh(7)=stresh(3)
      stresh(8)=stresh(6)

c     global sum of stress tensor
      
      if(mxnode.gt.1) call gdsum(stresh,9,buffer)
        
      return
      end subroutine kinstress

      subroutine kinstressg(ngrp,idnode,mxnode,stresh)

c*********************************************************************
c     
c     dl_poly routine to calculate kinetic contribution to the 
c     stress tensor
c     
c     copyright daresbury laboratory
c     author - m.leslie february 2003
c
c*********************************************************************
      
      integer ngrp,idnode,mxnode,igrp1,igrp2,ig,id
      real(8) stresh(9)
      
      igrp1=(idnode*ngrp)/mxnode+1
      igrp2=((idnode+1)*ngrp)/mxnode

c     initialise stress tensor

      do i=1,9
        stresh(i)=0.d0
      enddo

c     kinetic contribution to stress tensor
      
      do ig=igrp1,igrp2

        id=lstgtp(ig)
        stresh(1)=stresh(1)+gmass(id)*gvxx(ig)*gvxx(ig)
        stresh(2)=stresh(2)+gmass(id)*gvxx(ig)*gvyy(ig)
        stresh(3)=stresh(3)+gmass(id)*gvxx(ig)*gvzz(ig)
        stresh(5)=stresh(5)+gmass(id)*gvyy(ig)*gvyy(ig)
        stresh(6)=stresh(6)+gmass(id)*gvyy(ig)*gvzz(ig)
        stresh(9)=stresh(9)+gmass(id)*gvzz(ig)*gvzz(ig)
        
      enddo
      
      stresh(4)=stresh(2)
      stresh(7)=stresh(3)
      stresh(8)=stresh(6)

c     global sum of stress tensor
      
      if(mxnode.gt.1) call gdsum(stresh,9,buffer)
        
      return
      end subroutine kinstressg

      subroutine getkins(natms,idnode,mxnode,getkin)

c*********************************************************************
c
c     dl_poly routine to calculate system kinetic energy
c
c     copyright daresbury laboratory
c     author - w.smith october 2002
c     amended - w.smith january 2005: f90 conversion
c
c*********************************************************************
     
      implicit none

      integer natms,idnode,mxnode,iatm0,iatm1,i
      real(8) getkin,engke
      
      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode

      engke=0.d0
      
      do i=iatm0,iatm1
        engke=engke+weight(i)*(vxx(i)**2+vyy(i)**2+vzz(i)**2)
      enddo

      if(mxnode.gt.1)then
        
        buffer(1)=engke
        call gdsum(buffer(1),1,buffer(2))
        engke=buffer(1)
        
      endif

      getkin=0.5d0*engke

      return
      end subroutine getkins

      subroutine kinstressf(ntfree,idnode,mxnode,stresh)

c*********************************************************************
c     
c     dl_poly routine to calculate kinetic contribution to the 
c     stress tensor for unconstrained atoms
c     
c     copyright daresbury laboratory
c     author - m.leslie february 2003
c     amended - w.smith january 2005: f90 conversion
c
c*********************************************************************
      
      implicit none

      integer ntfree,idnode,mxnode,i,ifre1,ifre2,ifre
      real(8) stresh(9)
      
      ifre1=(idnode*ntfree)/mxnode+1
      ifre2=((idnode+1)*ntfree)/mxnode

c     initialise stress tensor

      do i=1,9
        stresh(i)=0.d0
      enddo

c     kinetic contribution to stress tensor
        
      do ifre=ifre1,ifre2
        
        i=lstfre(ifre)
        stresh(1)=stresh(1)+weight(i)*vxx(i)*vxx(i)
        stresh(2)=stresh(2)+weight(i)*vxx(i)*vyy(i)
        stresh(3)=stresh(3)+weight(i)*vxx(i)*vzz(i)
        stresh(5)=stresh(5)+weight(i)*vyy(i)*vyy(i)
        stresh(6)=stresh(6)+weight(i)*vyy(i)*vzz(i)
        stresh(9)=stresh(9)+weight(i)*vzz(i)*vzz(i)
        
      enddo

      stresh(4)=stresh(2)
      stresh(7)=stresh(3)
      stresh(8)=stresh(6)

c     global sum of stress tensor
      
      if(mxnode.gt.1) call gdsum(stresh,9,buffer)
        
      return
      end subroutine kinstressf
      
      subroutine nvtscale_shl
     x  (idnode,mxnode,ntshl,shlke,sigma_shl,tstep,qmass_shl,
     x   taut,chit_shl,conint)

c*********************************************************************
c
c     dl_poly routine to integrate and apply NVT thermostat
c     thermostats the core-shell relative motion
c
c     copyright daresbury laboratory
c     author - w.smith october 2002
c     amended - w.smith january 2005 : f90 conversion
c     adapted - d. quigley      2006 : core-shell motion
c
c*********************************************************************

      implicit none

      integer idnode,mxnode,ntshl,i,ishl1,ishl2,j,k,m
      real(8) shlke,sigma_shl,tstep,qmass_shl,chit_shl,conint
      real(8) dvx,dvy,dvz,tmx,tmy,tmz,rmu,scale,taut

      ishl1=(idnode*ntshl)/mxnode+1
      ishl2=((idnode+1)*ntshl)/mxnode

c     calculate kinetic energy
      
      call corshl(idnode,mxnode,ntshl,shlke)

c     update chit to 1/2 step
      
      chit_shl=chit_shl+tstep*(shlke-sigma_shl)/qmass_shl

c     thermostat the velocities
      
      scale=exp(-tstep*chit_shl)

      m=0
      do k=ishl1,ishl2
        
        m=m+1
        
        i=listshl(m,2)
        j=listshl(m,3)

        rmu=(weight(i)*weight(j))/(weight(i)+weight(j))
        
        if(rmu.gt.0.d0)then
          
          dvx=vxx(j)-vxx(i)
          dvy=vyy(j)-vyy(i)
          dvz=vzz(j)-vzz(i)
          
          tmx=weight(i)*vxx(i)+weight(j)*vxx(j)
          tmy=weight(i)*vyy(i)+weight(j)*vyy(j)
          tmz=weight(i)*vzz(i)+weight(j)*vzz(j)
          
          vxx(i)=tmx/(weight(i)+weight(j))-scale*rmu*dvx/weight(i)
          vxx(j)=tmx/(weight(i)+weight(j))+scale*rmu*dvx/weight(j)
          vyy(i)=tmy/(weight(i)+weight(j))-scale*rmu*dvy/weight(i)
          vyy(j)=tmy/(weight(i)+weight(j))+scale*rmu*dvy/weight(j)
          vzz(i)=tmz/(weight(i)+weight(j))-scale*rmu*dvz/weight(i)
          vzz(j)=tmz/(weight(i)+weight(j))+scale*rmu*dvz/weight(j)
          
        endif

      enddo

      shlke=shlke*scale**2

c     update chi to full step
      
      conint=conint+tstep*chit_shl*qmass_shl/taut**2

c     update chit to full step
      
      chit_shl=chit_shl+tstep*(shlke-sigma_shl)/qmass_shl

      if(mxnode.gt.1) call shlmerge(idnode,mxnode,ntshl)
      
      return
      end subroutine nvtscale_shl

      end module ensemble_tools_module

