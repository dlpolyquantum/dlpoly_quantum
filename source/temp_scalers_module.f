      module temp_scalers_module

c***********************************************************************
c     
c     dl_poly module for temperature scaling routines
c     
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2006
c     
c***********************************************************************
      
      use config_module
      use core_shell_module
      use error_module
      use rigid_body_module
      use setup_module
      use shake_module
      use utility_module
      
      contains

      subroutine quench(imcon,idnode,mxnode,natms,nscons,tolnce)

c*********************************************************************
c     
c     dl_poly subroutine for quenching the bond energies in the 
c     initial structure of a molecule defined by constraints
c     
c     copyright - daresbury laboratory 1992
c     author w.smith november 1992
c     
c*********************************************************************

      implicit none

      logical safe
      integer imcon,idnode,mxnode,natms,nscons,i,j,k,icyc
      integer fail
      real(8) tolnce,ddd,esig,vvv,ww1,ww2

      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: uxx(:),uyy(:),uzz(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)

      dimension fail(3)

      data fail/0,0,0/

c     allocate work arrays

      allocate (xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(1))
      allocate (uxx(mxatms),uyy(mxatms),uzz(mxatms),stat=fail(2))
      allocate (dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(3))

      do i=1,3
         if(fail(i).ne.0)call error(idnode,1770)
      enddo

c     calculate bond vectors

      do k=1,nscons
        
        i=listcon(k,2)
        j=listcon(k,3)
        
        dxt(k)=xxx(i)-xxx(j)
        dyt(k)=yyy(i)-yyy(j)
        dzt(k)=zzz(i)-zzz(j)
        
      enddo
      
      call images(imcon,0,1,nscons,cell,dxt,dyt,dzt)
      
c     normalise bond vectors
      
      do k=1,nscons
        
        ddd=sqrt(dxt(k)**2+dyt(k)**2+dzt(k)**2)
        
        dxt(k)=dxt(k)/ddd
        dyt(k)=dyt(k)/ddd
        dzt(k)=dzt(k)/ddd
        
      enddo
      
c     start of quenching cycle
      
      icyc=0
      safe=.false.
      do while(.not.safe.and.icyc.lt.mxshak)
        
        icyc=icyc+1

c     initialise velocity correction arrays
        
        do i=1,natms
          
          uxx(i)=0.d0
          uyy(i)=0.d0
          uzz(i)=0.d0
          
        enddo
        
c     calculate velocity corrections and error
        
        esig=0.d0

        do k=1,nscons
          
          i=listcon(k,2)
          j=listcon(k,3)
          
          vvv=dxt(k)*(vxx(i)-vxx(j))+dyt(k)*(vyy(i)-vyy(j))+
     x      dzt(k)*(vzz(i)-vzz(j))
          
          esig=max(esig,abs(vvv))
          
          ww1=weight(j)*vvv/(weight(i)+weight(j))
          ww2=weight(i)*vvv/(weight(i)+weight(j))
          uxx(i)=uxx(i)-ww1*dxt(k)
          uyy(i)=uyy(i)-ww1*dyt(k)
          uzz(i)=uzz(i)-ww1*dzt(k)
          uxx(j)=uxx(j)+ww2*dxt(k)
          uyy(j)=uyy(j)+ww2*dyt(k)
          uzz(j)=uzz(j)+ww2*dzt(k)
          
        enddo

        safe=(esig.lt.tolnce)
        
        if(mxnode.gt.1)call gstate(safe)
        
        if(.not.safe)then
          
c     transport velocity adjustments to other nodes
          
          if(mxnode.gt.1)then
            
            call shmove
     x        (idnode,mxnode,natms,lashap,lishap,uxx,uyy,uzz,
     x        xxt,yyt,zzt,buffer)
            
          endif
          
c     update velocities
          
          do k=1,nscons
            
            i=listcon(k,2)
            j=listcon(k,3)
            vxx(i)=vxx(i)+uxx(i)/dble(listme(i))
            vyy(i)=vyy(i)+uyy(i)/dble(listme(i))
            vzz(i)=vzz(i)+uzz(i)/dble(listme(i))
            vxx(j)=vxx(j)+uxx(j)/dble(listme(j))
            vyy(j)=vyy(j)+uyy(j)/dble(listme(j))
            vzz(j)=vzz(j)+uzz(j)/dble(listme(j))
            
          enddo
          
        endif
        
      enddo
      
c     error exit if quenching fails
      
      if(.not.safe)call error(idnode,70)
      
c     splice velocity arrays across nodes
      
      if(mxnode.gt.1) call splice
     x  (idnode,natms,listme,listot,vxx,vyy,vzz,buffer)

c     deallocate work arrays

      deallocate (xxt,yyt,zzt,stat=fail(1))
      deallocate (uxx,uyy,uzz,stat=fail(2))
      deallocate (dxt,dyt,dzt,stat=fail(3))
      
      return
      end subroutine quench

      subroutine quatqnch(idnode,imcon,mxnode,natms,ngrp)

c***********************************************************************
c     
c     dlpoly subroutine to convert atomic velocities to rigid body 
c     c.o.m. and angular velocity
c     
c     parallel replicated data version : block data
c     
c     copyright daresbury laboratory 1993.
c     author   - t.forester nov 1993.
c     amended  - t.forester dec 1994 : block data.
c     
c***********************************************************************

      implicit none
      
      integer idnode,imcon,mxnode,natms,ngrp,fail,ig,jr,id
      integer igrp1,igrp2,i,j
      real(8) rot,wxx,wyy,wzz
      
      dimension rot(9)
      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      
      data fail/0/
      
c     allocate work arrays 
      
      allocate (xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail)
      if(fail.ne.0)call error(idnode,1780)
      
c     block indices for groups
      
      igrp1=(idnode*ngrp)/mxnode+1
      igrp2=((idnode+1)*ngrp)/mxnode
      
c     translate atomic velocites to com velocity & angular velocity
      
      jr=0
      do ig=igrp1,igrp2
        
        gvxx(ig)=0.d0
        gvyy(ig)=0.d0
        gvzz(ig)=0.d0
        omx(ig)=0.d0
        omy(ig)=0.d0
        omz(ig)=0.d0
        
        id=lstgtp(ig)
        
        do j=1,numgsit(id)
          
          jr =jr+1
          i =lstrgd(jr)
          
c     centre of mass momentum
          
          gvxx(ig)=gvxx(ig)+weight(i)*vxx(i)
          gvyy(ig)=gvyy(ig)+weight(i)*vyy(i)
          gvzz(ig)=gvzz(ig)+weight(i)*vzz(i)
          
c     distance to c.o.m of molecule
          
          xxt(jr)=xxx(i)-gcmx(ig)
          yyt(jr)=yyy(i)-gcmy(ig)
          zzt(jr)=zzz(i)-gcmz(ig)
          
        enddo
        
c     centre of mass velocity
        
        gvxx(ig)=gvxx(ig)/gmass(id)
        gvyy(ig)=gvyy(ig)/gmass(id)
        gvzz(ig)=gvzz(ig)/gmass(id)
        
      enddo
      
      call images(imcon,0,1,jr,cell,xxt,yyt,zzt)
      
      jr=0
      do ig=igrp1,igrp2
        
c     rotational matrix
        
        rot(1)=q0(ig)**2+q1(ig)**2-q2(ig)**2-q3(ig)**2
        rot(2)=2.d0*(q1(ig)*q2(ig)-q0(ig)*q3(ig))
        rot(3)=2.d0*(q1(ig)*q3(ig)+q0(ig)*q2(ig))
        rot(4)=2.d0*(q1(ig)*q2(ig)+q0(ig)*q3(ig))
        rot(5)=q0(ig)**2-q1(ig)**2+q2(ig)**2-q3(ig)**2
        rot(6)=2.d0*(q2(ig)*q3(ig)-q0(ig)*q1(ig))
        rot(7)=2.d0*(q1(ig)*q3(ig)-q0(ig)*q2(ig))
        rot(8)=2.d0*(q2(ig)*q3(ig)+q0(ig)*q1(ig))
        rot(9)=q0(ig)**2-q1(ig)**2-q2(ig)**2+q3(ig)**2
        
c     angular momentum accumulators
        
        wxx=0.d0
        wyy=0.d0
        wzz=0.d0

        id=lstgtp(ig)

        do j=1,numgsit(id)
          
          jr =jr+1
          i =lstrgd(jr)
          
          wxx=wxx+weight(i)*(yyt(jr)*vzz(i)-zzt(jr)*vyy(i))
          wyy=wyy+weight(i)*(zzt(jr)*vxx(i)-xxt(jr)*vzz(i))
          wzz=wzz+weight(i)*(xxt(jr)*vyy(i)-yyt(jr)*vxx(i))
          
        enddo
        
c     angular velocity in body fixed frame
        
        omx(ig)=(rot(1)*wxx+rot(4)*wyy+rot(7)*wzz)*rotinx(id,2)
        omy(ig)=(rot(2)*wxx+rot(5)*wyy+rot(8)*wzz)*rotiny(id,2)
        omz(ig)=(rot(3)*wxx+rot(6)*wyy+rot(9)*wzz)*rotinz(id,2)
        
        jr=jr-numgsit(id)
        do j=1,numgsit(id)
          
          jr=jr +1
          i=lstrgd(jr)
          
c     site velocity in body frame 
          
          wxx=omy(ig)*gzz(id,j)-omz(ig)*gyy(id,j)
          wyy=omz(ig)*gxx(id,j)-omx(ig)*gzz(id,j)
          wzz=omx(ig)*gyy(id,j)-omy(ig)*gxx(id,j)
          
c     new atomic velocites in lab frame
          
          vxx(i)=rot(1)*wxx+rot(2)*wyy+rot(3)*wzz+gvxx(ig)
          vyy(i)=rot(4)*wxx+rot(5)*wyy+rot(6)*wzz+gvyy(ig)
          vzz(i)=rot(7)*wxx+rot(8)*wyy+rot(9)*wzz+gvzz(ig)
          
        enddo
        
      enddo
      
      if(mxnode.gt.1)then
        
        call merge(idnode,mxnode,ngrp,mxbuff,gvxx,gvyy,gvzz,buffer)
        call merge(idnode,mxnode,ngrp,mxbuff,omx,omy,omz,buffer)
        call merge1(idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)
        
      endif
      
c     deallocate work arrays
      
      deallocate (xxt,yyt,zzt,stat=fail)
      
      return
      end subroutine quatqnch

      subroutine vscaleg(idnode,mxnode,imcon,natms,ngrp,sigma)

c*********************************************************************
c     
c     dl_poly subroutine for scaling the velocity arrays to the
c     desired temperature
c     
c     zeroes angular momentum in non-periodic system.
c     
c     parallel replicated data version : block data
c     
c     copyright daresbury laboratory 1992.
c     author - w.smith july 1992
c     amended - t.forester oct 1993
c     amended - t.forester dec 1994 : block data
c     
c*********************************************************************
      
      implicit none
      
      integer idnode,mxnode,imcon,natms,ngrp,iatm1,iatm2,i
      real(8) sigma,roti,rotinv,cmx,cmy,cmz,cmvx,cmvy,cmvz,sysmas
      real(8) amx,amy,amz,det,scale,rsq,wxx,wyy,wzz,sumke
      
      dimension roti(9),rotinv(9)
      
c     block indices
      
      iatm1=(idnode*natms)/mxnode+1
      iatm2=((idnode+1)*natms)/mxnode
      
c     calculate centre of mass position and motion of the system
      
      cmx=0.d0
      cmy=0.d0
      cmz=0.d0
      cmvx=0.d0
      cmvy=0.d0
      cmvz=0.d0
      sysmas=0.d0
      
      do i=iatm1,iatm2
        
        if(lstfrz(i).eq.0.and.weight(i).gt.1.d-6)then
          
          cmx=cmx+weight(i)*xxx(i)
          cmy=cmy+weight(i)*yyy(i)
          cmz=cmz+weight(i)*zzz(i)
          sysmas=sysmas+weight(i)
          cmvx=cmvx+vxx(i)*weight(i)
          cmvy=cmvy+vyy(i)*weight(i)
          cmvz=cmvz+vzz(i)*weight(i)
          
        endif
        
      enddo
      
      if(mxnode.gt.1)then
        buffer(8)=sysmas
        buffer(9)=cmx
        buffer(10)=cmy
        buffer(11)=cmz
        buffer(12)=cmvx
        buffer(13)=cmvy
        buffer(14)=cmvz
        call gdsum(buffer(8),7,buffer(1))
        sysmas= buffer(8) 
        cmx=buffer(9) 
        cmy=buffer(10) 
        cmz=buffer(11) 
        cmvx=buffer(12) 
        cmvy=buffer(13) 
        cmvz=buffer(14) 
      endif
      
      cmx=cmx/sysmas
      cmy=cmy/sysmas
      cmz=cmz/sysmas
      
      cmvx=cmvx/sysmas
      cmvy=cmvy/sysmas
      cmvz=cmvz/sysmas
      
c     remove centre of mass motion  
      
      do i=1,natms
        
        if(lstfrz(i).eq.0.and.weight(i).gt.1.d-6)then
          
          vxx(i)=vxx(i)-cmvx
          vyy(i)=vyy(i)-cmvy
          vzz(i)=vzz(i)-cmvz
          
        else
          
          vxx(i)=0.d0
          vyy(i)=0.d0
          vzz(i)=0.d0
          
        endif
        
      enddo
      
c     zero angular momentum about centre of mass - non-periodic system
      
      if(imcon.eq.0)then
        
c     move to centre of mass origin
        
        do i=1,natms
          
          xxx(i)=xxx(i)-cmx
          yyy(i)=yyy(i)-cmy
          zzz(i)=zzz(i)-cmz
          
        enddo
        
c     angular momentum accumulators
        
        amx=0.d0
        amy=0.d0
        amz=0.d0
        
c     rotational inertia accumulators
        
        do i=1,9
          
          roti(i)=0.d0
          
        enddo
        
        do i=iatm1,iatm2
          
          amx=amx+weight(i)*(yyy(i)*vzz(i)-zzz(i)*vyy(i))
          amy=amy+weight(i)*(zzz(i)*vxx(i)-xxx(i)*vzz(i))
          amz=amz+weight(i)*(xxx(i)*vyy(i)-yyy(i)*vxx(i))
          
          rsq=xxx(i)**2+yyy(i)**2+zzz(i)**2
          roti(1)=roti(1)+weight(i)*(xxx(i)*xxx(i)-rsq)
          roti(2)=roti(2)+weight(i)* xxx(i)*yyy(i)
          roti(3)=roti(3)+weight(i)* xxx(i)*zzz(i)
          roti(5)=roti(5)+weight(i)*(yyy(i)*yyy(i)-rsq)
          roti(6)=roti(6)+weight(i)* yyy(i)*zzz(i)
          roti(9)=roti(9)+weight(i)*(zzz(i)*zzz(i)-rsq)
          
        enddo
        
c     complete rotational inertia matrix
        
        roti(4)=roti(2)
        roti(7)=roti(3)
        roti(8)=roti(6)
        
c     global sum
        
        if(mxnode.gt.1)then
          buffer(13)=amx
          buffer(14)=amy
          buffer(15)=amz
          do i=1,9
            buffer(15+i)=roti(i)
          enddo
          call gdsum(buffer(13),12,buffer(1))
          amx=buffer(13) 
          amy=buffer(14) 
          amz=buffer(15) 
          do i=1,9
            roti(i)=buffer(15+i)
          enddo
        endif
        
c     invert rotational inertia matrix
        
        call invert (roti,rotinv,det)
        
c     correction to angular velocity
        
        wxx=rotinv(1)*amx+rotinv(2)*amy+rotinv(3)*amz
        wyy=rotinv(4)*amx+rotinv(5)*amy+rotinv(6)*amz
        wzz=rotinv(7)*amx+rotinv(8)*amy+rotinv(9)*amz
        
c     correction to linear velocity
        
        do i=1,natms
          
          if(lstfrz(i).eq.0.and.weight(i).gt.1.d-6)then
            
            vxx(i)=vxx(i)+(wyy*zzz(i)-wzz*yyy(i))
            vyy(i)=vyy(i)+(wzz*xxx(i)-wxx*zzz(i))
            vzz(i)=vzz(i)+(wxx*yyy(i)-wyy*xxx(i))
            
          endif
          
        enddo
        
c     reset positions to original reference frame
        
        do i=1,natms
          
          xxx(i)=xxx(i)+cmx
          yyy(i)=yyy(i)+cmy
          zzz(i)=zzz(i)+cmz
          
        enddo
        
      endif
      
c     calculate temperature 
      
      sumke=0.d0
      
      do i=iatm1,iatm2
        
        sumke=sumke+weight(i)*
     x    (vxx(i)**2+vyy(i)**2+vzz(i)**2)
        
      enddo
      
      sumke=0.5d0*sumke
      if(mxnode.gt.1)then
        
        buffer(1)=sumke
        call gdsum(buffer(1),1,buffer(2))
        sumke=buffer(1)
        
      endif
      
c     apply temperature scaling
      
      if(sumke.gt.1.d-6)then

        scale=sqrt(sigma/sumke)

      else
        
        scale=1.d0
        
      endif
      
      do i=1,natms
        
        vxx(i)=scale*vxx(i)
        vyy(i)=scale*vyy(i)
        vzz(i)=scale*vzz(i)
        
      enddo

      if(ngrp.gt.0)then
        
        call quatqnch(idnode,imcon,mxnode,natms,ngrp)
        
      endif
      
      return
      end subroutine vscaleg

      subroutine shlqnch(idnode,mxnode,ntshl,temp)
      
c*********************************************************************
c     
c     dl_poly subroutine for quenching the internal bond energies
c     in ions defined by shell model
c     
c     copyright - daresbury laboratory 1994
c     author w.smith july  1994
c     
c*********************************************************************
      
      implicit none

      integer idnode,mxnode,ntshl,ishl1,ishl2,i,j,k,m
      real(8) temp,pke,rmu,dvx,dvy,dvz,tmx,tmy,tmz,scl

c     permitted core-shell internal kinetic energy 
      
      pke=boltz*temp*1.d-4

c     block indices

      ishl1 = (idnode*ntshl)/mxnode+1
      ishl2 = ((idnode+1)*ntshl)/mxnode

c     calculate core and shell velocities from total momentum
      
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
          
          scl=sqrt(pke/(rmu*(dvx*dvx+dvy*dvy+dvz*dvz)))
          
          tmx=weight(i)*vxx(i)+weight(j)*vxx(j)
          tmy=weight(i)*vyy(i)+weight(j)*vyy(j)
          tmz=weight(i)*vzz(i)+weight(j)*vzz(j)
          
          vxx(i)=tmx/(weight(i)+weight(j))-scl*rmu*dvx/weight(i)
          vxx(j)=tmx/(weight(i)+weight(j))+scl*rmu*dvx/weight(j)
          vyy(i)=tmy/(weight(i)+weight(j))-scl*rmu*dvy/weight(i)
          vyy(j)=tmy/(weight(i)+weight(j))+scl*rmu*dvy/weight(j)
          vzz(i)=tmz/(weight(i)+weight(j))-scl*rmu*dvz/weight(i)
          vzz(j)=tmz/(weight(i)+weight(j))+scl*rmu*dvz/weight(j)
          
        endif
        
      enddo

      if(mxnode.gt.1) call shlmerge(idnode,mxnode,ntshl)
      
      return
      end subroutine shlqnch
      
      subroutine regauss
     x  (idnode,imcon,mxnode,natms,ngrp,nscons,ntcons,
     x  ntshl,keyshl,sigma,temp,tolnce)
      
c***********************************************************************
c     
c     dl_poly subroutine for resetting the system velocities
c     
c     copyright - daresbury laboratory
c     author    - w. smith    may 2007
c     
c***********************************************************************

      implicit none

      integer idnode,imcon,mxnode,natms,ngrp,nscons
      integer ntcons,ntshl,i,k,keyshl
      real(8) temp,tolnce,sigma,rsq
      
c     set atomic velocities from gaussian distribution

      call gauss(natms,vxx,vyy,vzz)
      
      do i=1,natms
        
        rsq=sqrt(rmass(i))
        vxx(i)=vxx(i)*rsq
        vyy(i)=vyy(i)*rsq
        vzz(i)=vzz(i)*rsq
        
      enddo
      
      if(ntcons.gt.0)call quench
     x  (imcon,idnode,mxnode,natms,nscons,tolnce)
      
      if(ngrp.gt.0)call quatqnch
     x  (idnode,imcon,mxnode,natms,ngrp)
      
      if(keyshl.eq.1)then
        
        do k=1,4
          
          call vscaleg(idnode,mxnode,imcon,natms,ngrp,sigma)
          call shlqnch(idnode,mxnode,ntshl,temp)
          
        enddo
        
      else
        
        call vscaleg(idnode,mxnode,imcon,natms,ngrp,sigma)
        
      endif
      
      return
      end subroutine regauss
      
      subroutine impact(khit,natms,idnode,mxnode,ehit,xhit,yhit,zhit)
      
c*********************************************************************
c     
c     DLPOLY routinue for impacting a selected atom with a specified
c     recoil energy
c     
c     copyright daresbury laboratory
c     author w.smith september 2007
c     
c*********************************************************************
      
      use config_module
      use ensemble_tools_module
      
      implicit none
      
      integer i,khit,natms,idnode,mxnode,iatm0,iatm1
      real(8) ehit,vxo,vyo,vzo,xhit,yhit,zhit,fac,smass,vel
      
c     store original particle velocity
      
      vxo=vxx(khit)
      vyo=vyy(khit)
      vzo=vzz(khit)
      
c     determine recoil velocity
      
      vel=sqrt(2.d0*ehit/(weight(khit)*(xhit**2+yhit**2+zhit**2)))
      
c     reassign particle velocity
      
      vxx(khit)=vel*xhit
      vyy(khit)=vel*yhit
      vzz(khit)=vel*zhit
      
c     determine system mass
      
      smass=getmass(natms,idnode,mxnode)
      
c     calculate net system velocity
      
      vxo=(vxx(khit)-vxo)*weight(khit)/smass
      vyo=(vyy(khit)-vyo)*weight(khit)/smass
      vzo=(vzz(khit)-vzo)*weight(khit)/smass
      
c     reset system net velocity to zero
      
      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode
      do i=iatm0,iatm1
        
        vxx(i)=vxx(i)-vxo
        vyy(i)=vyy(i)-vyo
        vzz(i)=vzz(i)-vzo
        
      enddo
      
      if(mxnode.gt.1)
     x  call merge(idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)
      
      return
      end subroutine impact
      
      end module temp_scalers_module
