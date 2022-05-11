      module angles_module

c***********************************************************************
c     
c     dl_poly module for defining valence angle potentials
c     copyright - daresbury laboratory
c     author    - w. smith      sep 2003
c     modified  - p.-a.cazade   oct 2007 : solvation etc.
c     modified  - d. quigley        2010 : metadynamics
c     
c***********************************************************************

      use config_module
      use error_module
      use metafreeze_module
      use parse_module
      use property_module
      use setup_module
      use site_module
      use solvation_module
      use utility_module

      implicit none

      real(8), allocatable :: prmang(:,:)
      integer, allocatable :: listang(:,:)
      integer, allocatable :: numang(:),keyang(:),lstang(:,:)

      save prmang,listang,numang,keyang,lstang

      contains
      
      subroutine alloc_ang_arrays(idnode,mxnode)
      
c***********************************************************************
c     
c     dl_poly subroutine for defining valence angle potential arrays
c     copyright - daresbury laboratory
c     author    - w. smith    sep 2003
c     
c***********************************************************************
      
      implicit none

      logical safe
      integer i,fail,idnode,mxnode
      dimension fail(5)

      safe=.true.
      fail(:)=0

      allocate (prmang(mxtang,mxpang),stat=fail(1))
      allocate (numang(mxtmls),stat=fail(2))
      allocate (keyang(mxtang),stat=fail(3))
      allocate (lstang(mxtang,3),stat=fail(4))
      allocate (listang(mxangl,4),stat=fail(5))
      
      if(any(fail.gt.0))safe=.false.
      if(mxnode.gt.1)call gstate(safe)
      if(.not.safe)call error(idnode,1010)

c     initialise numang array
      
      do i=1,mxtmls
         numang(i)=0
      enddo
      
      end subroutine alloc_ang_arrays

      subroutine define_angles
     x  (safe,idnode,itmols,nangle,nsite,engunit)

c***********************************************************************
c     
c     dl_poly subroutine for defining bond angles
c     
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2003
c     
c***********************************************************************

      implicit none

      logical safe
      character*8 keyword
      character*1 message(80)
      integer idnode,itmols,nangle,nsite,ntmp,iang,j
      integer idum,iatm1,iatm2,iatm3,isite1,isite2,isite3,keytmp
      real(8) engunit,parpot(mxpang)

      ntmp=intstr(record,lenrec,idum)
      numang(itmols)=numang(itmols)+ntmp

      if(idnode.eq.0)then
        write(nrite,"(/,1x,'number of bond angles',
     x    10x,i10)")ntmp
        write(nrite,"(/,/,1x,'bond angle details:',
     x    /,/,18x,'unit',5x,'key',5x,'index',5x,'index',5x,
     x    'index',7x,'f-const',8x,'angle',/)")
      endif
      
      do iang=1,ntmp

        nangle=nangle+1
        if(nangle.gt.mxtang)call error(idnode,50)

c     read bond angle potential parameters
        
        call getrec(safe,idnode,nfield)
        if(.not.safe)return

        call copystring(record,message,80)
        call lowcase(record,4)
        call getword(keyword,record,4,lenrec)
        
        if(keyword(1:4).eq.'harm')then
          keytmp=1
        elseif(keyword(1:4).eq.'-hrm')then
          keytmp=-1
        elseif(keyword(1:4).eq.'quar')then
          keytmp=2
        elseif(keyword(1:4).eq.'-qur')then
          keytmp=-2
        elseif(keyword(1:4).eq.'thrm')then
          keytmp=3
        elseif(keyword(1:4).eq.'-thm')then
          keytmp=-3
        elseif(keyword(1:4).eq.'shrm')then
          keytmp=4
        elseif(keyword(1:4).eq.'-shm')then
          keytmp=-4
        elseif(keyword(1:4).eq.'bvs1')then
          keytmp=5
        elseif(keyword(1:4).eq.'-bv1')then
          keytmp=-5
        elseif(keyword(1:4).eq.'bvs2')then
          keytmp=6
        elseif(keyword(1:4).eq.'-bv2')then
          keytmp=-6
        elseif(keyword(1:4).eq.'hcos')then
          keytmp=7
        elseif(keyword(1:4).eq.'-hcs')then
          keytmp=-7
        elseif(keyword(1:4).eq.'cos ')then
          keytmp=8
        elseif(keyword(1:4).eq.'-cos')then
          keytmp=-8
        elseif(keyword(1:4).eq.'mmsb')then
          keytmp=9
        elseif(keyword(1:4).eq.'-msb')then
          keytmp=-9
        elseif(keyword(1:4).eq.'stst') then
          keytmp=10
        elseif(keyword(1:4).eq.'-sts') then
          keytmp=-10
        elseif(keyword(1:4).eq.'stbe') then
          keytmp=11
        elseif(keyword(1:4).eq.'-stb') then
          keytmp=-11
        elseif(keyword(1:4).eq.'cmps') then
          keytmp=12
        elseif(keyword(1:4).eq.'-cmp') then
          keytmp=-12
c *******************************************************************      
c     Dil Limbu
c     Method Development and Materials Simulation Laboratory
        elseif(keyword(1:4).eq.'mmbd')then
          keytmp=13
        elseif(keyword(1:4).eq.'-mbd')then
          keytmp=-13

        elseif(keyword(1:4).eq.'lcos')then
          keytmp=14
        elseif(keyword(1:4).eq.'-lcs')then
          keytmp=-14
c *******************************************************************      
        else
          if(idnode.eq.0)write(nrite,*)message
          call error(idnode,440)
        endif

        iatm1=intstr(record,lenrec,idum)
        iatm2=intstr(record,lenrec,idum)
        iatm3=intstr(record,lenrec,idum)
        parpot(1)=dblstr(record,lenrec,idum)
        parpot(2)=dblstr(record,lenrec,idum)
        parpot(3)=dblstr(record,lenrec,idum)
        parpot(4)=dblstr(record,lenrec,idum)
        parpot(5)=dblstr(record,lenrec,idum)
        parpot(6)=dblstr(record,lenrec,idum)
        
        isite1=nsite-numsit(itmols)+iatm1
        isite2=nsite-numsit(itmols)+iatm2
        isite3=nsite-numsit(itmols)+iatm3

c     test for frozen atom pairs

        if(lfzsit(isite1)*lfzsit(isite2)*
     x    lfzsit(isite3).ne.0)then
          
          numang(itmols)=numang(itmols)-1
          if(idnode.eq.0)
     x      write(nrite,"(a8,4x,i10,4x,a4,3i10,10f15.6)")
     x      '*frozen*',iang,keyword(1:4),iatm1,iatm2,iatm3,
     x      (parpot(j),j=1,mxpang)
          
        else
          
          if(idnode.eq.0)
     x      write(nrite,"(12x,i10,4x,a4,3i10,10f15.6)")
     x      iang,keyword,iatm1,iatm2,iatm3,(parpot(j),j=1,mxpang)
          
        endif
        
        keyang(nangle)=keytmp
        lstang(nangle,1)=iatm1
        lstang(nangle,2)=iatm2
        lstang(nangle,3)=iatm3
        prmang(nangle,1)=parpot(1)
        prmang(nangle,2)=parpot(2)
        prmang(nangle,3)=parpot(3)
        prmang(nangle,4)=parpot(4)
        prmang(nangle,5)=parpot(5)
        prmang(nangle,6)=parpot(6)
        
c     convert energies to internal units
        
        prmang(nangle,1)=prmang(nangle,1)*engunit
        if(abs(keytmp).eq.2)then
          prmang(nangle,3)=prmang(nangle,3)*engunit
          prmang(nangle,4)=prmang(nangle,4)*engunit
        elseif(abs(keytmp).eq.12)then
          prmang(nangle,2)=prmang(nangle,2)*engunit            
          prmang(nangle,3)=prmang(nangle,3)*engunit
        endif
        
c     convert angles to radians
        
        if(abs(keytmp).eq.12)then
          prmang(nangle,4)=prmang(nangle,4)*(pi/180.d0)
        elseif((abs(keytmp).ne.10).and.(abs(keytmp).ne.14))then
          prmang(nangle,2)=prmang(nangle,2)*(pi/180.d0) 
        endif
        
      enddo
      
      return
      end subroutine define_angles

      subroutine angfrc
     x  (lsolva,lfree,lexcite,idnode,imcon,mxnode,ntangl,engang,virang)

c***********************************************************************
c     
c     dl_poly subroutine for calculating bond angle energy and 
c     force terms in molecular dynamics.
c     
c     replicated data - blocked version
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith         may 1992
c     modified  - t. forester      feb 1993
c     modified  - t. forester      nov 1994 : block data
c     modified  - t. forester      may 1995 : stress tensor 
c     modified  - p.-a.cazade      oct 2007 : solvation etc.
c     modified  - d. quigley       nov 2010 : metadynamics
c     
c***********************************************************************
      
      implicit none
      
      logical safe,lsolva,lfree,lexcite,lselect
      logical idrive,jdrive,kdrive
      integer idnode,mxnode,imcon,ntangl,fail1,fail2
      integer ii,iang1,iang2,i,ia,ib,ic,kk,keya
      real(8) engang,virang,theta,fxc,fyc,fzc,rab,xab
      real(8) yab,zab,rbc,xbc,ybc,zbc,sint,cost,pterm,vterm
      real(8) gamma,gamsa,gamsc,rrbc,rrab,fxa,fya,fza
      real(8) strs(6),strs_loc(6)
      real(8), allocatable :: xdab(:),ydab(:),zdab(:)
      real(8), allocatable :: xdbc(:),ydbc(:),zdbc(:)

      real(8)  :: dtheta,xdeg
      
      safe=.true.
      
      allocate (xdab(msbad),ydab(msbad),zdab(msbad),stat=fail1)
      allocate (xdbc(msbad),ydbc(msbad),zdbc(msbad),stat=fail2)
      
      if(fail1.gt.0.or.fail2.gt.0)safe=.false.
      if(mxnode.gt.1)call gstate(safe)    
      if(.not.safe)call error(idnode,1020)
      
c     check size of work arrays

      if((ntangl-mxnode+1)/mxnode.gt.msbad)call error(idnode,419)

c     block indices
      
      iang1=(idnode*ntangl)/mxnode+1
      iang2=((idnode+1)*ntangl)/mxnode

c     zero accumulators
      
      engang=0.d0
      virang=0.d0
      ang_fre=0.d0
      ang_vir=0.d0
      strs(:)=0.d0
      strs_loc(:)=0.d0
      
      if(lsolva)then
        
        lcomp(2)=.true.
        ang_sol(:)=0.d0
        if(lexcite)ang_exc(:)=0.d0
        
      endif
      
c     calculate atom separation vectors
      
      ii=0
      do i=iang1,iang2
        
        ii=ii+1

c     indices of bonded atoms
        
        ia=listang(ii,2)
        ib=listang(ii,3)
        ic=listang(ii,4)

c     components of first bond vector
        
        xdab(ii)=xxx(ia)-xxx(ib)
        ydab(ii)=yyy(ia)-yyy(ib)
        zdab(ii)=zzz(ia)-zzz(ib)

c     components of second bond vector
        
        xdbc(ii)=xxx(ic)-xxx(ib)
        ydbc(ii)=yyy(ic)-yyy(ib)
        zdbc(ii)=zzz(ic)-zzz(ib)
        
      enddo

c     periodic boundary condition
      
      call images(imcon,0,1,ii,cell,xdab,ydab,zdab)
      call images(imcon,0,1,ii,cell,xdbc,ydbc,zdbc)

c     loop over all specified angle potentials
      
      ii=0
      do i=iang1,iang2
        
        ii=ii+1

c     define components of first bond vector
        
        rab=sqrt(xdab(ii)**2+ydab(ii)**2+zdab(ii)**2)
        rrab=1.d0/rab
        
        xab=xdab(ii)*rrab
        yab=ydab(ii)*rrab
        zab=zdab(ii)*rrab

c     define components of second bond vector
        
        rbc=sqrt(xdbc(ii)**2+ydbc(ii)**2+zdbc(ii)**2)
        rrbc=1.d0/rbc
        
        xbc=xdbc(ii)*rrbc
        ybc=ydbc(ii)*rrbc
        zbc=zdbc(ii)*rrbc

c     index of potential function parameters
        
        kk=listang(ii,1)

c     determine bond angle and calculate potential energy
        
        cost=(xab*xbc+yab*ybc+zab*zbc)
        theta=acos(cost)
        sint=max(1.d-8,sqrt(1.d0-cost**2))
        
        keya=abs(keyang(kk))
        
        if(keya.eq.1)then
          
c     harmonic potential
          
          pterm=0.5d0*prmang(kk,1)*(theta-prmang(kk,2))**2
          gamma=prmang(kk,1)*(theta-prmang(kk,2))/sint
          vterm=0.d0
          gamsa=0.d0
          gamsc=0.d0
          
        elseif(keya.eq.2)then

c     quartic potential
          
          pterm=0.5d0*prmang(kk,1)*(theta-prmang(kk,2))**2+
     x      1.d0/3.d0*prmang(kk,3)*(theta-prmang(kk,2))**3+
     x      0.25d0*prmang(kk,4)*(theta-prmang(kk,2))**4
          gamma=(prmang(kk,1)*(theta-prmang(kk,2))+
     x      prmang(kk,3)*(theta-prmang(kk,2))**2+
     x      prmang(kk,4)*(theta-prmang(kk,2))**3)/sint
          vterm=0.d0
          gamsa=0.d0
          gamsc=0.d0
          
        elseif(keya.eq.3)then

c     truncated harmonic potential
      
          pterm=0.5d0*prmang(kk,1)*(theta-prmang(kk,2))**2*
     x      exp(-(rab**8+rbc**8)/prmang(kk,3)**8)
          gamma=prmang(kk,1)*(theta-prmang(kk,2))*
     x      exp(-(rab**8+rbc**8)/prmang(kk,3)**8)/sint
          vterm=-8.d0*pterm*(rab**8+rbc**8)/prmang(kk,3)**8
          gamsa=(8.d0*pterm/prmang(kk,3)**8)*rab**7
          gamsc=(8.d0*pterm/prmang(kk,3)**8)*rbc**7
          
        elseif(keya.eq.4)then

c     screened harmonic potential
          
          pterm=0.5d0*prmang(kk,1)*(theta-prmang(kk,2))**2*
     x      exp(-(rab/prmang(kk,3)+rbc/prmang(kk,4)))
          gamma=prmang(kk,1)*(theta-prmang(kk,2))*
     x      exp(-(rab/prmang(kk,3)+rbc/prmang(kk,4)))/sint
          vterm=-pterm*(rab/prmang(kk,3)+rbc/prmang(kk,4))
          gamsa=(pterm/prmang(kk,3))
          gamsc=(pterm/prmang(kk,4))
          
        elseif(keya.eq.5)then

c     screened vessal potential (type 1)
      
          pterm=(prmang(kk,1)/(8.d0*(prmang(kk,2)-pi)**2)*
     x      (((prmang(kk,2)-pi)**2-(theta-pi)**2)**2))*
     x      exp(-(rab/prmang(kk,3)+rbc/prmang(kk,4)))
          gamma=(prmang(kk,1)/(2.d0*(prmang(kk,2)-pi)**2)*
     x      ((prmang(kk,2)-pi)**2-(theta-pi)**2)*(theta-pi))*
     x      exp(-(rab/prmang(kk,3)+rbc/prmang(kk,4)))/sint
          vterm=-pterm*(rab/prmang(kk,3)+rbc/prmang(kk,4))
          gamsa=(pterm/prmang(kk,3))
          gamsc=(pterm/prmang(kk,4))
          
        elseif(keya.eq.6)then

c     truncated vessal potential (type 2)
          
          pterm=prmang(kk,1)*(theta**prmang(kk,3)*(theta-prmang(kk,2))
     x      **2*(theta+prmang(kk,2)-2.d0*pi)**2-0.5d0*prmang(kk,3)*pi
     x      **(prmang(kk,3)-1.d0)*(theta-prmang(kk,2))**2*(pi-prmang(kk,
     x      2))**3)*exp(-(rab**8+rbc**8)/prmang(kk,4)**8)
          gamma=prmang(kk,1)*(theta**(prmang(kk,3)-1.d0)*(theta-prmang
     x      (kk,2))*(theta+prmang(kk,2)-2.d0*pi)*((prmang(kk,3)+4.d0)*
     x      theta**2-2.d0*pi*(prmang(kk,3)+2.d0)*theta+prmang(kk,3)*
     x      prmang(kk,2)*(2.d0*pi-prmang(kk,2)))-prmang(kk,3)*pi**
     x      (prmang(kk,3)-1.d0)*(theta-prmang(kk,2))*(pi-prmang(kk,2))
     x      **3)*exp(-(rab**8+rbc**8)/prmang(kk,4)**8)/sint
          vterm=-8.d0*pterm*(rab**8+rbc**8)/prmang(kk,4)**8
          gamsa=(8.d0*pterm/prmang(kk,4)**8)*rab**7
          gamsc=(8.d0*pterm/prmang(kk,4)**8)*rbc**7

        elseif(keya.eq.7)then
          
c     harmonic cosine potential
          
          pterm=0.5d0*prmang(kk,1)*(cos(theta)-cos(prmang(kk,2)))**2
          gamma=-prmang(kk,1)*(cos(theta)-cos(prmang(kk,2)))
          vterm=0.d0
          gamsa=0.d0
          gamsc=0.d0
          
        elseif(keya.eq.8)then
          
c     ordinary cosine potential
          
          pterm=prmang(kk,1)*(1+cos(prmang(kk,3)*theta-prmang(kk,2)))
          gamma=-prmang(kk,1)*prmang(kk,3)*sin(prmang(kk,3)*theta-
     x      prmang(kk,2))/sint
          vterm=0.d0
          gamsa=0.d0
          gamsc=0.d0
          
        elseif(keya.eq.9)then

c     mm3 stretch-bend potential

          pterm=prmang(kk,1)*(theta-prmang(kk,2))*
     x      (rab-prmang(kk,3))*(rbc-prmang(kk,4))
          gamma=prmang(kk,1)*(rab-prmang(kk,3))*(rbc-
     x      prmang(kk,4))/sint
          gamsa=-prmang(kk,1)*(theta-prmang(kk,2))*(rbc-prmang(kk,4))
          gamsc=-prmang(kk,1)*(theta-prmang(kk,2))*(rab-prmang(kk,3))
          vterm=-(gamsa*rab+gamsc*rbc)

        elseif(keya.eq.10)then

c     compass stretch-stretch potential

          pterm=prmang(kk,1)*(rab-prmang(kk,2))*(rbc-prmang(kk,3))
          gamma=0.d0
          gamsa=-prmang(kk,1)*(rbc-prmang(kk,3))
          gamsc=-prmang(kk,1)*(rab-prmang(kk,2))
          vterm=-(gamsa*rab+gamsc*rbc)

        elseif(keya.eq.11)then

c     compass stretch-bend potential

          pterm=prmang(kk,1)*(theta-prmang(kk,2))*(rab-prmang(kk,3))
          gamma=prmang(kk,1)*(rab-prmang(kk,3))/sint
          gamsa=-prmang(kk,1)*(theta-prmang(kk,2))
          gamsc=0.d0
          vterm=-gamsa*rab

        elseif(keya.eq.12)then

c     combined compass angle potential with 3 coupling terms

          pterm=prmang(kk,1)*(rab-prmang(kk,5))*(rbc-prmang(kk,6))+
     x          (theta-prmang(kk,4))*(prmang(kk,2)*(rab-prmang(kk,5))+
     x          prmang(kk,3)*(rbc-prmang(kk,6)))
          gamma=(prmang(kk,2)*(rab-prmang(kk,5))+
     x           prmang(kk,3)*(rbc-prmang(kk,6)))/sint
          gamsa=-prmang(kk,2)*(theta-prmang(kk,4))-
     x           prmang(kk,1)*(rbc-prmang(kk,6))
          gamsc=-prmang(kk,3)*(theta-prmang(kk,4))-
     x           prmang(kk,1)*(rab-prmang(kk,5))
          vterm=-(gamsa*rab+gamsc*rbc)

c *******************************************************************     
c     Dil Limbu 
c     Method Development and Materials Simulation Laboratory

        else if(keya.eq.13)then
c       enegry unit in mdyn and angle in degree
         xdeg = 180.d0/pi 
c     MM3 angle bend potential
         dtheta = (theta-prmang(kk,2))

        pterm=71.94d0*prmang(kk,1)*dtheta**2*(1.d0-0.014d0*(xdeg*dtheta)
     x         +5.6e-5*(xdeg*dtheta)**2-7.0e-7*(xdeg*dtheta)**3
     x         +2.2e-8*(xdeg*theta)**4)

         gamma=71.94d0*prmang(kk,1)*dtheta*(2.d0-0.042d0*(xdeg*dtheta)
     x         +2.24e-4*(xdeg*dtheta)**2-3.5e-6*(xdeg*dtheta)**3
     x         +1.32e-7*(xdeg*theta)**4)/sint

          vterm=0.d0
          gamsa=0.d0
          gamsc=0.d0


        else if(keya.eq.14)then
      
c     cosine periodic (S13-Smit2017_SI) angle potential

         pterm=prmang(kk,1)*(1.d0-prmang(kk,2)*(-1)**prmang(kk,3)
     x         *cos(prmang(kk,3)*theta))

        gamma=prmang(kk,1)*prmang(kk,3)*prmang(kk,2)*(-1)**prmang(kk,3)
     x         *sin(prmang(kk,3)*theta)/sint

c          pterm=prmang(kk,1)*(1.d0-cos(4.d0*theta))

c          gamma=prmang(kk,1)*4.d0*sin(4.d0*theta)/sint

          vterm=0.d0
          gamsa=0.d0
          gamsc=0.d0
c *******************************************************************     

        else

c     undefined potential

          safe=.false.
          pterm=0.d0
          vterm=0.d0
          gamma=0.d0
          gamsa=0.d0
          gamsc=0.d0
          
        endif
        
c     indices of bonded atoms
        
        ia=listang(ii,2)
        ib=listang(ii,3)
        ic=listang(ii,4)

c     metadynamics local definitions
        
        if(lmetadyn)then
          
          idrive=driven(ltype(ia))
          jdrive=driven(ltype(ib))
          kdrive=driven(ltype(ic))
          
        endif

c     set selection control
        
        lselect=.true.
        
        if(lexcite)then
          
c     selected excitation option
        
          if((atm_fre(ia).ne.1).and.(atm_fre(ib).ne.1).and.
     x      (atm_fre(ic).ne.1))then
            
c     reset selection control
            
            lselect=(atm_fre(ia)+atm_fre(ib)+atm_fre(ic).eq.0)
            
            if(lsolva)then
              ang_exc(atmolt(ia))=ang_exc(atmolt(ia))+pterm
            endif
            
          endif
          
        elseif(lfree)then
          
c     selected free energy option
          
          if((atm_fre(ia).eq.1).or.(atm_fre(ib).eq.1).or.
     x        (atm_fre(ic).eq.1))then
            
c     set hamiltonian mixing parameter

            ang_fre=ang_fre-pterm
            ang_vir=ang_vir-vterm
            pterm=lambda1*pterm
            vterm=lambda1*vterm
            gamma=lambda1*gamma
            gamsa=lambda1*gamsa
            gamsc=lambda1*gamsc
            
          elseif((atm_fre(ia).eq.2).or.(atm_fre(ib).eq.2).or.
     x        (atm_fre(ic).eq.2))then
            
c     set hamiltonian mixing parameter

            ang_fre=ang_fre+pterm
            ang_vir=ang_vir+vterm
            pterm=lambda2*pterm
            vterm=lambda2*vterm
            gamma=lambda2*gamma
            gamsa=lambda2*gamsa
            gamsc=lambda2*gamsc
                        
          endif
          
        endif
        
        if(lselect)then
          
c     sum potential energy and virial
          
          engang=engang+pterm
          virang=virang+vterm
          
c     calculate solvation energy
        
          if(lsolva)then
            ang_sol(atmolt(ia))=ang_sol(atmolt(ia))+pterm
          endif
          
c     calculate atomic forces
        
          fxa=gamma*(xbc-xab*cost)*rrab+gamsa*xab
          fya=gamma*(ybc-yab*cost)*rrab+gamsa*yab
          fza=gamma*(zbc-zab*cost)*rrab+gamsa*zab
          
          fxc=gamma*(xab-xbc*cost)*rrbc+gamsc*xbc
          fyc=gamma*(yab-ybc*cost)*rrbc+gamsc*ybc
          fzc=gamma*(zab-zbc*cost)*rrbc+gamsc*zbc
          
c     sum atomic forces

          fxx(ia)=fxx(ia)+fxa
          fyy(ia)=fyy(ia)+fya
          fzz(ia)=fzz(ia)+fza
          
          fxx(ib)=fxx(ib)-fxa-fxc
          fyy(ib)=fyy(ib)-fya-fyc
          fzz(ib)=fzz(ib)-fza-fzc
          
          fxx(ic)=fxx(ic)+fxc
          fyy(ic)=fyy(ic)+fyc
          fzz(ic)=fzz(ic)+fzc
          
c     calculate stress tensor
          
          strs(1)=strs(1)+rab*xab*fxa+rbc*xbc*fxc
          strs(2)=strs(2)+rab*xab*fya+rbc*xbc*fyc
          strs(3)=strs(3)+rab*xab*fza+rbc*xbc*fzc
          strs(4)=strs(4)+rab*yab*fya+rbc*ybc*fyc
          strs(5)=strs(5)+rab*yab*fza+rbc*ybc*fzc
          strs(6)=strs(6)+rab*zab*fza+rbc*zbc*fzc
          
        endif
        
c     metadynamics local parameters
        
        if(lmetadyn.and.(idrive.or.jdrive.or.kdrive))then
          
c     local energy and virial
          
          eng_loc=eng_loc+pterm
          vir_loc=vir_loc+vterm
          
c     local forces
          
          fxx_loc(ia)=fxx_loc(ia)+fxa
          fyy_loc(ia)=fyy_loc(ia)+fya
          fzz_loc(ia)=fzz_loc(ia)+fza
          
          fxx_loc(ib)=fxx_loc(ib)-fxa-fxc
          fyy_loc(ib)=fyy_loc(ib)-fya-fyc
          fzz_loc(ib)=fzz_loc(ib)-fza-fzc
          
          fxx_loc(ic)=fxx_loc(ic)+fxc
          fyy_loc(ic)=fyy_loc(ic)+fyc
          fzz_loc(ic)=fzz_loc(ic)+fzc
          
c     local stress tensor
          
          strs_loc(1)=strs_loc(1)+rab*xab*fxa+rbc*xbc*fxc
          strs_loc(2)=strs_loc(2)+rab*xab*fya+rbc*xbc*fyc
          strs_loc(3)=strs_loc(3)+rab*xab*fza+rbc*xbc*fzc
          strs_loc(4)=strs_loc(4)+rab*yab*fya+rbc*ybc*fyc
          strs_loc(5)=strs_loc(5)+rab*yab*fza+rbc*ybc*fzc
          strs_loc(6)=strs_loc(6)+rab*zab*fza+rbc*zbc*fzc

        endif
        
      enddo

c     check for undefined potentials

      if(mxnode.gt.1)call gstate(safe)
      if(.not.safe)call error(idnode,440)

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
      
c     sum up contributions to potential and virial
      
      if(mxnode.gt.1)then

        buffer(1)=engang
        buffer(2)=virang
        buffer(3)=ang_fre
        buffer(4)=ang_vir
        call gdsum(buffer(1),4,buffer(5))
        engang=buffer(1)
        virang=buffer(2)
        ang_fre=buffer(3)
        ang_vir=buffer(4)

c     sum up solvation energies
        
        if(lsolva)then

          call gdsum(ang_sol(1),mxtmls,buffer(1))
          if(lexcite)call gdsum(ang_exc(1),mxtmls,buffer(1))
          
        endif
        
      endif
      
      deallocate (xdab,ydab,zdab,stat=fail1)
      deallocate (xdbc,ydbc,zdbc,stat=fail2)
      
      return
      end subroutine angfrc
      
      end module angles_module
