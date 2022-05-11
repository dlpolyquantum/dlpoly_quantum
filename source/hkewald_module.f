      module hkewald_module

c***********************************************************************
c     
c     dl_poly module for defining hautman-klein ewald sum arrays
c     copyright - daresbury laboratory
c     author    - w. smith    nov 2003
c     
c***********************************************************************

      use config_module
      use exclude_module
      use error_module
      use pair_module
      use property_module
      use setup_module

      implicit none

      real(8), allocatable :: ahk(:),crn(:,:)
      real(8), allocatable :: elc(:,:),els(:,:)
      real(8), allocatable :: emc(:,:),ems(:,:)
      real(8), allocatable :: zzn(:),zzd(:)
      real(8), allocatable :: hon(:,:),znp(:,:)
      real(8), allocatable :: dhn(:,:),zgs(:)
      real(8), allocatable :: fon(:,:),zgc(:)
      real(8), allocatable :: ckc(:),cks(:)
      real(8), allocatable :: pp(:),sss(:)

      save ahk,crn,elc,els,emc,ems,zzn,zzd,hon,znp,dhn,zgs
      save fon,zgc,ckc,cks,pp,sss

      contains
      
      subroutine alloc_hke_arrays(idnode,mxnode)
      
      implicit none
      
      integer, parameter :: nnn=9

      logical safe
      integer i,fail,idnode,mxnode
      dimension fail(nnn)

      safe=.true.

c     allocate arrays

      fail(:)=0
      
      allocate (ahk(0:mxhko),crn(0:mxhko,0:mxhko),stat=fail(1))
      allocate (elc(mxewld,0:1),els(mxewld,0:1),stat=fail(2))
      allocate (emc(mxewld,0:kmaxb),ems(mxewld,0:kmaxb),stat=fail(3))
      allocate (zzn(mxxdf),zzd(mxxdf),stat=fail(4))
      allocate (hon(mxgrid,0:mxhko),znp(mxhke,0:2*mxhko),stat=fail(5))
      allocate (dhn(mxgrid,0:mxhko),zgs(0:2*mxhko),stat=fail(6))
      allocate (fon(mxegrd,0:7),zgc(0:2*mxhko),stat=fail(7))
      allocate (ckc(mxewld),cks(mxewld),stat=fail(8))
      allocate (pp(2*mxhko),sss(mxxdf),stat=fail(9))
      
      if(any(fail.gt.0))safe=.false.      
      if(mxnode.gt.1)call gstate(safe)    
      if(.not.safe)call error(idnode,1730)
      
      end subroutine alloc_hke_arrays

      subroutine hkgen(idnode,nhko,nlatt,alpha,drewd,rcut)

c***********************************************************************
c     
c     dl_poly subroutine for generating convergence function
c     arrays for hautman klein ewald method (up to order 3 only)
c     
c     copyright - daresbury laboratory 2000
c     author    - w. smith february 2000
c     
c***********************************************************************
      
      implicit none

      integer i,idnode,nhko,nlatt,k
      real(8) alpha,drewd,rcut,ecut,den,fac,ss1,aaa,ss2

      if(nhko.gt.mxhko)call error(idnode,332)

c     define effective cutoff

      ecut=rcut*dble(2*nlatt+1)

c     define grid resolution for potential arrays
      
      drewd=ecut/dble(mxegrd-4)

c     calculate HKE coefficients

      ahk(0)=1.d0

      do i=1,nhko

        ahk(i)=-0.25d0*ahk(i-1)*dble(2*i*(2*i-1))/dble(i*i)

      enddo

c     generate convergence function arrays

      do i=1,mxegrd

        hon(i,0)=0.d0
        hon(i,1)=dble(i-1)*drewd
        hon(i,2)=(2.d0*alpha/sqrpi)*exp(-(alpha*hon(i,1))**2)

      enddo

c     generate error function and derivatives by recursion

      do k=100,1,-1

        den=1.d0/dble(2*k-1)
        fac=(2.d0*alpha**2)**(k-1)

        do i=1,mxegrd

          hon(i,0)=den*(hon(i,0)*hon(i,1)**2+fac*hon(i,2))

        enddo

        if(k.le.2*nhko+2)then

          do i=1,mxegrd

            fon(i,k-1)=hon(i,0)

          enddo

        endif

      enddo
        
c     zeroth order function
c     note: hon(1,0)=2.d0*alpha/sqrpi

      do i=1,mxegrd

        hon(i,0)= fon(i,0)
        dhn(i,0)=-fon(i,1)

      enddo

      if(nhko.eq.0)then

        ss1=dble(mxegrd-1)*drewd
        aaa=abs(1.d0-hon(mxegrd,nhko)*ss1)
        if(aaa.gt.1.d-4)then

          call warning(idnode,100,aaa,0.d0,0.d0)

        endif
        
        return

      endif

c     first order function
c     note: hon(1,1)=8.d0*alpha**3/(3.d0*sqrpi)

      do i=1,mxegrd
        
        ss2=(dble(i-1)*drewd)**2
        
        hon(i,1)=-(2.d0*fon(i,1)-fon(i,2)*ss2)
        dhn(i,1)= (4.d0*fon(i,2)-fon(i,3)*ss2)
          
      enddo
        
      if(nhko.eq.1)then

        aaa=abs(1.d0-hon(mxegrd,nhko)*sqrt(ss2)**(2*nhko+1))
        if(aaa.gt.1.d-4)then

          call warning(idnode,100,aaa,0.d0,0.d0)

        endif
        
        return

      endif

c     second order function
c     note: hon(1,2)=64.d0*alpha**5/(45.d0*sqrpi)

      do i=1,mxegrd
        
        ss2=(dble(i-1)*drewd)**2
        
        hon(i,2)=(8.d0*fon(i,2)+ss2*(-8.d0*fon(i,3)+ss2*fon(i,4)))/9.d0
        dhn(i,2)=(-24.d0*fon(i,3)+ss2*(12.d0*fon(i,4)-ss2*fon(i,5)))
     x          /9.d0

      enddo

      if(nhko.eq.2)then

        aaa=abs(1.d0-hon(mxegrd,nhko)*sqrt(ss2)**(2*nhko+1))
        if(aaa.gt.1.d-4)then

          call warning(idnode,100,aaa,0.d0,0.d0)

        endif
        
        return
        
      endif

c     third order function (enough for anyone!)
c     note: hon(1,3)=768.d0*alpha**7/(14175.d0*sqrpi)

      do i=1,mxegrd
        
        ss2=(dble(i-1)*drewd)**2
        
        hon(i,3)=-(48.d0*fon(i,3)+ss2*(-72.d0*fon(i,4)+ss2*(
     x    18.d0*fon(i,5)-ss2*fon(i,6))))/225.d0
        dhn(i,3)= (192.d0*fon(i,4)+ss2*(-144.d0*fon(i,5)+ss2*(
     x    24.d0*fon(i,6)-ss2*fon(i,7))))/225.d0
        
      enddo

      if(nhko.eq.3)then

        aaa=abs(1.d0-hon(mxegrd,nhko)*sqrt(ss2)**(2*nhko+1))
        if(aaa.gt.1.d-4)then

          call warning(idnode,100,aaa,0.d0,0.d0)

        endif
        
        return
        
      endif

      return
      end subroutine hkgen

      subroutine hkewald1
     x  (idnode,mxnode,natms,imcon,nhko,kmax1,kmax2,
     x  engcpe,vircpe,alpha,epsq)
      
c***********************************************************************
c     
c     dl_poly subroutine for calculating coulombic forces in a
c     periodic system using Hautman Klein Ewald method
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 2000
c     author    - w. smith february 2000
c     
c     part 1 - reciprocal space terms (fourier part)
c     
c     note - in loop over all k vectors k=2pi(ll/cl,mm/cl)
c     the values of ll and mm are selected so that the symmetry of
c     reciprocal lattice is taken into account i.e. the following
c     rules apply.
c     
c     ll ranges over the values 0 to kmax1 only.
c     
c     mm ranges over 1 to kmax2 when ll=0 and over
c     -kmax2 to kmax2 otherwise.
c     
c***********************************************************************

      implicit none

      logical newjob
      integer idnode,mxnode,natms,imcon,nhko,kmax1,kmax2,i
      integer iatm0,iatm1,j,k,limit,mmin,l,ll,m,mm,n
      real(8) engcpe,vircpe,alpha,epsq,twopi,ralph,area,rarea
      real(8) det,rcpcut,rcpct2,aaa,engsic,pm1,pm2,term,ssx,ssy
      real(8) tmp,rkx1,rky1,rkx2,rky2,rksq,rkk,fac,eterm,fng,fn0,gaus
      real(8) bkk,force0,forcez,pterm,scale,cprop,omg,cs
c$$$      real(8) erfc

      dimension cprop(10),omg(9)

      save newjob,engsic

      data newjob/.true./

c     initialise coulombic potential energy

      engcpe=0.d0
      vircpe=0.d0
      if(alpha.lt.1.d-8)return

c     set working parameters

      twopi=2.d0*pi
      ralph=0.5d0/alpha
      call dcell(cell,cprop)
      area=cprop(1)*cprop(2)*sqrt(1.d0-cprop(4)**2)
      rarea=pi/area

c     set up atoms numbers for nodes

      iatm0 = (idnode*natms)/mxnode + 1
      iatm1 = ((idnode+1)*natms)/mxnode

c     initalize stress tensor working arrays

      do i = 1,9

        omg(i) = 0.d0

      enddo

c     construct reciprocal lattice vectors and set k vector range

      call invert(cell,rcell,det)
      if(abs(det).lt.1.d-6)call error(idnode,120)
      call dcell(rcell,cprop)
      rcpcut=min(dble(kmax1)*cprop(7),dble(kmax2)*cprop(8))*
     x  1.05d0*twopi
      rcpct2=rcpcut**2

c     compute quantities for first entry

      if(newjob)then

        newjob=.false.

c     pbc check and array bound checks

        if(imcon.ne.6)call error(idnode,66)
        if(mxhke.ne.msatms) call error(idnode,331)
        if(mxewld.ne.msatms) call error(idnode,330)

c     check hk screening function at cutoff

        aaa=cerfr(ralph,rcpcut)
c$$$        aaa=erfc(ralph*rcpcut)/rcpcut

        if(aaa.gt.1.d-4)then

          call warning(idnode,105,aaa,0.d0,0.d0)
c          call error(idnode,487)
          
        endif

c     calculate self interaction correction
        
        engsic=0.d0
        
        do i=iatm0,iatm1
          
          engsic=engsic+chge(i)**2
          
        enddo

        engsic=-r4pie0*alpha*engsic/(sqrpi*epsq)

c     binomial coefficients

        k=0
        crn(0,0)=0.5d0
        do i=1,2*nhko

          pp(i)=1.d0
          pm1=pp(1)

          do j=2,i

            pm2=pp(j)
            pp(j)=pm2+pm1
            pm1=pm2

          enddo

          if(mod(i,2).eq.0)then

            k=k+1
            do j=0,k

              term=pp(j+1)*(-1.d0)**j
              crn(j,k)=term
              crn(k,j)=term

            enddo

            crn(k,k)=0.5d0*crn(k,k)

          endif

        enddo

      endif

c     calculate and store powers of z_i
      
      i=0

      do j=iatm0,iatm1

        i=i+1
        znp(i,0)=1.d0
        znp(i,1)=zzz(j)

      enddo

      limit=i

      do k=2,2*nhko

        do i=1,limit

          znp(i,k)=znp(i,k-1)*znp(i,1)

        enddo

      enddo

c     calculate and store exponential factors

      i=0

      do j=iatm0,iatm1
        
        i=i+1
        elc(i,0)=1.d0
        emc(i,0)=1.d0
        els(i,0)=0.d0
        ems(i,0)=0.d0
        ssx=rcell(1)*xxx(j)+rcell(4)*yyy(j)
        ssy=rcell(2)*xxx(j)+rcell(5)*yyy(j)
        elc(i,1)=cos(twopi*ssx)
        emc(i,1)=cos(twopi*ssy)
        els(i,1)=sin(twopi*ssx)
        ems(i,1)=sin(twopi*ssy)
        
      enddo

      do l=2,kmax2
        
        do i=1,limit
          
          emc(i,l)=emc(i,l-1)*emc(i,1)-ems(i,l-1)*ems(i,1)
          ems(i,l)=ems(i,l-1)*emc(i,1)+emc(i,l-1)*ems(i,1)
          
        enddo
        
      enddo

c     start of main loop over k vectors

      mmin=1
      
      do ll=0,kmax1
        
        l=ll
        tmp = twopi*dble(ll)
        rkx1=tmp*rcell(1)
        rky1=tmp*rcell(4)

c     put cos(i,L) terms into cos(i,0) array

        if(l.eq.1) then
          
          do i=1,limit

            elc(i,0)=elc(i,1)
            els(i,0)=els(i,1)
            
          enddo

        elseif(l.gt.1) then

          do i=1,limit

            cs=elc(i,0)
            elc(i,0)=cs*elc(i,1)-els(i,0)*els(i,1)
            els(i,0)=els(i,0)*elc(i,1)+cs*els(i,1)
            
          enddo
          
        endif

        do mm=mmin,kmax2
          
          m=iabs(mm)
          tmp = twopi*dble(mm)
          rkx2=rkx1+tmp*rcell(2)
          rky2=rky1+tmp*rcell(5)

c     test on magnitude of k vector
          
          rksq=rkx2*rkx2+rky2*rky2

          if(rksq.le.rcpct2)then

c     calculate exp(ikr) terms and product with charges
            
            i=0

            if(mm.ge.0)then
              
              do j=iatm0,iatm1
                
                i=i+1
                ckc(i)=chge(j)*(elc(i,0)*emc(i,m)-els(i,0)*ems(i,m))
                cks(i)=chge(j)*(els(i,0)*emc(i,m)+ems(i,m)*elc(i,0))
                
              enddo
              
            else
              
              do j=iatm0,iatm1
                
                i=i+1
                
                ckc(i)=chge(j)*(elc(i,0)*emc(i,m)+els(i,0)*ems(i,m))
                cks(i)=chge(j)*(els(i,0)*emc(i,m)-ems(i,m)*elc(i,0))
                
              enddo
              
            endif

c     calculate sum of products of powers of z_i and q_i exp(ik.s_i)
            
            do k=0,2*nhko
              
              zgc(k)=0.d0
              zgs(k)=0.d0

              do i=1,limit

                zgc(k)=zgc(k)+ckc(i)*znp(i,k)
                zgs(k)=zgs(k)+cks(i)*znp(i,k)

              enddo

            enddo

c     perform global summation of zgc and zgs arrays
            
            if(mxnode.gt.1)then
              
              call gdsum(zgc(0),2*nhko+1,buffer)
              call gdsum(zgs(0),2*nhko+1,buffer)
              
            endif

c     calculate 0th order screening function

            rkk=sqrt(rksq)
            fn0=cerfr(ralph,rkk)
c$$$            fn0=erfc(ralph*rkk)/rkk
            gaus=exp(-(ralph*rkk)**2)/(alpha*sqrpi)

c     sum terms for orders of the screening function

            fac=1.d0

            do k=0,nhko

c     sum over z_i binomial contributions

              eterm=0.d0
              fng=fac*fn0
              do m=0,k

                n=2*k-m

c     sum energy terms

                eterm=eterm+crn(m,k)*(zgc(m)*zgc(n)+zgs(m)*zgs(n))

c     calculate force contribution to each site
                
                i=0
                bkk=-fng*crn(m,k)

                do j=iatm0,iatm1
                  
                  i=i+1
                  force0=bkk*(znp(i,n)*(zgs(m)*ckc(i)-zgc(m)*cks(i))+
     x              znp(i,m)*(zgs(n)*ckc(i)-zgc(n)*cks(i)))
                  fxx(j)=fxx(j)+rkx2*force0
                  fyy(j)=fyy(j)+rky2*force0

                  omg(3)=omg(3)+rkx2*force0*zzz(j)
                  omg(6)=omg(6)+rky2*force0*zzz(j)

                  if(k.gt.0)then

                    if(m.eq.0)then
                      
                      forcez=bkk*dble(n)*znp(i,n-1)*(zgc(m)*ckc(i)+
     x                  zgs(m)*cks(i))
                      
                    else
                      
                      forcez=bkk*(dble(m)*znp(i,m-1)*(zgc(n)*ckc(i)+
     x                  zgs(n)*cks(i))+dble(n)*znp(i,n-1)*(zgc(m)*
     x                  ckc(i)+zgs(m)*cks(i)))
                    
                    endif
                    
                    omg(9)=omg(9)+forcez*zzz(j)
                    fzz(j)=fzz(j)+forcez

                  endif

                enddo

              enddo

c     accumulate potential energy and stress tensor
              
              engcpe=engcpe+fng*eterm
              pterm=(dble(2*k-1)*fng-fac*gaus)/rksq
              omg(1)=omg(1)+eterm*(fng+pterm*rkx2*rkx2)
              omg(5)=omg(5)+eterm*(fng+pterm*rky2*rky2)
              omg(2)=omg(2)+eterm*pterm*rky2*rkx2
              fac=fac*rksq/(dble(2*(k+1))*dble(2*k+1))

c     end of loop over orders of screening function

            enddo

c     end of if-block for  rksq <  rcpct2

          endif

c     end of inner loop over reciprocal lattice vectors

        enddo
        
        mmin=-kmax2

c     end of outer loop over reciprocal lattice vectors

      enddo
      
      engcpe=engcpe/dble(mxnode)
      do i = 1,9

        omg(i) = omg(i)/dble(mxnode)

      enddo

c     add self interaction correction to potential

      scale=4.d0*rarea*r4pie0/epsq
      engcpe=scale*engcpe+engsic

c     virial term

      vircpe=vircpe-scale*(omg(1)+omg(5)+omg(9))

c     calculate final forces

      do i=iatm0,iatm1

        fxx(i)=scale*fxx(i)
        fyy(i)=scale*fyy(i)
        fzz(i)=scale*fzz(i)

      enddo

c     calculate stress tensor (symmetrical)

      stress(1) = stress(1)+scale*omg(1)
      stress(2) = stress(2)+scale*omg(2)
      stress(3) = stress(3)+scale*omg(3)
      stress(4) = stress(4)+scale*omg(2)
      stress(5) = stress(5)+scale*omg(5)
      stress(6) = stress(6)+scale*omg(6)
      stress(7) = stress(7)+scale*omg(3)
      stress(8) = stress(8)+scale*omg(6)
      stress(9) = stress(9)+scale*omg(9)
      
      return
      end subroutine hkewald1

      subroutine hkewald2
     x  (idnode,mxnode,nhko,nlatt,imcon,natms,engcpe,
     x  vircpe,drewd,rcut,epsq)

c***********************************************************************
c     
c     dl_poly subroutine for calculating real-space contributions to
c     the hautman-klein-ewald electrostatic method
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 2000
c     author    - w. smith  may 2000
c     
c***********************************************************************
            
      implicit none

      integer idnode,mxnode,nhko,nlatt,imcon,natms,nix,niy
      integer nboxes,i,j,k,n,m,ma,mpm2,npm2,ii,l0,l1,l2,last
      real(8) engcpe,vircpe,drewd,rcut,epsq
      real(8) step,rcsq,rdrewd,strs1,strs2,strs3
      real(8) strs5,strs6,strs9,dcx,dcy,udx,udy,fac,chgea,chgprd
      real(8) ddx,ddy,ssx,ssy,ssq,coul,fcoul,rrr,ppp,vk0,vk1,vk2
      real(8) eterm,t1,t2,egamma,fx,fy,fz,det

      dimension nix(25),niy(25)

      data nix/ 0, 1, 1, 0,-1,-1,-1, 0, 1, 2, 2,
     x  2, 1, 0,-1,-2,-2,-2,-2,-2,-1, 0, 1, 2, 2/
      data niy/ 0, 0, 1, 1, 1, 0,-1,-1,-1, 0, 1,
     x  2, 2, 2, 2, 2, 1, 0,-1,-2,-2,-2,-2,-2,-1/

CDIR$ CACHE_ALIGN fi
      
c     check boundary condition

      if(imcon.ne.6)call error(idnode,66)

c     number of neighbouring real space cells

      if(nlatt.gt.2)call error(idnode,488)
      step=dble(2*nlatt+1)
      nboxes=(2*nlatt+1)**2

c     initialise potential energy and virial
      
      engcpe=0.d0
      vircpe=0.d0

c     set cutoff condition for pair forces
      
      rcsq=(step*rcut)**2

c     reciprocal of interpolation interval
      
      rdrewd = 1.d0/drewd

c     reciprocal cell

      call invert(cell,rcell,det)
      do i=1,9

        rcell(i)=rcell(i)/step

      enddo

c     initialise stress tensor accumulators
      strs3 = 0.d0
      strs6 = 0.d0
      strs9 = 0.d0
      strs1 = 0.d0
      strs2 = 0.d0
      strs5 = 0.d0

c     loop over image cells, starting with central cell

      ma=1
      mpm2=natms/2
      npm2=(natms-1)/2

      do k=1,nboxes

        last=natms
        dcx=dble(nix(k))
        dcy=dble(niy(k))
        udx=cell(1)*dcx+cell(4)*dcy
        udy=cell(2)*dcx+cell(5)*dcy

c     outer loop over atoms

        do m=ma,mpm2

          fac=r4pie0/epsq
          if(m.eq.0)fac=fac*0.5d0
          if(m.gt.npm2)last=mpm2

c     set initial array values

          ii=0
          do i=idnode+1,last,mxnode
            
            ii=ii+1
            chgea=fac*chge(i)
            
            if(abs(chgea).gt.1.d-10)then
              
              j=i+m
              if(j.gt.natms)j=j-natms
              
              chgprd=chgea*chge(j)
              
              if(abs(chgprd).gt.1.d-10)then

                zzn(ii)=1.d0
                zzd(ii)=0.d0

c     calculate interatomic separation
                
                ddx=xxx(i)-xxx(j)+udx
                ddy=yyy(i)-yyy(j)+udy
                ssx=rcell(1)*ddx+rcell(4)*ddy
                ssy=rcell(2)*ddx+rcell(5)*ddy
                ssx=ssx-nint(ssx)
                ssy=ssy-nint(ssy)
                xdf(ii)=step*(ssx*cell(1)+ssy*cell(4))
                ydf(ii)=step*(ssx*cell(2)+ssy*cell(5))
                zdf(ii)=zzz(i)-zzz(j)
                rsqdf(ii)=xdf(ii)**2+ydf(ii)**2+zdf(ii)**2
                
              endif

            endif

          enddo

c     loop over HK orders

          do n=0,nhko

c     inner loop over atoms

            ii=0
            do i=idnode+1,last,mxnode

              ii=ii+1
              chgea = fac*chge(i)

              if(abs(chgea).gt.1.d-10)then
                
                j=i+m
                if(j.gt.natms)j=j-natms

                chgprd=chgea*chge(j)
                
                if(abs(chgprd).gt.1.d-10)then

c     apply truncation of potential
                  
                  ssq=rsqdf(ii)-zdf(ii)*zdf(ii)

                  if(rcsq.gt.ssq)then

c     calculate potential energy and virial
                    
                    coul=0.d0
                    fcoul=0.d0
                    rrr = sqrt(rsqdf(ii))
                    sss(ii)=sqrt(ssq)

                    if(n.eq.0)then

                      coul = chgprd/rrr
                      fcoul = coul/rsqdf(ii)

                    endif

c     interpolation parameters

                    l0=int(sss(ii)*rdrewd)
                    ppp=sss(ii)*rdrewd-dble(l0)
                    l0=l0+1
                    l1=l0+1
                    l2=l0+2

c     calculate interaction energy using 3-point interpolation
                    
                    vk0 = hon(l0,n)
                    vk1 = hon(l1,n)
                    vk2 = hon(l2,n)
                    t1 = vk0 + (vk1 - vk0)*ppp
                    t2 = vk1 + (vk2 - vk1)*(ppp - 1.0d0)
                    
                    eterm=(t1+(t2-t1)*ppp*0.5d0)*ahk(n)*chgprd
                    engcpe=engcpe+coul-eterm*zzn(ii)

c     calculate forces using 3pt interpolation
                    
                    vk0 = dhn(l0,n)
                    vk1 = dhn(l1,n)
                    vk2 = dhn(l2,n)
                    
                    t1 = vk0 + (vk1 - vk0)*ppp
                    t2 = vk1 + (vk2 - vk1)*(ppp - 1.0d0)

c     calculate in-plane forces
                    
                    egamma=fcoul+
     x                (t1+(t2-t1)*ppp*0.5d0)*chgprd*zzn(ii)*ahk(n)
                    fx=egamma*xdf(ii)
                    fy=egamma*ydf(ii)

c     calculate perpendicular forces
                    
                    fz=fcoul*zdf(ii)+2.d0*dble(n)*eterm*zzd(ii)

c     add to force accumulators

                    fxx(i)=fxx(i)+fx
                    fyy(i)=fyy(i)+fy
                    fzz(i)=fzz(i)+fz

                    fxx(j)=fxx(j)-fx
                    fyy(j)=fyy(j)-fy
                    fzz(j)=fzz(j)-fz

c     reset zzn array for next order of convergence function
              
                    zzd(ii)=zzn(ii)*zdf(ii)
                    zzn(ii)=zzd(ii)*zdf(ii)
              
c     calculate stress tensor
                    
                    strs1 = strs1 + xdf(ii)*fx
                    strs2 = strs2 + xdf(ii)*fy
                    strs3 = strs3 + xdf(ii)*fz
                    strs5 = strs5 + ydf(ii)*fy
                    strs6 = strs6 + ydf(ii)*fz
                    strs9 = strs9 + zdf(ii)*fz

                  endif
                  
                endif
                
              endif

            enddo

          enddo

        enddo

        ma=0

      enddo

c     calculate virial

      vircpe=vircpe-(strs1+strs5+strs9)

c     complete stress tensor
      
      stress(1) = stress(1) + strs1
      stress(2) = stress(2) + strs2
      stress(3) = stress(3) + strs3
      stress(4) = stress(4) + strs2
      stress(5) = stress(5) + strs5
      stress(6) = stress(6) + strs6
      stress(7) = stress(7) + strs3
      stress(8) = stress(8) + strs6
      stress(9) = stress(9) + strs9
      
      return
      end subroutine hkewald2

      subroutine hkewald3(iatm,ik,engcpe,vircpe,epsq)

c***********************************************************************
c     
c     dl_poly subroutine for calculating exclusion corrections to
c     the hautman-klein-ewald electrostatic method
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 2000
c     author    - w. smith  may 2000
c     
c***********************************************************************
      
      implicit none

      integer iatm,ik,m,jatm
      real(8) engcpe,vircpe,epsq,fx,fy,fz,strs1,strs2,strs3
      real(8) strs5,strs6,strs9,chgea,chgprd,rrr,rsq,coul,fcoul

CDIR$ CACHE_ALIGN fi
      
c     initialise stress tensor accumulators

      strs1 = 0.d0
      strs2 = 0.d0
      strs3 = 0.d0
      strs5 = 0.d0
      strs6 = 0.d0
      strs9 = 0.d0

c     initialise potential energy and virial
      
      engcpe=0.d0
      vircpe=0.d0     

c     start of primary loop for forces evaluation
      
      chgea=chge(iatm)/epsq*r4pie0
      
      if(abs(chgea).gt.1.d-10)then
        
        do m=1,nexatm(ik)

c     atomic index and charge product
          
          jatm=lexatm(ik,m)
          chgprd=chgea*chge(jatm)
          
          if(abs(chgprd).gt.1.d-10)then

c     calculate interatomic distance
            
            rsq=xdf(m)**2+ydf(m)**2+zdf(m)**2
            rrr = sqrt(rsq)

c     calculate potential energy and virial
            
            coul = chgprd/rrr
            engcpe = engcpe - coul
            
c     calculate forces
            
            fcoul = coul/rsq
            fx = fcoul*xdf(m)
            fy = fcoul*ydf(m)
            fz = fcoul*zdf(m)
            
            fxx(iatm) = fxx(iatm) - fx
            fyy(iatm) = fyy(iatm) - fy
            fzz(iatm) = fzz(iatm) - fz

            fxx(jatm) = fxx(jatm) + fx
            fyy(jatm) = fyy(jatm) + fy
            fzz(jatm) = fzz(jatm) + fz

c     calculate stress tensor
            
            strs1 = strs1 - xdf(m)*fx
            strs2 = strs2 - xdf(m)*fy
            strs3 = strs3 - xdf(m)*fz
            strs5 = strs5 - ydf(m)*fy
            strs6 = strs6 - ydf(m)*fz
            strs9 = strs9 - zdf(m)*fz
            
          endif
          
        enddo

c     virial
        
        vircpe=vircpe-engcpe

c     complete stress tensor
        
        stress(1) = stress(1) + strs1
        stress(2) = stress(2) + strs2
        stress(3) = stress(3) + strs3
        stress(4) = stress(4) + strs2
        stress(5) = stress(5) + strs5
        stress(6) = stress(6) + strs6
        stress(7) = stress(7) + strs3
        stress(8) = stress(8) + strs6
        stress(9) = stress(9) + strs9

      endif
      
      return
      end subroutine hkewald3

      subroutine hkewald4(iatm,ik,engcpe,vircpe,engcpl,vircpl,rcut,epsq)

c***********************************************************************
c     
c     dl_poly subroutine for calculating coulombic forces in a
c     periodic system using the hautman-klein-ewald method
c     
c     modified to allow direct calculation of primary (short-range)
c     interactions for multiple-time step corrections
      
c     primary neighbours are taken out of the Ewald sum
c     electrostatics are evaluated directly instead
c     
c     parallel replicated data version - real space terms
c     
c     copyright - daresbury laboratory 2000
c     author    - w. smith july 2000
c     
c***********************************************************************
      
      implicit none

      integer iatm,ik,m,jatm
      real(8) engcpe,vircpe,engcpl,vircpl,rcut,epsq
      real(8) fi,fli,rcsq,strs1,strs2,strs3,strs5,strs6
      real(8) strs9,strl1,strl2,strl3,strl5,strl6,strl9,chgea,chgprd
      real(8) rsq,rrr,coul,egamma,fx,fy,fz

      dimension fi(3),fli(3)

CDIR$ CACHE_ALIGN fi
CDIR$ CACHE_ALIGN fli

c     set cutoff condition for pair forces
      
      rcsq=rcut**2

c     initialise potential energy and virial
      
      engcpe=0.d0
      vircpe=0.d0
      
      engcpl=0.d0
      vircpl=0.d0

c     initialise stress tensor accumulators

      strs1 = 0.d0
      strs2 = 0.d0
      strs3 = 0.d0
      strs5 = 0.d0
      strs6 = 0.d0
      strs9 = 0.d0
      strl1 = 0.d0
      strl2 = 0.d0
      strl3 = 0.d0
      strl5 = 0.d0
      strl6 = 0.d0
      strl9 = 0.d0

c     start of primary loop for forces evaluation
      
      chgea=chge(iatm)/epsq*r4pie0
      if(abs(chgea).gt.1.d-10)then

c     temporary arrays for cache aligning
        
        fi(1) = fxx(iatm)
        fi(2) = fyy(iatm)
        fi(3) = fzz(iatm)
        
        fli(1) = flx(iatm)
        fli(2) = fly(iatm)
        fli(3) = flz(iatm)

        do m=1,ik

c     atomic index and charge product
          
          jatm=ilist(m)
          chgprd=chgea*chge(jatm)
          
          if(abs(chgprd).gt.1.d-10)then

c     calculate interatomic distance
            
            rsq=rsqdf(m)
            
            if(rcsq.gt.rsq)then

c     coulombic energy
              
              rrr = sqrt(rsq)
              coul = chgprd/rrr

c     sum contributions to the totals 
              
              engcpe = engcpe + coul
              engcpl = engcpl - coul
              vircpe = vircpe - coul
              vircpl = vircpl + coul

c     calculate coulombic forces
              
              egamma = coul/rsq

              fx = egamma*xdf(m)
              fy = egamma*ydf(m)
              fz = egamma*zdf(m)

c     add in contributions to instantaneous force
              
              fi(1) = fi(1) + fx
              fi(2) = fi(2) + fy
              fi(3) = fi(3) + fz
              
              fxx(jatm) = fxx(jatm) - fx
              fyy(jatm) = fyy(jatm) - fy
              fzz(jatm) = fzz(jatm) - fz

c     add in contributions to the long-range force

              fli(1) = fli(1) - fx
              fli(2) = fli(2) - fy
              fli(3) = fli(3) - fz
              
              flx(jatm) = flx(jatm) + fx
              fly(jatm) = fly(jatm) + fy
              flz(jatm) = flz(jatm) + fz

c     calculate long and short range stress tensors
              
              strs1 = strs1 + xdf(m)*fx
              strl1 = strl1 - xdf(m)*fx
              strs2 = strs2 + xdf(m)*fy
              strl2 = strl2 - xdf(m)*fy
              strs3 = strs3 + xdf(m)*fz
              strl3 = strl3 - xdf(m)*fz
              strs5 = strs5 + ydf(m)*fy
              strl5 = strl5 - ydf(m)*fy
              strs6 = strs6 + ydf(m)*fz
              strl6 = strl6 - ydf(m)*fz
              strs9 = strs9 + zdf(m)*fz
              strl9 = strl9 - zdf(m)*fz

            endif

          endif
          
        enddo

c     copy back temporaries

        fxx(iatm) = fi(1)
        fyy(iatm) = fi(2)
        fzz(iatm) = fi(3)

        flx(iatm) = fli(1)
        fly(iatm) = fli(2)
        flz(iatm) = fli(3)

c     complete stress tensor
        
        stresl(1) = stresl(1) + strl1
        stresl(2) = stresl(2) + strl2
        stresl(3) = stresl(3) + strl3
        stresl(4) = stresl(4) + strl2
        stresl(5) = stresl(5) + strl5
        stresl(6) = stresl(6) + strl6
        stresl(7) = stresl(7) + strl3
        stresl(8) = stresl(8) + strl6
        stresl(9) = stresl(9) + strl9
        
        stress(1) = stress(1) + strs1
        stress(2) = stress(2) + strs2
        stress(3) = stress(3) + strs3
        stress(4) = stress(4) + strs2
        stress(5) = stress(5) + strs5
        stress(6) = stress(6) + strs6
        stress(7) = stress(7) + strs3
        stress(8) = stress(8) + strs6
        stress(9) = stress(9) + strs9

      endif
      
      return
      end subroutine hkewald4

      function cerfr(alpha,rrr)

c***********************************************************************
c     
c     dl_poly function for generating complementary error function
c     divided by r
c     
c     copyright - daresbury laboratory 2001
c     author    - w. smith february 2001
c     
c***********************************************************************

      implicit none

      integer k
      real(8) cerfr,sqrpi,h0,h1,alpha,rrr,rr2,fac

      sqrpi=1.7724538509055159d0

c     starting values

      h0=0.d0
      h1=(2.d0*alpha/sqrpi)*exp(-(alpha*rrr)**2)

c     generate function by recursion

      rr2=rrr*rrr
      do k=100,1,-1

        fac=(2.d0*alpha**2)**(k-1)
        h0=(h0*rr2+fac*h1)/dble(2*k-1)

      enddo
        
      cerfr=1.d0/rrr-h0

      return
      end function cerfr
      
      end module hkewald_module

