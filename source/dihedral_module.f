      module dihedral_module
      
c***********************************************************************
c     
c     dl_poly module for defining dihedral potential arrays
c     copyright - daresbury laboratory
c     
c     author  - w. smith     sep 2003
c     adapted - p.-a. cazade oct 2007 : solvation, free energy, etc.
c     adapted - w.smith      jan 2011 : metadynamics
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
      use vdw_module
      
      implicit none
      
      real(8), allocatable :: prmdih(:,:)
      integer, allocatable :: listdih(:,:)
      integer, allocatable :: numdih(:),keydih(:),lstdih(:,:)
      
      save prmdih,listdih,numdih,keydih,lstdih
      
      contains
      
      subroutine alloc_dih_arrays(idnode,mxnode)
      
      implicit none

      logical safe
      integer i,fail,idnode,mxnode
      dimension fail(5)
      
      safe=.true.
      
c     allocate arrays
      
      fail(:)=0
      
      allocate (prmdih(mxtdih,mxpdih),stat=fail(1))
      allocate (numdih(mxtmls),stat=fail(2))
      allocate (keydih(mxtdih),stat=fail(3))
      allocate (lstdih(mxtdih,4),stat=fail(4))
      allocate (listdih(mxdihd,5),stat=fail(5))
      
      if(any(fail.gt.0))safe=.false.      
      if(mxnode.gt.1)call gstate(safe)    
      if(.not.safe)call error(idnode,1011)

c     initialise numdih array
      
      do i=1,mxtmls
        numdih(i)=0
      enddo
      
      end subroutine alloc_dih_arrays
      
      subroutine define_dihedrals
     x  (safe,idnode,itmols,ndihed,nsite,engunit)
      
c***********************************************************************
c     
c     dl_poly subroutine for defining dihedral angles
c     
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2003
c     
c***********************************************************************
      
      implicit none
      
      logical safe
      character*8 keyword
      character*1 message(80)
      integer idnode,itmols,ndihed,nsite,ntmp,idih,j,keytmp
      integer iatm1,iatm2,iatm3,iatm4,idum,isite1,isite2,isite3
      integer isite4
      real(8) engunit,parpot(mxpdih)
      
      ntmp=intstr(record,lenrec,idum)
      numdih(itmols)=numdih(itmols)+ntmp

      if(idnode.eq.0)then
        write(nrite,"(/,1x,'number of dihedral angles',
     x    6x,i10)")numdih(itmols)
        write(nrite,"(/,1x,'dihedral angle details:',
     x    /,/,12x,'unit',5x,'key',5x,'index',5x,'index',5x,'index',
     x    5x,'index',7x,'f-const',8x,'angle',9x,'trig',11x,
     x    '1-4 elec',7x,'1-4 vdw',/)")
      endif
      
      do idih=1,ntmp
        
        ndihed=ndihed+1
        if(ndihed.gt.mxtdih)call error(idnode,60)
        
c     read dihedral bond angle potential parameters
        
        call getrec(safe,idnode,nfield)
        if(.not.safe)return
        
        call copystring(record,message,80)
        call lowcase(record,4)
        call getword(keyword,record,4,lenrec)

        if(keyword(1:4).eq.'cos ')then
          keytmp=1
        elseif(keyword(1:4).eq.'harm')then
          keytmp=2
        elseif(keyword(1:4).eq.'hcos')then
          keytmp=3
        elseif(keyword(1:4).eq.'cos3')then
          keytmp=4
        elseif(keyword(1:4).eq.'ryck')then
          keytmp=5
        elseif(keyword(1:4).eq.'rbf')then 
          keytmp=6
        elseif(keyword(1:4).eq.'opls')then 
          keytmp=7
        else
          if(idnode.eq.0)write(nrite,*)message
          call error(idnode,448)
        endif

        iatm1=intstr(record,lenrec,idum)
        iatm2=intstr(record,lenrec,idum)
        iatm3=intstr(record,lenrec,idum)
        iatm4=intstr(record,lenrec,idum)
        parpot(1)=dblstr(record,lenrec,idum)
        parpot(2)=dblstr(record,lenrec,idum)
        parpot(3)=dblstr(record,lenrec,idum)
        parpot(4)=dblstr(record,lenrec,idum)
        parpot(5)=dblstr(record,lenrec,idum)
        
        isite1=nsite-numsit(itmols)+iatm1
        isite2=nsite-numsit(itmols)+iatm2
        isite3=nsite-numsit(itmols)+iatm3
        isite4=nsite-numsit(itmols)+iatm4
          
c     test for frozen atom pairs
          
        if(lfzsit(isite1)*lfzsit(isite2)*
     x    lfzsit(isite3)*lfzsit(isite4).ne.0)then
          
          numdih(itmols)=numdih(itmols)-1
          if(idnode.eq.0)
     x      write(nrite,"(4x,a8,i10,4x,a4,4i10,10f15.6)")
     x      '*frozen*',idih,keyword(1:4),iatm1,iatm2,iatm3,iatm4,
     x      (parpot(j),j=1,mxpdih)
          
        else
          
          if(idnode.eq.0)
     x      write(nrite,"(12x,i10,4x,a4,4i10,10f15.6)")
     x      idih,keyword(1:4),iatm1,iatm2,iatm3,iatm4,
     x      (parpot(j),j=1,mxpdih)
          
        endif
        
        keydih(ndihed)=keytmp
        lstdih(ndihed,1)=iatm1
        lstdih(ndihed,2)=iatm2
        lstdih(ndihed,3)=iatm3
        lstdih(ndihed,4)=iatm4
        prmdih(ndihed,1)=parpot(1)
        prmdih(ndihed,2)=parpot(2)
        prmdih(ndihed,3)=parpot(3)
        prmdih(ndihed,4)=parpot(4)
        prmdih(ndihed,5)=parpot(5)
          
c     convert energies to internal units and angles to radians
          
        prmdih(ndihed,1)=prmdih(ndihed,1)*engunit
        
        if(keytmp.eq.4)then
          
          prmdih(ndihed,2)=prmdih(ndihed,2)*engunit
          prmdih(ndihed,3)=prmdih(ndihed,3)*engunit
          
        elseif(keytmp.eq.7)then
          
          prmdih(ndihed,2)=prmdih(ndihed,2)*engunit
          prmdih(ndihed,3)=prmdih(ndihed,3)*engunit
          prmdih(ndihed,4)=prmdih(ndihed,4)*engunit
          prmdih(ndihed,5)=prmdih(ndihed,5)*(pi/180.d0)
          
        else
          
          prmdih(ndihed,2)=prmdih(ndihed,2)*(pi/180.d0)
          
        endif
        
      enddo
      
      return
      end subroutine define_dihedrals
      
      subroutine dihfrc
     x  (lsolva,lfree,lexcite,idnode,imcon,mxnode,ntdihd,keyfce,
     x  dlrpot,epsq,engcpe,engdih,engsrp,rcut,rvdw,alpha,vircpe,
     x  virdih,virsrp)
      
c***********************************************************************
c     
c     dl_poly subroutine for calculating dihedral energy and force 
c     terms in molecular dynamics.
c     
c     version 3: scale factors for reduces electrostatic and vdw
c     1-4 interactions.
c     
c     NOTE: assumes 1-4 interactions are in the exclude list
c     
c     block as opposed to stride version
c     
c     copyright - daresbury laboratory
c     author    - w. smith       mar 1992
c     modified  - t. forester    dec 1993
c     modified  - t. forester    jun 1995 - stress tensor added
c     modified  - a. smondyrev   may 2000 - ryckaert-bellemans potentials
c     modified  - p.-a. cazade   oct 2007 - solvation etc
c     
c***********************************************************************
      
      implicit none
      
      logical safe,lsolva,lfree,lexcite,lselect
      logical idrive,jdrive,kdrive,ldrive
      integer i,k,ii,kk,ntdihd,idnode,mxnode,idih1,idih2,ia,ib
      integer ic,id,imcon,ka,kb,l,keyfce,fail1,fail2,fail3,kkk
      real(8) phi,twopi,rtwopi,dterm,srpot
      real(8) engdih,virdih,engc14,engs14,virs14,rrbc,xab,yab,erc,fer
      real(8) zab,xbc,ybc,zbc,xcd,ycd,zcd,pbx,pby,pbz,pb2,rpb1,rpb2
      real(8) pcx,pcy,pcz,pc2,rpc1,rpc2,pbpc,cosp,sinp,rsinp,exp1
      real(8) gamma,fax,fay,faz,fcx,fcy,fcz,fb1x,fb1y,fb1z,fd1x,fd1y
      real(8) fd1z,scale,xad,yad,zad,rad,chgprd,coul,fcoul,fx,fy,fz
      real(8) ppp,dlrpot,t1,t2,vk0,vk1,vk2,gk0,gk1,gk2,epsq,engcpe
      real(8) vircpe,rcut,rvdw,engsrp,virsrp,xac,yac,zac,vcon,fcon
      real(8) virc14,b0,rfld0,rfld1,rfld2,alpha,a1,a2,a3,a4,a5,pp,tt
      real(8) cou14_vir,vdw14_vir,strs(6),strs_loc(6)
      real(8), allocatable :: xdab(:),ydab(:),zdab(:)
      real(8), allocatable :: xdbc(:),ydbc(:),zdbc(:)
      real(8), allocatable :: xdcd(:),ydcd(:),zdcd(:)
      data a1,a2,a3/0.254829592d0,-0.284496736d0,1.421413741d0/
      data a4,a5,pp/-1.453152027d0,1.061405429d0,0.3275911d0/
      data fail1,fail2,fail3/0,0,0/
      
      allocate (xdab(msbad),ydab(msbad),zdab(msbad),stat=fail1)
      allocate (xdbc(msbad),ydbc(msbad),zdbc(msbad),stat=fail2)
      allocate (xdcd(msbad),ydcd(msbad),zdcd(msbad),stat=fail3)
      if(fail1.ne.0.or.fail2.ne.0.or.fail3.ne.0)
     x  call error(idnode,1060)
      
      twopi=2.d0*pi
      rtwopi=1.d0/twopi
      safe=.true.
      
c     check size of work arrays
      
      if((ntdihd-mxnode+1)/mxnode.gt.msbad) call error(idnode,421)
      
c     block indices
      
      idih1=(idnode*ntdihd)/mxnode+1
      idih2=((idnode+1)*ntdihd)/mxnode
      
c     initialise accumulators
      
      engdih=0.d0
      virdih=0.d0
      dih_fre=0.d0
      dih_vir=0.d0
      strs(:)=0.d0
      strs_loc(:)=0.d0

      if(lsolva)then
        
        lcomp(3)=.true.
        dih_sol(:)=0.d0
        if(lexcite)dih_exc(:)=0.d0
        
      endif
      
      if(keyfce/2.eq.4)then
        
c     constant terms for shifted coulombic potential
        
        tt=1.d0/(1.d0+pp*alpha*rcut)
        exp1=exp(-(alpha*rcut)**2)
        vcon=tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1/rcut
        fcon=(vcon+2.d0*(alpha/sqrpi)*exp1)/rcut
        
      elseif(keyfce/2.eq.5)then
        
c     constant terms for reaction field potential
        
        b0=2.d0*(epsq-1.d0)/(2.d0*epsq+1.d0)
        rfld0=b0/rcut**3
        rfld1=(1.d0+b0*0.5d0)/rcut
        rfld2=rfld0*0.5d0
        tt=1.d0/(1.d0+pp*alpha*rcut)
        exp1=exp(-(alpha*rcut)**2)
        vcon=tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1/rcut
        fcon=(vcon+2.d0*(alpha/sqrpi)*exp1)/rcut-rfld0*rcut
        vcon=vcon+rfld2*rcut**2-rfld1
        
      endif
      
c     calculate bond vectors
      
      ii=0
      do i=idih1,idih2
        
        ii=ii+1
        
c     indices of bonded atoms
        
        ia=listdih(ii,2)
        ib=listdih(ii,3)
        ic=listdih(ii,4)
        id=listdih(ii,5)
        
c     define components of bond vectors
        
        xdab(ii)=xxx(ia)-xxx(ib)
        ydab(ii)=yyy(ia)-yyy(ib)
        zdab(ii)=zzz(ia)-zzz(ib)
        
        xdbc(ii)=xxx(ib)-xxx(ic)
        ydbc(ii)=yyy(ib)-yyy(ic)
        zdbc(ii)=zzz(ib)-zzz(ic)
        
        xdcd(ii)=xxx(ic)-xxx(id)
        ydcd(ii)=yyy(ic)-yyy(id)
        zdcd(ii)=zzz(ic)-zzz(id)
        
      enddo
      
c     periodic boundary condition
      
      call images(imcon,0,1,ii,cell,xdab,ydab,zdab)
      call images(imcon,0,1,ii,cell,xdbc,ydbc,zdbc)
      call images(imcon,0,1,ii,cell,xdcd,ydcd,zdcd)
      
c     zero dihedral energy accumulator
      
      engdih=0.d0
      virdih=0.d0
      
c     zero scaled 1-4 electrostatic and short range potential accumulators
      
      engc14=0.d0
      virc14=0.d0
      engs14=0.d0
      virs14=0.d0
      vdw14_vir=0.d0
      cou14_vir=0.d0
      
c     loop over all specified dihedrals
      
      ii=0
      do i=idih1,idih2
        
c     define components of bond vectors
        
        ii=ii+1
        
        xab=xdab(ii)
        yab=ydab(ii)
        zab=zdab(ii)
        
        xbc=xdbc(ii)
        ybc=ydbc(ii)
        zbc=zdbc(ii)
        rrbc=1.d0/sqrt(xbc*xbc+ybc*ybc+zbc*zbc)
        
        xcd=xdcd(ii)
        ycd=ydcd(ii)
        zcd=zdcd(ii)
        
        xac=xab+xbc
        yac=yab+ybc
        zac=zab+zbc
        
c     construct first dihedral vector
        
        pbx=yab*zbc-zab*ybc
        pby=zab*xbc-xab*zbc
        pbz=xab*ybc-yab*xbc
        pb2=pbx*pbx+pby*pby+pbz*pbz
        rpb1=1.d0/sqrt(pb2)
        rpb2=rpb1*rpb1
        
c     construct second dihedral vector
        
        pcx=ybc*zcd-zbc*ycd
        pcy=zbc*xcd-xbc*zcd
        pcz=xbc*ycd-ybc*xcd
        pc2=pcx*pcx+pcy*pcy+pcz*pcz
        rpc1=1.d0/sqrt(pc2)
        rpc2=rpc1*rpc1
        
c     determine dihedral angle 
        
        pbpc=pbx*pcx+pby*pcy+pbz*pcz
        cosp=pbpc*rpb1*rpc1
        sinp=(xbc*(pcy*pbz-pcz*pby)+ybc*(pbx*pcz-pbz*pcx)+
     x    zbc*(pcx*pby-pcy*pbx))*(rpb1*rpc1*rrbc)
        
        phi=atan2(sinp,cosp)
        
c     avoid singularity in sinp
        
        sinp=sign(max(1.d-8,abs(sinp)),sinp)
        rsinp=1.d0/sinp
        
c     selection of potential energy function type
        
        kk=listdih(ii,1)
        
c     calculate potential energy and scalar force term
        
        if(keydih(kk).eq.1)then
          
c     key=1 for torsion dihedral potential
          
          dterm=prmdih(kk,1)*(1.d0+cos(prmdih(kk,3)*phi-
     x      prmdih(kk,2)))
          gamma=-rpb1*rpc1*rsinp*prmdih(kk,1)*prmdih(kk,3)*
     x      sin(prmdih(kk,3)*phi-prmdih(kk,2))
          
        else if(keydih(kk).eq.2)then
          
c     key=2 for harmonic improper dihedral
          
          phi=phi-prmdih(kk,2)
          phi=phi-nint(phi*rtwopi)*twopi
          dterm=0.5d0*prmdih(kk,1)*(phi*phi)
          gamma=rpb1*rpc1*rsinp*prmdih(kk,1)*phi
          
        else if(keydih(kk).eq.3)then
          
c     key=3 for harmonic cosine dihedral
          
          dterm=0.5d0*prmdih(kk,1)*(cos(phi)-
     x      cos(prmdih(kk,2)))**2
          gamma=-rpb1*rpc1*prmdih(kk,1)*(cos(phi)-cos(prmdih(kk,2)))
          
        else if(keydih(kk).eq.4)then
          
c     key=4 for 3-term cosine dihedral
          
          dterm=0.5d0*(prmdih(kk,1)*(1.d0+cos(phi))+
     x      prmdih(kk,2)*(1.d0-cos(2.d0*phi))+prmdih(kk,3)*
     x      (1.d0+cos(3.d0*phi)))
          gamma=-rpb1*rpc1*rsinp*0.5d0*(prmdih(kk,1)*sin(phi)-
     x      2.d0*prmdih(kk,2)*sin(2.d0*phi)+3.d0*prmdih(kk,3)*
     x      sin(3.d0*phi))
          
        else if(keydih(kk).eq.5)then
          
c     key=5 for ryckaert-bellemans potential      
c     chem.phys.lett., vol.30, p.123, 1975.
c     ATTENTION !!! Modified to have trans configuration
c     correspond to phi=180 rather than 
c     phi=0 as in original form.
          
          dterm=prmdih(kk,1)*(1.116d0-1.462d0*cos(phi)-
     x      1.578d0*(cos(phi))**2+0.368d0*(cos(phi))**3+
     x      3.156d0*(cos(phi))**4+3.788d0*(cos(phi))**5)
          gamma=prmdih(kk,1)*(1.462d0+3.156d0*cos(phi)-
     x      1.104d0*(cos(phi))**2-12.624d0*(cos(phi))**3-
     x      18.94d0*(cos(phi))**4)*rpb1*rpc1
          
        else if(keydih(kk).eq.6)then
          
c     key=6 for fluorinated ryckaert-bellemans potential
c     Rice at al., JCP 104, 2101, (1996).
          
          dterm=prmdih(kk,1)*(3.55d0-2.78d0*cos(phi)-
     x      3.56d0*(cos(phi))**2-1.64d0*(cos(phi))**3+
     x      7.13d0*(cos(phi))**4+12.84d0*(cos(phi))**5+
     x      9.67d0*exp(-56.d0*(phi-pi)**2))
          gamma=(prmdih(kk,1)*(2.78d0+7.12d0*cos(phi)+
     x      4.92d0*(cos(phi))**2-28.52d0*(cos(phi))**3-
     x      64.2d0*(cos(phi))**4)-1083.04d0*(phi-pi)*
     x      exp(-56.0*(phi-pi)**2))*rpb1*rpc1
          
        else if(keydih(kk).eq.7)then
          
c     key=7 for opls cosine dihedral
          
          phi=phi-prmdih(kk,5)
          dterm=prmdih(kk,1)+0.5d0*(prmdih(kk,2)*
     x      (1.d0+cos(phi))+prmdih(kk,3)*(1.d0-cos(2.d0*phi))+
     x      prmdih(kk,4)*(1.d0+cos(3.d0*phi)))
          gamma=-0.5d0*(prmdih(kk,2)*sin(phi)-2.d0*prmdih(kk,3)*
     x      sin(2.d0*phi)+3.d0*prmdih(kk,4)*sin(3.d0*phi))*rpb1*
     x      rpc1*rsinp
          
        else
          
c     undefined potential
          
          safe=.false.
          dterm=0.d0
          gamma=0.d0
          
        endif
        
c     indices of bonded atoms
        
        ia=listdih(ii,2)
        ib=listdih(ii,3)
        ic=listdih(ii,4)
        id=listdih(ii,5)
        
c     metadynamics local definitions
        
        if(lmetadyn)then
          
          idrive=driven(ltype(ia))
          jdrive=driven(ltype(ib))
          kdrive=driven(ltype(ic))
          ldrive=driven(ltype(id))
          
        endif
        
c     set selection control for angle potential
        
        lselect=.true.
        
        if(lexcite)then
          
c     selected excitation option
          
          if((atm_fre(ia).ne.1).and.(atm_fre(ib).ne.1).and.
     x      (atm_fre(ic).ne.1).and.(atm_fre(id).ne.1))then
            
c     reset selection control
            
            lselect=(atm_fre(ia)+atm_fre(ib)+atm_fre(ic)+
     x        atm_fre(id).eq.0)
            
            if(lsolva)then
              dih_exc(atmolt(ia))=dih_exc(atmolt(ia))+dterm
            endif
            
          endif
          
        elseif(lfree)then
          
c     selected free energy option
          
          if((atm_fre(ia).eq.1).or.(atm_fre(ib).eq.1).or.
     x      (atm_fre(ic).eq.1).or.(atm_fre(id).eq.1))then
            
c     set hamiltonian mixing parameter
            
            dih_fre=dih_fre-dterm
            dterm=lambda1*dterm
            gamma=lambda1*gamma
            
          elseif((atm_fre(ia).eq.2).or.(atm_fre(ib).eq.2).or.
     x        (atm_fre(ic).eq.2).or.(atm_fre(id).eq.2))then
            
c     set hamiltonian mixing parameter
            
            dih_fre=dih_fre+dterm
            dterm=lambda2*dterm
            gamma=lambda2*gamma
            
          endif
          
        endif
        
        if(lselect)then
          
c     calculate potential energy
          
          engdih=engdih+dterm
          
c     calculate solvation energy for dihedral term
          
          if(lsolva)then
            dih_sol(atmolt(ia))=dih_sol(atmolt(ia))+dterm
          endif
          
c     calculate atomic forces for dihedral term
          
          fax=gamma*((-pcy*zbc+pcz*ybc)-pbpc*rpb2*(-pby*zbc+pbz*ybc))
          fay=gamma*(( pcx*zbc-pcz*xbc)-pbpc*rpb2*( pbx*zbc-pbz*xbc))
          faz=gamma*((-pcx*ybc+pcy*xbc)-pbpc*rpb2*(-pbx*ybc+pby*xbc))
          
          fcx=gamma*((-pcy*zab+pcz*yab)-pbpc*rpb2*(-pby*zab+pbz*yab))
          fcy=gamma*(( pcx*zab-pcz*xab)-pbpc*rpb2*( pbx*zab-pbz*xab))
          fcz=gamma*((-pcx*yab+pcy*xab)-pbpc*rpb2*(-pbx*yab+pby*xab))
          
          fb1x=gamma*((-pby*zcd+pbz*ycd)-pbpc*rpc2*(-pcy*zcd+pcz*ycd))
          fb1y=gamma*(( pbx*zcd-pbz*xcd)-pbpc*rpc2*( pcx*zcd-pcz*xcd))
          fb1z=gamma*((-pbx*ycd+pby*xcd)-pbpc*rpc2*(-pcx*ycd+pcy*xcd))
          
          fd1x=gamma*((-pby*zbc+pbz*ybc)-pbpc*rpc2*(-pcy*zbc+pcz*ybc))
          fd1y=gamma*(( pbx*zbc-pbz*xbc)-pbpc*rpc2*( pcx*zbc-pcz*xbc))
          fd1z=gamma*((-pbx*ybc+pby*xbc)-pbpc*rpc2*(-pcx*ybc+pcy*xbc))
          
          fxx(ia)=fxx(ia)+fax
          fyy(ia)=fyy(ia)+fay
          fzz(ia)=fzz(ia)+faz
          
          fxx(ib)=fxx(ib)-fax-fcx+fb1x
          fyy(ib)=fyy(ib)-fay-fcy+fb1y
          fzz(ib)=fzz(ib)-faz-fcz+fb1z
          
          fxx(ic)=fxx(ic)+fcx-fb1x-fd1x
          fyy(ic)=fyy(ic)+fcy-fb1y-fd1y
          fzz(ic)=fzz(ic)+fcz-fb1z-fd1z
          
          fxx(id)=fxx(id)+fd1x
          fyy(id)=fyy(id)+fd1y
          fzz(id)=fzz(id)+fd1z
          
c     stress tensor for dihedral term
          
          strs(1)=strs(1)+xab*fax+xbc*(fb1x-fcx)-xcd*fd1x 
          strs(2)=strs(2)+yab*fax+ybc*(fb1x-fcx)-ycd*fd1x 
          strs(3)=strs(3)+zab*fax+zbc*(fb1x-fcx)-zcd*fd1x 
          strs(4)=strs(4)+yab*fay+ybc*(fb1y-fcy)-ycd*fd1y 
          strs(5)=strs(5)+yab*faz+ybc*(fb1z-fcz)-ycd*fd1z 
          strs(6)=strs(6)+zab*faz+zbc*(fb1z-fcz)-zcd*fd1z 
          
        endif
        
c     metadynamics local parameters
        
        if(lmetadyn.and.(idrive.or.jdrive.or.kdrive.or.ldrive))then
          
c     local energy (no virial)
          
          eng_loc=eng_loc+dterm
          
c     local forces
          
          fxx_loc(ia)=fxx_loc(ia)+fax
          fyy_loc(ia)=fyy_loc(ia)+fay
          fzz_loc(ia)=fzz_loc(ia)+faz
          
          fxx_loc(ib)=fxx_loc(ib)-fax-fcx+fb1x
          fyy_loc(ib)=fyy_loc(ib)-fay-fcy+fb1y
          fzz_loc(ib)=fzz_loc(ib)-faz-fcz+fb1z
          
          fxx_loc(ic)=fxx_loc(ic)+fcx-fb1x-fd1x
          fyy_loc(ic)=fyy_loc(ic)+fcy-fb1y-fd1y
          fzz_loc(ic)=fzz_loc(ic)+fcz-fb1z-fd1z
          
          fxx_loc(id)=fxx_loc(id)+fd1x
          fyy_loc(id)=fyy_loc(id)+fd1y
          fzz_loc(id)=fzz_loc(id)+fd1z
          
c     local stress tensor
          
          strs_loc(1)=strs_loc(1)+xab*fax+xbc*(fb1x-fcx)-xcd*fd1x
          strs_loc(2)=strs_loc(2)+yab*fax+ybc*(fb1x-fcx)-ycd*fd1x
          strs_loc(3)=strs_loc(3)+zab*fax+zbc*(fb1x-fcx)-zcd*fd1x
          strs_loc(4)=strs_loc(4)+yab*fay+ybc*(fb1y-fcy)-ycd*fd1y
          strs_loc(5)=strs_loc(5)+yab*faz+ybc*(fb1z-fcz)-ycd*fd1z
          strs_loc(6)=strs_loc(6)+zab*faz+zbc*(fb1z-fcz)-zcd*fd1z
          
        endif
        
c     calculate 1-4 dihedral interactions (coulombic and short ranged)
c     assumes 1-4 interactions are in the exclude list
        
        kk=listdih(ii,1)
        
c     bypass OPLS 1-4 terms (not present)
        
        if(keydih(kk).ne.7)then
          
c     1-4 electrostatics : adjust by weighting factor
          
          scale=prmdih(kk,4)
          
          xad=xac+xcd
          yad=yac+ycd
          zad=zac+zcd
          
          rad=sqrt(xad**2+yad**2+zad**2)
          
c     scaled charge product*dielectric
          
          chgprd=scale*chge(ia)*chge(id)*r4pie0
          coul=0.d0
          fcoul=0.d0
          
c     truncation of potential for all schemes except ewald sum
          
          if(abs(chgprd).gt.1.d-10.and.keyfce.gt.0)then
            
c     electrostatics by ewald sum
            
            if(keyfce/2.eq.1.or.keyfce/2.eq.6.or.keyfce/2.eq.7)then
              
              coul=chgprd/(epsq*rad)
              fcoul=coul/(rad**2)
              
c     distance dependent dielectric
              
            elseif(rcut.gt.rad)then
              
              if(keyfce/2.eq.2)then
                
                coul=chgprd/(epsq*rad**2)
                fcoul=2.0d0*coul/(rad**2)
                
c     unmodified coulombic
                
              else if(keyfce/2.eq.3)then
                
                coul=chgprd/(epsq*rad)
                fcoul=coul/(rad**2)
                
c     truncated and shifted coulombic
                
              else if(keyfce/2.eq.4)then
                
                tt=1.d0/(1.d0+pp*alpha*rad)
                exp1=exp(-(alpha*rad)**2)
                erc=tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1/rad
                fer=(erc+2.d0*(alpha/sqrpi)*exp1)/rad**2
                
c     calculate potential energy and forces
                
                coul=chgprd*(erc-vcon+fcon*(rad-rcut))/epsq
                fcoul=chgprd*(fer-fcon/rad)/epsq
                
c     reaction field
                
              else if(keyfce/2.eq.5)then
                
                tt=1.d0/(1.d0+pp*alpha*rad)
                exp1=exp(-(alpha*rad)**2)
                erc=tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1/rad
                fer=(erc+2.d0*(alpha/sqrpi)*exp1)/rad**2
                coul=chgprd*(erc-vcon+fcon*(rad-rcut)+
     x            rfld2*rad*rad-rfld1)
                fcoul=chgprd*(fer-fcon/rad-rfld0)
                
              elseif(keyfce/2.eq.0)then
                
                coul=0.d0
                fcoul=0.d0
                
              else
                
                call error(idnode,446)
                
              endif
              
            endif
            
c     set selection control for coulombic 1-4 terms
            
            lselect=.true.
            
c     set double index
            
            if(lsolva)kkk=loc2(atmolt(ia),atmolt(id))
            
            if(lexcite)then
              
c     selected excitation option
              
              if((atm_fre(ia).ne.1).and.(atm_fre(id).ne.1))then
                
c     set selection control
                
                lselect=(atm_fre(ia)+atm_fre(id).eq.0)
                
                if(lsolva)cou_exc(kkk)=cou_exc(kkk)+coul
                
              endif
              
            elseif(lfree)then
              
c     selected free energy option
              
              if((atm_fre(ia).eq.1).or.(atm_fre(id).eq.1))then
                
c     set hamiltonian mixing parameter
                
                cou_fre=cou_fre-coul
                cou14_vir=cou14_vir+fcoul*rad**2
                coul=lambda1*coul
                fcoul=lambda1*fcoul
                
              elseif((atm_fre(ia).eq.2).or.(atm_fre(id).eq.2))then
                
c     set hamiltonian mixing parameter
                
                cou_fre=cou_fre+coul
                cou14_vir=cou14_vir-fcoul*rad**2
                coul=lambda2*coul
                fcoul=lambda2*fcoul
                
              endif
              
            endif
            
            if(lselect)then
              
c     electrostatic energy and virial for 1-4 term
              
              engc14=engc14+coul
              virc14=virc14-fcoul*rad**2
              
c     solvation energy for coulombic 1-4 term
              
              if(lsolva)cou_sol(kkk)=cou_sol(kkk)+coul
              
c     coulombic force for 1-4 term
              
              fx=fcoul*xad
              fy=fcoul*yad
              fz=fcoul*zad
              
              fxx(ia)=fxx(ia)+fx
              fyy(ia)=fyy(ia)+fy
              fzz(ia)=fzz(ia)+fz
              
              fxx(id)=fxx(id)-fx
              fyy(id)=fyy(id)-fy
              fzz(id)=fzz(id)-fz
              
c     stress tensor for coulombic 1-4 term
              
              strs(1)=strs(1)+xad*fx
              strs(2)=strs(2)+xad*fy
              strs(3)=strs(3)+xad*fz
              strs(4)=strs(4)+yad*fy
              strs(5)=strs(5)+yad*fz
              strs(6)=strs(6)+zad*fz
            
            endif
            
c     metadynamics local parameters
            
            if(lmetadyn.and.(idrive.or.ldrive))then
              
c     local energy and virial
              
              eng_loc=eng_loc+coul
              vir_loc=vir_loc-fcoul*rad**2
              
c     local forces
              
              fxx_loc(ia)=fxx_loc(ia)+fx
              fyy_loc(ia)=fyy_loc(ia)+fy
              fzz_loc(ia)=fzz_loc(ia)+fz
              
              fxx_loc(id)=fxx_loc(id)-fx
              fyy_loc(id)=fyy_loc(id)-fy
              fzz_loc(id)=fzz_loc(id)-fz
              
c     local stress tensor
              
              strs_loc(1)=strs_loc(1)+xad*fx
              strs_loc(2)=strs_loc(2)+xad*fy
              strs_loc(3)=strs_loc(3)+xad*fz
              strs_loc(4)=strs_loc(4)+yad*fy
              strs_loc(5)=strs_loc(5)+yad*fz
              strs_loc(6)=strs_loc(6)+zad*fz
              
            endif
            
          endif
          
c     1-4 short ranged : adjust by weighting factor
          
          scale=prmdih(kk,5)
          gamma=0.d0
          
          if(mod(keyfce,2).eq.1)then
            
c     atomic and potential function indices
            
            ka=max(ltype(ia),ltype(id))
            kb=min(ltype(ia),ltype(id))
            k=lstvdw((ka*(ka-1))/2+kb)
            
            if(abs(scale*vvv(1,k)).gt.1.d-10)then
              
c     apply truncation of potential
              
              if(rvdw.gt.rad)then
                
c     determine interpolation panel for force arrays
                
                l=int(rad/dlrpot)
                ppp=rad/dlrpot-dble(l)
                
c     calculate interaction energy using 3-point interpolation
                
                vk0=vvv(l,k)
                vk1=vvv(l+1,k)
                vk2=vvv(l+2,k)
                
                t1=vk0+(vk1-vk0)*ppp
                t2=vk1+(vk2-vk1)*(ppp-1.0d0)
                
                srpot=scale*(t1+(t2-t1)*ppp*0.5d0)
                
c     calculate forces using 3-point interpolation
                
                gk0=ggg(l,k)
                gk1=ggg(l+1,k)
                gk2=ggg(l+2,k)
                
                t1=gk0+(gk1-gk0)*ppp
                t2=gk1+(gk2-gk1)*(ppp-1.0d0)
                
                gamma=scale*(t1+(t2-t1)*ppp*0.5d0)/(rad**2)
                
c     set selection control for short ranged 1-4 terms
                
                lselect=.true.
                
c     set double index
                
                if(lsolva)kkk=loc2(atmolt(ia),atmolt(id))
                
                if(lexcite)then
                  
c     selected excitation option
                  
                  if((atm_fre(ia).ne.1).and.(atm_fre(id).ne.1))then
                    
c     set selection control
                    
                    lselect=(atm_fre(ia)+atm_fre(id).eq.0)
                    
                    if(lsolva)vdw_exc(kkk)=vdw_exc(kkk)+srpot
                    
                  endif
                  
                elseif(lfree)then
                  
c     selected free energy option
                  
                  if((atm_fre(ia).eq.1).or.(atm_fre(id).eq.1))then
                    
c     set hamiltonian mixing parameter
                    
                    vdw_fre=vdw_fre-srpot
                    vdw14_vir=vdw14_vir+gamma*rad**2
                    srpot=lambda1*srpot
                    gamma=lambda1*gamma
                    
                  elseif((atm_fre(ia).eq.2).or.(atm_fre(id).eq.2))then
                    
c     set hamiltonian mixing parameter
                    
                    vdw_fre=vdw_fre+srpot
                    vdw14_vir=vdw14_vir-gamma*rad**2
                    srpot=lambda2*srpot
                    gamma=lambda2*gamma
                    
                  endif
                  
                endif
                
                if(lselect)then
                  
c     short ranged energy and virial for 1-4 term
                  
                  engs14=engs14+srpot
                  virs14=virs14-gamma*rad**2
                  
c     solvation energy for short ranged 1-4 term
                  
                  if(lsolva)vdw_sol(kkk)=vdw_sol(kkk)+srpot
                  
c     short ranged forces for 1-4 term
                  
                  fx=gamma*xad
                  fy=gamma*yad
                  fz=gamma*zad
                  
                  fxx(ia)=fxx(ia)+fx
                  fyy(ia)=fyy(ia)+fy
                  fzz(ia)=fzz(ia)+fz
                  
                  fxx(id)=fxx(id)-fx
                  fyy(id)=fyy(id)-fy
                  fzz(id)=fzz(id)-fz
                  
c     stress tensor for short ranged 1-4 term
                  
                  strs(1)=strs(1)+xad*fx
                  strs(2)=strs(2)+xad*fy
                  strs(3)=strs(3)+xad*fz
                  strs(4)=strs(4)+yad*fy
                  strs(5)=strs(5)+yad*fz
                  strs(6)=strs(6)+zad*fz
                  
                endif
                
c     metadynamics local parameters
                  
                if(lmetadyn.and.(idrive.or.ldrive))then
                  
c     local energy and virial
                  
                  eng_loc=eng_loc+srpot
                  vir_loc=vir_loc-gamma*rad**2
                  
c     local forces
                  
                  fxx_loc(ia)=fxx_loc(ia)+fx
                  fyy_loc(ia)=fyy_loc(ia)+fy
                  fzz_loc(ia)=fzz_loc(ia)+fz
                  
                  fxx_loc(id)=fxx_loc(id)-fx
                  fyy_loc(id)=fyy_loc(id)-fy
                  fzz_loc(id)=fzz_loc(id)-fz
                  
c     stress tensor for short ranged 1-4 term
                  
                  strs_loc(1)=strs_loc(1)+xad*fx
                  strs_loc(2)=strs_loc(2)+xad*fy
                  strs_loc(3)=strs_loc(3)+xad*fz
                  strs_loc(4)=strs_loc(4)+yad*fy
                  strs_loc(5)=strs_loc(5)+yad*fz
                  strs_loc(6)=strs_loc(6)+zad*fz
                  
                endif
                
              endif
              
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
      
c     sum contributions to potentials
      
      if(mxnode.gt.1)then
        
        buffer(1)=engdih
        buffer(2)=engc14
        buffer(3)=virc14
        buffer(4)=engs14
        buffer(5)=virs14
        buffer(6)=dih_fre
        buffer(7)=dih_vir
        buffer(8)=vdw14_vir
        buffer(9)=cou14_vir
        call gdsum(buffer(1),9,buffer(10))
        engdih=buffer(1)
        engc14=buffer(2)
        virc14=buffer(3)
        engs14=buffer(4)
        virs14=buffer(5)
        dih_fre=buffer(6)
        dih_vir=buffer(7)
        vdw14_vir=buffer(8)
        cou14_vir=buffer(9)
        
        if(lsolva)then
          
          call gdsum(dih_sol(1),mxtmls,buffer(1))
          if(lexcite)call gdsum(dih_exc(1),mxtmls,buffer(1))
          
        endif
        
      endif
      
      engcpe=engcpe+engc14
      vircpe=vircpe+virc14
      engsrp=engsrp+engs14
      virsrp=virsrp+virs14
      vdw_vir=vdw_vir+vdw14_vir
      cou_vir=cou_vir+cou14_vir
      
c     check for undefined potentials
      
      if(mxnode.gt.1)call gstate(safe)
      if(.not.safe)call error(idnode,448)
      
      deallocate (xdab,ydab,zdab,stat=fail1)
      deallocate (xdbc,ydbc,zdbc,stat=fail2)
      deallocate (xdcd,ydcd,zdcd,stat=fail3)
      
      return
      end subroutine dihfrc
      
      end module dihedral_module
