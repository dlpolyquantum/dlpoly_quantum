      module property_module

c***********************************************************************
c     
c     dl_poly module for defining simulation property data
c     copyright - daresbury laboratory
c     author    - w. smith    nov 2003
c     
c***********************************************************************

      use config_module
      use error_module
      use pair_module
      use setup_module
      use site_module
      use tether_module
      use utility_module
      use vdw_module

      implicit none

      real(8), allocatable :: rdf(:,:),zdens(:,:)
      real(8), allocatable :: stpval(:),sumval(:)
      real(8), allocatable :: ssqval(:),zumval(:)
      real(8), allocatable :: ravval(:),stkval(:,:)
      real(8), allocatable :: xx0(:),yy0(:),zz0(:)
      real(8), allocatable :: amsd(:)

      save rdf,zdens,stpval,sumval,ssqval,xx0,yy0,zz0
      save zumval,ravval,stkval

      contains
      
      subroutine alloc_prp_arrays(idnode,mxnode)

      implicit none

      integer, parameter :: nnn=6

      logical safe
      integer i,fail,idnode,mxnode,numatm
      dimension fail(nnn)

      safe=.true.
      numatm=nbeads*mxatms

c     allocate arrays

      fail(:)=0
      
      allocate (zdens(mxzdn,mxatyp),stat=fail(1))
      allocate (rdf(mxrdf,mxxtyp),amsd(mxatyp),stat=fail(2))
      allocate (stpval(mxnstk),sumval(mxnstk),stat=fail(3))
      allocate (ssqval(mxnstk),zumval(mxnstk),stat=fail(4))
      allocate (ravval(mxnstk),stkval(mxstak,mxnstk),stat=fail(5))
      allocate (xx0(numatm),yy0(numatm),zz0(numatm),stat=fail(6))
      
      if(any(fail.gt.0))safe=.false.      
      if(mxnode.gt.1)call gstate(safe)    
      if(.not.safe)call error(idnode,1740)
      
      end subroutine alloc_prp_arrays

      subroutine result
     x  (ltad,lbpd,lgofr,lpgr,lzden,lpimd,idnode,imcon,keyens,mxnode,
     x  natms,nbeads,levcfg,nzden,nstep,ntpatm,numacc,numrdf,keybpd,
     x  chip,chit,conint,rcut,tstep,engcfg,volm,virtot,vircom,zlen,
     x  tboost,chit_shl,gaumom)

c***********************************************************************
c     
c     dl_poly subroutine for writing simulation summary and
c     saving the restart data
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith dec 1992.
c     
c***********************************************************************

      implicit none

      character*1 hms,dec
      logical lgofr,lpgr,lzden,check,ltad,lbpd,lpimd,goprint
      
      integer idnode,imcon,keyens,mxnode,natms,nbeads,nzden,nstep
      integer ntpatm,levcfg,numacc,numrdf,keybpd,i,iadd,io,j
      real(8) chip,chit,conint,rcut,tstep,volm,timelp,avvol,zlen,dc
      real(8) engcfg,virtot,vircom,prntim,simtim,tboost,chit_shl
      real(8) rgau,gaumom(0:5)

c     sum the gaussian moments
      
      if(mxnode.gt.1)then
        
        buffer(1)=gaumom(0)
        call gdsum(buffer(1),1,buffer(2))
        
        if(nint(buffer(1)).eq.0)then
          rgau=0.d0
        else
          rgau=1.d0/buffer(1)
        endif
        
        do i=1,5
          
          gaumom(i)=(gaumom(0)*rgau)*gaumom(i)
          
        enddo
        
        call gdsum(gaumom(0),6,buffer(1))
        
      endif

c     save restart data
      
      call revive
     x  (lgofr,lzden,lpimd,idnode,imcon,mxnode,natms,nbeads,levcfg,
     x  nstep,nzden,numacc,numrdf,chip,chit,conint,tstep,engcfg,virtot,
     x  vircom,tboost,chit_shl,gaumom)

c     for TAD and BPD system averages not generally meaningful 
c     useful only for BPD in configurational sampling mode
      
      goprint=.not.(ltad.or.(lbpd.and.keybpd.gt.1))
      
      if(goprint)then
        
c     calculate final fluctuations
        
        do i=1,mxnstk
          ssqval(i)=sqrt(max(0.d0,ssqval(i)))
        enddo
        
c     final averages and fluctuations
        
        call timchk(0,timelp)
        
        if(idnode.eq.0)then
          
          write(nrite,
     x      "(/,/,1x,'run terminated after',i8,' steps.',
     x      ' final averages calculated over',i8,' steps.',/,/)") 
     x      nstep,numacc
          write(nrite,"(1x,120('-'),
     x      /,/,1x,'    step',5x,'eng_tot',4x,'temp_tot',5x,'eng_cfg',
     x      5x,'eng_vdw',5x,'eng_cou',5x,'eng_bnd',5x,'eng_ang',5x,
     x      'eng_dih',5x,'eng_tet',/,1x,'time    ',5x,' eng_pv',4x,
     x      'temp_rot',5x,'vir_cfg',5x,'vir_vdw',5x,'vir_cou',5x,
     x      'vir_bnd',5x,'vir_ang',5x,'vir_con',5x,'vir_tet',/,
     x      1x,'cpu time',6x,'volume',4x,'temp_shl',5x,'eng_shl',
     x      5x,'vir_shl',7x,'alpha',8x,'beta',7x,'gamma',5x,'vir_pmf',
     x    7x,'press')")
          if(lpimd)write(nrite,"(1x,'pimd    ',5x,'eng_qpi',5x,
     x      'eng_qvr',5x,'eng_rng',5x,'vir_rng',5x,'qms_rgr',5x,
     x      'qms_bnd',5x,'eng_the')")
          write(nrite,"(/,/,1x,120('-'))")
          
          call get_prntime(hms,timelp,prntim)
          call get_simtime(dec,nstep,tstep,simtim)
          write(nrite,'(1x,i8,1p,9e12.4,/,1x,0p,f7.3,a1,1p,9e12.4,
     x      /,1x,0p,f7.3,a1,1p,9e12.4)') 
     x      nstep,(sumval(i),i=1,9),
     x      simtim,dec,(sumval(i),i=10,18),
     x      prntim,hms,(sumval(i),i=19,27)
          iadd=mxatyp+45
          if(lpimd)then
            write(nrite,'(9x,1p,9e12.4)')(sumval(iadd+i),i=1,7)
          endif
          write(nrite,"(/,1x,' r.m.s. ',1p,9e12.4,/,1x,'fluctn. ',
     x      1p,9e12.4,/,9x,9e12.4)") (ssqval(i),i=1,27)
          if(lpimd)then
            write(nrite,'(9x,1p,9e12.4)')(ssqval(iadd+i),i=1,7)
          endif
          write(nrite,"(1x,120('-'))")
          
c     write out bias potential boost factor
          
          if(lbpd)write(nrite,"(/,/,1x,
     x      'calculated bias potential boost factor',1p,e16.8)")tboost
          
          if(numacc.gt.0)then
            
            iadd=27
            
c     write out estimated diffusion coefficients
            
            if(numacc.gt.0)then
              
              write(nrite,"(/,/,12x,'Approximate 3D Diffusion',
     x          '  coefficients (10^-9 m^2 / s)',/,/,12x,'atom',7x,
     x          ' D ')")
              
              do i=1,ntpatm
                
                dc=(ravval(iadd)-sumval(iadd))/
     x            (3.d0*dble(numacc-min(mxnstk,numacc-1))*tstep)*10.d0
                if(dc.lt.1d-10) dc=0.d0
                if(lbpd)dc=dc/tboost
                write(nrite,'(12x,a8,1p,e13.4)') unqatm(i),dc
                
              enddo
              
            endif

            iadd=iadd+mxatyp
            
c     print out average pressure tensor
            
            write(nrite,"(/,/,16x,'Average pressure tensor',
     x        39x,'r.m.s. fluctuations ',/)")
            
            do i=iadd,iadd+6,3
              write(nrite,'(9x,1p,3e12.4,24x,3e12.4)')
     x          (sumval(i+j),j=1,3),(ssqval(i+j),j=1,3)
            enddo
            
            iadd=iadd+9
            
            write(nrite,'(/,12x,a,1p,e12.4)') 'trace/3. ',
     x        (sumval(iadd)+sumval(iadd-4)+sumval(iadd-8))/3.d0
            
c     write out mean cell vectors for npt 
            
            if(keyens.gt.3.and.(keyens.le.7))then
              
              write(nrite,"(/,/,17x,'Average cell vectors',
     x          41x,'r.m.s. fluctuations ',/)")
              
              do i=iadd,iadd+6,3
                write(nrite,'(3f20.10,9x,1p,3e12.4)')
     x            (sumval(i+j),j=1,3),(ssqval(i+j),j=1,3)
              enddo
              
            endif
            
            iadd=iadd+9
            if(lpimd)iadd=iadd+7
            
c     write out remaining nonzero registers (if any)
            
            check=.false.
            do i=iadd+1,mxnstk
              
              if((abs(sumval(i)).gt.1.d-10).or.
     x          (abs(ssqval(i)).gt.1.d-10)) check=.true.
              
            enddo
            
            if(check)then
              
              write(nrite,"(/,/,12x,
     x          'Remaining non-zero statistics registers ',/,/,12x,
     x          'Register',7x,'Average value',8x,'r.m.s. fluc.')")
              
              do i=iadd+1,mxnstk
                
                if((abs(sumval(i)).gt.1.d-10).or.
     x            (abs(ssqval(i)).gt.1.d-10))
     x            write(nrite,'(10x,i10,2f20.10)') i,sumval(i),ssqval(i)
                
              enddo
              
            endif
            
          endif
          
c     print out gaussian moments of momenta
          
          write(nrite,'(/,/,1x,a34,1p,e13.6,a8)')
     x      "gaussian moments of momenta. over ",
     x      gaumom(0)," samples"
          write(nrite,'(1p,5e16.6)')(gaumom(i),i=1,5)
          
        endif
        
      endif
      
c     print out sample of final configuration 
      
      if(idnode.eq.0)then
        
        write(nrite,"(/,/,1x,'sample of final configuration',/)")
        write(nrite,"(6x,'i',7x,'x(i)',8x,'y(i)',8x,'z(i)',
     x7x,'vx(i)',7x,'vy(i)',7x,'vz(i)',7x,'fx(i)',7x,
     x'fy(i)',7x,'fz(i)',/,/)")
        io=(natms*nbeads+19)/20
        
        do i=1,natms*nbeads,io
          
          write(nrite,"(1x,i6,1p,3e12.4,3e12.4,3e12.4)") 
     x      i,xxx(i),yyy(i),zzz(i),vxx(i),vyy(i),vzz(i),
     x      fxx(i),fyy(i),fzz(i)
          
        enddo

      endif
      
c     bypass printing averages for certain tad and bpd options
      
      if(goprint)then
        
c     average volume
        
        avvol=sumval(19)
        if(imcon.eq.0.or.imcon.eq.6)then
          avvol=4.d0*pi/3.d0*rcut**3
          volm=avvol
        endif
        
c     calculate and print radial distribution functions
        
        if(lgofr.and.lpgr.and.(numrdf.gt.0))then 
          
c     scale densities for average volume
          
          do i=1,ntpatm
            dens(i)=dens(i)*(volm/avvol)
          enddo
          
          call rdf1
     x      (lpgr,idnode,mxnode,ntpatm,numrdf,avvol,rcut)
          
        endif
        
        if(lzden.and.lpgr.and.(nzden.gt.0))then 
          call zden1(lpgr,idnode,mxnode,ntpatm,nzden,avvol,zlen)
        endif
        
        if(imcon.eq.0)volm=0.d0
        
      endif
      
c     print final time check
        
      call timchk(1,timelp)
        
      return
      end subroutine result

      subroutine diffsn0(idnode,natms,mxnode,tstep)

c***********************************************************************
c     
c     DL_POLY routine for calculating displacements of sites from
c     t=0 positions
c     
c     use diffsn1 for mean squared displacements
c     
c     parallel version - replicated data.
c     
c     copyright daresbury laboratory 1993
c     
c     author - t. forester      june 1993
c     
c***********************************************************************

      implicit none

      logical newjob 
      integer idnode,natms,mxnode,iatm1,iatm2,i
      real(8) tstep

      save newjob,iatm1,iatm2
      data newjob/.true./

      if(newjob)then

        newjob=.false.
        iatm1=(idnode*natms)/mxnode+1
        iatm2=((idnode+1)*natms)/mxnode

      endif

      do i=iatm1,iatm2

        xx0(i)=xx0(i)+vxx(i)*tstep
        yy0(i)=yy0(i)+vyy(i)*tstep
        zz0(i)=zz0(i)+vzz(i)*tstep
        
      enddo
      
      return
      end subroutine diffsn0

      subroutine diffsn1(idnode,natms,ntpatm,mxnode)
      
c***********************************************************************
c     
c     DL_POLY routine for calculating mean squared displacements
c     
c     displacements calculated in diffsn0
c     
c     parallel version - replicated data.
c     
c     copyright daresbury laboratory 1993
c     
c     author - t. forester      june 1993
c     
c***********************************************************************
      
      implicit none

      logical newjob
      integer idnode,natms,ntpatm,mxnode,iatm1,iatm2,k,i

      save newjob,iatm1,iatm2

      data newjob/.true./

      if(newjob)then

        newjob=.false.
        iatm1=(idnode*natms)/mxnode+1
        iatm2=((idnode+1)*natms)/mxnode

      endif

c     running sum of squared displacements
      
      do k=1,ntpatm
        
        amsd(k)=0.d0
        
      enddo
      
c     calculate square of displacements for each atom type

      do i=iatm1,iatm2
        
        k=ltype(i)
        amsd(k)=amsd(k)+xx0(i)**2+yy0(i)**2+zz0(i)**2
        
      enddo
      
c     global sum - replicated data strategy
      
      if(mxnode.gt.1)then
        
        do k=1,ntpatm
          
          buffer(k+ntpatm)=amsd(k)
          
        enddo
        
        call  gdsum(buffer(1+ntpatm),ntpatm,buffer(1))
        
        do k=1,ntpatm
          
          amsd(k)=buffer(k+ntpatm)
          
        enddo
        
      endif
      
c     mean squared displacement
      
      do k=1,ntpatm
        
        amsd(k)=amsd(k)/dble(max(numtyp(k),1))
        
      enddo
      
      return
      end subroutine diffsn1

      subroutine rdf0(iatm,ik,rcut)

c***********************************************************************
c     
c     dl_poly subroutine for accumulating statistic for radial
c     distribution functions.
c     double precision accumulators
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1994
c     author    - t. forester    march 1994    
c     
c***********************************************************************
      
      implicit none

      integer iatm,ik,m,jatm,ll,k
      real(8) rcut,rcsq,rdelr,ai,aj,rsq,rrr
      
c     set cutoff condition for pair forces
      
      rcsq=rcut*rcut

c     grid interval for rdf tables

      rdelr=dble(mxrdf)/rcut

c     set up atom iatm type

      ai=ltype(iatm) 

c     start of primary loop for rdf accumulation

      do m=1,ik

c     atomic and potential function indices
        
        jatm=ilist(m)

        aj=ltype(jatm)
        if(ai.gt.aj)then
          k=int(ai*(ai-1.d0)*0.5d0+aj+0.5d0)
        else
          k=int(aj*(aj-1.d0)*0.5d0+ai+0.5d0)
        endif

c     apply truncation of potential
        
        rsq=rsqdf(m)
        
        if(rcsq.gt.rsq)then

          rrr=sqrt(rsq)
          ll=Min(1+Int(rrr*rdelr),mxrdf)

c     accumulate statistics

          rdf(ll,k)=rdf(ll,k)+1.d0

        endif
        
      enddo
      
      return
      end subroutine rdf0

      subroutine rdf1
     x   (lpgr,idnode,mxnode,ntpatm,numrdf,volm,rcut)

c***********************************************************************
c     
c     dl_poly subroutine for calculating radial distribution functions
c     from accumulated data.
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1994
c     author    - t. forester    march 1994    
c     
c***********************************************************************
      
      implicit none

      integer idnode,mxnode,ntpatm,numrdf,ia,ib,k,j
      real(8) volm,factor,sum,rrr,dvol,gofr,rcut,delrdf
      
      logical lpgr,zero

      if(idnode.eq.0) write(nrite,
     x   "(/,/,12X,'RADIAL DISTRIBUTION FUNCTIONS',/,/,
     x   'calculated using ',i10,' configurations')") numrdf

      if(lpgr)then
        
c     open RDF file and write headers
        
        if(idnode.eq.0)then
          
          open(nrdfdt,file='RDFDAT')
          
          write(nrdfdt,'(80a1)')cfgname
          write(nrdfdt,'(2i10)')mxxtyp,mxrdf
          
        endif

c     default bin width
        
        delrdf=rcut/dble(mxrdf)
        
c     construct rdf tables
        
        do ia=1,ntpatm
          
          do ib=ia,ntpatm
            
            k=(ib*(ib-1))/2+ia
            
            if(idnode.eq.0)then
                
              write(nrite,
     x         "(/,'g(r)  :',2a8,/,/,8x,'r',6x,'g(r)',9x,'n(r)',/)") 
     x           unqatm(ia),unqatm(ib)
              write(nrdfdt,'(2a8)')unqatm(ia),unqatm(ib)
              
            endif
              
c     global sum of data on all nodes
            
            if(mxnode.gt.1) call gdsum(rdf(1,k),mxrdf,buffer)
            
c     normalisation factor 
            
            factor=volm*dens(ia)*dens(ib)*dble(numrdf)
            if((ia.eq.ib).and.(volm*dens(ia).gt.1.d0)) 
     x         factor=factor*0.5d0
            
c     running integration of rdf
            
            sum=0.d0
              
c     loop over distances
              
            zero=.true.
            
            do j=1,mxrdf
              
              if(zero.and.(j.lt.mxrdf-3)) 
     x          zero=(rdf(j+2,k).le.0.d0)
              
              rrr=(dble(j)-0.5d0)*delrdf
              dvol=4.d0*pi*(delrdf*rrr*rrr+(delrdf**3)/12.d0)
              
              gofr=rdf(j,k)/(factor*dvol)
              sum=sum+gofr*dvol*dens(ib)
              
c     print out information
                
              if(idnode.eq.0)then
                
                write(nrdfdt,"(1p,2e14.6)")rrr,gofr
                if(.not.zero)
     x             write(nrite,"(f10.4,1p,2e14.6)")rrr,gofr,sum
                
              endif
              
            enddo
              
          enddo
          
        enddo
        
        if(idnode.eq.0)close (nrdfdt)
        
      endif
      
      return
      end subroutine rdf1

      subroutine static
     x  (lbpd,lzeql,lpimd,idnode,intsta,imcon,keyens,natms,nstack,
     x  nstep,nsteql,ntpatm,numacc,mxnode,nblock,keybpd,numbpd,
     x  consv,degfre,degrot,engang,engbnd,engcpe,engdih,enginv,
     x  engke,engrot,engsrp,engunit,engcfg,stpeng,stpeth,stpprs,
     x  stptmp,stpvir,stpvol,tstep,virbnd,engfbp,vircom,vircon,
     x  vircpe,virsrp,engfld,virfld,engtbp,virtbp,virpmf,virshl,
     x  engshl,engtet,virtet,degshl,shlke,virang,width,engmet,
     x  virmet,engter,virter,boost,tboost,engqpi,engqvr,engrng,
     x  virrng,qmsrgr,qmsbnd,engthe)

c***********************************************************************
c     
c     dl_poly subroutine for accumulating periodic data during the
c     molecular dynamics simulation and computing the rolling averages
c     
c     copyright daresbury laboratory 1992
c     
c     author - w. smith       august 1992
c     
c***********************************************************************

c      use vv_pimd_module, only:press2

      implicit none

      logical lbpd,lzeql,lpimd,newjob
      integer idnode,intsta,imcon,keyens,natms,nstack,nstep,j
      integer nsteql,ntpatm,numacc,mxnode,i,iadd,k,kstak
      integer nblock,keybpd,numbpd
      real(8) consv,degfre,degrot,engang,engbnd,engcpe,engdih
      real(8) enginv,engke,engrot,engsrp,engunit,engcfg,stpeng
      real(8) stpeth,stpprs,stptmp,stpvir,stpvol,tstep,virbnd
      real(8) engfbp,vircom,vircon,vircpe,virsrp,engfld,virfld
      real(8) engtbp,virtbp,virpmf,virshl,engshl,engtet,virtet
      real(8) degshl,shlke,virang,width,sclnv1,sclnv2,stprot
      real(8) stpcns,stpshl,zistk,engmet,virmet,engter,virter
      real(8) tbold,aterm,bterm,cterm,boost,tboost,engqpi,engqvr
      real(8) engrng,virrng,qmsrgr,qmsbnd,engthe

      save newjob
      
      data newjob/.true./

c     open statistics file for append
      
      if(newjob.and.idnode.eq.0.and.intsta.gt.0)then

        open(nstats,file='STATIS',position='append')
        newjob=.false.
        
      endif
      
      if(idnode.eq.0.and.nstep.eq.intsta.and.intsta.gt.0)then
        
        write(nstats,'(80a1)') cfgname
        if(abs(engunit-9648.530821d0).le.1.d-10) write(nstats,
     x    "(' ENERGY UNITS=electron Volts ')")
        if(abs(engunit-9648530.821d0).le.1.d-10) write(nstats,
     x    "(' ENERGY UNITS=kilo electron Volts ')")
        if(abs(engunit-418.4d0).le.1.d-10)       write(nstats,
     x    "(' ENERGY UNITS=kcal/mol ')")
        if(abs(engunit-1.d2).le.1.d-10)          write(nstats,
     x    "(' ENERGY UNITS=kjoule/mol ')")
        if(abs(engunit-boltz).lt.1.d-10)         write(nstats,
     x    "(' ENERGY UNITS=kelvin ')")
        if(abs(engunit-1.d0).lt.1.d-10)          write(nstats,
     x    "(' ENERGY UNITS=DL_POLY Internal Units ')")
        
      endif
      
c     calculate cell volume and minimum cell half-width
      
      if(imcon.eq.0)then
        
        width=0.d0
        
        stpvol=0.d0
        do i=1,10
           celprp(i)=0.d0
        enddo
        
      else
        
        call dcell(cell,celprp)
        stpvol=celprp(10)
        width=min(celprp(7),celprp(8),celprp(9))/2.d0
        
        if(imcon.eq.4)then

          stpvol=0.5d0*celprp(10)
          width=sqrt(3.d0)*cell(1)/4.d0

        elseif(imcon.eq.5)then
        
          stpvol=0.5d0*celprp(10)
          width=cell(1)/2.d0

        elseif(imcon.eq.6)then

          width=min(celprp(7),celprp(8))/2.d0

        elseif(imcon.eq.7)then
        
          stpvol=0.5d0*celprp(10)

        endif
        
      endif

c     energetic properties of system
      
      stpvir=virsrp+vircpe+virbnd+vircon+vircom+virtbp+virang
     x  +virshl+virtet+virter+virmet+virfld+virrng
c     x  +virshl+virtet+virter+virmet+virfld
      stpeng=engcfg+engke+engrot
      stprot=2.d0*engrot/(boltz*max(1.d0,degrot))
      stpshl=2.d0*shlke/(boltz*max(1.d0,degshl))
      stptmp=2.d0*(engke+engrot)/(boltz*degfre)
      stpprs=0.d0
      if(imcon.gt.0)stpprs=(2.d0*engke-stpvir)/(3.d0*stpvol)
      stpeth=stpeng+stpprs*stpvol
      stpcns=stpeng+consv+engthe

c     convert pressure to units of katm
      
      stpprs=stpprs*prsunt
      
c     calculate mean squared displacements 
c     atomic displacements from origin of production run
      
      if((.not.lzeql).or.(nstep.gt.nsteql))then

        call diffsn0(idnode,natms,mxnode,tstep)
        call diffsn1(idnode,natms,ntpatm,mxnode)
        
      endif

c     zero statistics arrays
      
      if((nstep.le.0).or.(numacc.eq.0))then
        
        numacc=0
        
        do i=1,mxnstk
          
          stpval(i)=0.d0
          sumval(i)=0.d0
          ssqval(i)=0.d0
          
        enddo
        
        do i=1,mxatms
          
          xx0(i)=0.d0
          yy0(i)=0.d0
          zz0(i)=0.d0
          
        enddo
        
      endif
      
c     store current values in statistics array
      
      stpval(1) =stpcns/engunit
      stpval(2) =stptmp
      stpval(3) =engcfg/engunit
      stpval(4) =(engsrp+engmet+engter)/engunit
      stpval(5) =engcpe/engunit
      stpval(6) =engbnd/engunit
      stpval(7) =(engang+engtbp)/engunit
      stpval(8) =(engdih+enginv+engfbp)/engunit
      stpval(9) =engtet/engunit
      stpval(10)=stpeth/engunit
      stpval(11)=stprot
      stpval(12)=stpvir/engunit
      stpval(13)=(virsrp+virmet+virter)/engunit
      stpval(14)=vircpe/engunit
      stpval(15)=virbnd/engunit
      stpval(16)=(virtbp+virang)/engunit
      stpval(17)=vircon/engunit
      stpval(18)=virtet/engunit
      stpval(19)=stpvol
      stpval(20)=stpshl
      stpval(21)=engshl/engunit
      stpval(22)=virshl/engunit
      stpval(23)=acos(celprp(6))*180.d0/pi
      stpval(24)=acos(celprp(5))*180.d0/pi
      stpval(25)=acos(celprp(4))*180.d0/pi
      stpval(26)=virpmf/engunit
      stpval(27)=stpprs

      iadd=27

c     mean squared displacements 
      
      if((.not.lzeql).or.(nstep.gt.nsteql))then
        
        do k=1,ntpatm
          
          stpval(iadd+k)=amsd(k)
          
        enddo
        
      endif

      iadd=iadd+mxatyp

c     stress tensor

      if(abs(stpvol).le.1.d-10) stpvol=1.d0
      
      do i=1,9
        stpval(iadd+i)=stress(i)*prsunt/(stpvol)
      enddo
      
      iadd=iadd+9

c     cell vectors
      
      do i=1,9
        stpval(iadd+i)=cell(i)
      enddo
      
      iadd=iadd+9

c     store pimd variables
      
      if(lpimd)then
        
         stpval(iadd+1)=engqpi/engunit
         stpval(iadd+2)=engqvr/engunit
         stpval(iadd+3)=engrng/engunit
         stpval(iadd+4)=virrng/engunit
         stpval(iadd+5)=qmsrgr
         stpval(iadd+6)=qmsbnd
         stpval(iadd+7)=engthe/engunit
         
      endif
      
      iadd=iadd+7
      
c     check on number of variables for stack - 
      
      if(iadd.gt.mxnstk) call error(idnode,170)

c     accumulate totals over steps
      
      numacc=numacc+1
      sclnv2=1.d0/dble(numacc)
      sclnv1=dble(numacc-1)/dble(numacc)
      
      if(lbpd.and.keybpd.eq.1)then
        
c     calculate true thermodynamic averages in bias potential system
c     note integers numacc and numbpd should be equal in this case
        
        tbold=tboost*dble(numbpd)/dble(numbpd-1)-boost/dble(numbpd-1)
        cterm=0.d0
        do i=1,mxnstk
          
          aterm=sumval(i)*tbold
          bterm=ssqval(i)*tbold**2
          if(tbold.gt.0.d0)cterm=(bterm+aterm**2)/tbold
          ssqval(i)=(sclnv1*(sclnv1*bterm+boost*sclnv2*(cterm+
     x      (tbold*stpval(i)-2.d0*aterm)*stpval(i))))/tboost**2
          sumval(i)=(sclnv1*aterm+boost*sclnv2*stpval(i))/tboost
          
        enddo
        
      else
        
c     calculate true thermodynamic averages in normal system

        do i=1,mxnstk
          
          ssqval(i)=sclnv1*(ssqval(i)+sclnv2*(stpval(i)-sumval(i))**2)
          sumval(i)=sclnv1*sumval(i)+sclnv2*stpval(i)
          
        enddo
        
      endif
      
c     write statistics file
      
      if(idnode.eq.0.and.intsta.gt.0)then

        if(mod(nstep,intsta).eq.0)then
        
          write(nstats,'(i10,1p,e14.6,0p,i10,/,(1p,5e14.6))')
     x       nstep,nstep*tstep,iadd,(stpval(k),k=1,iadd)
          call flush(nstats)
c$$$c     write option for Excel spreadsheet
c$$$          write(nstats,'(i10,1p,e14.6,0p,i10,300(1p,5e14.6))')
c$$$     x      nstep,nstep*tstep,iadd,(stpval(k),k=1,iadd)

        endif
        
      endif
      
c     zero rolling average accumulators
      
      if(nstep.le.0)then
        
        numacc=0
        
        do i=1,mxnstk
          
          zumval(i)=0.d0
          
          do j=1,mxstak
            
            stkval(j,i)=0.d0
            
          enddo
          
        enddo
        
      endif
      
c     store quantities in stack
      
      kstak=mod(nstep-1,nstack)+1
      
      if(nstep.gt.nstack)then
        
        do i=1,mxnstk
          
          zumval(i)=zumval(i)-stkval(kstak,i)
          
        enddo
        
      endif
      
      do i=1,mxnstk
        
        stkval(kstak,i)=stpval(i)
        zumval(i)=zumval(i)+stpval(i)
        
      enddo
      
c     calculate rolling averages
      
      zistk=min(nstack,nstep)
      
      do i=1,mxnstk
        
        ravval(i)=zumval(i)/zistk
        
      enddo
      
c     zero accumulators during equilibration period
      
      if(lzeql.and.nstep.le.nsteql)then
        
        numacc=0
        do i=1,mxnstk
          
          sumval(i)=0.d0
          ssqval(i)=0.d0
          
        enddo
        
      endif

c     close statistics file at regular intervals
      
      if(.not.newjob.and.mod(nstep,ndump).eq.0)then
        
        if(idnode.eq.0)close (nstats)
        newjob=.true.
        
      endif
      
      return
      end subroutine static

      subroutine revive
     x  (lgofr,lzden,lpimd,idnode,imcon,mxnode,natms,nbeads,levcfg,
     x  nstep,nzden,numacc,numrdf,chip,chit,conint,tstep,engcfg,virtot,
     x  vircom,tboost,chit_shl,gaumom)

c***********************************************************************
c     
c     dl_poly subroutine for writing restart files at job termination
c     or at selected intervals in simulation
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith dec 1992.
c     
c***********************************************************************
     
      implicit none
      
      logical lgofr,lzden,lpimd
      integer idnode,imcon,mxnode,natms,nstep,nbeads,nzden,numacc
      integer numrdf,levcfg,nsum,nbuff,i,j,numatm
      real(8) chip,chit,conint,tstep,engcfg,rmxnode,virtot,vircom
      real(8) tboost,chit_shl
      real(8) gaumom(0:5)

      numatm=natms*nbeads
      
      if(mxnode.gt.1)then

c     merge displacement data

        call merge(idnode,mxnode,numatm,mxbuff,xx0,yy0,zz0,buffer)

c     globally sum rdf information before saving
        
        if(lgofr)then

c     maximum rdfs that can be summed in each step
          
          nsum=mxbuff/mxrdf
          if(nsum.eq.0) call error(idnode,200)
          
          nbuff=nsum*mxrdf
          
          do i=1,mxxtyp,nsum
            
            if((mxxtyp+1-i).lt.nsum) nbuff=(mxxtyp+1-i)*mxrdf
            call  gdsum(rdf(1,i),nbuff,buffer)
            
          enddo
          
        endif

c     globally sum zden information before saving
        
        if(lzden)then

c     maximum zdfs that can be summed in each step
          
          nsum=mxbuff/mxzdn
          if(nsum.eq.0) call error(idnode,200)
          
          nbuff=nsum*mxzdn
          
          do i =1,mxatyp,nsum
            
            if((mxatyp+1-i).lt.nsum) nbuff=(mxatyp+1-i)*mxzdn
            call  gdsum(zdens(1,i),nbuff,buffer)
            
          enddo
          
        endif
        
      endif

c     node 0 handles i/o

      if(idnode.eq.0)then

c     write configuration data to new configuration file

        call config_write('REVCON',levcfg,imcon,numatm,engcfg)

c     write accumulator data to dump file
        
        open(nrest,file='REVIVE',form='unformatted')
        
        write(nrest) dble(nstep),dble(numacc),dble(numrdf),chit,
     x    chip,conint,dble(nzden),tboost,chit_shl
        write(nrest) virtot,vircom,eta,strcns,strbod
        write(nrest) stpval
        write(nrest) sumval
        write(nrest) ssqval
        write(nrest) zumval
        write(nrest) ravval
        write(nrest) stkval
        write(nrest) xx0,yy0,zz0
        write(nrest) xxs,yys,zzs
        if(lgofr) write(nrest) rdf
        if(lzden) write(nrest) zdens
        write(nrest) gaumom
        
        close (nrest)

      endif

c     divide rdf data between nodes
      
      rmxnode=1.d0/dble(mxnode)
      
      if(lgofr)then
        
        do i=1,mxxtyp
          
          do j=1,mxrdf
            
            rdf(j,i)=rdf(j,i)*rmxnode
            
          enddo
          
        enddo
        
      endif

c     divide zdensity data between nodes

      if(lzden)then
        
        do i=1,mxatyp
          
          do j=1,mxzdn
            
            zdens(j,i)=zdens(j,i)*rmxnode
            
          enddo
          
        enddo
        
      endif
      
      return
      end subroutine revive

      subroutine zden0(idnode,natms,mxnode,nzden,zlen)

c***********************************************************************
c     
c     dl_poly subroutine for accumulating statistic for density profile
c     zlen=length of cell in z direction
c     
c     double precision accumulators
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1994
c     author    - t. forester    march 1994    
c     
c***********************************************************************

      implicit none

      integer idnode,mxnode,natms,nzden,iatm,ll,k
      real(8) zlen,zleno2,rzdn

c     accumulator

      nzden=nzden+1

c     half of z length

      zleno2=zlen*0.5d0

c     grid interval for density profiles

      rzdn=dble(mxzdn)/zlen

c     set up atom iatm type

      do iatm=idnode+1,natms,mxnode

        k =ltype(iatm) 

        ll=int((zzz(iatm)+zleno2)*rzdn+1.0d0)

c     accumulate statistic

        if(ll.gt.0.and.ll.le.mxzdn)zdens(ll,k)=zdens(ll,k)+1.d0

      enddo
      
      return
      end subroutine zden0

      subroutine zden1
     x  (lpgr,idnode,mxnode,ntpatm,nzden,volm,zlen)

c***********************************************************************
c     
c     dl_poly subroutine for calculating Z density profile
c     from accumulated data.
c     double precision version
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1994
c     author    - t. forester    march 1994    
c     
c***********************************************************************
      
      implicit none

      logical lpgr
      integer idnode,mxnode,ntpatm,nzden,k,j
      real(8) volm,zlen,delzdn,dvolz,factor,sum,rrr,rho

      if(idnode.eq.0) write(nrite,
     x  "(/,/,12X,'Z DENSITY PROFILES',/,/,
     x  'calculated using ',i10,' configurations')") nzden
      
      if(lpgr)then

c     open Z density file and write headers

        if(idnode.eq.0)then

          open(nzdndt,file='ZDNDAT')

          write(nzdndt,'(80a1)')cfgname
          write(nzdndt,'(2i10)')ntpatm,mxzdn

        endif

c     volume of z strip (arbitrary)

        delzdn=zlen/dble(mxzdn)
        dvolz=(volm/zlen)*delzdn
        
c     normalisation factor 
        
        nzden=max(nzden,1)
        factor=1.d0/(dble(nzden)*dvolz)

        do k=1,ntpatm
          
          if(idnode.eq.0)then

             write(nrite,
     x      "(/,'rho(r)  :',a8,/,/,8x,'r',6x,'rho',9x,'n(r)',/)")
     x      unqatm(k)
             write(nzdndt,'(a8)')unqatm(k)

          endif

c     global sum of data on all nodes
          
          if(mxnode.gt.1)call gdsum(zdens(1,k),mxzdn,buffer)

c     running integration of z-density
          
          sum=0.d0

c     loop over distances
          
          do j=1,mxzdn
            
            rrr=(dble(j)-0.5d0)*delzdn-zlen*0.5d0
            rho=zdens(j,k)*factor
            sum=sum+rho*dvolz

c     print out information
            
            if(idnode.eq.0)then

              write(nrite,"(f10.4,1p,2e14.6)") rrr,rho,sum
              write(nzdndt,"(1p,2e14.6)") rrr,rho

            endif
            
          enddo
          
        enddo
        
        if(idnode.eq.0)close (nzdndt)

      endif
      
      return
      end subroutine zden1

      subroutine rdf0neu(ik,rcut)

c***********************************************************************
c     
c     dl_poly subroutine for accumulating statistic for radial
c     distribution functions.
c     neutral group implementation
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1994
c     author    - t. forester    march 1994    
c     amended     t. forester    april 1994
c     
c***********************************************************************
      
      implicit none

      integer ik,m,iatm,jatm,ll,k
      real(8) rcut,a0,a1,a2,a3,a4,a5,rcsq,rrcsq,sqlim,rdelr
      real(8) ai,aj,rsq,rrr,sss
      
      data a0,a1,a2,a3,a4,a5/.0837557783d0,2.9399054d0,-7.8475201d0,
     x  14.1328992d0,-12.6228528d0,4.32084948d0/
      
c     set cutoff condition for pair forces
      
      rcsq=rcut**2
      rrcsq=1.d0/rcsq
      sqlim=0.01d0*rcsq
      
c     grid interval for rdf tables
      
      rdelr=dble(mxrdf)/rcut

c     start of primary loop for rdf accumulation
      
      do m=1,ik
        
c     atomic and potential function indices
        
        iatm=ilist(m)
        ai=ltype(iatm)
        
        jatm=jlist(m) 
        aj=ltype(jatm)
        
        if(ai.gt.aj)then
          ll=int(ai*(ai-1.d0)*0.5d0+aj+0.5d0)
          k=lstvdw(ll)
        else
          ll=int(aj*(aj-1.d0)*0.5d0+ai+0.5d0)
          k=lstvdw(ll)
        endif
        
        rsq=rsqdf(m)
        
        if(rcsq.gt.rsq)then
          
c     determine interpolation panel for rdf table
          
          if(rsq.lt.sqlim)then
            
            rrr=sqrt(rsq)
            
          else

c     interpolate square-root by polynomial plus newton-raphson
            
            sss=rsq*rrcsq
            rrr=1.d0/
     x        (a0 +sss*(a1+sss*(a2+sss*(a3+sss*(a4+sss*a5)))))
            rrr=0.5d0*rrr*(3.d0-sss*rrr*rrr)
            rrr=0.5d0*rrr*(3.d0-sss*rrr*rrr)
            rrr=0.5d0*rrr*(3.d0-sss*rrr*rrr)*sss*rcut
            
          endif
          
          ll=int(rrr*rdelr+0.999999d0)

c     accumulate statistics
          
          rdf(ll,k)=rdf(ll,k)+1.d0
          
        endif
        
      enddo
      
      return
      end subroutine rdf0neu
      
      end module property_module
