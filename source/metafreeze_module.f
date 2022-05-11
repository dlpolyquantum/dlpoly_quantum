      module metafreeze_module

c---------------------------------------------------------------------
c     
c     Metafreeze module for metadynamics
c     Author D. Quigley - University of Warwick
c     Copyright D. Quigley
c     Adapted for dl_poly classic w. smith dec 2010
c
c---------------------------------------------------------------------
      
      implicit none
      
c     Minimise usage by making everything private
      
      private

c---------------------------------------------------------------------
c     P u b l i c   R o u t i n e s
c     ...  unless exposed here.
c---------------------------------------------------------------------
      
      public :: define_metadynamics ! Called to intitialise this module
      public :: metafreeze_driver ! Called at every MD time-step

c---------------------------------------------------------------------
c     P u b l i c   V a r i a b l e s
c---------------------------------------------------------------------
      
      public :: lmetadyn  ! Is this a metadynamics calculation?
      public :: driven    ! Is this atom site involved in metadynamics
      public :: lstein,ltet,lglobpe,llocpe,ncolvar,nq4,nq6,ntet,hkey
      public :: meta_step_int,globpe_scale,locpe_scale,ref_W_aug
      public :: h_aug,wt_Dt

c     Local copies of exluded atom arrays note that these are
c     not indexed by atom number, but by an index which loops
c     over the atoms handles by each MPI rank.
      public :: mtd_nexatm,mtd_lexatm

c------------------------------------
c     Populated from CONTROL file   
c------------------------------------
      
      logical,save :: lmetadyn =.false. ! Master metadynamics flag
      
c-------------------------------------------------------
c     Energy, virial, stress and forces from 'local' pe
c     Populated by dlpoly energy routines  
c-------------------------------------------------------
      
      
      public :: eng_loc,vir_loc,stress_loc  
      public :: fxx_loc,fyy_loc,fzz_loc

c----------------------------------------------------------------------
c     Data accumulated for local potential energy (computed elsewhere) 
c----------------------------------------------------------------------
      
      real(8),save                           :: eng_loc,vir_loc
      real(8),dimension(9),save              :: stress_loc
      real(8),allocatable,dimension(:),save  :: fxx_loc,fyy_loc,fzz_loc
      
c----------------------------------------------------------------------
c     Arrays holding information on excluded interactions. Replicated 
c     here to avoid a compilation dependency loop which would occur
c     if simply using the arrays already in exclude_module.f
c----------------------------------------------------------------------
      integer,allocatable,dimension(:),save   :: mtd_nexatm
      integer,allocatable,dimension(:,:),save :: mtd_lexatm

c---------------------------------------------------------------------
c     P r i v a t e   V a r i a b l e s 
c---------------------------------------------------------------------

c----------------------------------------------------------------
c     Collective variables and derivatives of v_aug w.r.t colvars
c----------------------------------------------------------------

      integer,parameter :: maxhis = 150000 ! Size of history arrays      
      real(8),allocatable,dimension(:),save   :: colvar,dcolvar
      real(8),allocatable,dimension(:),save   :: colvar_scale

c----------------------------------------------------------------
c     Positions and heights of previous Gaussians in colvar space
c----------------------------------------------------------------
      
      real(8),allocatable,dimension(:,:),save :: colvar_his
      real(8),allocatable,dimension(:),save   :: w_aug
      
c------------------------------
c     Read from CONTROL file
c------------------------------
      
      integer,save :: ncolvar  = 0 ! Total number of collvars
      logical,save :: lstein   =.false. ! Q4/Q6 collective variables
      logical,save :: ltet     =.false. ! Tetrahedral order parameter
      logical,save :: lglobpe  =.false. ! Global potential energy
      logical,save :: llocpe   =.false. ! Local potential energy
      
      integer,save :: nq4   = 0 ! Number of Q4 pair types
      integer,save :: nq6   = 0 ! Number of Q6 pair types
      integer,save :: ntet  = 0 ! Number of zeta triplets
      
      real(8),save :: globpe_scale = 1.0d0 ! Scaling factors for local 
      real(8),save :: locpe_scale  = 1.0d0 ! and global pe colvars
      
      real(8),save :: ref_W_aug=1.0d0 ! Reference Gaussian height
      real(8),save :: h_aug=1.0d0     ! Gaussian width
      integer,save :: hkey=0          ! Height control scheme
      real(8),save :: wt_Dt=100.0d0   ! "Well-tempered" parameter
      integer,save :: meta_step_int=5 ! interval between depositions
      
c----------------------------------------
c     Read from STEINHARDT or TETRAHEDRAL
c----------------------------------------
      
c     Global Steinhardt order parameters
      
      real(8),allocatable,dimension(:),save   :: q4_global 
      real(8),allocatable,dimension(:),save   :: q6_global 
      
c     Global Tetrahedral order parameters
      
      real(8),allocatable,dimension(:),save   :: zeta_global 
      
c     Bookkeeping arrays for order parameter computation
      
      character(8),allocatable,dimension(:,:),save  :: q4label
      character(8),allocatable,dimension(:,:),save  :: q6label
      character(8),allocatable,dimension(:),save    :: zetalabel
      
c     Inner and outer cutoffs
      
      real(8),allocatable,dimension(:,:),save :: q4cutoff
      real(8),allocatable,dimension(:,:),save :: q6cutoff
      real(8),allocatable,dimension(:,:),save :: zetacutoff
      
c     Scaling factors for q4 and q6
      
      real(8),allocatable,dimension(:),save   :: q4scale
      real(8),allocatable,dimension(:),save   :: q6scale
      real(8),allocatable,dimension(:),save   :: zetascale
      
c     Number of nearest neighbours for q4, q6 and zeta
      
      integer,allocatable,dimension(:),save   :: q4nn,q6nn
      integer,allocatable,dimension(:),save   :: zetann

c------------------------------------------------------------
c     Arrays holding data for computation of order parameters
c------------------------------------------------------------
      
c     Steinhardt site-site interaction arrays
      
      integer,allocatable,dimension(:,:),save  :: q4site
      integer,allocatable,dimension(:,:),save  :: q6site
      integer,allocatable,dimension(:)  ,save  :: zetasite
      
c     Number of included sites
      
      integer,allocatable,dimension(:),save    :: q4ninc
      integer,allocatable,dimension(:),save    :: q6ninc
      integer,allocatable,dimension(:),save    :: zetaninc
      
c     Real and imaginary parts of q4bar and q6bar
      
      real(8),allocatable,dimension(:,:),save :: ReQ6bar,ImQ6bar
      real(8),allocatable,dimension(:,:),save :: ReQ4bar,ImQ4bar
      
c     Max number of entries in co-ordination shell
      
      integer,parameter :: mxflist = 50 
      integer :: mxninc
      
c     Full neighbour list for Tetrahedral order parameter
      
      integer,allocatable,dimension(:)  ,save   :: nflist 
      integer,allocatable,dimension(:,:),save   :: flist
      
c-------------------------------------
c     Internal bookkeeping
c-------------------------------------
      
      logical,allocatable,dimension(:),save :: driven ! Metadynamics option
      integer,save :: meta_step=1 ! Current metadynamics step number
      real(8),save :: meta_energy ! Value of metadynamics bias potential
      
      integer,save :: wl_nbins=30    ! Number of bins for WL recursion
      integer,save :: wl_cycle=0     ! Current WL cycle
      real(8),save :: wl_range=0.175 ! range of WL 
      real(8),allocatable,dimension(:),save :: wl_bin ! WL bins

c--------------------------------------
c     Miscellaneous internal variables
c--------------------------------------

      integer,allocatable,dimension(:) :: buff ! Comms buffer

c     File units
      
      integer,save :: stn = 91  ! STEINHARDT
      integer,save :: mtd = 92  ! METADYNAMICS
      integer,save :: zta = 93  ! ZETA
      integer,save :: wlb = 94  ! WL_BINS.DAT
      
c     Error flag
      
      integer,dimension(100) :: ierr = 0
      
c     Local store of comms variables
c     Assuming no task farming, comms will require changing if farmed
      
      integer, save :: myrank,commsize 
      logical, save :: onroot
      real(8),save :: kt

      contains

      Subroutine Metafreeze_Driver
     x  (imcon,natms,temp,nstep,engcfg,virtot,engord,virord)
      
c---------------------------------------------------------------------
c     Top level metadynamics routine called after evaluation of all 
c     other energetic and force terms within the main molecular 
c     dynamics loop. 
c     
c     1. Computes the ncolvar order parameters
c     2. Deposits a new Gaussian at the current collective variables 
c        as the current number of steps reaches meta_step_int
c     3. Computed the bias potential and its derivative w.r.t. the 
c        ncolvar collective variables.
c     4. Computes the forces stresses and virial resulting from the 
c        bias
c     
c     Author  D. Quigley - University of Warwick
c     Copyright D. Quigley
c
c---------------------------------------------------------------------

      use setup_module,    only : boltz
      use config_module,   only : fxx,fyy,fzz,stress
      
      implicit none
      
      integer,intent(in)    :: nstep,imcon,natms
      real(8),intent(in)    :: engcfg,virtot,temp
      real(8),intent(out)   :: engord,virord
      
c     Local variables
      
      integer       :: k,iq,itet,ibin,nfail,my_meta_step
      integer,save  :: nlastg = 0
      real(8)       :: height,buff1,wl_mean
      logical       :: flat,safe
      
c------------------------------------------------------
c     Compute order parameters / collective variables
c------------------------------------------------------

c     Steinhardt order parameters
      
      if ( nq4>0.or.nq6>0 ) call compute_steinhardt(imcon,natms)
      
      k = 1
      do iq = 1,nq4
        colvar(k) = q4_global(iq)
        k = k + 1
      end do
      do iq = 1,nq6
        colvar(k) = q6_global(iq)
        k = k + 1
      end do
      
c     Tetrahedral order parameters
      
      if ( ntet > 0 ) then
        call compute_tet_nlist(imcon,natms)
        call compute_tetrahedral(imcon,natms)
      end if  
      
      do itet = 1,ntet
        colvar(k) = zeta_global(itet)
        k = k + 1
      end do
      
c     Energy order parameters
      
      if (lglobpe) then
        colvar(k) = engcfg
        k = k + 1
      end if
      if (llocpe) then
        
c     Global reduction of local virial and energy
        
        if ( commsize > 1 ) call gdsum(eng_loc,1,buff1)
        if ( commsize > 1 ) call gdsum(vir_loc,1,buff1)
        colvar(k)   = eng_loc
        k = k + 1
      end if
      
      if ( k-1/=ncolvar ) call Mfrz_Error(2500,0.d0) 
      
      if ( hkey==1 ) then
        k = int(dble(wl_nbins)*colvar(1)/wl_range) + 1
        if ( k < wl_nbins) wl_bin(k) = wl_bin(k) + 1.0d0
      end if
      
c--------------------------------------------------------
c     Deposit a new Gaussian if now is the correct time
c--------------------------------------------------------
      
      if ( (mod(nstep,meta_step_int)==0).and.(nstep>nlastg) ) then
        nlastg = nstep          ! Avoid multiple depositions at the
                                ! same timestep (relaxed shell model)
        
        select case (hkey)
          
        case(0)
          
c     Always deposit Gaussians of the same height
          
          height = ref_W_aug
          
        case(1)
          
c     Wang-Landau style recursion
          
          open(unit=wlb,file='WL_BINS.DAT',status='replace')
          
          do ibin = 1,wl_nbins
            write(wlb,*)ibin,wl_bin(ibin)
          end do
          
          close(wlb)
          
          if ( ncolvar/=1 ) then
            call Mfrz_Error(2501,0.d0)
          else
            
            height    = ref_W_aug*(0.5d0**dble(wl_cycle))
            
            nfail = 0
            wl_mean = 0.d0
            do ibin = 6,wl_nbins-5
              wl_mean = wl_mean + wl_bin(ibin)
              nfail = nfail + 1
            end do
            wl_mean = wl_mean/dble(nfail)
            
            nfail = 0
            flat = .true.
            do ibin = 6,wl_nbins-5
              if ( wl_bin(ibin) < 0.8d0*wl_mean ) then
                if ( nfail > 2 ) flat = .false.
                nfail = nfail + 1
              end if
            end do
            
            if ( flat.and.(sum(wl_bin)>50.0d0) ) then
              wl_cycle = wl_cycle + 1
              wl_bin   = 0.0d0
            end if
            
            height = ref_W_aug*(0.5d0**dble(wl_cycle))
            
          end if   
          
        case(2)
          
c     Well-tempered metadynamics
          
          meta_energy = 0.0d0
          call compute_bias_potential()
          
          height = ref_W_aug*exp(-meta_energy/wt_Dt)
          
        case default
          
          call Mfrz_Error(2502,0.d0)
          
        end select
        
        call deposit_gaussian(height,temp)
        my_meta_step = (meta_step-1)/commsize + 1
        safe = ( maxhis >= my_meta_step )
        call gstate(safe)
        if ( .not.safe ) call Mfrz_Error(2503,0.d0)
        
      end if
      
c-----------------------------------------------------------
c     Compute the bias potential and its derivatives w.r.t.
c     to the ncolvar collective variables.                 
c-----------------------------------------------------------
      
      call compute_bias_potential()

c-----------------------------------------------------------
c     Add in the forces, stresses and virial contributions
c     from this derivative.
c-----------------------------------------------------------
      
      virord = 0.0d0            ! Zero the virial
      
c     Must compute contributions from pe order parameters
c     first before we change any forces.
      
      k = nq4+nq6+ntet+1
      
c     Energy order parameters
      
      if (lglobpe) then
        
        fxx(:) = fxx(:)*(1.0d0+dcolvar(k))
        fyy(:) = fyy(:)*(1.0d0+dcolvar(k))  
        fzz(:) = fzz(:)*(1.0d0+dcolvar(k))
        
c     correct for later summation:
        
        virord = virord+dcolvar(k)*virtot/dble(commsize)
        stress = stress*(1.0d0+dcolvar(k))
        
      end if
      if (llocpe) then
        
        fxx(:) = fxx(:) + fxx_loc(:)*dcolvar(k)
        fyy(:) = fyy(:) + fyy_loc(:)*dcolvar(k)
        fzz(:) = fzz(:) + fzz_loc(:)*dcolvar(k)
        
c     correct for later summation:
        
        virord = virord + dcolvar(k)*vir_loc/dble(commsize)
        stress = stress + stress_loc*dcolvar(k)
        
      end if
      
c     Steinhardt order parameters
      
      if ( nq4>0.or.nq6>0 ) call 
     x  compute_steinhardt_forces(imcon,natms,engord,virord)
      
c     Tetrahedral order parameters
      
      if ( ntet > 0 ) call
     x  compute_tetrahedral_forces(imcon,natms,engord,virord)
      
c     global reduction of virord

      if ( commsize > 1 ) call gdsum(virord,1,buff1)
      
      engord = meta_energy

c     write(0,'("DEBUG : engord = ",F12.6)')engord/(temp*boltz)
      
      return
      
      end Subroutine Metafreeze_Driver
      
      Subroutine Deposit_Gaussian(height,temp)

c---------------------------------------------------------------------
c     
c     Deposits a new Gaussian at the current collective variables and 
c     appends to the METADYNAMICs file.
c     
c     Author D. Quigley - University of Warwick
c     Copyright D. Quigley
c
c---------------------------------------------------------------------
      
      use setup_module, only : boltz
      
      implicit none
      
      real(8),intent(in) :: height,temp
      integer       :: my_meta_step
      character(11) :: fmtstring
      
c     store current order parameters and Gaussian height
      
      if ( mod(meta_step-1,commsize) == myrank ) then
        
        my_meta_step = (meta_step-1)/commsize + 1
        w_aug(my_meta_step)        = height
        colvar_his(:,my_meta_step) = colvar(:)
        
      end if
      
      if (onroot) then
        
c     Create format string
        
        write(fmtstring,'("(I8,",I1,"E15.6)")')ncolvar+1
        
c     write METADYNAMICS file
        
        open(unit=mtd,file='METADYNAMICS',status='old',position=
     x    'append',iostat=ierr(1))
        write(unit=mtd,fmt=fmtstring)meta_step,colvar(:),
     x    height/(temp*boltz)
        close(unit=mtd)
        
      end if
      
      meta_step = meta_step+1
      
      return
      
      end Subroutine Deposit_Gaussian

      Subroutine Compute_Bias_Potential()
      
c---------------------------------------------------------------------
c
c     Computes the augmenting bias potential as a function of the 
c     collective variables. Also computes the derivative of the bias 
c     potential w.r.t. the collective variables required to compute 
c     the metadynamics forces.
c
c     Author D. Quigley - University of Warwick
c     Copyright D. Quigley
c     
c---------------------------------------------------------------------
      
      implicit none
      
      integer :: istep,ic,k,my_meta_step
      real(8) :: vsq,exp1,dWpsq
      real(8),allocatable,dimension(:) :: buff1,buff2
      
      allocate(buff1(1:2*(ncolvar+1)),stat=ierr(1))
      allocate(buff2(1:2*(ncolvar+1)),stat=ierr(2))
      
      if (any(ierr/=0)) call Mfrz_Error(2504,0.d0)

c     Set squared-width of gaussians
      
      dWpsq = 1.0d0/h_aug**2
      meta_energy  = 0.0d0
      
c     Zero accumulators of derivative w.r.t. each order parameter
      
      dcolvar(:) = 0.0d0
      my_meta_step = (meta_step-1)/commsize + 1
      do istep=1,my_meta_step
        
        vsq = 0.0d0
        do ic = 1,ncolvar
          vsq = vsq + ( colvar_scale(ic)*(colvar(ic) - 
     x          colvar_his(ic,istep)) )**2
        end do
        exp1 = w_aug(istep)*exp(-0.5d0*vsq*dWpsq)
        do ic = 1,ncolvar
          dcolvar(ic) = dcolvar(ic) - (colvar_scale(ic)**2)*exp1*
     x                  (colvar(ic) - colvar_his(ic,istep))*dWpsq
        end do
        
        meta_energy = meta_energy + exp1
        
      end do
      
      buff1(1) = meta_energy
      k = 2
      do ic = 1,ncolvar
        buff1(k) = dcolvar(ic)
        k = k + 1
      end do
      
      if ( commsize > 1 ) then
        call gdsum(buff1,ncolvar+1,buff2)
      end if
      
      meta_energy = buff1(1)
      
      k = 2
      do ic=1,ncolvar
        dcolvar(ic) = buff1(k)
        k = k + 1
      end do

c     write(0,'("DEBUG : CV derivs = ",6F15.6)')dcolvar

      deallocate(buff1,buff2,stat=ierr(1))
      
      return
      
      end Subroutine Compute_Bias_Potential
      
      Subroutine Define_Metadynamics(tm,ts,natms,ntpatm,temp)
      
c---------------------------------------------------------------------
c     Processes the metadynamics input file. This is done in several 
c     stages.
c     1. Process the metadynamics control data read from the CONTROL
c        file, which defines the number of collective variables and
c        indicates if we need to read from auxilliary input files 
c        (e.g. STEINHARDT or ZETA) which define order parameters.
c     2. Read and process these auxilliary files.
c     3. Process the information obtained from the CONTROL file which
c        controls the properties of the Gaussians used to build the 
c        bias potential.
c     
c     Author D. Quigley - University of Warwick
c     Copyright D. Quigley
c     Adapted w. smith  - jan 2011
c     
c---------------------------------------------------------------------
      
      use site_module
      use config_module, only : ltype
      use parse_module
      
      implicit none
      
      integer,intent(in) :: tm,ts,natms,ntpatm
      real(8),intent(in) :: temp
      
c     Local variables
      
      integer :: isite,jsite,ilin,i,iq,iatm0,iatm1,istrd,iatm,k,n,idum
      logical :: lexist,go,safe
      real(8) :: waug,dummy(1)
      
c     Allocate the driven list. Modifications elsewhere in DLPOLY always
c     check if sites are driven (for accumulation of local PE) and hence
c     this should be allocated if this is a metadynamics run or not.
      
      allocate(driven(1:size(unqatm)),stat=ierr(1))
      if (ierr(1)/=0) call Mfrz_Error(2505,0.d0)
      driven = .false.

c     Do nothing else if the metafreeze flag has not been set in CONTROL
      
      if (.not.lmetadyn) then
        return
      end if
      
      myrank=tm 
      commsize = ts
      onroot = (myrank==0)

    
c$$$     DEBUG
cc$$$    if (onroot) write(0,'("================================")')
cc$$$    if (onroot) write(0,'("Available sites from site module")')
cc$$$    if (onroot) write(0,'("================================")')
cc$$$    if (onroot) then
cc$$$       do isite = 1,size(unqatm)
cc$$$          write(0,'("Site index ",i4,": ",a8)')isite,unqatm(isite)
cc$$$       end do
cc$$$    end if
    
c     Cannot bias the global and local PE in the same run.
      
      if ( lglobpe.and.llocpe ) call Mfrz_Error(2509,0.d0)

c     Allocate local force arrays if anything is driven.
      
      allocate(fxx_loc(1:mxatms),stat=ierr(1))
      allocate(fyy_loc(1:mxatms),stat=ierr(2))
      allocate(fzz_loc(1:mxatms),stat=ierr(3))
      if (any(ierr/=0)) call Mfrz_Error(2510,0.d0)
   
c     Allocate arrays to hold collective variables 
      
      allocate( colvar(1:ncolvar),stat=ierr(1))
      allocate(dcolvar(1:ncolvar),stat=ierr(2))
      allocate(colvar_his(1:ncolvar,1:maxhis),stat=ierr(3))
      allocate(colvar_scale(1:ncolvar),stat=ierr(4))
      allocate(w_aug(1:maxhis),stat=ierr(5))
      if (any(ierr/=0)) call Mfrz_Error(2511,0.d0)
      
c     Allocate Wang-Landau bin array

      allocate(wl_bin(1:wl_nbins),stat=ierr(1))
      if (any(ierr/=0)) call Mfrz_Error(2512,0.d0)
      wl_bin = 0.0d0
      
c-------------------------------------------------------------
c     Process Steinhardt order parameter settings if present
c-------------------------------------------------------------
      
      if ( lstein ) then
        
        allocate(q4label(1:2,1:nq4),stat=ierr(1))
        allocate(q6label(1:2,1:nq6),stat=ierr(2))
        allocate(q4cutoff(1:2,1:nq4),stat=ierr(3))
        allocate(q6cutoff(1:2,1:nq6),stat=ierr(4))
        allocate(q4scale(1:nq4),stat=ierr(5))
        allocate(q6scale(1:nq6),stat=ierr(6))
        allocate(q4nn(1:nq4),stat=ierr(7))
        allocate(q6nn(1:nq6),stat=ierr(8))
        allocate(q4ninc(1:nq4),stat=ierr(9))
        allocate(q6ninc(1:nq6),stat=ierr(10))
        allocate(buff(1:max(nq4,nq6)),stat=ierr(11))
        allocate(q4_global(1:nq4),stat=ierr(12))
        allocate(q6_global(1:nq6),stat=ierr(13))
        allocate(ReQ4Bar(-4:+4,1:nq4),stat=ierr(14))
        allocate(ImQ4Bar(-4:+4,1:nq4),stat=ierr(15))
        allocate(ReQ6Bar(-6:+6,1:nq6),stat=ierr(16))
        allocate(ImQ6Bar(-6:+6,1:nq6),stat=ierr(17))       
        if (any(ierr/=0)) call Mfrz_Error(2515,0.d0)
        
c     Open STEINHARDT file and process
        
        if (onroot) then
          open(unit=stn,file='STEINHARDT',status='old',iostat=ierr(1))
        else
          ierr(1)=0
        endif
        call gisum(ierr(1),1,ierr(2))
        if ( ierr(1)/=0 ) call Mfrz_Error(2516,0.d0)
        
        ilin = 1
        safe=.true.
        if (nq4>0) then
          call getrec(safe,myrank,stn) ! Ignore q4 comment line
          ilin = ilin + 1
          do i = 1,nq4
            call getrec(safe,myrank,stn)
            if (safe) then
              call getword(q4label(1,i),record,8,lenrec)
              call getword(q4label(2,i),record,8,lenrec)
              q4cutoff(1,i)=dblstr(record,lenrec,idum)
              q4cutoff(2,i)=dblstr(record,lenrec,idum)
              q4scale(i)=dblstr(record,lenrec,idum)
              q4nn(i)=intstr(record,lenrec,idum)
              ierr(ilin)=0
            else
              ierr(ilin)=1
            endif
            ilin = ilin + 1
          end do
        end if
        if (nq6>0) then
          call getrec(safe,myrank,stn) ! Ignore q6 comment line
          ilin = ilin + 1
          do i = 1,nq6
            call getrec(safe,myrank,stn)
            if (safe) then
              call getword(q6label(1,i),record,8,lenrec)
              call getword(q6label(2,i),record,8,lenrec)
              q6cutoff(1,i)=dblstr(record,lenrec,idum)
              q6cutoff(2,i)=dblstr(record,lenrec,idum)
              q6scale(i)=dblstr(record,lenrec,idum)
              q6nn(i)=intstr(record,lenrec,idum)
              ierr(ilin)=0
            else
              ierr(ilin)=1
            endif
            ilin = ilin + 1
          end do
        end if
        if (onroot) close(unit=stn)
        
        call gisum(ierr(1),ilin-1,ierr(ilin))
        do i = 1,ilin-1
          if (ierr(i)/=0) then
            call Mfrz_Error(2521,dble(i))
          end if
        end do
        
c     Create array indicating which site-site connections use
c     which set of q4 cut-offs, scaling factors and num neighbours.
       
       allocate(q4site(1:size(unqatm),1:size(unqatm)),stat=ierr(1))
       if (ierr(1)/=0) call Mfrz_Error(2517,0.d0)
       q4site(:,:) = 0
       
       do isite = 1,ntpatm
         do jsite = isite,ntpatm
           do iq = 1,nq4
             if ((q4label(1,iq)==unqatm(isite)).and.
     x         (q4label(2,iq)==unqatm(jsite))) then
               q4site(jsite,isite) = iq
               q4site(isite,jsite) = iq
               driven(jsite) = .true.
               driven(isite) = .true.
             end if
           end do
         end do
       end do

       allocate(q6site(1:size(unqatm),1:size(unqatm)),stat=ierr(1))
       if (ierr(1)/=0) call Mfrz_Error(2518,0.d0)
       q6site(:,:) = 0
       do isite = 1,ntpatm
         do jsite = isite,ntpatm
           do iq = 1,nq6
             if ((q6label(1,iq)==unqatm(isite)).and.
     x         q6label(2,iq)==unqatm(jsite)) then
               q6site(jsite,isite) = iq
               q6site(isite,jsite) = iq
               driven(jsite) = .true.
               driven(isite) = .true.
             end if
           end do
         end do
       end do
       
c     Count number of included sites
       
       iatm0 = myrank+1
       iatm1 = natms
       istrd = commsize
       
       q4ninc = 0 
       q6ninc = 0
       do iatm = iatm0,iatm1,istrd
         
         isite = ltype(iatm)
         
         do iq = 1,nq4
           if (unqatm(isite)==q4label(1,iq)) q4ninc(iq) = q4ninc(iq) + 1
         end do
         do iq = 1,nq6
           if (unqatm(isite)==q6label(1,iq)) q6ninc(iq) = q6ninc(iq) + 1
         end do
         
       end do
       
       if ( commsize > 0 ) then
         if (nq4>0) call gisum(q4ninc,nq4,buff(1:nq4))
         if (nq6>0) call gisum(q6ninc,nq6,buff(1:nq6))
       end if
       
       deallocate(buff,stat=ierr(1))
       if (ierr(1)/=0) call Mfrz_Error(2519,0.d0)
       
      end if                    ! end if steinhardt order parameters
      
      if ( ltet ) then
        
        allocate(zetacutoff(1:2,1:ntet),stat=ierr(1))
        allocate(zeta_global(1:ntet),stat=ierr(2))
        allocate(zetascale(1:ntet),stat=ierr(3))
        allocate(zetalabel(1:ntet),stat=ierr(4))
        allocate(zetann(1:ntet),stat=ierr(5))
        allocate(zetaninc(1:ntet),stat=ierr(6))
        allocate(buff(1:ntet),stat=ierr(7))
        if (any(ierr/=0)) call Mfrz_Error(2522,0.d0)
        
c     Open ZETA file and process

        if (onroot) then
          open(unit=zta,file='ZETA',status='old',iostat=ierr(1))
        else
          ierr(1)=0
        endif
        call gisum(ierr(1),1,ierr(2))
        if ( ierr(1)/=0 ) call Mfrz_Error(2523,0.d0)
        
        ilin = 1
        safe=.true.
        if (ntet>0) then
          call getrec(safe,myrank,zta) ! Ignore comment line
          ilin = ilin + 1
          do i = 1,ntet
            call getrec(safe,myrank,zta)
            if (safe) then
              call getword(zetalabel(i),record,8,lenrec)
              zetacutoff(1,i)=dblstr(record,lenrec,idum)
              zetacutoff(2,i)=dblstr(record,lenrec,idum)
              zetascale(i)=dblstr(record,lenrec,idum)
              zetann(i)=intstr(record,lenrec,idum)
              ierr(ilin)=0
            else
              ierr(ilin)=1
            endif
            ilin = ilin + 1
          end do
        end if
        if (onroot) close(unit=zta)

        call gisum(ierr(1),ilin-1,ierr(ilin))
        do i = 1,ilin-1
          if (ierr(i)/=0) then
            call Mfrz_Error(2529,dble(i))
          end if
        end do
        
c     Create array indicating which site-site connections use
c     which set of q4 cut-offs, scaling factors and num neighbours.
        
        allocate(zetasite(1:size(unqatm)),stat=ierr(1))
        if (ierr(1)/=0) call Mfrz_Error(2524,0.d0)
        zetasite(:) = 0
        
        do isite = 1,size(unqatm)
          do iq = 1,ntet
            if (zetalabel(iq)==unqatm(isite)) then
              zetasite(isite) = iq
              driven(isite) = .true.
            end if
          end do
        end do
        
c     Count number of included sites
        
        iatm0 = myrank+1
        iatm1 = natms
        istrd = commsize
        
        zetaninc(:) = 0
        do iatm = iatm0,iatm1,istrd
          
          isite = ltype(iatm)
          
          do iq = 1,ntet
            if (unqatm(isite)==zetalabel(iq)) 
     x        zetaninc(iq) = zetaninc(iq) + 1
          end do
        end do
        
        if (commsize>1) then
          if (ntet>0) call gisum(zetaninc,ntet,buff)
        end if
        
c$$$        do iq = 1,ntet
c$$$          write(0,'("Number of sites for zeta type ",I5," : ",I5)')
c$$$     x      iq,zetaninc(iq)
c$$$        end do
        
        mxninc = max(100,4*maxval(zetaninc)/commsize)
        allocate(nflist(1:mxninc),stat=ierr(1))
        allocate(flist(1:mxflist,1:mxninc),stat=ierr(2))
        if (any(ierr/=0)) call Mfrz_Error(2525,0.d0)
        
        deallocate(buff,stat=ierr(1))
        if (ierr(1)/=0) call Mfrz_Error(2519,0.d0)
        
      end if                    ! end if tetrahedral order parameters
      
c     Check total number of collective variables (ncolvar) matches total
c     number specified by nq4, nq6, ntet and potential energy flags.
      
      k = 0
      if (llocpe  ) k = k + 1
      if (lglobpe ) k = k + 1
      k = k + ntet + nq4 + nq6
      if ( k /= ncolvar ) call Mfrz_Error(2527,0.d0)
    
c     populate colvar_scale
      
      k = 1
      do iq = 1,nq4
        colvar_scale(k) = q4scale(iq)
        k = k + 1
      end do
      do iq = 1,nq6
        colvar_scale(k) = q6scale(iq)
        k = k + 1
      end do
      do iq = 1,ntet
        colvar_scale(k) = zetascale(iq)
        k = k + 1
      end do
      if (lglobpe) then
        colvar_scale(k) = globpe_scale
        k = k + 1
      end if
      if (llocpe) then
        colvar_scale(k) = locpe_scale
        k = k + 1
      end if
      
c     write(0,*)lglobpe,llocpe
c     write(0,'("DEBUG : CV Scaling factors : ",6F15.6)')colvar_scale(:)


c     Convert into internal units
      
      wt_Dt = wt_Dt*temp*boltz
      ref_W_aug = ref_W_aug*temp*boltz
      kt        = temp*boltz
      
c---------------------------------------------------------------------
c     Purge the METADYNAMICS file or re-open and read if this is a 
c     restart. N.B. we assume a restart if REVOLD is present and 
c     ignore keyres.
c---------------------------------------------------------------------
      
      if (onroot) then
        inquire(file='REVOLD',exist=lexist)
      else
        lexist=.true.
      endif
      call gstate(lexist)
      
      if (lexist) then
        
c     read contents of METADYNAMICS file
        
        if (onroot) then
          
          open(unit=mtd,file='METADYNAMICS',status='old',iostat=ierr(1))
          
          k = 0
          do
            read(unit=mtd,fmt=*,end=10)meta_step,colvar(:),waug
            waug = waug*temp*boltz
            if (k == 0) then
              n = (meta_step-1)/commsize + 1
              colvar_his(:,n)=colvar(:)
              w_aug(n)=waug
            else
              dummy(1)=dble(meta_step)
              call csend(17947,dummy,1,k,idum)
              call csend(17948,colvar,ncolvar,k,ierr(3))
              dummy(1)=waug
              call csend(17949,dummy,1,k,ierr(4))
            end if
            
            k = k + 1
            if (k == commsize) k = 0
          end do
          
   10     close(unit=mtd)
          
          do k=1,commsize-1
            dummy(1)=-dble(meta_step)
            call csend(17947,dummy,1,k,ierr(2))
          end do
          
        else
          
          go = .true.
          do while(go)
            
            call crecv(17947,dummy,1)
            meta_step=nint(dummy(1))
            ierr(2)=0
            
            if ( meta_step < 0 ) then
              meta_step = -meta_step
              go = .false.
            else
              call crecv(17948,colvar,ncolvar)
              ierr(3)=0
              call crecv(17949,dummy,1)
              waug=dummy(1)
              ierr(4)=0
              n = (meta_step-1)/commsize + 1
              colvar_his(:,n)=colvar(:)
              w_aug(n)=waug
            end if
          
          enddo
          
        end if
        call gisum(ierr(1),4,ierr(5))
        do i=1,4
          if (ierr(i)/=0) call Mfrz_Error(2531,0.d0)
        enddo
        meta_step = meta_step + 1
        
      else
        
c     purge any existing METADYNAMICS file
        
        if (onroot) then
          
          open(unit=mtd,file='METADYNAMICS',status='replace',
     x      iostat=ierr(1))
          close(unit=mtd)
          
        end if
        
      end if
      
      return
      
      end Subroutine Define_Metadynamics
      
      Function Fc(r,inner_cut,outer_cut)
      
c---------------------------------------------------------------------
c
c     Computes the smooth cut-off function used when computing an order
c     parameter as a function of pair separation.
c     
c     Author D. Quigley - University of Warwick
c     Copyright D. Quigley
c     
c---------------------------------------------------------------------
      
      implicit none
      
      real(8),parameter :: Pi=3.141592653589793238462643383279502884d0
      real(8),intent(in) :: r,inner_cut,outer_cut
      real(8) :: fc
      
      if ( r > outer_cut ) then
        fc = 0.0d0
      elseif ( r > inner_cut ) then
        fc = 0.5d0*cos((r-inner_cut)*Pi/(outer_cut-inner_cut))+0.5d0
      elseif ( r <= inner_cut ) then
        fc = 1.0d0
      else
        call Mfrz_Error(2532,r)
      end if
      
      return
      
      end Function Fc
      
      Function Dfc(r,inner_cut,outer_cut)
      
c---------------------------------------------------------------------
c     Computes the derivative of the smooth cut-off function used when
c     computing an order parameter as a function of pair separation.
c     
c     Author D. Quigley - University of Warwick
c     Copyright D. Quigley
c     
c---------------------------------------------------------------------
      
      implicit none
      
      real(8),parameter :: Pi=3.141592653589793238462643383279502884d0
      real(8),intent(in) :: r,inner_cut,outer_cut
      real(8) :: dfc
      
      if ( r > outer_cut ) then
        dfc = 0.0d0
      elseif ( r > inner_cut ) then
        dfc = -0.5d0*sin((r-inner_cut)*Pi/(outer_cut-inner_cut))
     x    *Pi/(outer_cut-inner_cut)
      else
        dfc = 0.0d0
      end if
      
      return
      
      end Function Dfc
  
      subroutine compute_steinhardt(imcon,natms)
      
c---------------------------------------------------------------------
c     
c     Computes nq4 Q4 and nq6 Q6 global order parameters.
c
c     Author D. Quigley - University of Warwick
c     Copyright D. Quigley October 2008
c     
c---------------------------------------------------------------------
      
      use config_module
      use site_module
      use setup_module,   only : mxlist
      use utility_module, only : images
      
      implicit none
      
c     Arguments
      
      integer,intent(in) :: imcon,natms
      
c     sqrt(pi/13) , sqrt(pi/9) , 1/3
      
      real(8),parameter :: rpibtt = 0.491590249d0
      real(8),parameter :: rpibn  = 0.590817950d0
      real(8),parameter :: third  = 1.0d0/3.0d0
      
c     Prefactors for spherical harmonics with l = 6
      
      real(8),parameter :: ypre6m6 =  0.48308411358006625446748d0
      real(8),parameter :: ypre6m5 =  1.67345245810009801701312d0
      real(8),parameter :: ypre6m4 =  0.35678126285399802686271d0
      real(8),parameter :: ypre6m3 =  0.65139048586771575166665d0
      real(8),parameter :: ypre6m2 =  0.32569524293385787583333d0
      real(8),parameter :: ypre6m1 =  0.41197551630114082055201d0
      real(8),parameter :: ypre6m0 =  0.06356920226762842462964d0
      real(8),parameter :: ypre6p1 = -0.41197551630114082055201d0
      real(8),parameter :: ypre6p2 =  0.32569524293385787583333d0
      real(8),parameter :: ypre6p3 = -0.65139048586771575166665d0
      real(8),parameter :: ypre6p4 =  0.35678126285399802686271d0
      real(8),parameter :: ypre6p5 = -1.67345245810009801701312d0
      real(8),parameter :: ypre6p6 =  0.48308411358006625446748d0
      
c     Prefactors for spherical harmonics with l = 4
      
      real(8),parameter :: ypre4m4 =  0.44253269244498261159038d0
      real(8),parameter :: ypre4m3 =  1.25167147089835228968013d0
      real(8),parameter :: ypre4m2 =  0.33452327177864460416856d0
      real(8),parameter :: ypre4m1 =  0.47308734787878004013351d0
      real(8),parameter :: ypre4m0 =  0.10578554691520430930396d0
      real(8),parameter :: ypre4p1 = -0.47308734787878004013351d0
      real(8),parameter :: ypre4p2 =  0.33452327177864460416856d0
      real(8),parameter :: ypre4p3 = -1.25167147089835228968013d0
      real(8),parameter :: ypre4p4 =  0.44253269244498261159038d0
      
c     Maximum no. of entries in solvation list
      
      integer :: maxneigh
      
c     Solvation shell information - Q4
      
      real(8),allocatable,dimension(:) :: solvx4,solvy4,solvz4
      real(8),allocatable,dimension(:) :: solvrmag4,solvimag4,solvrsq4
      integer,allocatable,dimension(:) :: solvlist4,solvtype4
      integer :: isolvmax4
      
c     Solvation shell information - Q6
      
      real(8),allocatable,dimension(:) :: solvx6,solvy6,solvz6
      real(8),allocatable,dimension(:) :: solvrmag6,solvimag6,solvrsq6
      integer,allocatable,dimension(:) :: solvlist6,solvtype6
      integer :: isolvmax6
      
c     separation vectors and powers thereof
      
      real(8),allocatable,dimension(:) :: xdf,ydf,zdf
      real(8) :: x,y,z
      real(8) :: x2,y2,z2,x3,y3,z3
      real(8) :: x4,y4,z4,x5,y5,z5
      real(8) :: x6,y6,z6
      real(8) :: invrc,invrs

c     list of separation vectors
      integer :: numdst
      integer,allocatable,dimension(:) :: dstlst

      
c     Comms buffers
      
      real(8),allocatable,dimension(:) :: buff1,buff2
      
c     Temporaries
      
      real(8) :: tmpsq,f_ij,df_ij,ReYlm,ImYlm,tmpvar
      
c     Loop counters
      
      integer :: iatm1,iatm0,iatm,isite,istrd,ii,isolv4,isolv6,isolv
      integer :: idi,idj,limit,nn,k,jatm,jsite,q4type,q6type
      integer :: itype,jtype,l,m,iq
      
      maxneigh = 100            ! Max number of atoms in coordination shell

      ierr = 0                  ! Error flags
      

c     DQ - modified 10/12/11, arrays now big enough
c     to hold maximum number of neighbours plus
c     maximum number of excluded atoms.
      allocate(xdf(1:mxlist+mxexcl),stat=ierr(1))
      allocate(ydf(1:mxlist+mxexcl),stat=ierr(2))
      allocate(zdf(1:mxlist+mxexcl),stat=ierr(3))

c     DQ - modified 10/12/11, array to hold a list of
c     all atom entries in the above three arrays
      allocate(dstlst(1:mxlist+mxexcl),stat=ierr(4))
      
      allocate(solvx4(1:maxneigh),stat=ierr(5))
      allocate(solvy4(1:maxneigh),stat=ierr(6))
      allocate(solvz4(1:maxneigh),stat=ierr(7))
      allocate(solvrmag4(1:maxneigh),stat=ierr(8))
      allocate(solvimag4(1:maxneigh),stat=ierr(9))
      allocate(solvrsq4 (1:maxneigh),stat=ierr(10))
      allocate(solvlist4(1:maxneigh),stat=ierr(11)) 
      allocate(solvtype4(1:maxneigh),stat=ierr(12))
      
      allocate(solvx6(1:maxneigh),stat=ierr(13))
      allocate(solvy6(1:maxneigh),stat=ierr(14))
      allocate(solvz6(1:maxneigh),stat=ierr(15))
      allocate(solvrmag6(1:maxneigh),stat=ierr(16))
      allocate(solvimag6(1:maxneigh),stat=ierr(17))
      allocate(solvrsq6 (1:maxneigh),stat=ierr(18))
      allocate(solvlist6(1:maxneigh),stat=ierr(19)) 
      allocate(solvtype6(1:maxneigh),stat=ierr(20))
      if (any(ierr/=0)) call Mfrz_Error(2533,0.d0) 
      
      allocate(buff1(1:18*nq4+26*nq6),stat=ierr(1))
      allocate(buff2(1:18*nq4+26*nq6),stat=ierr(2))    
      if (any(ierr/=0)) call Mfrz_Error(2534,0.d0)
      
c     Zero accumulators used in Steinhardt order parameters
      
      ReQ6bar = 0.0d0
      ImQ6bar = 0.0d0
      ReQ4bar = 0.0d0
      ImQ4bar = 0.0d0
      
c     Set atoms looped over by current rank
      
      iatm0 = myrank+1
      iatm1 = natms
      istrd = commsize
      
      ii = 0
      do iatm = iatm0,iatm1,istrd
        
c --------------------------------------------------------------
c     Build a list of the required connections to iatm. This  
c     differs depending on the version of DLPOLY we are using.
c     First we loop over atoms in the neighbour list of iatm.
c---------------------------------------------------------------
      
      ii = ii + 1
      isite=ltype(iatm)
      limit=lentry(ii)
      nn = 0
      do k = 1,limit
        
        jatm  = list(ii,k)
        jsite = ltype(jatm)
        
        if ( q4site(jsite,isite)+q6site(jsite,isite)==0 ) cycle
        
        nn = nn + 1

        dstlst(nn) = jatm
        
        xdf(nn)=xxx(jatm)-xxx(iatm)
        ydf(nn)=yyy(jatm)-yyy(iatm)
        zdf(nn)=zzz(jatm)-zzz(iatm) 
        
      end do

c --------------------------------------------------------------
c     Next we loop over the excluded atom list of iatm and add 
c     and pairs needed for computation of the current OP.
c---------------------------------------------------------------

ccc   DEBUG
ccc      write(0,'("atom ",I5," has ",I5," excluded interactions")')
ccc     x iatm,mtd_nexatm(iatm)

      do k = 1,mtd_nexatm(ii)

         jatm  = mtd_lexatm(ii,k)
         jsite = ltype(jatm)

ccc   DEBUG
ccc         write(0,'("Interaction with atom ",I5," is excluded. ")')jatm

         if ( q4site(jsite,isite)+q6site(jsite,isite)==0 ) cycle
        
         nn = nn + 1
        
         dstlst(nn) = jatm

         xdf(nn)=xxx(jatm)-xxx(iatm)
         ydf(nn)=yyy(jatm)-yyy(iatm)
         zdf(nn)=zzz(jatm)-zzz(iatm) 

      end do

ccc   DEBUG
ccc      write(0,'("Num neighbours to consider for atom ",I5," : ",I5)')
ccc     x iatm,nn

      numdst = nn

      call images(imcon,0,1,nn,cell,xdf,ydf,zdf)
      nn = 0
      isolvmax4 = 0
      isolvmax6 = 0
      isolv4 = 0
      isolv6 = 0

      do k = 1,numdst
        jatm  = dstlst(k)
        jsite = ltype(jatm)
        
        if ( q4site(jsite,isite)+q6site(jsite,isite)==0 ) cycle
        
        nn = nn + 1
        
        q4type = q4site(jsite,isite)
        q6type = q6site(jsite,isite)
        
        tmpsq = xdf(nn)*xdf(nn)+ydf(nn)*ydf(nn)+zdf(nn)*zdf(nn)
        
        if (nq4>0) then
          
c     Add to solvation lists if within cut-off 
          
          if (  tmpsq < q4cutoff(2,q4type)**2 ) then
            isolv4 = isolv4 + 1
            solvlist4(isolv4) = jatm
            solvrsq4(isolv4)  = tmpsq
            solvrmag4(isolv4) = sqrt(tmpsq)
            solvimag4(isolv4) = 1.0d0/solvrmag4(isolv4)
            solvx4(isolv4)    = xdf(nn)
            solvy4(isolv4)    = ydf(nn)
            solvz4(isolv4)    = zdf(nn)
            solvtype4(isolv4) = q4type
          end if
        end if
        
        if (nq6>0) then
          
c     Add to solvation lists if within cut-off 
          
          if (  tmpsq < q6cutoff(2,q6type)**2 ) then
            isolv6 = isolv6 + 1
            solvlist6(isolv6) = jatm
            solvrsq6(isolv6)  = tmpsq
            solvrmag6(isolv6) = sqrt(tmpsq)
            solvimag6(isolv6) = 1.0d0/solvrmag6(isolv6)
            solvx6(isolv6)    = xdf(nn)
            solvy6(isolv6)    = ydf(nn)
            solvz6(isolv6)    = zdf(nn)
            solvtype6(isolv6) = q6type
          end if
        end if
        isolvmax4 = isolv4
        isolvmax6 = isolv6
        if ((isolv4>maxneigh) .or. (isolv6>maxneigh))
     x    call Mfrz_Error(2535,0.d0)
        
      end do                    ! end loop over k

ccc      write(0,'("Num in range for OPs on atom ",I5," : ",I5)')
ccc     x iatm,isolvmax4

      
c---------------------------------------------------------
c     Compute Q4 Steinhardt order parameters              
c---------------------------------------------------------
      
      if ( (nq4>0).and.isolvmax4>0 ) then
        
        do isolv4 = 1,isolvmax4
          
          jatm = solvlist4(isolv4)
          itype = solvtype4(isolv4)
          
          invrc = solvimag4(isolv4)**6
          invrs = solvimag4(isolv4)**4
          
          x     = solvx4(isolv4)
          y     = solvy4(isolv4)
          z     = solvz4(isolv4)
          
          f_ij  = fc(solvrmag4(isolv4),q4cutoff(1,itype),
     x              q4cutoff(2,itype))
          
          x2 = x*x
          y2 = y*y
          z2 = z*z
          
          x3 = x2*x
          y3 = y2*y
          z3 = z2*z
          
          x4 = x2*x2
          y4 = y2*y2
          z4 = z2*z2
          
          x5 = x4*x
          y5 = y4*y
          z5 = z4*z
          
          x6 = x4*x2
          y6 = y4*y2
          z6 = z4*z2
          
c----------------------------------------------------------
c     Real and imaginary contribution to Q4bar(-4/+4)
c----------------------------------------------------------
          
          ReYlm = ypre4m4*invrs*(x4-6.d0*x2*y2+y4)
          ImYlm = ypre4m4*invrs*(-4.d0*x3*y+4.d0*x*y3)
          
          ReQ4bar(-4,itype) = ReQ4bar(-4,itype) + f_ij*ReYlm
          ImQ4bar(-4,itype) = ImQ4bar(-4,itype) + f_ij*ImYlm

          ReQ4bar(+4,itype) = ReQ4bar(+4,itype) + f_ij*ReYlm
          ImQ4bar(+4,itype) = ImQ4bar(+4,itype) - f_ij*ImYlm
          
c----------------------------------------------------------
c     Real and imaginary contribution to Q4bar(-3/+3)
c----------------------------------------------------------
          
          ReYlm = ypre4m3*invrs*z*(x3-3.d0*x*y2)
          ImYlm = ypre4m3*invrs*z*(-3.d0*x2*y+y3)
          
          ReQ4bar(-3,itype) = ReQ4bar(-3,itype) + f_ij*ReYlm
          ImQ4bar(-3,itype) = ImQ4bar(-3,itype) + f_ij*ImYlm

          ReQ4bar(+3,itype) = ReQ4bar(+3,itype) - f_ij*ReYlm
          ImQ4bar(+3,itype) = ImQ4bar(+3,itype) + f_ij*ImYlm
          
c----------------------------------------------------------
c     Real and imaginary contribution to Q4bar(-2/+2)
c----------------------------------------------------------
          
          ReYlm = -ypre4m2*invrs*(x2-y2)*(-6.d0*z2+x2+y2)
          ImYlm = ypre4m2*invrs*2.d0*(-6.d0*z2+x2+y2)*x*y
          
          ReQ4bar(-2,itype) = ReQ4bar(-2,itype) + f_ij*ReYlm
          ImQ4bar(-2,itype) = ImQ4bar(-2,itype) + f_ij*ImYlm

          ReQ4bar(+2,itype) = ReQ4bar(+2,itype) + f_ij*ReYlm
          ImQ4bar(+2,itype) = ImQ4bar(+2,itype) - f_ij*ImYlm
          
c----------------------------------------------------------
c     Real and imaginary contribution to Q4bar(-1/+1)
c----------------------------------------------------------
          
          ReYlm = -ypre4m1*invrs*z*(-4.d0*z2+3.d0*x2+3.d0*y2)*x
          ImYlm = ypre4m1*invrs*z*(-4.d0*z2+3.d0*x2+3.d0*y2)*y
          
          ReQ4bar(-1,itype) = ReQ4bar(-1,itype) + f_ij*ReYlm
          ImQ4bar(-1,itype) = ImQ4bar(-1,itype) + f_ij*ImYlm

          ReQ4bar(+1,itype) = ReQ4bar(+1,itype) - f_ij*ReYlm
          ImQ4bar(+1,itype) = ImQ4bar(+1,itype) + f_ij*ImYlm
          
c----------------------------------------------------------
c     Real and imaginary contribution to Q4bar(0)
c----------------------------------------------------------
          
          ReYlm = ypre4m0*invrs*(8.d0*z4-24.d0*z2*x2-24.d0*z2*y2+
     x      3.d0*x4+6.d0*x2*y2+3.d0*y4)
          
          ReQ4bar(0,itype)  = ReQ4bar(0,itype) + f_ij*ReYlm
          
        end do                  ! end loop over connection list for iatm
        
      end if                    ! end if computing Q4
      
c------------------------------------------------
c     Compute Q6 Steinhardt order parameters     
c------------------------------------------------
      
      if ( (nq6>0).and.isolvmax6>0 ) then
        
        do isolv6 = 1,isolvmax6
          
          jatm = solvlist6(isolv6)
          itype = solvtype6(isolv6)
          
          invrc = solvimag6(isolv6)**6
          invrs = solvimag6(isolv6)**4
          
          x     = solvx6(isolv6)
          y     = solvy6(isolv6)
          z     = solvz6(isolv6)
          
          f_ij  =  fc(solvrmag6(isolv6),q6cutoff(1,itype),
     x               q6cutoff(2,itype))
          
          x2 = x*x
          y2 = y*y
          z2 = z*z
          
          x3 = x2*x
          y3 = y2*y
          z3 = z2*z
          
          x4 = x2*x2
          y4 = y2*y2
          z4 = z2*z2
          
          x5 = x4*x
          y5 = y4*y
          z5 = z4*z
          
          x6 = x4*x2
          y6 = y4*y2
          z6 = z4*z2
          
c----------------------------------------------------------
c     Real and imaginary conribution to Q6bar(-6/+6)
c----------------------------------------------------------
          
          ReYlm = ypre6m6*invrc*(x6-15.0d0*x4*y2+15.0d0*x2*y4-y6)
          ImYlm = ypre6m6*invrc*(-6.0d0*x5*y+20.0d0*x3*y3-6.0d0*x*y5)
          
          ReQ6bar(-6,itype) = ReQ6bar(-6,itype) + f_ij*ReYlm
          ImQ6bar(-6,itype) = ImQ6bar(-6,itype) + f_ij*ImYlm

          ReQ6bar(+6,itype) = ReQ6bar(+6,itype) + f_ij*ReYlm
          ImQ6bar(+6,itype) = ImQ6bar(+6,itype) - f_ij*ImYlm
          
c----------------------------------------------------------
c     Real and imaginary conribution to Q6bar(-5/+5)
c----------------------------------------------------------
          
          ReYlm = -ypre6m5*invrc*z*(-x5+10.0d0*x3*y2-5.0d0*x*y4)
          ImYlm = -ypre6m5*invrc*z*(5.0d0*x4*y-10.0d0*x2*y3+y5)
          
          ReQ6bar(-5,itype) = ReQ6bar(-5,itype) + f_ij*ReYlm
          ImQ6bar(-5,itype) = ImQ6bar(-5,itype) + f_ij*ImYlm

          ReQ6bar(+5,itype) = ReQ6bar(+5,itype) - f_ij*ReYlm
          ImQ6bar(+5,itype) = ImQ6bar(+5,itype) + f_ij*ImYlm
          
c----------------------------------------------------------
c     Real and imaginary conribution to Q6bar(-4/+4)
c----------------------------------------------------------
          
          ReYlm = ypre6m4*invrc*(10.0d0*z2-x2-y2)*(x4-6.0d0*x2*y2+y4)
          ImYlm = ypre6m4*invrc*(10.0d0*z2-x2-y2)*(-4.0d0*x3*y+
     x            4.0d0*x*y3)
          
          ReQ6bar(-4,itype) = ReQ6bar(-4,itype) + f_ij*ReYlm
          ImQ6bar(-4,itype) = ImQ6bar(-4,itype) + f_ij*ImYlm

          ReQ6bar(+4,itype) = ReQ6bar(+4,itype) + f_ij*ReYlm
          ImQ6bar(+4,itype) = ImQ6bar(+4,itype) - f_ij*ImYlm
          
c----------------------------------------------------------
c     Real and imaginary conribution to Q6bar(-3/+3)
c----------------------------------------------------------
          
          ReYlm = -ypre6m3*invrc*z*(8.0d0*z2-3.0d0*x2-3.0d0*y2)*
     x             (-x3+3.0d0*x*y2)
          ImYlm = -ypre6m3*invrc*z*(8.0d0*z2-3.0d0*x2-3.0d0*y2)*
     x             (3.0d0*x2*y-y3)
          
          ReQ6bar(-3,itype) = ReQ6bar(-3,itype) + f_ij*ReYlm
          ImQ6bar(-3,itype) = ImQ6bar(-3,itype) + f_ij*ImYlm

          ReQ6bar(+3,itype) = ReQ6bar(+3,itype) - f_ij*ReYlm
          ImQ6bar(+3,itype) = ImQ6bar(+3,itype) + f_ij*ImYlm
          
c----------------------------------------------------------
c     Real and imaginary conribution to Q6bar(-2/+2)
c----------------------------------------------------------

          ReYlm =  ypre6m2*invrc*(16.0d0*z4-16.0d0*z2*x2-16.0d0*z2*y2+
     x             x4+2.0d0*x2*y2+y4)*(x2-y2)
          ImYlm = -ypre6m2*invrc*2.0d0*(16.0d0*z4-16.0d0*z2*x2-16.0d0*
     x             z2*y2+x4+2.0d0*x2*y2+y4)*x*y
          
          ReQ6bar(-2,itype) = ReQ6bar(-2,itype) + f_ij*ReYlm
          ImQ6bar(-2,itype) = ImQ6bar(-2,itype) + f_ij*ImYlm

          ReQ6bar(+2,itype) = ReQ6bar(+2,itype) + f_ij*ReYlm
          ImQ6bar(+2,itype) = ImQ6bar(+2,itype) - f_ij*ImYlm
          
c----------------------------------------------------------
c     Real and imaginary conribution to Q6bar(-1/+1)
c----------------------------------------------------------
          
          ReYlm =  ypre6m1*z*invrc*(8.0d0*z4-20.0d0*z2*x2-20.0d0*z2*y2+
     x             5.0d0*x4+10.0d0*x2*y2+5.0d0*y4)*x
          ImYlm = -ypre6m1*z*invrc*(8.0d0*z4-20.0d0*z2*x2-20.0d0*z2*y2+
     x             5.0d0*x4+10.0d0*x2*y2+5.0d0*y4)*y
          
          ReQ6bar(-1,itype) = ReQ6bar(-1,itype) + f_ij*ReYlm
          ImQ6bar(-1,itype) = ImQ6bar(-1,itype) + f_ij*ImYlm

          ReQ6bar(+1,itype) = ReQ6bar(+1,itype) - f_ij*ReYlm
          ImQ6bar(+1,itype) = ImQ6bar(+1,itype) + f_ij*ImYlm
          
c----------------------------------------------------------
c     Real and imaginary conribution to Q6bar(0)
c----------------------------------------------------------
          
          ReYlm =  ypre6m0*invrc*(16.0d0*z6-120.0d0*z4*x2-120.0d0*z4*
     x             y2+90.0d0*z2*x4+180.0d0*z2*x2*y2+90.0d0*z2*y4-5.0d0
     x             *x6-15.0d0*x4*y2-15.0d0*x2*y4-5.0d0*y6)
          
          ReQ6bar(0,itype) = ReQ6bar(0,itype) + f_ij*ReYlm
          
        end do                  ! end loop over connection list for iatm
        
      end if                    ! end if computing Q6
      
      end do                    ! end loop over iatm
      
c-----------------------------------------------
c     Global summation of order parameters      
c-----------------------------------------------
      
      l = 1
      do itype = 1,nq4
        do m = -4,4
          buff1(l) = ReQ4bar(m,itype)
          l = l + 1
        end do
        do m = -4,4
          buff1(l) = ImQ4bar(m,itype)
          l = l + 1
        end do
      end do
      do itype = 1,nq6
        do m = -6,6
          buff1(l) = ReQ6bar(m,itype)
          l = l + 1
        end do
        do m = -6,6
          buff1(l) = ImQ6bar(m,itype)
          l = l + 1
        end do
      end do
      
      if (commsize>1)   call gdsum(buff1,18*nq4+26*nq6,buff2)
      
      l = 1
      do itype = 1,nq4
        do m = -4,4
          ReQ4bar(m,itype) = buff1(l)
          l = l + 1
        end do
        do m = -4,4
          ImQ4bar(m,itype) = buff1(l)
          l = l + 1
        end do
      end do
      do itype = 1,nq6
        do m = -6,6
          ReQ6bar(m,itype) = buff1(l)
          l = l + 1
        end do
        do m = -6,6
          ImQ6bar(m,itype) = buff1(l)
          l = l + 1
        end do
      end do
      
c---------------------------------------------------
c     Final computation of global order parameters  
c---------------------------------------------------
      
      l = 1
      do iq = 1,nq4
        tmpvar = 0.0d0
        do m = 1,18
          tmpvar = tmpvar + buff1(l)**2
          l = l + 1
        end do
        q4_global(iq) = 4.0d0*rpibn*sqrt(tmpvar)/
     x                  dble(q4ninc(iq)*q4nn(iq))
      end do
      
      do iq = 1,nq6
        tmpvar = 0.0d0
        do m = 1,26
          tmpvar = tmpvar + buff1(l)**2
          l = l + 1
        end do
        q6_global(iq) = 4.0d0*rpibtt*sqrt(tmpvar)/
     x                  dble(q6ninc(iq)*q6nn(iq))
      end do
      
c     Tidy up
      
      deallocate(xdf,stat=ierr(1))
      deallocate(ydf,stat=ierr(2))
      deallocate(zdf,stat=ierr(3))

      deallocate(dstlst,stat=ierr(4))
      
      deallocate(solvx4,stat=ierr(5))
      deallocate(solvy4,stat=ierr(6))
      deallocate(solvz4,stat=ierr(7))
      deallocate(solvrmag4,stat=ierr(8))
      deallocate(solvimag4,stat=ierr(9))
      deallocate(solvrsq4 ,stat=ierr(10))
      deallocate(solvlist4,stat=ierr(11)) 
      deallocate(solvtype4,stat=ierr(12))
      
      deallocate(solvx6,stat=ierr(13))
      deallocate(solvy6,stat=ierr(14))
      deallocate(solvz6,stat=ierr(15))
      deallocate(solvrmag6,stat=ierr(16))
      deallocate(solvimag6,stat=ierr(17))
      deallocate(solvrsq6 ,stat=ierr(18))
      deallocate(solvlist6,stat=ierr(19)) 
      deallocate(solvtype6,stat=ierr(20))
      if (any(ierr/=0)) call Mfrz_Error(2536,0.d0) 
      
      deallocate(buff1,stat=ierr(1))
      deallocate(buff2,stat=ierr(2))    
      if (any(ierr/=0)) call Mfrz_Error(2537,0.d0)

      return
      
      end Subroutine Compute_Steinhardt
      
      Subroutine Compute_Steinhardt_Forces(imcon,natms,engord,virord)

c---------------------------------------------------------------------
c     
c     Computes forces from nq4 Q4 and nq6 Q6 global order parameters.
c
c     Author D. Quigley - University of Warwick
c     Copyright D. Quigley October 2008
c     
c---------------------------------------------------------------------
      
      use config_module
      use site_module
      use setup_module,   only : mxlist
      use utility_module, only : images
      implicit none
      
c     Arguments
      
      integer,intent(in) :: imcon,natms
      
      real(8),intent(inout) :: engord,virord
      
c     sqrt(pi/13) , sqrt(pi/9) , 1/3
      
      real(8),parameter :: rpibtt = 0.491590249d0
      real(8),parameter :: rpibn  = 0.590817950d0
      real(8),parameter :: third  = 1.0d0/3.0d0
      
c     Prefactors for spherical harmonics with l = 6
      
      real(8),parameter :: ypre6m6 =  0.48308411358006625446748d0
      real(8),parameter :: ypre6m5 =  1.67345245810009801701312d0
      real(8),parameter :: ypre6m4 =  0.35678126285399802686271d0
      real(8),parameter :: ypre6m3 =  0.65139048586771575166665d0
      real(8),parameter :: ypre6m2 =  0.32569524293385787583333d0
      real(8),parameter :: ypre6m1 =  0.41197551630114082055201d0
      real(8),parameter :: ypre6m0 =  0.06356920226762842462964d0
      real(8),parameter :: ypre6p1 = -0.41197551630114082055201d0
      real(8),parameter :: ypre6p2 =  0.32569524293385787583333d0
      real(8),parameter :: ypre6p3 = -0.65139048586771575166665d0
      real(8),parameter :: ypre6p4 =  0.35678126285399802686271d0
      real(8),parameter :: ypre6p5 = -1.67345245810009801701312d0
      real(8),parameter :: ypre6p6 =  0.48308411358006625446748d0
      
c     Prefactors for speherical harmonics with l = 4
      
      real(8),parameter :: ypre4m4 =  0.44253269244498261159038d0
      real(8),parameter :: ypre4m3 =  1.25167147089835228968013d0
      real(8),parameter :: ypre4m2 =  0.33452327177864460416856d0
      real(8),parameter :: ypre4m1 =  0.47308734787878004013351d0
      real(8),parameter :: ypre4m0 =  0.10578554691520430930396d0
      real(8),parameter :: ypre4p1 = -0.47308734787878004013351d0
      real(8),parameter :: ypre4p2 =  0.33452327177864460416856d0
      real(8),parameter :: ypre4p3 = -1.25167147089835228968013d0
      real(8),parameter :: ypre4p4 =  0.44253269244498261159038d0
      
c     Maximum no. of entries in solvation list
      
      integer :: maxneigh
      
c     Solvation shell information - Q4
      
      real(8),allocatable,dimension(:) :: solvx4,solvy4,solvz4
      real(8),allocatable,dimension(:) :: solvrmag4,solvimag4,solvrsq4
      integer,allocatable,dimension(:) :: solvlist4,solvtype4
      integer :: isolvmax4
      
c     Solvation shell information - Q6
      
      real(8),allocatable,dimension(:) :: solvx6,solvy6,solvz6
      real(8),allocatable,dimension(:) :: solvrmag6,solvimag6,solvrsq6
      integer,allocatable,dimension(:) :: solvlist6,solvtype6
      integer :: isolvmax6
      
c     Prefactors arising from derivative of bias potential

      real(8),allocatable,dimension(:) :: q4prefactor,q6prefactor
      
c     Separation vectors and powers thereof
      
      real(8),allocatable,dimension(:) :: xdf,ydf,zdf
      real(8) :: x,y,z
      real(8) :: x2,y2,z2,x3,y3,z3
      real(8) :: x4,y4,z4,x5,y5,z5
      real(8) :: x6,y6,z6
      real(8) :: invrc,invrs,invrq

c     list of separation vectors
      integer :: numdst
      integer,allocatable,dimension(:) :: dstlst
      
c     Comms buffers
      
      real(8),allocatable,dimension(:) :: buff1,buff2
      
c     Temporaries

      real(8) :: tmpsq,f_ij,df_ij,ReYlm,ImYlm,tmpvar,invrN
      real(8) :: fx,fy,fz,fx2,fy2,fz2,prefactor2,fx1,fy1,fz1
      real(8) :: strs1,strs2,strs3,strs4,strs5,strs6,strs7,strs8,strs9
      
      integer :: iatm1,iatm0,iatm,isite,istrd,ii,isolv4,isolv6,isolv
      integer :: idi,idj,limit,nn,k,jatm,jsite,q4type,q6type
      integer :: itype,jtype,l,m,iq
      
      maxneigh = 500         ! Max number of atoms in coordination shell
      
      ierr = 0               ! Error flags
      
c     DQ - modified 10/12/11, arrays now big enough
c     to hold maximum number of neighbours plus
c     maximum number of excluded atoms.
      allocate(xdf(1:mxlist+mxexcl),stat=ierr(1))
      allocate(ydf(1:mxlist+mxexcl),stat=ierr(2))
      allocate(zdf(1:mxlist+mxexcl),stat=ierr(3))

c     DQ - modified 10/12/11, array to hold a list of
c     all atom entries in the above three arrays
      allocate(dstlst(1:mxlist+mxexcl),stat=ierr(4))
      
      allocate(solvx4(1:maxneigh),stat=ierr(5))
      allocate(solvy4(1:maxneigh),stat=ierr(6))
      allocate(solvz4(1:maxneigh),stat=ierr(7))
      allocate(solvrmag4(1:maxneigh),stat=ierr(8))
      allocate(solvimag4(1:maxneigh),stat=ierr(9))
      allocate(solvrsq4 (1:maxneigh),stat=ierr(10))
      allocate(solvlist4(1:maxneigh),stat=ierr(11)) 
      allocate(solvtype4(1:maxneigh),stat=ierr(12))
      
      allocate(solvx6(1:maxneigh),stat=ierr(13))
      allocate(solvy6(1:maxneigh),stat=ierr(14))
      allocate(solvz6(1:maxneigh),stat=ierr(15))
      allocate(solvrmag6(1:maxneigh),stat=ierr(16))
      allocate(solvimag6(1:maxneigh),stat=ierr(17))
      allocate(solvrsq6 (1:maxneigh),stat=ierr(18))
      allocate(solvlist6(1:maxneigh),stat=ierr(19)) 
      allocate(solvtype6(1:maxneigh),stat=ierr(20))
      if (any(ierr/=0)) call Mfrz_Error(2538,0.d0)
      
      allocate(buff1(1:18*nq4+26*nq6),stat=ierr(1))
      allocate(buff2(1:18*nq4+26*nq6),stat=ierr(2))    
      if (any(ierr/=0)) call Mfrz_Error(2534,0.d0)
      
      allocate(q4prefactor(1:nq4),stat=ierr(1))
      allocate(q6prefactor(1:nq6),stat=ierr(2))
      if (any(ierr/=0)) call Mfrz_Error(2540,0.d0)
      
c     Compute the prefactors associated from dV_aug/d_q4
      
      k = 1
      do iq = 1,nq4
        invrN = 1.0d0/dble(q4ninc(iq)*q4nn(iq))
        q4prefactor(iq) = -16.0d0*(rpibn**2)*(invrN**2)*dcolvar(k)/
     x                    Q4_global(iq)
        k = k + 1
      end do
      
c     Compute the prefactors associated from dV_aug/d_q6
      
      do iq = 1,nq6
        invrN = 1.0d0/dble(q6ninc(iq)*q6nn(iq))
        q6prefactor(iq) = -16.0d0*(rpibtt**2)*(invrN**2)*dcolvar(k)/
     x                    Q6_global(iq)
        k = k + 1
      end do

c     write(0,'("DEBUG : q4prefactors = ",5F15.6)')q4prefactor
      
c     Set atoms looper over by current rank
      
      iatm0 = myrank+1
      iatm1 = natms
      istrd = commsize
      
      strs1 = 0.0d0
      strs2 = 0.0d0
      strs3 = 0.0d0
      strs4 = 0.0d0
      strs5 = 0.0d0
      strs6 = 0.0d0
      strs7 = 0.0d0
      strs8 = 0.0d0
      strs9 = 0.0d0
      
      ii = 0
      do iatm = iatm0,iatm1,istrd
        
c --------------------------------------------------------------
c     Build a list of the required connections to iatm. This  
c     differs depending on the version of DLPOLY we are using.
c     First we loop over atoms in the neighbour list of iatm.
c---------------------------------------------------------------
        
        ii = ii + 1
        isite=ltype(iatm)
        limit=lentry(ii)
        
        nn = 0
        do k = 1,limit
          
          jatm  = list(ii,k)
          jsite = ltype(jatm)
          
          if ( q4site(jsite,isite)+q6site(jsite,isite)==0 ) cycle
          
          nn = nn + 1

          dstlst(nn) = jatm
          
          xdf(nn)=xxx(jatm)-xxx(iatm)
          ydf(nn)=yyy(jatm)-yyy(iatm)
          zdf(nn)=zzz(jatm)-zzz(iatm) 
          
        end do

c --------------------------------------------------------------
c     Next we loop over the excluded atom list of iatm and add 
c     and pairs needed for computation of the current OP.
c---------------------------------------------------------------

ccc   DEBUG
ccc        write(0,'("atom ",I5," has ",I5," excluded interactions")')
ccc     x       iatm,mtd_nexatm(iatm)

        do k = 1,mtd_nexatm(ii)
           
           jatm  = mtd_lexatm(ii,k)
           jsite = ltype(jatm)
           
ccc   DEBUG
ccc           write(0,'("Interaction with atom ",I5," is excluded. ")')jatm
           
           if ( q4site(jsite,isite)+q6site(jsite,isite)==0 ) cycle
        
           nn = nn + 1
           
           dstlst(nn) = jatm
           
           xdf(nn)=xxx(jatm)-xxx(iatm)
           ydf(nn)=yyy(jatm)-yyy(iatm)
           zdf(nn)=zzz(jatm)-zzz(iatm) 
           
        end do
        
ccc   DEBUG
ccc        write(0,*)
        
        numdst = nn
        
        call images(imcon,0,1,nn,cell,xdf,ydf,zdf)
        
        nn = 0
        isolvmax4 = 0
        isolvmax6 = 0
        isolv4 = 0
        isolv6 = 0
        do k = 1,numdst
          jatm  = dstlst(k)
          jsite = ltype(jatm)
          
          if ( q4site(jsite,isite)+q6site(jsite,isite)==0 ) cycle
          
          nn = nn + 1
          
          q4type = q4site(jsite,isite)
          q6type = q6site(jsite,isite)
          
          tmpsq = xdf(nn)*xdf(nn)+ydf(nn)*ydf(nn)+zdf(nn)*zdf(nn)
          if (nq4>0) then
            
c     Add to solvation lists if within cut-off 
            
            if (  tmpsq < q4cutoff(2,q4type)**2 ) then
              isolv4 = isolv4 + 1
              solvlist4(isolv4) = jatm
              solvrsq4(isolv4)  = tmpsq
              solvrmag4(isolv4) = sqrt(tmpsq)
              solvimag4(isolv4) = 1.0d0/solvrmag4(isolv4)
              solvx4(isolv4)    = xdf(nn)
              solvy4(isolv4)    = ydf(nn)
              solvz4(isolv4)    = zdf(nn)
              solvtype4(isolv4) = q4type
            end if
          end if
          
          if (nq6>0) then
            
c     Add to solvation lists if within cut-off 
            
            if (  tmpsq < q6cutoff(2,q6type)**2 ) then
              isolv6 = isolv6 + 1
              solvlist6(isolv6) = jatm
              solvrsq6(isolv6)  = tmpsq
              solvrmag6(isolv6) = sqrt(tmpsq)
              solvimag6(isolv6) = 1.0d0/solvrmag6(isolv6)
              solvx6(isolv6)    = xdf(nn)
              solvy6(isolv6)    = ydf(nn)
              solvz6(isolv6)    = zdf(nn)
              solvtype6(isolv6) = q6type
            end if
          end if
          isolvmax4 = isolv4
          isolvmax6 = isolv6
          if ((isolv4>maxneigh) .or. (isolv6>maxneigh))
     x      call Mfrz_Error(2535,0.d0)
          
        end do                  ! end loop over k
        
c---------------------------------------------------------
c---------------------------------------------------------
c     Compute forces arising from  Q4 order parameters    
c---------------------------------------------------------
c---------------------------------------------------------
        
        if ( (nq4>0).and.isolvmax4>0 ) then
          
          do isolv4 = 1,isolvmax4
            
            jatm = solvlist4(isolv4)
            itype = solvtype4(isolv4)
            
            invrc = solvimag4(isolv4)**6
            invrq = solvimag4(isolv4)**8
            invrs = solvimag4(isolv4)**4
            
            x     = solvx4(isolv4)
            y     = solvy4(isolv4)
            z     = solvz4(isolv4)
            
            f_ij  =  fc(solvrmag4(isolv4),q4cutoff(1,itype),
     x                 q4cutoff(2,itype))
            df_ij = dfc(solvrmag4(isolv4),q4cutoff(1,itype),
     x                 q4cutoff(2,itype))
            
            x2 = x*x
            y2 = y*y
            z2 = z*z
            
            x3 = x2*x
            y3 = y2*y
            z3 = z2*z
            
            x4 = x2*x2
            y4 = y2*y2
            z4 = z2*z2
            
            x5 = x4*x
            y5 = y4*y
            z5 = z4*z
            
            x6 = x4*x2
            y6 = y4*y2
            z6 = z4*z2
            
            fx = 0.0d0 
            fy = 0.0d0 
            fz = 0.0d0
            
c-------------------------------------
c     Gradient of f_ij w.r.t. r_{j}   
c-------------------------------------
            
            fx2 = df_ij*x*solvimag4(isolv4)
            fy2 = df_ij*y*solvimag4(isolv4)
            fz2 = df_ij*z*solvimag4(isolv4)
            
c--------------------------------------------------------
c     Real and imaginary spherical harmonics for m = -4  
c--------------------------------------------------------

            ReYlm = ypre4m4*invrs*(x4-6.d0*x2*y2+y4)
            ImYlm = ypre4m4*invrs*(-4.d0*x3*y+4.d0*x*y3)
            
c--------------------------------------------------
c     Force contributions from m = -4 (real part)
c--------------------------------------------------
            
            prefactor2 = 2.0d0*q4prefactor(itype)*ReQ4bar(-4,itype)
            
c--------------------------------------------
c     Gradient of Re(Y_{4,-4}) w.r.t r_{j}
c--------------------------------------------
            
            fx1 =  invrc*ypre4m4*4.d0*x*(4.d0*x2*y2-4.d0*y4+z2*x2-
     x             3.d0*z2*y2)
            fy1 = -invrc*ypre4m4*4.d0*y*(4.d0*x4-4.d0*x2*y2+3.d0*
     x             z2*x2-z2*y2)
            fz1 = -invrc*ypre4m4*4.d0*z*(x4-6.d0*x2*y2+y4)
            
            fx = fx + prefactor2*(f_ij*fx1 + fx2*ReYlm)
            fy = fy + prefactor2*(f_ij*fy1 + fy2*ReYlm)
            fz = fz + prefactor2*(f_ij*fz1 + fz2*ReYlm)
            
c--------------------------------------------------------
c     Force contributions from m = -4 (imaginary part) 
c--------------------------------------------------------
            
            prefactor2 = 2.0d0*q4prefactor(itype)*ImQ4bar(-4,itype)
            
c----------------------------------------------
c     Gradient of Im(Y_{4,-4}) w.r.t r_{j}
c----------------------------------------------

            fx1 =  invrc*ypre4m4*4.d0*y*(x4-6.d0*x2*y2-3.d0*z2*x2+
     x             y4+z2*y2)
            fy1 = -invrc*ypre4m4*4.d0*x*(-6.d0*x2*y2+y4+x4+z2*x2-
     x             3.d0*z2*y2)
            fz1 =  invrc*ypre4m4*16.d0*x*y*z*(x2-y2)
            
            fx = fx + prefactor2*(f_ij*fx1 + fx2*ImYlm)
            fy = fy + prefactor2*(f_ij*fy1 + fy2*ImYlm)
            fz = fz + prefactor2*(f_ij*fz1 + fz2*ImYlm)
            
c--------------------------------------------------------
c     Real and imaginary spherical harmonics for m = -3  
c--------------------------------------------------------
            
            ReYlm = ypre4m3*invrs*z*(x3-3.d0*x*y2)
            ImYlm = ypre4m3*invrs*z*(-3.d0*x2*y+y3)
            
c--------------------------------------------------
c     Force contributions from m = -3 (real part)  
c--------------------------------------------------
            
            prefactor2 = 2.0d0*q4prefactor(itype)*ReQ4bar(-3,itype)
            
c---------------------------------------------
c     Gradient of Re(Y_{4,-3}) w.r.t r_{j}
c---------------------------------------------
            
            fx1 = -invrc*ypre4m3*z*(x4-12.d0*x2*y2-3.d0*z2*x2+
     x             3.d0*y4+3.d0*z2*y2)
            fy1 = -invrc*ypre4m3*2.d0*z*x*y*(5.0d0*x2-3.d0*y2+
     x             3.d0*z2)
            fz1 =  invrc*ypre4m3*x*(x2-3.d0*y2)*(x2+y2-3.d0*z2)

            fx = fx + prefactor2*(f_ij*fx1 + fx2*ReYlm)
            fy = fy + prefactor2*(f_ij*fy1 + fy2*ReYlm)
            fz = fz + prefactor2*(f_ij*fz1 + fz2*ReYlm)
            
c-------------------------------------------------------
c     Force contributions from m = -3 (imaginary part) 
c-------------------------------------------------------

            prefactor2 = 2.0d0*q4prefactor(itype)*ImQ4bar(-3,itype)
            
c--------------------------------------------------
c     Gradient of Im(Y_{4,-3}) w.r.t r_{j}
c--------------------------------------------------

            fx1 =  invrc*ypre4m3*2.d0*z*x*y*(3.d0*x2-5.d0*y2-3.0d0*z2)
            fy1 = -invrc*ypre4m3*z*(-12.d0*x2*y2+y4+3.d0*x4+3.d0*z2*x2
     x            -3.d0*z2*y2)
            fz1 = -invrc*ypre4m3*y*(3.d0*x2-y2)*(x2+y2-3.d0*z2)
            
            fx = fx + prefactor2*(f_ij*fx1 + fx2*ImYlm)
            fy = fy + prefactor2*(f_ij*fy1 + fy2*ImYlm)
            fz = fz + prefactor2*(f_ij*fz1 + fz2*ImYlm)
            
c--------------------------------------------------------
c     Real and imaginary spherical harmonics for m = -2  
c--------------------------------------------------------
            
            ReYlm = -ypre4m2*invrs*(x2-y2)*(-6.d0*z2+x2+y2)
            ImYlm = ypre4m2*invrs*2.d0*(-6.d0*z2+x2+y2)*x*y
            
c--------------------------------------------------
c     Force contributions from m = -2 (real part)  
c--------------------------------------------------
            
            prefactor2 = 2.0d0*q4prefactor(itype)*ReQ4bar(-2,itype)
            
c---------------------------------------------
c     Gradient of Re(Y_{4,-2}) w.r.t r_{j}
c---------------------------------------------

            fx1 = -invrc*ypre4m2*4.d0*x*(4.d0*z2*x2+y4-9.d0*z2*y2-
     x             3.d0*z4+x2*y2)
            fy1 =  invrc*ypre4m2*4.d0*y*(x4-9.d0*z2*x2+4.d0*z2*y2-
     x             3.d0*z4+x2*y2)
            fz1 =  invrc*ypre4m2*4.d0*z*(x2-y2)*(4.d0*x2+4.d0*y2-
     x             3.d0*z2)
            
            fx = fx + prefactor2*(f_ij*fx1 + fx2*ReYlm)
            fy = fy + prefactor2*(f_ij*fy1 + fy2*ReYlm)
            fz = fz + prefactor2*(f_ij*fz1 + fz2*ReYlm)
            
c-------------------------------------------------------
c     Force contributions from m = -2 (imaginary part) 
c-------------------------------------------------------
            
            prefactor2 = 2.0d0*q4prefactor(itype)*ImQ4bar(-2,itype)
            
c-----------------------------------------
c Gradient of Im(Y_{4,-2}) w.r.t r_{j}
c-----------------------------------------
            
            fx1 = -invrc*ypre4m2*2.d0*y*(x4-21.d0*z2*x2+5.d0*z2*y2+
     x             6.d0*z4-y4)
            fy1 =  invrc*ypre4m2*2.d0*x*(-y4+21.d0*z2*y2-5.d0*z2*x2
     x            -6.d0*z4+x4)
            fz1 = -invrc*ypre4m2*8.d0*z*x*y*(4.d0*x2+4.d0*y2-3.d0*z2)
            
            fx = fx + prefactor2*(f_ij*fx1 + fx2*ImYlm)
            fy = fy + prefactor2*(f_ij*fy1 + fy2*ImYlm)
            fz = fz + prefactor2*(f_ij*fz1 + fz2*ImYlm)
            
c--------------------------------------------------------
c     Real and imaginary spherical harmonics for m = -1
c--------------------------------------------------------
            
            ReYlm = -ypre4m1*invrs*z*(-4.d0*z2+3.d0*x2+3.d0*y2)*x
            ImYlm = ypre4m1*invrs*z*(-4.d0*z2+3.d0*x2+3.d0*y2)*y
            
c--------------------------------------------------
c     Force contributions from m = -1 (real part)  
c--------------------------------------------------
            
            prefactor2 = 2.0d0*q4prefactor(itype)*ReQ4bar(-1,itype)
            
c----------------------------------------
c Gradient of Re(Y_{4,-1}) w.r.t r_{j}
c----------------------------------------
            
            fx1 =  invrc*ypre4m1*z*(3.d0*x4-21.d0*z2*x2+z2*y2+4.d0*z4-
     x             3.d0*y4)
            fy1 =  invrc*ypre4m1*2.d0*z*x*y*(3.d0*x2+3.d0*y2-11.d0*z2)
            fz1 = -invrc*ypre4m1*x*(-21.d0*x2*z2-21.d0*z2*y2+4.d0*z4+
     x             3.d0*x4+6.d0*x2*y2+3.d0*y4)
            
            fx = fx + prefactor2*(f_ij*fx1 + fx2*ReYlm)
            fy = fy + prefactor2*(f_ij*fy1 + fy2*ReYlm)
            fz = fz + prefactor2*(f_ij*fz1 + fz2*ReYlm)
            
c-------------------------------------------------------
c     Force contributions from m = -1 (imaginary part) 
c-------------------------------------------------------
            
            prefactor2 = 2.0d0*q4prefactor(itype)*ImQ4bar(-1,itype)
            
c--------------------------------------------
c     Gradient of Im(Y_{4,-1}) w.r.t r_{j}
c--------------------------------------------
            
            fx1 = -invrc*ypre4m1*2.d0*z*x*y*(3.d0*x2+3.d0*y2-11.d0*z2)
            fy1 =  invrc*ypre4m1*z*(-3.d0*y4+21.d0*z2*y2-z2*x2-4.d0*z4+
     x             3.d0*x4)
            fz1 =  invrc*ypre4m1*y*(-21.d0*z2*x2-21.d0*z2*y2+4.d0*z4+
     x             3.d0*x4+6.d0*x2*y2+3.d0*y4)
            
            fx = fx + prefactor2*(f_ij*fx1 + fx2*ImYlm)
            fy = fy + prefactor2*(f_ij*fy1 + fy2*ImYlm)
            fz = fz + prefactor2*(f_ij*fz1 + fz2*ImYlm)
            
c--------------------------------------------------------
c     Real spherical harmonics for m = 0                 
c--------------------------------------------------------
            
            ReYlm =  ypre4m0*invrs*(8.d0*z4-24.d0*z2*x2-24.d0*z2*y2+
     x        3.d0*x4+6.d0*x2*y2+3.d0*y4)
            
c--------------------------------------------------
c     Force contributions from m = 0 (real part)   
c--------------------------------------------------
            
            prefactor2 = q4prefactor(itype)*ReQ4bar(0,itype)
            
c-------------------------------------------
c     Gradient of Re(Y_{4,0}) w.r.t r_{j}
c-------------------------------------------
            
            fx1 =  20.d0*ypre4m0*invrc*z2*(-4.d0*z2+3.d0*x2+3.d0*y2)*x
            fy1 =  20.d0*ypre4m0*invrc*z2*(-4.d0*z2+3.d0*x2+3.d0*y2)*y
            fz1 = -20.d0*ypre4m0*invrc*z*(-4.d0*z2*x2-4.d0*z2*y2+3.d0*
     x             x4+6.d0*x2*y2+3.d0*y4)
            
            fx = fx + prefactor2*(f_ij*fx1 + fx2*ReYlm)
            fy = fy + prefactor2*(f_ij*fy1 + fy2*ReYlm)
            fz = fz + prefactor2*(f_ij*fz1 + fz2*ReYlm)
            
c     Add into global force and stress arrays
            
            fxx(jatm) = fxx(jatm) + fx
            fyy(jatm) = fyy(jatm) + fy
            fzz(jatm) = fzz(jatm) + fz
            fxx(iatm) = fxx(iatm) - fx
            fyy(iatm) = fyy(iatm) - fy
            fzz(iatm) = fzz(iatm) - fz
            
c     Virial
            
            virord = virord - (fx*x + fy*y + fz*z)
            
c     Stress
            
            strs1=strs1+x*fx
            strs2=strs2+x*fy
            strs3=strs3+x*fz
            strs5=strs5+y*fy
            strs6=strs6+y*fz
            strs9=strs9+z*fz
            
          end do                ! end loop over connection list for iatm
          
        end if                  ! end of computing Q4
        
c---------------------------------------------------------
c---------------------------------------------------------
c     Compute forces arising from  Q6 order parameters        
c---------------------------------------------------------
c---------------------------------------------------------
        
        if ( (nq6>0).and.isolvmax6>0 ) then
          
          do isolv6 = 1,isolvmax6
            
            jatm = solvlist6(isolv6)
            itype = solvtype6(isolv6)
            
            invrc = solvimag6(isolv6)**6
            invrq = solvimag6(isolv6)**8
            invrs = solvimag6(isolv6)**4
            
            x     = solvx6(isolv6)
            y     = solvy6(isolv6)
            z     = solvz6(isolv6)
            
            f_ij  =  fc(solvrmag6(isolv6),q6cutoff(1,itype),
     x                 q6cutoff(2,itype))
            df_ij = dfc(solvrmag6(isolv6),q6cutoff(1,itype),
     x                 q6cutoff(2,itype))
            
            x2 = x*x
            y2 = y*y
            z2 = z*z
            
            x3 = x2*x
            y3 = y2*y
            z3 = z2*z
            
            x4 = x2*x2
            y4 = y2*y2
            z4 = z2*z2
            
            x5 = x4*x
            y5 = y4*y
            z5 = z4*z
            
            x6 = x4*x2
            y6 = y4*y2
            z6 = z4*z2
            
            fx = 0.0d0
            fy = 0.0d0
            fz = 0.0d0
            
c----------------------------------------
c     Gradient of f_ij w.r.t. r_{j}  
c----------------------------------------
            
            fx2 = df_ij*x*solvimag6(isolv6)
            fy2 = df_ij*y*solvimag6(isolv6)
            fz2 = df_ij*z*solvimag6(isolv6)
            
c-----------------------------------------------------------
c     Real and imaginary spherical harmonics for m = -6     
c-----------------------------------------------------------
            
            ReYlm      = ypre6m6*invrc*(x6-15.0d0*x4*y2+15.0d0*x2*y4-y6)
            ImYlm      = ypre6m6*invrc*(-6.0d0*x5*y+20.0d0*x3*y3-6.0d0*
     x                   x*y5)
            
c-----------------------------------------------------
c     Force contributions from m = -6 (real part)     
c-----------------------------------------------------
            
            prefactor2 = 2.0d0*q6prefactor(itype)*ReQ6bar(-6,itype)
            
c--------------------------------------------
c     Gradient of Re(Y_{6,-6}) w.r.t r_{j}
c--------------------------------------------
            
            fx1 =  invrq*ypre6m6*6.0d0*x*(6.0d0*x4*y2-20.0d0*x2*y4+
     x             6.0d0*y6+z2*x4-10.0d0*z2*x2*y2+5.0d0*z2*y4)
            fy1 = -invrq*ypre6m6*6.0d0*y*(6.0d0*x6-20.0d0*x4*y2+6.0d0
     x             *x2*y4+5.0d0*z2*x4-10.0d0*z2*x2*y2+z2*y4)
            fz1 = -invrq*ypre6m6*6.0d0*z*(x6-15.0d0*x4*y2+15.0d0*x2*
     x             y4-y6)
            
            fx = fx + prefactor2*(f_ij*fx1 + fx2*ReYlm)
            fy = fy + prefactor2*(f_ij*fy1 + fy2*ReYlm)
            fz = fz + prefactor2*(f_ij*fz1 + fz2*ReYlm)
            
c----------------------------------------------------------
c     Force contributions from m = -6 (Imaginary part)     
c----------------------------------------------------------
            
            prefactor2 = 2.0d0*q6prefactor(itype)*ImQ6bar(-6,itype)
            
c--------------------------------------------
c     Gradient of Im(Y_{6,-6}) w.r.t r_{j}
c--------------------------------------------
            
            fx1 = -invrq*ypre6m6*6.0d0*y*(-x6+15.0d0*x4*y2-15.0d0*x2*
     x             y4+5.0d0*z2*x4-10.0d0*z2*x2*y2+y6+z2*y4)
            fy1 = -invrq*ypre6m6*6.0d0*x*(-15.0d0*x4*y2+15.0d0*x2*y4-
     x             y6+x6+z2*x4-10.0d0*z2*x2*y2+5.0d0*z2*y4)
            fz1 =  invrq*ypre6m6*12.0d0*x*y*z*(3.0d0*x4-10.0d0*x2*y2+
     x             3.0d0*y4)
            
            fx = fx + prefactor2*(f_ij*fx1 + fx2*ImYlm)
            fy = fy + prefactor2*(f_ij*fy1 + fy2*ImYlm)
            fz = fz + prefactor2*(f_ij*fz1 + fz2*ImYlm)
            
c---------------------------------------------------------
c     Real and imaginary spherical harmonics for m = -5   
c---------------------------------------------------------
            
            ReYlm = -ypre6m5*invrc*z*(-x5+10.0d0*x3*y2-5.0d0*x*y4)
            ImYlm = -ypre6m5*invrc*z*(5.0d0*x4*y-10.0d0*x2*y3+y5)
            
c--------------------------------------------------
c     Force contributions from m = -5 (real part)  
c--------------------------------------------------
            
            prefactor2 = 2.0d0*q6prefactor(itype)*ReQ6bar(-5,itype)
            
c--------------------------------------------
c     Gradient of Re(Y_{6,-5}) w.r.t r_{j}
c--------------------------------------------
            
            fx1 = ypre6m5*invrq*z*(-x6+35.0d0*x4*y2-55.d0*x2*y4+
     x            5.0d0*z2*x4-30.0d0*z2*x2*y2+5.0d0*y6+5.0d0*z2*y4)
            fy1 = -2.0d0*ypre6m5*invrq*x*y*z*(13.0d0*x4-30.0d0*x2
     x        *y2+5.0d0*y4+10.0d0*z2*x2-10.0d0*z2*y2) 
            fz1 = -ypre6m5*invrq*x*(x4-10.0d0*x2*y2+5.0d0*y4)*
     x        (-x2-y2+5.0d0*z2)
            
            fx = fx + prefactor2*(f_ij*fx1 + fx2*ReYlm)
            fy = fy + prefactor2*(f_ij*fy1 + fy2*ReYlm)
            fz = fz + prefactor2*(f_ij*fz1 + fz2*ReYlm)
            
c--------------------------------------------------------
c     Force contributions from m = -5 (Imaginary part)   
c--------------------------------------------------------

            prefactor2 = 2.0d0*q6prefactor(itype)*ImQ6bar(-5,itype)
            
c--------------------------------------------
c     Gradient of Im(Y_{6,-5}) w.r.t r_{j}
c--------------------------------------------
            
            fx1 = -2.0d0*ypre6m5*invrq*x*y*z*(-5.0d0*x4+30.0d0*x2*y2
     x            -13.0d0*y4+10.0d0*z2*x2-10.0d0*z2*y2)
            fy1 = -ypre6m5*invrq*z*(-55.0d0*x4*y2+35.0d0*x2*y4-y6+
     x             5.0d0*x6+5.0d0*z2*x4-30.0d0*z2*x2*y2+5.0d0*z2*y4)
            fz1 =  ypre6m5*invrq*y*(5.0d0*x4-10.0d0*x2*y2+y4)*
     x             (-x2-y2+5.0d0*z2)
            
            fx = fx + prefactor2*(f_ij*fx1 + fx2*ImYlm)
            fy = fy + prefactor2*(f_ij*fy1 + fy2*ImYlm)
            fz = fz + prefactor2*(f_ij*fz1 + fz2*ImYlm)
            
c--------------------------------------------------------
c     Real and imaginary spherical harmonics for m = -4  
c--------------------------------------------------------
            
            ReYlm = ypre6m4*invrc*(10.0d0*z2-x2-y2)*(x4-6.0d0*x2*y2+y4)
            ImYlm = ypre6m4*invrc*(10.0d0*z2-x2-y2)*(-4.0d0*x3*y+4.0d0*
     x              x*y3)
            
c--------------------------------------------------
c     Force contributions from m = -4 (real part)  
c--------------------------------------------------
            
            prefactor2 = 2.0d0*q6prefactor(itype)*ReQ6bar(-4,itype)
            
c--------------------------------------------
c     Gradient of Re(Y_{6,-4}) w.r.t r_{j}
c--------------------------------------------
            
            fx1 = 2.0d0*ypre6m4*invrq*x*(-8.0d0*x4*y2-13.0d0*z2*x4+
     x            150.0d0*z2*x2*y2+8.0d0*y6-85.0d0*z2*y4+20.0d0*z4*
     x            x2-60.0d0*z4*y2)
            fy1 =-2.0d0*ypre6m4*invrq*y*(-8.0d0*x6+85.0d0*z2*x4+8.0d0
     x           *x2*y4-150.0d0*z2*x2*y2+13.0d0*z2*y4+60.0d0*z4*x2-
     x            20.0d0*z4*y2)
            fz1 =-2.0d0*ypre6m4*invrq*z*(x4-6.0d0*x2*y2+y4)*(-13.0d0
     x           *x2-13.0d0*y2+20.0d0*z2)
            
            fx = fx + prefactor2*(f_ij*fx1 + fx2*ReYlm)
            fy = fy + prefactor2*(f_ij*fy1 + fy2*ReYlm)
            fz = fz + prefactor2*(f_ij*fz1 + fz2*ReYlm)
            
c--------------------------------------------------------
c     Force contributions from m = -4 (Imaginary part)   
c--------------------------------------------------------
            
            prefactor2 = 2.0d0*q6prefactor(itype)*ImQ6bar(-4,itype)
            
c--------------------------------------------
c     Gradient of Im(Y_{6,-4}) w.r.t r_{j}
c--------------------------------------------
            
            fx1 = -4.0d0*ypre6m4*invrq*y*(x6-35.0d0*z2*x4-5.0d0*x2*y4
     x             +80.0d0*z2*x2*y2+30*z4*x2-9.0d0*z2*y4-10*z4*y2-5.0d0
     x             *x4*y2+y6)
            fy1 = -4.0d0*ypre6m4*invrq*x*(5.0d0*x4*y2-80.0d0*z2*x2*y2
     x            -y6+35.0d0*z2*y4+9.0d0*z2*x4+10.0d0*z4*x2-30.0d0*z4
     x            *y2-x6+5.0d0*x2*y4)
            fz1 =  8.0d0*ypre6m4*invrq*z*x*y*(x2-y2)*(-13.0d0*x2-
     x             13.0d0*y2+20.0d0*z2)
            
            fx = fx + prefactor2*(f_ij*fx1 + fx2*ImYlm)
            fy = fy + prefactor2*(f_ij*fy1 + fy2*ImYlm)
            fz = fz + prefactor2*(f_ij*fz1 + fz2*ImYlm)
            
c--------------------------------------------------------
c     Real and imaginary spherical harmonics for m = -3  
c--------------------------------------------------------

            ReYlm = -ypre6m3*invrc*z*(8.0d0*z2-3.0d0*x2-3.0d0*y2)*
     x               (-x3+3.0d0*x*y2)
            ImYlm = -ypre6m3*invrc*z*(8.0d0*z2-3.0d0*x2-3.0d0*y2)*
     x               (3.0d0*x2*y-y3)
            
c--------------------------------------------------
c     Force contributions from m = -3 (real part)  
c--------------------------------------------------
            
            prefactor2 = 2.0d0*q6prefactor(itype)*ReQ6bar(-3,itype)
            
c--------------------------------------------
c     Gradient of Re(Y_{6,-3}) w.r.t r_{j}
c--------------------------------------------
            
            fx1 =  3.d0*ypre6m3*invrq*z*(x6-11.d0*x4*y2-13.d0*z2*x4-
     x             9.d0*x2*y4+54.d0*z2*x2*y2+8.d0*z4*x2-5.d0*z2*y4-
     x             8.d0*z4*y2+3.d0*y6)
            fy1 = -6.d0*ypre6m3*invrq*z*x*y*(-5.d0*x4-2.d0*x2*y2+14.d0
     x            *z2*x2+3.d0*y4-22.d0*z2*y2+8.d0*z4)
            fz1 = -3.d0*ypre6m3*invrq*x*(x2-3.d0*y2)*(-13.d0*z2*x2-
     x             13.d0*z2*y2+8.d0*z4+x4+2.d0*x2*y2+y4)
            
            fx = fx + prefactor2*(f_ij*fx1 + fx2*ReYlm)
            fy = fy + prefactor2*(f_ij*fy1 + fy2*ReYlm)
            fz = fz + prefactor2*(f_ij*fz1 + fz2*ReYlm)
            
c--------------------------------------------------------
c     Force contributions from m = -3 (Imaginary part)   
c--------------------------------------------------------
            
            prefactor2 = 2.0d0*q6prefactor(itype)*ImQ6bar(-3,itype)
            
c--------------------------------------------
c     Gradient of Im(Y_{6,-3}) w.r.t r_{j}
c--------------------------------------------
            
            fx1 = -6.d0*ypre6m3*invrq*z*x*y*(3.d0*x4-2.d0*x2*y2-22.d0
     x            *z2*x2-5.d0*y4+14.d0*z2*y2+8.d0*z4)
            fy1 = -3.d0*ypre6m3*invrq*z*(9.d0*x4*y2+11.d0*x2*y4-54.d0
     x            *z2*x2*y2-y6+13.d0*z2*y4+5.d0*z2*x4+8.d0*z4*x2-8.d0
     x            *z4*y2-3.d0*x6)
            fz1 =  3.d0*ypre6m3*invrq*y*(3.d0*x2-y2)*(-13.d0*z2*x2-
     x             13.d0*z2*y2+8.d0*z4+x4+2.d0*x2*y2+y4)
            
            fx = fx + prefactor2*(f_ij*fx1 + fx2*ImYlm)
            fy = fy + prefactor2*(f_ij*fy1 + fy2*ImYlm)
            fz = fz + prefactor2*(f_ij*fz1 + fz2*ImYlm)
            
c--------------------------------------------------------
c     Real and imaginary spherical harmonics for m = -2  
c--------------------------------------------------------
            
            ReYlm =  ypre6m2*invrc*(16.0d0*z4-16.0d0*z2*x2-16.0d0*z2
     x              *y2+x4+2.0d0*x2*y2+y4)*(x2-y2)
            ImYlm = -ypre6m2*invrc*2.0d0*(16.0d0*z4-16.0d0*z2*x2-
     x               16.0d0*z2*y2+x4+2.0d0*x2*y2+y4)*x*y
            
c--------------------------------------------------
c     Force contributions from m = -2 (real part)  
c--------------------------------------------------
            
            prefactor2 = 2.0d0*q6prefactor(itype)*ReQ6bar(-2,itype)
            
c--------------------------------------------
c     Gradient of Re(Y_{6,-2}) w.r.t r_{j}
c--------------------------------------------
            
            fx1 =  2.d0*ypre6m2*invrq*x*( 19.d0*z2*x4-64.d0*z4*x2-
     x             49.d0*z2*y4+64.d0*z4*y2+2.d0*x4*y2+4.d0*x2*y4+
     x             2.d0*y6+16.d0*z6-30.d0*z2*x2*y2)
            fy1 = -2.d0*ypre6m2*invrq*y*(-49.d0*z2*x4+64.d0*z4*x2+
     x             19.d0*z2*y4-64.d0*z4*y2+2.d0*x6+4.d0*x4*y2+2.d0*
     x             x2*y4+16.d0*z6-30.d0*z2*x2*y2)
            fz1 = -2.d0*ypre6m2*invrq*z*(x2-y2)*(-64.d0*z2*x2-64.d0
     x            *z2*y2+16.d0*z4+19.d0*x4+38.d0*x2*y2+19.d0*y4)
            
            fx = fx + prefactor2*(f_ij*fx1 + fx2*ReYlm)
            fy = fy + prefactor2*(f_ij*fy1 + fy2*ReYlm)
            fz = fz + prefactor2*(f_ij*fz1 + fz2*ReYlm)
            
c--------------------------------------------------------
c     Force contributions from m = -2 (Imaginary part)   
c--------------------------------------------------------
            
            prefactor2 = 2.0d0*q6prefactor(itype)*ImQ6bar(-2,itype)
            
c--------------------------------------------
c     Gradient of Im(Y_{6,-2}) w.r.t r_{j}
c--------------------------------------------
            
            fx1 = -2.d0*ypre6m2*invrq*y*(53.d0*z2*x4+38.d0*z2*x2*
     x             y2-128.d0*z4*x2-x6-x4*y2+x2*y4+16.d0*z6-15.d0*
     x             z2*y4+y6)
            fy1 = -2.d0*ypre6m2*invrq*x*(38.d0*z2*x2*y2+53.d0*z2*
     x             y4-128.d0*z4*y2+x4*y2-x2*y4-y6+16.d0*z6-15.d0*
     x             z2*x4+x6)
            fz1 =  4.d0*ypre6m2*invrq*z*x*y*(-64.d0*z2*x2-64.d0*
     x             z2*y2+16.d0*z4+19.d0*x4+38.d0*x2*y2+19.d0*y4)
            
            fx = fx + prefactor2*(f_ij*fx1 + fx2*ImYlm)
            fy = fy + prefactor2*(f_ij*fy1 + fy2*ImYlm)
            fz = fz + prefactor2*(f_ij*fz1 + fz2*ImYlm)
            
c--------------------------------------------------------
c     Real and imaginary spherical harmonics for m = -1  
c--------------------------------------------------------
            
            ReYlm =  ypre6m1*z*invrc*(8.0d0*z4-20.0d0*z2*x2-20.0d0*z2*y2
     x               +5.0d0*x4+10.0d0*x2*y2+5.0d0*y4)*x
            ImYlm = -ypre6m1*z*invrc*(8.0d0*z4-20.0d0*z2*x2-20.0d0*z2*y2
     x               +5.0d0*x4+10.0d0*x2*y2+5.0d0*y4)*y
            
c--------------------------------------------------
c     Force contributions from m = -1 (real part)  
c--------------------------------------------------
            
            prefactor2 = 2.0d0*q6prefactor(itype)*ReQ6bar(-1,itype)
            
c--------------------------------------------
c     Gradient of Re(Y_{6,-1}) w.r.t r_{j}
c--------------------------------------------
            
            fx1 =       ypre6m1*invrq*z*(85.d0*z2*x4+70.d0*z2*x2*y2
     x            -100.d0*z4*x2-5.d0*x6-5.d0*x4*y2+5.d0*x2*y4-12.d0
     x            *z4*y2+8.d0*z6-15.d0*z2*y4+5.d0*y6)
            fy1 = -2.d0*ypre6m1*invrq*z*x*y*(-50.d0*z2*x2-50.d0*z2*
     x             y2+44.d0*z4+5.d0*x4+10.d0*x2*y2+5.d0*y4)
            fz1 =      -ypre6m1*invrq*x*(-100.d0*z4*x2-100.d0*z4*y2+
     x             8.d0*z6+85.d0*z2*x4+170.d0*z2*x2*y2+85.d0*z2*y4-
     x             5.d0*x6-15.d0*x4*y2-15.d0*x2*y4-5.d0*y6)
            
            fx = fx + prefactor2*(f_ij*fx1 + fx2*ReYlm)
            fy = fy + prefactor2*(f_ij*fy1 + fy2*ReYlm)
            fz = fz + prefactor2*(f_ij*fz1 + fz2*ReYlm)
            
c--------------------------------------------------------
c     Force contributions from m = -1 (Imaginary part)   
c--------------------------------------------------------
            
            prefactor2 = 2.0d0*q6prefactor(itype)*ImQ6bar(-1,itype)
            
c--------------------------------------------
c     Gradient of Im(Y_{6,-1}) w.r.t r_{j}
c--------------------------------------------
            
            fx1 =  2.d0*ypre6m1*invrq*z*x*y*(-50.d0*z2*x2-50.d0*
     x             z2*y2+44.d0*z4+5.d0*x4+10.d0*x2*y2+5.d0*y4)
            fy1 =      -ypre6m1*invrq*z*(70.d0*z2*x2*y2+85.d0*z2
     x            *y4-100.d0*z4*y2+5.d0*x4*y2-5.d0*x2*y4-5.d0*y6
     x            -12.d0*z4*x2+8.d0*z6-15.d0*z2*x4+5.d0*x6)
            fz1 =       ypre6m1*invrq*y*(-100.d0*z4*x2-100.d0*z4
     x            *y2+8.d0*z6+85.d0*z2*x4+170.d0*z2*x2*y2+85.d0*
     x             z2*y4-5.d0*x6-15.d0*x4*y2-15.d0*x2*y4-5.d0*y6)
            
            fx = fx + prefactor2*(f_ij*fx1 + fx2*ImYlm)
            fy = fy + prefactor2*(f_ij*fy1 + fy2*ImYlm)
            fz = fz + prefactor2*(f_ij*fz1 + fz2*ImYlm)
            
c--------------------------------------------------------
c     Real spherical harmonics for m = 0                 
c--------------------------------------------------------
            
            ReYlm = ypre6m0*invrc*(16.0d0*z6-120.0d0*z4*x2-120.0d0
     x             *z4*y2+90.0d0*z2*x4+180.0d0*z2*x2*y2+90.0d0*z2*y4
     x             -5.0d0*x6-15.0d0*x4*y2-15.0d0*x2*y4-5.0d0*y6)
            
c--------------------------------------------------
c     Force contributions from m = 0 (real part)       
c--------------------------------------------------
            
            prefactor2 = q6prefactor(itype)*ReQ6bar(0,itype)
            
c--------------------------------------------
c     Gradient of Re(Y_{6,0}) w.r.t r_{j}
c--------------------------------------------
            
            fx1 = -42.d0*ypre6m0*invrq*z2*(8.d0*z4-20.d0*z2*x2-20.d0
     x            *z2*y2+5.d0*x4+10.d0*x2*y2+5.d0*y4)*x
            fy1 = -42.d0*ypre6m0*invrq*z2*(8.d0*z4-20.d0*z2*x2-20.d0
     x            *z2*y2+5.d0*x4+10.d0*x2*y2+5.d0*y4)*y
            fz1 =  42.d0*ypre6m0*invrq*z*(8.d0*z4*x2+8.d0*z4*y2-20.d0
     x            *z2*x4-40.d0*z2*x2*y2-20.d0*z2*y4+5.d0*x6+15.d0*x4*
     x             y2+15.d0*x2*y4+5.d0*y6)
            
            fx = fx + prefactor2*(f_ij*fx1 + fx2*ReYlm)
            fy = fy + prefactor2*(f_ij*fy1 + fy2*ReYlm)
            fz = fz + prefactor2*(f_ij*fz1 + fz2*ReYlm)
            
            
c     Add into global force and stress arrays
            
            fxx(jatm) = fxx(jatm) + fx
            fyy(jatm) = fyy(jatm) + fy
            fzz(jatm) = fzz(jatm) + fz
            fxx(iatm) = fxx(iatm) - fx
            fyy(iatm) = fyy(iatm) - fy
            fzz(iatm) = fzz(iatm) - fz
            
c     Virial
            
            virord = virord - (fx*x + fy*y + fz*z)
            
c     Stress
            
            strs1=strs1+x*fx
            strs2=strs2+x*fy
            strs3=strs3+x*fz
            strs5=strs5+y*fy
            strs6=strs6+y*fz
            strs9=strs9+z*fz
            
          end do                ! end loop over connection list for iatm
          
        end if                  ! end of computing Q6
        
        
      end do                    ! end loop over iatm
      
c     complete stress tensor
      
      stress(1)=stress(1)+strs1
      stress(2)=stress(2)+strs2
      stress(3)=stress(3)+strs3
      stress(4)=stress(4)+strs2
      stress(5)=stress(5)+strs5
      stress(6)=stress(6)+strs6
      stress(7)=stress(7)+strs3
      stress(8)=stress(8)+strs6
      stress(9)=stress(9)+strs9
      
c     tidy up
      
      ierr = 0
      
      deallocate(xdf,stat=ierr(1))
      deallocate(ydf,stat=ierr(2))
      deallocate(zdf,stat=ierr(3))

      deallocate(dstlst,stat=ierr(4))
      
      deallocate(solvx4,stat=ierr(5))
      deallocate(solvy4,stat=ierr(6))
      deallocate(solvz4,stat=ierr(7))
      deallocate(solvrmag4,stat=ierr(8))
      deallocate(solvimag4,stat=ierr(9))
      deallocate(solvrsq4 ,stat=ierr(10))
      deallocate(solvlist4,stat=ierr(11)) 
      deallocate(solvtype4,stat=ierr(12))
      
      deallocate(solvx6,stat=ierr(13))
      deallocate(solvy6,stat=ierr(14))
      deallocate(solvz6,stat=ierr(15))
      deallocate(solvrmag6,stat=ierr(16))
      deallocate(solvimag6,stat=ierr(17))
      deallocate(solvrsq6 ,stat=ierr(18))
      deallocate(solvlist6,stat=ierr(19)) 
      deallocate(solvtype6,stat=ierr(20))
      if (any(ierr/=0)) call Mfrz_Error(2536,0.d0)
      
      deallocate(buff1,stat=ierr(1))
      deallocate(buff2,stat=ierr(2))    
      if (any(ierr/=0)) call Mfrz_Error(2537,0.d0)
      
      deallocate(q4prefactor,stat=ierr(1))
      deallocate(q6prefactor,stat=ierr(2))
      if (any(ierr/=0)) call Mfrz_Error(2540,0.d0)
      
      return
      
      end Subroutine Compute_Steinhardt_Forces
      
      Subroutine Compute_Tet_Nlist(imcon,natms)
      
c---------------------------------------------------------------------
c
c     The existing neighbour list is not known by all nodes and 
c     therefore we compute a new one from scratch rather than trying
c     to merge a full neighbour list across all MPI tasks or restore
c     symmetry.
c     
c     Author D. Quigley - University of Warwick
c     Copyright D. Quigley October 2008
c
c---------------------------------------------------------------------
      
      use config_module
      use site_module
      use utility_module, only : images
      
      implicit none
      integer,intent(in) :: imcon,natms
      integer :: nlast,megatm
      
      integer :: iatm0,iatm1,istrd,iatm,jatm,nn,itype,jtype
      integer :: ninclude,ifi,ii,ztype,nnn,k
      real(8) :: rsq,rangesq
      real(8),allocatable,dimension(:) :: xdf,ydf,zdf
      
      iatm0 = myrank+1
      iatm1 = natms
      istrd = commsize
      
      ninclude = maxval(zetaninc)
      
      nlast  = natms
      megatm = natms
      
      nnn = int(dble(ninclude*nlast*1.2)/dble(megatm))
      
      allocate(xdf(1:nnn),stat=ierr(1))
      allocate(ydf(1:nnn),stat=ierr(2))
      allocate(zdf(1:nnn),stat=ierr(3))
      if (any(ierr/=0)) call Mfrz_Error(2541,0.d0) 
      
      ii = 0
      do iatm = iatm0,iatm1,istrd
        
        itype = ltype(iatm)
        
        ztype = zetasite(itype)
        if ( ztype==0 ) cycle
        
        nn  = 0                 ! Number of images to compute
        ii  = ii + 1            ! index for this list
        ifi = 0                 ! index for entries in this list
        
        if (ii>mxninc) call Mfrz_Error(2542,0.d0)
        
        do jatm = 1,nlast
          
          jtype = ltype(jatm)
          
          if ( itype/=jtype ) cycle
          if ( iatm == jatm ) cycle
          
          nn = nn + 1
          xdf(nn) = xxx(iatm) - xxx(jatm) ! separation vector
          ydf(nn) = yyy(iatm) - yyy(jatm)
          zdf(nn) = zzz(iatm) - zzz(jatm)
          
        end do
        
        if ( nn > nnn ) call Mfrz_Error(2543,0.d0)
        
        call images(imcon,0,1,nn,cell,xdf,ydf,zdf)
        
        nn = 0
        do jatm = 1,nlast
          
          jtype = ltype(jatm)
          
          if ( itype/=jtype ) cycle
          if ( iatm == jatm ) cycle
          
          nn = nn + 1
          rsq = xdf(nn)*xdf(nn) + ydf(nn)*ydf(nn) + zdf(nn)*zdf(nn)
          
          rangesq =  zetacutoff(2,ztype)**2
          
          if ( rsq < rangesq ) then
            
            ifi = ifi + 1
            flist(ifi,ii) = jatm      
            
          end if
          
        end do
        if ( ifi > mxflist ) call Mfrz_Error(2544,0.d0)
        nflist(ii) = ifi
        
      end do                    ! end loop over iatm
      
      deallocate(xdf,stat=ierr(1))
      deallocate(ydf,stat=ierr(2))
      deallocate(zdf,stat=ierr(3))
      if (any(ierr/=0))call Mfrz_Error(2545,0.d0) 
      
      return
      
      end Subroutine Compute_Tet_Nlist
      
      Subroutine Compute_Tetrahedral(imcon,natms)
      
c---------------------------------------------------------------------
c     
c     Author D. Quigley - University of Warwick
c     Copyright D. Quigley October 2008
c     
c---------------------------------------------------------------------
      
      use config_module
      use site_module
      use utility_module, only : images
      
      implicit none
      integer,intent(in) :: imcon,natms
      
      integer :: iatm0,iatm1,istrd,iatm,jatm,nn,itype,jtype
      integer :: ninclude,ifi,ii,k,ztype,mm,katm,ktype,it
      real(8) :: f_ij,f_ik,r_ij,r_ik,x,y,z,dot
      real(8),parameter :: third=1.0d0/3.0d0
      real(8),allocatable,dimension(:) :: xdf,ydf,zdf
      real(8),allocatable,dimension(:) :: buff1
      
      iatm0 = myrank+1
      iatm1 = natms
      istrd = commsize
      
      allocate(xdf(1:mxflist),stat=ierr(1))
      allocate(ydf(1:mxflist),stat=ierr(2))
      allocate(zdf(1:mxflist),stat=ierr(3))
      allocate(buff1(1:ntet),stat=ierr(4))
      if (any(ierr/=0)) call Mfrz_Error(2546,0.d0)
      
      zeta_global(:) = 0.0d0
      
      ii = 0
      do iatm = iatm0,iatm1,istrd
        
        itype = ltype(iatm)
        
c     no tetrahedral contributions from this atom type?
        
        ztype = zetasite(itype)
        if ( ztype == 0 ) cycle
        
        ii  = ii + 1            ! index for this list
        
        do k = 1,nflist(ii)
          
          jatm = flist(k,ii)
          
          xdf(k) = xxx(jatm) - xxx(iatm) ! separation vector
          ydf(k) = yyy(jatm) - yyy(iatm)
          zdf(k) = zzz(jatm) - zzz(iatm)
          
        end do
        
        nn = nflist(ii)
        call images(imcon,0,1,nn,cell,xdf,ydf,zdf)
        
        do k = 1,nflist(ii)
          
          jatm = flist(k,ii)
          
          r_ij = sqrt(xdf(k)*xdf(k) + ydf(k)*ydf(k) + zdf(k)*zdf(k))
          f_ij = fc(r_ij,zetacutoff(1,ztype),zetacutoff(2,ztype)) 
          
          x = xdf(k) ; y = ydf(k) ; z = zdf(k)
          
c     loop over all other entries katm
          
          do mm = k+1,nflist(ii)
            
c     katm also in solvation shell of iatm
            
            katm  = flist(mm,ii)
            
            r_ik = sqrt(xdf(mm)*xdf(mm) + ydf(mm)*ydf(mm) + 
     x        zdf(mm)*zdf(mm))
            f_ik = fc(r_ik,zetacutoff(1,ztype),zetacutoff(2,ztype)) 
            
c     The node holding the central atom keeps the contrib
            
            dot = (x*xdf(mm) + y*ydf(mm) + z*zdf(mm)) / (r_ij*r_ik) +
     x            third
            zeta_global(ztype) = zeta_global(ztype) + f_ij*f_ik*dot*dot
          end do       
        end do
        
      end do
      
      call gdsum(zeta_global(1),ntet,buff1(1:ntet))
      
      do it = 1,ntet
        zeta_global(it) = 1.0d0 - zeta_global(it)/dble(zetaninc(it)*
     x    zetann(it))
      end do
      
      deallocate(xdf ,stat=ierr(1))
      deallocate(ydf ,stat=ierr(2))
      deallocate(zdf ,stat=ierr(3))
      deallocate(buff1,stat=ierr(4))
      if (any(ierr/=0)) call Mfrz_Error(2547,0.d0)
      
      return
      
      end Subroutine Compute_Tetrahedral
      
      Subroutine Compute_Tetrahedral_Forces(imcon,natms,engord,virord)
      
c---------------------------------------------------------------------
c     
c     Author D. Quigley - University of Warwick
c     Copyright D. Quigley October 2008
c     
c---------------------------------------------------------------------
      
      use config_module
      use site_module
      use utility_module, only : images
      
      implicit none
      integer,intent(in) :: imcon,natms
      real(8),intent(inout) :: engord,virord
      
      integer :: iatm,jatm,katm,iatm0,iatm1,istrd,nn
      integer :: ii,k,m,itet,itype,jtype,ztype,it
      
      real(8),parameter :: third = 1.0d0/3.0d0
      real(8) :: strs1,strs2,strs3,strs5,strs6,strs9
      real(8) :: xj,yj,zj,xk,yk,zk,tmpvar,tmpvar2,dot
      real(8) :: r_ij,r_ik,f_ij,f_ik,df_ij,df_ik
      real(8) :: invrij,invrik,ctheta
      real(8) :: fxj,fyj,fzj,fxk,fyk,fzk
      
      real(8),dimension(3) :: rij_hat,rik_hat
      
      real(8),allocatable,dimension(:) :: xdf,ydf,zdf
      real(8),allocatable,dimension(:) :: tetprefactor
      
      iatm0 = myrank+1
      iatm1 = natms
      istrd = commsize
      
      allocate(xdf(1:mxflist),stat=ierr(1))
      allocate(ydf(1:mxflist),stat=ierr(2))
      allocate(zdf(1:mxflist),stat=ierr(3))
      allocate(tetprefactor(1:ntet),stat=ierr(4))
      if (any(ierr/=0)) call Mfrz_Error(2548,0.d0)
      
c     Compute the prefactor
      
      k = nq4+nq6+1
      do it = 1,ntet
        tetprefactor = dcolvar(k)/dble(zetaninc(it)*zetann(it))
        k = k + 1
      end do
      
c     zero contribution to the stress tensor
      
      strs1=0.0d0
      strs2=0.0d0
      strs3=0.0d0
      strs5=0.0d0
      strs6=0.0d0
      strs9=0.0d0
      
      ii = 0
      do iatm = iatm0,iatm1,istrd
        
        itype = ltype(iatm)
        
c     no tetrahedral contributions from this atom type?
        
        ztype = zetasite(itype)
        if ( ztype == 0 ) cycle
        
        ii  = ii + 1            ! index for this list
        
        do k = 1,nflist(ii)
          
          jatm = flist(k,ii)
          
          xdf(k) = xxx(jatm) - xxx(iatm) ! separation vector
          ydf(k) = yyy(jatm) - yyy(iatm)
          zdf(k) = zzz(jatm) - zzz(iatm)
          
        end do
        
        nn = nflist(ii)
        call images(imcon,0,1,nn,cell,xdf,ydf,zdf)
        
        do k = 1,nflist(ii)
          
          jatm  = flist(k,ii)
          
          r_ij  = sqrt(xdf(k)*xdf(k) + ydf(k)*ydf(k) + zdf(k)*zdf(k))
          f_ij  =  fc(r_ij,zetacutoff(1,ztype),zetacutoff(2,ztype))
          df_ij = dfc(r_ij,zetacutoff(1,ztype),zetacutoff(2,ztype))
          
          xj = xdf(k) ; yj = ydf(k) ; zj = zdf(k)
          
          invrij = 1.0d0/r_ij
          
          rij_hat(1) = xj*invrij
          rij_hat(2) = yj*invrij
          rij_hat(3) = zj*invrij
          
          do m = k + 1,nflist(ii)
            
            r_ik = sqrt(xdf(m)*xdf(m) + ydf(m)*ydf(m) + zdf(m)*zdf(m))
            f_ik  =  fc(r_ik,zetacutoff(1,ztype),zetacutoff(2,ztype))
            df_ik = dfc(r_ik,zetacutoff(1,ztype),zetacutoff(2,ztype))
            
            xk = xdf(m) ; yk = ydf(m) ; zk = zdf(m)
            
            invrik = 1.0d0/r_ik
            
            rik_hat(1) = xk*invrik
            rik_hat(2) = yk*invrik
            rik_hat(3) = zk*invrik
            
            ctheta = dot_product(rij_hat,rik_hat)
            dot    = ctheta + third
            
            tmpvar  = 2.0d0*dot*f_ij*f_ik*tetprefactor(ztype)*invrij
            tmpvar2 = tetprefactor(ztype)*dot*dot*df_ij*f_ik
            
c     force between atom i and atom j due to second term 
c     i.e. ( f_ij*f_ik*dot*dot )
            
            fxj =  tmpvar*(rik_hat(1) - rij_hat(1)*ctheta) + 
     x             tmpvar2*rij_hat(1)
            fyj =  tmpvar*(rik_hat(2) - rij_hat(2)*ctheta) + 
     x             tmpvar2*rij_hat(2)
            fzj =  tmpvar*(rik_hat(3) - rij_hat(3)*ctheta) + 
     x             tmpvar2*rij_hat(3)
            
            tmpvar  = 2.0d0*dot*f_ij*f_ik*tetprefactor(ztype)*invrik
            tmpvar2 = tetprefactor(ztype)*dot*dot*df_ik*f_ij
            
c     force between atom i and atom k due to second term 
c     i.e ( f_ij*f_ik*dot*dot )
            
            fxk =  tmpvar*(rij_hat(1) - rik_hat(1)*ctheta) + 
     x             tmpvar2*rik_hat(1)
            fyk =  tmpvar*(rij_hat(2) - rik_hat(2)*ctheta) + 
     x             tmpvar2*rik_hat(2)
            fzk =  tmpvar*(rij_hat(3) - rik_hat(3)*ctheta) + 
     x             tmpvar2*rik_hat(3)        
            
c     Add in to forces, virial and stress tensor 
            
            katm = flist(m,ii)
            
            fxx(iatm) = fxx(iatm) + fxj + fxk
            fyy(iatm) = fyy(iatm) + fyj + fyk
            fzz(iatm) = fzz(iatm) + fzj + fzk
            
            fxx(jatm) = fxx(jatm) - fxj
            fyy(jatm) = fyy(jatm) - fyj
            fzz(jatm) = fzz(jatm) - fzj
            
            fxx(katm) = fxx(katm) - fxk
            fyy(katm) = fyy(katm) - fyk
            fzz(katm) = fzz(katm) - fzk         
            
            virord = virord - fxj*xj - fyj*yj - fzj*zj
            virord = virord - fxk*xk - fyk*yk - fzk*zk
            
            strs1 = strs1 + xk*fxk + xj*fxj
            strs2 = strs2 + xk*fyk + xj*fyj
            strs3 = strs3 + xk*fzk + xj*fzj
            strs5 = strs5 + yk*fyk + yj*fyj
            strs6 = strs6 + yk*fzk + yj*fzj
            strs9 = strs9 + zk*fzk + zj*fzj
            
          end do
        end do
      end do
      
c     Complete stress tensor
      
      stress(1)=stress(1)+strs1
      stress(2)=stress(2)+strs2
      stress(3)=stress(3)+strs3
      stress(4)=stress(4)+strs2
      stress(5)=stress(5)+strs5
      stress(6)=stress(6)+strs6
      stress(7)=stress(7)+strs3
      stress(8)=stress(8)+strs6
      stress(9)=stress(9)+strs9
      
      
c     Tidy up
      
      deallocate(xdf,stat=ierr(1))
      deallocate(ydf,stat=ierr(2))
      deallocate(zdf,stat=ierr(3))
      if (any(ierr/=0)) call Mfrz_Error(2547,0.d0)
      
      return
      
      end Subroutine Compute_Tetrahedral_Forces
      
      Subroutine Mfrz_Error(kode,arg)
      
c---------------------------------------------------------------------
c     
c     Author W. Smith Daresbury Laboratory January 2011
c     Adapted from D. Quigley - University of Warwick
c     Copyright D. Quigley October 2008
c     
c---------------------------------------------------------------------
      
      use setup_module,   only : nrite,nhist,nread,nconf,nstats,
     x                           nrest,nfield,ntable,nevnt
      
      implicit none
      integer,intent(in) :: kode
      real(8),intent(in) :: arg
      
      if(onroot)then
        
        if(kode.eq.2500)then
          
          write(nrite,"(
     x      'Error in number of collective variables - '//
     x      'ncolvar too small?'
     x      )")
          
        elseif(kode.eq.2501)then
          
          write(nrite,"(
     x      'Wang-Landau style recursion not yet implemented'//
     x      'for ncolvar > 1'             
     x      )")
          
        elseif(kode.eq.2502)then
          
          write(nrite,"('Unrecognised Gaussian height scheme')")
          
        elseif(kode.eq.2503)then
          
          write(nrite,"('Error maxhis exceeded in metadynamics')")
          
        elseif(kode.eq.2504)then
          
          write(nrite,"(
     x      'Error allocating comms buffer in compute_bias_potential'
     x      )")
          
        elseif(kode.eq.2505)then
          
          write(nrite,"('Error allocating driven array')")
          
        elseif(kode.eq.2506)then
          
          write(nrite,"('Could not open METACONTROL')")
          
        elseif(kode.eq.2508)then
          
          write(nrite,"('Comms error in metadynamics setup')")
          
        elseif(kode.eq.2509)then
          
          write(nrite,"(
     x      'Cannot bias local and global PE in same run'
     x      )")
          
        elseif(kode.eq.2510)then
          
          write(nrite,"('Error allocating local force arrays')")
          
        elseif(kode.eq.2511)then
          
          write(nrite,"(
     x      'Error allocating collective variables arrays'
     x      )")
          
        elseif(kode.eq.2512)then
          
          write(nrite,"('Error allocating Wang-Landau bins')")
          
        elseif(kode.eq.2515)then
          
          write(nrite,"(
     x      'Error allocating Steinhardt parameter arrays'       
     x      )")
          
        elseif(kode.eq.2516)then
          
          write(nrite,"('Could not open STEINHARDT')")
          
        elseif(kode.eq.2517)then
          
          write(nrite,"('Error allocating q4site')")
          
        elseif(kode.eq.2518)then
          
          write(nrite,"('Error allocating q6site')")
          
        elseif(kode.eq.2519)then
          
          write(nrite,"('Error deallocating buff')")
          
        elseif(kode.eq.2521)then
          
          write(nrite,"('Error reading line ',i5,' of STEINHARDT'
     x      )")nint(arg)
          
        elseif(kode.eq.2522)then
          
          write(nrite,"(
     x      'Error allocating Steinhardt parameter arrays'       
     x      )")
          
        elseif(kode.eq.2523)then
          
          write(nrite,"('Could not open ZETA')")
          
        elseif(kode.eq.2524)then
          
          write(nrite,"('Error allocating zetasite')")
          
        elseif(kode.eq.2525)then
          
          write(nrite,"('Error allocating full neighbour list')")
          
        elseif(kode.eq.2527)then
          
          write(nrite,"(
     x      'Number of collective variables incorrect  for specified'//
     x      'order parameters'
     x      )")
          
        elseif(kode.eq.2529)then
          
          write(nrite,"('Error reading line ',i5,' of ZETA'
     x      )")nint(arg)
          
        elseif(kode.eq.2531)then
          
          write(nrite,"('Comms error on reading METADYNAMICS')")
          
        elseif(kode.eq.2532)then
          
          write(nrite,"('Error in fc function - out of range')")
          write(nrite,"('Value of r was ',1p,e14.6)")arg
          
        elseif(kode.eq.2533)then
          
          write(nrite,"(
     x      'Error allocating solvation arrays for metadynamics'
     x      )")
          
        elseif(kode.eq.2534)then
          
          write(nrite,"('Error allocating comms buffer arrays')")
          
        elseif(kode.eq.2535)then
          
          write(nrite,"('Solvation list overrun')")
          
        elseif(kode.eq.2536)then
          
          write(nrite,"(
     x      'Error deallocating solvation arrays for metadynamics'
     x      )")
          
        elseif(kode.eq.2537)then
          
          write(nrite,"('Error deallocating comms buffer arrays')")
          
        elseif(kode.eq.2538)then
          
          write(nrite,"(
     x      'Error allocating solvation arrays for metadynamics'
     x      )")
          
        elseif(kode.eq.2540)then
          
          write(nrite,"('Error allocating force prefactor arrays')")
          
        elseif(kode.eq.2541)then
          
          write(nrite,"(
     x      'Memory allocation error in compute_tet_nlist'
     x      )")
          
        elseif(kode.eq.2542)then
          
          write(nrite,"(
     x      'Error in metafreeze_module.f90 mxninc too small'
     x      )")
          
        elseif(kode.eq.2543)then
          
          write(nrite,"('nnn too small in compute_tet_nlist')")
          
        elseif(kode.eq.2544)then
          
          write(nrite,"('mxflist too small in metafreeze_module')")
          
        elseif(kode.eq.2545)then
          
          write(nrite,"(
     x      'Memory deallocation error in compute_tet_nlist'
     x      )")
          
        elseif(kode.eq.2546)then
          
          write(nrite,"(
     x      'Memory allocation error in compute_tet_nlist'
     x      )")
          
        elseif(kode.eq.2547)then
          
          write(nrite,"(
     x      'Memory deallocation error in compute_tet_nlist'
     x      )")
          
        elseif(kode.eq.2548)then
          
          write(nrite,"(
     x      'Memory allocation error in compute_tet_nlist'
     x      )")
          
        endif
        
c     close all i/o channels
        
        close (nrite)
        close (nhist)
        close (nread)
        close (nconf)
        close (nstats)
        close (nrest)
        close (nfield)
        close (ntable)
        close (nevnt)
        close (stn)
        close (zta)
        close (mtd)
        
      endif
      
c     shut down communications
      
      call gsync()
      call exitcomms()
        
      end Subroutine Mfrz_Error
      
      end module metafreeze_module
      
