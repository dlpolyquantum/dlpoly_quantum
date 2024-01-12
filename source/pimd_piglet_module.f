      module pimd_piglet_module

c**********************************************************************
c     
c     dl_poly_quantum module to implement PIGLET thermostat for path 
c     integral molecular dynamics
c     
c     reference: ceriotti, manolopoulos,
c                phys. rev. lett. 109, 100604 (2012)  
c                uhl, marx, ceriotti,  
c                j. chem. phys. 145, 054101 (2016)     
c     
c     author   - f. uhl
c
c     adapted/implemented in dl_poly_quantum by
c
c              - Dil Limbu and Nathan London 2023
c
c**********************************************************************
      
      use setup_module,   only : pi,boltz,hbar,mspimd,nrite
      use pimd_module,    only : zmass,rzmass,pxx,pyy,pzz,nbeads,nsp1
      use error_module,   only : error
      use utility_module, only : puni

      implicit none

      real(8),allocatable,save :: a_mat(:,:,:),c_mat(:,:,:)
      real(8),allocatable,save :: gle_T(:,:,:),gle_S(:,:,:)
      real(8),allocatable,save :: smallsx(:,:),smallsy(:,:),smallsz(:,:)
      real(8),allocatable,save :: zmass2(:),rzmass2(:)

      real(8),allocatable      :: tempx1(:,:),tempy1(:,:),tempz1(:,:)
      real(8),allocatable      :: tempx2(:,:),tempy2(:,:),tempz2(:,:)


      public alloc_piglet_arrays,dealloc_piglet_arrays
      public piglet_init,piglet_thermo_step


      contains

      subroutine alloc_piglet_arrays(idnode,mxnode)

c**********************************************************************
c     
c     dl_poly_quantum routine to allocate arrays for PIGLET thermostat
c     in normal mode for path integral molecular dynamics
c     
c     copyright - 
c     authors   - Dil Limbu and Nathan London 2023
c     
c**********************************************************************

      implicit none
      
      logical                  :: safe
      integer, intent(in)      :: idnode,mxnode
      integer, dimension(1:8) :: fail

      safe=.true.

c     Allocate arrays

      fail(:)=0

      allocate(a_mat(1:nsp1,1:nsp1,1:nbeads),      stat=fail(1))
      allocate(c_mat(1:nsp1,1:nsp1,1:nbeads),      stat=fail(2))
      allocate(gle_T(1:nsp1,1:nsp1,1:nbeads),      stat=fail(3))
      allocate(gle_S(1:nsp1,1:nsp1,1:nbeads),      stat=fail(4))
      allocate(smallsx(1:nsp1,1:mspimd),           stat=fail(5))
      allocate(smallsy(1:nsp1,1:mspimd),           stat=fail(6))
      allocate(smallsz(1:nsp1,1:mspimd),           stat=fail(7))
      allocate(zmass2(1:mspimd),rzmass2(1:mspimd), stat=fail(8))

      if(any(fail.gt.0))safe=.false.
      if(mxnode.gt.1)call gstate(safe)
      if(.not.safe)call error(idnode,4001)

      end subroutine alloc_piglet_arrays


      subroutine dealloc_piglet_arrays(idnode,mxnode)

c**********************************************************************
c     
c     dl_poly_quantum routine to deallocate array for PIGLET thermostat
c     in normal mode for path integral molecular dynamics
c     
c     copyright - 
c     authors   - Dil Limbu and Nathan London 2023
c     
c**********************************************************************

      implicit none
      
      logical                 :: safe
      integer, intent(in)     :: idnode,mxnode
      integer, dimension(1:3) :: fail

      fail(:)=0
      safe=.true.

      deallocate(a_mat,c_mat,gle_T,gle_S,  stat=fail(1))
      deallocate(smallsx,smallsy,smallsz,  stat=fail(2))
      deallocate(zmass2,rzmass2,           stat=fail(3))

      if(any(fail.gt.0))safe=.false.
      if(mxnode.gt.1)call gstate(safe)
      if(.not.safe)call error(idnode,4002)

      end subroutine dealloc_piglet_arrays


      subroutine read_gle_matrix(idnode,mxnode,natms,temp)

c**********************************************************************
c     
c     dl_poly_quantum routine to read gle matrix to get A MATRIX and
c     C MATRIX to run PIGLET thermostat  in normal mode for path 
c     integral molecular dynamics
c     
c     copyright - 
c     authors   - Dil Limbu and Nathan London 2023
c
c**********************************************************************

      implicit none

      logical             :: safe,file_exist
      integer, intent(in) :: idnode,mxnode,natms
      real(8), intent(in) :: temp
      
      integer           :: i,j,k
      integer           :: iatm0,iatm1,ndim
      integer           :: iunit,obrac,cbrac,read_err
      real(8)           :: tTemp

      character(len=30) :: mat_fname='GLEMAT'
      character(len=40) :: temp_input
      character(len=30) :: read_unit

      safe=.true.

      iatm0=((idnode*natms)/mxnode)
      iatm1=(((idnode+1)*natms)/mxnode)

      ndim=(iatm1-iatm0)

c      read A-MATRIX & C-MATRIX

      inquire(file=trim(mat_fname),exist=file_exist)
      if(.not.file_exist)then
        safe=.false.
        if(mxnode.gt.1)call gstate(safe)
        if(.not.safe)call error(idnode,4003)
      endif

      iunit=1
      open(unit=iunit,file=trim(mat_fname),action="READ",status="OLD")
      read_err = 0

      if(idnode.eq.0)write(nrite,'(/,1x,a)')
     x        'Reading GLE MATRIX for PIGLET thermostat:'

      DO WHILE(read_err == 0)
        READ (iunit, '(A40)', iostat=read_err) temp_input

        IF (read_err /= 0) EXIT

c       !Parse comment section
        IF(INDEX(temp_input, "#") /= 0)THEN
          IF(INDEX(temp_input, "T=") /= 0)THEN

            obrac=INDEX(temp_input, "T=") + 2
            cbrac=INDEX(temp_input, "K") - 2
            read_unit=temp_input(obrac:cbrac)
            read(read_unit,'(F4.0)')tTemp
            if(tTemp.ne.temp)then
               safe=.false.
               if(mxnode.gt.1)call gstate(safe)
               if(.not.safe)call error(idnode,4004)
            endif

          ELSEIF(INDEX(temp_input, "A MATRIX") /= 0)THEN

            obrac=INDEX(temp_input, "(") + 1
            cbrac=INDEX(temp_input, ")") - 1
            read_unit=temp_input(obrac:cbrac)

            if(idnode.eq.0)write(nrite,*)'# A MATRIX Unit: ',read_unit

            DO k=1,nbeads
              READ(iunit,'(A40)') temp_input
              if(idnode.eq.0)write(nrite,*) temp_input
              DO i=1,nsp1
                READ(iunit,*,iostat=read_err)
     x                       (a_mat(i,j,k),j=1,nsp1)
                if(idnode.eq.0)write(nrite,'(9e13.4)')
     x                              (a_mat(i,j,k),j=1,nsp1)
                IF(read_err /= 0)THEN
                  if(idnode.eq.0)WRITE(nrite,FMT=*)"Invalid PIGLET 
     x                              A-matrix Nr, imode: ",i-1,k
                  EXIT
                END IF
              END DO
            END DO

c           !convert to internal system units
            IF(read_err == 0)THEN
              CALL a_mat_to_sysunit(a_mat,read_unit)
            END IF

          ELSEIF(INDEX(temp_input, "C MATRIX") /= 0)THEN

            obrac=INDEX(temp_input, "(") + 1
            cbrac=INDEX(temp_input, ")") - 1
            read_unit=temp_input(obrac:cbrac)

            if(idnode.eq.0)write(nrite,*)'# C MATRIX Unit: ',read_unit

            DO k=1,nbeads
              READ(iunit,'(A40)') temp_input
              if(idnode.eq.0)write(nrite,*) temp_input
              DO i=1,nsp1
                READ(iunit,*,iostat=read_err)
     x                       (c_mat(i,j,k), j=1,nsp1)
                if(idnode.eq.0)write(6,'(9e13.4)')
     x                              (c_mat(i,j,k),j=1,nsp1)
                IF(read_err /= 0)THEN
                  if(idnode.eq.0) WRITE(nrite,FMT=*)"Invalid PIGLET 
     x                              C-matrix Nr, imode :",i-1,k
                  EXIT
                END IF
              END DO
            END DO

c           !convert to internal system units
            IF(read_err == 0)THEN
              CALL c_mat_to_sysunit(c_mat,read_unit)
            END IF

          ENDIF
        ENDIF

      ENDDO
      if(idnode.eq.0)write(nrite,'(/,1x,a)')
     x        'End of GLE MATRIX for PIGLET thermostat'

      close(iunit)

      end subroutine read_gle_matrix


      subroutine piglet_init(idnode,mxnode,natms,tstep,temp,uuu)

c**********************************************************************
c     
c     dl_poly_quantum routine to initialize to run PIGLET thermostat 
c     in normal mode for path integral molecular dynamics
c     
c     adapted/implemented in dl_poly_quantum by
c
c              - Dil Limbu and Nathan London 2023
c
c**********************************************************************

      implicit none

      logical             :: safe,noskip
      integer, intent(in) :: idnode,mxnode,natms
      real(8), intent(in) :: tstep,temp,uuu(102)
      
      integer  :: i,j,k,l
      real(8)  :: tmpMat(nsp1,nsp1),tmpSS(nsp1,nsp1)
      real(8)  :: dt
      integer  :: iatm0,iatm1,ndim
      integer  :: fail(1:6)

      iatm0=((idnode*natms)/mxnode)
      iatm1=(((idnode+1)*natms)/mxnode)

      ndim=(iatm1-iatm0)

      dt = 0.5d0*tstep

c      allocate temporary arrays

      safe=.true.
      fail(:)=0
      allocate(tempx1(1:nsp1,1:ndim), stat=fail(1))
      allocate(tempy1(1:nsp1,1:ndim), stat=fail(2))
      allocate(tempz1(1:nsp1,1:ndim), stat=fail(3))
      allocate(tempx2(1:nsp1,1:ndim), stat=fail(4))
      allocate(tempy2(1:nsp1,1:ndim), stat=fail(5))
      allocate(tempz2(1:nsp1,1:ndim), stat=fail(6))
      if(any(fail.gt.0))safe=.false.
      if(mxnode.gt.1)call gstate(safe)
      if(.not.safe)call error(idnode,4005)

c      read A-MATRIX & C-MATRIX

      call read_gle_matrix(idnode,mxnode,natms,temp)

c     compute T and S matrices on every process

      do k=1,nbeads

c        determine the deterministic part of the propagator
c        !T = EXP(-A_mat*dt) = EXP(-A_mat*dt/2)
c        values for j and k = 15 are way to high, but to be sure.
c        (its only executed once anyway)

        call matrix_exp(-dt*a_mat(:,:,k),nsp1,15,15,gle_T(:,:,k))

c        !S*TRANSPOSE(S) = C-T*C*TRANSPOSE(T)
c         tempS=c_mat(:,:,i)-matmul(gle_T(:,:,i),matmul(c_mat(:,:,i),
c     x          transpose(gle_T(:,:,i))))
c        !T*C:

        call dgemm('N','N',nsp1,nsp1,nsp1,1.d0,gle_T(:,:,k),nsp1,
     x            c_mat(:,:,k),nsp1,0.d0,tmpMat,nsp1)

c        !T*C*TRANSPOSE(T):

        call dgemm('N','T',nsp1,nsp1,nsp1,1.d0,tmpMat,nsp1,
     x            gle_T(:,:,k),nsp1,0.d0,tmpSS,nsp1)

        !C - T*C*TRANSPOSE(T):

        tmpSS(:,:)=c_mat(:,:,k)-tmpSS(:,:)

c        determine the stochastic part of the propagator
c        get S by cholesky decomposition of tmpSS

        call cholesky(tmpSS,gle_S(:,:,k),nsp1)

      enddo

c      ! initialize extra degrees of freedom for Markovian Dynamics
c      ! as a cholesky decomposition of C-matrix multiplied by a random
c      ! number vector ( or from restart - NOT implemented yet)

      do k=1,nbeads

        call cholesky(c_mat(:,:,k),tmpSS,nsp1)

c       fill a vector with random numbers

        noskip=.true.

        do i=1,ndim

          do j=1,nsp1

            tempx2(j,i)=gssrnd(noskip,uuu)
            tempy2(j,i)=gssrnd(noskip,uuu)
            tempz2(j,i)=gssrnd(noskip,uuu)

          enddo

        enddo

        call dgemm('N','N',nsp1,ndim,nsp1,1.d0,tmpSS,nsp1,
     x            tempx2,nsp1,0.d0,tempx1,nsp1)

        call dgemm('N','N',nsp1,ndim,nsp1,1.d0,tmpSS,nsp1,
     x            tempy2,nsp1,0.d0,tempy1,nsp1)

        call dgemm('N','N',nsp1,ndim,nsp1,1.d0,tmpSS,nsp1,
     x            tempz2,nsp1,0.d0,tempz1,nsp1)

        do i=1,ndim

          l=(k-1)*ndim+i

          do j=1,nsp1

            smallsx(j,l)=tempx1(j,i)
            smallsy(j,l)=tempy1(j,i)
            smallsz(j,l)=tempz1(j,i)

          enddo

        enddo
        
      enddo

c     fill the array for the sqrt of the masses
  
      do i=1,nbeads*ndim

         zmass2(i) = sqrt(zmass(i))
         rzmass2(i) = sqrt(rzmass(i))

      enddo

      end subroutine piglet_init


      subroutine piglet_thermo_step(idnode,mxnode,natms,tstep,temp,uuu)

c**********************************************************************
c     
c     dl_poly_quantum routine to apply PIGLET thermostat to momentum
c     in normal mode for path integral molecular dynamics
c     
c     adapted/implemented in dl_poly_quantum by
c
c              - Dil Limbu and Nathan London 2023
c
c**********************************************************************

      implicit none

      logical             :: noskip,newjob
      integer, intent(in) :: idnode,mxnode,natms
      real(8), intent(in) :: temp,uuu(102)
      real(8)             :: tstep
      integer             :: i,j,k
      integer             :: iatm0,iatm1,ndim
      integer             :: init

      save newjob
      data newjob/.true./

      if(newjob)then
        call piglet_init(idnode,mxnode,natms,tstep,temp,uuu)
        newjob=.false.
      endif

      iatm0=((idnode*natms)/mxnode)
      iatm1=(((idnode+1)*natms)/mxnode)

      ndim=(iatm1-iatm0)

      do k=1,nbeads

c        copy mass scaled momenta to temp1 matrix
c        p/sqrt(m) and so velocity v*sqrt(m)

        do i=1,ndim

          tempx1(1,i)=pxx((i-1)*nbeads+k)*rzmass2((i-1)*nbeads+k)
          tempy1(1,i)=pyy((i-1)*nbeads+k)*rzmass2((i-1)*nbeads+k)
          tempz1(1,i)=pzz((i-1)*nbeads+k)*rzmass2((i-1)*nbeads+k)

        enddo

c        copy extra degrees of freedom to the temp1 matrix

        do i=1,ndim

          do j=2,nsp1

            tempx1(j,i)=smallsx(j,(k-1)*ndim+i)
            tempy1(j,i)=smallsy(j,(k-1)*ndim+i)
            tempz1(j,i)=smallsz(j,(k-1)*ndim+i)

          enddo

        enddo

c        fill temp2 with gaussian random number

        noskip=.true.

        do j=1,nsp1

          do i=1,ndim

          tempx2(j,i)=gssrnd(noskip,uuu)
          tempy2(j,i)=gssrnd(noskip,uuu)
          tempz2(j,i)=gssrnd(noskip,uuu)

          enddo

        enddo

        i=(k-1)*ndim+1
c        smalls(:,i)=1*S*temp2 + 0*smalls

        call dgemm('N','N',nsp1,ndim,nsp1,1.d0,gle_S(:,:,k),nsp1,tempx2,
     x             nsp1,0.d0,smallsx(:,i),nsp1)

        call dgemm('N','N',nsp1,ndim,nsp1,1.d0,gle_S(:,:,k),nsp1,tempy2,
     x             nsp1,0.d0,smallsy(:,i),nsp1)

        call dgemm('N','N',nsp1,ndim,nsp1,1.d0,gle_S(:,:,k),nsp1,tempz2,
     x             nsp1,0.d0,smallsz(:,i),nsp1)

c        now add the product of T-matrix * old smalls vector
c        samlls(:,i) = 1*T*temp1 + 1*smalls

        call dgemm('N','N',nsp1,ndim,nsp1,1.d0,gle_T(:,:,k),nsp1,tempx1,
     x             nsp1,1.d0,smallsx(:,i),nsp1)

        call dgemm('N','N',nsp1,ndim,nsp1,1.d0,gle_T(:,:,k),nsp1,tempy1,
     x             nsp1,1.d0,smallsy(:,i),nsp1)

        call dgemm('N','N',nsp1,ndim,nsp1,1.d0,gle_T(:,:,k),nsp1,tempz1,
     x             nsp1,1.d0,smallsz(:,i),nsp1)

      enddo

c      ! Copy the mass scales momenta to the outgoing velocities

      do i=1,ndim

        do k=1,nbeads

            pxx((i-1)*nbeads+k)=smallsx(1,(k-1)*ndim+i)*
     x                             zmass2((i-1)*nbeads+k)

            pyy((i-1)*nbeads+k)=smallsy(1,(k-1)*ndim+i)*
     x                             zmass2((i-1)*nbeads+k)

            pzz((i-1)*nbeads+k)=smallsz(1,(k-1)*ndim+i)*
     x                             zmass2((i-1)*nbeads+k)

        enddo

      enddo

      end subroutine piglet_thermo_step

      
      subroutine a_mat_to_sysunit(a_mat,myunit)
      
c**********************************************************************
c     
c     dl_poly_quantum routine to convert A-MATRIX (unit of time) into
c     internal system unit (picosecond^-1) for PIGLET thermostat to 
c     momentum in normal mode for path integral molecular dynamics
c     
c     copyright - 
c     authors   - Dil Limbu and Nathan London 2023
c
c**********************************************************************

      implicit none

      integer                  :: i,j,k
      real(8), intent(inout)   :: a_mat(nsp1,nsp1,nbeads)
      real(8)                  :: sys2fac
      character(len=30), intent(inout) :: myunit

      select case(TRIM(myunit))
c      fs^-1 to internal unit of ps^-1
        case("femtoseconds^-1")
          myunit="fs^-1"
          sys2fac=1000.d0

        case("picoseconds^-1")
          myunit="ps^-1"
          sys2fac=1.d0

c      s^-1 to internal unit of ps^-1
        case("seconds^-1")
          myunit="s^-1"
          sys2fac=1.0e-12

c      a.u.^-1 to internal unit of ps^-1
        case("atomic time units^-1")
          myunit="a.u.^-1"
          sys2fac=41341.37458d0

        case default
          myunit="ps^-1"
          sys2fac=1.d0

      end select

      do  k=1,nbeads
        do j = 1, nsp1
          do i = 1, nsp1
            a_mat(i,j,k)=a_mat(i,j,k)*sys2fac
          enddo
        enddo
      enddo

      end subroutine


      subroutine c_mat_to_sysunit(c_mat,myunit)

c**********************************************************************
c     
c     dl_poly_quantum routine to convert C-MATRIX (unit of energy) into
c     internal system unit (K) for PIGLET thermostat to 
c     momentum in normal mode for path integral molecular dynamics
c     
c     C-MATRIX is scaled by nbeads to simulate at the physical
c     temperature 
c
c     copyright - 
c     authors   - Dil Limbu and Nathan London 2023
c
c**********************************************************************

      implicit none

      integer                  :: i,j,k
      real(8), intent(inout)   :: c_mat(nsp1,nsp1,nbeads)
      real(8)                  :: sys2fac
      character(len=30), intent(inout) :: myunit

      select case(TRIM(myunit))
c      eV to internal unit
        case("eV")
          myunit="eV"
          sys2fac=11604.524d0*boltz/dble(nbeads)

        case("K")
c      K to internal unit
          myunit="K"
          sys2fac=boltz/dble(nbeads)

        case("atomic energy units")
c      a.u. to internal unit
          myunit="a.u."
          sys2fac=315774.d0*boltz/dble(nbeads)

        case default
c      assuming K to internal unit
          myunit="K"
          sys2fac=boltz/dble(nbeads)
      end select

      do  k=1,nbeads
        do j = 1, nsp1
          do i = 1, nsp1
            c_mat(i,j,k)=c_mat(i,j,k)*sys2fac
          enddo
        enddo
      enddo

      end subroutine


      subroutine matrix_exp(M,n,j,k,EM)

c**********************************************************************
c     
c     routine to compute the expoential of a (square) matrix using 
c     scale and square algorithm.
c     Parameters:
c     M - The square matrix to exponentiate
c     n - The size of the matrix
c     j - The number of terms in the exponential to evaluate
c     k - The exponent of the power of two used to scale the matrix
c     Returns:
c     EM - The expoential of the given matrix
c     
c**********************************************************************

      implicit none
      integer, intent(in)  :: n, j, k
      real(8), intent(in)  :: M(n,n)
      real(8), intent(out) :: EM(n,n)

      integer :: i,p
      real(8) :: tc(j+1),SM(n,n)

      tc(1)=1.d0
      do i=1,j
        tc(i+1)=tc(i)/dble(i)
      enddo

c    !scale
      SM=M*(1.d0/2.d0**k)
      EM=0.d0
      do i=1,n
        EM(i,i)=tc(j+1)
      enddo

c    !taylor exp of scaled matrix
      do p=j,1,-1
        EM=matmul(SM,EM);
        do i=1,n
          EM(i,i)=EM(i,i)+tc(p)
        enddo
      enddo

c    !square
      do p=1,k
        EM=matmul(EM,EM)
      enddo

      end subroutine matrix_exp


      subroutine cholesky(SST,S,n)

c**********************************************************************
c     
c     routine to compute the brute-force stabilized Cholesky 
c     decomposition of a square matrix.
c     Parameters:
c     SST - The matrix to determine the Cholesky decomposition of
c     n - The size of the matrix
c     Returns:
c     S - The computed Cholesky decomposition
c     
c**********************************************************************

      implicit none
      integer, intent(in)  :: n
      real(8), intent(in)  :: SST(n,n)
      real(8), intent(out) :: S(n,n)
    
      integer  :: i,j,k
      real(8)  :: D(n), L(n,n)
   
      D=0.d0
      L=0.d0
      S=0.d0

      do i=1,n
        L(i,i)=1.d0
        D(i)=SST(i,i)
        do j=1,i-1
          L(i,j)=SST(i,j)
          do k=1,j-1
            L(i,j)=L(i,j)-L(i,k)*L(j,k)*D(k)
          end do
          if (D(j).ne.0.) L(i,j)=L(i,j)/D(j)
        end do
        do k=1,i-1
          D(i)=D(i)-L(i,k)*L(i,k)*D(k)
        end do
      end do

      do i=1,n
        do j=1,i
          if (D(j)>0.d0) S(i,j)=S(i,j)+L(i,j)*sqrt(D(j))
        end do
      end do

      end subroutine cholesky


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
      
c**********************************************************************

      end module pimd_piglet_module
