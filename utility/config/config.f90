
!************************************************************************************
! DL_POLY Quantum additional program for generating multiple CONFIG files from a single
! HISTORY file. Used to generate intial configuration files for real-time dynamics 
! simulations after performing a sampling simulation
!
! Authors: Nathan London and Dil Limbu
!
!************************************************************************************
program config

  character(80) :: header
  character(20) :: filenm
  character(5)  :: filenum
  character(8)  :: timestep
  character(8),allocatable :: atmnam(:)
  integer       :: i,j,k
  integer       :: keytrj,imcon,natms,nstep,iatm,stepread,stepwrite
  real(8)       :: tstep,hold(9)
  real(8),allocatable :: xxx(:,:),yyy(:,:),zzz(:,:),vxx(:,:),vyy(:,:),vzz(:,:)
  real(8),allocatable :: fxx(:,:),fyy(:,:),fzz(:,:),cell(:,:)
  real(8),allocatable :: weight(:),charge(:)
 
  open(1,file="CONTROL",status='old')

  read(1,*) stepread
  read(1,*) stepwrite
  
  close(1)
  open(2,file='HISTORY',status='old')

  read(2,'(a)') header
  read(2,'(3I10)') keytrj,imcon,natms

  allocate(atmnam(natms))
  allocate(weight(natms))
  allocate(charge(natms))
  allocate(xxx(stepread,natms))
  allocate(yyy(stepread,natms))
  allocate(zzz(stepread,natms))
  allocate(cell(stepread,9))

  if(keytrj.gt.0)then
    allocate(vxx(stepread,natms))
    allocate(vyy(stepread,natms))
    allocate(vzz(stepread,natms))
  endif
  if(keytrj.gt.1)then
    allocate(fxx(stepread,natms))
    allocate(fyy(stepread,natms))
    allocate(fzz(stepread,natms))
  endif

  do i=1,stepread
    read(2,'(a8,4i10,f12.6)') timestep,nstep,natms,keytrj,imcon,tstep
    if(imcon.gt.0)then
      read(2,'(3g12.4)') cell(i,1:3) !hold(1:3)
      read(2,'(3g12.4)') cell(i,4:6) !hold(4:6)
      read(2,'(3g12.4)') cell(i,7:9) !hold(7:9)
    endif
    do j=1,natms
      read(2,'(a8,i10,2f12.6)') atmnam(j),iatm,weight(j),charge(j)
      read(2,'(3e12.4)') xxx(i,j),yyy(i,j),zzz(i,j)
      if(keytrj.gt.0) read(2,'(3e12.4)') vxx(i,j),vyy(i,j),vzz(i,j)
      if(keytrj.gt.1) read(2,'(3e12.4)') fxx(i,j),fyy(i,j),fzz(i,j)
    enddo
  enddo

  close(2)

  
  do k=1,stepwrite
    write(filenum,'(i0.3)') k
    filenm ='CONFIG' // trim(filenum)
    open(3,file=filenm,status='replace')  
    
    i=(k-1)*stepread/stepwrite+1
    write(3,'(a80)') filenm
    write(3,'(3i10)') keytrj,imcon,natms
    write(3,'(3f20.12)') cell(i,1:3)
    write(3,'(3f20.12)') cell(i,4:6)
    write(3,'(3f20.12)') cell(i,7:9)
    do j=1,natms
      write(3,'(a8,i10)') atmnam(j),j
      write(3,'(3g20.10)') xxx(i,j),yyy(i,j),zzz(i,j)
      if(keytrj.gt.0) write(3,'(3g20.12)') vxx(i,j),vyy(i,j),vzz(i,j)
      if(keytrj.gt.1) write(3,'(3g20.12)') fxx(i,j),fyy(i,j),fzz(i,j)
    enddo
    close(3)
  enddo

end program config
