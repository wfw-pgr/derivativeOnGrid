program main
  use derivGridMod
  implicit none
  integer         , parameter   :: cLen      = 300
  integer         , parameter   :: lun       = 50
  integer         , parameter   :: x_=1, y_=2, z_=3
  character(cLen) , parameter   :: fieldFile = 'dat/field.dat'
  character(cLen) , parameter   :: diffeFile = 'dat/diffe.dat'
  integer                       :: nData, LI, LJ, LK, i, j, k, m
  character(cLen)               :: cmtMark
  double precision, allocatable :: xyz(:,:,:,:), phi(:,:,:,:), dphi(:,:,:,:)

  ! ------------------------------------------------------ !
  ! --- [1] preparation / Field Fetching               --- !
  ! ------------------------------------------------------ !
  open(lun,file=trim(fieldFile),status='old',form='formatted')
  read(lun,*)
  read(lun,*) cmtMark, nData, LI, LJ, LK
  allocate( xyz(3,LI,LJ,LK), phi(nData,LI,LJ,LK), dphi(3*3*nData,LI,LJ,LK) )
  write(6,*)
  write(6,*) '[main] Load Data from fieldFile       :: ', trim( fieldFile )
  write(6,*) "[main] Data Size :: nData, LI, LJ, LK :: ", nData, LI, LJ, LK
  write(6,*)
  do k=1, LK
     do j=1, LJ
        do i=1, LI
           read(lun,*) xyz(x_:z_,i,j,k), phi(:,i,j,k)
        enddo
     enddo
  enddo
  close(lun)
  do k=1, LK
     do j=1, LJ
        do i=1, LI
           do m=1, nData*3*3
              dphi(m,i,j,k) = 0.d0
           enddo
        enddo
     enddo
  enddo

  ! ------------------------------------------------------ !
  ! --- [2] Difference Calculation                     --- !
  ! ------------------------------------------------------ !
  call derivativeOnGrid3D( xyz, phi, dphi, LI, LJ, LK, nData )
  
  ! ------------------------------------------------------ !
  ! --- [3] save Difference results                    --- !
  ! ------------------------------------------------------ !
  open(lun,file=trim(diffeFile),status='replace',form='formatted')
  write(lun,*) '# xp yp zp phi dphidx dphidy dphidz'
  write(lun,*) '#', nData, LI, LJ, LK
  do k=1, LK
     do j=1, LJ
        do i=1, LI
           write(lun,*) xyz(x_:z_,i,j,k), phi(:,i,j,k), dphi(:,i,j,k)
        enddo
     enddo
  enddo
  close(lun)
  write(6,*) '[main] Save Data into diffeFile       :: ', trim( diffeFile )
  write(6,*)
  
end program main
