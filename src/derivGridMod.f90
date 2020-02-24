module derivGridMod
contains

  ! ====================================================== !
  ! === derivative calculation on Cartesian Grid       === !
  ! ====================================================== !
  subroutine derivativeOnGrid3D( xyz, phi, dphi, LI, LJ, LK, nData )
    implicit none
    integer         , intent(in)  :: LI, LJ, LK, nData
    double precision, intent(in)  :: xyz(3,LI,LJ,LK), phi(nData,LI,LJ,LK)
    double precision, intent(out) :: dphi(3*3*nData,LI,LJ,LK)
    double precision              :: dx, dy, dz, dxInv, dyInv, dzInv, coef(9)
    integer                       :: i, j, k, m, of_
    integer         , parameter   :: x_ =1, y_ =2, z_ =3
    integer         , parameter   :: x1_=1, y1_=2, z1_=3
    integer         , parameter   :: x2_=4, y2_=5, z2_=6
    integer         , parameter   :: x3_=7, y3_=8, z3_=9

    ! ------------------------------------------------------ !
    ! --- [1] preparation1 Delta                         --- !
    ! ------------------------------------------------------ !
    dx    = xyz(x_,2,1,1) - xyz(x_,1,1,1)
    dy    = xyz(y_,1,2,1) - xyz(y_,1,1,1)
    dz    = xyz(z_,1,1,2) - xyz(z_,1,1,1)
    dxInv = 1.d0 / dx
    dyInv = 1.d0 / dy
    dzInv = 1.d0 / dz
    ! ------------------------------------------------------ !
    ! --- [2] preparation2 coef of derivative            --- !
    ! ------------------------------------------------------ !
    coef(x1_) = 0.5d0 * dxInv
    coef(y1_) = 0.5d0 * dyInv
    coef(z1_) = 0.5d0 * dzInv
    coef(x2_) =         dxInv**2
    coef(y2_) =         dyInv**2
    coef(z2_) =         dzInv**2
    coef(x3_) = 0.5d0 * dxInv**3
    coef(y3_) = 0.5d0 * dyInv**3
    coef(z3_) = 0.5d0 * dzInv**3
    do k=3, LK-2
       do j=3, LJ-2
          do i=3, LI-2
             do m=1, nData
                of_                 = 3*3*( m-1 )
                dphi(x1_+of_,i,j,k) = ( phi(m,i+1,j  ,k  )                           - phi(m,i-1,j  ,k  ) ) * coef(x1_)
                dphi(y1_+of_,i,j,k) = ( phi(m,i  ,j+1,k  )                           - phi(m,i  ,j-1,k  ) ) * coef(y1_)
                dphi(z1_+of_,i,j,k) = ( phi(m,i  ,j  ,k+1)                           - phi(m,i  ,j  ,k-1) ) * coef(z1_)
                dphi(x2_+of_,i,j,k) = ( phi(m,i+1,j  ,k  ) - 2.d0*phi(m,i  ,j  ,k  ) + phi(m,i-1,j  ,k  ) ) * coef(x2_)
                dphi(y2_+of_,i,j,k) = ( phi(m,i  ,j+1,k  ) - 2.d0*phi(m,i  ,j  ,k  ) + phi(m,i  ,j-1,k  ) ) * coef(y2_)
                dphi(z2_+of_,i,j,k) = ( phi(m,i  ,j  ,k+1) - 2.d0*phi(m,i  ,j  ,k  ) + phi(m,i  ,j  ,k-1) ) * coef(z2_)
                dphi(x3_+of_,i,j,k) = ( phi(m,i+2,j  ,k  ) - 2.d0*phi(m,i+1,j  ,k  ) &
                     &                                     + 2.d0*phi(m,i-1,j  ,k  ) - phi(m,i-2,j  ,k  ) ) * coef(x3_)
                dphi(y3_+of_,i,j,k) = ( phi(m,i  ,j+2,k  ) - 2.d0*phi(m,i  ,j+1,k  ) &
                     &                                     + 2.d0*phi(m,i  ,j-1,k  ) - phi(m,i  ,j-2,k  ) ) * coef(y3_)
                dphi(z3_+of_,i,j,k) = ( phi(m,i  ,j  ,k+2) - 2.d0*phi(m,i  ,j  ,k+1) &
                     &                                     + 2.d0*phi(m,i  ,j  ,k-1) - phi(m,i  ,j  ,k-2) ) * coef(z3_)
             enddo
          enddo
       enddo
    enddo
    
    return
  end subroutine derivativeOnGrid3D

end module derivGridMod
