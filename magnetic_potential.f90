    subroutine nablaAsp
    use setup2d
    implicit none
    include 'mpif.h'

    ! calculate dAzdx and dAzdy for the constrained transport
    ! equation for magnetic potential

    integer :: is,js,ic,rfp
    double precision :: dApsi,dAeta
    
    do ic=1,NCELL
        do js=1,N
        do is=1,N
            dApsi = 0.d0
            do rfp=1,N+1
                dApsi = dApsi + Azfi(rfp,js,ic) * Mmat(rfp,is) * S1(1,1,rfp,js,ic)
                ! ksi_x
            end do
            dAeta = 0.d0
            do rfp=1,N+1
                dAeta = dAeta + Azfj(is,rfp,ic) * Mmat(rfp,js) * S2(2,1,is,rfp,ic)
                ! eta_x
            end do
            nablaAs(1,is,js,ic) = ( dApsi + dAeta ) / Jac(is,js,ic)
    !**********************************************************************************************
            dApsi = 0.d0
            do rfp=1,N+1
                dApsi = dApsi + Azfi(rfp,js,ic) * Mmat(rfp,is) * S1(1,2,rfp,js,ic)
            end do
      
            dAeta = 0.d0
            do rfp=1,N+1
                dAeta = dAeta + Azfj(is,rfp,ic) * Mmat(rfp,js) * S2(2,2,is,rfp,ic)
            end do
            nablaAs(2,is,js,ic) = ( dApsi + dAeta ) / Jac(is,js,ic)
    !**********************************************************************************************
        end do
        end do
    end do

    end subroutine nablaAsp

    subroutine update_B
    use setup2d
    implicit none
    include 'mpif.h'

    ! compute the updated magnetic field via the curl of Az

    integer :: ic,is,js
    double precision :: Bx,By,Bz

    do ic = 1,NCELL
        do js = 1,N
        do is = 1,N
            Bx = nablaQs(9,2,is,js,ic)
            By = -nablaQs(9,1,is,js,ic)
            Bz = Q(8,is,js,ic)
            Q(6,is,js,ic) = Bx
            Q(7,is,js,ic) = By
            Q(8,is,js,ic) = Bz
        end do
        end do
    end do
            
    end subroutine update_B
