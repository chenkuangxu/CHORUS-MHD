    SUBROUTINE SETINITIALCOND
    use setup2d
    implicit none
    include 'mpif.h'

    integer :: ic,is,js
    double precision :: Qvi(numv),x,y,u,v,w,rho,Bx,By,Bz,&
    NORM_U,NORM_B,Az,pr,xx(2),error_B(3)

    do ic=1,NCELL
        do js=1,N
        do is=1,N
            xx(1:2) = XXsolu(1:2,is,js,ic)
            x = xx(1)
            y = xx(2)
            rho = gam**2
            u = -sin(y)
            v = sin(x)
            w = 0.d0
            pr = gam
            Bx = -sin(y)
            By = sin(x*2)
            Bz = 0.d0
            Az = 0.5d0*cos(x*2) + cos(y)
            NORM_U = u**2+v**2+w**2
            NORM_B = Bx**2+By**2+Bz**2
            Qvi(1) = rho
            Qvi(2) = rho*u
            Qvi(3) = rho*v
            Qvi(4) = rho*w
            Qvi(5) = pr/(gam-1) + 0.5d0*rho*NORM_U + 0.5d0*NORM_B
            Qvi(6) = Bx
            Qvi(7) = By
            Qvi(8) = Bz
            Qvi(9) = Az
            Q(1:numv,is,js,ic) = Qvi(1:numv)
        end do
        end do
    end do

    ctime = 0.d0
    iter = 0

    call correct_B_via_curl_Az

    do ic=1,NCELL
        do js=1,N
        do is=1,N
            xx(1:2) = XXsolu(1:2,is,js,ic)
            x = xx(1)
            y = xx(2)
            rho = gam**2
            u = -sin(y)
            v = sin(x)
            w = 0.d0
            pr = gam
            Bx = Q(6,is,js,ic)
            By = Q(7,is,js,ic)
            Bz = Q(8,is,js,ic)
            error_B(1) = Bx + sin(y)
            error_B(2) = By - sin(x*2)
            ! print *,error_B(1:2)
            Az = 0.5d0*cos(x*2) + cos(y)
            NORM_U = u**2+v**2+w**2
            NORM_B = Bx**2+By**2+Bz**2
            Qvi(1) = rho
            Qvi(2) = rho*u
            Qvi(3) = rho*v
            Qvi(4) = rho*w
            Qvi(5) = pr/(gam-1) + 0.5d0*rho*NORM_U + 0.5d0*NORM_B
            Qvi(6) = Bx
            Qvi(7) = By
            Qvi(8) = Bz
            Qvi(9) = Az
            Q(1:numv,is,js,ic) = Qvi(1:numv)
        end do
        end do
    end do

    END SUBROUTINE SETINITIALCOND

    SUBROUTINE correct_B_via_curl_Az
    use setup2d
    implicit none

    call comp_inviscid_flux

    call nablaQsp(9,9)

    call update_B

    END SUBROUTINE correct_B_via_curl_Az