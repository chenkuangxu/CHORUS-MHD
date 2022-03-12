    subroutine compresid
    use setup2d
    implicit none
    include 'mpif.h'

    integer :: is,js,ic,rfp
    double precision :: dFpsi(8),dGeta(8),dFvpsi(9),dGveta(9),Qs(9),&
    advect_u,advect_v

    call compflux

    do ic=1,NCELL

        ! need to find the flux derivatives at each solution pt
   
        do js=1,N
        do is=1,N

            Qs(1:9) = Q(1:9,is,js,ic)
            advect_u = Qs(2)/Qs(1)
            advect_v = Qs(3)/Qs(1)
   
            ! compute the dFdpsi derivative
            dFpsi(1:8) = 0.d0
            dFvpsi(1:9) = 0.d0
            do rfp=1,N+1
                dFpsi(1:8) = dFpsi(1:8) + F1(1:8,rfp,js,ic)*Mmat(rfp,is)
                if (vismode==1) then
                dFvpsi(1:9) = dFvpsi(1:9) + Fv1(1:9,rfp,js,ic)*Mmat(rfp,is)
                end if
            end do
   
            ! compute the dGdeta derivative
            dGeta(1:8) = 0.d0
            dGveta(1:9) = 0.d0
            do rfp=1,N+1
                dGeta(1:8) = dGeta(1:8) + G2(1:8,is,rfp,ic)*Mmat(rfp,js)
                if (vismode==1) then
                dGveta(1:9) = dGveta(1:9) + Gv2(1:9,is,rfp,ic)*Mmat(rfp,js)
                end if
            end do
   
            resid(1:8,is,js,ic) = -( dFpsi(1:8) + dGeta(1:8) )

            resid(9,is,js,ic) = -advect_u * nablaAs(1,is,js,ic) - &
            advect_v * nablaAs(2,is,js,ic)
            resid(9,is,js,ic) = resid(9,is,js,ic) * Jac(is,js,ic)
   
            if (vismode==1) then
                resid(1:9,is,js,ic) = resid(1:9,is,js,ic) + &
                (dFvpsi(1:9)+dGveta(1:9))
            end if
   
        end do
        end do
   
    end do
    
    do ic=1,NCELL
        do js=1,N
        do is=1,N
            resid(1:9,is,js,ic)=resid(1:9,is,js,ic)/Jac(is,js,ic)
        end do
        end do
    end do

    end subroutine compresid
    