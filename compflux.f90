    subroutine compflux
    use setup2d
    implicit none
    include 'mpif.h'

    eigvmax = 0.d0

    call comp_inviscid_flux

    call nablaAsp

    if (vismode==1) then

        call nablaQsp(1,9)
        call nablaQ_interior_FP
        call interfaceVISflux
        call PROCINTVISFLUX
        call CYCREMVISFLUX
        call BCVISFLUX
        call CYCLOCVISFLUX

        if (ARTIV==1) then
            call comp_artificial_viscosity
            call compAVflux
        end if

        call comp_viscous_flux

    end if

    end subroutine compflux

    subroutine comp_inviscid_flux
    use setup2d
    implicit none

    call inviscid_flux_interior

    call interfaceflux

    call PROCINTFLUX

    call CYCREMFLUX

    call CYCLOCFLUX

    call BCFLUX

    end subroutine comp_inviscid_flux

    subroutine inviscid_flux_interior
    use setup2d
    implicit none

    integer :: ic,ifp,jfp,sp
    double precision,dimension(9) :: Qf,Qs
    double precision,dimension(8) :: F,G

    do ic=1,NCELL
        ! firstly the i-dir flux points (only interior)
        ! FOR INTERIOR QF is the one which will be used by
        ! Qvfi and Qvfj store them temporarily
        do jfp=1,N
        do ifp=2,N
            ! we need to find the Q at the flux points
            Qf(1:9) = 0.d0
            do sp=1,N
                Qs(1:9) = Q(1:9,sp,jfp,ic)
                Qf(1:9) = Qf(1:9) + Qs(1:9)*Lmat(ifp,sp)
            end do

            if (vismode==1) then
                Qvfi(1:9,ifp,jfp,ic) = Qf(1:9)
                Azfi(ifp,jfp,ic) = Qf(9)
            end if
   
            CALL getIfluxvectors(Qf(1:8),F(1:8),G(1:8))
            F1(1:8,ifp,jfp,ic) = F(1:8)*S1(1,1,ifp,jfp,ic) + G(1:8)*S1(1,2,ifp,jfp,ic) 
        end do
        end do

        ! secondly the j-dir flux points (only interior)
        do jfp=2,N
        do ifp=1,N
            ! we need to find the Q at the flux points
            Qf(1:9) = 0.d0
            do sp=1,N
                Qs(1:9) = Q(1:9,ifp,sp,ic)
                Qf(1:9) = Qf(1:9) + Qs(1:9)*Lmat(jfp,sp)
            end do
            if (vismode==1) then
                Qvfj(1:9,ifp,jfp,ic) = Qf(1:9)
                Azfj(ifp,jfp,ic) = Qf(9)
            end if

            call getIfluxvectors(Qf(1:8),F(1:8),G(1:8))
            G2(1:8,ifp,jfp,ic) = F(1:8)*S2(2,1,ifp,jfp,ic) + G(1:8)*S2(2,2,ifp,jfp,ic)

        end do
        end do

    end do ! end of loop over cells

    end subroutine inviscid_flux_interior


    Subroutine nablaQsp(v_sta,v_end)
    !*********************************************************************************************
    !get nablaQ at the solution points
    ! using the formula (2.15) in Sun et al. (2007) Communications in Comput. Physics
    !**********************************************************************************************
    use setup2d
    IMPLICIT NONE
    
    integer,intent(in) :: v_sta,v_end
    integer :: is,js,ic,rfp
    double precision,allocatable,dimension(:) :: dQSpsi,dQSeta

    allocate(dQSpsi(numv),dQSeta(numv))
    
    do ic=1,NCELL
        do js=1,N
        do is=1,N
            dQSpsi(v_sta:v_end) = 0.d0
            do rfp=1,N+1
                dQSpsi(v_sta:v_end) = dQSpsi(v_sta:v_end) + &
                Qvfi(v_sta:v_end,rfp,js,ic) * Mmat(rfp,is) * S1(1,1,rfp,js,ic)
                ! ksi_x
            end do
            dQSeta(v_sta:v_end) = 0.d0
            do rfp=1,N+1
                dQSeta(v_sta:v_end) = dQSeta(v_sta:v_end) + &
                Qvfj(v_sta:v_end,is,rfp,ic) * Mmat(rfp,js) * S2(2,1,is,rfp,ic)
                ! eta_x
            end do
            nablaQs(v_sta:v_end,1,is,js,ic) = &
            ( dQSpsi(v_sta:v_end) + dQSeta(v_sta:v_end) ) / Jac(is,js,ic)
    !**********************************************************************************************
            dQSpsi(v_sta:v_end) = 0.d0
            do rfp=1,N+1
                dQSpsi(v_sta:v_end)= dQSpsi(v_sta:v_end) + &
                Qvfi(v_sta:v_end,rfp,js,ic) * Mmat(rfp,is) * S1(1,2,rfp,js,ic)
            end do
      
            dQSeta(v_sta:v_end) = 0.d0
            do rfp=1,N+1
                dQSeta(v_sta:v_end) = dQSeta(v_sta:v_end) + &
                Qvfj(v_sta:v_end,is,rfp,ic) * Mmat(rfp,js) * S2(2,2,is,rfp,ic)
            end do
            nablaQs(v_sta:v_end,2,is,js,ic) = &
            ( dQSpsi(v_sta:v_end) + dQSeta(v_sta:v_end) ) / Jac(is,js,ic)
    !**********************************************************************************************
        end do
        end do
    end do
    
    end Subroutine nablaQsp

    subroutine nablaQ_interior_FP
    use setup2d
    implicit none

    integer :: ic,ifp,jfp,sp
    double precision :: Qf(9),Qf2(9)

    do ic=1,NCELL
        do jfp=1,N
            do ifp=2,N ! for interior flux points only
                Qf(1:9) = 0.d0
                Qf2(1:9) = 0.d0
                do sp=1,N
                Qf(1:9) = Qf(1:9) + nablaQs(1:9,1,sp,jfp,ic)*Lmat(ifp,sp)
                Qf2(1:9) = Qf2(1:9) + nablaQs(1:9,2,sp,jfp,ic)*Lmat(ifp,sp)
                end do
                nablaQvfi(1:9,1,ifp,jfp,ic) = Qf(1:9)
                nablaQvfi(1:9,2,ifp,jfp,ic) = Qf2(1:9)
            end do
        end do

        ! secondly the j-dir flux points (only interior)

        do jfp=2,N
            do ifp=1,N
                Qf(1:9) = 0.d0
                Qf2(1:9) = 0.d0
                do sp=1,N
                Qf(1:9) = Qf(1:9) + nablaQs(1:9,1,ifp,sp,ic)*Lmat(jfp,sp)
                Qf2(1:9) = Qf2(1:9) + nablaQs(1:9,2,ifp,sp,ic)*Lmat(jfp,sp)
                end do
                nablaQvfj(1:9,1,ifp,jfp,ic) = Qf(1:9)
                nablaQvfj(1:9,2,ifp,jfp,ic) = Qf2(1:9)
            end do
        end do
    end do		! end of loop over cells for nablaQf

    end subroutine nablaQ_interior_FP

    subroutine comp_viscous_flux
    use setup2d
    implicit none

    integer :: ic,ifp,jfp
    double precision :: Qs(9),nablaQ(9,2),Fv(9),Gv(9),Fv_art(9),Gv_art(9)

    do ic = 1,NCELL

        do jfp = 1,N
        do ifp = 1,N+1

            Qs(1:9) = Qvfi(1:9,ifp,jfp,ic)
            nablaQ(1:9,1) = nablaQvfi(1:9,1,ifp,jfp,ic)
            nablaQ(1:9,2) = nablaQvfi(1:9,2,ifp,jfp,ic)

            call getVfluxvectors(Qs,nablaQ,Fv,Gv)

            Fv1(1:9,ifp,jfp,ic) = Fv(1:9)*S1(1,1,ifp,jfp,ic) &
            + Gv(1:9)*S1(1,2,ifp,jfp,ic)

            if (ARTIV==1) then

                call getVfluxvectors_AV(Qs,nablaQ,Fv_art,Gv_art,1,ifp,jfp,ic)

                Fv1(1:9,ifp,jfp,ic) = Fv1(1:9,ifp,jfp,ic) + &
                Fv_art(1:9)*S1(1,1,ifp,jfp,ic) &
                + Gv_art(1:9)*S1(1,2,ifp,jfp,ic)

            end if     

        end do
        end do

        do jfp = 1,N+1
        do ifp = 1,N

            Qs(1:9) = Qvfj(1:9,ifp,jfp,ic)
            nablaQ(1:9,1) = nablaQvfj(1:9,1,ifp,jfp,ic)
            nablaQ(1:9,2) = nablaQvfj(1:9,2,ifp,jfp,ic)

            call getVfluxvectors(Qs,nablaQ,Fv,Gv)

            Gv2(1:9,ifp,jfp,ic) = Fv(1:9)*S2(2,1,ifp,jfp,ic) &
            + Gv(1:9)*S2(2,2,ifp,jfp,ic)

            if (ARTIV==1) then

                call getVfluxvectors_AV(Qs,nablaQ,Fv_art,Gv_art,2,ifp,jfp,ic)

                Gv2(1:9,ifp,jfp,ic) = Gv2(1:9,ifp,jfp,ic) + &
                Fv_art(1:9)*S2(2,1,ifp,jfp,ic) &
                + Gv_art(1:9)*S2(2,2,ifp,jfp,ic)

            end if
            
        end do
        end do
    end do

    end subroutine comp_viscous_flux