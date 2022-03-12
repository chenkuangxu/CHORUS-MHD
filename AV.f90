    SUBROUTINE comp_artificial_viscosity
    use setup2d
    implicit none
    include 'mpif.h'

    double precision,allocatable :: Qloc_low(:,:,:),Qsm_high(:,:,:)
    integer :: is,js,isl,jsl,ic,ivar
    double precision :: xsi,eta,hhval,hhval_low,Se(numv),Se_max(numv),Se_min(numv),&
    Se_max_g(numv),Se_min_g(numv),k_range(numv),s0(numv),se_low(numv),se_high(numv),&
    h_msh
    double precision,parameter :: mu_ratio = 1.d0

    allocate(Qloc_low(numv,N-1,N-1),Qsm_high(numv,N,N))

    do ic = 1,NCELL

        Qloc_low = 0.d0

        do jsl = 1,N-1
        do isl = 1,N-1

            xsi = Xs_low(isl)
            eta = Xs_low(jsl)

            do js = 1,N
            do is = 1,N

                Qloc_low(1:numv,isl,jsl) = Qloc_low(1:numv,isl,jsl) + &
                Q(1:numv,is,js,ic) * hhval(is,xsi) * hhval(js,eta)
                
            end do
            end do

        end do
        end do

! ====================================================================

        Qsm_high = 0.d0

        do js = 1,N
        do is = 1,N

            xsi = Xs(is)
            eta = Xs(js)

            do jsl = 1,N-1
            do isl = 1,N-1

                Qsm_high(1:numv,is,js) = Qsm_high(1:numv,is,js) + &
                Qloc_low(1:numv,isl,jsl) * &
                hhval_low(isl,xsi) * hhval_low(jsl,eta)

            end do
            end do

        end do
        end do

        Se = 0.d0
        do js = 1,N
        do is = 1,N

            do ivar = 1,numv
                if (ivar.ne.4.and.ivar.ne.8) then
                Se(ivar) = Se(ivar) + &
                (Qsm_high(ivar,is,js)-Q(ivar,is,js,ic))**2 / &
                (Q(ivar,is,js,ic)**2)
                end if
            end do

        end do
        end do

        Se(4) = 1.d0
        Se(8) = 1.d0

        do ivar = 1,numv
            Se(ivar) = log10(Se(ivar))
        end do

        SM_IND(1:numv,ic) = Se(1:numv)

        if (ic.eq.1) then
            do ivar = 1,numv
                Se_min(ivar) = Se(ivar)
                Se_max(ivar) = Se(ivar)
            end do
        else
            do ivar = 1,numv
                Se_min(ivar) = min(Se(ivar),Se_min(ivar))
                Se_max(ivar) = max(Se(ivar),Se_max(ivar))
            end do
        end if
            
    end do

    deallocate(Qloc_low)

    CALL MPI_ALLREDUCE(Se_min(1),Se_min_g(1),1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,IERROR)
    CALL MPI_ALLREDUCE(Se_min(2),Se_min_g(2),1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,IERROR)
    CALL MPI_ALLREDUCE(Se_min(3),Se_min_g(3),1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,IERROR)
    CALL MPI_ALLREDUCE(Se_min(4),Se_min_g(4),1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,IERROR)
    CALL MPI_ALLREDUCE(Se_min(5),Se_min_g(5),1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,IERROR)
    CALL MPI_ALLREDUCE(Se_min(6),Se_min_g(6),1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,IERROR)
    CALL MPI_ALLREDUCE(Se_min(7),Se_min_g(7),1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,IERROR)
    CALL MPI_ALLREDUCE(Se_min(8),Se_min_g(8),1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,IERROR)
    CALL MPI_ALLREDUCE(Se_min(9),Se_min_g(9),1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,IERROR)

    CALL MPI_ALLREDUCE(Se_max(1),Se_max_g(1),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERROR)
    CALL MPI_ALLREDUCE(Se_max(2),Se_max_g(2),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERROR)
    CALL MPI_ALLREDUCE(Se_max(3),Se_max_g(3),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERROR)
    CALL MPI_ALLREDUCE(Se_max(4),Se_max_g(4),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERROR)
    CALL MPI_ALLREDUCE(Se_max(5),Se_max_g(5),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERROR)
    CALL MPI_ALLREDUCE(Se_max(6),Se_max_g(6),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERROR)
    CALL MPI_ALLREDUCE(Se_max(7),Se_max_g(7),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERROR)
    CALL MPI_ALLREDUCE(Se_max(8),Se_max_g(8),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERROR)
    CALL MPI_ALLREDUCE(Se_max(9),Se_max_g(9),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERROR)

    do ivar = 1,numv
        se_low(ivar) = -6.d0
        se_high(ivar) = -1.d0
    end do

    do ivar = 1,numv
        if (ivar.eq.4.or.ivar.eq.7.or.ivar.eq.8.or.ivar.eq.9) then
            se_low(ivar) = 0.d0
            se_high(ivar) = 0.d0
        end if
    end do

    do ivar = 1,numv
        if (se_low(ivar).gt.se_high(ivar)) then
            print *,'AV lower bound greater than higher bound! Wrong!'
            if (rank.eq.0) then
                print *,ivar,se_low(ivar),se_high(ivar)
            end if
            stop
        end if
    end do

    k_range(1:numv) = (se_high(1:numv)-se_low(1:numv)) * 0.5d0
    s0(1:numv) = (se_high(1:numv)+se_low(1:numv)) * 0.5d0

    CALL MPI_ALLREDUCE(eigvmax,eigvmax_g,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERROR)

    h_msh = dsqrt(DXCYCL*DYCYCL/(NCELLGLOB))
    epsilon0 = eigvmax_g * h_msh / (N-1) * mu_ratio

    ! if (rank.eq.0) then
    !     ivar = 1
    !     print *,Se_min_g(ivar),se_low(ivar),se_high(ivar),Se_max_g(ivar)
    !     ivar = 6
    !     print *,Se_min_g(ivar),se_low(ivar),se_high(ivar),Se_max_g(ivar)
    ! end if

    do ic = 1,NCELL
        do ivar = 1,numv
            if (ivar.eq.4.or.ivar.eq.7.or.ivar.eq.8.or.ivar.eq.9) then
                muAV(ivar,ic) = 0.d0
            else
                if (SM_IND(ivar,ic).lt.se_low(ivar)) then
                    muAV(ivar,ic) = 0.d0
                else if (SM_IND(ivar,ic).gt.se_high(ivar)) then
                    muAV(ivar,ic) = epsilon0
                else
                    muAV(ivar,ic) = epsilon0 / 2.d0 * &
                    ( 1.d0 + sin( pi*(SM_IND(ivar,ic)-s0(ivar))/2.d0/k_range(ivar) ) )
                end if
            end if
        end do
    end do

    END SUBROUTINE comp_artificial_viscosity

    function hhval_low(i,xval)
    use setup2d
    implicit none

    integer,intent(in) :: i
    double precision,intent(in) :: xval
    double precision :: hhval_low, hvaln, hvald
    integer :: s, np

    hvaln = 1.0d0
    hvald = 1.0d0
    np = N-1
    do s = 1, np
        if (s/= i)then
            hvaln = hvaln * (xval - Xs_low(s))
            hvald = hvald * (Xs_low(i) - Xs_low(s))
        end if
    end do
    
    hhval_low = hvaln/hvald

    end function hhval_low

    SUBROUTINE compAVflux
    use setup2d
    implicit none
    include 'mpif.h'

    integer :: ic,ffp

    do ic = 1,NCELL
        do ffp = 1,2
            muAVfi(1:numv,ffp,ic) = muAV(1:numv,ic)
        end do
        do ffp = 1,2
            muAVfj(1:numv,ffp,ic) = muAV(1:numv,ic)
        end do
    end do

    call interfaceAVflux

    call procintAVflux

    call cycremAVflux

    call cyclocAVflux

    END SUBROUTINE compAVflux


    SUBROUTINE interfaceAVflux
    use setup2d
    implicit none
    
    double precision :: muAV_left(numv),muAV_right(numv)
    integer :: ifinter,iface,icleft,icright,ifacelc,ifacerc,ic

    ! calculate the interface fluxes at all interior faces (non-boundary faces)
    do ifinter=1,NINTER
      
        iface 	= IBFINTER(ifinter)
      
        icleft  = IF2C(iface,1)
        icright = IF2C(iface,2)
      
        ifacelc = IF2C(iface,3)
        ifacerc = IF2C(iface,4)

        ic = icleft
        select case(ifacelc)
        case(1)
            muAV_left(1:numv) = muAVfj(1:numv,1,ic)
        case(2)
            muAV_left(1:numv) = muAVfi(1:numv,2,ic)
        case(3)
            muAV_left(1:numv) = muAVfj(1:numv,2,ic)
        case(4)
            muAV_left(1:numv) = muAVfi(1:numv,1,ic)
        end select

        ic = icright
        select case(ifacerc)
        case(1)
            muAV_right(1:numv) = muAVfj(1:numv,1,ic)
        case(2)
            muAV_right(1:numv) = muAVfi(1:numv,2,ic)
        case(3)
            muAV_right(1:numv) = muAVfj(1:numv,2,ic)
        case(4)
            muAV_right(1:numv) = muAVfi(1:numv,1,ic)
        end select

        ! averaging
        ic = icleft
        select case(ifacelc)
        case(1)
            muAVfj(1:numv,1,ic) = 0.5d0 * (muAV_left(1:numv)+muAV_right(1:numv))
        case(2)
            muAVfi(1:numv,2,ic) = 0.5d0 * (muAV_left(1:numv)+muAV_right(1:numv))
        case(3)
            muAVfj(1:numv,2,ic) = 0.5d0 * (muAV_left(1:numv)+muAV_right(1:numv))
        case(4)
            muAVfi(1:numv,1,ic) = 0.5d0 * (muAV_left(1:numv)+muAV_right(1:numv))
        end select

        ic = icright
        select case(ifacerc)
        case(1)
            muAVfj(1:numv,1,ic) = 0.5d0 * (muAV_left(1:numv)+muAV_right(1:numv))
        case(2)
            muAVfi(1:numv,2,ic) = 0.5d0 * (muAV_left(1:numv)+muAV_right(1:numv))
        case(3)
            muAVfj(1:numv,2,ic) = 0.5d0 * (muAV_left(1:numv)+muAV_right(1:numv))
        case(4)
            muAVfi(1:numv,1,ic) = 0.5d0 * (muAV_left(1:numv)+muAV_right(1:numv))
        end select

        ! if (icright==NCELL/2) then
        !     print *,icleft,icright
        !     print *,muAV_left
        !     print *,'--'
        !     print *,muAV_right
        !     print *,'---'
        ! end if
        
    end do

    END SUBROUTINE interfaceAVflux

    SUBROUTINE procintAVflux
    use setup2d
    IMPLICIT NONE
    include 'mpif.h'
    
    integer :: ifacelc,ic
    integer :: k,MLENGTH,IB
    INTEGER :: TAG,DEST,SOURCE
    INTEGER :: REQHANDLE(NPROCINT),ISTAT(MPI_STATUS_SIZE,NPROCINT)
      
    double precision,dimension(numv) :: muAV_left,muAV_right
    double precision, dimension(:,:),allocatable :: RBUFC, SBUFC
      
    MLENGTH = numv
      
    allocate(RBUFC(MLENGTH,NPROCINT))
    allocate(SBUFC(MLENGTH,NPROCINT))
      
    ! FIRST POST NON-BLOCKING RECEIVES
      
    do IB=1,NPROCINT
        K = IBFPROC(IB)
        TAG = K
        SOURCE=PROCINT2PROC(IB)
        CALL MPI_IRECV(RBUFC(1,IB),MLENGTH,MPI_DOUBLE_PRECISION,SOURCE,TAG,&
        MPI_COMM_WORLD,REQHANDLE(IB),IERROR)
    end do
      
    ! Then loop over faces and send muAV on the left side
      
    do IB=1,NPROCINT

        K= IBFPROC(IB)
        ic = IF2C(K,1)
        ifacelc = IF2C(K,3)
      
        TAG = PROCINT2F_PROC(IB)
        DEST = PROCINT2PROC(IB)
        
        select case(ifacelc)
        case(1)
            muAV_left(1:numv) = muAVfj(1:numv,1,ic)
        case(2)
            muAV_left(1:numv) = muAVfi(1:numv,2,ic)
        case(3)
            muAV_left(1:numv) = muAVfj(1:numv,2,ic)
        case(4)
            muAV_left(1:numv) = muAVfi(1:numv,1,ic)
        end select

        SBUFC(1:numv,IB) = muAV_left(1:numv)
        
        ! Now send the buffer
        CALL MPI_SEND(SBUFC(1,IB),MLENGTH,MPI_DOUBLE_PRECISION,DEST,TAG,&
        MPI_COMM_WORLD,IERROR)
      
    end do  ! do loop over processor interface faces
      
    ! ----------------------------------------------------------
    CALL MPI_WAITALL(NPROCINT,REQHANDLE,ISTAT,IERROR)
    ! ------------------------------------------------------------

    ! -----------------------------------------------------
    ! Now that we have muAV_left and muAV_right, we average them
    ! -----------------------------------------------------
      
    do IB=1,NPROCINT
        K= IBFPROC(IB)
        ic  = IF2C(K,1)
      
        ifacelc = IF2C(K,3)

        muAV_right(1:numv) = RBUFC(1:numv,IB)

        muAV_left(1:numv) = SBUFC(1:numv,IB)

        select case(ifacelc)
        case(1)
            muAVfj(1:numv,1,ic) = 0.5d0 * (muAV_left(1:numv)+muAV_right(1:numv))
        case(2)
            muAVfi(1:numv,2,ic) = 0.5d0 * (muAV_left(1:numv)+muAV_right(1:numv))
        case(3)
            muAVfj(1:numv,2,ic) = 0.5d0 * (muAV_left(1:numv)+muAV_right(1:numv))
        case(4)
            muAVfi(1:numv,1,ic) = 0.5d0 * (muAV_left(1:numv)+muAV_right(1:numv))
        end select
      
    end do  ! do loop over processor interface faces
        
    deallocate(RBUFC)
    deallocate(SBUFC)
      
    END SUBROUTINE procintAVflux

    SUBROUTINE cycremAVflux
    use setup2d
    IMPLICIT NONE
    include 'mpif.h'
        
    integer :: ifacelc,ic
    integer :: k,MLENGTH,IB
    INTEGER :: TAG,DEST,SOURCE
    INTEGER :: REQHANDLE(NCYCREM),ISTAT(MPI_STATUS_SIZE,NCYCREM)
          
    double precision,dimension(numv) :: muAV_left,muAV_right
    double precision, dimension(:,:),allocatable :: RBUFC, SBUFC
          
    MLENGTH = numv
          
    allocate(RBUFC(MLENGTH,NCYCREM))
    allocate(SBUFC(MLENGTH,NCYCREM))
          
    ! FIRST POST NON-BLOCKING RECEIVES
          
    do IB=1,NCYCREM
        K = IBFCYCREM(IB)
        TAG = K
        SOURCE=CYCREM2PROC(IB)
        CALL MPI_IRECV(RBUFC(1,IB),MLENGTH,MPI_DOUBLE_PRECISION,SOURCE,TAG,&
        MPI_COMM_WORLD,REQHANDLE(IB),IERROR)
    end do
          
    ! Then loop over faces and send muAV on the left side
          
    do IB=1,NCYCREM
    
        K= IBFCYCREM(IB)
        ic = IF2C(K,1)
        ifacelc = IF2C(K,3)
          
        TAG = CYCREM2F_PROC(IB)
        DEST = CYCREM2PROC(IB)
            
        select case(ifacelc)
        case(1)
            muAV_left(1:numv) = muAVfj(1:numv,1,ic)
        case(2)
            muAV_left(1:numv) = muAVfi(1:numv,2,ic)
        case(3)
            muAV_left(1:numv) = muAVfj(1:numv,2,ic)
        case(4)
            muAV_left(1:numv) = muAVfi(1:numv,1,ic)
        end select
    
        SBUFC(1:numv,IB) = muAV_left(1:numv)
            
        ! Now send the buffer
        CALL MPI_SEND(SBUFC(1,IB),MLENGTH,MPI_DOUBLE_PRECISION,DEST,TAG,&
        MPI_COMM_WORLD,IERROR)
          
    end do  ! do loop over processor interface faces
          
    ! ----------------------------------------------------------
    CALL MPI_WAITALL(NCYCREM,REQHANDLE,ISTAT,IERROR)
    ! ------------------------------------------------------------
    
    ! -----------------------------------------------------
    ! Now that we have muAV_left and muAV_right, we average them
    ! -----------------------------------------------------
          
    do IB=1,NCYCREM

        K= IBFCYCREM(IB)
        ic  = IF2C(K,1)
          
        ifacelc = IF2C(K,3)
    
        muAV_right(1:numv) = RBUFC(1:numv,IB)
    
        muAV_left(1:numv) = SBUFC(1:numv,IB)
    
        select case(ifacelc)
        case(1)
            muAVfj(1:numv,1,ic) = 0.5d0 * (muAV_left(1:numv)+muAV_right(1:numv))
        case(2)
            muAVfi(1:numv,2,ic) = 0.5d0 * (muAV_left(1:numv)+muAV_right(1:numv))
        case(3)
            muAVfj(1:numv,2,ic) = 0.5d0 * (muAV_left(1:numv)+muAV_right(1:numv))
        case(4)
            muAVfi(1:numv,1,ic) = 0.5d0 * (muAV_left(1:numv)+muAV_right(1:numv))
        end select
          
    end do  ! do loop over processor interface faces
            
    deallocate(RBUFC)
    deallocate(SBUFC)
          
    END SUBROUTINE cycremAVflux
    
    SUBROUTINE cyclocAVflux
    use setup2d
    implicit none

    integer :: ipair,K1,K2,icleft,icright,ifacelc,ifacerc,ic
    double precision,dimension(numv) :: muAV_left,muAV_right

    ! calculate the interface AV at all interior faces (non-boundary faces)
    do ipair=1,NPAIR
      
        K1 = CYC2FPAIR(ipair,1)
        K2 = CYC2FPAIR(ipair,2)
      
        icleft  = IF2C(K1,1)
        icright = IF2C(K2,1)
      
        ifacelc = IF2C(K1,3)
        ifacerc = IF2C(K2,3)
      
        ic = icleft
        select case(ifacelc)
        case(1)
            muAV_left(1:numv) = muAVfj(1:numv,1,ic)
        case(2)
            muAV_left(1:numv) = muAVfi(1:numv,2,ic)
        case(3)
            muAV_left(1:numv) = muAVfj(1:numv,2,ic)
        case(4)
            muAV_left(1:numv) = muAVfi(1:numv,1,ic)
        end select

        ic = icright
        select case(ifacerc)
        case(1)
            muAV_right(1:numv) = muAVfj(1:numv,1,ic)
        case(2)
            muAV_right(1:numv) = muAVfi(1:numv,2,ic)
        case(3)
            muAV_right(1:numv) = muAVfj(1:numv,2,ic)
        case(4)
            muAV_right(1:numv) = muAVfi(1:numv,1,ic)
        end select

        ! averaging
        ic = icleft
        select case(ifacelc)
        case(1)
            muAVfj(1:numv,1,ic) = 0.5d0 * (muAV_left(1:numv)+muAV_right(1:numv))
        case(2)
            muAVfi(1:numv,2,ic) = 0.5d0 * (muAV_left(1:numv)+muAV_right(1:numv))
        case(3)
            muAVfj(1:numv,2,ic) = 0.5d0 * (muAV_left(1:numv)+muAV_right(1:numv))
        case(4)
            muAVfi(1:numv,1,ic) = 0.5d0 * (muAV_left(1:numv)+muAV_right(1:numv))
        end select

        ic = icright
        select case(ifacerc)
        case(1)
            muAVfj(1:numv,1,ic) = 0.5d0 * (muAV_left(1:numv)+muAV_right(1:numv))
        case(2)
            muAVfi(1:numv,2,ic) = 0.5d0 * (muAV_left(1:numv)+muAV_right(1:numv))
        case(3)
            muAVfj(1:numv,2,ic) = 0.5d0 * (muAV_left(1:numv)+muAV_right(1:numv))
        case(4)
            muAVfi(1:numv,1,ic) = 0.5d0 * (muAV_left(1:numv)+muAV_right(1:numv))
        end select
      
    end do	! do loop over interior faces

    END SUBROUTINE cyclocAVflux