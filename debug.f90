    subroutine debug
    use setup2d
    implicit none

    integer :: ic,is,js,K,IFA,IB,aa,bb
    double precision :: Bx(N,N,3,3),By(N,N,3,3)

    call RK_S5P3

    call correct_B_via_curl_Az

    call comp_inviscid_flux

    call nablaQsp(1,9)
    
    if (rank.eq.0) then
        ic = 10
        print *,abs(nablaQs(6,1,1,1,ic)+nablaQs(7,2,1,1,ic))
        ! ic = 9
        ! IFA = 1
        ! K = IC2F(ic,IFA)
        ! IB = IF2IBCYCREM(K)
        ! print *,IB,NCYCREM
        ! print *,CYCREM2F_PROC(IB),CYCREM2PROC(IB)
        ! print *,Q(9,1,1,9),Q(9,1,1,10),Q(9,1,1,11)
        ! print *,Qvfj(9,1,1:2,9)
        ! print *,Q(6,1,1,9),Q(6,1,1,10),Q(6,1,1,11)
        ! do IFA = 1,4
        !     K = IC2F(ic,IFA)
        !     print *,K,IFA
        !     print *,'IF2C=',IF2C(K,1:2)
        ! end do
        ! print *,Q(9,1,1,10)
        ! print *,Q(9,1,1,6862)
        ! print *,Q(9,1,1,11)
        ! print *,Q(9,1,1,82)
        ! print *,Q(9,1,1,9)
        ! print *,(Q(9,1,1,82)-Q(9,1,1,6862)) / (2*pi*2/200)
        ! print *,Q(9,1,1,82),Q(9,1,1,10),Qvfj(9,1,2,10)
        ! print *,Q(9,1,1,6862),Q(9,1,1,10),Qvfj(9,1,1,10)
        ! print *,Q(6,1,1,10)
    end if

    if (rank.eq.3) then
        ! print *,IF2C(13358,1:2),NCELL
        ! print *,IF2C(13360,1:2)
        ! print *,IF2C(13362,1:2)
        ! print *,Q(9,1,1,6593)
        ! print *,Q(9,1,1,6594)
        ! print *,Q(9,1,1,6595)
    end if

    ! if (rank.eq.0) then
    !     ic = 9
    !     print *,Q(9,1,1,ic),Q(9,2,1,ic),Q(9,2,2,ic),Q(9,1,2,ic)
    !     ic = 10
    !     print *,Q(9,1,1,ic),Q(9,2,1,ic),Q(9,2,2,ic),Q(9,1,2,ic)
    !     ic = 11
    !     print *,Q(9,1,1,ic),Q(9,2,1,ic),Q(9,2,2,ic),Q(9,1,2,ic)
    !     ic = 81
    !     print *,Q(9,1,1,ic),Q(9,2,1,ic),Q(9,2,2,ic),Q(9,1,2,ic)
    !     ic = 82
    !     print *,Q(9,1,1,ic),Q(9,2,1,ic),Q(9,2,2,ic),Q(9,1,2,ic)
    !     ic = 83
    !     print *,Q(9,1,1,ic),Q(9,2,1,ic),Q(9,2,2,ic),Q(9,1,2,ic)
    ! end if

    ! if (rank.eq.3) then
    !     ic = 6593
    !     print *,Q(9,1,1,ic),Q(9,2,1,ic),Q(9,2,2,ic),Q(9,1,2,ic)
    !     ic = 6594
    !     print *,Q(9,1,1,ic),Q(9,2,1,ic),Q(9,2,2,ic),Q(9,1,2,ic)
    !     ic = 6595
    !     print *,Q(9,1,1,ic),Q(9,2,1,ic),Q(9,2,2,ic),Q(9,1,2,ic)
    ! end if

    ! if (rank.eq.0) then
    !     print *,Q(9,1,1,10),Q(9,1,2,10),Q(9,1,1,82),Q(9,1,2,82)
    ! end if

    if (rank.eq.0) then
        ic = 9
        aa = 1
        bb = 2
        do js = 1,N
        do is = 1,N
            Bx(is,js,aa,bb) = Q(6,is,js,ic)
            By(is,js,aa,bb) = Q(7,is,js,ic)
        end do
        end do
        !!!!!!!!!!!!
        ic = 10
        aa = 2
        bb = 2
        do js = 1,N
        do is = 1,N
            Bx(is,js,aa,bb) = Q(6,is,js,ic)
            By(is,js,aa,bb) = Q(7,is,js,ic)
        end do
        end do
        !!!!!!!!!!!!
        ic = 11
        aa = 3
        bb = 2
        do js = 1,N
        do is = 1,N
            Bx(is,js,aa,bb) = Q(6,is,js,ic)
            By(is,js,aa,bb) = Q(7,is,js,ic)
        end do
        end do
        !!!!!!!!!!!!
        ic = 81
        aa = 1
        bb = 3
        do js = 1,N
        do is = 1,N
            Bx(is,js,aa,bb) = Q(6,is,js,ic)
            By(is,js,aa,bb) = Q(7,is,js,ic)
        end do
        end do
        !!!!!!!!!!!!
        ic = 82
        aa = 2
        bb = 3
        do js = 1,N
        do is = 1,N
            Bx(is,js,aa,bb) = Q(6,is,js,ic)
            By(is,js,aa,bb) = Q(7,is,js,ic)
        end do
        end do
        !!!!!!!!!!!!
        ic = 83
        aa = 3
        bb = 3
        do js = 1,N
        do is = 1,N
            Bx(is,js,aa,bb) = Q(6,is,js,ic)
            By(is,js,aa,bb) = Q(7,is,js,ic)
        end do
        end do
    end if

    if (rank.eq.3) then
        !!!!!!!!!!!!
        ic = 6593
        aa = 1
        bb = 1
        do js = 1,N
        do is = 1,N
            Bx(is,js,aa,bb) = Q(6,is,js,ic)
            By(is,js,aa,bb) = Q(7,is,js,ic)
        end do
        end do
        !!!!!!!!!!!!
        ic = 6594
        aa = 2
        bb = 1
        do js = 1,N
        do is = 1,N
            Bx(is,js,aa,bb) = Q(6,is,js,ic)
            By(is,js,aa,bb) = Q(7,is,js,ic)
        end do
        end do
        !!!!!!!!!!!!
        ic = 6595
        aa = 3
        bb = 1
        do js = 1,N
        do is = 1,N
            Bx(is,js,aa,bb) = Q(6,is,js,ic)
            By(is,js,aa,bb) = Q(7,is,js,ic)
        end do
        end do
    end if

    if (rank.eq.0) then
        ! print *,Q(9,1,1,9),Q(9,2,1,9),Q(9,1,1,10),Q(9,2,1,10),Q(9,1,1,11),Q(9,2,1,11),'bb'
        ! print *,Q(9,1,2,9),Q(9,2,2,9),Q(9,1,2,10),Q(9,2,2,10),Q(9,1,2,11),Q(9,2,2,11),'bb'
        ! print *,Q(9,1,1,81),Q(9,2,1,81),Q(9,1,1,82),Q(9,2,1,82),Q(9,1,1,83),Q(9,2,1,83),'bb'
        ! print *,Q(9,1,2,81),Q(9,2,2,81),Q(9,1,2,82),Q(9,2,2,82),Q(9,1,2,83),Q(9,2,2,83),'bb'
    end if

    if (rank.eq.3) then
        ! print *,Q(9,1,1,6593),Q(9,2,1,6593),Q(9,1,1,6594),Q(9,2,1,6594),Q(9,1,1,6595),&
        ! Q(9,2,1,6595),'bb'
        ! print *,Q(9,1,2,6593),Q(9,2,2,6593),Q(9,1,2,6594),Q(9,2,2,6594),Q(9,1,2,6595),&
        ! Q(9,2,2,6595),'bb'
    end if

    if (rank.eq.0) then
        print *,Q(6,1,1,9),Q(6,2,1,9),Q(6,1,1,10),Q(6,2,1,10),Q(6,1,1,11),Q(6,2,1,11),'b'
    end if


    if (rank.eq.3) then
        print *,Q(9,1,1,6594),Q(9,1,2,6594),'rank==3'
    end if

    if (rank.eq.0) then
        print *,Q(9,1,1,10),Q(9,1,2,10),Q(9,1,1,82),Q(9,1,2,82),'rank==0'
    end if

    ! if (rank.eq.0) then
    !     print *,Qvfj(9,1,1,10),Qvfj(9,1,2,10),Qvfj(9,1,3,10),'Qvfj'
    ! end if

    if (rank.eq.0) then
        print *,Qvfi(6,3,1,9),'Qvfi'
    end if

    if (rank.eq.0) then
        print *,nablaQs(9,2,1,1,10),'nablaQs'
        print *,Q(6,1,1,10),'Bx'
        print *,nablaQs(6,1,1,1,10),nablaQs(7,2,1,1,10),'der'
    end if

    end subroutine debug

    subroutine check_divergence_B
    use setup2d
    implicit none
    include 'mpif.h'

    double precision :: divB,dbuf(NPROC),divBglob
    integer :: ic,is,js,i

    call compflux

    divB = 0.d0
    do ic=1,NCELL
    do js=1,N
    do is=1,N
        divB = divB + (nablaQs(6,1,is,js,ic)+nablaQs(7,2,is,js,ic))**2
    end do
    end do 
    end do

    ! Collect the residual norm from each processor
    CALL MPI_GATHER(divB, 1, MPI_DOUBLE_PRECISION, dbuf, 1, &
    MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERROR)

    if (rank.eq.0) then
        divBglob = 0.d0
        do i=1,NPROC
            divBglob = divBglob + dbuf(i)
        enddo
        divBglob = sqrt(divBglob/(NCELLGLOB*N*N))
        write (*,*) 'Rank:',rank,'Iteration:',  &
        iter,'divB: ',divBglob
    end if

    end subroutine check_divergence_B