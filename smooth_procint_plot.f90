    subroutine smooth_procint_plot
    use setup2d
    implicit none
    include 'mpif.h'

    integer,allocatable :: iface2tecp(:,:),jface2tecp(:,:)
    integer :: iface,Np,pt,ptr,K,ii,jj,gnid,ivar
    integer :: TAG,DEST,M,SOURCE,ifacelc,faml,ic,MLENGTH
    INTEGER :: REQHANDLE(NPROCINT),ISTAT(MPI_STATUS_SIZE,NPROCINT)
    double precision,allocatable :: SBUFP(:,:),RBUFP(:,:)
    double precision :: factor_right,xx(2),dist

    Np = NRE
    MLENGTH = (Np+1)*12

    allocate(SBUFP(MLENGTH,NPROCINT),RBUFP(MLENGTH,NPROCINT))

    allocate(iface2tecp(Np+1,4),jface2tecp(Np+1,4))

    ! for face 1
    iface = 1
    do pt=1,Np+1
        iface2tecp(pt,iface) = pt
        jface2tecp(pt,iface) = 1
    end do
  
    ! for face 2
    iface = 2
    do pt=1,Np+1
        iface2tecp(pt,iface) = Np+1
        jface2tecp(pt,iface) = pt
    end do
  
    ! for face 3
    iface = 3
    do pt=1,Np+1
        iface2tecp(pt,iface) = Np - pt + 2
        jface2tecp(pt,iface) = Np+1
    end do
  
    ! for face 4
    iface = 4
    do pt=1,Np+1
        iface2tecp(pt,iface) = 1
        jface2tecp(pt,iface) = Np - pt + 2
    end do

    ! FIRST POST NON-BLOCKING RECEIVES
    do K = 1,NPROCINT
        iface = IBFPROC(K)
        TAG = iface
        SOURCE = PROCINT2PROC(K)
        CALL MPI_IRECV(RBUFP(1,K),MLENGTH,MPI_DOUBLE_PRECISION,SOURCE,TAG,&
        MPI_COMM_WORLD,REQHANDLE(K),IERROR)
    end do

    do K = 1,NPROCINT
        iface = IBFPROC(K)
        ic = IF2C(iface,1)

        TAG = PROCINT2F_PROC(K)
        DEST = PROCINT2PROC(K)
        M = 1

        ifacelc = IF2C(iface,3)
        faml = mod(ifacelc,2) + 1

        do pt = 1,Np+1

            ii = iface2tecp(pt,ifacelc)
            jj = jface2tecp(pt,ifacelc)

            gnid = gnumnode(ii,jj,ic)

            SBUFP(M,K)   = dble(gnode_factor(gnid))
            do ivar = 1,9
                SBUFP(M+ivar,K) = plotQ(ivar,gnid)
            end do
            SBUFP(M+10,K) = plotX(1,gnid)
            SBUFP(M+11,K) = plotX(2,gnid)

            M = M+12

        end do ! do loop over points on processor interface

        call MPI_SEND(SBUFP(1,K),MLENGTH,MPI_DOUBLE_PRECISION,DEST,TAG,&
        MPI_COMM_WORLD,IERROR)

    end do ! do loop over processor interface faces

    call MPI_WAITALL(NPROCINT,REQHANDLE,ISTAT,IERROR)

    do K = 1,NPROCINT
        iface = IBFPROC(K)
        ic = IF2C(iface,1)

        ifacelc = IF2C(iface,3)
        faml = mod(ifacelc,2) + 1

        M = 1

        do ptr = 1,Np+1
            pt = Np+2-ptr
            ii = iface2tecp(pt,ifacelc)
            jj = jface2tecp(pt,ifacelc)
            gnid = gnumnode(ii,jj,ic)

            factor_right = RBUFP(M,K)
            
            gnode_factor(gnid) = gnode_factor(gnid)+nint(factor_right)

            do ivar = 1,9
            plotQ(ivar,gnid) = plotQ(ivar,gnid) + RBUFP(M+ivar,K)
            end do

            xx(1) = RBUFP(M+10,K)
            xx(2) = RBUFP(M+11,K)
            
            dist = sqrt( (plotX(1,gnid)-xx(1))**2 + (plotX(2,gnid)-xx(2))**2 )

            if (dist.gt.1d-5) then
                print *,'something wrong in tec smoothing on the processor interfaces!'
            end if
            
            M = M+12
            
        end do
        
    end do

    deallocate(SBUFP,RBUFP)

    deallocate(iface2tecp,jface2tecp)

    end subroutine smooth_procint_plot