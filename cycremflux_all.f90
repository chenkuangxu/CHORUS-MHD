    SUBROUTINE CYCREMFLUX
    use setup2d
    IMPLICIT NONE
    include 'mpif.h'
  
    integer :: icright,faml
    integer :: ifacelc
    integer :: ifpl,jfpl,ifpr,jfpr,ic,sp,ifinter,sign_l,sign_r
    integer :: nfp,mfp,count
    integer :: k,MLENGTH,mfprm,nfprm,M,IB
    INTEGER :: TAG,DEST,SOURCE,IC2,ivar
    INTEGER :: REQHANDLE(NCYCREM),ISTAT(MPI_STATUS_SIZE,NCYCREM)
  
    double precision,dimension(9)  :: Qs,Qfl,Qfr
    double precision,dimension(8) :: Fnl,Fnr
    double precision, dimension(2) :: normf
    double precision                :: eigv,Az
    double precision, dimension(:,:),allocatable :: RBUFC, SBUFC
  
  
    MLENGTH = N*numv
  
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
  
    ! Then loop over faces and send Qfl
  
    do IB=1,NCYCREM
        K= IBFCYCREM(IB)
        ic = IF2C(K,1)
        ifacelc = IF2C(K,3)
        faml = mod(ifacelc,2) + 1
  
        TAG = CYCREM2F_PROC(IB)
        DEST = CYCREM2PROC(IB)
        M = 1
  
        do nfp=1,N
            ifpl = iface2fp(nfp,ifacelc)
            jfpl = jface2fp(nfp,ifacelc)
   
            Qfl_c(1:9,nfp,IB) = 0.d0
  
            if(faml==1) then
                do sp=1,N
                Qs(1:9) = Q(1:9,sp,jfpl,ic)
                Qfl_c(1:9,nfp,IB) = Qfl_c(1:9,nfp,IB) + Qs(1:9)*Lmat(ifpl,sp)
                end do
            else
                do sp=1,N
                Qs(1:9) = Q(1:9,ifpl,sp,ic)
                Qfl_c(1:9,nfp,IB) = Qfl_c(1:9,nfp,IB) + Qs(1:9)*Lmat(jfpl,sp)
                end do
            end if

            do ivar = 1,9
                SBUFC(M+ivar-1,IB) = Qfl_c(ivar,nfp,IB)
            end do
  
            M=M+9
        end do  ! do loop over points on interface
  
        ! Now send the buffer
        CALL MPI_SEND(SBUFC(1,IB),MLENGTH,MPI_DOUBLE_PRECISION,DEST,TAG,&
        MPI_COMM_WORLD,IERROR)
  
    end do  ! do loop over processor interface faces
  
    ! ----------------------------------------------------------
    CALL MPI_WAITALL(NCYCREM,REQHANDLE,ISTAT,IERROR)
    ! ------------------------------------------------------------
    
    ! Now receive the information and store in appropriate Qfr
  
    do IB=1,NCYCREM
        M = 1
        ! Loop over remote values of and nfp
        do nfprm=1,N
  
            nfp = N+1-nfprm
            ! Now we have all the left information
            ! Need to receive Qfr
            do ivar = 1,9
                Qfr_c(ivar,nfp,IB) = RBUFC(M+ivar-1,IB)
            end do
  
            M=M+9
        end do
    end do
    ! -----------------------------------------------------
    ! Now that we have Qfl and Qfr, we can call getrusanov
    ! -----------------------------------------------------
  
    do IB=1,NCYCREM
        K= IBFCYCREM(IB)
        ic  = IF2C(K,1)
  
        ifacelc = IF2C(K,3)
        faml = mod(ifacelc,2) + 1
  
        do nfp=1,N
  
            ifpl = iface2fp(nfp, ifacelc)
            jfpl = jface2fp(nfp, ifacelc)  
            sign_l =1
  
            if(faml==1) then
                normf(1:2) = S1(1,1:2,ifpl,jfpl,ic)
                if(ifpl==1) sign_l=-1
            else
                normf(1:2) = S2(2,1:2,ifpl,jfpl,ic)
                if(jfpl==1) sign_l=-1
            end if
  
            ! Now we have all the left information
  
            sign_r = 1 ! don't care the right side

            Qfl(1:9) = Qfl_c(1:9,nfp,IB)
            Qfr(1:9) = Qfr_c(1:9,nfp,IB)
  
            CALL getrusanovflux(Qfl,Qfr,Fnl,Fnr,Az,&
            normf,sign_l,sign_r,eigv)
  
            if(faml==1) then
                F1(1:8,ifpl,jfpl,ic)  = Fnl(1:8)
                Azfi(ifpl,jfpl,ic) = Az
            else
                G2(1:8,ifpl,jfpl,ic)  = Fnl(1:8)
                Azfj(ifpl,jfpl,ic) = Az
            end if
  
            IF (vismode==1) Then
                if(faml==1) then
                    Qvfi(1:9,ifpl,jfpl,ic) = 0.5d0*Qfl(1:9) +0.5d0*Qfr(1:9)
                else
                    Qvfj(1:9,ifpl,jfpl,ic) = 0.5d0*Qfl(1:9) +0.5d0*Qfr(1:9)
                end if
            END IF
        end do  ! do loop over points on interface
  
    end do  ! do loop over processor interface faces
    
    deallocate(RBUFC)
    deallocate(SBUFC)
  
    END SUBROUTINE CYCREMFLUX

    SUBROUTINE CYCREMVISFLUX

    use setup2d
          
    IMPLICIT NONE
    include 'mpif.h'
          
    integer :: ifacelc,faml
    integer :: ifpl,jfpl,ic,iface,sp
    integer :: nfp,mfp
    integer :: k,MLENGTH,mfprm,nfprm,M,ivar
    INTEGER :: TAG,DEST,SOURCE,IC2
    INTEGER :: REQHANDLE(NCYCREM),ISTAT(MPI_STATUS_SIZE,NCYCREM)
          
    double precision, dimension(9)  :: Qs,Qfl,Qfr,Qfl2,Qfr2
    double precision, dimension(:,:),allocatable :: RBUF2C, SBUF2C
    
          
    MLENGTH = N*18
          
    allocate(RBUF2C(MLENGTH,NCYCREM))
    allocate(SBUF2C(MLENGTH,NCYCREM))
          
    do k=1,NCYCREM
        iface = IBFCYCREM(k)
        TAG = iface
        SOURCE=CYCREM2PROC(k)
        CALL MPI_IRECV(RBUF2C(1,K),MLENGTH,MPI_DOUBLE_PRECISION,SOURCE,TAG,&
        MPI_COMM_WORLD,REQHANDLE(k),IERROR)
    end do
          
    ! Then loop over faces and send Qfl,Qfl2
          
    do k=1,NCYCREM
        iface = IBFCYCREM(k)
        ic    = IF2C(iface,1)
          
        TAG = CYCREM2F_PROC(k)
        DEST = CYCREM2PROC(k)
        M = 1
          
        ifacelc = IF2C(iface,3)  
          
        faml = mod(ifacelc,2) + 1   
          
        do nfp=1,N
          
            ifpl = iface2fp(nfp,ifacelc)  
            jfpl = jface2fp(nfp,ifacelc) 
          
            Qfl_c(1:9,nfp,K) = 0.d0
            Qfl_c2(1:9,nfp,K) = 0.d0
          
            if(faml==1) then
          
                do sp=1,N
                    Qfl_c(1:9,nfp,K) = Qfl_c(1:9,nfp,K) &
                    + nablaQs(1:9,1,sp,jfpl,ic)*Lmat(ifpl,sp)
                    Qfl_c2(1:9,nfp,K) = Qfl_c2(1:9,nfp,K) &
                    + nablaQs(1:9,2,sp,jfpl,ic)*Lmat(ifpl,sp)
                end do
          
            else
          
                do sp=1,N
                    Qfl_c(1:9,nfp,K) = Qfl_c(1:9,nfp,K) &
                    + nablaQs(1:9,1,ifpl,sp,ic)*Lmat(jfpl,sp)
                    Qfl_c2(1:9,nfp,K) = Qfl_c2(1:9,nfp,K) &
                    + nablaQs(1:9,2,ifpl,sp,ic)*Lmat(jfpl,sp)
                end do
          
            end if
    
            do ivar = 1,9
                SBUF2C(M+ivar-1,K) = Qfl_c(ivar,nfp,K)
            end do
          
            M = M+9
    
            do ivar = 1,9
                SBUF2C(M+ivar-1,K) = Qfl_c2(ivar,nfp,K)
            end do
    
            M = M+9
    
        end do  ! do loop over points on interface,nfp
          
        ! Now send the buffer
        CALL MPI_SEND(SBUF2C(1,K),MLENGTH,MPI_DOUBLE_PRECISION,DEST,TAG,&
        MPI_COMM_WORLD,IERROR)
          
    end do  ! do loop over processor interface faces
          
    ! ----------------------------------------------------------
    CALL MPI_WAITALL(NCYCREM,REQHANDLE,ISTAT,IERROR)
    ! ------------------------------------------------------------
          
    ! Now receive the information and store in appropriate Qfr
          
    do k=1,NCYCREM
         M = 1
        ! Loop over remote values of mfp and nfp
        do nfprm=1,N
          
            nfp = N+1-nfprm
            ! Now we have all the left information
            ! Need to receive Qfr
    
            do ivar = 1,9
                Qfr_c(ivar,nfp,K) = RBUF2C(M+ivar-1,K)
            end do
            M = M+9
    
            do ivar = 1,9
                Qfr_c2(ivar,nfp,K) = RBUF2C(M+ivar-1,K)
            end do
            M = M+9
          
        end do
    end do
          
    ! -----------------------------------------------------
    ! Now that we have Qfl,Qfl2 and Qfr,Qfr2 we can calculate nablaQvf
    ! -----------------------------------------------------
          
    do k=1,NCYCREM
        iface   = IBFCYCREM(k)
        ic  = IF2C(iface,1)
    
        ifacelc = IF2C(iface,3)
            
        faml = mod(ifacelc,2) + 1
          
        do nfp=1,N
          
            ifpl = iface2fp(nfp,ifacelc)
            jfpl = jface2fp(nfp,ifacelc) 
          
            Qfl(1:9)  = Qfl_c(1:9,nfp,K)
            Qfl2(1:9) = Qfl_c2(1:9,nfp,K)
          
            Qfr(1:9)  = Qfr_c(1:9,nfp,K)
            Qfr2(1:9) = Qfr_c2(1:9,nfp,K)
    
            if(faml==1) then
          
                nablaQvfi(1:9,1,ifpl,jfpl,ic) = &
                0.5*Qfl(1:9)+0.5*Qfr(1:9)
                nablaQvfi(1:9,2,ifpl,jfpl,ic) = &
                0.5*Qfl2(1:9)+0.5*Qfr2(1:9)
          
            else
          
                nablaQvfj(1:9,1,ifpl,jfpl,ic) = &
                0.5*Qfl(1:9)+0.5*Qfr(1:9)
                nablaQvfj(1:9,2,ifpl,jfpl,ic) = &
                0.5*Qfl2(1:9)+0.5*Qfr2(1:9)
          
            end if
          
        end do  ! do loop over points on interface, nfp
          
    end do  ! do loop over processor interface faces     
        
    deallocate(RBUF2C)
    deallocate(SBUF2C)
      
    END SUBROUTINE CYCREMVISFLUX