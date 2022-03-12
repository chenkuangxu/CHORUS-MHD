    SUBROUTINE BCFLUX
    use setup2d
    implicit none

    integer :: ifn,iface,icleft,ifacelc,faml,ifpl,jfpl,nfp,ic,sp
    integer :: sign_l,sign_r
    double precision,dimension(9) :: Qfl,Qfr,Qs
    double precision,dimension(8) :: Fnl,Fnr
    double precision,dimension(2) :: normf
    double precision :: eigv,Az

    do ifn=1,NFREE

        iface = IBFREE(ifn)
        icleft  = IF2C(iface,1)
        ifacelc = IF2C(iface,3)
        faml = mod(ifacelc,2) + 1

        do nfp=1,N

            ! lets get the left solution vectors at flux point
            ifpl = iface2fp(nfp,ifacelc)
            jfpl = jface2fp(nfp,ifacelc)
            ic = icleft
            sign_l = 1

            if(faml==1) then
                Qfl(1:9) = 0.d0
                do sp=1,N
                    Qs(1:9)  = Q(1:9,sp,jfpl,ic)
                    Qfl(1:9) = Qfl(1:9) + Qs(1:9)*Lmat(ifpl,sp)
                end do
                if(ifpl==1) sign_l=-1
                normf(1:2) = S1(1,1:2,ifpl,jfpl,ic)
           else
                Qfl(1:9) = 0.d0
                do sp=1,N
                    Qs(1:9) = Q(1:9,ifpl,sp,ic)
                    Qfl(1:9) = Qfl(1:9) + Qs(1:9)*Lmat(jfpl,sp)
                end do
                if(jfpl==1) sign_l=-1
                normf(1:2) = S2(2,1:2,ifpl,jfpl,ic)
            end if
            
            ! Set the BC on the right sol vector
            Qfr(1:9) = Qfl(1:9)
   
            if (vismode==1) then
                if(faml==1) then
                    Qvfi(1:9,ifpl,jfpl,icleft)  = 0.5d0*Qfl(1:9)+0.5d0*Qfr(1:9)
                else
                    Qvfj(1:9,ifpl,jfpl,icleft)  = 0.5d0*Qfl(1:9)+0.5d0*Qfr(1:9)
                end if
            end if
   
            ! now we have the Q on either side of face
            CALL getrusanovflux(Qfl,Qfr,Fnl,Fnr,Az,normf,sign_l,sign_r,eigv)
   
            if(faml==1) then
                F1(1:8,ifpl,jfpl,icleft)  = Fnl(1:8)
                Azfi(ifpl,jfpl,icleft) = Az
            else
                G2(1:8,ifpl,jfpl,icleft)  = Fnl(1:8)
                Azfj(ifpl,jfpl,icleft) = Az
            end if
        end do ! loop over flux point on exit boundary face
   
    end do ! loop over exit boundary faces

    END SUBROUTINE BCFLUX


    SUBROUTINE BCVISFLUX
    use setup2d
    implicit none

    integer :: ifn,iface,icleft,ifacelc,faml,ifpl,jfpl,nfp,ic,sp
    double precision,dimension(9) :: Qfl,Qfl2,Qfr,Qfr2

    DO ifn=1,NFREE

        iface = IBFREE(ifn)
        icleft  = IF2C(iface,1)
        ifacelc = IF2C(iface,3)
        faml = mod(ifacelc,2) + 1

        do nfp=1,N

            ifpl = iface2fp(nfp,ifacelc)
            jfpl = jface2fp(nfp,ifacelc)
            ic = icleft
            Qfl(1:9) = 0.d0
            Qfl2(1:9) = 0.d0

            if(faml==1) then
                do sp=1,N
                    Qfl(1:9) = Qfl(1:9) + nablaQs(1:9,1,sp,jfpl,ic)*Lmat(ifpl,sp)
                    Qfl2(1:9) = Qfl2(1:9) + nablaQs(1:9,2,sp,jfpl,ic)*Lmat(ifpl,sp)
                end do
            else
                do sp=1,N
                    Qfl(1:9) = Qfl(1:9) + nablaQs(1:9,1,ifpl,sp,ic)*Lmat(jfpl,sp)
                    Qfl2(1:9) = Qfl2(1:9) + nablaQs(1:9,2,ifpl,sp,ic)*Lmat(jfpl,sp)
                end do
            end if
   
            Qfr(1:9) = Qfl(1:9)
            Qfr2(1:9) = Qfl2(1:9)
           
            if(faml==1) then
                nablaQvfi(1:9,1,ifpl,jfpl,icleft) = 0.5d0*Qfl(1:9)+0.5d0*Qfr(1:9)
                nablaQvfi(1:9,2,ifpl,jfpl,icleft) = 0.5d0*Qfl2(1:9)+0.5d0*Qfr2(1:9)
            else
                nablaQvfj(1:9,1,ifpl,jfpl,icleft) = 0.5d0*Qfl(1:9)+0.5d0*Qfr(1:9)
                nablaQvfj(1:9,2,ifpl,jfpl,icleft) = 0.5d0*Qfl2(1:9)+0.5d0*Qfr2(1:9)
            end if
        end do
   
    END DO  

    END SUBROUTINE BCVISFLUX