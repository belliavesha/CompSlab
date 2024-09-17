
subroutine cbf(res,arg,cr)
      external cr
      real(4), intent(out)::res
      real(4), intent(in)::arg

      res=cr(arg)
      write(*,*) 'from fortran', arg,res
end subroutine cbf


subroutine fill(RX,sigma,NE,NZ,mw,xw,x,m,Theta,Compton_redistribution_aa)
!f2py intent(callback) Compton_redistribution_aa
      interface
            subroutine Compton_redistribution_aa(CR,x1,x2,m1,m2)
                  real(4), intent(in) :: x1,x2,m1,m2
                  real(4),dimension(2,2), intent(out) :: CR
            end subroutine
      end interface
!f2py external Compton_redistribution_aa
!f2py real(4) :: x1,x2,m1,m2
!f2py real(4) :: r(2,2)
!f2py real(4) :: call Compton_redistribution_aa(r,x1,x2,m1,m2)
      
      integer, intent(in):: NE,NZ
      real(4), dimension(NE,NE,NZ,NZ,2,2), intent(out) :: RX
      real(4), dimension(2,2) :: r,rm,rt,rmt,rf,rmf,rtf,rmtf 
      real(4), dimension(NZ), intent(in) :: mw,m
      real(4), dimension(NE), intent(in) :: x
      real(4), dimension(NE), intent(out) :: sigma
      real(4), intent(in) :: xw, Theta
      integer :: e,e1,d,d1,md,md1
      real(4) :: w,mu
      write(*,*) 'imin'
      do e=1,NE; do e1=e,NE;
            do d=1,NZ/2; do d1=d,NZ/2
                  ! write(*,*) 'well, the indeces are as following : ',e,e1,d,d1
                  md=NZ-d-1 ! # -mu
                  md1=NZ-d1-1 !# -mu1
                  w=mw(d1)*xw*mw(d)
                  !t= d1 > d
                  !f= e1 > e

                  if (.TRUE.) then 
                        write(*,*) 'well, the truth is true'
                        call Compton_redistribution_aa(r,x(e),x(e1),m(d),m(d1))
                        call Compton_redistribution_aa(rm,x(e),x(e1),m(d),m(md1))

                        ! write(*,*) 'well, you do not see me'
                        sigma(e1)=sigma(e1)+(r(1,1)+rm(1,1))*w
                        RX(e,e1,d,d1,:,:)=r
                        RX(e,e1,md,md1,:,:)=r
                        RX(e,e1,d,md1,:,:)=rm
                        RX(e,e1,md,d1,:,:)=rm
                  endif
                  if (d1 > d) then 
                        rt=transpose(r) 
                        rmt=transpose(rm) 
                        sigma(e1)=sigma(e1)+(rt(1,1)+rmt(1,1))*w
                        RX(e,e1,d1,d,:,:)=rt
                        RX(e,e1,md1,md,:,:)=rt
                        RX(e,e1,d1,md,:,:)=rmt
                        RX(e,e1,md1,d,:,:)=rmt
                  endif
                        
                  if (e1 > e) then

                        mu=exp((x(e)-x(e1))/Theta) *x(e1)**3/x(e)**3
                        rf=r*mu 
                        rmf=rm*mu  
                        sigma(e)=sigma(e)+(rf(1,1)+rmf(1,1))*w
                        RX(e1,e,d,d1,:,:)=rf
                        RX(e1,e,md,md1,:,:)=rf
                        RX(e1,e,d,md1,:,:)=rmf
                        RX(e1,e,md,d1,:,:)=rmf
                  endif
                  if (d1 > d .and. e1 > e) then 
                        rtf=rt*mu
                        rmtf=rmt*mu
                        sigma(e)=sigma(e)+(rtf(1,1)+rmtf(1,1))*w
                        RX(e1,e,d1,d,:,:)=rtf
                        RX(e1,e,md1,md,:,:)=rtf
                        RX(e1,e,d1,md,:,:)=rmtf
                        RX(e1,e,md1,d,:,:)=rmtf
                  endif
            ! write(*,*) 'well, itiration just ended, it is gonna continue up to',NE,NE,NZ/2,NZ/2
            enddo;enddo;
      enddo;enddo
      write(*,*) 'end'
end subroutine fill




! subroutine hello(sigma,x,ne,f)
!       integer, intent(in) :: ne
!       real(4), dimension(NE), intent(in) :: x
!       real(4), dimension(NE), intent(out) :: sigma
!       external f
!       write(*,*) f(3d0)
!       write(*,*) 'hi'
!       sigma=x
!       write(*,*) 'by'
! end subroutine hello
