!-----------------------------------------------------------------------------!
! MAE 5093 DOCUMENTATION | Engineering Numerical Analysis
!-----------------------------------------------------------------------------!
! >>> Adaptive Runge Kutta Method for system of equations (Lorenz system)
!     Bogacki-Shampine Method(ode23) for solving du/dt=rhs
!-----------------------------------------------------------------------------!
! References: 
! * Fundamentals of Engineering Numerical Analysis by P. Moin (2012) -- Ch.1
! * Numerical Recipes: The Art of Scientific Computing, 2nd Edition (1992) 
!-----------------------------------------------------------------------------!
! Written by Omer San
!            CFDLab, Oklahoma State University, cfdlab.osu@gmail.com
!            www.cfdlab.org
! 
! Last updated: Oct. 1, 2015
!-----------------------------------------------------------------------------!

program adaptiveRK
implicit none
integer,parameter::ne=3
real*8 ::a,b,hmin,hmax,sf,tol
real*8 ::u(ne)

call RKTABLE


!minimum and maximum available time steps
hmin= 1.0d-5
hmax= 1.0d-1

!tolerance and safety factor
tol = 1.0d-4
sf  = 0.5d0

!time interval
a=0.0d0
b=100.0d0

!initial guess
u(1) = 1.0d0 
u(2) = 1.0d0 
u(3) = 1.0d0 

open(20,file='num.plt')
write(20,*) 'variables ="t","u1","u2","u3"'
write(20,*) a,u(1),u(2),u(3)

open(30,file='dt.plt')
write(30,*) 'variables ="t","dt"'


!Bogacki-Shampine Method(ode23) for solving du/dt=f(u)
call BS23(a,b,hmin,hmax,sf,tol,ne,u)  !ode23


end


!-----------------------------------------------------------!
!Adaptive time stepping
!Bogacki-Shampine Method(ode23) for solving du/dt=rhs
!-----------------------------------------------------------!
SUBROUTINE BS23(a,b,hmin,hmax,sf,tol,ne,u)
implicit none
integer::ne
real*8 ::a,b,hmin,hmax,sf,tol
real*8 ::u(ne)
integer::i,j,k
real*8 ::abs23(4,4),cbs23(4),bbs2(4),bbs3(4)
real*8 ::r1(ne),r2(ne),r3(ne),r4(ne),uu(ne)
real*8 ::dt,t,tt,est,ratio,err(ne)

common /bs23table/abs23,cbs23,bbs2,bbs3

!initial step size (start optimistically)
dt=hmax

!t holds the current time 
t=a

!u holds the current u(t)
!do k=1,ne 
!u(k)=u(k)
!end do

do k=1,ne 

end do
    
j=0
i=0 !step counter
do while (t.le.b-hmin/2.0d0)

	!compute r1,r2,r3,r4
	call rhs(ne,t,u,r1)
    
		do k=1,ne 
		uu(k)=u(k)+dt*abs23(2,1)*r1(k)
		end do
		tt=t+cbs23(2)*dt
        
	call rhs(ne,tt,uu,r2)
    
    	do k=1,ne 
		uu(k)=u(k)+dt*(abs23(3,1)*r1(k)+abs23(3,2)*r2(k))
        end do
		tt=t+cbs23(3)*dt
        
	call rhs(ne,tt,uu,r3)
    
    	do k=1,ne 
		uu(k)=u(k)+dt*(abs23(4,1)*r1(k)+abs23(4,2)*r2(k)+abs23(4,3)*r3(k))
        end do
		tt=t+cbs23(4)*dt
        
	call rhs(ne,tt,uu,r4)

	!compute second-order solution
    do k=1,ne 
	uu(k) = u(k) + dt*(bbs2(1)*r1(k)+bbs2(2)*r2(k)+bbs2(3)*r3(k)+bbs2(4)*r4(k))
    end do

	!compute third-order solution
    do k=1,ne 
	r3(k) = u(k) + dt*(bbs3(1)*r1(k)+bbs3(2)*r2(k)+bbs3(3)*r3(k)+bbs3(4)*r4(k))
    end do


!estimated error per unit time, should be at most tol
do k=1,ne
err(k)=dabs(r3(k)-uu(k))/dt
end do

!compute maximum error
est = maxval(err)

!check
if (est.le.tol.or.dt.le.hmin) then !accept
t=t+dt
	do k=1,ne
	u(k)=r3(k)
	end do
i=i+1
print*,t
write(20,*) t,u(1),u(2),u(3)
write(30,*) t,dt
end if

ratio =(sf*tol/(est+1.0d-8))**(1.0d0/2.0d0) 

!adjust step size for next time step
dt=dt*ratio

if(dt.lt.hmin) dt=hmin
if(dt.gt.hmax) dt=hmax
if((t+dt).gt.b) dt=b-t

j=j+1
end do

print*,"number of time steps taken",i
print*,"total number of time steps taken",j
print*,"number of null steps",j-i

return
end





!---------------------------------------------------!
!Butcher Table for the Runge-Kutta Method 
!---------------------------------------------------!
SUBROUTINE RKTABLE
implicit none
integer::i,j
real*8::abs23(4,4),cbs23(4),bbs2(4),bbs3(4)

common /bs23table/abs23,cbs23,bbs2,bbs3

!Ititialize Butcher table
do i=1,4
do j=1,4
abs23(i,j)=0.0d0
end do
cbs23(i)=0.0d0
bbs2(i) =0.0d0
bbs3(i) =0.0d0
end do

!Butcher table for Bogacki-Shampine method (ode23)
abs23(2,1)=1.0d0/2.0d0
abs23(3,1)=0.0d0
abs23(3,2)=3.0d0/4.0d0
abs23(4,1)=2.0d0/9.0d0
abs23(4,2)=3.0d0/9.0d0
abs23(4,3)=4.0d0/9.0d0

cbs23(2)=1.0d0/2.0d0
cbs23(3)=3.0d0/4.0d0
cbs23(4)=1.0d0

bbs2(1)=7.0d0/24.0d0
bbs2(2)=1.0d0/4.0d0
bbs2(3)=1.0d0/3.0d0
bbs2(4)=1.0d0/8.0d0


bbs3(1)=2.0d0/9.0d0
bbs3(2)=3.0d0/9.0d0
bbs3(3)=4.0d0/9.0d0
bbs3(4)=0.0d0

Return 
End


!---------------------------------------------------!
!rhs function (Lorenz system = chaotic attractor)
!---------------------------------------------------!
subroutine rhs(ne,t,u,f)
integer::ne
real*8 ::t,u(ne),f(ne)
real*8 ::r,s,b

r = 28.0d0
s = 10.0d0
b = 8.0d0/3.0d0

f(1) = s*(u(2)-u(1))
f(2) = r*u(1) - u(2) - u(1)*u(3)
f(3) = u(1)*u(2) - b*u(3)

return
end

