
!-----------------------------------------------------------------------------!
! MAE 5093 DOCUMENTATION | Engineering Numerical Analysis
!-----------------------------------------------------------------------------!
! >>> Multigrid Method for solving 1D Poisson equation
!     d2u/dx2 = f(x)
!     Drichlet b.c.代码非常清晰
!     
!-----------------------------------------------------------------------------!
! References: 
! * Fundamentals of Engineering Numerical Analysis by P. Moin (2012) 
! * Numerical Recipes: The Art of Scientific Computing, 2nd Edition (1992) 
!-----------------------------------------------------------------------------!
! Written by Omer San
!            CFDLab, Oklahoma State University, cfdlab.osu@gmail.com
!            www.cfdlab.org
! 
! Last updated: Oct. 29, 2015
! update by Chunye Gong, 2024.02.21
!-----------------------------------------------------------------------------!



program poisson1d
implicit none
integer::i,nitr,nx
real*8 ::dx,x0,xL,tol,pi
real*8,allocatable ::f(:),u(:),x(:)

!Domain
x0 = 0.0d0 !left
xL = 1.0d0 !right

!number of points
nx = 64

!grid spacing (spatial)
dx = (xL-x0)/dfloat(nx)


!spatial coordinates 
allocate(x(0:nx))
do i=0,nx
x(i) = x0 + dfloat(i)*dx
end do

!given source term 
allocate(f(0:nx))

pi = 4.0d0*datan(1.0d0)

!右端项,sin(pi*x)+sin(16pi*x)
do i=0,nx
	f(i) = 0.5d0*(dsin(pi*x(i))+dsin(16.0d0*pi*x(i)))
end do


open(19,file='residual.plt')
write(19,*) 'variables ="n","r"'


allocate(u(0:nx))

!initial guess which satisfies boundary condtions
do i=0,nx
u(i) = 0.0d0
end do

!Tolerance
tol = 1.0d-6

 

call GuassSeidel(nx, f,u,dx, tol)  

return

!initial guess which satisfies boundary condtions
do i=0,nx
u(i) = 0.0d0
end do



!Gauss Seidel scheme
!call GS(nx,dx,f,u,tol,nitr) 

!V-cycle Multigrid Scheme (4-level)
call MG4(nx,dx,f,u,tol,nitr) 

close(19)


print*,"--------------------------------------------"
print*,"MG iterations:", nitr
print*,"--------------------------------------------"


open(11,file='field.plt')
write(11,*) 'variables ="x","u"'
do i=0,nx
write(11,*)x(i),u(i)
end do


end

!-----------------------------------------------------------------------------!
!Gauss-Seidel iteration
!Homogeneous Drichlet B.C.
!i.e., boundary points are not updated, GS迭代仍然有问题 2024.02.21
!-----------------------------------------------------------------------------!
subroutine GuassSeidel(nx,f,u, dx, tol)
implicit none
integer::nx,i, it
real*8 ::dx, res, res1, res2, tol
real*8 ::f(0:nx) 
real*8 ::u(0:nx),u1(0:nx)
res=100

do i=0,nx
	u1(i) = 1d0
end do

do it = 1, 1000000
	res1 = 0
	res1 =0
	do i=1,nx-1
		u1(i) = 0.5d0*(u(i+1)+u(i-1)-dx*dx*f(i))
	end do
	u1(0) = 0.0d0
	u1(nx)= 0.0d0
	

	do i=1,nx-1
		!res1 = max(abs(u1(i)-u(i)), res1)
		res2 = max(abs(f(i) - (u(i+1)-2.0d0*u(i)+u(i-1))/(dx*dx)), res1)		
		u(i) = u1(i)
	end do
	
	res = min(res,res2)
	!print *, res, res1
	if(mod(it,100) .eq. 0) print *, it, res
	if (res .lt. tol) then
		print*,"--------------------------------------------"
		print*,"GS iterations:", it, res
		print*,"--------------------------------------------"
		return
	end if
end do
return

end 


!-----------------------------------------------------------------------------!
!Gauss-Seidel iterative scheme
!-----------------------------------------------------------------------------!
subroutine GS(nx,dx,f,u,tol,nitr)  
implicit none
integer::nx,k,nI,nitr
real*8 ::dx,tol,res
real*8 ::f(0:nx),u(0:nx),r(0:nx)

nI = 90000000 

	!compute initial residual
	call residual(nx,dx,f,u,r)
    res = maxval(dabs(r))
    write(19,*) 0, res
    write(*,*) 0, res
    
!Iterations:
do k=1,nI

	!GS relaxation scheme
	call relax(nx,dx,f,u) 

	!compute residual
	call residual(nx,dx,f,u,r)

    !compute max residual and check with tolerance
    res = maxval(dabs(r))

    write(19,*) k, res
    write(*,*) k, res
     
    if(res.le.tol) goto 100
       
end do
100 continue

nitr = k

end 


!-----------------------------------------------------------------------------!
!V-cycle multigrid scheme (4 level)
!-----------------------------------------------------------------------------!
subroutine MG4(nx,dx,f,u,tol,nitr)  
implicit none
integer::nx,k,nI,nitr,n,i
real*8 ::dx,tol,res
real*8 ::f(0:nx),u(0:nx)
integer::nx2,nx3,nx4
real*8 ::dx2,dx3,dx4
integer::v1,v2,v3
real*8,allocatable::r(:),r2(:),r3(:),r4(:)
real*8,allocatable::p(:),p2(:),p3(:)
real*8,allocatable::u2(:),u3(:),u4(:)
real*8,allocatable::f2(:),f3(:),f4(:)

nI = 100000 !maximum number of outer iteration
v1 = 2   	!number of relaxation for restriction in V-cycle
v2 = 2   	!number of relaxation for prolongation in V-cycle
v3 = 100 	!number of relaxation at coarsest level

nx2 = nx/2
nx3 = nx/4
nx4 = nx/8

dx2 = dx*2.0d0
dx3 = dx*4.0d0
dx4 = dx*8.0d0

allocate(r(0:nx))
allocate(r2(0:nx2))
allocate(r3(0:nx3))
allocate(r4(0:nx4))

allocate(p(0:nx))
allocate(p2(0:nx2))
allocate(p3(0:nx3))

allocate(u2(0:nx2))
allocate(u3(0:nx3))
allocate(u4(0:nx4))

allocate(f2(0:nx2))
allocate(f3(0:nx3))
allocate(f4(0:nx4))



!compute initial residual
call residual(nx,dx,f,u,r)
res = maxval(dabs(r))
write(19,*) 0, res
write(*,*) 0, res
	
    
!Iterations:
do k=1,nI

	!relax v1 times on level-h
    do n=1,v1
   		!GS relaxation scheme，松驰,其实就是GS迭代
		call relax(nx,dx,f,u)  
    end do
    
	!compute residual残差 on level-h
	call residual(nx,dx,f,u,r)

	!compute max residual and check with tolerance
	res = maxval(dabs(r))
	write(19,*) k, res
	!第几次迭代及残差，简单一维问题会收敛很快
	write(*,*) k, res
	!满足就终止循环
	if(res.le.tol) goto 100


    !perform restriction to level-2h (h --> 2h)，细网格到粗网格
	!把残差r，映射到新网络f2上。不是把值映射。类似高阶计算方法？
	!一是粗化，二是对残差进行计算
    call restriction(nx,r,nx2,f2)
    
    !initial guess for level-2h
	do i=0,nx2  
    u2(i) = 0.0d0
    end do

    !relax v1 times on level-2h
    do n=1,v1
   		!GS relaxation scheme，f2是右端项，u2是GS迭代结果
		call relax(nx2,dx2,f2,u2)  
    end do
    
    !compute residual on level-2h,r2是残差
	call residual(nx2,dx2,f2,u2,r2)

    !perform restriction to level-4h (2h --> 4h)，双映射到第3层粗网格
	!因为这里的网格完善，可以直接映射，如果不整除怎么弄？
    call restriction(nx2,r2,nx3,f3)  

    !initial guess for level-4h
    do i=0,nx3 
    u3(i) = 0.0d0 
    end do

    !relax v1 times on level-4h
    do n=1,v1
   		!GS relaxation scheme
		call relax(nx3,dx3,f3,u3)  
    end do
    
    !compute residual on level-4h
	call residual(nx3,dx3,f3,u3,r3)

    !perform restriction to level-8h (4h --> 8h)，映射到第4层网格
    call restriction(nx3,r3,nx4,f4)  

    !initial guess for level-8h
    do i=0,nx4 
    u4(i) = 0.0d0 
    end do
        
    !relax v3 times on level-8h (coarsest grid)
    do n=1,v3
   		!GS relaxation scheme，解是u4，f4是右端项。并不是我们想象的把解在做插值
		call relax(nx4,dx4,f4,u4)  
    end do

	!perform prolongation (8h --> 4h)，u4->p3
    call prolongation(nx4,u4,nx3,p3)   
   
    !correct on level-4h,用u3更新p3. p3其实是误差，根据f3得到u3
    do i=0,nx3 
    u3(i) = u3(i) + p3(i)
    end do

    !relax v2 times on level-4h,双光滑两次
    do n=1,v2
   		!GS relaxation scheme
		call relax(nx3,dx3,f3,u3)  
    end do      
 
	!perform prolongation延长;延伸; (4h --> 2h),u3扩展到p2,p2为u2的误差
    call prolongation(nx3,u3,nx2,p2)   
   
    !correct on level-2h，更新u2
    do i=0,nx2 
    u2(i) = u2(i) + p2(i)
    end do

    !relax v2 times on level-2h,,双光滑两次
    do n=1,v2
   		!GS relaxation scheme
		call relax(nx2,dx2,f2,u2)  
    end do   
     
 	!perform prolongation (2h --> h)，把u2扩展到p
    call prolongation(nx2,u2,nx,p)   
 
    !correct on level-h，更新u
    do i=0,nx
    u(i) = u(i) + p(i)
    end do

    !relax v2 times on level-h ,对u再光滑二次
    do n=1,v2
   		!GS relaxation scheme
		call relax(nx,dx,f,u)  
    end do  
          
end do !k

100 continue

nitr = k

end 




!-----------------------------------------------------------------------------!
!Gauss-Seidel relaxations to update u
!Homogeneous Drichlet B.C.
!i.e., boundary points are not updated, 也可以叫光滑
!-----------------------------------------------------------------------------!
subroutine relax(nx,dx,f,u)  
implicit none
integer::nx,i
real*8 ::dx
real*8 ::f(0:nx),u(0:nx)

!就是那个离散公式, Gauss-Seidel一次
do i=1,nx-1
u(i) = 0.5d0*(u(i+1)+u(i-1)-dx*dx*f(i))
end do

return
end 

!-----------------------------------------------------------------------------!
!Compute residual，计算残差，这个残差是根据方程在计算，并不是根据Ax=b中的迭代误差在计算，
!因为需要重新计算，而GS根本没有残差更新的过程
!-----------------------------------------------------------------------------!
subroutine residual(nx,dx,f,u,r)
implicit none
integer::nx,i
real*8 ::dx
real*8 ::f(0:nx),u(0:nx),r(0:nx)

!误差，右端项减去左边的离散值。
do i=1,nx-1
	r(i) = f(i) - (u(i+1)-2.0d0*u(i)+u(i-1))/(dx*dx)
end do
r(0) = 0.0d0
r(nx)= 0.0d0

return
end 


!-----------------------------------------------------------------------------!
!restriction operator
!nxf: number of points on finel level 
!nxc: number of points on coarse level
!-----------------------------------------------------------------------------!
subroutine restriction(nxf,rf,nxc,fc)  
implicit none
integer::nxf,nxc,i
real*8 ::rf(0:nxf),fc(0:nxc)

! rf(1,2,3)-> fc(1); rf(3,4,5)-> fc(2); rf(5,6,7)->fc(3)
do i=1,nxc-1
fc(i) = 0.25d0*(rf(2*i-1)+2.0d0*rf(2*i)+rf(2*i+1))
end do

fc(0)   = rf(0)
fc(nxc) = rf(nxf)

return
end 

!-----------------------------------------------------------------------------!
!restriction operator
!nxf: number of points on finel level 
!nxc: number of points on coarse level，把rc按照规则分配给pf
!-----------------------------------------------------------------------------!
subroutine prolongation(nxc,rc,nxf,pf)  
implicit none
integer::nxf,nxc,i
real*8 ::pf(0:nxf),rc(0:nxc)

do i=0,nxc-1
pf(2*i)   = rc(i)
pf(2*i+1) = 0.5d0*(rc(i)+rc(i+1))
end do

pf(nxf)   = rc(nxc)


return
end 






