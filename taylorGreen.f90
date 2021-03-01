program taylorGreen
implicit none
real,parameter::pi=4*ATAN(1.0)
integer,parameter::nx=100 !integer multiples of pi
integer,parameter::ny=100
real,parameter::n=100 !pi factor
!arrays
integer::i,j,iunit,tcount,tn
real::ii,ij,k
real,dimension(0:nx,0:ny)::u,v,p,x,y,vres
real::Vi,z
real::t,nu,diss,dt
character*14::filename,filename2
character*50::filename1
Vi=10
k=1
open(1,file="x_y_V.dat",status='replace')
open(2,file="velocity.plt",status='replace')
open(3,file="x_y_p.dat",status='replace')
open(4,file="pressure.plt",status='replace')
open(5,file="velocity.csv",status='replace')
do i=0,nx
    do j=0,ny
        x(i,j)=(2*pi/nx)*i
        y(i,j)=(2*pi/ny)*j
        ij=j
        u(i,j)=Vi*sin(K*x(i,j))*cos(k*y(i,j))
        v(i,j)=-1*Vi*cos(x(i,j))*sin(y(i,j))
        p(i,j)=-0.25*(cos(2*k*x(i,j))+cos(2*k*y(i,j)))
        vres(i,j)=sqrt((u(i,j)**2)+(v(i,j)**2))
        write(3,*)x(i,j),y(i,j),p(i,j)
        write(1,*)x(i,j),y(i,j),vres(i,j)
        write(5,*)x(i,j),',',y(i,j),',',vres(i,j)
        write(*,*)x(i,j)
    end do
end do
write(2,*)'set xlabel "x"'
write(2,*)'set ylabel "y"'
write(2,*)'set style function pm3d'
write(2,*)'set palette model RGB'
write(2,*)'set palette rgb 33,13,10'
write(2,*)'set zlabel "U velocity"'
write(2,*)'set autoscale xy'
write(2,*)'set title "U plot"'
write(2,*)'plot "x_y_V.dat" using 1:2:3 with image'
write(4,*)'set xlabel "x"'
write(4,*)'set ylabel "y"'
write(4,*)'set style function pm3d'
write(4,*)'set palette model RGB'
write(4,*)'set palette rgb 33,13,10'
write(4,*)'set zlabel "Pressure"'
write(4,*)'set autoscale xy'
write(4,*)'set title "Pressure Field"'
write(4,*)'plot "x_y_p.dat" using 1:2:3 with image'
CALL SYSTEM('gnuplot -p pressure.plt')
CALL SYSTEM('gnuplot -p velocity.plt')
close(1)
close(2)
close(3)
close(4)
close(5)
t=20 !in seconds
dt=0.1 !in seconds
tcount=t/dt
nu=0.001 !kinematic viscosity
diss=(-2*nu*(k**2))
do tn=0,tcount
    open(16,file="tempV.dat",status='replace')    
    z=modulo(tn,10)
    if (z.EQ.0) then
        write(filename,'("time",I0,".dat")')tn 
        open(unit=tn,file=filename,status='replace')
    end if 
    do i=0,nx
        do j=0,ny
            vres(i,j)=vres(i,j)*exp(diss*tn*dt)
            p(i,j)=p(i,j)*exp(diss*tn*dt)
            z=modulo(tn,10)
            if (z.EQ.0) then
                write(tn,*)x(i,j),y(i,j),vres(i,j)
                write(16,*)x(i,j),y(i,j),vres(i,j)
            end if
        end do
    end do
    close(unit=tn)
    
        if (z.EQ.0) then
        open(15,file="timeplot.plt",status='replace')
        write(15,*)'set xlabel "x"'
        write(15,*)'set ylabel "y"'
        write(15,*)'set style function pm3d'
        write(15,*)'set cbrange [-1:10]'
        write(15,*)'set palette model RGB'
        write(15,*)'set palette rgb 33,13,10'
        write(15,*)'set zlabel "V velocity at ',tn*dt,'"'
!        write(15,*)'set autoscale xy'
       ! write(15,*)'set title "V plot"'
        write(15,*)'plot "tempV.dat" using 1:2:3 with image title "Taylor Green Vortex Velocity Field in Domain [0,2*pi] (2D)"'
        write(15,*)'set term png'
        write(filename2,'("time",I0,".png")')tn
        WRITE(15,'(A,A,A)')'set output "',filename2,'"'
        write(15,*)'replot'
        call SYSTEM('gnuplot -p timeplot.plt')
        end if 
    close(15) 
    close(16)  
end do

end program
