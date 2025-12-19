        program sg2
        implicit double precision (a-h,o-z)
        do i=1,100
        xk=rand()
        if(xk.lt.0.5) then 
        z1=sg1(-1.0d0)
        else 
        z1=sg1(+1.0d0)
        end if  
        end do

        
        stop 
        end
        double precision function sg1(spin)
        implicit double precision (a-h,o-z)
        parameter (N=1000)
c------ This is an attempt to code the stern gerlach experiment in 3d 
c------ here x(n,1) =x(n), x(n,2)=y(n),x(n,3)=z(n)
        dimension x(n,3),f(n,3),v(n,3),b(n)
c----- let's say the particle is launched from origin and a mass=1 
        vy=100.0d0 
        L_magnet_y=10.0
        t=L_magnet_y/vy 
        dt=t/dble(N)
        dy=L_magnet_y/dble(N)
        dx=dy 
        dz=dy 
        v(1,1)=0.0
        v(1,2)=Vy 
        v(1,3)=0.0
        x(1,1)=0.0
        x(1,2)=0.0
        x(1,3)=0.0
        do i=2,N 
        call force(f(i,1),f(i,2),f(i,3),x(i-1,1),
     & x(i-1,2),x(i-1,3),dx,dy,dz,spin)
        v(i,1)=v(i-1,1)+f(i,1)*dt
        v(i,2)=V(i-1,2)+f(i,2)*dt
        v(i,3)=v(i-1,3)+f(i,3)*dt
        x(i,1)=x(i-1,1)+v(i,1)*dt
        x(i,2)=x(i-1,2)+v(i,2)*dt
        x(i,3)=x(i-1,3)+v(i,3)*dt
        call Bfield(b(i),x(i,1),x(i,2),x(i,3))
        end do 
        sg1=0.0 
        open(10,file="sg1.txt",status="unknown")
        do i=1,N
        write(10,*)x(i,1),x(i,2),x(i,3),v(i,1),v(i,2),v(i,3),b(i)
        end do 

        return
        end 
c-----------
        subroutine force(Fx,Fy,Fz,x,y,z,dx,dy,dz,spin)
        implicit double precision (a-h,o-z)
        bz1=0.0
        bz2=0.0
        call Bfield(bz1,x,y,z+dz)
        call Bfield(bz2,x,y,z-dz)
        Fz=0.50d0*spin*(bz1-bz2)/(2*dz)+rand()*5.0d0
        call Bfield(bz1,x,y+dy,z)
        call Bfield(bz2,x,y-dy,z)
        Fx=0.50d0*spin*(bz1-bz2)/(2*dx)
        call Bfield(bz1,x+dx,y,z)
        call Bfield(bz2,x-dx,y,z)
        Fy=0.50d0*spin*(bz1-bz2)/(2*dx)
        return 
        end 
c ----------
        subroutine Bfield(Bz,x,y,z)
        implicit double precision (a-h,o-z)
        Bz=10.0d0+50.0d0*z
        return
        end 

        