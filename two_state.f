c--- this code solves for a two state system, 
        program timeevo
        implicit double precision(a-g,q-z)
        implicit double complex (h,o-p)
        Parameter (N=1000000)
        dimension Psi1(N),Psi2(N),H(2,2)
        dimension A1(N),A2(N)
        dt=0.001d0
        Psi1(1)=(1.0d0,0.0d0)
        Psi2(1)=(0.0d0,0.0d0)
        dB=0.1d0
        open(10,file="te.txt",status="unknown")
        open(20,file="te2.txt",status="unknown")
        do k=1,100
        H(1,1)=-0.5d0*(0.0d0,0.0d0)
        H(1,2)=-0.5d0*(1.0d0,0.0d0)*k*dB
        H(2,1)=-0.5d0*(1.0d0,0.0d0)*k*dB
        H(2,2)=-0.5d0*(0.0d0,0.0d0)
        do i=2,N
        Psi1(i)=Psi1(i-1)+dt*((-1.0d0*(0.0d0,1.0d0))*
     &  (H(1,1)*Psi1(i-1)+H(1,2)*Psi2(i-1)))
        Psi2(i)=Psi2(i-1)+dt*((-1.0d0*(0.0d0,1.0d0))*
     &  (H(2,1)*Psi1(i-1)+H(2,2)*Psi2(i-1)))
        end do 
        do i=1,N
        A1(i)=abs(psi1(i))**2
        A2(i)=abs(psi2(i))**2
        end do 
        if(k.eq.50) then 
        do i =1,N/100
        write(20,*)i*dt,A1(i),A2(i)
        end do 
        end if
c----- As we ARE WORKING WITH PROBABILITY , WE HAVE MAXIMUM PROBABILITY OF 1 AND MINIMUM OF 0.
c----- this portion is to confirm that the frequency of any of the states is proportional to the magnetic 
c----  feild applied to it. we find the frequency but analysing the postion of consecutive crests in the data
c---- as we know the solutions are sinusoidal  
        flag=0
        t1=0
        t2=0
        j=0
        DO  I=1,N
        if(A1(i-1).lt.A1(i).and.A1(i+1).lt.A1(i)) then
        if(flag.eq.0) then 
        t1=i*dt
        flag=1
        else
        t2=i*dt
        exit
        end if 
        end if 
        END DO 
        f=1.0d0/(t2-t1)
        write(10,*)f,k*db
        end do 
        stop
        end
