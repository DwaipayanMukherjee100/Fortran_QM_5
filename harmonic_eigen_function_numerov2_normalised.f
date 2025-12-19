        program Harmonic_with_numerov
c--- this one is same as previous one , but we shall normalise before the plot 

        implicit double precision (a-h,o-z)
        parameter (N=10000)
        dimension x(N),psi(N)
        E=0.5
        h=1.0d-3
        a=-5.0d0
        b=5.0d0
        h=(b-a)/real(N)
        do i=1,N
                x(i)=i*h+a 
                psi(i)=0.0d0 
        end do 
        psi(2)=h
        open(10,file="numer2.txt",status="replace")
        
        do i=3,N
                z1=2.0d0*(1.0d0-(5.0d0*h*h*g(x(i-1),E)/12.0d0))*psi(i-1)
                z2=(1+(h*h*g(x(i-2),E)/12.0d0))*psi(i-2)
                z3= (1+(h*h*g(x(i),E)/12.0d0))
                if (z3.lt.1.0d-5) then
                write(*,*)"er",i 
                end if 
                psi(i)=(z1-z2)/z3
        end do 
        A=0.0d0 
        do i=1,N
c---- simple reimann sum
                A=A+psi(i)*psi(i)*h 
        end do 
        A=sqrt(A)
        do i =1,N
                write(10,*) x(i),psi(i)/A 
        end do 
        stop 
        end
        double precision function g(x,E)
        implicit double precision (a-h,o-z)
        g=2.0d0*(E-V(x))
        return 
        end 
        double precision function V(x)
        implicit double precision (a-h,o-z)
        g1=1.0d0
        V=0.50d0*g1*x*x
        return 
        end  
