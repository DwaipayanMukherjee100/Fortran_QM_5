        program Harmonic_with_numerov
c-- this code solved the Time-independent schrodinger equation using numerov method 
c--- the next version of this code will aim at finding the eigen vales of a given harmonic oscilator 
c---- in this entire repo , we shall assume \bar{h} =1, m=1 and \omega=1

        implicit double precision (a-h,o-z)
        parameter (N=100000)
        dimension x(N),psi(N)
c----- Checked for ground state convergence , Warning: if you change the values of a and b , 
c----- Numerov accumulate errors and diverges , work with only the given value here.
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
        open(10,file="numer1.txt",status="replace")
        write(10,*)x(1),psi(1)
        write(10,*)x(2),psi(2) 
        do i=3,N
                z1=2.0d0*(1.0d0-(5.0d0*h*h*g(x(i-1),E)/12.0d0))*psi(i-1)
                z2=(1+(h*h*g(x(i-2),E)/12.0d0))*psi(i-2)
                z3= (1+(h*h*g(x(i),E)/12.0d0))
                if (z3.lt.1.0d-5) then
                write(*,*)"er",i 
                end if 
                psi(i)=(z1-z2)/z3
                write(10,*)x(i),psi(i)
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
