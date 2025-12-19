        program Harmonic_with_numerov 
        implicit double precision (a-h,o-z)
c---- now we have created a function out of the last code which returns the end point of the 
c ---- wavefunction , such that we can use bisection method to find it's eigen values 
        tol=1.0d-12
        Emin=1.0d0 
        Emax=2.4d0  
        do i =1,1000 
        Emid=(Emin+Emax)/2.0d0
        z1=xnumerend(Emin,0)
        z2=xnumerend(Emid,0)
        z3=xnumerend(Emax,0)
        if(abs(Emax-Emin).lt.tol) then 
        write(*,*)"Eigenvalue is : ",Emid
        z2=xnumerend(Emid,1)
        exit 
        
        else if (z1*z2 .lt. 0.0d0) then
        Emax = Emid
        else
        Emin = Emid
        end if  
        end do 
        stop 
        end 
        double precision function xnumerend(E,nflag)
        implicit double precision (a-h,o-z)
        parameter (N=10000)
        dimension x(N),psi(N)
        h=1.0d-3
        a=-5.0d0
        b=5.0d0
        h=(b-a)/real(N)
        do i=1,N
                x(i)=i*h+a 
                psi(i)=0.0d0 
        end do 
        psi(2)=h
        
        
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
        if(nflag.eq.1) then 
        open(10,file="numer3.txt",status="replace")
        do i =1,N
                write(10,*) x(i),psi(i)/A 
        end do 
        close(10)
        end if 
        xnumerend=psi(N)
        return 
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