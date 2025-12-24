        program Harmonic_with_numerov 
        implicit double precision (a-h,o-z)
c---- now we have created a function out of the last code which returns the end point of the 
c ---- wavefunction , such that we can use Newton-raphson method to find it's eigen values 
c-----  THis works better than i even expected honestly 
c----- the numerical accuracy is much better than bisection.
        do k=0,1
        E_guess=-0.2-k
        dE=1.0d-5
        Tol=10d-12
        do i=1,50
        fde=xnumerend(E_guess+dE,0)
        fe=xnumerend(E_guess,0)
        fde=(fde-fe)/dE 
        z=fe/fde
        if(abs(z).lt.Tol)  then
        write(*,*) "Eigenvalue =",E_guess 
        p=xnumerend(E_guess,1+k)
        exit  
        endif 
        E_guess=E_guess-z
        End do 
        end do 
        stop 
        end 
        double precision function xnumerend(E,nflag)
        implicit double precision (a-h,o-z)
        parameter (N=10000)
        dimension x(N),psi(N)
        h=1.0d-3
        a=0.0001d0
        b=10.00d0
        h=(b-a)/dble(N)
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
        open(10,file="raphson1.txt",status="replace")
        do i =1,N
               write(10,*) x(i),psi(i)/A 
        end do 
        close(10)
        end if 
          if(nflag.eq.2) then 
        open(20,file="raphson2.txt",status="replace")
        do i =1,N
               write(20,*) x(i),psi(i)/A 
        end do 
        close(20)
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
        lambda=0.0005d0
        g1=exp(-lambda*x)
        V=-g1/x 
        return 
        end  