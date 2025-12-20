c--- this code will attempt to make a transmission probability vs Energy plot for a finite potential barrier problem 
c--- mathematics of this is stored in a corresponding seperate pdf different from main  QM_notes 
	        program transmission
        implicit double precision (a-h,r-z)
        z=rk4(2.0d0,2.0d0)
        stop
        end

        double precision function v(x,a)
        implicit double precision (a-h,r-z)
        V0=2.0d0
        if(abs(x).le.a)then
        v=V0
        else
        v=0.0d0
        endif
        return
        end

        double precision function rk4(a,E)
        implicit double precision (a-h,r-z)
        implicit double complex (p,q)
        parameter(N=200000)
        dimension psi(N),psid(N),x(N)

        x_min=-5.0d0
        x_max=5.0d0
        h=(x_max-x_min)/dble(N)

        do i=1,N
        x(i)=x_min+(i-1)*h
        psi(i)=(0.0d0,0.0d0)
        psid(i)=(0.0d0,0.0d0)
        enddo

        k=sqrt(2.0d0*E)
        psi(1)=exp((0.0d0,1.0d0)*k*x(1))
        psid(1)=(0.0d0,1.0d0)*k*psi(1)

        open(10,file='trans1.txt',status='replace')
        write(10,*)x(1),abs(psi(1))**2

        do i=1,N-1
        k1_psi=psid(i)
        k1_psid=2.0d0*(v(x(i),a)-E)*psi(i)

        k2_psi=psid(i)+0.5d0*h*k1_psid
        k2_psid=2.0d0*(v(x(i)+0.5d0*h,a)-E)*
     &           (psi(i)+0.5d0*h*k1_psi)

        k3_psi=psid(i)+0.5d0*h*k2_psid
        k3_psid=2.0d0*(v(x(i)+0.5d0*h,a)-E)*
     &           (psi(i)+0.5d0*h*k2_psi)

        k4_psi=psid(i)+h*k3_psid
        k4_psid=2.0d0*(v(x(i)+h,a)-E)*
     &           (psi(i)+h*k3_psi)

        psi(i+1)=psi(i)+(h/6.0d0)*
     &          (k1_psi+2.0d0*k2_psi+
     &           2.0d0*k3_psi+k4_psi)

        psid(i+1)=psid(i)+(h/6.0d0)*
     &           (k1_psid+2.0d0*k2_psid+
     &            2.0d0*k3_psid+k4_psid)

        write(10,*)x(i+1),abs(psi(i+1))**2
        enddo

        close(10)
        rk4=0.0d0
        return
        end
