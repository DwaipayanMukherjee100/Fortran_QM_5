        program sg 
        implicit double precision (a-h,o-z)
        parameter (N=1000)
        dimension z(n)
        external final 
c---- We shall assume the particle was launched with a constant velocity vy and then we figure out y(t) as 
c a function of time, now the particle enters the non-uniform magnetic feild , for practical purposes we shall assume
c Bz= z**2 , and take the numerical derivative (although an overkill but might help if the functional form of Bz changes)
c , then we simply integrate twice to find the displacement in z and finally we obtain the stern gerlach experiment result in the 
c form of the histogram of z axis final points. 
        s=0
        open(10,file="sg.txt",status="unknown")
        do i=1,n
        xk=abs(rand())
        if(xk.ge.0.5) then 
        s=1
        else 
        s=-1
        end if 
        write(10,*)final(s)        
        end do 
        stop 
        end 


        double precision function final(spin)
        implicit double precision (a-h,o-z)
        parameter (N=10000)
        dimension y(n),z1(n),fz(n),bz(n),s(n),t(n),z(n),vz(n)
c-- this function will solve for the final pozition of each particle given the spin 
c--- orientation of + or - 1. 
        xL_b=1.0d0 
        xL_m=2.0d0
        xL_a=1.0d0
        vy=300.0d0
c--- we are working in SI units, where we assume the total span of experiment is 4 meters
        dy=4.0d0/real(N)
        do i =1,N 
        y(i)=(i-1)*dy
        bz(i)=0
        t(i)=y(i)/vy 
        end do 
        dz = 1.0d-4
        vz(1) = 0.0d0
        z(1)  = 0.0d0

        do i = 2, N-1

        if ( y(i).gt.xL_b .and. y(i).lt.(xL_b+xL_m) ) then
        bzp = B(z(i-1)+dz, y(i), xL_b, xL_a, xL_m)
        bzm = B(z(i-1)-dz, y(i), xL_b, xL_a, xL_m)
        fz(i) = spin*(bzp - bzm)/(2.0d0*dz)
        else
        fz(i) = 0.0d0
        end if

        vz(i) = vz(i-1) + fz(i)*(t(i)-t(i-1))
        z(i)  = z(i-1)  + vz(i)*(t(i)-t(i-1))

        end do
        final=z(N-1)
        return
        end 

        double precision function B(z,y,xL_b,xL_a,xL_m)
        double precision z,y,xL_b,xL_a,xL_m
        if(y.gt.xL_b.and.(y.lt.(xL_b+xL_m))) then 
        B=z*1.0d0/(z**2+1.0d0)
        else 
        B=0.0d0
        end if
                 
        return 
        end 
