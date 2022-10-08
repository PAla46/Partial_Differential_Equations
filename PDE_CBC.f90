! Solving partial differential equations using Cauchy Boundary equations

Program laplace
    IMPLICIT NONE
    Integer,parameter :: lx=34,ly=34
    double precision :: old_temp(1:lx,1:ly), temp(1:lx,1:ly)
    Integer :: i,j,ii,jj,kk
    double precision :: bound_temp, increment_temp, dx,dy,prefactor

    Character(len=30) :: filename
    Integer :: iunit,ci,cii,test,counter

    bound_temp = 0.0d0; increment_temp = 0.1d0

    old_temp = 0.0d0
    
    do j=1,ly
        old_temp(1,j) = 7.4d0
        old_temp(34,j) = 0.8d0
    end do
    do i=1,lx
        old_temp(i,1) = 7.5d0 - dfloat(i)*increment_temp
        old_temp(i,34) = 7.5d0 - dfloat(i)*increment_temp
    end do

    temp = old_temp

    iunit = 71; ci=0
    write(filename,'("initialize_",i0,".dat")')ci 
    open(newunit = iunit,file=filename)
    do ii=1,lx
        do jj=1,ly
            write(iunit,*) ii,jj,old_temp(ii,jj)
        end do
    end do 
    close(iunit)

    test=0;counter=0
    do
        counter=counter+1
        test=0
        do jj=2,ly-1
            do ii=2,lx-1
                temp(ii,jj) = 0.25d0*(old_temp(ii-1,jj) + old_temp(ii+1,jj) + old_temp(ii,jj-1) + old_temp(ii,jj+1))
            end do
        end do

        do jj=2,ly-1
            do ii=2,lx-1
                if((abs(temp(jj,ii) - old_temp(jj,ii))).gt.0.0001) test=1
            end do
        end do

        if(test.eq.0) exit
        old_temp = temp
    end do

    write(*,*) 'counter',counter


    iunit=60;ci=counter
    write(filename,'("initialize_",i0,".dat")') ci
    open(newunit=iunit,file=filename)
    do ii=1,lx
        do jj=1,ly
            write(iunit,*) ii,jj,temp(ii,jj)
            if(ii.eq.20.and.jj.eq.20)then
                print*,'T at (20,20) after convergence is', temp(ii,jj)
            end if
            if(ii.eq.14.and.jj.eq.11)then
                print*,'T at (14,11) after convergence is', temp(ii,jj)
            end if
            if(ii.eq.12.and.jj.eq.9)then
                print*,'T at (12,9) after convergence is', temp(ii,jj)
            end if
            if(ii.eq.17.and.jj.eq.17)then
                print*,'T at (17,17) after convergence is', temp(ii,jj)
            end if
        end do
    end do
    close(iunit)

End Program laplace
