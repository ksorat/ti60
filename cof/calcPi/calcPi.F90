program calcPi

	implicit none
	integer, parameter :: N = 1000000

	real, dimension(N) :: x,y, r
	real :: piApprox
	integer :: myID, Np, i, totIn
	integer, codimension[*] :: locIn

	myID = this_image()
	Np = num_images()

	!Create proc-dependent seed
	call random_seed(put=(/myID/))

	!Generate random numbers
	call random_number(x)
	call random_number(y)
	

	r = sqrt( x**2.0 + y**2.0 )
	locIn = count( r <= 1 )
	sync all

	if (myID == 1) then
		totIn = 0
		do i=1,Np
			totIn = totIn + locIn[i]
		end do
		piApprox = 4.0*totIn/(Np*N)
		write(*,'(a,f)') 'Pi = ', piApprox
		write(*,'(a,i)') 'Np = ', Np
		write(*,'(a,i)') 'N  = ', N

	endif

end program calcPi


