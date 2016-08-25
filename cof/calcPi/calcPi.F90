program calcPi

	implicit none
	integer, parameter :: N = 1000

	real, dimension(N) :: x,y
	real :: rMin,rMax
	real, dimension(N) :: r
	integer :: myID, Np, i,locIn[*]
	rMin = -1.0
	rMax = 1.0

	myID = this_image()
	Np = num_images()

	call random_seed(put=(/myID/))

	!Generate random numbers
	call random_number(r)
	x = rMin + (rMax-rMin)*r
	call random_number(r)
	y = rMin + (rMax-rMin)*r

	r = sqrt( x**2.0 + y**2.0 )
	locIn = count( r <= 1 )
	if (myID == 1) then
		do concurrent(i=1:Np)
			write(*,*) locIn[i]
		end do
	endif

end program calcPi


