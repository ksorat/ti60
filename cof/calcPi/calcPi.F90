program calcPi

	implicit none
	integer, parameter :: N = 1000

	real, dimension(N) :: x,y
	real :: rMin,rMax, piApprox
	real, dimension(N) :: r
	integer :: myID, Np, i, totIn
	integer, codimension[*] :: locIn
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
	sync all

	if (myID == 1) then
		totIn = 0
		do i=1,Np
			totIn = totIn + locIn[i]
		end do
		piApprox = 4.0*totIn/(Np*N)
		write(*,'(a,f)') 'Pi = ', piApprox
		write(*,'(a,i)') 'Np = ', Np

	endif

end program calcPi


