program hello_coarray
	
	use ifport ! Used to get hostnames
	implicit none

	character(MAX_HOSTNAM_LENGTH + 1) :: hostname
	integer :: istat
	integer, parameter :: max_str_len = 100
	character(max_str_len), allocatable :: greeting[:]
	integer :: image

	! get the host we are running on
	istat = HOSTNAM(hostname)

	! allocate and write greeting message
	allocate(greeting[*])
	write(greeting,"(2(a,i2),a,a,a)") "Image ", this_image(),&
	" of ", num_images(), " (", trim(hostname), ")"

	! wait for all images to write their message
	sync all
	
	! have image 1 write all the messages
	if (this_image() == 1) then
		do concurrent(image=1:num_images())
			write(*,'(a)') trim(greeting[image])
		end do
	endif

end program hello_coarray