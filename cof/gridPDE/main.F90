!Main param module

module params
    integer, parameter :: dp = kind(1.0D0)
    integer, parameter :: cp = dp  !cp=current precision,"cp" or "dp"
    integer, parameter :: Nxp = 50, Nyp = 30 !Number of cells per core in each dimension
    integer, parameter :: Px = 3 !3 Procs in x dimension
    integer, parameter :: Ng = 1 !Number of ghost cells
    real(cp), parameter :: xMin=-1.0, xMax=1.0, yMin=-1.0,yMax=1.0 !Grid domain

    integer :: isd = 1-Ng, ied = Nxp+Ng
    integer :: jsd = 1-Ng, jed = Nyp+Ng
    integer :: is = 1, ie = Nxp
    integer :: js = 1, je = Nyp

end module params

module gridOps
    use params
    real(cp) :: dx,dy !Grid spacing, assumed uniform
    real(cp), allocatable :: Q(:,:)[:,:]
    real(cp), allocatable :: xi(:)[:,:], yi(:)[:,:] !Local grid info

    subroutine initGrid()
        allocate( Q(isd:ied,jsd:jed)[Px,*] )
        allocate(xi(is:ie)[Px,*])
        allocate(yi(js:je)[Px,*])
        !Implicit sync b/c allocate
    end subroutine initGrid

    subroutine destroyGrid()
        deallocate(Q,xi,yi)
    end subroutine destroyGrid
end module gridOps

!module pdeOps
!end module pdeOps

program Main
    use params

    integer :: myID, myIDx, myIDy, NumP
    integer :: gridShape(2)
    myID = this_image() !1D rank
    NumP = num_images()
    if (myID == 1) then
        write(*,*) 'Dimensions per core: ', Nxp, Nyp
        write(*,*) 'Number of cores = ', NumP
    endif

    gridShape = this_image(Q)

    write(*,*) 'My rank is ', gridShape(1), gridShape(2)
    call initGrid()
end program Main
