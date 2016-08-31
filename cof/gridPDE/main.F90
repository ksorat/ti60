!Main param module
!ifort -coarray=shared -coarray-num-images=12 main.F90
module params
    implicit none
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
    implicit none

    real(cp) :: dx,dy,Pdx,Pdy !Grid spacing, assumed uniform
    real(cp) :: xMinP, xMaxP, yMinP, yMaxP

    real(cp), allocatable :: Q(:,:)[:,:]
    real(cp), allocatable :: xc(:)[:,:], yc(:)[:,:] !Local grid info
    integer :: myIDx, myIDy, NumX, NumY
    contains

    subroutine initGrid()
        integer :: n
        integer :: gridID(2), gridShape(2)

        allocate( Q(isd:ied,jsd:jed)[Px,*] )
        allocate(xc(isd:ied)[Px,*])
        allocate(yc(jsd:jed)[Px,*])
        !Implicit sync b/c allocate

        !Construct cell-centers for grid
        gridID = this_image(Q)
        do n=1,2
            gridShape(n) = ucobound(Q,n) - lcobound(Q,n) + 1
        end do

        myIDx = gridID(1)
        myIDy = gridID(2)
        NumX  = gridShape(1)
        NumY  = gridShape(2)
        write(*,*) 'My rank is ', gridID(1), gridID(2)
        write(*,*) '   of ', gridShape(1), gridShape(2)
        
        !Grid spacing
        dx = (xMax-xMin)/(NumX*Nxp)
        dy = (yMax-yMin)/(NumY*Nyp)

        !Processor spacing
        dxP = dx*Nxp 
        dyP = dy*Nyp

        xMinP = xMin+dxP*(myIDx-1)
        xMaxP = xMinP + dxP
        yMinP = yMin+dyP*(myIDy-1)
        yMaxP = yMinP + dyP

        critical
            write(*,'(a,I,a,I)') 'My rank is (' myIDx, ',', myIDy, ')'
        end critical 
    end subroutine initGrid

    subroutine destroyGrid()
        deallocate(Q,xc,yc)
    end subroutine destroyGrid
end module gridOps

!module pdeOps
!end module pdeOps

program Main
    use params
    use gridOps

    implicit none

    integer :: myID, NumP
    integer :: gridShape(2)
    myID = this_image() !1D rank
    NumP = num_images()
    if (myID == 1) then
        write(*,*) 'Dimensions per core: ', Nxp, Nyp
        write(*,*) 'Number of cores = ', NumP
    endif
    call initGrid()
    
    call destroyGrid()

end program Main
