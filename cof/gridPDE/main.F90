!Main param module
!ifort -coarray=shared -coarray-num-images=12 main.F90
module params
    implicit none
    integer, parameter :: dp = kind(1.0D0)
    integer, parameter :: cp = dp  !cp=current precision,"cp" or "dp"
    integer, parameter :: Nxp = 5, Nyp = 4 !Number of cells per core in each dimension
    integer, parameter :: Px = 3 !3 Procs in x dimension
    integer, parameter :: Ng = 1 !Number of ghost cells
    real(cp), parameter :: xMin=-1.0, xMax=1.0, yMin=-1.0,yMax=1.0 !Grid domain

    real(cp), parameter :: qNorth = 0.0, qSouth = 0.0, qWest = 0.0, qEast = 0.0

    integer :: isd = 1-Ng, ied = Nxp+Ng
    integer :: jsd = 1-Ng, jed = Nyp+Ng
    integer :: is = 1, ie = Nxp
    integer :: js = 1, je = Nyp

end module params

module gridOps
    use params
    implicit none

    real(cp) :: dx,dy,dxP,dyP !Grid spacing, assumed uniform
    real(cp) :: xMinP, xMaxP, yMinP, yMaxP

    real(cp), allocatable :: Q(:,:)[:,:]
    real(cp), allocatable :: xc(:)[:,:], yc(:)[:,:] !Local grid info
    integer :: myIDx, myIDy, NumX, NumY
    contains

    subroutine initGrid()
        integer :: n, gID
        integer :: gridID(2), gridShape(2)

        allocate( Q(isd:ied,jsd:jed)[Px,*] )
        allocate(xc(isd:ied)[Px,*])
        allocate(yc(jsd:jed)[Px,*])
        !Implicit sync b/c allocate

        !Construct cell-centers for grid
        gID = this_image()
        gridID = this_image(Q)
        do n=1,2
            gridShape(n) = ucobound(Q,n) - lcobound(Q,n) + 1
        end do

        myIDx = gridID(1)
        myIDy = gridID(2)
        NumX  = gridShape(1)
        NumY  = gridShape(2)
        
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

        !This procs chunk
        xc = xMinP + 0.5*dx + dx*(/isd-1:ied-1/)
        yc = yMinP + 0.5*dy + dy*(/jsd-1:jed-1/)

        critical !Coarray throttle
            write(*,'(a,I4,a,I4,a)') 'My rank is (', myIDx, ',', myIDy, ')'
            write(*,'(a,f7.3)') '    xMin = ', xMinP
            write(*,'(a,f7.3)') '    xMax = ', xMaxP
            write(*,'(a,f7.3)') '    yMin = ', yMinP
            write(*,'(a,f7.3)') '    yMax = ', yMaxP
            !write(*,*) 'xc = ', xc
            !write(*,*) 'yc = ', yc
        end critical 

        Q(:,:) = 1.0*gID

        sync all
    end subroutine initGrid

    subroutine destroyGrid()
        deallocate(Q,xc,yc)
    end subroutine destroyGrid

    !Does halo swap
    subroutine Halo()

        sync all !Assuming need to sync at beginning
        !West boundary
        if (myIDx == 1) then
            Q(isd:is-1,:) = qWest
        else
            Q(isd:is-1,:) = Q(ie-Ng+1:ie,:)[myIDx-1,myIDy]
        endif
        !East boundary
        if (myIDx == NumX) then
            Q(ie+1:ied,:) = qEast
        else
            Q(ie+1:ied,:) = Q(is:is+Ng-1,:)[myIDx+1,myIDy]
        endif
        !South boundary
        if (myIDy == 1) then
            Q(:,jsd:js-1) = qSouth
        else
            Q(:,jsd:js-1) = Q(:,je-Ng+1:je)[myIDx,myIDy-1]
        endif
        !North boundary
        if (myIDy == 1) then
            Q(:,je+1:jed) = qNorth
        else
            Q(:,je+1:jed) = Q(:,js:js+Ng-1)[myIDx,myIDy+1]
        endif

        !Final sync
        sync all
    end subroutine Halo
end module gridOps

!module pdeOps
!end module pdeOps

program Main
    use params
    use gridOps

    implicit none

    integer :: myID, NumP, i,j
    integer :: gridShape(2)
    myID = this_image() !1D rank
    NumP = num_images()
    if (myID == 1) then
        write(*,*) 'Dimensions per core: ', Nxp, Nyp
        write(*,*) 'Number of cores = ', NumP
    endif

    call initGrid()
    call Halo()
    critical
        write(*,*) 'I am rank ', myID
        do j=jsd,jed
            do i=isd,ied
                write(*,*) 'i,j,Q(i,j) = ', i, ',', j, ',', Q(i,j)
            end do
        end do
    end critical
    call destroyGrid()

end program Main
