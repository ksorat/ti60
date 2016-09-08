!Main param module
!ifort -coarray=shared -coarray-num-images=16 main.F90
module params
    implicit none
    integer, parameter :: dp = kind(1.0D0)
    integer, parameter :: cp = dp  !cp=current precision,"cp" or "dp"
    integer, parameter :: NxTot = 96, NyTot = 96
    
    integer, parameter :: Ng = 1 !Number of ghost cells
    
    real(cp), parameter :: qEdge = -10, qInterior = 10.0
    real(cp), parameter :: qNorth = qEdge, qSouth = qEdge, qWest = qEdge, qEast = qEdge
    real(cp), parameter :: resTol = 1.0e-1
    real(cp) :: k0 = 1.0

    integer :: Px[*]  !Procs in x dimension
    integer :: Nxp, Nyp !Number of cells per core in each dimension

end module params

module gridOps
    use params
    implicit none
    integer :: isd, ied, is, ie
    integer :: jsd,jed,js,je


    real(cp) :: xMin,xMax,yMin,yMax !Grid domain
    real(cp) :: dt, dx,dy,dxP,dyP !Grid spacing, assumed uniform
    real(cp) :: xMinP, xMaxP, yMinP, yMaxP

    real(cp), allocatable :: Q(:,:)[:,:]
    real(cp), allocatable :: xc(:)[:,:], yc(:)[:,:] !Local grid info
    integer :: myIDx, myIDy, NumX, NumY
    contains

    subroutine initGrid()
        integer :: n, gID
        integer :: gridID(2), gridShape(2)

        isd = 1-Ng
        ied = Nxp+Ng
        jsd = 1-Ng
        jed = Nyp+Ng
        is = 1
        ie = Nxp
        js = 1
        je = Nyp

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

        !Timestep
        dt = 0.25*( min(dx,dy)**2.0/k0 )

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
            write(*,'(a,f9.3)') '    xMin = ', xMinP
            write(*,'(a,f9.3)') '    xMax = ', xMaxP
            write(*,'(a,f9.3)') '    yMin = ', yMinP
            write(*,'(a,f9.3)') '    yMax = ', yMaxP
            !write(*,*) 'xc = ', xc
            !write(*,*) 'yc = ', yc
        end critical 

        Q(:,:) = qInterior

        sync all
        call Halo()

    end subroutine initGrid

    subroutine destroyGrid()
        deallocate(Q,xc,yc)
    end subroutine destroyGrid

    !Does halo swap
    subroutine Halo()

        sync all !Need to sync at beginning
        !East/West
        if ( (myIDx > 1) .and. (myIDx < NumX) ) then
            !Interior, talk to neighbors
            Q(isd :is-1 , js:je) = Q(ie-Ng+1:ie      , js:je)[myIDx-1,myIDy]
            Q(ie+1:ied  , js:je) = Q(is     :is+Ng-1 , js:je)[myIDx+1,myIDy]
        else if (myIDx == 1) then 
            !West boundary
            Q(isd:is-1,js:je) = qWest
        else 
            !East boundary
            Q(ie+1:ied,js:je) = qEast
        endif

        !North/South
        if ( (myIDy > 1) .and. (myIDy < NumY) ) then
            !Interior, talk to neighbors
            Q(is:ie,jsd : js-1) = Q(is:ie , je-Ng+1:je)     [myIDx,myIDy-1]
            Q(is:ie,je+1: jed)  = Q(is:ie , js     :js+Ng-1)[myIDx,myIDy+1]
        else if (myIDy == 1) then
            !South boundary
            Q(is:ie,jsd:js-1) = qSouth
        else
            !North boundary
            Q(is:ie,je+1:jed) = qNorth
        endif

        !Final sync
        sync all
    end subroutine Halo

    !Evolve local grid using relaxation
    !Assume halos up to date
    function Relax()
        integer :: i,j
        real(cp) :: Relax
        real(cp) :: qSwap(is:ie,js:je)
        real(cp) :: wx,wy,wxx,wyy,rk

        wx = k0*dt/(dy*dy)
        wy = k0*dt/(dx*dx)
        wxx = -2 + dy*dy/(2*k0*dt)
        wyy = -2 + dx*dx/(2*k0*dt)


        do j=js,je
            do i=is,ie
                qSwap(i,j) = wx* ( Q(i-1,j) + Q(i+1,j) + wxx*Q(i,j) ) &
                           + wy* ( Q(i,j-1) + Q(i,j+1) + wyy*Q(i,j) ) 
            end do
        end do

        Relax = 0.0
        do j=js,je
            do i=is,ie
                rk = Q(i,j) - qSwap(i,j)
                Relax = Relax + rk*rk
                Q(i,j) = qSwap(i,j)
            end do
        end do

    end function Relax
    
end module gridOps

!module pdeOps
!end module pdeOps

program Main
    use params
    use gridOps

    implicit none

    integer :: myID, NumP, i,j, ts
    integer :: gridShape(2)
    character(len=500) :: inpArg
    real(cp), codimension[*] :: locRes
    real(cp) :: totRes !Total residual
    logical :: doIter = .true.
    myID = this_image() !1D rank
    NumP = num_images()

    !Initialize grid info
    if (myID == 1) then
        call get_command_argument(1,inpArg)
        write(*,*) inpArg
        read(inpArg,*) Px
        do i=1,NumP
            Px[i] = Px
        enddo
    endif
    sync all

    Nxp = NxTot/Px
    Nyp = NyTot/(NumP/Px)

    xMin = 0.0
    xMax = dble(NxTot)
    yMin = 0.0
    yMax = dble(NyTot)

    call initGrid()
    if (myID == 1) then
        write(*,*) 'Physical Grid dimension: ', NxTot,NyTot
        write(*,*) 'Processor Grid dimension: ', Px, NumP/Px
        write(*,*) 'Dimensions per core: ', Nxp, Nyp
        write(*,*) 'Number of cores = ', NumP
        write(*,*) 'Timestep = ', dt
    endif

    totRes = 1.0e+8
    ts = 0
    
    !Grid initialized, halos initialized, sync complete
    do while(doIter)

        locRes[myID] = Relax()
        if (myID == 1) then
            totRes = 0.0
            do i=1,NumP
                totRes = totRes + locRes[i]
            enddo
            totRes = sqrt(totRes)
            do i=1,NumP
                locRes[i] = totRes
            enddo
            !locRes[:] = totRes
            write(*,*) 'TS = ',ts
            write(*,*) '   Residual = ', totRes
        endif
        !call co_broadcast(totRes,source_image=1)
        sync all
        if (locRes<=resTol) then
            doIter = .false.
        else
            call Halo()
        endif
        ts = ts+1
    enddo

    call destroyGrid()

end program Main
