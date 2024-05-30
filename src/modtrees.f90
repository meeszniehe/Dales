!> \file modtrees.f90
!! By Meesz Niehe, email: meesz@niehe.com, TU Delft, section Atmospheric Physics 

module modtrees
    use modtreesdata, only : lapply_treedrag, lreadfile_trees, lapply_sourceSGS, ltree, Cd, A_pad
    use modprecision        
    implicit none
    save
    public :: inittrees, exittrees, applytrees 

    contains
        subroutine inittrees
        ! Global variables,itot=xdim, imax=itot/nprocx, kmax=maxheightvaluedimensions(given)
        ! zh=thickness half-level, zf=thickness at full level
        ! ifinput/ifnamopt are integer number to organise file by numbers
        ! i1,j1,k1 = imax+1
        ! ih, jh = 3, kh = 1 ghost cells for extra conditions at, for example, boundaries
        ! dx, dy = gridspacing in x and y direction
        use modglobal,  only :  zh, zf, itot, jtot, i1, j1, k1, ih, jh, imax, jmax, kmax, dx, dy, &
                                ifinput, ifnamopt, fname_options, cexpnr, cu, cv, checknamelisterror
        use modmpi,     only :  myid, comm3d, mpierr,  myidx, myidy, d_mpi_bcast, excjs, &
                                D_MPI_ALLREDUCE, mpi_max, MPI_SUM

        real(field_r), allocatable :: tree_height(:,:)              !< 2D array to store tree heights at each grid point (x,y)             
        integer                 :: i, j, k, ierr, ii, jj, kk, n     !< initialize integer to loop over
        integer                 :: startIdx, endIdx, di, dj         !< initialize integer to loop over for tree_crown creation
        integer, allocatable    :: kindex_tree(:,:)                 !< index of tree height
        integer                 :: tempi, tempj                     !< temporary index for creating tree shape
        character(100)          :: readstring                       !< read files as text
            
        namelist/NAMTREES/ lapply_treedrag, lreadfile_trees, lapply_sourceSGS, Cd, A_pad 

        if (myid==0) then 
            open(ifnamopt,file=fname_options,status='old',iostat=ierr) ! fname_options='namoptions', iostat=0 if operation is successful, otherwise non-zero value
            read(ifnamopt,NAMTREES,iostat=ierr)
            call checknamelisterror(ierr, ifnamopt, 'NAMTREES') ! print 'problem' + invalid lines in namoptions file
            write(6 ,NAMTREES)
            close(ifnamopt)
        endif

        ! broadcast lapply_treedrag from process with rank 0 to all processes in comm3d communicator, for consistent setting everywhere
        call D_MPI_BCAST(lapply_treedrag,1,0,comm3d,mpierr) ! 1=onevalue, mpierr=checkforerrors
        if (.not. (lapply_treedrag)) return
        call D_MPI_BCAST(lreadfile_trees,1,0,comm3d,mpierr)
        call D_MPI_BCAST(lapply_sourceSGS,1,0,comm3d,mpierr)
        call D_MPI_BCAST(Cd, 1, 0, comm3d, mpierr)
        call D_MPI_BCAST(A_pad, 1, 0, comm3d, mpierr)
        
        if (abs(cu)>1e-15 .or. abs(cv)>1e-15) then
            if(myid==0) print *, 'Problem in namoptions'
            if(myid==0) print *, 'cu or cv cannot be nonzero while using trees'
            if(myid==0) print *, 'The trees would move in that case'
            if(myid==0) print *, 'Set cu and cv to 0. to solve this problem or simulate without trees'
            stop 'ERROR: Problem in namoptions NAMTREES with cu and cv'
        endif
 
        allocate(tree_height(itot+1,jtot+1))                ! +1=extra room to store bc or staggered variables
        allocate(ltree(2-ih:i1+ih,2-jh:j1+jh,k1))      
        allocate(kindex_tree(2-ih:i1+ih,2-jh:j1+jh))

        ! Set default values to zero and false
        tree_height(:,:) = 0
        ltree(:,:,:) = .false.
        kindex_tree(:,:) = 0
        
        !!! 1D input -> 2D tree_height map !!!
        if(myid==0) then
            if (lreadfile_trees) then
                write(6,*) 'Reading inputfile in modtrees'
                open(ifinput, file='trees.inp.'//cexpnr) 
                    do k=1,7 
                        read (ifinput,'(a100)') readstring  ! read in a100 = alphanumeric, otherwise use * to let the code decide its format
                        write (6,*) readstring              ! output the first 100 characters read in standard output format (6)
                    enddo
        
                    !SvdL, 20231218: er wordt hier omgekeerd geloopt, omdat de data-ordening van lezen in fortran anders is dan de standaard ordening van onze wiskundige assen. 
                    do j=jtot+1,2,-1                                ! loop backwarts and keep first value open for bc 
                        do i=2,itot+1    
                            read(ifinput,'(F6.1)') tree_height(i,j) ! F=floating point number, 6.1=6 characters wide with one digit behind decimal point
                        enddo
                    enddo 
                close(ifinput) 
                
                tree_height(1,:)=tree_height(itot+1,:)              !< boundary conditions
                tree_height(:,1)=tree_height(:,jtot+1)
                write(6,*) 'Succesfully read inputfile in modtrees'      
            else 
                write(6,*) 'No trees.inp file found. Stopping modtrees.'
                stop
            endif
        endif !myid==0
        
        call D_MPI_BCAST(tree_height,(itot+1)*(jtot+1),0,comm3d,mpierr)

        do i=2,i1 ! i1=imax+1
            do j=2,j1
                do k=1,kmax !
                    if(zf(k).LE.tree_height(i+myidx*imax,j+myidy*jmax)) then  ! obstacle height is above mid point of vertical grid                        
                        ! Part below created to grow tree shapes from 2D tree_height map to model overhanging parts
                        startIdx = 0 
                        endIdx = 0
                        tempi = i
                        tempj = j
                        do di = startIdx, endIdx
                            do dj = startIdx, endIdx
                                ltree(tempi+di,tempj+dj,k) = .true.
                                write(6,*) 'ltree', ltree(tempi+di,tempj+dj,k), tempi+di+myidx*imax,tempj+di+myidy*jmax, &
                                            tempi, tempj, k, tempi+di, tempj+dj, tree_height(tempi+myidx*imax,tempj+myidy*jmax), zh(k)
                            end do
                        end do
                    
                    endif
                    
                end do  !k
            end do      !j
        end do          !i

        call excjs(ltree,2,i1,2,j1,1,k1,ih,jh)       
        deallocate(tree_height)
        deallocate(kindex_tree)

        return 
    end subroutine inittrees


    subroutine exittrees
        implicit none
    
        if (.not. (lapply_treedrag)) return
        deallocate(ltree)

        return
    end subroutine exittrees
  
    subroutine applytrees
        use modglobal,      only:   kmax, i1, j1, k1, ih, jh, dx, dy, dzh, dzf, rdt, timee     ! rdt=timeintegrationinterval, timee=elapsed time since start
        use modfields,      only:   um, vm, wm, e12m, &   !t-1
                                    u0, v0, w0, e120, &   !t
                                    up, vp, wp, e12p    !tendency of ..m
        use modtreesdata,   only:   lapply_treedrag, lapply_sourceSGS, Cd, A_pad
        use modmpi,         only:   excjs    
        use modprecision,   only:   field_r
    
        ! Declare local variables
        integer :: i, j, k
        real :: treedrag_u, treedrag_v, source_SGS

        if (.not. lapply_treedrag) return
        do i=2,i1
            do j=2,j1
                do k=1,kmax                 
                    if(ltree(i,j,k)) then               
                        !!! Drag on resolved TKE due to trees !!!
                        treedrag_u = 0 
                        treedrag_v = 0
                        call F_treedrag(Cd, A_pad, u0(i-1,j,k), v0(i,j-1,k), u0(i,j,k), v0(i,j,k), treedrag_u, treedrag_v)
                        up(i-1,j,k) = up(i-1,j,k) + treedrag_u/2        ! drag is calculated in the center of the cell, divided back to the faces here 
                        up(i,j,k) = up(i,j,k) + treedrag_u/2  
                        vp(i,j-1,k) = vp(i,j-1,k) + treedrag_v/2      
                        vp(i,j,k) = vp(i,j,k) + treedrag_v/2

                        !!! Source for SGS-TKE !!!
                        if (lapply_sourceSGS) then
                            source_SGS = 0
                            call source_SGS_TKE(Cd, A_pad, u0(i-1, j, k), v0(i, j-1, k), u0(i, j, k), v0(i, j, k), e120(i, j, k), source_SGS)
                            write(6,*) e12p(i, j, k)
                            e12p(i, j, k) = e12p(i, j, k) + source_SGS
                            write(6,*) source_SGS, e12p(i, j, k), e120(i, j, k)
                        endif    
                    endif
                end do
            end do
        end do
       
        call excjs(up,2,i1,2,j1,1,k1,ih,jh)
        call excjs(vp,2,i1,2,j1,1,k1,ih,jh)
        call excjs(e12p,2,i1,2,j1,1,k1,ih,jh)
        ! write(6,* ) 'applytrees succesfull'
        return
    end subroutine applytrees


    subroutine source_SGS_TKE(Cd, A_pad, u1, v1, u2, v2, e120, source_SGS) ! Drag is overestimated, so sourceterm in SGS to compensate (Patton et al. 2015)
        implicit none

        ! input variables
        real, intent(in) :: Cd              ! Drag coefficient for tree
        real, intent(in) :: A_pad           ! Plant area density
        real, intent(in) :: u1,v1,u2,v2     ! velocity component
        real, intent(in) :: e120            ! SGS-TKE (scalar at cell-center)

        ! output variables
        real, intent(out) :: source_SGS     ! set to zero before call in applytrees

        ! Local variables 
        real :: u_mag                       ! magnitude of the velocity vector

        u_mag = sqrt((0.5*(u1+u2))**2 + (0.5*(v1+v2))**2)   ! Magnitude of the velocity vector at centre of gridcell
        source_SGS = (8/3)*Cd * A_pad * u_mag * e120       

    end subroutine source_SGS_TKE

    
    subroutine F_treedrag(Cd, A_pad, u1, v1, u2, v2, treedrag_u, treedrag_v) ! Calculate drag in centre of cell in u and v direction
        implicit none
        
        ! Input variables
        real, intent(in) :: Cd              ! Drag coefficient for tree
        real, intent(in) :: A_pad           ! Plant area density
        real, intent(in) :: u1, v1, u2, v2  ! Velocity components

        ! Output variables1
        real, intent(out) :: treedrag_u, treedrag_v 

        ! Local variables
        real :: u_mag 

        u_mag = sqrt((0.5*(u1+u2))**2 + (0.5*(v1+v2))**2) ! Magnitude of the velocity vector at centre of gridcell

        ! Drag force component calculated at cell centre, later divided to the faces
        treedrag_u = - Cd * A_pad * 0.5 * (u1+u2) * u_mag
        treedrag_v = - Cd * A_pad * 0.5 * (v1+v2) * u_mag
        ! influence vertical direction very small so neglected
    end subroutine F_treedrag

end module modtrees




    

  