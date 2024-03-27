!> \file modtrees.f90
!! By Meesz Niehe, email: meesz@niehe.com, TU Delft, section Atmospheric Physics, date 

module modtrees
    use modtreesdata, only : lapply_trees, lreadfile_trees, ltree_stem, C_stem, A_stem !, ltree_leaves
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

        real(field_r), allocatable :: tree_height(:,:)              !< 2D array to store tree heights at each grid point (x,y), field_r precision?
        integer                 :: i, j, k, ierr, ii, jj, kk, n     !< initialize integer to loop over
        integer, allocatable    :: kindex_stem(:,:)                      !< index of stem height
        character(100)          :: readstring                       !< read files as text
            
        namelist/NAMTREES/ lapply_trees, lreadfile_trees, C_stem, A_stem 

        if(myid==0) then 
            open(ifnamopt,file=fname_options,status='old',iostat=ierr) ! fname_options='namoptions', iostat=0 if operation is successful, otherwise non-zero value
            read(ifnamopt,NAMTREES,iostat=ierr)
            call checknamelisterror(ierr, ifnamopt, 'NAMTREES') ! print 'problem' + invalid lines in namoptions file
            write(6 ,NAMTREES)
            close(ifnamopt)
        endif

        !!! BROADCASTING !!!
        ! broadcast lapply_trees from process with rank 0 to all processes in comm3d communicator, for consistent setting everywhere
        call D_MPI_BCAST(lapply_trees,1,0,comm3d,mpierr) ! 1=onevalue, mpierr=checkforerrors
        if (.not. (lapply_trees)) return
        call D_MPI_BCAST(lreadfile_trees,1,0,comm3d,mpierr)
        call D_MPI_BCAST(C_stem, 1, 0, comm3d, mpierr)
        call D_MPI_BCAST(A_stem, 1, 0, comm3d, mpierr)
        
        !!! MOVING TREES !!!
        if (abs(cu)>1e-15 .or. abs(cv)>1e-15) then
            if(myid==0) print *, 'Problem in namoptions' ! only print when running the mainprocess myid=0
            if(myid==0) print *, 'cu or cv cannot be nonzero while using trees'
            if(myid==0) print *, 'The trees would move in that case'
            if(myid==0) print *, 'Set cu and cv to 0. to solve this problem or simulate without trees'
            stop 'ERROR: Problem in namoptions NAMTREES with cu and cv'
        endif

        !!! ALLOCATE VARIABLES !!!  
        allocate(tree_height(itot+1,jtot+1))                ! +1=extra room to store bc or staggered variables, velocity on face/pressure in center?
        allocate(ltree_stem(2-ih:i1+ih,2-jh:j1+jh,k1))      ! 'true' means there is a stem
        allocate(kindex_stem(2-ih:i1+ih,2-jh:j1+jh))
        !allocate(ltree_leaves(2-ih:i1+ih,2-jh:j1+jh,k1))   ! 'true' means leaves

        ! Set default values to zero and false
        tree_height(:,:) = 0
        ltree_stem(:,:,:) = .false.
        kindex_stem(:,:) = 0
        !ltree_leaves (:,:,:) = .false.
        
        !!! 1D input -> 2D tree_height map !!!
        if(myid==0) then
            if (lreadfile_trees) then
                write(6,*) 'Reading inputfile in modtrees'
                open(ifinput, file='trees.inp.'//cexpnr) ! is cexpnr working?
                    do k=1,7 
                        read (ifinput,'(a100)') readstring ! read in a100 = alphanumeric, otherwise use * to let the code decide its format
                        write (6,*) readstring ! output the first 100 characters read in standard output format (6)
                    enddo
        
                    !SvdL, 20231218: er wordt hier omgekeerd geloopt, omdat de data-ordening van lezen in fortran anders is dan de standaard ordening van onze wiskundige assen. 
                    !SvdL, 20231218: moet je ook rekening mee houden met het maken van de case: waardes in xy-vlak omgekeerd wegschrijven (zie voorbeeldscript, moet ik nog sturen..)
                    do j=jtot+1,2,-1        ! loop backwarts? keep first value open for bc 
                        do i=2,itot+1    
                            read(ifinput,'(F6.1)') tree_height(i,j) ! F=floating point number, 6.1=6 characters wide with one digit behind decimal point
                        enddo
                    enddo 
                close(ifinput) 
                
                tree_height(1,:)=tree_height(itot+1,:) !< boundary conditinos
                tree_height(:,1)=tree_height(:,jtot+1)
                
                write(6,*) 'Succesfully read inputfile in modtrees'      
            else 
                write(6,*) 'No trees.inp file found. Stopping modtrees.'
                stop
            endif
        endif !myid==0
        
        call D_MPI_BCAST(tree_height,(itot+1)*(jtot+1),0,comm3d,mpierr)

        !!! INDICATE STEM & LEAF CELLS !!!
        do i=2,i1 ! i1=imax+1
            do j=2,j1
                do k=1,kmax !
                    if(zf(k).LE.tree_height(i+myidx*imax,j+myidy*jmax)) then  ! obstacle height is above mid point of vertical grid
                        !!!! Tree is represented by one straight line going straight up !!!!
                        ltree_stem(i,j,k) = .true.     ! true/false array to indicate stem cells
                        kindex_stem(i,j)   = k + 1     ! werkt niet voor overhangende bladeren/takken
                        write(6,*) 'ltree_stem',i+myidx*imax,j+myidy*jmax,i,j,k,ltree_stem(i,j,k),tree_height(i+myidx*imax,j+myidy*jmax),zh(kindex_stem(i,j))
                        
                        ! !!!! Divide tree in smaller stem and thicker crown, same C_d value !!!!
                        ! if (k <= 2) then ! For the lowest two levels, thickness is 1 grid cell
                        !     startIdx = 0
                        !     endIdx = 0
                        ! else ! For higher levels, thickness is 3 grid cells
                        !     startIdx = -1
                        !     endIdx = 1
                        ! endif 
                        ! do di = startIdx, endIdx
                        !     do dj = startIdx, endIdx
                        !         ! Check boundaries to prevent out-of-bound errors
                        !         if (i+di >= 2 .AND. i+di <= imax .AND. j+dj >= 2 .AND. j+dj <= jmax) then ! for i=2 en i+di=1, moet die 1 niet de laatste cel van het naastgelegen gebied zijn?
                        !             ltree_stem(i+di,j+dj,k) = .true.
                        !             ! Only update kindex_stem for the central cell
                        !             if (di == 0 .AND. dj == 0) then
                        !                 kindex_stem(i,j) = k + 1
                        !             endif
                        !         endif
                        !     end do
                        ! end do

                    endif
                end do  !k
            end do      !j
        end do          !i

        ! from ltree_stem data exchange in x,y,z dimension, 2->i1, 2->j1, 1->k1, 
        call excjs(ltree_stem,2,i1,2,j1,1,k1,ih,jh)       
        deallocate(tree_height)
        deallocate(kindex_stem)

        return ! why? !SvdL, 20231218: ik weet het niet zeker, de fortran beschrijving op internet is er ook niet heel duidelijk over. In feite sluit je hiermee de subroutine af en geef je controle terug aan de routine erboven, maar het statement END SUBROUTINE zou in principe hetzelfde al moeten doen. Dus het lijkt me dubbelop. 
    end subroutine inittrees


    subroutine exittrees
        implicit none
    
        if (.not. (lapply_trees)) return
        deallocate(ltree_stem)
        !deallocate(ltree_leaves)

        return
    end subroutine exittrees
  
    subroutine applytrees
        ! no use of: thl=liquid potential temperature, qt=total specific humidity, e12=sqrt(tke), sv=scalar variance
        use modglobal,      only:   kmax, i1, j1, k1, ih, jh, dx, dy, dzh, dzf, rdt, timee     ! rdt=timeintegrationinterval, timee=elapsed time since start
        use modfields,      only:   um, vm, wm, e12m&   !t-1
                                    u0, v0, w0, e120&   !t
                                    up, vp, wp, e12p    !tendency of ..m
        use modtreesdata,   only:   lapply_trees, C_stem, A_stem
        use modmpi,         only:   excjs    
        use modprecision,   only:   field_r
    
        ! Declare local variables
        integer :: i, j, k
        real :: drag_stem_u, drag_stem_v 

        if (.not. lapply_trees) return
    
        !!! IMPLEMENTATION DRAG FORCE !!!
        do i=2,i1
            do j=2,j1
                do k=1,kmax                 
                    if(ltree_stem(i,j,k)) then   ! could be faster by limiting k to highest tree value?
                        write(6,*) 'ltree is true for index', i, j, k                        
                        ! Drag on resolved TKE
                        drag_stem_u = 0 
                        drag_stem_v = 0
                        call drag_force_stem(C_stem, A_stem, u0(i-1,j,k), v0(i,j-1,k), u0(i,j,k), v0(i,j,k), drag_stem_u, drag_stem_v)
                        ! Reassign the velocity value at the faces adjusted for drag in u and v direction
                        up(i-1,j,k) = up(i-1,j,k) - drag_stem_u/2        ! both sides get half the drag calculated from the middle. 
                        up(i,j,k) = up(i,j,k) - drag_stem_u/2  
                        vp(i,j-1,k) = vp(i,j-1,k) - drag_stem_v/2      
                        vp(i,j,k) = vp(i,j,k) - drag_stem_v/2
                        write(6,*) 'resolved drag force applied'
                        ! Drag on SFS-TKE
                        drag_SFS = 0
                        call drag_force_SFS_TKE(C_stem, A_stem, u0(i-1,j,k), v0(i,j-1,k), u0(i,j,k), v0(i,j,k), e120, drag_SFS) ! e120?? 
                        write(6,*) 'e12p beforehand: ', e12p(i,j,k)
                        e12p(i,j,k) = e12p(i,j,k) - drag_SFS
                        write(6,*) 'SFS drag force applied, drag: ', drag_SFS, 'and e12p afterwards: ', e12p(i,j,k)
                        !wp(i,j,k-1) = 0
                        !wp(i,j,k) = 0
                    !elseif (ltree_leaves(i,j,k)) then   ! Drag force due to leaves
                        !call drag_force_leaves(C_leaves, A_leaves, um(i,j,k), vm(i,j,k), wm(i,j,k), drag_leaves_u, drag_leaves_v, drag_leaves_w)
                        !up(i,j,k) = up(i,j,k) - (drag_leaves_u*rdt)/(dx*rho_air)         ! averaged for gridspacing dx
                        !vp(i,j,k) = vp(i,j,k) - (drag_leaves_v*rdt)/(dy*rho_air)         ! averaged for gridspacing dy
                        !wp(i,j,k) = wp(i,j,k) - (drag_leaves_w*rdt)/(dzf*rho_air)
                    endif
                end do
            end do
        end do
        ! E.g., synchronizing data across processors, call excjs, also for e12,thl,qt,svp possibly?
        call excjs(up,2,i1,2,j1,1,k1,ih,jh)
        call excjs(vp,2,i1,2,j1,1,k1,ih,jh)
        !call excjs(wp,2,i1,2,j1,1,k1,ih,jh)
        write(6,* ) 'applytrees succesfull'
        return
    end subroutine applytrees


    subroutine drag_force_SFS_TKE(C_stem, A_stem, u1, v1, u2, v2, e120, drag_SFS) 
        ! Drag force on SFS-TKE, take into account direction of surrounding windspeeds as SFS-TKE is a scalar at the cell-center
        implicit none

        ! input variables
        real, intent(in) :: C_stem   ! Drag coefficient for stem
        real, intent(in) :: A_stem   ! Area of stem
        real, intent(in) :: u1,v1,u2,v2 ! velocity component
        
        ! output variables
        real, intent(out) :: drag_SFS ! set to zero before call in applytrees

        ! Local variables 
        real :: u_mag ! magnitude of the velocity vector

        ! Magnitude of the velocity vector at centre of gridcell
        u_mag = 0.25*sqrt((u1+u2)**2 + (v1+v2)**2) !+ (w1+w2)**2)
        
        ! work performed by SFS motions against canopy drag (Patton et al. 2015)
        drag_SFS = (8/3)*C_stem * A_stem * u_mag * e120

    end subroutine drag_force_SFS_TKE


    !SvdL, 20231218: 3 dingen, (1) ik was verrast dat dit zo mocht met intent(out) drag_stem_u en niet declareren in enige outer scope? In C++ zou dat nooit mogen, dus niet geheel zeker over Fortran
    ! en (2) nu bereken je per i,j,k iteratie de kracht, maar dat houdt ook N^3 functiecalls in. Dat gebeurt in modibm ook, al vraag ik me af of dat inderdaad de beste optie is.. voor nu gewoon later zou ik zeggen.
    ! en (3) de snelheden zijn gegeven op de wanden van de cellen, dus de snelheid in het midden van de zelf is een middeling van deze snelheden. Dat zou eventueel nog geimplemteerd moeten worden (zie onder) 
    !SvdL, 20231218: heb de indentatie deels aangepast.
    subroutine drag_force_stem(C_stem, A_stem, u1, v1, u2, v2, drag_stem_u, drag_stem_v) ! Calculate drag in centre of cell in u and v direction
        implicit none
        
        ! Input variables
        real, intent(in) :: C_stem   ! Drag coefficient for stem
        real, intent(in) :: A_stem   ! Cross-sectional area of the stem
        real, intent(in) :: u1, v1, u2, v2  ! Velocity components

        ! Output variables
        real, intent(out) :: drag_stem_u, drag_stem_v 

        ! Local variables
        real :: u_mag   ! Magnitude of the velocity vector

        ! Magnitude of the velocity vector at centre of gridcell
        u_mag = 0.25*sqrt((u1+u2)**2 + (v1+v2)**2) !+ (w1+w2)**2)

        ! Calculate the drag force components
        drag_stem_u = C_stem * A_stem * 0.5 * (u1+u2) * u_mag ! 0.5(u1+u2) -> avg in cell center
        drag_stem_v = C_stem * A_stem * 0.5 * (v1+v2) * u_mag

        !SvdL, 20231218: deze heb ik uitgecommend: vanaf bovenaf gekeken is A_stem niet relevant, maar waarschijnlijk een veel kleiner oppervlak. Ook zal w zelf erg klein zijn.
        !drag_stem_w = -C_stem * A_stem * w * u_mag
    end subroutine drag_force_stem

end module modtrees




    

  