!> \file modtrees.f90
!! By Meesz Niehe, email: meesz@niehe.com, TU Delft, section Atmospheric Physics, date 

! First feedback Steven
! Github test
module modtrees
    ! Using tree data module
    use modtreesdata, only : lapply_trees, lreadfile_trees, ltree_stem, C_stem, A_stem !, ltree_leaves
    use modprecision        
    implicit none
    save
    public :: inittrees, exittrees, applytrees 

    contains
        subroutine inittrees
        !!! INITIALIZATION !!!
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
        
        !< Number of grid points in a slab excluding obstacles, and the number of obstacle points
        !integer, allocatable    :: Nairl_trees(:)
        !integer, allocatable    :: Nair_trees(:) !SvdL, 20231218: veranderd in Nair_trees omdat naam Nair in conflict zou kunnen komen met Nair (modibm). Uiteindelijk moet afhankelijk van gebruik modules Nair = som(Nai_apart) genomen worden.
        !< for calculating free spaces of air 

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
        
        !!! MOVING TREES !!!
        if (abs(cu)>1e-15 .or. abs(cv)>1e-15) then
            if(myid==0) print *, 'Problem in namoptions' ! only print when running the mainprocess myid=0
            if(myid==0) print *, 'cu or cv cannot be nonzero while using trees'
            if(myid==0) print *, 'The trees would move in that case'
            if(myid==0) print *, 'Set cu and cv to 0. to solve this problem or simulate without trees'
            stop 'ERROR: Problem in namoptions NAMTREES with cu and cv'
        endif

        !!! ALLOCATE VARIABLES !!!  
        !write(6,*) 'allocating fields in modtrees'
        allocate(tree_height(itot+1,jtot+1))                ! +1=extra room to store bc or staggered variables, velocity on face/pressure in center?
        allocate(ltree_stem(2-ih:i1+ih,2-jh:j1+jh,k1))      ! 'true' means there is a stem
        allocate(kindex_stem(2-ih:i1+ih,2-jh:j1+jh))
        !allocate(ltree_leaves(2-ih:i1+ih,2-jh:j1+jh,k1))   ! 'true' means leaves

        !allocate(Nair_trees(k1))                                  ! k1=kmax+1
        !allocate(Nairl_trees(k1))

        ! Set standard value to zero
        tree_height(:,:) = 0 ! choose value higher than ground and lower than lowest tree
 
        ltree_stem(:,:,:) = .false.
        kindex_stem(:,:) = 0
        !ltree_leaves (:,:,:) = .false.
        
        !!! 1D input -> 2D Array !!!
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
        !write(6,*) 'broadcasted tree_height'

        !!! INDICATE STEM & LEAF CELLS !!!
        do i=2,i1
            do j=2,j1
                do k=1,kmax
                    if(zf(k).LE.tree_height(i+myidx*imax,j+myidy*jmax)) then  ! obstacle height is above mid point of vertical grid
                        ltree_stem(i,j,k) = .true.         ! true/false array to indicate stem cells
                        kindex_stem(i,j)   = k + 1     	! werkt niet voor overhangende bladeren/takken
                        ! 
                        write(6,*) 'ltree_stem',i+myidx*imax,j+myidy*jmax,i,j,k,ltree_stem(i,j,k),tree_height(i+myidx*imax,j+myidy*jmax),zh(kindex_stem(i,j))
                        !write(6,*) 'indicated stem and leaf cells'
                    endif
                end do  !k
            end do      !j
        end do          !i

        ! from ltree_stem data exchange in x,y,z dimension, 2->i1, 2->j1, 1->k1, 
        call excjs(ltree_stem,2,i1,2,j1,1,k1,ih,jh)       
        !write(6,*) 'excjs ltree_stem succesfull'
        
        !!! COUNT NUMBER OF GRIDCELLS WITHOUT TREE PER VERTICAL LEVEL !!!
        ! Voor statistiek gebruikt, niet voor simuleren
        !Nair_trees(:) = 0 
        !Nairl_trees(:) = 0    ! local counter
        !do i=2,i1
        !    do j=2,j1
        !        do k=kindex_stem(i,j),k1          
        !            Nairl_trees(k) = Nairl_trees(k)+1 
        !        enddo
        !    enddo
        !enddo
        !write(6,*) 'Statistics done' 
         
        ! In d_mpi_allreduce k1 is the count which differs per i,j position due to a changing kindex_stem
        !call D_MPI_ALLREDUCE(Nairl_trees,Nair_trees,k1,MPI_SUM,comm3d,mpierr)
        !write(6,*) 'MPI_ALLREDUCE succesfull'

        !write(6,* ) 'start deallocate'
        deallocate(tree_height)
        !write(6,* ) 'deallocated tree_height'
        !deallocate(Nairl_trees)
        !write(6,* ) 'deallocated Nairl_trees'
        !deallocate(Nair_trees)
        !write(6,* ) 'deallocated Nair_trees'
        deallocate(kindex_stem)
        !write(6,* ) 'deallocated kindex_stem'
        !write(6,*) 'exit initrees'

        return ! why? !SvdL, 20231218: ik weet het niet zeker, de fortran beschrijving op internet is er ook niet heel duidelijk over. In feite sluit je hiermee de subroutine af en geef je controle terug aan de routine erboven, maar het statement END SUBROUTINE zou in principe hetzelfde al moeten doen. Dus het lijkt me dubbelop. 
    end subroutine inittrees


    subroutine exittrees
        implicit none
    
        if (.not. (lapply_trees)) return
        deallocate(ltree_stem)
        !write(6,* ) 'deallocated ltree_stem'
        !deallocate(kindex_stem) ! deze is niet gealloceerd in hoofd routine, dus kan niet hier gedeallocate worden, maar onderaan inittrees?
        !deallocate(ltree_leaves)

        return
    end subroutine exittrees
  
    subroutine applytrees
        ! no use of: thl=liquid potential temperature, qt=total specific humidity, e12=sqrt(tke), sv=scalar variance
        use modglobal,      only:   kmax, i1, j1, k1, ih, jh, dx, dy, dzh, dzf, rdt, timee     ! rdt=timeintegrationinterval, timee=elapsed time since start
        use modfields,      only:   um, vm, wm, &   !t-1
                                    u0, v0, w0, &   !t
                                    up, vp, wp   !tendency of ..m
        use modtreesdata,   only:   lapply_trees, C_stem, A_stem
        use modmpi,         only:   excjs    
        use modprecision,   only:   field_r
    
        ! Declare local variables
        integer :: i, j, k
        real :: drag_stem_u, drag_stem_v

        !real(field_r) :: ...   

        !write(6,* ) 'starting with applytrees after module and variable declarations'
        if (.not. lapply_trees) return
    
        !!! IMPLEMENTATION DRAG FORCE !!!
        ! Runge Kutta coefficients? (rk3step, rk3coef, rk3coefi)
    
        !SvdL, 20231218: ik snap je berekening/formule hier niet helemaal (gebruik rdt, etc). Morgen bespreken. Code technisch werkt het waarschijnlijk wel.

        do i=2,i1
            do j=2,j1
                do k=2,kmax                  
                    if(ltree_stem(i,j,k)) then   ! could be faster by limiting k to highest tree value?
                        ! Calculate drag in centre of cell in u and v direction
                        call drag_force_stem(C_stem, A_stem, um(i-1,j,k), vm(i,j-1,k), wm(i,j,k-1), um(i,j,k), vm(i,j,k), wm(i,j,k), drag_stem_u, drag_stem_v)
                        ! Reassign the velocity value at the faces adjusted for drag in u and v direction
                        up(i-1,j,k) = up(i-1,j,k) - drag_stem_u/2        ! averaged for gridspacing dx, up = in kracht uitgedrukt, dx weg?
                        up(i,j,k) = up(i,j,k) - drag_stem_u/2  
                        vp(i,j-1,k) = vp(i,j-1,k) - drag_stem_v/2        ! averaged for gridspacing dy
                        vp(i,j,k) = vp(i,j,k) - drag_stem_v/2
                        wp(i,j,k-1) = 0
                        wp(i,j,k) = 0
                        write(6,*) 'tendencies',up(i-1,j,k), up(i,j,k), vp(i,j-1,k), vp(i,j,k), wp(i,j,k-1), wp(i,j,k)
                    !elseif (ltree_leaves(i,j,k)) then   ! Drag force due to leaves
                        !call drag_force_leaves(C_leaves, A_leaves, um(i,j,k), vm(i,j,k), wm(i,j,k), drag_leaves_u, drag_leaves_v, drag_leaves_w)
                        !up(i,j,k) = up(i,j,k) - (drag_leaves_u*rdt)/(dx*rho_air)         ! averaged for gridspacing dx
                        !vp(i,j,k) = vp(i,j,k) - (drag_leaves_v*rdt)/(dy*rho_air)         ! averaged for gridspacing dy
                        !wp(i,j,k) = wp(i,j,k) - (drag_leaves_w*rdt)/(dzf*rho_air)
                    endif
                end do
            end do
        end do
        !write(6,* ) 'dragforce applied, move to excjs variables'
        ! E.g., synchronizing data across processors, call excjs, also for e12,thl,qt,svp possibly?
        call excjs(up,2,i1,2,j1,1,k1,ih,jh)
        call excjs(vp,2,i1,2,j1,1,k1,ih,jh)
        call excjs(wp,2,i1,2,j1,1,k1,ih,jh)
        write(6,* ) 'applytrees succesfull'
        return
    end subroutine applytrees

    !SvdL, 20231218: 3 dingen, (1) ik was verrast dat dit zo mocht met intent(out) drag_stem_u en niet declareren in enige outer scope? In C++ zou dat nooit mogen, dus niet geheel zeker over Fortran
    ! en (2) nu bereken je per i,j,k iteratie de kracht, maar dat houdt ook N^3 functiecalls in. Dat gebeurt in modibm ook, al vraag ik me af of dat inderdaad de beste optie is.. voor nu gewoon later zou ik zeggen.
    ! en (3) de snelheden zijn gegeven op de wanden van de cellen, dus de snelheid in het midden van de zelf is een middeling van deze snelheden. Dat zou eventueel nog geimplemteerd moeten worden (zie onder) 
    !SvdL, 20231218: heb de indentatie deels aangepast.
    subroutine drag_force_stem(C_stem, A_stem, u1, v1, w1, u2, v2, w2, drag_stem_u, drag_stem_v)
        implicit none
        
        ! Input variables
        real, intent(in) :: C_stem   ! Drag coefficient for stem
        real, intent(in) :: A_stem   ! Cross-sectional area of the stem
        real, intent(in) :: u1, v1, w1, u2, v2, w2  ! Velocity components

        ! Output variables
        real, intent(out) :: drag_stem_u, drag_stem_v 

        ! Local variables
        real :: u_mag   ! Magnitude of the velocity vector

        ! Magnitude of the velocity vector at centre of gridcell
        u_mag = 0.25*sqrt((u1+u2)**2 + (v1+v2)**2 + (w1+w2)**2)

        ! Calculate the drag force components
        drag_stem_u = C_stem * A_stem * 0.5 * (u1+u2) * u_mag
        drag_stem_v = C_stem * A_stem * 0.5 * (v1+v2) * u_mag

        write(6,*) 'drag_force', u_mag, drag_stem_u, drag_stem_v, u1, u2, v1, v2, w1, w2
        !SvdL, 20231218: deze heb ik uitgecommend: vanaf bovenaf gekeken is A_stem niet relevant, maar waarschijnlijk een veel kleiner oppervlak. Ook zal w zelf erg klein zijn.
        ! drag_stem_w = -C_stem * A_stem * w * u_mag
        !write(6,* ) 'dragforce calculated'
     end subroutine drag_force_stem
    
    
end module modtrees




    

  