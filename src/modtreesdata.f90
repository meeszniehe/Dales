!< modtrees_testcase_data.f90 file

module modtreesdata
    implicit none
    save
  
    ! Global settings
    logical :: lapply_trees = .true.       !< Switch to enable tree method
    logical :: lreadfile_trees = .true.   !< Switch to read tree height data from a file
    
  
    !< Number of grid points in a slab excluding obstacles, and the number of obstacle points
    integer, allocatable :: Nair_trees(:) !SvdL, 20231218: veranderd in Nair_trees omdat naam Nair in conflict zou kunnen komen met Nair (modibm). Uiteindelijk moet afhankelijk van gebruik modules Nair = som(Nai_apart) genomen worden.
    
    real    :: C_stem        = 0.15            !< Drag coefficient for stem, based on modcanopy
    real    :: A_stem        = 1            !< Cross-sectional area of the stem

    ! real    :: C_leaves        = ???            !< Drag coefficient for leaves
    ! real    :: A_leaves        = ???            !< Cross-sectional area of the leaves, value between 0.1-2.0 m2 m−3 depending on season (grylls 2021)

    logical, allocatable    :: ltree_stem(:,:,:)                !< true/false array to indicate stem cells, !SvdL, 20231218: als het een true/false array moet zijn, moet je hem als logical declareren
    !logical, allocatable    :: ltree_leaves(:,:,:)
    
end module modtreesdata

