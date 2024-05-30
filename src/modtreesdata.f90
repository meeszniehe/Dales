!< modtrees_testcase_data.f90 file

module modtreesdata
    implicit none
    save
  
    ! Global settings
    logical :: lapply_treedrag = .false.    !< Switch to enable tree method
    logical :: lreadfile_trees = .false.    !< Switch to read tree height data from a file
    logical :: lapply_sourceSGS = .false.   !< Switch to apply source SGS model
    real        :: Cd        = 0.264        !< Drag coefficient tree (taken from modcanopy?)
    real        :: A_pad        = 1         !< Cross-sectional area

    logical, allocatable    :: ltree(:,:,:)                !< true/false array to indicate stem cells
end module modtreesdata

