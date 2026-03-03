module global_parameters
    implicit none
    save

    !==============================
    ! File paths
    !==============================
    character(len=*), parameter :: input_file = 'datasets/input_catalog.txt'
    character(len=*), parameter :: rescaled_dist = 'results/rescaled_distances.txt'
    character(len=*), parameter :: susc_index_file_name = 'results/susc_index.txt'
    character(len=*), parameter :: output_cat = 'results/output_catalog.txt'
    character(len=*), parameter :: final_results = 'results/final_results.txt'

    !==============================
    ! Catalog parameters
    !==============================
    integer, parameter :: max_events = 2000000
    integer, parameter :: max_thr = 6000
    real(8), parameter :: mcut = 3.0d0

    !==============================
    ! Susceptibility index parameters
    !==============================
    real(8), parameter :: pthmin = 1d-18
    real(8), parameter :: pthmax = 1d15
    real(8), parameter :: bin0   = 0.05d0

    !==============================
    ! Metric parameters (Baiesi–Paczuski)
    !==============================
    real(8), parameter :: df   = 1.6d0     ! fractal dimension
    real(8), parameter :: bval = 1.0d0     ! Gutenberg–Richter b-value

    !==============================
    ! Smoothing parameters
    !==============================
    integer, parameter :: max_points = 10000

end module global_parameters
