program susceptibility_index
  use iso_fortran_env, only: real64, int64
  use ieee_arithmetic, only: ieee_is_finite
  use global_parameters
  implicit none

  integer, parameter :: dp = real64, derivative_step = 10
  real(dp), parameter :: year_seconds = 365.25_dp*86400.0_dp
  real(dp), parameter :: radians = acos(-1.0_dp)/180.0_dp
  real(dp), parameter :: smoothing_fraction = 0.10_dp
  real(dp), parameter :: tail_fraction = 0.02_dp
  real(dp), parameter :: shift_fraction = 0.30_dp

  ! True: shift used in both supplied scripts.
  ! False: linear shift stated in the paper.
  logical, parameter :: logarithmic_shift = .true.

  real(dp), allocatable :: time(:), lat(:), lon(:), mag(:)
  real(dp), allocatable :: phi(:), lambda(:), cos_phi(:), weight(:)
  integer, allocatable :: parent_bin(:)

  integer(int64), save :: links(-max_thr:max_thr)
  integer(int64) :: delta_links
  integer :: nn(-max_thr-1:max_thr), background(-max_thr:max_thr)
  real(dp), save :: si(max_points), threshold(max_points)
  integer :: k(max_points), bins(max_points)

  integer :: n, np, ipmin, ipmax, ip, u, i, flag, nbg
  integer :: best, valley, left_peak, right_peak

  call read_catalog()

  print *, 'Retained events:', n
  print *, 'Metric b-value:', bval, '  df:', df

  if (bin0 <= 0 .or. pthmin <= 0 .or. pthmax <= pthmin) &
    error stop 'Invalid threshold settings'

  ipmin = ceiling(log(pthmin)/bin0)
  ipmax = floor(log(pthmax)/bin0)

  if (ipmin < -max_thr .or. ipmax > max_thr) &
    error stop 'Increase max_thr'
  if (ipmax-ipmin <= derivative_step) &
    error stop 'Threshold range too narrow'

  call count_links()

  ! Output: K, SI, NN, similarity threshold, bin index, link count.
  open(newunit=u, file=susc_index_file_name, status='replace')

  np = 0
  do ip = ipmax, ipmin+derivative_step, -1
    if (background(ip) <= 1) exit

    delta_links = links(ip-derivative_step)-links(ip)
    if (delta_links <= 0) cycle

    if (np == max_points) error stop 'Increase max_points'
    np = np+1

    k(np) = background(ip)
    bins(np) = ip
    threshold(np) = exp(ip*bin0)

    si(np) = real(background(ip)-background(ip-derivative_step), dp) &
           / real(delta_links, dp)*real(links(ip), dp)/real(n, dp)

    write(u, *) k(np), si(np), nn(ip), threshold(np), ip, links(ip)
  end do
  close(u)

  call select_threshold()
  ip = bins(best)

  ! Output: time, latitude, longitude, magnitude, classification.
  ! Classification: 0 = background, 1 = triggered.
  open(newunit=u, file=output_cat, status='replace')

  nbg = 0
  do i = 1, n
    flag = 0
    if (i > 1) then
      if (parent_bin(i) >= ip) flag = 1
    end if

    if (flag == 0) nbg = nbg+1
    write(u, *) time(i), lat(i), lon(i), mag(i), flag
  end do
  close(u)

  ! Preserve the repository's ten summary columns.
  open(newunit=u, file=final_results, status='replace')
  write(u, *) k(best), si(best), threshold(best), &
              k(valley), si(valley), threshold(valley), &
              k(left_peak), threshold(left_peak), &
              k(right_peak), threshold(right_peak)
  close(u)

  print *, 'Background events:', nbg
  print *, 'Similarity threshold:', threshold(best)

contains

  subroutine read_catalog()
    integer :: unit, ios, pass, count
    real(dp) :: row(4)

    open(newunit=unit, file=input_file, status='old', action='read')

    ! Count first, then allocate only as much memory as needed.
    do pass = 1, 2
      count = 0

      do
        read(unit, *, iostat=ios) row
        if (ios < 0) exit
        if (ios > 0) &
          error stop 'Expected: time latitude longitude magnitude'

        if (.not. all(ieee_is_finite(row))) &
          error stop 'Non-finite input value'

        if (row(4) < mcut) cycle

        count = count+1
        if (count > max_events) error stop 'Increase max_events'
        if (pass == 1) cycle

        time(count) = row(1)
        lat(count) = row(2)
        lon(count) = row(3)
        mag(count) = row(4)

        if (abs(lat(count)) > 90) error stop 'Invalid latitude'

        if (count > 1) then
          if (time(count) < time(count-1)) &
            error stop 'Sort catalog by time'
        end if
      end do

      if (pass == 1) then
        n = count
        if (n < 3) error stop 'Too few retained events for SI'
        allocate(time(n), lat(n), lon(n), mag(n), parent_bin(n))
        rewind(unit)
      end if
    end do

    close(unit)
  end subroutine read_catalog


  subroutine count_links()
    integer :: event, parent, b, bmax, unit
    real(dp) :: a, distance, dt, tr, rr, eta
    real(dp) :: best_eta, best_t, best_r

    ! Precompute quantities reused in many pair comparisons.
    allocate(phi(n), lambda(n), cos_phi(n), weight(n))
    phi = lat*radians
    lambda = lon*radians
    cos_phi = cos(phi)
    weight = 10.0_dp**(-0.5_dp*bval*mag)

    links = 0_int64
    nn = 0
    parent_bin = ipmin-1

    open(newunit=unit, file=rescaled_dist, status='replace')

    do event = 2, n
      best_eta = huge(1.0_dp)
      best_t = 0
      best_r = 0
      bmax = ipmin-1

      do parent = event-1, 1, -1
        ! Haversine epicentral distance in kilometres.
        a = sin(0.5_dp*(phi(event)-phi(parent)))**2 &
          + cos_phi(event)*cos_phi(parent) &
          * sin(0.5_dp*(lambda(event)-lambda(parent)))**2

        distance = 2.0_dp*6370.0_dp &
                 * asin(sqrt(min(1.0_dp, max(0.0_dp, a))))

        dt = time(event)-time(parent)

        ! Regularize equal times and identical locations independently.
        if (dt <= 0) dt = 1.0_dp
        if (distance <= 0) distance = 0.01_dp

        ! T is in years; R is in km**df; eta = T*R.
        tr = (dt/year_seconds)*weight(parent)
        rr = distance**df*weight(parent)
        eta = tr*rr

        ! Similarity = 1/eta.
        ! FLOOR is necessary when the logarithm is negative.
        b = floor(-log(eta)/bin0)
        b = max(ipmin-1, min(ipmax, b))

        if (b >= ipmin) links(b) = links(b)+1_int64
        bmax = max(bmax, b)

        ! Keep the actual nearest neighbour, including within-bin ties.
        if (eta < best_eta) then
          best_eta = eta
          best_t = tr
          best_r = rr
        end if
      end do

      parent_bin(event) = bmax
      nn(bmax) = nn(bmax)+1
      write(unit, *) best_t, best_r
    end do

    close(unit)
    deallocate(phi, lambda, cos_phi, weight)

    ! Accumulate numbers of pairs above each similarity threshold.
    do b = ipmax-1, ipmin, -1
      links(b) = links(b)+links(b+1)
    end do

    ! The first event has no parent and is always background.
    background(ipmin) = 1+nn(ipmin-1)
    do b = ipmin+1, ipmax
      background(b) = background(b-1)+nn(b-1)
    end do

    if (links(ipmax) > 0 .or. background(ipmin) > 1) &
      print *, 'Note: consider widening pthmin/pthmax to cover all links.'
  end subroutine count_links


  subroutine select_threshold()
    real(dp) :: prefix(0:np), smoothed(np), target
    integer :: maxima(np), minima(np)
    integer :: width, half, limit, lo, hi, j, nmax, nmin
    logical :: found

    if (np < 5) &
      error stop 'Too few SI points; inspect the SI output'

    ! Prefix sums make each moving-average evaluation inexpensive.
    prefix(0) = 0
    do j = 1, np
      prefix(j) = prefix(j-1)+si(j)
    end do

    limit = min(np-2, max(3, int(smoothing_fraction*np)))
    found = .false.

    do width = 3, limit, 2
      half = width/2
      lo = half+1
      hi = np-half

      do j = lo, hi
        smoothed(j) = max(0.0_dp, &
          (prefix(j+half)-prefix(j-half-1))/real(width, dp))
      end do

      nmax = 0
      nmin = 0

      ! Both neighbours must have computed smoothing values.
      do j = lo+1, hi-1
        if (k(j) <= tail_fraction*n .or. &
            k(j) >= (1-tail_fraction)*n) cycle

        if (smoothed(j) > smoothed(j-1) .and. &
            smoothed(j) > smoothed(j+1)) then
          nmax = nmax+1
          maxima(nmax) = j
        end if

        if (smoothed(j) < smoothed(j-1) .and. &
            smoothed(j) < smoothed(j+1)) then
          nmin = nmin+1
          minima(nmin) = j
        end if
      end do

      if (nmax == 2 .and. nmin == 1) then
        if (maxima(1) < minima(1) .and. &
            minima(1) < maxima(2)) then
          valley = minima(1)
          right_peak = maxima(1)
          left_peak = maxima(2)
          found = .true.
          exit
        end if
      end if
    end do

    if (.not. found) then
      ! Fallback: lowest detected interior local minimum.
      if (nmin == 0) &
        error stop 'No interior SI minimum; inspect the curve'

      valley = minima(1)
      do j = 2, nmin
        if (smoothed(minima(j)) < smoothed(valley)) &
          valley = minima(j)
      end do

      right_peak = lo-1+maxloc(smoothed(lo:valley-1), dim=1)
      left_peak = valley+maxloc(smoothed(valley+1:hi), dim=1)

      if (smoothed(right_peak) <= smoothed(valley) .or. &
          smoothed(left_peak) <= smoothed(valley)) &
        error stop 'No bracketed SI minimum; inspect the curve'

      print *, 'Using the interior-minimum fallback.'
    end if

    if (logarithmic_shift) then
      if (smoothed(valley) <= 0) &
        error stop 'Logarithmic shift needs positive SI'

      target = exp((1-shift_fraction)*log(smoothed(valley)) &
                 + shift_fraction*log(smoothed(left_peak)))
    else
      target = smoothed(valley) &
             + shift_fraction*(smoothed(left_peak)-smoothed(valley))
    end if

    ! K decreases with array index: move toward the smaller-K peak.
    best = left_peak
    do j = valley, left_peak
      if (smoothed(j) >= target) then
        best = j
        exit
      end if
    end do
  end subroutine select_threshold

end program susceptibility_index