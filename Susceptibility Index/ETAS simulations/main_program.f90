program ssm_simulations
  use global_parameters
  implicit none

  ! ===== File handling =====
  character(len=256) :: name

  ! ===== ETAS variables =====
  real(4) :: kk,alpha,beta,cc_days,pp,DD,gamma,qdec,bg_rate,n_branch,bb
  real(8) :: final_time, t_start
  real(8) :: xxxx, yyyy
  integer :: seed1, nevent, jjj, lll
  integer :: ZBQLPOI

  ! Background coordinate pools
  integer :: nev1, nevc, jj, jk
  real(8) :: xxe(max_bg_in), yye(max_bg_in)

  ! Simulation arrays
  real(8) :: q(q_buf_size), qq, rr, ran2
  real(8) :: muaa, deltar, qq0, ll, time, tt0, theta, x
  real(8) :: xx(max_events), yy(max_events), xx0, yy0, itime(max_events)
  integer :: event_id(max_events), parent_id(max_events)
  integer :: nmain, n_main, n_off, nn, j, i
  integer :: y_true, n_c

  ! Target arrays (declustering input)
  real(8) :: time_c(max_events), q_c(max_events), lat_c(max_events), lon_c(max_events)
  integer :: event_id_c(max_events), parent_id_c(max_events), icluster_c(max_events)
  integer :: icluster(max_events)

  ! SI variables
  integer :: ipmin, ipmax, ip
  integer :: kr(-max_thr:max_thr), nlink(-max_thr:max_thr), nd(-max_thr:max_thr)
  integer :: ipmaxt, ipf, nd_bg(-max_thr:max_thr), nd_tr(-max_thr:max_thr), k_max, pflag
  integer :: thres_ind(max_points)
  real(8) :: ixxx, iyyy, l1, zmexp, der(max_points)
  real(8) :: sindex(max_points), thresh(max_events)

  ! Classification metrics arrays
  real(8) :: recall(max_points), preci(max_points), f1(max_points), acc(max_points), f1_true
  real(8) :: prec_shft, rec_shft, f1_shft, f1_min, acc_shft
  integer :: tn(-max_thr:max_thr), fn(-max_thr:max_thr), fp(-max_thr:max_thr), tp(-max_thr:max_thr)
  integer :: iflag1, iflag2
  real(8) :: f1_max, ratio

  ! Metric outputs (rescaled)
  real(8) :: dt2, dr2

  ! Minimum extraction outputs
  real(8) :: zmin, zmin_shft, z_max1, z_max2, z_min, zmax1, zmax2
  integer :: smw, num_points, kmin, kmax, kmax2, kmin_shft
  integer :: k(max_points)

  ! Export classification
  integer :: id(max_events), id_pr(max_events), ip_thr
  integer :: unit_summary, unit_comb
  integer :: unit_ssm, unit_cat

  ! Timing
  integer :: count_rate, start_clock, finish_clock
  real(8) :: elapsed_min

  seed1 = seed_default

  call system_clock(count_rate=count_rate)
  call system_clock(start_clock)

  ! Ensure output directories exist (best effort; harmless if already there)
  call ensure_dir(results_dir)
  call ensure_dir(ssm_index_dir)
  call ensure_dir(sim_catalog_dir)

  unit_summary = 94
  open(unit=unit_summary, file=summary_file, status="replace", action="write")

  unit_comb = 30
  open(unit=unit_comb, file=combinations_file, status="replace", action="write")

  ! -------------------- Load background coordinate pools --------------------
  call load_bg_coords(bg_coords_m4_file, bg_coords_m2_file, xxe, yye, nev1, nevc)

  do lll = 1, num_sets
    write(name, "('comb_',I0)") lll
    print *, "Simulation", lll, "of", num_sets

    ! --- Generate ETAS parameters ---
    call param_gen(alpha, kk, pp, cc_days, DD, gamma, qdec, qmin, qsup, beta, bg_rate, n_branch)

    write(unit_comb,*) alpha,kk,pp,cc_days,DD,gamma,qdec,bg_rate,beta,n_branch
    call flush(unit_comb)

    ! --- Transform parameters ---
    bb = beta / log(10.0)                ! b-value GR
    cc_days = cc_days * seconds_per_day  ! now cc in seconds
    t_start  = seconds_per_year * tc_years
    final_time = seconds_per_year * tlast_years

    ! expected number of parents in auxiliary window
    nmain = int((bg_rate / seconds_per_day) * final_time * (x0max-x0min) * (y0max-y0min))

    ! init q buffer
    do i=1,q_buf_size
      q(i) = -100d0
    enddo

    ! -------------------- Generate nmain parent events --------------------
    nevent = 0
    n_main = 0
    jj = 0

    do jjj = 1, q_buf_size
      time = ran2(seed1) * final_time

      rr = ran2(seed1)
      qq = qmin - (1.0 / bb) * log10(1.0 - rr)

      if (qq >= 4.0) then
        jk = int(1 + ran2(seed1) * nev1)
      else
        jk = int(nev1 + 1 + ran2(seed1) * (nevc - nev1))
      endif

      xxxx = xxe(jk) + 2d0*coord_jitter_deg*ran2(seed1) - coord_jitter_deg
      yyyy = yye(jk) + 2d0*coord_jitter_deg*ran2(seed1) - coord_jitter_deg

      if (xxxx > x0min .and. xxxx < x0max .and. yyyy > y0min .and. yyyy < y0max) then
        jj = jj + 1
        nevent = nevent + 1

        xx(nevent) = xxxx
        yy(nevent) = yyyy
        itime(nevent) = time
        event_id(nevent) = nevent
        parent_id(nevent) = 0
        q(nevent) = min(qq, qsup)

        if (q(nevent) >= qc .and. time >= t_start .and. xx(nevent) >= xmin .and. &
            xx(nevent) <= xmax .and. yy(nevent) >= ymin .and. yy(nevent) <= ymax) then
          n_main = n_main + 1
        endif
      endif

      if (jj >= nmain) exit
    enddo

    ! -------------------- Generate offsprings --------------------
    n_off = 0
    do j=1,q_buf_size
      if (q(j) <= -10d0) exit

      muaa = kk * exp(alpha * (q(j) - qmin))
      nn = ZBQLPOI(muaa)

      tt0 = itime(j)
      xx0 = xx(j)
      yy0 = yy(j)
      qq0 = q(j)

      do i=1,nn
        time = cc_days * (ran2(seed1)**(1d0/(1d0-pp))) - cc_days
        if (tt0 + time > final_time) cycle

        ll = DD * exp(gamma * (qq0 - qmin))
        x = ran2(seed1)
        deltar = sqrt(ll * (x**(1d0/(1d0-qdec))) - ll)
        theta = ran2(seed1) * 2d0 * 3.141592653589793d0

        xxxx = xx0 + deltar * sin(theta)
        yyyy = yy0 + deltar * cos(theta)

        if (xxxx > x0min .and. xxxx < x0max .and. yyyy > y0min .and. yyyy < y0max) then
          nevent = nevent + 1

34        rr = ran2(seed1)
          qq = qmin - (1.0 / bb) * log10(1.0 - rr)
          if (qq > qsup) goto 34

          q(nevent) = min(qq, qsup)
          xx(nevent) = xxxx
          yy(nevent) = yyyy
          itime(nevent) = tt0 + time
          event_id(nevent) = nevent
          parent_id(nevent) = j

          if (q(nevent) >= qc .and. itime(nevent) >= t_start .and. xx(nevent) >= xmin .and. &
              xx(nevent) <= xmax .and. yy(nevent) >= ymin .and. yy(nevent) <= ymax) then
            n_off = n_off + 1
          endif
        endif
      enddo
    enddo

    ! Sort by time
    call HPSORT(nevent, itime, q, xx, yy, event_id, parent_id)

    ! Cluster ids (truth)
    call cluster_index(nevent, event_id, parent_id, icluster)

    ! Build target arrays
    y_true = 0
    n_c = 0
    do i=1,nevent
      if (q(i) >= qc .and. itime(i) >= t_start .and. xx(i) >= xmin .and. xx(i) <= xmax .and. &
          yy(i) >= ymin .and. yy(i) <= ymax) then

        n_c = n_c + 1
        time_c(n_c) = itime(i)
        q_c(n_c) = q(i)
        lon_c(n_c) = xx(i)
        lat_c(n_c) = yy(i)
        event_id_c(n_c) = event_id(i)
        parent_id_c(n_c) = parent_id(i)
        icluster_c(n_c) = icluster(i)
        if (parent_id(i) == 0) y_true = y_true + 1
      endif
    enddo

    ! -------------------- Susceptibility Index --------------------
    ipmin = int(dlog(pthmin)/bin0)
    ipmax = int(dlog(pthmax)/bin0)

    do ip=ipmin,ipmax
      nlink(ip)=0; kr(ip)=0; nd(ip)=0
      nd_tr(ip)=0; nd_bg(ip)=0
      tn(ip)=0; fn(ip)=0; fp(ip)=0; tp(ip)=0
    enddo

    do i=2,n_c
      pflag = 0
      if (parent_id_c(i) == 0) pflag = 1

      ixxx = lat_c(i)
      iyyy = lon_c(i)
      ipmaxt = ipmin

      do j=i-1,1,-1
        call metric(l1, ixxx, iyyy, time_c(i), lat_c(j), lon_c(j), time_c(j), q_c(j), dt2, dr2)

        ipf = max(ipmin, int(log(l1)/bin0))

        if (int(dlog(l1)/bin0) < ipmin .or. int(dlog(l1)/bin0) > ipmax) then
          print *, "WARNING: change limits (ip out of range). sim=", lll
        endif

        ipmaxt = max(ipf, ipmaxt)

        do ip=ipmin,ipf
          nlink(ip)=nlink(ip)+1
        enddo
      enddo

      id_pr(i) = ipmaxt
      nd(ipmaxt) = nd(ipmaxt) + 1
      nd_bg(ipmaxt) = nd_bg(ipmaxt) + pflag
      nd_tr(ipmaxt) = nd_tr(ipmaxt) + (1-pflag)

      do ip=ipmaxt+1,ipmax
        kr(ip) = kr(ip) + 1
        tn(ip) = tn(ip) + pflag
        fn(ip) = fn(ip) + (1-pflag)
      enddo

      do ip=ipmin,ipmaxt
        tp(ip) = tp(ip) + (1-pflag)
        fp(ip) = fp(ip) + pflag
      enddo
    enddo

    ! Compute SI curve + metrics
    num_points = 0
    f1_max = 0d0
    k_max  = 0
    iflag1 = 0
    iflag2 = 0

    do ip=ipmax, ipmin+smooth_step, -1
      if (kr(ip) == 1) exit
      if (nlink(ip-smooth_step) > nlink(ip)) then
        num_points = num_points + 1
        k(num_points) = kr(ip)

        zmexp = (kr(ip)-kr(ip-smooth_step))*1d0 / (nlink(ip-smooth_step)-nlink(ip))
        sindex(num_points) = zmexp * nlink(ip) * 1d0 / n_c
        der(num_points) = zmexp
        thresh(kr(ip)) = exp(ip*bin0)
        thres_ind(num_points) = ip

        if ((tn(ip)+fn(ip)) /= 0) iflag1 = 1
        if ((tn(ip)+fp(ip)) /= 0) iflag2 = 1

        preci(num_points) = (tn(ip)*1d0/(tn(ip)+fn(ip))) * iflag1
        recall(num_points)= (tn(ip)*1d0/(tn(ip)+fp(ip))) * iflag2

        f1(num_points) = (2d0*preci(num_points)*recall(num_points) / &
                          max(1d-30, (recall(num_points)+preci(num_points)))) * (iflag1*iflag2)

        acc(num_points)= (tp(ip)+tn(ip))*1d0 / max(1, (tn(ip)+tp(ip)+fp(ip)+fn(ip)))

        if (f1(num_points) > f1_max) then
          f1_max = f1(num_points)
          k_max = k(num_points)
        endif
      endif
    enddo

    f1_true = 0d0
    do i=2,num_points-1
      if (k(i-1) >= y_true .and. k(i+1) < y_true) f1_true = f1(i)
    enddo

    call min_ssm_shft(num_points, sindex, k, preci, recall, f1, acc, &
                      kmin_shft, zmin_shft, prec_shft, rec_shft, f1_shft, acc_shft, &
                      smw, z_max1, z_max2, z_min, kmin, f1_min, zmax1, kmax, zmax2, kmax2, zmin)

    write(unit_summary, *) n_c, y_true, f1_true, k_max, f1_max, kmin, kmin_shft, f1_min, f1_shft, &
                           prec_shft, rec_shft, acc_shft, ((n_c-y_true)*1d0)/(n_c*1d0), &
                           ((n_c-kmin_shft)*1d0)/(n_c*1d0), kmin_shft*1d0/(y_true*1d0), &
                           kmin*1d0/(y_true*1d0), k_max*1d0/(y_true*1d0)
    call flush(unit_summary)

    ratio = kmin_shft*1d0/(y_true*1d0)

    if (ratio > export_ratio_hi .or. ratio < export_ratio_lo) then
      unit_ssm = 93
      open(unit=unit_ssm, file=trim(ssm_index_dir)//"/ssm_index_"//trim(name)//".txt", status="replace", action="write")
      do i=1,num_points
        write(unit_ssm,*) k(i), sindex(i), nd(thres_ind(i)), nd_bg(thres_ind(i)), nd_tr(thres_ind(i)), &
                          nd(thres_ind(i))*1d0/n_c, nlink(thres_ind(i))*1d0/n_c, der(i), exp(thres_ind(i)*bin0), &
                          thres_ind(i), nlink(thres_ind(i)), preci(i), recall(i), f1(i), acc(i)
      enddo
      close(unit_ssm)

      unit_cat = 88
      open(unit=unit_cat, file=trim(sim_catalog_dir)//"/sim_cat_"//trim(name)//".txt", status="replace", action="write")

      ip_thr = int(log(thresh(kmin_shft))/bin0)
      id(1)=0

      write(unit_cat,*) time_c(1), lon_c(1), lat_c(1), q_c(1), event_id_c(1), parent_id_c(1), icluster_c(1), id(1)
      do i=2,n_c
        if (id_pr(i) >= ip_thr) then
          id(i)=1
        else
          id(i)=0
        endif
        write(unit_cat,*) time_c(i), lon_c(i), lat_c(i), q_c(i), event_id_c(i), parent_id_c(i), icluster_c(i), id(i)
      enddo
      close(unit_cat)
    endif

  enddo

  call system_clock(finish_clock)
  elapsed_min = (finish_clock - start_clock) * 1d0 / count_rate / 60d0
  print *, "Elapsed time:", elapsed_min, "minutes"

  close(unit_comb)
  close(unit_summary)

contains

  subroutine ensure_dir(dir)
    character(len=*), intent(in) :: dir
    character(len=512) :: cmd
    cmd = 'mkdir -p "'//trim(dir)//'"'
    call execute_command_line(cmd, wait=.true.)
  end subroutine ensure_dir

  subroutine load_bg_coords(file_m4, file_m2, xxe, yye, nev1, nevc)
    character(len=*), intent(in) :: file_m4, file_m2
    real(8), intent(out) :: xxe(max_bg_in), yye(max_bg_in)
    integer, intent(out) :: nev1, nevc
    integer :: ii, jj
    real(8) :: z, t, qe, x, y

    jj = 0
    open(36, file=file_m4, status="old")
    do ii=1,max_bg_in
      read(36,*,end=199) z, t, qe, x, y, z
      jj=jj+1
      xxe(jj)=x
      yye(jj)=y
    enddo
199 close(36)
    nev1 = jj

    open(35, file=file_m2, status="old")
    do ii=1,max_bg_in
      read(35,*,end=99) z, t, qe, x, y, z
      if (qe >= 1) then
        jj=jj+1
        xxe(jj)=x
        yye(jj)=y
      endif
    enddo
99  close(35)
    nevc = jj
  end subroutine load_bg_coords

  subroutine param_gen(a,kk,pp,cc,DD,gamma,qdec,mo,magn_max,beta_out,bg_rate_out,n_branch)
    real(4), intent(out) :: pp, cc, a, kk, gamma, dd, qdec, beta_out, bg_rate_out, n_branch
    real(4), intent(in)  :: mo, magn_max
    real(4) :: u

    beta_out = beta_val

20  call random_number(u); pp = pp_min + (pp_max-pp_min)*u
    call random_number(u); cc = cc_min_days + (cc_max_days-cc_min_days)*u
    call random_number(u); a  = alpha_min + (alpha_max-alpha_min)*u
    call random_number(u); kk = kk_min + (kk_max-kk_min)*u
    call random_number(u); gamma = gamma_min + (gamma_max-gamma_min)*u
    call random_number(u); dd = dd_min + (dd_max-dd_min)*u
    call random_number(u); qdec = qdec_min + (qdec_max-qdec_min)*u

    if (a == beta_out) then
      n_branch = (kk * beta_out * (magn_max - mo)) / (1 - exp(-beta_out * (magn_max - mo)))
    else
      n_branch = (kk * beta_out * (1 - exp(-(magn_max - mo) * (beta_out - a)))) / &
                 ((beta_out - a) * (1 - exp(-beta_out * (magn_max - mo))))
    endif
    if (n_branch <= nbranch_min .or. n_branch >= nbranch_max) goto 20

    call random_number(u)
    bg_rate_out = bg_rate_min + (bg_rate_max-bg_rate_min)*u
  end subroutine param_gen

  subroutine HPSORT(N,time,mag,lon,lat,event_id,parent_id)
    integer, intent(in) :: N
    real*8, intent(inout) :: mag(N), time(N), lon(N), lat(N)
    integer, intent(inout) :: event_id(N), parent_id(N)
    real*8 :: RRA, lonRA, latRA, qRA
    integer :: IR, L, I, J, event_idRA, parent_idRA

    L = N/2 + 1
    IR = N

10  continue
    if (L > 1) then
      L = L - 1
      RRA = time(L); qRA = mag(L); lonRA = lon(L); latRA = lat(L)
      event_idRA = event_id(L); parent_idRA = parent_id(L)
    else
      RRA = time(IR); qRA = mag(IR); lonRA = lon(IR); latRA = lat(IR)
      event_idRA = event_id(IR); parent_idRA = parent_id(IR)
      time(IR)=time(1); mag(IR)=mag(1); lon(IR)=lon(1); lat(IR)=lat(1)
      event_id(IR)=event_id(1); parent_id(IR)=parent_id(1)
      IR = IR - 1
      if (IR == 1) then
        time(1)=RRA; mag(1)=qRA; lon(1)=lonRA; lat(1)=latRA
        event_id(1)=event_idRA; parent_id(1)=parent_idRA
        return
      endif
    endif

    I = L
    J = L + L
20  if (J <= IR) then
      if (J < IR) then
        if (time(J) < time(J+1)) J = J + 1
      endif
      if (RRA < time(J)) then
        time(I)=time(J); mag(I)=mag(J); lon(I)=lon(J); lat(I)=lat(J)
        event_id(I)=event_id(J); parent_id(I)=parent_id(J)
        I=J; J=J+J
      else
        J = IR + 1
      endif
      goto 20
    endif

    time(I)=RRA; mag(I)=qRA; lon(I)=lonRA; lat(I)=latRA
    event_id(I)=event_idRA; parent_id(I)=parent_idRA
    goto 10
  end subroutine HPSORT

  subroutine cluster_index(N,event_id,parent_id,icluster)
    integer, intent(in) :: N
    integer, intent(in) :: event_id(N), parent_id(N)
    integer, intent(out) :: icluster(N)
    integer :: ncluster, i
    integer :: ixcluster(N)

    ncluster = 0
    do i=1,N-1
      if (parent_id(i) == 0) then
        ncluster = ncluster + 1
        icluster(i) = ncluster
      else
        icluster(i) = ixcluster(parent_id(i))
      endif
      ixcluster(event_id(i)) = icluster(i)
    enddo
  end subroutine cluster_index

  subroutine metric(aaa,x1,y1,t1,x3,y3,t3,q3,dt2,dr2)
    use global_parameters, only: df, bval, seconds_per_year
    implicit none
    real(8), intent(out) :: aaa, dt2, dr2
    real(8), intent(in)  :: x1,y1,t1,x3,y3,t3,q3

    real(8), parameter :: pr=3.141592653589793238d0/180d0
    real(8) :: dr, prob, dt
    real(8) :: lat1,lat2,lon1,lon2,phi1,phi2,lambda1,lambda2
    real(8) :: dt0, dr0
    integer :: flag1

    dt0 = 1d0     ! 1 second
    dr0 = 0.01d0  ! km (small shift if identical)

    lat1=x1; lat2=x3; lon1=y1; lon2=y3
    phi1=lat1*pr; phi2=lat2*pr; lambda1=lon1*pr; lambda2=lon2*pr

    dr = acos(sin(phi1)*sin(phi2)+cos(phi1)*cos(phi2)*cos(lambda1-lambda2))
    dr = dr * 6370d0  ! km

    dt = (t1 - t3)

    flag1 = 0
    if (dr == 0d0 .or. dt <= 0d0) flag1=1
    dt = dt + dt0*flag1
    dr = dr + dr0*flag1
    ! Similarity in 1/(second*km^df) units, rescaled by the magnitude term
    prob = 1d0/(dt*(dr**df)*10d0**(-bval*q3))
    dt2  = (dt/seconds_per_year) * (10d0**(-bval*q3*0.5d0))
    dr2  = (dr**df) * (10d0**(-bval*q3*0.5d0))
    aaa  = prob
  end subroutine metric
  
  subroutine min_ssm_shft(num_points,sindex,k,preci,recall,f1,acc,kmin_shft,zmin_shft,prec_shft,rec_shft&
                            ,f1_shft,acc_shft,sw,z_max1,z_max2,z_min,kmin,f1_min,zmax1,kmax,zmax2,kmax2,zmin)
  !!!!!
  ! Compute the minimum of the SI and the shifted one by 30% based on the left peak of the curve
  !
  !! INPUT !! 
  ! num_points: Number of thresholds,     sindex:     Susceptibility Index/threshold
  ! k:          Number clusters/threshold preci:      Precision for each threshold
  ! recall:     Recall/threshold          f1:         F-1 score/threshold
  
  !! OUTPUT !!
  ! kmin_shft:  Number clusters-Shft min  zmin_shft:  Shifted SI
  ! prec_shft:  Precision-Shft min        rec_shft:   Recall-shft min f1_shft: F-1 score-shft min 
  ! sw:         The final smoothing factor
  ! z_max1:     SI value of 1st max       z_max2:     SI value of 2nd max
  ! z_min:      SI value of the minimum   kmin:       K value of the minimum
  ! zmax1:      Smoothed SI of 1st max    kmax:       K value of the 1st max
  ! zmax2:      Smoothed SI of 2nd max    kmax2:      K value of the 2nd max
  ! zmin:       Smoothed SI of minimum
  !!!!!
    use global_parameters
    implicit none
    real*8  :: sindex(max_points),smoothed_data(max_points) 
    real*8  :: zmin,z_min,zmin_shft,zbest,zmax1,z_max1,zmax2,z_max2
    integer :: maxima_indices(max_points),k(max_points)
    integer :: minima_indices(max_points),sw
    integer :: num_maxima,num_minima,imax,imax2,imin,ibest
    integer :: num_points,kmin,kmax,kmax2,kmin_shft,tot_max
    real*8 :: preci(max_points),recall(max_points),f1(max_points),acc(max_points)
    real*8 :: prec_shft,rec_shft,f1_shft,max_max,min_min
    real*8 :: prec_min,rec_min,f1_min,smooth_limit,acc_shft

      smooth_limit = 0.10 ! 10% of the total number of points
      tot_max = k(1)      ! maximum number of clusters
      ! loop over smoothing step until 2 maximum 1 minimum are reached or end of iterations
      do sw=2,int(num_points*smooth_limit) 
        num_maxima=0
        num_minima=0
        
        ! Apply the smoothing window from right to left of the SI curve
        do i = 1+sw/2, num_points-sw/2
           smoothed_data(i) = 0.0
           do j=max(1,i-sw/2),min(num_points,i+sw/2)
              smoothed_data(i) = smoothed_data(i)+sindex(j)
           end do
           smoothed_data(i)=smoothed_data(i)*1./sw
        end do

        ! Find all local maxima and minima in the smoothed data
        do i = 1+sw/2, num_points-sw/2
           ! neglect max's and min's that lie in the first and last 2% of total events
           if(smoothed_data(i)>smoothed_data(i-1).and.&
              smoothed_data(i)>smoothed_data(i + 1).and.&
              k(i).gt.int(tot_max*0.02).and.k(i).lt.int(tot_max*0.98)&
              )then
              num_maxima = num_maxima + 1
              maxima_indices(num_maxima)=i
           elseif(smoothed_data(i)<smoothed_data(i-1).and.&
              smoothed_data(i)<smoothed_data(i + 1).and.&
              k(i).gt.int(tot_max*0.05).and.k(i).lt.int(tot_max*0.95)&
              )then
              num_minima = num_minima + 1
              minima_indices(num_minima)=i
           end if
        end do

        ! Increase smoothing factor until the condition is reached
        if(num_minima==1.and.num_maxima==2)then

           imin=minima_indices(1)       ! index of minimum
           zmin=smoothed_data(imin)     ! smoothed ssm_index
           z_min = sindex(imin)         ! original ssm_index
           kmin=k(imin)                 ! number of clusters
                      
           imax=maxima_indices(2)       ! index of max1 - left peak
           kmax=k(imax)                 ! number of clusters for max1
           zmax1=smoothed_data(imax)    ! smoothed ssm_index
           z_max1 = sindex(imax)        ! original ssm_index
           
           imax2=maxima_indices(1)      ! index of max2 - right peak
           kmax2=k(imax2)               ! num clusters for max2
           zmax2=smoothed_data(imax2)   ! smoothed ssm_index
           z_max2 = sindex(imax2)       ! original ssm_index
           
           prec_min=preci(imin)         ! precision of minimum
           rec_min=recall(imin)         ! recall of minimum  
           f1_min=f1(imin)              ! f1 of minimum
          
           exit ! exit loop over smoothing factor
        endif

        ! if smoothing limit is reached
        if(sw.eq.int(num_points*smooth_limit))then
            
            ! find the global minimum
            min_min=1E10
            do i = 1, num_minima
               if(smoothed_data(minima_indices(i)).lt.min_min)then
                  min_min=smoothed_data(minima_indices(i))
                  imin=minima_indices(i)
                  prec_min=preci(imin)
                  rec_min=recall(imin)
                  f1_min=f1(imin)
               endif
            end do
            zmin=smoothed_data(imin)  ! smoothed ssm_index
            kmin=k(imin)              ! number of clusters
            z_min=sindex(imin)        ! original ssm_index
            
            ! find the right peak
            max_max=-1
            do i=1+sw/2,imin-1
               if(smoothed_data(i).gt.max_max)then
                  kmax2=k(i)
                  max_max=smoothed_data(i)
                  imax2=i
               endif
            enddo
            zmax2=max_max             ! smoothed ssm_index - right peak
            z_max2=sindex(imax2)      ! origianl ssm_index - right peak

            ! find left peak
            max_max=-1
            do i=imin+1,num_points-sw/2
               if(smoothed_data(i).gt.max_max)then
                  kmax=k(i)
                  max_max=smoothed_data(i)
                  imax=i
               endif
            enddo
            zmax1=max_max            ! smoothed ssm_index - left peak
            z_max1=sindex(imax)      ! original ssm_index - left peak

            exit ! exit loop over smoothing factor
        
          endif
      enddo
        
      ! Shift minimum 30% higher with respect to left peak - zmax1
      zbest=exp(log(zmin)+(log(zmax1)-log(zmin))*0.3)
      do i=imin,num_points-sw/2
        if(smoothed_data(i).ge.zbest)then
           ibest=i
           prec_shft=preci(ibest)   ! precision of shifted minimum
           rec_shft=recall(ibest)   ! recall of shifted minimum
           f1_shft=f1(ibest)        ! f1 of shifted minimum
           acc_shft=acc(ibest)
           exit
        endif
     enddo

     kmin_shft = k(ibest)          ! num of clusters for shifted min
     zmin_shft = sindex(ibest)     ! original ssm_index
     
  end subroutine min_ssm_shft
  
end program ssm_simulations
  
  FUNCTION ran2(idum)
    INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
    REAL*8 ran2,AM,EPS,RNMX
    PARAMETER (IM1=2147483563,IM2=2147483399,AM=1d0/IM1,IMM1=IM1-1,&
    IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,&
    NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
    INTEGER idum2,j,k,iv(NTAB),iy
    SAVE iv,iy,idum2
    DATA idum2/123456789/, iv/NTAB*0/, iy/0/
    
    if (idum.le.0) then
      idum=max(-idum,1)
      idum2=idum
      do 11 j=NTAB+8,1,-1
        k=idum/IQ1
        idum=IA1*(idum-k*IQ1)-k*IR1
        if (idum.lt.0) idum=idum+IM1
        if (j.le.NTAB) iv(j)=idum
  11      continue
      iy=iv(1)
    endif
    k=idum/IQ1
    idum=IA1*(idum-k*IQ1)-k*IR1
    if (idum.lt.0) idum=idum+IM1
    k=idum2/IQ2
    idum2=IA2*(idum2-k*IQ2)-k*IR2
    if (idum2.lt.0) idum2=idum2+IM2
    j=1+iy/NDIV
    iy=iv(j)-idum2
    iv(j)=idum
    if(iy.lt.1)iy=iy+IMM1
    ran2=min(AM*iy,RNMX)
    return
  END FUNCTION ran2
