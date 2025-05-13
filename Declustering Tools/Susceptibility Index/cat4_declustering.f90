!*************************************************************************************************!
!    Program with steps                                                                           !
!    1. Import data with format [time, longitude, latitude, magnitude]                            !
!    2. Compute Susceptibility index                                                              !
!    3. Compute the minimum with smoothing algorithm                                              !                
!    4. Exports i) The catalog ii) Susceptibility index iii) True and predicted num of clusters   !
!               iv) rescaled distances                                                            !
!*************************************************************************************************!
program susceptibility_index  

    implicit none
    
    !!!!! File Processing variables !!!!!
    !! Spatial boundaries need to be defined a-priori - Magnitude is not filtered
    character(len=*), parameter :: input_file = 'datasets/input_data.txt'
    character(len=500), parameter :: rescaled_dist = 'results/rescaled_distances.txt'
    character(len=256), parameter :: susc_index_file_name = 'results/susc_index.txt'
    character(len=256), parameter :: output_cat='results/output_catalog.txt'
    character(len=256), parameter ::final_results='results/final_results.txt'
    
    !!!!! Catalog variables !!!!!
    integer, parameter :: max_events = 2000000,max_thr = 6000
    real*8, parameter :: mcut = 3.0
    integer :: i,ii,j,jj,n_event,cl_id(max_events),p_id(max_events)
    real*8 :: time(max_events),mag(max_events),xxe(max_events),yye(max_events)
    real*8 :: t,x,y,qe

    !!!!! Susceptibility Index variables !!!!!
    real*8, parameter:: pthmin=1E-18,pthmax=1E15,bin0=0.05
    integer :: ipmin, ipmax,ip,kr(-max_thr:max_thr),nlink(-max_thr:max_thr),nd(-max_thr:max_thr)
    integer :: ipmaxt,ipf,smooth_step
    real*8 :: ixxx,iyyy,l1
    real*8 :: zmexp
    
    !!!!! Minimum value variables !!!!!
    integer, parameter :: max_points = 10000
    real*8 :: zmin, sindex(max_points),zmin_shft,z_max1,z_max2 
    real*8 :: z_min, zmax1,zmax2,thresh(max_events)
    integer :: sw, num_points
    integer :: k(max_points),kmin_shft,kmin,kmax,kmax2
  
    !!!!! METRICS !!!!
    real*8 :: dt2,dr2,maxdt(max_events),maxdr(max_events),prc
    integer :: id(max_events),id_pr(max_events),ip_thr

    !!!!! Counter variables !!!!
    integer :: unit_number
    integer :: start_clock, end_clock, rate
    real(8) :: elapsed_time,total_time
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                        IMPORT DATA                             !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    open(35,file=input_file,status='old')
    
    jj=0
    do ii=1,max_events
      ! elapsed times since origin (seconds),longitude (deg), latitude (deg), magnitude
      read(35,*,end=199)t,x,y,qe
      if(t.gt.0) then
         if(qe.ge.mcut)then
            jj=jj+1
            xxe(jj)=x
            yye(jj)=y
            time(jj)=t
            mag(jj)=qe
         endif
      end if
    end do
    199 close(35)
    n_event = jj-1
    print*, "Number of events",n_event

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                        !
    !                 SUSCEPTIBILITY INDEX                   ! 
    !                                                        !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!! Compute the proximities among elements !!!!!!!!!!

    ipmin=int(dlog(pthmin)/bin0) ! Set minimum threshold index
    ipmax=int(dlog(pthmax)/bin0) ! Set maximum threshold index
    do ip=ipmin,ipmax
       nlink(ip)=0
       kr(ip)=0
       nd(ip)=0
    end do
    unit_number = 85
    open (unit_number,file=rescaled_dist,status="replace",action="write") 
    prc = 0.025
    ! Get the clock rate and start time
    elapsed_time = 0.
    call SYSTEM_CLOCK(start_clock, rate)
    do i=2,n_event  
      if (i / real(n_event) >= prc)then
         ! Get the end time
         call SYSTEM_CLOCK(end_clock)
         ! Calculate elapsed time
         elapsed_time = real(end_clock - start_clock, 8) / rate + elapsed_time
         total_time = elapsed_time / prc 
         write(15,*) prc, i, n_event, total_time - elapsed_time
         prc = prc + 0.025
         print*,prc,i
         call flush(15)
         ! Get the clock rate and start time
         call SYSTEM_CLOCK(start_clock, rate)
      end if
      ixxx=xxe(i)
      iyyy=yye(i)
      ipmaxt = ipmin
      do j=i-1,1,-1
            ! Compute similarity sij
            call metric(l1,ixxx,iyyy,time(i),xxe(j),yye(j),time(j),mag(j),dt2,dr2)
            if(int(dlog(l1)/bin0).lt.ipmin.or.int(dlog(l1)/bin0).gt.ipmax)then 
               print*,int(dlog(l1)/bin0),'WARNING change limits'
               print*,l1,ixxx,iyyy,time(i),xxe(j),yye(j),time(j),mag(j),dt2,dr2,i,j
               pause
            endif
            ipf=max(ipmin,int(log(l1)/bin0)) ! find m1: sth(m1)>sij>sth(m1+1)    
            if(ipf.gt.ipmaxt)then
                maxdt(i)=dt2
                maxdr(i)=dr2
            endif
            ipmaxt=max(ipf,ipmaxt)           ! find the maximum similarity s_i(max)  
            do ip=ipmin,ipf
                nlink(ip)=nlink(ip)+1
            end do 
       end do
       
       write(unit_number,*)maxdt(i),maxdr(i)
       nd(ipmaxt) = nd(ipmaxt)+1
       id_pr(i) = ipmaxt
       do ip=ipmaxt+1,ipmax
        kr(ip)=kr(ip)+1                      ! no links for sth>s_i(max)
       end do
    end do
    close(unit_number)
    close(15)
   pause

   !**********Compute and export the Susceptibility Index *************!
   unit_number = 86
   open (unit_number,file=susc_index_file_name,status="replace",action="write")
   num_points = 0
   smooth_step = 10
   do ip=ipmax,ipmin+smooth_step,-1
      if(kr(ip).eq.1)exit
      if(nlink(ip-smooth_step).gt.nlink(ip))then
         num_points = num_points + 1
         k(num_points) = kr(ip)                         ! Number of clusters
         zmexp=(kr(ip)-kr(ip-smooth_step))*1d0/(nlink(ip-smooth_step)-nlink(ip))
         sindex(num_points) = zmexp*nlink(ip)*1./n_event ! SSM index
         thresh(kr(ip)) =  exp(ip*bin0)
         !                 "Clusters",    "Susc index",      "NN" ,  "Threshold", "Thr index", "No links"
         write(unit_number,*)kr(ip)  ,  sindex(num_points), nd(ip),  exp(ip*bin0),    ip     ,  nlink(ip)
      endif    
   end do
   close(unit_number)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !                                                        !
   !                 FIND MINIMUM VALUE AND EXPORT          ! 
   !                                                        !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   !! Smoothing algorithm !!
   call min_ssm_shft(num_points,sindex,k,kmin_shft,zmin_shft,sw,z_max1,z_max2,z_min,kmin,zmax1,kmax,zmax2,&
                     kmax2,zmin)

   ! Open the file with the desired name
   unit_number = 87
   open(unit=unit_number, file=final_results, status='replace', action='write')

   ! Write the data to the file
   write(unit_number, *) kmin_shft, zmin_shft, thresh(kmin_shft), kmin, z_min, thresh(kmin), &
                           kmax, thresh(kmax), kmax2, thresh(kmax2)

   ! Close the file
   close(unit_number)

   !!!!!! Export catalog with id !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   unit_number = 88
   open(unit=unit_number, file=output_cat, status='replace', action='write')
   ip_thr = int(log(thresh(kmin_shft))/bin0)
   id(1)=0
   write(unit_number,*)time(1),xxe(1),yye(1),mag(1),id(1)
   do i=2,n_event
      if(id_pr(i).ge.ip_thr)then
         id(i) = 1  ! triggered
      elseif(id_pr(i).lt.ip_thr)then
         id(i) = 0  ! background
      endif
      write(unit_number,*)time(i),xxe(i),yye(i),mag(i),id(i)
   enddo
   close(unit_number)

   contains

   subroutine metric(aaa,x1,y1,t1,x3,y3,t3,q3,dt2,dr2)
        implicit none
                
        real*8 :: AAA,dr,prob,dt,dt0,dt2,dr2
        real*8 :: x1,y1,t1,x3,y3,t3,q3,pr,dr0,df,bval
        parameter(pr=3.14159/180.)
        integer :: flag1
        !real*8 :: cc,pp,K,a,gamma,dd,q,qc,espo,fdr,dm,l1

        dt0 = 1    ! 1 sec
        dr0 = 1E-3 ! in degrees 
        
        if(x1.eq.x3.and.y1.eq.y3)then
           x3=x3+dr0
           y3=y3+dr0
        endif

        !! Convert degrees to radians
        lat1 = x1
        lat2 = x3
        lon1 = y1
        lon2 = y3
        phi1 = lat1 * pr
        phi2 = lat2 * pr
        lambda1 = lon1 * pr
        lambda2 = lon2 * pr
        
        !! Harvesine distance !!
         ! ! Calculate differences
         ! dphi = phi2 - phi1
         ! dlambda = lambda2 - lambda1
         ! ! compute distance
         ! a = sin(dphi / 2.0d0)**2 + cos(phi1) * cos(phi2) * sin(dlambda / 2.0d0)**2
         ! dr = 2.0d0 * asin(sqrt(a))   ! angular distance in radians
         ! ! dr = dr * (1 / pr)         ! in degrees
         ! dr=dr*6370                   ! in km
        
        !! Distance !!
         dr = acos(sin(phi1)*sin(phi2)+cos(phi1)*cos(phi2)&
               *cos(lambda1-lambda2))
         dr = dr * 6370. ! in km
         dr = dr * 1000. ! in mt
        
        !! Euclidean distance !!
        !dr=sqrt((x1-x3)**2+(y1-y3)**2) ! euclidean distance
        
        ! Temporal distance
        dt=(t1-t3)

        flag1 = 0
        dr0 = 10 ! mt
        if(dr.eq.0.or.dt.eq.0)flag1=1
        dt = dt+dt0*flag1
        dr = dr+dr0*flag1
        
        !!!! Baiesi-Paszuski metric !!!!
        df = 1.6   ! fractal dimension
        bval = 1.0 ! b-value
        prob=1/(dt*(dr**df)*10**(-bval*q3)) ! similarity in sec*mt
        !dt2 = dt*(10**(-bval*q3*0.5))                     ! in seconds
        dt2 = (dt/(365.25*24.*3600.))*(10**(-bval*q3*0.5)) ! in years
        !dr2 = (dr**df)*(10**(-bval*q3*0.5))               ! in mt
        dr = dr / 1000.
        dr2 = (dr**df)*(10**(-bval*q3*0.5))                ! in km
        aaa=prob
   end subroutine metric

   subroutine min_ssm_shft(num_points,sindex,k,kmin_shft,zmin_shft,sw,z_max1,z_max2,z_min,kmin,zmax1,kmax,zmax2,&
                              kmax2,zmin)
      ! Compute the minimum of the SI and the shifted one by 30% based on the left peak of the curve
      !
      !! INPUT !! 
      ! num_points: Number of thresholds,     sindex:     Susceptibility Index/threshold
      ! k:          Number clusters/threshold
      
      !! OUTPUT !!
      ! kmin_shft:  Number clusters-Shft min  zmin_shft:  Shifted SI
      ! sw:         The final smoothing factor
      ! z_max1:     SI value of 1st max       z_max2:     SI value of 2nd max
      ! z_min:      SI value of the minimum   kmin:       K value of the minimum
      ! zmax1:      Smoothed SI of 1st max    kmax:       K value of the 1st max
      ! zmax2:      Smoothed SI of 2nd max    kmax2:      K value of the 2nd max
      ! zmin:       Smoothed SI of minimum
      !!!!!

      implicit none
      integer, parameter :: max_points = 10000 ! Adjust this to the maximum number of data points
      real*8  :: sindex(max_points),smoothed_data(max_points) ! Adjust the data type (real or integer) as needed
      real*8  :: zmin,zbest,zmax2,smooth_limit
      integer :: maxima_indices(max_points),k(max_points)
      integer :: minima_indices(max_points),kmax,kmax2,sw
      integer :: num_maxima,num_minima,imax,imax2,imin,ibest
      integer :: num_points, kmin,kmin_shft,tot_max
      real*8 :: z_max1,z_max2,z_min
      real*8 :: zmax1,min_min,max_max,zmin_shft
      
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
           exit
        endif
     enddo

     kmin_shft = k(ibest)          ! num of clusters for shifted min
     zmin_shft = sindex(ibest)     ! original ssm_index

      
   end subroutine min_ssm_shft

end program
