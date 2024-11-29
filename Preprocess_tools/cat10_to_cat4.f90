!*************************************************************************************************!
!    Program that converts 10 column data to 4 column data with steps                                                                           !
!    1. Import data with format [Year, Month, Day, Hour, Min, Sec, Longitude, Latitude, Dep, Mag] !                            !
!    2. Transform to julian dates (seconds)                                                       !
!    3. Compute the elapsed time since first event                                                !                
!    4. Export catalog with format [Elapsed_time, Longitude, Latitude, Mag]                       !                                     !
!*************************************************************************************************!
program cat10_to_cat4  

    implicit none
    
    !!!!! File Processing variables !!!!!
    character(len=*), parameter :: input_file = 'datasets/input_catalog.txt'
    character(len=256), parameter :: output_cat='results/output_catalog.txt'
    !!!!! Catalog variables !!!!!
    integer, parameter :: max_events = 2000000
    real*8 :: first_event_time, times(max_events)
    integer :: ios,ii,t1,t2,t3,t4,t5
    real*8 :: t6,x,y,d,qe
    integer :: n_event,i,jj,unit_number
    real*8 :: time(max_events),mag(max_events),xxe(max_events),yye(max_events)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                        IMPORT DATA                             !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Open the input file
    open(unit=35, file=input_file, status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print*, "Error opening file:", input_file
        stop
    end if

    ! Read data and parse time
    jj=0
    do ii=1,2000000
      ! year,month,day,hour,min,sec, longitude (deg), latitude (deg), depth, magnitude
      read(35,*,end=199)t1,t2,t3,t4,t5,t6,x,y,d,qe
            jj=jj+1
            xxe(jj)=x
            yye(jj)=y
            times(jj)=date_to_seconds(t1,t2,t3,t4,t5,t6)
            mag(jj)=qe
    enddo
    199 close(35)
    n_event = jj-1
    print*, "Number of events",n_event

    ! Calculate the time of the first event
    first_event_time = times(1)

    ! Write the transformed data to output file
    do i = 1, n_event
        time(i) = times(i) - first_event_time
    end do

    unit_number = 88
    open(unit=unit_number, file=output_cat, status='replace', action='write')
    write(unit_number,*)time(1),xxe(1),yye(1),mag(1)
    do i=2,n_event
      write(unit_number,*)time(i),xxe(i),yye(i),mag(i)
   enddo
   close(unit_number)

   contains

   ! Function to convert date and time to seconds
   real(8) function date_to_seconds(yy, mm, dd, hh, mins, ss)
   integer, intent(in) :: yy, mm, dd, hh, mins
   real(8), intent(in) :: ss
   integer :: days_in_months(12) = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
   integer :: total_days, i

   ! Check for leap year
   if ((mod(yy, 4) == 0 .and. mod(yy, 100) /= 0) .or. mod(yy, 400) == 0) then
       days_in_months(2) = 29
   end if

   ! Calculate total days up to this date in the year
   total_days = (yy - 1970) * 365 + (yy - 1969) / 4 - (yy - 1901) / 100 + (yy - 1601) / 400
   do i = 1, mm - 1
       total_days = total_days + days_in_months(i)
   end do
   total_days = total_days + dd - 1

   ! Convert to seconds
   date_to_seconds = total_days * 86400.0 + hh * 3600.0 + mins * 60.0 + ss
  end function date_to_seconds

end program
