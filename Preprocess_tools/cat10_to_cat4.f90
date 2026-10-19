! Convert [year, month, day, hour, minute, second, longitude, latitude, depth, magnitude]
! to [elapsed_seconds, latitude, longitude, magnitude] for the SI tool.
! Input rows must be chronological; elapsed time starts at the first event.
program cat10_to_cat4
    implicit none

    character(len=*), parameter :: input_file = 'datasets/input_catalog.txt'
    character(len=*), parameter :: output_cat = 'results/output_catalog.txt'
    integer, parameter :: max_events = 2000000
    integer :: ios, i, n_event, year, month, day, hour, minute
    real(8) :: second, x, y, depth, magnitude, first_event_time
    real(8) :: times(max_events), lon(max_events), lat(max_events), mag(max_events)

    open(35, file=input_file, status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print *, 'Error opening file: ', input_file
        stop 1
    end if

    n_event = 0
    do i = 1, max_events
        read(35, *, end=100) year, month, day, hour, minute, second, x, y, depth, magnitude
        n_event = n_event + 1
        times(n_event) = date_to_seconds(year, month, day, hour, minute, second)
        lon(n_event) = x
        lat(n_event) = y
        mag(n_event) = magnitude
    end do
100 close(35)

    if (n_event == 0) error stop 'Input catalog is empty'
    print *, 'Number of events', n_event
    first_event_time = times(1)

    open(88, file=output_cat, status='replace', action='write')
    do i = 1, n_event
        write(88, *) times(i) - first_event_time, lat(i), lon(i), mag(i)
    end do
    close(88)

contains

    real(8) function date_to_seconds(yy, mm, dd, hh, mins, ss)
        integer, intent(in) :: yy, mm, dd, hh, mins
        real(8), intent(in) :: ss
        integer :: days_in_months(12), total_days, i

        ! Reset every call, including after a leap year.
        days_in_months = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
        if ((mod(yy, 4) == 0 .and. mod(yy, 100) /= 0) .or. mod(yy, 400) == 0) then
            days_in_months(2) = 29
        end if

        ! Days since 1970-01-01, counting leap years through the previous year.
        total_days = (yy - 1970)*365 + (yy - 1)/4 - 1969/4 &
                     - (yy - 1)/100 + 1969/100 + (yy - 1)/400 - 1969/400
        do i = 1, mm - 1
            total_days = total_days + days_in_months(i)
        end do
        total_days = total_days + dd - 1

        date_to_seconds = total_days*86400d0 + hh*3600d0 + mins*60d0 + ss
    end function date_to_seconds

end program cat10_to_cat4
