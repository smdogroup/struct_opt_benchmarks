module quicksort
  ! A module for quick-sorting things

  contains 
  subroutine quicksortPartition(array, mid)
    ! Partition an array of values into two for the quicksort algorithm
    ! This will be called separately on the left and right segments of
    ! the array.
    !
    ! Input:
    ! array:  the array of integers
    !
    ! Output:
    ! mid:    the position in the array immediatley after the partition

    implicit none

    integer, intent(inout), dimension(:) :: array
    integer, intent(out) :: mid
    integer :: i, j, temp, pivot

    ! Select the pivot
    pivot = array(1)

    ! Set i to the le
    i = 0
    j = size(array) + 1

    do
       ! Search from the last element of the array backwards until we
       ! have an element violating the pivot point
       j = j-1
       do
          if (array(j) <= pivot) exit
          j = j-1
       end do

       ! Search from the first element of the array until we have a
       ! point violating the pivot
       i = i+1
       do
          if (array(i) >= pivot) exit
          i = i+1
       end do

       ! Now we can either exchange entries or exit if no
       ! pivots were found
       if (i < j) then
          ! Exchange the entries to partition the array
          temp = array(i)
          array(i) = array(j)
          array(j) = temp
       else if (i == j) then
          mid = i+1
          return
       else
          mid = i
          return
       end if
    end do

  end subroutine quicksortPartition

  recursive subroutine quicksortArray(array)
    ! Sort an array of integers using a recursive quick-sort
    ! implementation
    !
    ! Input:
    ! array:   An array of integers to be sorted

    implicit none

    integer, intent(inout), dimension(:) :: array
    integer :: mid

    ! If the array is larger than 1, partition it 
    if (size(array) > 1) then
       call quicksortPartition(array, mid)

       ! Sort the left- and right-hand partitions of the array
       call quicksortArray(array(:mid-1))
       call quicksortArray(array(mid:))
    endif

  end subroutine quicksortArray

end module quicksort
