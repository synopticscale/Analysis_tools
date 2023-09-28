module evaluate_fss_fortran 
  implicit none

contains

subroutine evaluate_fss(model_probabilities, obs_probabilities, roi, fss, bias, false_alarm, pod)
  implicit none
  real(8), dimension(:,:) :: model_probabilities, obs_probabilities
  integer, dimension(3) :: roi
  real(8), dimension(:), allocatable :: FN, FP, TN, TP, O, M
  real(8), dimension(:,:,:) :: possible
  integer :: i, j, r, ri
  real(8) :: a_in(:,:), b_in(:,:)
  real(8) :: obcount(:), modcount(:), fneg(:), fpos(:), tneg(:), tpos(:)
  real(8) :: MSE, MSE_ref
  logical :: mask_valid
  real(8), dimension(3) :: fss, bias, false_alarm, pod

  FN = []  ! False negatives
  FP = []  ! False positives
  TN = []  ! True negatives
  TP = []  ! True positives
  O = []   ! Observed fractions
  M = []   ! Modeled fractions

  allocate(possible(0:len(roi), size(model_probabilities, 1), size(model_probabilities, 2)))

  ! Check if FSS_MASK is not 'none'
  mask_valid = .true.
  if (FSS_MASK /= 'none') then
    ! Load MASK from FSS_MASK (assuming FSS_MASK is a character variable)
    call load_mask(MASK, FSS_MASK)
  end if

  do i = 35, 264  ! Change to allow variability in domain size
    do j = 35, 264
      if (FSS_MASK /= 'none') then
        ! Check if the mask value is less than 0.5
        mask_valid = MASK(j, i) >= 0.5
        if (.not. mask_valid) then
          cycle  ! Skip to the next iteration if the mask value is invalid
        end if
      end if

      ! Loop through roi's
      obcount = 0.0d0
      modcount = 0.0d0
      fneg = 0.0d0
      fpos = 0.0d0
      tneg = 0.0d0
      tpos = 0.0d0
      ri = 1

      do r = 1, size(roi)
        ! Take subset of points around the current one
        a_in = obs_probabilities(max(j - roi(r), 1):min(j + roi(r), 299), max(i - roi(r), 1):min(i + roi(r), 299))
        b_in = model_probabilities(max(j - roi(r), 1):min(j + roi(r), 299), max(i - roi(r), 1):min(i + roi(r), 299))

        ! Calculate the number of valid elements
        possible(ri, j, i) = count_non_nan(a_in)

        ! If > 50% of FSS window is valid, record fractions
        if (possible(ri, j, i) >= 0.5d0 * ((2 * real(roi(r)) + 1.0d0)**2)) then
          obcount(ri) = sum(a_in) / possible(ri, j, i)
          modcount(ri) = sum(b_in) / possible(ri, j, i)

          fneg(ri) = sum(merge(1.0d0 - b_in, 0.0d0, (a_in == 1.0d0) .and. (b_in < 1.0d0)))
          fpos(ri) = sum(merge(b_in, 0.0d0, (a_in == 0.0d0) .and. (b_in > 0.0d0)))
          tneg(ri) = sum(merge(1.0d0 - b_in, 0.0d0, (a_in == 0.0d0) .and. (b_in < 1.0d0)))
          tpos(ri) = sum(merge(b_in, 0.0d0, (a_in == 1.0d0) .and. (b_in > 0.0d0)))
        else
          ! If < 50% of FSS window is valid, record NAN's instead of fractions
          obcount(ri) = NaN
          modcount(ri) = NaN
          fneg(ri) = NaN
          fpos(ri) = NaN
          tneg(ri) = NaN
          tpos(ri) = NaN
        end if

        ri = ri + 1
      end do

      ! Append the results to the arrays
      O = [O, obcount]
      M = [M, modcount]
      FN = [FN, fneg]
      FP = [FP, fpos]
      TN = [TN, tneg]
      TP = [TP, tpos]
    end do
  end do

  ! Calculate FSS, bias, pod, and false_alarm
  fss = 0.0d0
  bias = 0.0d0
  pod = 0.0d0
  false_alarm = 0.0d0

  do i = 1, size(roi)
    pod(i) = 1.0d0 - mean(FN(:, i)) / mean(FN(:, i) + TP(:, i))
    false_alarm(i) = mean(FP(:, i)) / mean(FP(:, i) + TN(:, i))

    MSE = sum((O(:, i) - M(:, i))**2) / count_non_nan(O

end module evaluate_fss_fortran 

