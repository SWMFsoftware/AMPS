! These fixtures use the same by-reference scalar/array calling convention as
! the legacy AMPS C++ declarations.  The first routine also pins the historical
! lower-case/trailing-underscore linker contract.  The second uses ISO_C_BINDING
! to provide a standards-defined control case on the same compiler.

subroutine amps_fortran_legacy_echo(tag, scalar, vector, matrix, &
                                    tag_out, scalar_out, vector_out, probes_out)
  use, intrinsic :: iso_c_binding, only: c_int, c_double
  implicit none
  integer(c_int), intent(in) :: tag
  real(c_double), intent(in) :: scalar
  real(c_double), intent(in) :: vector(3)
  real(c_double), intent(in) :: matrix(3, 3)
  integer(c_int), intent(out) :: tag_out
  real(c_double), intent(out) :: scalar_out
  real(c_double), intent(out) :: vector_out(3)
  real(c_double), intent(out) :: probes_out(3)

  tag_out = tag + 101_c_int
  scalar_out = 4.0_c_double * scalar - real(tag, c_double)
  vector_out(1) = vector(3) + real(tag, c_double)
  vector_out(2) = vector(1) - 2.0_c_double * vector(2)
  vector_out(3) = vector(2) + scalar
  probes_out(1) = matrix(1, 1) + 100.0_c_double * matrix(3, 2)
  probes_out(2) = matrix(2, 3) - matrix(1, 2)
  probes_out(3) = matrix(3, 1) + matrix(2, 2)
end subroutine amps_fortran_legacy_echo

subroutine amps_fortran_bindc_echo(tag, scalar, vector, matrix, &
                                   tag_out, scalar_out, vector_out, probes_out) &
    bind(C, name="amps_fortran_bindc_echo")
  use, intrinsic :: iso_c_binding, only: c_int, c_double
  implicit none
  integer(c_int), intent(in) :: tag
  real(c_double), intent(in) :: scalar
  real(c_double), intent(in) :: vector(3)
  real(c_double), intent(in) :: matrix(3, 3)
  integer(c_int), intent(out) :: tag_out
  real(c_double), intent(out) :: scalar_out
  real(c_double), intent(out) :: vector_out(3)
  real(c_double), intent(out) :: probes_out(3)

  tag_out = tag - 9_c_int
  scalar_out = scalar + 0.5_c_double * real(tag, c_double)
  vector_out(1) = vector(2) - scalar
  vector_out(2) = vector(3) + vector(1)
  vector_out(3) = vector(1) + vector(2) + vector(3)
  probes_out(1) = matrix(3, 3) - matrix(1, 1)
  probes_out(2) = matrix(2, 1) + 10.0_c_double * matrix(1, 3)
  probes_out(3) = matrix(3, 2) - matrix(2, 2)
end subroutine amps_fortran_bindc_echo
