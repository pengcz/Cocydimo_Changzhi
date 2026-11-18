program test_gas_refinement
  use m_solver_lib
  implicit none
  type(box_t) :: b
  integer :: flags(8,8)
  b%n_cell = 8
  b%dr = 0.001_dp
  b%r_min = [0.0_dp, 0.0_dp]
  b%tag = 0
  allocate(b%cc(0:9,0:9,1:5))
  i_E_norm = 1
  i_gas_vars = [2,3,4,5]
  b%cc = 0
  b%cc(:,:,2) = 1
  b%cc(:,:,5) = 1e5_dp/(gas_gamma-1)
  refine_max_dx = 1
  refine_min_dx = 1e-4_dp
  refine_field_threshold = 3e6_dp
  derefine_field_threshold = 2e6_dp
  derefine_dx = 0.01_dp
  global_time = 1e-7_dp
  refine_max_distance_head = 0.005_dp
  global_n_streamers = 0
  gas_refine_dx = 0.00025_dp
  gas_refine_radius = 0.002_dp
  call refinement_criterion(b, flags)
  if (any(flags /= af_do_ref)) error stop 'Axis region did not refine'
  b%r_min = [0.1_dp, 0.1_dp]
  call refinement_criterion(b, flags)
  if (any(flags /= af_rm_ref)) error stop 'Uniform distant gas did not coarsen'
  gas_inv_molecular_weight = 1
  gas_refine_rtol = 0.05_dp
  b%cc(0,:,5) = 1.2_dp*b%cc(0,:,5)
  call refinement_criterion(b, flags)
  if (any(flags /= af_do_ref)) error stop 'Gas feature across box boundary missed'
  b%dr = gas_refine_dx
  call refinement_criterion(b, flags)
  if (any(flags /= af_keep_ref)) error stop 'Resolved gas feature was coarsened'
  b%cc(:,:,5) = 1e5_dp/(gas_gamma-1)
  b%cc(0,:,5) = 1.04_dp*b%cc(0,:,5)
  call refinement_criterion(b, flags)
  if (any(flags /= af_keep_ref)) error stop 'Hysteresis failed'
  ! A nearly uniform warm region must refine even without a large gradient.
  gas_refine_rtol = 0
  gas_refine_temperature = 301
  gas_inv_molecular_weight = 1/(28.97_dp * Da)
  b%dr = 0.001_dp
  b%cc(:,:,5) = 302 * k_b * gas_inv_molecular_weight / (gas_gamma-1)
  call refinement_criterion(b, flags)
  if (any(flags /= af_do_ref)) error stop 'Warm region did not refine'
  b%cc(:,:,5) = 300.75_dp * k_b * gas_inv_molecular_weight / (gas_gamma-1)
  call refinement_criterion(b, flags)
  if (any(flags /= af_keep_ref)) error stop 'Temperature hysteresis failed'
  b%cc(:,:,5) = 300 * k_b * gas_inv_molecular_weight / (gas_gamma-1)
  call refinement_criterion(b, flags)
  if (any(flags /= af_rm_ref)) error stop 'Ambient gas did not coarsen'
  b%cc(0,:,5) = 302 * k_b * gas_inv_molecular_weight / (gas_gamma-1)
  call refinement_criterion(b, flags)
  if (any(flags /= af_do_ref)) error stop 'Warm boundary buffer did not refine'
  gas_refine_dx = 0
  focus_dx = 1e-4_dp
  focus_lo = [0.0_dp, 0.055_dp]
  focus_hi = [0.001_dp, 0.065_dp]
  b%r_min = [0.0_dp, 0.056_dp]
  b%dr = 0.001_dp
  call refinement_criterion(b, flags)
  if (any(flags /= af_do_ref)) error stop 'Fixed tip region did not refine'
  b%dr = focus_dx
  call refinement_criterion(b, flags)
  if (any(flags /= af_keep_ref)) error stop 'Fixed tip region coarsened'
  b%r_min = [0.01_dp, 0.056_dp]
  call refinement_criterion(b, flags)
  if (any(flags /= af_rm_ref)) error stop 'Distant box incorrectly protected'
  print *, 'Gas refinement checks passed'
end program
