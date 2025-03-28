program main
    use test_globals, only: test_mmult
    use test_cheby, only: test_scalar_cheby, test_spice_subset
    use test_tensorops, only: test_basic_tensorops_d, test_basic_tensorops_q
    use test_quat, only: rot_hist_test, quat_ops_test
    use test_sh
    use test_frkmin
    use test_frkmin_q
    use test_genqist, only: test_make_spice_subset, &
                            test_make_rotation, &
                            test_grav_load, &
                            test_gq_init, &
                            test_kernel_write, &
                            test_make_qist
    use test_makemodel, only: test_kep_grav, &
                              test_threebody_grav, &
                              test_SH_grav, &
                              end_to_end_integration_test
    logical testpass

    write (*,*) "TEST START matrix multiplication functions. . ."
    call test_mmult(testpass)
    write (*,*) "TEST START DOUBLE tensor operations. . ."
    call test_basic_tensorops_d(testpass)
    write (*,*) "TEST START QUAD tensor operations. . ."
    call test_basic_tensorops_q(testpass)
    write (*,*) "TEST START Chebyshev interpolation operations. . ."
    call test_scalar_cheby(testpass)
    write (*,*) "TEST START rotation interpolation operations. . ."
    call rot_hist_test(testpass)
    write (*,*) "TEST START quaternion operations. . ."
    call quat_ops_test(testpass)
    write (*,*) "TEST START spherical harmonics operations. . ."
    call sh_test(testpass)
    write (*,*) "TEST START SPICE kernel resampling operations. . ."
    call test_spice_subset("../../kernels/mk_test.tf",testpass)
    ! write (*,*) "TEST START DOUBLE DOP853 (Runge Kutta integrator). . ."
    ! call run_frk_tests(testpass)
    ! write (*,*) "TEST START QUAD DOP853 (Runge Kutta integrator). . ."
    ! call run_frk_tests_q(testpass)

    write (*,*) "TEST SPICE kernel resampling operations. . ."
    call test_make_spice_subset(testpass)
    write (*,*) "TEST rotation kernel resampling operations. . ."
    call test_make_rotation(testpass)
    
    write (*,*) "TEST gravity model loading operations. . ."
    call test_grav_load(testpass)
    write (*,*) "TEST QIST model initialization operations. . ."
    call test_gq_init(testpass)
    write (*,*) "TEST keplerian gravity partials. . ."
    call test_kep_grav(testpass)
    write (*,*) "TEST third body perturbing gravity partials. . ."
    call test_threebody_grav(testpass)
    write (*,*) "TEST Spherical Harmonics perturbing gravity partials. . ."
    call test_SH_grav(testpass)
    write (*,*) "TEST SPICE kernel writing. . ."
    call test_kernel_write(testpass)
    write (*,*) "TEST kernel integration end-to-end. . ."
    call end_to_end_integration_test(testpass)
    write (*,*) "TEST QIST model end-to-end check. . ."
    call test_make_qist(testpass)
end program main
