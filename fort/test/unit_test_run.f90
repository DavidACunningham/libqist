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

    ! call test_mmult(testpass)
    ! call test_basic_tensorops_d(testpass)
    ! call test_basic_tensorops_q(testpass)
    ! call test_scalar_cheby(testpass)
    ! call rot_hist_test(testpass)
    ! call quat_ops_test(testpass)
    ! call sh_test(testpass)
    ! call test_spice_subset("/home/david/wrk/nstgro/qist/kernels/mk_test.tf",testpass)
    ! call run_frk_tests(testpass)
    ! call run_frk_tests_q(testpass)
    ! call test_make_spice_subset(testpass)
    ! call test_make_rotation(testpass)
    ! call test_grav_load(testpass)
    ! call test_gq_init(testpass)
    ! call test_kep_grav(testpass)
    ! call test_threebody_grav(testpass)
    ! call test_SH_grav(testpass)
    ! call test_kernel_write(testpass)
    ! call end_to_end_integration_test(testpass)
    call test_make_qist(testpass)
end program main
