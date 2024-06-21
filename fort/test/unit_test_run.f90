program main
    use test_globals, only: test_mmult
    use test_cheby, only: test_scalar_cheby, test_spice_subset
    use test_tensorops, only: test_basic_tensorops_d, test_basic_tensorops_q
    use test_quat, only: rot_hist_test
    logical testpass

    call test_mmult(testpass)
    call test_basic_tensorops_d(testpass)
    call test_basic_tensorops_q(testpass)
    call test_scalar_cheby(testpass)
    call rot_hist_test(testpass)
    call test_spice_subset("/home/david/wrk/nstgro/qist/kernels/mk.tf",testpass)

end program main
