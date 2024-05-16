program main
    use genqist, only: make_qist_model
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp=>real128
    implicit none
    character(len=1000) :: namelist_filepath
    namelist_filepath = "/home/david/wrk/nstgro/qist/libqist/fort/clean_building/qist_config.nml"
    call make_qist_model(namelist_filepath)
end program main
