! Title: frkmin_q.f90 
! Description:
!     Quad precision implementation of DOP853 integrator tailored for use
!     with QIST. Based on SciPy implementation of DOP853.
!
!     PORTIONS OF THIS SOURCE FILE WERE PORTED FROM SCIPY. THE FOLLOWING TEXT 
!     (BSD 3-CLAUSE LICENSE) IS INCLUDED FOR LICENSING COMPLIANCE AND APPLIES
!     ONLY TO PORTIONS PORTED FROM SCIPY. THE REST OF THE SOURCE CODE IN THIS
!     FILE IS COVERED BY THE COPYRIGHT AND LICENSE PROVIDED IN THE ROOT 
!     DIRECTORY OF THE ASSOCIATED REPOSITORY.
!     START LICENSE FOR SCIPY PORTIONS ~~~~~~~~~~~~~~
!          Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
!          All rights reserved.
!          
!          Redistribution and use in source and binary forms, with or without
!          modification, are permitted provided that the following conditions
!          are met:
!          
!          1. Redistributions of source code must retain the above copyright
!             notice, this list of conditions and the following disclaimer.
!          
!          2. Redistributions in binary form must reproduce the above
!             copyright notice, this list of conditions and the following
!             disclaimer in the documentation and/or other materials provided
!             with the distribution.
!          
!          3. Neither the name of the copyright holder nor the names of its
!             contributors may be used to endorse or promote products derived
!             from this software without specific prior written permission.
!          
!          THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!          "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!          LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
!          A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
!          OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
!          SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
!          LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
!          DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
!          THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!          (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!          OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!     END LICENSE FOR SCIPY PORTIONS ~~~~~~~~~~~~~~
!
! References:
!     Hairer, E.: Solving Ordinary Differential Equations I. Springer Series 
!         in Computational Mathematics, vol. 8. Springer, Berlin, 
!         Heidelberg (1993)
!     Virtanen, P.. et. al.: Fundamental algorithms for scientific computing
!         in Python. Nat. Methods 17, 261–272 (2020). 
!         https://doi.org/10.1038/s41592-019-0686-2
! 
! author: David Cunningham
! Last edited: See git log
module frkmin_q
    use, intrinsic :: iso_fortran_env, only: real64, real128
    implicit none
    ! Basic parameters for shapes and precision
    integer, parameter :: DP=real64, &
                    & QP=real128, &
                    & WP=QP, & ! CHOOSE PRECISION HERE
                    & N_STAGES = 12, &
                    & N_STAGES_EXTENDED = 16, &
                    & INTERPOLATOR_POWER = 7
    ! Machine epsilon needed for step sizing
    real(WP), parameter :: de= epsilon(1._wp), &
                       nan = huge(nan)
    real(WP), parameter :: &
    !~~~BUTCHER TABLEAU VALUES BEGIN HERE~~~
    ! C: Values prescribing intermediate independent variable values
    C(0:N_STAGES_EXTENDED-1) = [0.0_wp, &
                             0.526001519587677318785587544488e-01_wp, &
                             0.789002279381515978178381316732e-01_wp, &
                             0.118350341907227396726757197510_wp, &
                             0.281649658092772603273242802490_wp, &
                             0.333333333333333333333333333333_wp, &
                             0.25_wp, &
                             0.307692307692307692307692307692_wp, &
                             0.651282051282051282051282051282_wp, &
                             0.6_wp, &
                             0.857142857142857142857142857142_wp, &
                             1.0_wp, &
                             1.0_wp, &
                             0.1_wp, &
                             0.2_wp, &
                             0.777777777777777777777777777778_wp], &
    ! A: Values weighting intermediate function calls 
    Ar0(0:N_STAGES_EXTENDED-1) =   0._WP, &
    Ar1(0:N_STAGES_EXTENDED-1) = [ &
                                5.26001519587677318785587544488e-2_wp, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP  &
                                ], &
    Ar2(0:N_STAGES_EXTENDED-1) = [ &
                                1.97250569845378994544595329183e-2_WP, &
                                5.91751709536136983633785987549e-2_WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP  &
                                ], &
    Ar3(0:N_STAGES_EXTENDED-1) = [ &
                                2.95875854768068491816892993775e-2_WP, &
                                0._WP, &
                                8.87627564304205475450678981324e-2_WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP  &
                                ], &
    Ar4(0:N_STAGES_EXTENDED-1) = [ &
                                2.41365134159266685502369798665e-1_WP, &
                                0._WP, &
                               -8.84549479328286085344864962717e-1_WP, &
                                9.24834003261792003115737966543e-1_WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP  &
                                ], &
    Ar5(0:N_STAGES_EXTENDED-1) = [ &
                                3.7037037037037037037037037037e-2_WP, &
                                0._WP, &
                                0._WP, &
                                1.70828608729473871279604482173e-1_WP, &
                                1.25467687566822425016691814123e-1_WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP  &
                                ], &
    Ar6(0:N_STAGES_EXTENDED-1) = [ &
                                3.7109375e-2_WP, &
                                0._WP, &
                                0._WP, &
                                1.70252211019544039314978060272e-1_WP, &
                                6.02165389804559606850219397283e-2_WP, &
                               -1.7578125e-2_WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP  &
                                ], &
    Ar7(0:N_STAGES_EXTENDED-1) = [ &
                                3.70920001185047927108779319836e-2_WP, &
                                0._WP, &
                                0._WP, &
                                1.70383925712239993810214054705e-1_WP, &
                                1.07262030446373284651809199168e-1_WP, &
                               -1.53194377486244017527936158236e-2_WP, &
                                8.27378916381402288758473766002e-3_WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP  &
                                ], &
    Ar8(0:N_STAGES_EXTENDED-1) = [ &
                                6.24110958716075717114429577812e-1_WP, &
                                0._WP, &
                                0._WP, &
                               -3.36089262944694129406857109825_WP, &
                               -8.68219346841726006818189891453e-1_WP, &
                                2.75920996994467083049415600797e1_WP, &
                                2.01540675504778934086186788979e1_WP, &
                               -4.34898841810699588477366255144e1_WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP  &
                                ], &
    Ar9(0:N_STAGES_EXTENDED-1) = [ &
                                4.77662536438264365890433908527e-1_WP, &
                                0._WP, &
                                0._WP, &
                               -2.48811461997166764192642586468_WP, &
                               -5.90290826836842996371446475743e-1_WP, &
                                2.12300514481811942347288949897e1_WP, &
                                1.52792336328824235832596922938e1_WP, &
                               -3.32882109689848629194453265587e1_WP, &
                               -2.03312017085086261358222928593e-2_WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP  &
                                ], &
    Ar10(0:N_STAGES_EXTENDED-1) = [&
                               -9.3714243008598732571704021658e-1_WP, &
                                0._WP, &
                                0._WP, &
                                5.18637242884406370830023853209_WP, &
                                1.09143734899672957818500254654_WP, &
                               -8.14978701074692612513997267357_WP, &
                               -1.85200656599969598641566180701e1_WP, &
                                2.27394870993505042818970056734e1_WP, &
                                2.49360555267965238987089396762_WP, &
                               -3.0467644718982195003823669022_WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP  &
                                ], &
    Ar11(0:N_STAGES_EXTENDED-1) = [&
                                2.27331014751653820792359768449_WP, &
                                0._WP, &
                                0._WP, &
                               -1.05344954667372501984066689879e1_WP, &
                               -2.00087205822486249909675718444_WP, &
                               -1.79589318631187989172765950534e1_WP, &
                                2.79488845294199600508499808837e1_WP, &
                               -2.85899827713502369474065508674_WP, &
                               -8.87285693353062954433549289258_WP, &
                                1.23605671757943030647266201528e1_WP, &
                                6.43392746015763530355970484046e-1_WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP  &
                                ], &
    Ar12(0:N_STAGES_EXTENDED-1) = [&
                                5.42937341165687622380535766363e-2_WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                4.45031289275240888144113950566_WP, &
                                1.89151789931450038304281599044_WP, &
                               -5.8012039600105847814672114227_WP, &
                                3.1116436695781989440891606237e-1_WP, &
                               -1.52160949662516078556178806805e-1_WP, &
                                2.01365400804030348374776537501e-1_WP, &
                                4.47106157277725905176885569043e-2_WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP  &
                                ], &
    Ar13(0:N_STAGES_EXTENDED-1) = [&
                                5.61675022830479523392909219681e-2_WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                2.53500210216624811088794765333e-1_WP, &
                               -2.46239037470802489917441475441e-1_WP, &
                               -1.24191423263816360469010140626e-1_WP, &
                                1.5329179827876569731206322685e-1_WP, &
                                8.20105229563468988491666602057e-3_WP, &
                                7.56789766054569976138603589584e-3_WP, &
                               -8.298e-3_WP, &
                                0._WP, &
                                0._WP, &
                                0._WP  &
                                ], &
    Ar14(0:N_STAGES_EXTENDED-1) = [&
                                3.18346481635021405060768473261e-2_WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                2.83009096723667755288322961402e-2_WP, &
                                5.35419883074385676223797384372e-2_WP, &
                               -5.49237485713909884646569340306e-2_WP, &
                                0._WP, &
                                0._WP, &
                               -1.08347328697249322858509316994e-4_WP, &
                                3.82571090835658412954920192323e-4_WP, &
                               -3.40465008687404560802977114492e-4_WP, &
                                1.41312443674632500278074618366e-1_WP, &
                                0._WP, &
                                0._WP  &
                                ], &
    Ar15(0:N_STAGES_EXTENDED-1) = [&
                               -4.28896301583791923408573538692e-1_WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                               -4.69762141536116384314449447206_WP, &
                                7.68342119606259904184240953878_WP, &
                                4.06898981839711007970213554331_WP, &
                                3.56727187455281109270669543021e-1_WP, &
                                0._WP, &
                                0._WP, &
                                0._WP, &
                               -1.39902416515901462129418009734e-3_WP, &
                                2.9475147891527723389556272149_WP, &
                               -9.15095847217987001081870187138_WP, &
                                0._WP  &
                                ], &
    A(0:N_STAGES_EXTENDED-1, 0:N_STAGES_EXTENDED-1) = reshape(&
     [ Ar0, Ar1, Ar2, Ar3, Ar4, Ar5, Ar6, Ar7, Ar8, & 
     & Ar9, Ar10, Ar11, Ar12, Ar13, Ar14, Ar15], &
     & [N_STAGES_EXTENDED,N_STAGES_EXTENDED],order=[2,1]), &
    ! B: Values weighting final contribution of computed stages
    B(0:N_STAGES-1) = A(N_STAGES, :N_STAGES-1), &
    ! Error stage coefficients
    E3(0:N_STAGES) = [ &
                  B(0) - 0.244094488188976377952755905512_WP, &
                  B(1:7), &
                  B(8) - 0.733846688281611857341361741547_WP, &
                  B(9:10), &
                  B(11) - 0.220588235294117647058823529412e-1_WP, &
                  0._WP &
                  ], &
    E5(0:N_STAGES) = [ &
                  0.1312004499419488073250102996e-1_WP, &
                  0._WP, &
                  0._WP, &
                  0._WP, &
                  0._WP, &
                 -0.1225156446376204440720569753e+1_WP, &
                 -0.4957589496572501915214079952_WP, &
                  0.1664377182454986536961530415e+1_WP, &
                 -0.3503288487499736816886487290_WP, &
                  0.3341791187130174790297318841_WP, &
                  0.8192320648511571246570742613e-1_WP, &
                 -0.2235530786388629525884427845e-1_WP, &
                  0._WP &
                  ], &
    ! Interpolation stage coefficients
    ! First 3 coefficients are computed separately.
    dr0(0:15) = [ &
             -0.84289382761090128651353491142e+1_WP, &
              0._WP, &
              0._WP, &
              0._WP, &
              0._WP, &
              0.56671495351937776962531783590_WP, &
             -0.30689499459498916912797304727e+1_WP, &
              0.23846676565120698287728149680e+1_WP, &
              0.21170345824450282767155149946e+1_WP, &
             -0.87139158377797299206789907490_WP, &
              0.22404374302607882758541771650e+1_WP, &
              0.63157877876946881815570249290_WP, &
             -0.88990336451333310820698117400e-1_WP, &
              0.18148505520854727256656404962e+2_WP, &
             -0.91946323924783554000451984436e+1_WP, &
             -0.44360363875948939664310572000e+1_WP &
             ],&
    dr1(0:15) = [ &
             0.10427508642579134603413151009e+2_WP, &
             0._WP, &
             0._WP, &
             0._WP, &
             0._WP, &
             0.24228349177525818288430175319e+3_WP, &
             0.16520045171727028198505394887e+3_WP, &
            -0.37454675472269020279518312152e+3_WP, &
            -0.22113666853125306036270938578e+2_WP, &
             0.77334326684722638389603898808e+1_WP, &
            -0.30674084731089398182061213626e+2_WP, &
            -0.93321305264302278729567221706e+1_WP, &
             0.15697238121770843886131091075e+2_WP, &
            -0.31139403219565177677282850411e+2_WP, &
            -0.93529243588444783865713862664e+1_WP, &
             0.35816841486394083752465898540e+2_WP &
             ],&
    dr2(0:15) = [ &
             0.19985053242002433820987653617e+2_WP, &
             0._WP, &
             0._WP, &
             0._WP, &
             0._WP, &
            -0.38703730874935176555105901742e+3_WP, &
            -0.18917813819516756882830838328e+3_WP, &
             0.52780815920542364900561016686e+3_WP, &
            -0.11573902539959630126141871134e+2_WP, &
             0.68812326946963000169666922661e+1_WP, &
            -0.10006050966910838403183860980e+1_WP, &
             0.77771377980534432092869265740_WP, &
            -0.27782057523535084065932004339e+1_WP, &
            -0.60196695231264120758267380846e+2_WP, &
             0.84320405506677161018159903784e+2_WP, &
             0.11992291136182789328035130030e+2_WP &
             ],&
    dr3(0:15) = [ &
            -0.25693933462703749003312586129e+2_WP, &
             0._WP, &
             0._WP, &
             0._WP, &
             0._WP, &
            -0.15418974869023643374053993627e+3_WP, &
            -0.23152937917604549567536039109e+3_WP, &
             0.35763911791061412378285349910e+3_WP, &
             0.93405324183624310003907691704e+2_WP, &
            -0.37458323136451633156875139351e+2_WP, &
             0.10409964950896230045147246184e+3_WP, &
             0.29840293426660503123344363579e+2_WP, &
            -0.43533456590011143754432175058e+2_WP, &
             0.96324553959188282948394950600e+2_WP, &
            -0.39177261675615439165231486172e+2_WP, &
            -0.14972683625798562581422125276e+3_WP &
             ], &
    D(0:INTERPOLATOR_POWER-4, 0:N_STAGES_EXTENDED-1) = reshape(&
     [ dr0, dr1, dr2, dr3 ], & 
     & [INTERPOLATOR_POWER-3,N_STAGES_EXTENDED],order=[2,1])
    real(WP), parameter :: SAFETY = 0.9_WP, &
                           MIN_FACTOR = 0.2_WP, &  ! Minimum allowed decrease in a step size.
                           MAX_FACTOR = 10_WP  ! Maximum allowed increase in a step size.
    character(len=15)  :: MESSAGES(0:1) = ["Solver success.", &
                                           "Solver failure."]
    type :: RungeKutta
        ! Type storing DOP853 initialization, step, error estimation, etc.
        ! Parameters
        ! ----------
        ! fun : function
        !     Right-hand side of the system. The calling signature is
        !     found below in the interface ``eomsig''
        !     corresponds to a single column in ``y``
        ! t0 : real
        !     Initial time.  
        ! y0 : real, shape (n)
        !     Initial state.
        ! t_bound : real
        ! max_step : real, optional
        !     Maximum allowed step size. Default is unbounded.
        ! rtol, atol : real, optional
        !     Relative and absolute tolerances. The solver keeps the local error
        !     estimates less than ``atol + rtol * abs(y)``. 
        !     1e-5 for `rtol` and 1e-8 for `atol`.

        ! Attributes
        ! ----------
        ! n : integer
        !     Number of equations.
        ! status : character
        !     Current status of the solver: 'running', 'finished' or 'failed'.
        ! t_bound : real
        !     Boundary time.
        ! direction : real
        !     Integration direction: +1 or -1.
        ! t : real
        !     Current time.
        ! y : real(n)
        !     Current state.
        ! t_old : real
        !     Previous time.
        ! step_size : real
        !     Size of the last successful step.
        ! nfev : integer
        !     Number evaluations of the system's right-hand side.
        character(len=56)     :: TOO_SMALL_STEP = "Required step size is less than spacing between numbers."
        character(len=20)     :: statmsg
        character(len=20)     :: message
        logical               :: success
        real(WP)              :: t_old, t, t_bound, &
                                 A(0:n_stages,0:n_stages) = A(:n_stages, :n_stages), &
                                 B(0:n_stages-1) = B, &
                                 C(0:n_stages) = C(:n_stages), &
                                 E3(0:n_stages) = E3, &
                                 E5(0:n_stages) = E5, &
                                 D(0:Interpolator_power-4, 0:N_stages_extended-1) = D, &
                                 A_EXTRA(0:n_stages_extended-n_stages - 2, 0:n_stages_extended -1) = A(n_stages + 1:,:), &
                                 C_EXTRA(0:n_stages_extended-n_stages - 2) = C(n_stages + 1:), &
                                 h_previous, error_exponent
        real(WP), allocatable :: y(:), y_old(:), rtol(:), atol(:), f(:)
        real(WP)              :: direction, max_step, h_abs
        integer               :: n
        integer               :: nfev
        procedure(eomsig), pointer :: eom => null()
        real(WP), allocatable :: K(:,:)
        real(WP), allocatable :: K_extended(:,:)
        integer               :: order = 8, &
                                 error_estimator_order = 7, &
                                 n_stages = N_STAGES, &
                                 n_stages_extended = N_STAGES_EXTENDED
        contains 
            procedure :: init => rkinit
            procedure estimate_error_norm
            procedure step_impl
            procedure :: step_size
            procedure :: step
            procedure :: dense_output => DOP853DenseOutputImpl
            procedure :: fn
    end type
    type  :: Dop853DenseOutput
        ! Type storing DOP853 solution
        real(WP) :: t_old, t, t_min, t_max, h
        real(WP), allocatable :: y_old(:), Frev(:,:)
        integer  :: n
        contains
            generic :: dcall => scall, mcall
            procedure :: scall
            procedure :: mcall
            procedure :: init => DOP853denseinit
            procedure :: call_impl => DOP853densecall
            procedure :: write => DOwrite
            procedure :: read => DOread
    end type
    type Odesolution
        !       Continuous ODE solution.
        !       A collection of `DenseOutput` objects which represent
        !       local interpolants. 
        !       The correct interpolant to use is determined by binary search.
        !       The interpolants cover the range between `t_min` and `t_max` (see
        !       Attributes below). Evaluation outside this interval is allowed, but
        !       the accuracy is not guaranteed.
        !       When evaluating at a breakpoint (one of the values in `ts`) a segment with
        !       the lower index is selected.

        !       Parameters
        !       ----------
        !       ts : real (n_segments + 1)
        !           Time instants between which local interpolants are defined. Must
        !           be strictly increasing or decreasing (zero segment with two points is
        !           also allowed).
        !       interpolants (n_segments) : 
        !           List of DenseOutput with n_segments elements
        !           Local interpolants. The i-th interpolant is assumed to be defined
        !           between ``ts(i)'' and ``ts(i + 1)''.

        !       Attributes
        !       ----------
        !       t_min, t_max : real
        !           Time range of the interpolation.
        integer                              :: n_segments, nfev, ndim
        real(WP),                allocatable :: ts(:), ts_sorted(:), ys(:,:)
        real(WP)                             :: t_min, t_max 
        logical                              :: ascending, dense
        character(20)                        :: statmsg
        type(DOP853DenseOutput), allocatable :: interpolants(:)
        contains
            procedure :: init => init_odesol
            generic   :: call => call_single, call_
            procedure :: call_single
            procedure :: call_
            procedure :: write =>solwrite
            procedure :: read => solread
    end type
    type ExtensibleDOArray
        ! A dynamic-length array for storing dense output objects
        ! Necessary for computing variable step length solutions,
        ! When you don't know beforehand how many steps the integrator
        ! Will take.
        ! Once the integration is done, the extensible array is 
        ! ``sealed'' and can be stored in a normal Fortran array.
        integer               :: length, allocatedSize
        type(DOP853DenseOutput), allocatable :: dataArray(:)
        contains
            procedure :: create => DOcreateflat ! start an array with an initial element
            procedure :: yield => DOyield
            procedure :: seal => DOsealflat
            procedure :: append => DOappend
            procedure :: doublesize => DOdoublesize
    end type ExtensibleDOArray
    type ExtensibleRealArray
        ! A dynamic-length array for storing real numbers.
        ! Necessary for computing variable step length solutions,
        ! When you don't know beforehand how many steps the integrator
        ! Will take.
        ! Once the integration is done, the extensible array is 
        ! ``sealed'' and can be stored in a normal Fortran array.
        integer               :: length, allocatedSize, vecdim
        real(WP), allocatable :: dataArray(:,:)

        contains

            generic :: create => createflat, createmat ! start an array with an initial column element
            procedure :: createflat
            procedure :: createmat
            procedure :: yield
            generic   :: seal => sealflat, sealmat
            procedure :: sealflat
            procedure :: sealmat
            procedure :: append
            procedure :: appendMultiple
            procedure :: doublesize
    end type ExtensibleRealArray
    interface 
        function eomsig(self,t,y) result(res)
            ! EOMs for the solve_ivp function MUST
            ! follow this function signature.
            import
            class(RungeKutta), intent(inout) :: self
            real(WP),            intent(in)    :: t
            real(WP),            intent(in)    :: y(:)
            real(WP)                           :: res(size(y))
        end function eomsig
    end interface
    contains
    ! ODESolution Procedures
    subroutine init_odesol(self, ts, ys, nfev, statmsg, interpolants)
        ! Initialize ode solution
        class(Odesolution),      intent(inout)        :: self
        real(WP),                intent(in)           :: ts(:), ys(:,:)
        integer,                 intent(in)           :: nfev
        character(len=20),       intent(in)           :: statmsg
        type(DOP853DenseOutput), intent(in), optional :: interpolants(:)
        real(WP)                                      :: d(size(ts)-1)

        d = diff(ts)
        self%statmsg = statmsg
        self%nfev = nfev
        ! The first case covers integration on zero segment.
        if (.not. ((size(ts) == 2 .and. ts(1) == ts(size(ts))) &
                .or. (all(d > 0) .or. all(d < 0)))) then
            print *, "`ts` must be strictly increasing or decreasing."
            error stop
        end if
        if (present(interpolants)) then
            self%n_segments = size(interpolants)
        else
            self%n_segments = size(ts)
        endif
        if (present(interpolants).and. (size(ts) .ne. (self%n_segments + 1))) then
            print *, "Numbers of time stamps and interpolants don't match."
            error stop
        endif

        allocate(self%ts(self%n_segments))
        allocate(self%ts_sorted(self%n_segments))
        allocate(self%ys,mold=ys)
        self%ts = 0._wp
        self%ys = 0._wp
        self%ts = ts
        self%ys = ys
        self%ndim = size(ys,1)
        if (present(interpolants)) then
            allocate(self%interpolants(self%n_segments))
            self%interpolants = interpolants
            self%dense = .true.
        else
            self%dense = .false.
        endif
        if (ts(self%n_segments) >= ts(1)) then
            self%t_min = ts(1)
            self%t_max = ts(self%n_segments)
            self%ascending = .True.
            self%ts_sorted = ts
            
        else
            self%t_min = ts(self%n_segments)
            self%t_max = ts(1)
            self%ascending = .False.
            self%ts_sorted = ts(self%n_segments:1:-1)
        endif
    end subroutine
    function call_single(self, t) result(res)
        ! Call ODE solution at time t
        class(Odesolution), intent(in) :: self
        real(WP),           intent(in) :: t
        real(WP)                       :: res(self%interpolants(1)%n)
        type(DOP853DenseOutput)        :: seg
        integer                        :: ind, segment
        ! Here we preserve a certain symmetry that when t is in ts,
        ! then we prioritize a segment with a lower index.
        ind = binsearch(self%ts_sorted, t)

        segment = min(max(ind, 1), self%n_segments)
        if (.not. self%ascending) then
            segment = self%n_segments - 1 - segment
        end if

        seg = self%interpolants(segment)
        res = seg%dcall(t)
    end function

    subroutine solwrite(self,unit_num,prec)
        ! Write ODE solution to disk. 
        ! the file handle given by unit_num must be opened
        ! with ``stream'' access and writable.
        class(Odesolution), intent(inout) :: self
        integer, intent(in)           :: unit_num
        integer, intent(in), optional :: prec
        real(WP)                      :: ts(size(self%ts)), &
                                       & ts_sorted(size(self%ts_sorted)), &
                                         ys(size(self%ys,1),size(self%ys,2))
        integer i
        integer  :: p
        if(present(prec)) then
            p = prec
        else
            p = WP
        endif
        ts = self%ts
        ts_sorted = self%ts_sorted
        ys = self%ys
        write(unit_num) self%ascending
        write(unit_num) self%dense
        write(unit_num) self%statmsg
        write(unit_num) self%ndim
        write(unit_num) self%n_segments
        write(unit_num) self%nfev
        select case (p)
        case (DP)
            write(unit_num) real(self%t_min,DP)
            write(unit_num) real(self%t_max,DP)
            write(unit_num) real(ts,DP)
            write(unit_num) real(ts_sorted,DP)
            write(unit_num) real(ys,DP)
        case (QP)
            write(unit_num) real(self%t_min,QP)
            write(unit_num) real(self%t_max,QP)
            write(unit_num) real(ts,QP)
            write(unit_num) real(ts_sorted,QP)
            write(unit_num) real(ys,QP)
        end select
        do i = 1,self%n_segments
            call self%interpolants(i)%write(unit_num, p)
        end do
    end subroutine solwrite
    subroutine solread(self,unit_num)
        ! Write ODE solution to disk. 
        ! the file handle given by unit_num must be opened
        ! with ``stream'' access and readable.
        class(Odesolution), intent(inout) :: self
        integer, intent(in)           :: unit_num
        integer i
        read(unit_num) self%ascending
        read(unit_num) self%dense
        read(unit_num) self%statmsg
        read(unit_num) self%ndim
        read(unit_num) self%n_segments
        read(unit_num) self%nfev
        read(unit_num) self%t_min
        read(unit_num) self%t_max
        if (allocated(self%ts)) deallocate(self%ts)
        if (allocated(self%ts_sorted)) deallocate(self%ts_sorted)
        if (allocated(self%ys)) deallocate(self%ys)
        if (allocated(self%interpolants)) deallocate(self%interpolants)
        allocate(self%ts(self%n_segments+1), &
                 self%ts_sorted(self%n_segments+1), &
                 self%ys(self%ndim, self%n_segments+1), &
                 self%interpolants(self%n_segments))
        read(unit_num) self%ts
        read(unit_num) self%ts_sorted
        read(unit_num) self%ys
        do i = 1,self%n_segments
            call self%interpolants(i)%read(unit_num)
        end do
    end subroutine solread

    function call_(self, t) result(res)
        !Evaluate the solution.
        !Parameters
        !----------
        !t : real array with shape (n_points)
        !    Points to evaluate at.
 
        !Returns
        !-------
        !y : real (n_states, n_points)

        ! t must be sorted, unlike the scipy version
        class(Odesolution), intent(in) :: self
        real(WP),           intent(in) :: t(:)
        real(WP)                       :: res(self%interpolants(1)%n,size(t))
        integer                        :: segments(size(t))
        type(DOP853DenseOutput)        :: seg(size(t))
        integer iter
        !$OMP PARALLEL DO
        do iter=1,size(t)
            segments(iter) = binsearch(self%ts_sorted, t(iter))
        end do
        !$OMP END PARALLEL DO
        where (segments < 0._WP) segments = 1
        where (segments > self%n_segments - 1) segments  = self%n_segments
        if (.not. self%ascending) then
            segments = self%n_segments - 1 - segments
        endif

        !$OMP PARALLEL DO
        do iter=1,size(t)
            seg(iter) = self%interpolants(segments(iter))
            res(:,iter) = seg(iter)%dcall(t(iter))
        end do
        !$OMP END PARALLEL DO
    end function call_

    function step_size(self) result(res)
        ! Compute the DOP853 step size
        ! given the current solver state
        class(RungeKutta), intent(in) :: self
        real(WP)                      :: res
        if (self%t_old==nan) then
            res = nan
        else
            res  = abs(self%t - self%t_old)
        end if
    end function step_size

    function step(self) result(res)
        ! Take a runge-kutta step given the current solver state
        class(RungeKutta), intent(inout) :: self
        character(len=20) :: res
        logical           :: success
        real(WP)          :: t

        if (self%statmsg.ne.'running') then
            print *, "Attempt to step on a failed or finished solver."
            error stop
        end if 
        if (self%n == 0 .or. self%t == self%t_bound) then
            ! Handle corner cases of empty solver or no integration.
            self%t_old = self%t
            self%t = self%t_bound
            res = "None"
            self%statmsg = "finished"
            success = .true.
        else
            t = self%t
            success = self%step_impl()
            res = self%statmsg

            if (.not.success) then
                self%statmsg = 'failed'
            else
                self%t_old = t
                if (self%direction * (self%t - self%t_bound) >= 0) then
                    self%statmsg = 'finished'
                endif
            end if
        endif
    end function step
    function fn(self, t, y) result(res)
        ! Wrapper for using the EOMs that increments
        ! the solver nfev counter
        class(RungeKutta), intent(inout) :: self
        real(WP),         intent(in)    :: t, y(:)
        real(WP)                        :: res(size(y))

        self%nfev = self%nfev+1
        res = self%eom(t,y)
    end function 
    ! ExtensibleDOArray Procedures
    subroutine DOcreateflat(self, element)
        class(ExtensibleDOArray), intent(inout) :: self
        type(DOP853DenseOutput),  optional,  intent(in)    :: element
        self%allocatedSize = 2
        if (allocated(self%dataArray)) deallocate(self%dataArray)
        allocate(self%dataArray(self%allocatedSize))
        if (present(element)) then
            self%length = 1
            self%dataArray(1) = element
        else
            self%length = 0
        endif
    end subroutine DOcreateflat
    subroutine DOsealflat(self,tgt, destroy)
        class(ExtensibleDOArray), intent(inout)  :: self
        type(DOP853DenseOutput), allocatable,      intent(out)    :: tgt(:)
        logical, optional,          intent(in)     :: destroy
        logical                                    :: dest
        if (allocated(tgt)) deallocate(tgt)
        allocate(tgt(self%length))
        tgt = self%dataArray(:self%length)
        dest = .true.
        if (present(destroy)) dest = destroy
        if (dest) deallocate(self%dataArray)
    end subroutine
    function DOyield(self,item) result(res)
        class(ExtensibleDOArray), intent(inout) :: self
        integer                   , intent(in)    :: item
        type(DOP853DenseOUtput)                   :: res
        if (item < self%length) then
            print *, "No data in element ", item
        else
            res = self%dataArray(item)
        end if
    end function
    subroutine DOappend(self, item)
        class(ExtensibleDOArray), intent(inout) :: self
        type(DOP853DenseOutput),  intent(in)    :: item
        if (self%length + 1 == self%allocatedSize) then
            call self%doublesize()
        endif
            self%dataArray(self%length+1) = item
            self%length = self%length+1
    end subroutine
    subroutine DOdoublesize(self)
        class(ExtensibleDOArray), intent(inout) :: self
        type(DOP853DenseOutput)                 :: dummyArray(self%allocatedSize)

        dummyArray = self%dataArray
        deallocate(self%dataArray)
        allocate(self%dataArray(2*self%allocatedSize))
        self%dataArray(:self%allocatedSize) = dummyArray
        self%allocatedSize = 2*self%allocatedSize
    end subroutine
    ! ExtensibleRealArray Procedures
    subroutine createflat(self, element, moldonly)
        class(ExtensibleRealArray), intent(inout) :: self
        real(WP),                   intent(in)    :: element
        logical, optional,          intent(in)    :: moldonly
        self%vecdim = 1
        if (present(moldonly)) then
            if (moldonly) then 
                self%length = 0
            else
                self%length = 1
            endif
        else
            self%length = 1
        endif
        self%allocatedSize = 2
        if (allocated(self%dataArray)) deallocate(self%dataArray)
        allocate(self%dataArray(self%vecdim, self%allocatedSize))
        self%dataArray = 0._wp
        self%dataArray(:,1) = element
    end subroutine createflat
    subroutine createmat(self, array, moldonly)
        class(ExtensibleRealArray), intent(inout) :: self
        real(WP),                   intent(in)    :: array(:)
        logical, optional,          intent(in)    :: moldonly
        self%vecdim = size(array)
        if (present(moldonly)) then
            if (moldonly) then 
                self%length = 0
            else
                self%length = 1
            endif
        else
            self%length = 1
        endif

        self%allocatedSize = 2
        if (allocated(self%dataArray)) deallocate(self%dataArray)
        allocate(self%dataArray(self%vecdim, self%allocatedSize))
        self%dataArray = 0._wp
        self%dataArray(:,1) = array
    end subroutine createmat
    subroutine sealflat(self,tgt, destroy)
        class(ExtensibleRealArray), intent(inout)  :: self
        real(WP), allocatable,      intent(out)    :: tgt(:)
        logical, optional,          intent(in)     :: destroy
        logical                                    :: dest
        if (self%vecdim.ne.1) then
            print *, "Error: must provide 2-d array to seal 2-d extensible array"
            error stop
        endif
        if (allocated(tgt)) deallocate(tgt)
        allocate(tgt(self%length))
        tgt = self%dataArray(1,:self%length)
        dest = .true.
        if (present(destroy)) dest = destroy
        if (dest) deallocate(self%dataArray)
    end subroutine
    subroutine sealmat(self,tgt, destroy)
        class(ExtensibleRealArray), intent(inout)  :: self
        real(WP), allocatable,      intent(out)    :: tgt(:,:)
        logical, optional,          intent(in)     :: destroy
        logical                                    :: dest
        if (allocated(tgt)) deallocate(tgt)
        allocate(tgt(self%vecdim,self%length))
        tgt = self%dataArray(:,:self%length)
        dest = .true.
        if (present(destroy)) dest = destroy
        if (dest) deallocate(self%dataArray)
    end subroutine
    function yield(self,item) result(res)
        class(ExtensibleRealArray), intent(inout) :: self
        integer                   , intent(in)    :: item
        real(WP)                                  :: res(self%vecdim)
        if (item < self%length) then
            print *, "No data in element ", item
        else
            res = self%dataArray(:,item)
        end if
    end function
    subroutine append(self, item)
        class(ExtensibleRealArray), intent(inout) :: self
        real(WP)                  , intent(in)    :: item(self%vecdim)
        if (self%length + 1 == self%allocatedSize) then
            call self%doublesize()
        endif
            self%dataArray(:,self%length+1) = item
            self%length = self%length+1
    end subroutine
    subroutine appendMultiple(self,items)
        class(ExtensibleRealArray), intent(inout) :: self
        real(WP)                  , intent(in)    :: items(:,:)
        integer i
        if (size(items,1).ne.self%vecdim) then
            print *, "to append, array must have ", self%vecdim, " rows, but an array with ", size(items,1), " rows was provided."
        else
            ! Make sure dataArray is long enough to handle all the data
            ! if not, double size enough times to accomodate additional length
            do while (self%allocatedSize.lt.self%length+size(items,2))
                call self%doublesize()
            end do
            ! Add all the data in parallel
            !$OMP PARALLEL DO
            do i=self%length, self%length+size(items,2)
                self%dataArray(:,i) = items(:,i)
            end do
            !$OMP END PARALLEL DO
            self%length = self%length + size(items,2)
        endif
    end subroutine
    subroutine doublesize(self)
        class(ExtensibleRealArray), intent(inout) :: self
        real(WP), allocatable                     :: dummyArray(:,:)

        allocate(dummyArray(self%vecdim,self%allocatedSize))
        dummyArray= self%dataArray
        deallocate(self%dataArray)
        allocate(self%dataArray(self%vecdim,2*self%allocatedSize))
        self%dataArray = 0._wp
        self%dataArray(:,:self%allocatedSize) = dummyArray
        self%allocatedSize = 2*self%allocatedSize
    end subroutine
    subroutine rk_step(solver, t, y, f, h, A, B, C, K, y_new, f_new)
        ! Perform a single Runge-Kutta step.

        ! This function computes a prediction of an explicit Runge-Kutta method and
        ! also estimates the error of a less accurate method.

        ! Parameters
        ! ----------
        ! solver : RungeKutta type
        ! t : real
        !     Current time.
        ! y : real (n)
        !     Current state.
        ! f : real (n)
        !     Current value of the derivative, i.e., ``fun(x, y)``.
        ! h : real
        !     Step to use.
        ! A : real(n_stages, n_stages)
        !     Coefficients for combining previous RK stages to compute the next
        !     stage. For explicit methods the coefficients at and above the main
        !     diagonal are zeros.
        ! B : real(n_stages,)
        !     Coefficients for combining RK stages for computing the final
        !     prediction.
        ! C : real(n_stages,)
        !     Coefficients for incrementing time for consecutive RK stages.
        !     The value for the first stage is always zero.
        ! K : real(n_stages + 1, n)
        !     Storage array for putting RK stages here. Stages are stored in rows.
        !     The last row is a linear combination of the previous rows with
        !     coefficients

        ! Returns
        ! -------
        ! y_new : real(n)
        !     Solution at t + h computed with a higher accuracy.
        ! f_new : real(n)
        !     Derivative ``fun(t + h, y_new)``.
        class(RungeKutta),   intent(inout) :: solver
        real(WP),           intent(in)    :: t, &
                                           & y(:), &
                                           & f(:), &
                                           & h, &
                                           & A(0:,0:), &
                                           & B(0:), &
                                           & C(0:) 
        real(WP),           intent(out)   :: y_new(size(y)), &
                                           & f_new(size(f)), &
                                           & K(0:,0:)
        real(WP)                          :: dy(size(y))
        integer s
        
        dy = 0._wp
        K = 0._wp
        K(0,:) = f
        do s=1,solver%N_STAGES-1
            dy = matmul(transpose(K(:s,:)), A(s,:s)) * h
            K(s,:) = solver%fn(t + C(s) * h, y + dy)
        end do
        y_new = y + h * matmul(transpose(K(:N_STAGES-1,:)), B)
        f_new = solver%fn(t + h, y_new)
        K(N_STAGES,:) = f_new
    end subroutine
    ! procedures for base class RungeKutta
    subroutine rkinit(self, fun, t0, y0, t_bound, max_step_opt, &
                 rtol_opt, atol_opt)
        ! initialize the solver
        class(RungeKutta),  intent(inout) :: self
        procedure(eomsig)                 :: fun
        real(WP),           intent(in)    :: t0, y0(:), t_bound
        real(WP), optional, intent(in)    :: max_step_opt, &
                                             rtol_opt(:), &
                                             atol_opt(:)
        real(WP)                          :: max_step, rtol(size(y0)), atol(size(y0))

        if (present(max_step_opt)) then
            max_step = max_step_opt
        else
            max_step = huge(max_step_opt)
        end if
        if (present(rtol_opt)) then
            if (size(rtol_opt).eq.size(y0)) then
                rtol = rtol_opt
            else
                rtol = rtol_opt(1)
            endif
        else
            rtol = 1.e-5_wp
        endif
        if (present(atol_opt)) then
            if (size(atol_opt).eq.size(y0)) then
                atol = atol_opt
            else
                atol = atol_opt(1)
            endif
        else
            atol = 1.e-8_wp
        endif

        if (t0 <= t_bound) then
            self%direction=1._WP
        else
            self%direction = -1._WP
        end if
        allocate(self%y, self%y_old, self%rtol, self%atol, self%f, mold=y0)
        self%y = 0._wp
        self%y_old = 0._wp
        self%rtol = 0._wp
        self%atol = 0._wp
        self%f = 0._wp
        self%t_old = nan
        self%t = nan
        self%t_bound = nan
        self%h_previous = 0._wp
        self%error_exponent = 0._wp
        self%max_step = 0._wp
        self%h_abs = 0._wp
        self%n = 0

        self%n = size(y0)
        self%y_old = huge(self%y_old)
        self%y = y0
        self%max_step = validate_max_step(max_step)
        self%rtol = rtol
        self%atol = atol
        self%t = t0
        self%t_bound = t_bound
        call validate_tol(self%rtol, self%atol, self%n)
        self%eom =>fun
        self%f = self%fn(self%t, self%y)
        self%h_abs = select_initial_step( &
            self, t0, y0, self%f, self%direction, &
            self%error_estimator_order, self%rtol, self%atol)
        allocate(self%K(self%n_stages+1, self%n))
        self%K = 0._wp
        if (allocated(self%K_extended)) deallocate(self%K_extended)
        if (allocated(self%K)) deallocate(self%K)
        allocate(self%K_extended(0:N_STAGES_EXTENDED-1, 0:self%n-1), &
               & self%K(0:self%n_stages, 0:self%n-1))
        self%error_exponent = -1._WP / (real(self%error_estimator_order,WP) + 1._WP)
        self%h_previous = huge(self%h_previous)
        self%statmsg = "running"
    end subroutine
    !DOP853 Type Procedures
    function estimate_error_norm(self, K, h, scaleval) result(res)
        class(RungeKutta), intent(in) :: self
        real(WP),      intent(in) :: K(:,:), h, scaleval(self%n)
        real(WP)                  :: res, err5(self%n), err3(self%n), err5_norm_2, err3_norm_2, denom
        err5 = 0._wp
        err3 = 0._wp
        where (abs(scaleval) > 0._WP) 
            err3 =  matmul(transpose(K), self%E3) / scaleval
            err5 =  matmul(transpose(K), self%E5) / scaleval
        endwhere
        err5_norm_2 = sum(err5**2)
        err3_norm_2 = sum(err3**2)
        if (err5_norm_2 == 0._wp .and. err3_norm_2 == 0._wp) then
            res =  0._wp
        end if
        denom = err5_norm_2 + 0.01_WP * err3_norm_2
        res = abs(h) * err5_norm_2 / sqrt(denom * size(scaleval))
    end function

    function step_impl(self) result(res)
        ! Carry out bookkeeping for accepting/rejecting an RK step
        ! including error estimation and step size selection
        class(RungeKutta), intent(inout) :: self
        real(WP)                         :: t, y(self%n), &
                                          & rtol(self%n), &
                                          & atol(self%n), &
                                          & max_step, min_step, & 
                                          & h_abs, h, t_new, y_new(self%n), &
                                          & scaleval(self%n), &
                                          & error_norm, f_new(self%n), &
                                          & factor
        logical                          :: res, step_accepted, step_rejected
        integer nrejected
        t = self%t
        y = self%y

        max_step = self%max_step
        rtol = self%rtol
        atol = self%atol

        min_step = 10._WP *abs(spacing(t))

        if (self%h_abs > max_step) then
            h_abs = max_step
        else if (self%h_abs < min_step) then
            h_abs = min_step
        else
            h_abs = self%h_abs
        end if

        step_accepted = .False.
        step_rejected = .False.

        nrejected = 0
        do while (.not. step_accepted)
            if (h_abs < min_step) then
                print *, self%TOO_SMALL_STEP," ", nrejected
                res =  .False.
                return
            endif

            h = h_abs * self%direction
            t_new = t + h

            if (self%direction * (t_new - self%t_bound) > 0._WP) then
                t_new = self%t_bound
            endif
            h = t_new - t
            h_abs = abs(h)
            call rk_step(self, t, y, self%f, h, self%A, &
                               self%B, self%C, self%K, y_new, f_new)
            ! sk = atol + rtol*max(abs(y(i)), abs(k5(i)))
            scaleval = atol + max(abs(y), abs(y_new)) * rtol
            error_norm = self%estimate_error_norm(self%K, h, scaleval)

            if (error_norm < 1) then
                if (error_norm == 0) then
                    factor = MAX_FACTOR
                else
                    factor = min(MAX_FACTOR, &
                                 SAFETY * error_norm ** self%error_exponent)
                endif

                if (step_rejected) then
                    factor = min(1._WP, factor)
                endif

                h_abs =  h_abs*factor

                step_accepted = .TRUE.
            else
                h_abs = h_abs *max(MIN_FACTOR, &
                             SAFETY * error_norm ** self%error_exponent)
                step_rejected = .TRUE.
                nrejected = nrejected + 1
            endif
        end do

        self%h_previous = h
        self%y_old = y

        self%t = t_new
        self%y = y_new

        self%h_abs = h_abs
        self%f = f_new

        res = .True. 
    end function
    subroutine DOP853DenseOutputImpl(self,tgt) 
        ! Copmute the extra stages required for dense output
        class(RungeKutta), intent(inout)          :: self
        type(DOP853DenseOutput), intent(out)      :: tgt
        real(WP), allocatable        :: testmat(:,:)
        real(WP)                     :: accum(0:size(self%D,1)-1,0:size(self%K,2)-1)
        real(WP)                     :: K(0:size(self%K_extended,1)-1,0:size(self%K_extended,2)-1), h, &
                                      & dy(self%n), F(0:INTERPOLATOR_POWER-1, 0:self%n-1), &
                                      & f_old(self%n), delta_y(self%n)
        integer s, iter2, i, j
        K = 0._WP
        testmat = self%K
        K(0:self%n_stages, 0:self%n-1)= self%K
        h = self%h_previous
        iter2 = 0
        do s= self%n_stages+1, self%N_STAGES_EXTENDED-1
            dy = matmul(transpose(K(:s,:)),self%A_EXTRA(iter2,:s))*h
            K(s,:) = self%fn(self%t_old + self%C_EXTRA(iter2) * h, self%y_old + dy)
            iter2 = iter2 + 1
        end do
        f_old = K(0,:)
        delta_y = self%y - self%y_old

        F(0,:) = delta_y
        F(1,:) = h * f_old - delta_y
        F(2,:) = 2 * delta_y - h * (self%f + f_old)
        do i=0,size(self%D,1)-1
        do j=0,size(self%K,2)-1
            accum(i,j) = dot_product(self%D(i,:),K(:,j))
            end do
        end do
        F(3:,:) = h * accum
        testmat=self%D

        call tgt%init(self%t_old, self%t, self%y_old, F)
    end subroutine


    ! Procedures for DOP853DenseOutput type
    subroutine DOP853denseinit(self, t_old, t, y_old, F)
        ! Initializer
        class(DOP853DenseOutput), intent(inout) :: self
        real(WP), intent(in)                    :: t_old, t, y_old(:), F(0:,0:)
        self%t_old = t_old
        self%t = t
        self%t_min = min(t, t_old)
        self%t_max = max(t, t_old)
        self%h = t-t_old
        allocate(self%y_old,mold=y_old)
        self%y_old = y_old
        self%n = size(y_old)
        allocate(self%Frev(0:size(F,1)-1,0:size(F,2)-1))
        self%Frev = F(size(F,1)-1:0:-1,:)
    end subroutine
    subroutine DOwrite(self, unit_num,prec)
        ! Write to disk
        class(DOP853DenseOutput), intent(inout) :: self
        integer, intent(in)                     :: unit_num
        integer, intent(in), optional           :: prec
        integer                                 :: yos, frs1, frs2
        real(WP)                                :: y_old(size(self%y_old)), &
                                                   Frev(0:size(self%Frev,1)-1,0:size(self%Frev,2)-1)
        integer  :: p
        if(present(prec)) then
            p = prec
        else
            p = WP
        endif
        y_old = self%y_old
        yos = size(self%y_old)
        Frev = self%Frev
        frs1 = size(self%Frev,1)
        frs2 = size(self%Frev,2)
        write(unit_num) yos
        write(unit_num) frs1
        write(unit_num) frs2
        write(unit_num) self%n
        select case (p)
        case (DP)
            write(unit_num) real(self%t_old,DP)
            write(unit_num) real(self%t    ,DP)
            write(unit_num) real(self%t_min,DP)
            write(unit_num) real(self%t_max,DP)
            write(unit_num) real(self%h    ,DP)
            write(unit_num) real(y_old     ,DP)
            write(unit_num) real(Frev      ,DP)
        case (QP)
            write(unit_num) real(self%t_old,QP)
            write(unit_num) real(self%t    ,QP)
            write(unit_num) real(self%t_min,QP)
            write(unit_num) real(self%t_max,QP)
            write(unit_num) real(self%h    ,QP)
            write(unit_num) real(y_old     ,QP)
            write(unit_num) real(Frev      ,QP)
        end select
    end subroutine
    subroutine DOread(self, unit_num)
        ! Read from disk
        class(DOP853DenseOutput), intent(inout) :: self
        integer, intent(in)                     :: unit_num
        integer yos, frs1, frs2
        read(unit_num) yos
        read(unit_num) frs1
        read(unit_num) frs2
        read(unit_num) self%n

        read(unit_num) self%t_old
        read(unit_num) self%t 
        read(unit_num) self%t_min
        read(unit_num) self%t_max
        read(unit_num) self%h
        if (allocated(self%y_old)) deallocate(self%y_old)
        if (allocated(self%Frev)) deallocate(self%Frev)
        allocate(self%y_old(yos), self%Frev(0:frs1-1,0:frs2-1))
        read(unit_num) self%y_old
        read(unit_num) self%Frev
    end subroutine

    function DOP853densecall(self, t) result(res)
        ! compute the value of the stored interpolant
        class(DOP853DenseOutput), intent(in) :: self
        real(WP),                 intent(in)    :: t(:)
        real(WP)                                :: res(self%n,size(t)), &
                                                 & x(size(t)), &
                                                 & y(self%n,size(t))
        integer i, j
        x = (t - self%t_old) / self%h

        y = 0._WP
        do i=0,size(self%Frev,1)-1
            forall (j=1:size(t)) y(:,j) = y(:,j) + self%Frev(i,:)
            if (mod(i,2)== 0) then
                forall (j=1:size(t)) y(:,j) = y(:,j)*x(j)
            else
                forall (j=1:size(t)) y(:,j) = y(:,j)*(1._WP - x(j))
            endif
        end do
        !$OMP PARALLEL DO
        do i=1,size(t)
            y(:,i) = y(:,i) + self%y_old
        end do
        !$OMP END PARALLEL DO
        res = y
    end function
    function scall(self, t) result(res)
        ! Call the stored interpolant once
        class(DOP853DenseOutput), intent(inout) :: self
        real(WP),           intent(in) :: t
        real(WP) :: resarray(self%n,1), res(self%n)
        resarray = self%call_impl([t])
        res = resarray(:,1)
    end function scall
    function mcall(self, t) result(res)
        ! Call the stored interpolant at multiple values
        class(DOP853DenseOutput), intent(inout) :: self
        real(WP),           intent(in) :: t(:)
        real(WP)                       :: res(self%n,size(t))
        res = self%call_impl(t)
    end function mcall

    function select_initial_step(solver, t0, y0, f0, direction, order, rtol, atol) result(res)
        ! Empirically select a good initial step.

        ! Parameters
        ! ----------
        ! solver : RungeKutta type
        ! t0 : real
        !     Initial value of the independent variable.
        ! y0 : real (n)
        !     Initial state
        ! f0 : real (n)
        !     Initial value of the derivative, i.e., ``fun(t0, y0)``.
        ! direction : real
        !     Integration direction.
        ! order : integer
        !     Error estimator order. It means that the error controlled by the
        !     algorithm is proportional to ``step_size ** (order + 1)`.
        ! rtol : real
        !     Desired relative tolerance.
        ! atol : real
        !     Desired absolute tolerance.
        ! Returns
        ! -------
        ! h_abs : real
        !     Absolute value of the suggested initial step.
        class(RungeKutta)       :: solver
        real(WP), intent(in)    :: t0, y0(:), f0(:), direction, rtol(:), atol(:)
        integer                 :: order
        real(WP)                :: res
        real(WP)                :: d0, d1, d2, h0, h1, y1(size(y0)), &
                                 & f1(size(y0)), normvec(size(y0)), &
                                 & sc(size(y0))
        sc = atol + abs(y0) * rtol
        normvec = y0/sc
        d0 = norm(normvec)
        normvec = f0/sc
        d1 = norm(normvec)
        if (d0 < 1.e-5_WP .or. d1 < 1.e-5_WP) then
            h0 = 1.e-6_WP
        else
            h0 = 0.01_WP * d0 / d1
        endif
        y1 = y0 + h0 * direction * f0
        f1 = solver%fn(t0 + h0 * direction, y1)
        normvec = (f1-f0)/sc
        d2 = norm(normvec) / h0

        if (d1 <= 1.e-15_WP .and. d2 <= 1.e-15_WP) then
            h1 = max(1.e-6_WP, h0 * 1.e-3_WP)
        else
            h1 = (0.01_WP / max(d1, d2)) ** (1._WP / (real(order,WP) + 1._WP))
        endif 
        res = min(100 * h0, h1)
    end function
    function solve_ivp(fun, t_span, y0, &
                       t_eval_in, dense_output, &
                       rtol, atol, istep) result(endsol)
        procedure(eomsig)                      :: fun
        real(WP),         intent(in)           :: t_span(2), &
                                                & y0(:)
        real(WP),         intent(in), optional :: t_eval_in(:)
        real(WP),         intent(in), optional :: rtol, atol, istep
        logical,          intent(in), optional :: dense_output
        character(20)                          :: message
        logical                                :: dovar, loud
        real(WP)                               :: t0, tf, &
                                                & t_old, &
                                                & t, &
                                                & y(size(y0)), &
                                                & rtolact, atolact
        integer                                :: t_eval_i, &
                                                & solverstat
        type(ExtensibleRealArray)              :: ts, ys, ti
        real(WP),                  allocatable :: t_eval(:), &
                                                & ts_final(:), &
                                                & ys_final(:,:)
        type(RungeKutta)                           :: solver
        type(DOP853DenseOutput)                :: sol
        type(ExtensibleDOArray)       :: interpolants
        type(DOP853DenseOutput),        allocatable  :: interpolants_final(:)
        type(Odesolution)                      :: endsol
        if (present(dense_output)) then
            dovar=dense_output
        else
            dovar=.false.
        endif
        loud = .true.
        if (present(rtol)) then
            rtolact = rtol
        else
            rtolact=1.e-8_WP
        endif

        if (present(atol)) then
            atolact = atol
        else
            atolact=1.e-10_WP
        endif

        t0 = t_span(1); tf = t_span(2)

        if (present(t_eval_in)) then
            if (tf > t0) then
                t_eval = t_eval_in
            else
                ! Make order of t_eval increasing to use np.searchsorted.
                t_eval = t_eval(size(t_eval_in):1:-1)
                ! This will be an upper bound for slices.
                t_eval = size(t_eval_in)
            endif
            t_eval_i = 0
        endif

        call solver%init(fun, t0, y0, tf, &
                         rtol_opt=[rtolact], &
                         atol_opt=[atolact])
        ! Allow for manual specification of initial step
        if (present(istep)) then
            solver%h_abs = istep
        endif

        if (.not. present(t_eval_in)) then
            call ts%create(t0)
            call ys%create(y0)
            if (dovar) call interpolants%create()
        else if (present(t_eval_in) .and. dovar) then
            call ts%create(t0,moldonly=.true.)
            call ti%create(t0)
            call ys%create(y0,moldonly=.true.)
            call interpolants%create()
        else
            call ts%create(t0,moldonly=.true.)
            call ys%create(y0,moldonly=.true.)
        endif
        solverstat = -50
        do while (solverstat.eq.-50)
            message = solver%step()

            if (trim(solver%statmsg) == 'finished') then
                solverstat = 0
            else if (trim(solver%statmsg) == 'failed') then
                solverstat = 1
                if (loud) then
                else
                    exit
                endif
            end if 
            t_old = solver%t_old
            t = solver%t
            y = solver%y

            if (dovar) then
                call solver%dense_output(sol)
                call interpolants%append(sol)
            end if
            if (.not.present(t_eval_in)) then
                call ts%append([t])
                call ys%append(y)
            else
                if (.not.dovar) then
                    call solver%dense_output(sol)
                end if
                do while (t_eval(t_eval_i) .le. solver%t)
                    call ts%append(t_eval(t_eval_i))
                    call ys%append(sol%dcall(t_eval(t_eval_i)))
                    t_eval_i = t_eval_i + 1
                end do
            end if

            if (present(t_eval_in).and.dovar) then 
                call ti%append([t])
            end if
        end do

        message = MESSAGES(solverstat)

        call ts%seal(ts_final)
        call ys%seal(ys_final)

        if (dovar) then
            if (.not.present(t_eval_in)) then
                call interpolants%seal(interpolants_final)
                call endsol%init(ts_final, ys_final, solver%nfev, message, interpolants_final)
            else
                call interpolants%seal(interpolants_final)
                call endsol%init(ts_final, ys_final, solver%nfev, message, interpolants_final)
            end if 
        else 
            call endsol%init(ts_final, ys_final, solver%nfev, message)
        endif
    end function
    ! utility functions 
    function validate_first_step(first_step, t0, tf) result(res)
        ! Assert that first_step is valid and return it.
        real(WP), intent(in) :: first_step, t0, tf
        real(WP)             :: res
        if (first_step <= 0._WP) then
         print *, "`first_step` must be positive."
         error stop
        endif 
        if (first_step > abs(tf - t0)) then
         print *, "`first_step` exceeds bounds."
         error stop
        endif
        res = first_step
    end function
    function validate_max_step(max_step) result(res)
        real(WP), intent(in) :: max_step
        real(WP)             :: res
        !Assert that max_Step is valid and return it.
        if (max_step <= 0._WP) then
         print*, "`max_step` must be positive."
         error stop
        endif
        res = max_step
    end function
    subroutine validate_tol(rtol, atol, n)
        !Validate tolerance values
        real(WP), intent(inout) :: rtol(:), atol(:)
        integer, intent(in)     :: n
        integer i
        if (any((rtol < 100._WP * de))) then
         print *, "At least one element of rtol is too small. Setting small elements to ", 100._WP * de
         forall (i=1:n) rtol(i) = max(rtol(i), 100 * de)
        endif
        if (size(atol) > 0 .and. size(atol).ne.n) then
         print *, "`atol` has wrong shape."
         error stop
        endif 
        if (any(atol < 0._WP)) then
         print *, "`atol` must be positive."
         error stop
        endif 
    end subroutine
    function norm(x) result(res)
        real(WP), intent(in) :: x(:)
        real(WP)             :: res
        ! Compute RMS norm.
        res = norm2(x) / size(x) ** 0.5_WP
    end function
    function diff(x) result(res)
        real(WP), intent(in) :: x(:)
        real(WP)             :: res(size(x)-1)
        integer i
        !$OMP PARALLEL DO
        do i=1,size(x)-1
         res(i) = x(i+1) - x(i)
        end do
        !$OMP END PARALLEL DO
    end function
    function binsearch(array,val) result(res)
        real(WP), intent(in) :: array(:), val
        real(WP)             :: workarray(size(array))
        logical              :: flipped
        integer :: j, top, bot, n, res
        n = size(array)
        if (array(1).lt.array(n)) then
         workarray = array
         flipped = .false.
        else
         workarray = array(n:1:-1)
         flipped = .true.
        end if
        res = -1
        bot= 1
        top= n
        j = (bot+top)/2
        do while (bot.lt.top)
         if (val.lt.workarray(j)) then
             top= j
         else if (workarray(j+1).lt.val) then
             bot= j+1
         else  
             res = j 
             if (val.eq.workarray(j+1)) res = j+1;
             if (flipped) res = n-res + 1
             exit
         end if 
         j = (bot+ top)/2
        end do
        if (bot==n) res = n-1
    end function binsearch
end module frkmin_q
