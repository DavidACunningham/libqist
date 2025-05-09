% Matlab example driver for QIST
% Author: David Cunningham
% Date: 1 May 2025
% currently must be run from folder /home/david/wrk/nstgro/qist
% See comments inline for usage info

% Environment setup
clear all;
close all;
clc;
qistdir="../fort/lib/";																 % This works if you're running the script from the $LIBQIST/examples folder
cspice_furnsh('../kernels/mk_example_with_traj.tf')
mu_d = 0.0000985; % deimos
nmlstring = "./curve_deimos_config_namelist_2026Nov26120000002026Nov2700000000.nml"; % the QIST configuration namelist
charnml = convertStringsToChars(nmlstring);											 % all text passed to the Fortran library must be of character type
loadlibrary(strcat(qistdir,"mqist.a"),strcat(qistdir,"mqist.h"));					 % mqist.a must be on your system.
calllib('mqist','m_init_n', charnml);

% Times: These initial and final values are the ones set in the configuration namelist.
t0_utc = convertStringsToChars("2026 Nov 26 12:00:00.00");
tf_utc = convertStringsToChars("2026 Nov 27 00:00:00.00");
% QIST, like MICE, computes times in seconds past J2000.
t0  =  cspice_str2et(t0_utc);
tf  =  cspice_str2et(tf_utc);

% End of setup section


%%%%% REFERENCE TRAJECTORY %%%%%
% Now we will retrieve the reference trajectory from the kernel using MICE 
% purely so we can plot it. You don't have to use MICE when you use the QIST runtime
% but it can be nice to have around.

times = linspace(t0,tf,1000);
states = zeros(6,1000);
for i =1:1000
    [this_state,~] = cspice_spkgeo(-31415, times(i), 'J2000', 402);
	states(:,i) = this_state;
end

%%%%% PROPAGATING THE RELATIVE TRAJECTORY WITH QIST %%%%%
% We use QIST to propagate the relative trajectory. The specified initial conditions should
% result in a spiral when viewed in the reference-centered IJK frame.
% QIST natively propagates in the J2000 frame.
dx0 = [0, 0.1, 0., 0.1*sqrt(mu_d/10.^3)/2., 0., 0.]';
dxs = zeros(6,1000);
for i =1:1000
	% Initialize the output state to NaN
	this_state = NaN(8,1);

	% Here is the actual call to the QIST library, specifically the "prop_once" function.
	% Note that the physical state is the first six elements of the initial vector, in this case
	% [dx0; 0.; 0.]. the last two elements are for a change in t0 and a change in TOF and should
	% almost always be zero. These are % supplied so that free-time optimization can be conducted 
	% using QIST. A MATLAB convenience wrapper for QIST could hide this, but is not performed here.
	[~,~,~,~,this_state] = calllib('mqist','m_prop_once', t0, times(i), [dx0;0.;0.], 2, this_state);

	% Assign the output state to the actual relative state variable
	dxs(:,i) = this_state(1:6);
end

ref_fig=figure('Units','inches', 'Position',[0.5,2, 6,6]);
ref_state_plot = plot3(states(1,:), states(2,:), states(3,:));
hold on;
ref_state_plot.LineWidth=2;
ellipsoid(0,0,0,8.04,5.89,5.11, 'facecolor','#D3D3D3','edgecolor','none');
grid on;
axis equal;
xlabel("X (km)");
ylabel("Y (km)");
zlabel("Z (km)");
title("Deimos with Reference Trajectory");
rel_fig=figure('Units','inches', 'Position',[0.6,4, 6,6]);

rel_state_plot = plot3(dxs(1,:),dxs(2,:),dxs(3,:));
hold on;
rel_state_plot.LineWidth=2;
ref_point = scatter3(0,0,0,'x');
grid on;
axis equal;
xlabel("\delta x (km)");
ylabel("\delta y (km)");
zlabel("\delta z (km)");
legend("Relative Trajectory", "Ref")
title("Relative Trajectory");
