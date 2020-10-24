%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%                                                       %%
 %%      This script runs  FullNyeFowlerForEnviro_1.m     %%
 %%                                                       %%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script is desinged to run the function FullNyeFowlerForEnviro_1.m
% while varying crucial parameters like the water supply to the channel MR,
% the input to the lake Qin and the basic hydraulic gradient parameter a.
%
% As it is currently setup defining MRVarying as a vector will cause the
% model to carry out multiple simulations using each value of MR. Uncommenting the 
% save call at the bottom of the script enables you to save the results of each simulation.  
%
% The following parameters recreates one simulation from the parameter
% search shown in Figure 3.7 in Chapter 3 of Kingslake (2013)
% http://etheses.whiterose.ac.uk/4630/1/J_Kingslake_ThesisFINAL.pdf:
% MRdim = 0.0005 m^2 s^-1, a = 8, Qin_dim = 5 m^3 s^-1, using the numerical
% method 'Finite'.
%
% See section 2.2.5 of Kingslake (2013) for a decription of the numercial
% methods. 
%
% THe parameter RunInfo.plots allows the user to choose which plots they
% want to see as the simulation runs. Define it as an array with either
% none, one or multiple entries:
%    
%      Entry     Plot
%      1         Along-channel profile of discharge, Q(s) 
%      2         Along-channel profile of effective pressure, N(s) 
%      3         Time series of discharge at the lake, Q(0,t)
%      4         Along-channel profile of channel cross-sectional area, S(s)
%      5         Time series of effective pressure at the terminus, N(s0,t) 
%      6         Time series of the minumum area of the channel,    min(S(s))(t)
%      7         Along-channel profile of disharge and effective pressure together
%      8         Phase-space plot of lake depth vs discharge at the lake
%      9         Rate of change of discharge at the lake, dQR/dt
%      10        Along-channel profile of disharge, effective pressure and discharge together
%      11        Rate of change of discharge with lake height, dQR/dh 
%      12        Along-channel profile of channel radius at the lake, hr(s) 
%      13        Time series of lake input, Qin(t) 
% 
% e.g. if you want to plot Q(s,t) and N(s0,t) use:
% RunInfo.plots    = [1 5];

%
% Written by J. Kingslake, Department of Geography, University of Sheffield
% Modified by A. Jenson and J. Amundson, University of Alaska Southeast

global rho_i rho_w g
rho_i = 917; % ice density [kg m^-3]
rho_w = 1000; % water density [km m^-3]
g = 9.81; % gravitational acceleration [m s^-2]

% This model assumes that the glacier flow model has already been run.
% Here, specify the year(s) of the glacier geometry to use in the outburst
% flood simulations, keeping in mind that terminus retreat doesn't start
% until year 100.
t_years = [6:2:400];
%t_years = [12:12:400];
%t_years = 200;
%t_years = 006;

% Choose which slope, climate, basin location, basin-shape, and basin with/without ice

%slope = 2;
slope = 4;
%slope = 6;

%climate ='maritime';
climate = 'continental';

basin = 'ice';
%basin = 'noice';
%basin = 'ice_thickness_defined' ; 

%pL = 1; %box-shaped basin;
pL = 2 ; %wedge-shaped basin
%pL = 3 ;    % cone-shaped basin

s0 = 0.75; % location of basin in terms of percentage of glacier from the divide
   
for i = 1:length(t_years)

% load glacier model for specified years; note that the location of the
% model output is hard-coded and may vary for your computer and for
% different simulations...
yearlabel = num2str(round(t_years(i)), '%03.f');
%load(['./output.' yearlabel '.mat'])
load(['glacier_model/output/ELA1500/RCP8.5/slope' num2str(slope) '/' climate '/s0=' num2str(s0) '/output.' yearlabel '.mat'])


% compute basin geometry for each time step

if pL ==1
    basin_shape = 'box-shaped'; 
    a = 8.5e5; % a = W*L
end
    
if pL == 2
    W = 1000; %width of basin in meters    
    theta_b = 15; % basin slope [deg]
    basin_shape = 'wedge-shaped';
    a = W*cotd(theta_b); %for a wedge shaped basin with a 15 degree slope
end

if pL == 3
    theta_b = 19.75;
    basin_shape = 'cone-shaped';
    hLi = rho_i/rho_w*(sol.H_basin(end));
    a = pi/2*(cotd(19.75))^2;% for a cone-shaped basin
end


basin_geometry = LakeWithIce_with_GC(sol,a,pL,basin); % returns initial water level and volume, depends on if there is ice in the basin

% specify run parameters
RunInfo.DateAndTime            = datestr(now);
RunInfo.Qin_dim                = 10;  % Input to the lake m^3 s^-^1
RunInfo.InitialLakeDepthDim    = (rho_i/rho_w*basin_geometry.H_basin)-rho_i/rho_w*basin_geometry.h_ice;  % Inital lake level
RunInfo.plots                  = [3];  % choose what plots to display (see table above)
RunInfo.PlotPeriod             = 10;   % The interval between plots in time steps
RunInfo.a                      = 0;    % this is the basic hydraulic gradient (psi) parameter, when a>1, psi < 0 at the lake, larger a --> more -ve at lake and region of -ve psi extends downglacier see 
RunInfo.Notes                  = 'None';
RunInfo.system                 = 'Typical';   % this can be ignored
RunInfo.InitialSGuess          = 1;  % in m^2
RunInfo.DPing                  = 0.00 ;% you can try adjusting this numerical (\mu in equation 2.74 of Kingslake, 2013) it can speed up convergence when using 'Relax'.
RunInfo.PlotRelax              = 0;    % an option to plot effective pressure and discharge profiles as they 'Relax' to a solution using the Relaxation method
% This is where you choose the numerical method to use. 'Relax' is the so-called Relaxation
% method described in section 2.2.5 of Kingslake (2013). It solves the full version of the model.
% 'Finite' solves a reduced version of the equations (see Folwer, 1999)
% also described in section 2.2.5 of Kingslake (2013). The advantage of
% 'finite' is that it is much fster than 'Relax'  

% Note: Because the equations solved and the boundary conditions used in each
% case are different, results obtained using the two methods are
% quantitatively different. 
RunInfo.method                 = 'Relax';           
%RunInfo.method                 = 'Finite';
RunInfo.MRdim                  = 0; %10^-5; %input into channel along length

modelParams.RunInfo = RunInfo; %specific run parameters

output = FullNyeFowlerForEnviro_updated(RunInfo,basin_geometry,glacier,pL);

% add time and basin geometry at that time to the model parameters array; this is for saving purposes 
modelParams.t_years = t_years(i);
modelParams.basin_geometry = basin_geometry;
modelParams.RunInfo;

mkdir(['outburst_flood_model/FullModel/ELA1500/RCP8.5/slope' num2str(slope) '/' climate '/s0=' num2str(s0) '/' basin_shape '/' basin '/hL=floatation' ]);           
save(['outburst_flood_model/FullModel/ELA1500/RCP8.5/slope' num2str(slope) '/' climate '/s0=' num2str(s0) '/' basin_shape '/' basin '/hL=floatation/FullModelt_years=' num2str(t_years(i),'%03.f') '___.mat'],'output','modelParams');
%mkdir(['outburst_flood_model/FullModel/ELA1500/RCP8.5/slope' num2str(slope) '/' climate '/s0=' num2str(s0) '/' basin_shape '/' basin ]);           
%save(['outburst_flood_model/FullModel/ELA1500/RCP8.5/slope' num2str(slope) '/' climate '/s0=' num2str(s0) '/' basin_shape '/' basin '/' 'h_ice=' num2str(basin_geometry.h_ice,'%03.f') '.mat'],'output','modelParams');
end