function output = FullNyeFowlerForEnviro_updated(RunInfo, basin_geometry, glacier, pL, varargin)


%function [TSamp,QRLakeSamp,xTposSamp,hrLakeSamp,hLSamp,VLSamp,QRendSamp,meanNRSamp,TDays,t0,QR0,s0,hL0,VL0,dt,ds,MRdim,s,exitflag,QRSamp,SRSamp,SR0,Peak,PeakTimeSecs,PeakYearFrac,LakeEmptied,ChannelClosed,Highstand,Lowstand] ...
%    = FullNyeFowlerForEnviro_updated(RunInfo,basin_geometry, varargin)

%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                                          %%
%%         This function solves equations describing channelised subglacial water flow      %%
%%            coupled to a marginal lake in terms of both pressure and discharge.           %%
%%                                                                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
% Function can be called by RunningNyeFowlerCode.m  
%
% If using the relaxation method it requires NewRaphToFindQtop3.m
%
% 
%
% exitflag = 1       model reached the end of the run succesfully
% exitflag = -1      relaxation method failed at some point in the run
% exitflag = NaN     model didnt reach t=T but dont know what happened, exitflag not defined again after model setupWholeDays = NaN;
%
% Written by J. Kingslake, 2011, Department of Geography, University of
% Sheffield, UK.
% Modified by A. Jenson and J. Amundson, University of Alaska Southeast

global rho_i rho_w g

% parse model input
Qin_dim = RunInfo.Qin_dim; % input to the lake
InitialLakeDepthDim = RunInfo.InitialLakeDepthDim; % initial lake depth, set so that ice at bottom of the dam is at flotation
plots = RunInfo.plots; % just specify which plots to display while this is running
PlotFreq = RunInfo.PlotPeriod; % how frequently to plot model output while the simulation is running
system = RunInfo.system; %???
%InitialSGuess = RunInfo.InitialSGuess; % initial channel area
a = RunInfo.a; % basic hydraulic gradient, psi???
DPing = RunInfo.DPing;
PlotRelax = RunInfo.PlotRelax;
method = RunInfo.method;
MRdim = RunInfo.MRdim;



Q_old = NaN;
hL_old = NaN;
P = 1;
Llo = 1;
Lhi = 1;

% if ishandle(999) ==0; 
%     figure(999);
%     text(0.3,0.5,'Control-click to end simulation','FontSize',18)
%     set(999,'WindowButtonDownFcn',@EndRun);
% end

global UserReturn
UserReturn = 0;


LakeEmptied = 0;
ChannelClosed = 0;
Peak = NaN;
PeakTimeSecs = NaN;
PeakYearFrac = NaN;
Highstand = NaN;
Lowstand = NaN;

if size(varargin,1)==1      % if a values has been entered for Ice thickness assign it
    IceThickness = varargin{1};
elseif size(varargin,1)==0 % if nothing was entered default to 100m
    IceThickness = basin_geometry.hLi*rho_w/rho_i; 
end


global paused
paused=0;
exitflag = NaN;

%%%%%% MODEL PHYSICAL PARAMETERS %%%%%%%%%%%%
% R model
L = 3.34e5;                                          % latent heat of fusion of water [J kg^-1]
ni = 0.0600;                                        % roughness of ice wall [m^-1/3 s]
nb = 0.1629;                                        % roughness of bed material [m^-1/3 s]
n_full = ni * 0.611 + (1 - 0.611) * nb;             % roughness of channel when full [m^-1/3 s]
f = 6.6 * n_full^2;                                 % friction factor from Fowler(1999) when channel is full See page 80 of book 6 of notes
sqg = sqrt(g);                                      % square root of gravity
K = (2.4e-24)*2/9; %10^-24;                         % ice flow constant
Slope = -(glacier.hb(end)-glacier.hb(1))/glacier.x(end); % conduit/glacier bed slope [rad]
n = 3;                                              % Glen's flow law exponent


% lake reference geometry, assumed at lake full
hLi = basin_geometry.hLi;    % lake reference depth (when ice dam is at hydrostatic equilibrium)
VLi = basin_geometry.VLi;    % lake reference volume (when ice dam is at hydrostatic equilibrium)   

IceDamHeight = basin_geometry.H_basin;

% scales for non-dimensionalization
s0 = (glacier.x(end) - glacier.x_basin)/cos(Slope); % length of channel [m] (length scale)

% lake depth and volume scales (use lake depth and volume at flotation)
%hL0 = hLi;
%VL0 = VLi;
hL0 = basin_geometry.H_basin*rho_i/rho_w;
VL0 = (hL0/hLi)^pL*VLi;

QR0   = 1000;                                   % discharge scale [m^3 s^-1]
psi0 = rho_w*g*sin(Slope) + rho_i*g*IceDamHeight/s0;            % basic hydraulic gradient scale [Pa m^-1]; note that we are assuming that this is constant
SR0   = (f*rho_w*g*QR0^2/psi0)^(3/8);          % area of channel scale [m^2]
m0   = psi0*QR0/L;                             % melting of walls scale
t0   = rho_i*SR0/m0;                           % time scale [s]
N0   = (K * t0) ^(-1/3);                       % effective pressure scale [NR m^-2]
hr0  = sqrt(2*SR0/pi);                         % channel roof height scale [m]
%hw0  = sqrt(2*SR0/pi);                         % water in channel depth scale [m]
%A0   = (pi/2)*hr0^2;                           % Area of flow scale [m^2]
%P0   = hr0*(pi+2);                             % wetted perimeter scale [m]
%n0   = (2+4/pi)^-(2/3) * hr0^(1/6) * sqrt(Slope/(2*g)); %roughness scale [m^-1/3 s]
%Qin0 = QR0; % lake input scale the same as discharge scale.

InitialSGuess = RunInfo.InitialSGuess/SR0; %make sure that initialSguess is the same for all floods

% dimensionless parameters of model
% F model

epsilonR = s0 * psi0 /rho_i / L;
lambda = t0 * hLi^pL * QR0/(   pL * VLi * hL0^pL);
delta = N0/(s0*psi0);
% kapR = N0 * s0 * delta * con / QR0;
beta = rho_w*g*hL0/N0;
r = rho_i/rho_w;


%%%%%%%%%%%%%%% SET UP TIME %%%%%%%%%%%%
dt=0.01;                                          % time step default 0.01
T = 500;                                              % dimensionless simulation time default 500
t=0:dt:T;                                           % time vector
Lt = length(t);                                     % number of time steps
TSamp = t;%(1:10:Lt);                                 % set up a sampling time vector
TDays = TSamp*t0/3600/24;
tDays = t*t0/3600/24;

% set up space grid
S_end = 1;
ds = 0.01;                                            % space step    usually 0.01
s = 0:ds:S_end;                                         % space vector
Ls = length(s);                                     % number of space steps
% 1 = ceil(Ls/2);                                   % position of start of cavities
% Cpos = 1;
%%%%%%%%%%% PREALLOCATE ARRAYS %%%%%%%%
disp('Pre-allocating Sampling Arrays....')

hLSamp = NaN(length(TSamp),1);                    % array for sampling values of lake depth [m]
VLSamp = NaN(length(TSamp),1);
hrLakeSamp = NaN(length(TSamp),1);                % array for sampling values of channel roof height [m]
% hwLakeSamp= NaN(length(TSamp),1);               % array for sampling values of water flow depth    [m]
xTposSamp = NaN(length(TSamp),1);                 % array for sampling values of transition position [x0] [m]
QRLakeSamp = NaN(length(TSamp),1);                % array for sampling values of the discharge at the lake [QR0]
QRendSamp  = NaN(length(TSamp),1);                % array for sampling values of the discharge at the end of the tunnel [QR0]
ubSamp = NaN(length(TSamp),Ls);                     % array for sampling values of the sliding velocity profile
QRSamp = NaN(length(TSamp),Ls);                     % array for sampling values of the channel discharge profile
SRSamp = NaN(length(TSamp),Ls);                     % array for sampling values of the channel area profile
meanNRSamp = NaN(length(TSamp),1);                % array for sampling mean values of NR
disp('Done')

disp('Pre-allocating Variable Arrays....')
% Lake variables
hL  = zeros(2,1);                                   % lake level  [m]
VL = zeros(2,1);                                    % lake volume

% R Channel variables
hr = zeros(2,Ls);                                   % channel roof height [m]
SR = zeros(2,Ls);                                   % channel cross-section
NR = zeros(2,Ls);                                   % channel effective pressure
QR  = zeros(2,Ls);                                  % discharge in channel  [m^3 s^-1]
% Xs = zeros(2,1);                                  % divide position


% Newtons method arrays
aa = zeros(2,2,Ls);
bb = zeros(2,2,Ls);
cc = zeros(2,2,Ls);
RR = zeros(2,1,Ls);


% Bed supply variable
% omega = zeros(Lt,1);
% MRdim = zeros(Lt,1);
% MCdim = zeros(Lt,1);
disp('Done')


%geometry
%hb = zeros(1,Ls);                                   % glacier bed profile [m]
%hs = zeros(1,Ls);                                   % glacier surface profile [m]
%H = zeros(1,Ls);
% disp('Getting Inputs from GMInputs1.m....')

% non-dimensioalize the lake input
Qin = Qin_dim/QR0.*ones(Lt,1);                                            
disp('Done')

% Boundary condition on N at termius (top boundary condition is supplied by
% the lake level.
% NbottomVariable = 0.1*sin(2*pi*t*t0/6/3600)    % Allows you to impose a time varying pressure boundary condition at the terminus. 
NbottomVariable = zeros(length(t),1);


% non-dimensional input to channel
MR = MRdim*s0 /QR0.*ones(Lt,1) ;                     % R channel input
% MRVariable = MRdimVariable*s0 /QR0 ;     % variable R channel input

% MR(1:21) = [0 1:20]./20.*MR(end);    % This ramps up the input to the channel over the first 20 time time steps because the inital guess at S is uniform and this is unrealistic, so this allows a bit of time for S to evolve towards a  realistic profile while MR is ramped up - maybe not needed.

%%%%%%%%%%%%%%%%%%%%%%%%
%%% Glacier Geometry %%%
%%%%%%%%%%%%%%%%%%%%%%%%
% thickness gradient along the entire length of the glacier
dHdx = gradient(glacier.H,glacier.x);

% interpolate to find the thickness gradient along the length of the
% conduit, which has Ls grid points
dHdx = interp1(glacier.x,dHdx,linspace(glacier.x_basin,glacier.x(end),Ls));

% calculate and nondimensionalize psi
%%% define hydraulic potential
psi = rho_w*g*sin(Slope) - rho_i*g*dHdx;
b = 20;
psi = (psi/psi0).*(1-a*exp(-b*s)); % idealised hydraulic gradient 

% Sampling Counter
SampNumber = 2;
SpinUpFin = 1;
% Qin = QinAfterSpinUp/QR0;
% Qin(:) = 0;

%disp('Done')

%%%%%%% initial conditions %%%%%%%%%

disp('Starting Initial Conditions...')

SR_temp(1,:) = ones(1,Ls).*InitialSGuess;   % initial channel radius as function of position
hL(1,1) = InitialLakeDepthDim/hL0;          % set initial water depth equal to 1
VL(1,1) = VLi*(hL(1,1)/hLi)^pL;             % initial volume in scaled units
NL = beta*(1-hL(1,1))-rho_i*g*basin_geometry.h_ice/N0;    % initial effective pressure (should equal 0 if starting at flotation)
Nbottom = 0;                                % effective pressure at the end of the glacier
% DELETE? NTopFixed = 0.2;
% Define N at top end
NR(1,1) = NL;
% define the end boundary condition on NR which will not change
NR(1,Ls) = Nbottom;

% Use a made up function to define a smooth NR profile which matchs these
% two end conditions
% NR_temp(1,:) = sin((Ls-[1:Ls])./(Ls-1).*(pi-asin(NR(1,1)))) +[(s(:)-s(1))./(1-s(1))*NR(1,Ls)]';
NR_temp = (1-s./1).*NR(1,1); % assume that effective pressure varies linearly along the length of the channel
Xs = NaN(2,1); %???


QR_temp = SR_temp(1,:).^(4/3); % channel discharge depends on channel radius to 4/3 power?


switch method
    case 'Finite'
        % (3) initial conditions

        % S at far end of channel fixes Q there
        QR_temp(1,Ls) = sqrt((SR_temp(1,Ls)^(8/3))*psi(Ls));
        % fixes Q for the whole channel
        QR_temp(1,:) = QR_temp(1,Ls) - MR(1)*(1-s./S_end);

        for k = 1:Ls-1
            NR_temp(1,k+1) = NR_temp(1,k) + ds * ( QR_temp(1,k)*abs(QR_temp(1,k))/(SR_temp(1,k).^(8/3)) - psi(k) );
        end
    case 'Relax'

        % Use the function NewRaphToFindQtop3.m to get the effective pressure at the         InitialQLakeGuess.*ones(1,Ls)
        % lake a guess a the discharge there
        [QR_temp,NR_temp,exitflagfzero,exitflagRelax] = NewRaphToFindQtop3(Nbottom,SR_temp(1,:),psi,ds,s,epsilonR,r,delta,MR(1),NL,QR_temp,DPing,PlotRelax);


        if exitflagfzero ~=1 || exitflagRelax ~=1
            warning 'Model Reached a point where the Relaxation Method Failed'
            exitflag = -1;
            save ('OuterSaveOfFail')
            return
        end

    otherwise
        error 'Unknown method at initial conditions'
end
disp('...')
% SR_temp = QR_temp.^(3/4);
% end

% redefine channel properties (in dimensionless units?)
QR(1,:) = QR_temp; % discharge
NR(1,:) = NR_temp; % effective pressure
SR(1,:) = SR_temp; % channel area



%try loading the post spin-up conditions
% SpinUpFin=1;
% load 'TempCondtitions1'
% Qin(:) = 50/QR0;
disp('Done')


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% MAIN LOOP %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Starting Main Loop...')

output.hL = [];
for i = 2:Lt


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Step Channel Crossection Forward %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update the R channel cross-sectional area
    rhs = abs(QR(1,:).^3) ./ (SR(1,:).^(8/3)) - SR(1,:).*(NR(1,:)).^3; % equation 2.59 in Kingslake thesis
    SR(2,:) = SR(1,:)+ dt * rhs; % equation 2.59 in Kingslake thesis

    if any(SR(2,:)<=0)
        disp 'SR has gone to zero - reducing the time step or increasing the lake input Qin can help to prevent this.'
        
        ChannelClosed = 1;
        return
    end
    
    % calculate channel height as function of distance (assuming half
    % circle cross-section?)
    hr(2,:) = sqrt(SR(2,:)*SR0); % [m]???

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Step Lake and lake effective pressure forward %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        hL(2,1) = InitialLakeDepthDim/hL0;
    hL(2,1) = max(0,hL(1,1) + dt * lambda/hL(1,1)^(pL-1) * (Qin(i)- QR(1,1) )); % equation 2-63 in Kingslake thesis
    VL(2,1) = VLi/VL0*(hL0*hL(2,1)/hLi)^pL; 
 
 
    if hL(2,1)*hL0 <= hr(2,:)
        
       disp 'Lake level dropped below the height of the channel roof! - This version of the model cannot deal with this. A later version described in Ch. 6 of Kingslake (2013) [thesis] implements an open-channel model to simulate what happens when the lake level drops below the level of the channel roof. '
       LakeEmptied = 1;
       return
                 error 'Lake Emptied!!!!'
    end 
    
    % place maximum water level constraint
    if hL(2,1)*hL0 > basin_geometry.H_basin
        %disp 'Overflowed basin'
        %return
       hL(2,1) = basin_geometry.H_basin/hL0;
    end
    
    %if  hL(2,1) == hL(1,1)
 %    if  abs(QR(2,1) - QR(1,1))*QR0 <= 0.15 && QR(2,1) < QR(1,1) ||  hL(2,1) == hL(1,1)
  %     disp 'One flood cycle complete'
 %       return
 %   end
    
    
    % Define N at top end
    
    NL = beta*(1-hL(2,1))-rho_i*g*basin_geometry.h_ice/N0;
    
    %%% What is this????
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Find NR and QR which fit this channel shape and BC's %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % define the end boundary condition on NR which will not change

    % Use a made up function to define a smooth NR profile which matchs these
    % two end conditions
    %         NR_temp(1,:) = sin((Ls-[1:Ls])./(Ls-1).*(pi-asin(NR(1,1)))) +[(s(:)-s(1))./(1-s(1))*NR(2,Ls)]';
    % channels

    switch method
        case 'Finite'
%             Xs(2,1) = S_end - ( (psi(Ls)*(SR(2,Ls)^(8/3)))/MR(i)^2 )^0.5;
%             % This also fixes the Q profile along channel, including inlet Q.
%             QR(2,:) = MR(i)*(s-Xs(2,1));
            QR(2,:) = sqrt(psi(Ls)*(SR(2,Ls)^(8/3))) - MR(i)*(S_end-s);


            % starting from boundary condition N_L, find N profile
            NR(2,1) = NL;
            vec1 = cumsum(ds * ( QR(2,1:Ls-1).*abs(QR(2,1:Ls-1))./(SR(2,1:Ls-1).^(8/3)) - psi(1:Ls-1) ));
            NR(2,2:Ls) = NR(2,1) + vec1;
%             for k = 1:Ls-1,
%                 NR(2,k+1) = NR(2,k) + ds * ( QR(2,k)*abs(QR(2,k))/(SR(2,k).^(8/3)) - psi(k) );
%             end
        case 'Relax'
            %%%%%%%%% Use Relaxation method %%%%%%%%%
            try

                [QR(2,:),NR(2,:),exitflagfzero,exitflagRelax] = NewRaphToFindQtop3(NbottomVariable(i),SR(2,:),psi,ds,s,epsilonR,r,delta,MR(i),NL,QR(1,:),DPing,PlotRelax);


            catch
                warning 'Model Reached a point where the Relaxation Method Failed';
                exitflag = -1;
                %save ('ModelFailTempSave')
                return
            end
       
    
            

    end
    


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%% Plotting %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if rem(i,PlotFreq)==0
        
%         ddsQ_ext = (QR(2,2:end) - QR(2,1:end-1))/ds;   
%   if nnz(ddsQ_ext<0)~=0
%       pause(0.001)
%   end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Roughly calc divide positions%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %         XR = interp1(QR(2,:),s0*s,0);
        sprintf('%3.3f %% Complete  ',t(i)/T * 100)
        %          pause
        TotalSeconds = t(i)*t0;
        TotalDays = TotalSeconds/(24*3600);
        WholeDays = fix(TotalDays)
        %         disp(WholeDays)
        disp(datestr(datenum(0,0,0,0,0,TotalSeconds),'HH:MM:SS'))
        %         disp(t(i))
        assignin('base','TotalDays', TotalDays);
        %output.TotalDays = TotalDays;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%           Plotting           %%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % QR       
        if nnz(plots==1)~=0
            if ishandle(1) ==0; figure(1);set(1,'WindowButtonDownFcn',@PauseSim);end
            set(0,'CurrentFigure',1)
            plot(s*s0,QR(2,:)*QR0,'b')

            title 'Q profiles'
            pause(0.0001)
        end
        
        % effective pressure profile
        if nnz(plots==2)~=0
            if ishandle(2) ==0; figure(2);set(2,'WindowButtonDownFcn',@PauseSim);end
            set(0,'CurrentFigure',2)
            plot(s*s0,NR(2,:)*N0,'b')
            %       axis([0 1 -3 1])
            title 'Effective Pressure Profiles'
            pause(0.0001)
        end
        
        % QR time series
           if nnz(plots==3)~=0
            if ishandle(3) ==0; figure(3);set(3,'WindowButtonDownFcn',@PauseSim);end
            set(0,'CurrentFigure',3)
            %colors = parula
            %set(groot,'defaultaxescolororder',colors)
            hold on            
            plot(tDays(i),QR(2,1)*QR0,'.r', 'MarkerSize',5)
            title 'QR time series'
            pause(0.001)
           end
           
        % QR(end) and SR(end)   time series
        if nnz(plots==4)~=0
            if ishandle(4) ==0; figure(4);set(4,'WindowButtonDownFcn',@PauseSim);end
            set(0,'CurrentFigure',4)
            hold off
            plot(s*s0,SR(2,:)*SR0,'b')
            xlabel s
            ylabel SR
            title 'Cross-Sectional Area along channel'
            pause(0.001)
        end
        
        % Effective Pressure at Terminus
        if nnz(plots==5)~=0
            if ishandle(5) ==0; figure(5);set(5,'WindowButtonDownFcn',@PauseSim);end
            set(0,'CurrentFigure',5)
            hold on
            if i>1
            plot(tDays(i),NR(2,end)*N0,'r.')
            end
            title 'N at the terminus'
            pause(0.00001)
        end
        
        %%%  Time series of the minumum area of the channel
        if nnz(plots==6)~=0
            if ishandle(6) ==0; figure(6);set(6,'WindowButtonDownFcn',@PauseSim);end
            set(0,'CurrentFigure',6)
            hold on
            if i>1
            plot(tDays(i),min(SR(2,:))*SR0,'r.')
            end
            title 'min(SR)'
            pause(0.00001)
        end
        
        % QR, NR and QC and NC profiles
        if nnz(plots==7)~=0
            if ishandle(7) ==0; figure(7);set(7,'WindowButtonDownFcn',@PauseSim);end
            set(0,'CurrentFigure',7)
            hold off
            plot(s*s0,QR(2,:)*QR0,'b',s*s0,NR(2,:)*N0,'--b')
            hold on
%             Line2 = line([XR XR],[-1 1]);
%             set(Line2,'Color','b')
            axis([0 s0 -0.1 1])
            title 'Profiles of Q and N'
            pause(0.0001)
        end
        % QR-hL phase space
        if nnz(plots==8)~=0
            if ishandle(8) ==0; figure(8);set(8,'WindowButtonDownFcn',@PauseSim);end
            set(0,'CurrentFigure',8)
            hold on
            plot(hL(2,1)*hL0,QR(2,1)*QR0,'.r','MarkerSize',5)
            xlabel hL
            ylabel QR
            grid on
            pause(0.00001)
        end
        % ubdim
        if nnz(plots==9)~=0
            if ishandle(9) ==0; figure(9);set(9,'WindowButtonDownFcn',@PauseSim);end
            set(0,'CurrentFigure',9)
            hold on
            if i>1
            plot(t(i),((QR(2,1)-QR(1,1))*QR0)/(tDays(i)-tDays(i-1)),'r.')
            end
            title 'dQR/dt'
            pause(0.00001)
        end

        if nnz(plots==10)~=0
            if ishandle(10) ==0; figure(10);set(10,'WindowButtonDownFcn',@PauseSim);end
            set(0,'CurrentFigure',10)          
            plot(s*s0,NR(2,:)*N0,'r',s*s0,SR(2,:)*SR0,s*s0,QR(2,:)*QR0)          
            title 'NR, SR and QR profiles'
            hold off
            pause(0.00001)
        end
        
        %Rate of change of discharge with lake height, dQR/dh
        if nnz(plots==11)~=0
            if ishandle(11) ==0; figure(11);set(11,'WindowButtonDownFcn',@PauseSim);end
            set(0,'CurrentFigure',11)
            hold on
            plot(t(i),((QR(2,1)-QR(1,1))*QR0)/((hL(2,1)-hL(1,1))*hL0),'r.')
            title 'dQR/dh'
            pause(0.00001)
        end

        if nnz(plots==12)~=0
            if ishandle(12) ==0; figure(12);set(12,'WindowButtonDownFcn',@PauseSim);end
            set(0,'CurrentFigure',12)
            hold off
            %plot(s*s0,sqrt(SR(2,:)*SR0)/pi)
            plot(s*s0,sqrt(2*SR(2,:)*SR0/pi))
            title 'Channel radius profile'
            xlim([0 100])
            pause(0.00001)
        end
        if nnz(plots==13)~=0
            if ishandle(13) ==0; figure(13);set(13,'WindowButtonDownFcn',@PauseSim);end
            set(0,'CurrentFigure',13)
            hold on
            plot(t(i),Qin(i),'*')
            title 'Lake Input'
            pause(0.00001)
        end

    end

    %%%%%%%%%%%%%%%%%%
    %%%  SAMPLING  %%%
    %%%%%%%%%%%%%%%%%%

    if TSamp(SampNumber)==t(i)
        hrLakeSamp(SampNumber,1) = hr(2,1);
        %         hwLakeSamp(SampNumber,1) = hw(2,1);
        %         xTposSamp(SampNumber,1) = xT(2,1);
        QRLakeSamp(SampNumber,1) = QR(2,1);
        QRendSamp(SampNumber,1) = QR(2,end);
        hLSamp(SampNumber,1) = hL(2,1);
        VLSamp(SampNumber,1) = VL(2,1);
        meanNRSamp(SampNumber,1) = mean(NR(2,:));
        QRSamp(SampNumber,:) = QR(2,:);
        SRSamp(SampNumber,:) = SR(2,:);

        SampNumber = SampNumber + 1;

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Detect Peak of floods %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if QR(2,1)<QR(1,1) && QR(1,1)>Q_old
    Peak(P) = QR(1,1);
    PeakTimeSecs(P) = t(i-1)*t0;          % the seconds since the beginning
    PeakTimeDays(P) = PeakTimeSecs(P)/3600/24;
    PeakYearFrac(P) = PeakTimeDays(P)/365-floor(PeakTimeDays(P)/365);
    P=P+1;

    %save (['FullModelt_years=' num2str(t_years(i)) '___3.mat'], 'Peak(P)')
    %Save2ws('QPeak', P); 
    assignin('base','Peak',Peak);
    assignin('base','PeakTimeDays',PeakTimeDays);
    
end
Q_old = QR(1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Detect Lake Highstand %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if hL(2,1)<hL(1,1) && hL(1,1)>hL_old
    Highstand(Lhi) = hL(1,1);
    HighstandTimeSecs(Lhi) = t(i-1)*t0;          % the seconds since the beginning
    HighstandTimeDays(Lhi) = HighstandTimeSecs(Lhi)/3600/24;
    HighstandYearFrac(Lhi) = HighstandTimeDays(Lhi)/365-floor(HighstandTimeDays(Lhi)/365);
    Lhi=Lhi+1;
    
    assignin('base','Highstand', Highstand);
    assignin('base','HighstandTimeDays', HighstandTimeDays);
    

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Detect Lake Lowstand %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if hL(2,1) >hL(1,1) && hL(1,1)<hL_old
    Lowstand(Llo) = hL(1,1);
    LowstandTimeSecs(Llo) = t(i-1)*t0;          % the seconds since the beginning
    LowstandTimeDays(Llo) = LowstandTimeSecs(Llo)/3600/24;
    LowstandYearFrac(Llo) = LowstandTimeDays(Llo)/365-floor(LowstandTimeDays(Llo)/365);
    Llo=Llo+1;
    assignin('base','Lowstand', Lowstand);
    

end
hL_old = hL(1,1);

%%% create output variable %%%

output.QR0 = QR0;
output.hL0 = hL0;
output.VL0 = VL0;
output.SR0 = SR0;
output.t = tDays;
output.QR = QRLakeSamp;
output.QRend = QRendSamp;
output.hL = hLSamp;
output.SR = SRSamp;
output.hrLake = hrLakeSamp;
output.hr0 = hr0;
output.Ls = Ls;
%output.TotalDays = TotalDays;

% assign to output variable
%output.Highstand = Highstand;
%output.HighstandTimeDays = HighstandTimeDays;
%output.Peak = Peak;
%output.PeakTimeSecs = PeakTimeSecs;
%output.PeakTimeDays = PeakTimeDays;
%output.Lowstand = Lowstand;
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Detect Convergence to a limit cycle %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Converged = P>=4 && Lhi>=4 && Llo>=4 ...
    && abs((1-Peak(P-1)/Peak(P-2))) <0.0005 && abs((1-Peak(P-1)/Peak(P-3))) <0.0005...
    && abs((1-Highstand(Lhi-1)/Highstand(Lhi-2))) <0.0005 && abs((1-Highstand(Lhi-1)/Highstand(Lhi-3))) <0.0005...
    && abs((1-Lowstand(Llo-1)/Lowstand(Llo-2))) <0.0005 && abs((1-Lowstand(Llo-1)/Lowstand(Llo-3))) <0.0005;
if Converged
    disp(['Limit Cycle Reached - change line code around line 640 in ' mfilename  '.m to change this behaivour'])

    return

end

%%%%%%%%%%%%%%%%%%
%%%%% Saving data for plotting figure 9 for the enviro chapter
%%%%%%%%%%%%%%%%%%

% if WholeDays>941
%     save DataForFigure9EnviroChap
%     
%     return
% end




    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% Loop values round %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % copy values into position 1

    hL(1,1) = hL(2,1);
    VL(1,1) = VL(2,1);
    hr(1,:) = hr(2,:);
    SR(1,:)  = SR(2,:);
    QR(1,:) = QR(2,:);
    NR(1,:) = NR(2,:);
    Xs(1,1) = Xs(2,1);


    % wipe old values
    hL(2,1) = 0;
    VL(2,1) = 0;
    hr(2,:) = 0;
    SR(2,:)  = 0;
    QR(2,:) = 0;
    NR(2,:) = 0;
    Xs(2,1) = 0;

    %%%%%%%%%%%%%%%%%%%%
    %%%  PAUSE CODE  %%%
    %%%%%%%%%%%%%%%%%%%%

    if paused ==true
        for l=1:inf
            if paused ==false
                break
            end
            pause(0.01)
        end
    end

    if UserReturn == 1
        button = questdlg('Return?');
        if strcmp(button,'Yes')
            return
        else
            UserReturn = 0;
        end
    end




end


%%
disp('Model reached t=T')
exitflag = 1




end




%%% Pause function %%%

function [paused] = PauseSim(src,evt);

global paused
if paused == false;
    paused =true;
    disp('Paused')
    return
end
if paused == true;
    paused =false;
    disp('UnPaused')
    return
end
end

function [UserReturn] = EndRun(src,evnt)
   global UserReturn
if strcmp(get(src,'SelectionType'),'alt')
      disp Return?
      UserReturn = 1;
   else
      disp('Use control-click return')
   end
end
