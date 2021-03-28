function [precip,Qb,runoff,glacier,domain,sol]=glacier_flow_model(model,ELA_0,RCP,theta,SMB,precip,s0,plots)
% Original code from Ellyn M. Enderlin, ellyn.enderlin@gmail.com
% Modified 2014 Jason M. Amundson, jason.amundson@uas.alaska.edu
% Modified for runoff modeling of land terminating glaciers, 2017
% elcarnahan@alaska.edu
%
%Formerly v9, has varying initialization structure and constant basal stress 1e5
%has ability for gradual increases in ELA as well as stepchange 
% and has varying nyears based on slope and model run 
%
%SMB array includ .grad .max; precip array includes .grad 
%
%model: specify whether the model run is a 'spinup' or for terminus
% 'advance' or 'retreat'. To for an advance with a moraine, use 'moraine'.
%
% 
warning off; % turn off warnings (velocity coefficient matrix is close to singular)


%% GLOBAL CONSTANTS
    global rho_i rho_w g m n A dx slope phi mu tau_b

    % Ice and water densities and gravitational acceleration
    rho_i = 917; % ice density (kg m^-3)
    rho_w = 1000; % water density (kg m^-3)
    g = 9.81; % acceleration (m s^-2)

    % Stress parameters
    m = 1; % basal sliding exponent
    n = 3; % flow law exponent *NOTE: cannot use nthroot for Rxx calculation if n is even
    A = 2.4e-24; % flow rate factor (Pa^-n s^-1); use constant flow rate factor for temperate ice
    
    mu = theta;
    slope = -theta*pi/180;
    phi = 0*pi/180;
    
    tau_b = 2e5; % max basal shear stress [Pa]
    
%% CONSTANTS THAT ARE NOT GLOBAL IN SCOPE
    % Desired grid spacing (m)
    dx0 = 100; 
    dx = dx0; % dx initially set to dx0; it will be adjusted later
    
    % Time step (s)
%     dt = dx0/(5000/(86400*365)); % assume max speed = 2 km/yr; make sure this satisfies the CFL condition (dt <=dx/U)
    dt = 0.05*86400*365;
     if strcmp(model,'spinup') || strcmp(model,'stepchange')
        nyears = 500; % number of years to run the simulation for
     elseif strcmp(model,'gradual')
         nyears = 500;
     end
     
    
%% LOAD INITIAL DATA, INITIALIZE THE SOLUTION VECTORS, AND SET COUNTERS
    % load the starting file
        
    if strcmp(model,'spinup') % spin-up the model
        [domain,glacier]=initialization_data; % for initial spin-up
        %load spinup.mat
        %dx = dx0; 
        
    elseif strcmp(model,'gradual') && SMB.max >=3 || strcmp(model,'stepchange') && SMB.max >=3
        load(['spinup.ELA' num2str(ELA_0) '.slope' num2str(theta) '.maritime.v8.mat'],'glacier','domain'); % for model runs
    
    elseif strcmp(model,'gradual') && SMB.max < 3 || strcmp(model,'stepchange') && SMB.max < 3
        load(['spinup.ELA' num2str(ELA_0) '.slope' num2str(theta) '.continental.v8.mat'],'glacier','domain'); % for model runs
    
    end

    L=glacier.x(end); % initial glacier length;
    
      if strcmp(model,'spinup') %set basin at divide for spinup
         
       x_basin = 0;
        
      elseif strcmp(model,'gradual') %then get basin distance from the divide after setting length from the terminus
      
          %set the basin some specific distance from the terminus 
          if s0 > 1
            x_basin = glacier.x(end)-(s0*cos(slope));  
       
          %set basin some percentage of the glacier length from the divide 
          elseif s0 < 1
            x_basin = glacier.x(end)*s0/cos(slope);
        
           
          elseif s0 == 1 %set basin where the glacier thickness is consistent
            x_basin = interp1(glacier.H, glacier.x, 300); 
      
          end
       end 

    % Initialize solution structure for time series plots
    sol.L = [];
    sol.SMB = [];
    sol.V = [];
    sol.t = [];    
    sol.dx = [];
    sol.AAR = [];
    sol.H_basin = [];
    sol.SMB_basin = [];
    sol.dhdx_basin = [];
    
    % Initialize Catchment Basin Precipitation and Runoff structures
    precip.F = []; %amount of rainfall on glacier in m we. per year
    
    Qb.F = []; %amount of runoff from glacier in m w.e. per year
    
    runoff.F = []; % sum of precip on glacier and negative mass balance in m we. per year
    
    % Set the timers
    year = 0; %keeps track of the number of model years 
    year_previous = 0; %used to plot profiles once per model year
    year_start = year; %tracks the time elapsed
    year_end = year_start+nyears; %stops time for the model


%% RUN FLOWLINE MODEL
j=1;
while j
    %% SOLVE THE STRESS BALANCE EQUATIONS TO FIND GLACIER VELOCITIES
    
        [glacier] = U_convergence(glacier,dx,L,model,year);
        
    %% CALCULATE THE ICE FLUX AT EACH GRID POINT
        F = glacier.U.*glacier.H.*glacier.W; % ice flux (m^2 s^-1)
        
    %% CALCULATE SURFACE MASS BALANCE AND USE MASS CONTINUITY TO FIND NEW GLACIER THICKNESS DISTRIBUTION
    
        % Assume constant balance gradient with elevation, and place maximum
        % bound on accumulation


       if strcmp(model,'gradual')
            if RCP == 2.6
                if year < 100
                ELA_year = ELA_0+158*(1-exp(-1/28*year));
                end
                
            elseif RCP == 8.5         
                ELA_year = ELA_0+1*(year);
            end
            
       elseif strcmp(model,'spinup') % spin-up the model
            ELA_year = ELA_0; % for initial spin-up
       end
                
        glacier.SMB = SMB.grad*(glacier.h-ELA_year);
        glacier.SMB(glacier.SMB>SMB.max) = SMB.max;
        glacier.SMB = glacier.SMB/(86400*365); % put this in m/s
        
        % Calculate the change in ice thickness from mass continuity
        dHdt = glacier.SMB-gradient(F,glacier.x)./glacier.W;

        % New glacier thickness.
        glacier.H = glacier.H + dHdt*dt;
           
        

    %% FIND NEW TERMINUS LOCATION

        dLdt = glacier.U(end);
        L = L+dLdt*dt; % find new glacier length
            
    %% ADJUST GRID TO ACCOUNT FOR NEW GLACIER LENGTH
       
        [glacier,dx] = track_glacier_terminus(domain,glacier,L);
        
        
        if glacier.H(end) < 1 % if terminus is less than 1 m thick, find the location where it is 0.1 m thick...
            L = interp1(glacier.H,glacier.x,1,'linear');
            [glacier,dx] = track_glacier_terminus(domain,glacier,L);
        end
        
  

    %% ADVANCE THE MODEL TIME    
        year = year + dt/31536000;
        
         %% Calculate Catchment Basin Precipitation and Runoff
  
    % recalculate the SMB in case adjusting the grid resulted in a
        % different number of grid points
        glacier.SMB = SMB.grad*(glacier.h-ELA_year);
        glacier.SMB(glacier.SMB>SMB.max) = SMB.max;
        glacier.SMB = glacier.SMB/(86400*365); % put this in m/s
        
            %precip.grad = (SMB.max*rho_i/rho_w-precip.sealevel)/(SMB.max/SMB.grad+1000); %allows gradient to stay constant
            %and keeps value of precipitation above SMB rate by allowing the
            %gradient to vary and the sea level precipitation rate to be input
    precip.sealevel = SMB.max*(1-precip.grad/SMB.grad)-precip.grad*1000; % allows sea level precip to vary depending 
            %gradient prescribed, 1000 represents the minimum ELA to allow for
            %precip to remain greater then SMB
            
    precip.sealevel = precip.sealevel*(rho_i/rho_w); %keeps gradient constant and prescribes sea level precip rate again 
   
   
    precip.elev = precip.sealevel + precip.grad * glacier.h; %precipitation in m per year only on glacier
       %precip.elev(precip.elev>4.580) = 4.580; %max precip converted to m we. to find max precip rate
   
    precip.flux = trapz(precip.elev.*glacier.W)*dx; % net precipitation over glacier surface
        %in m^3 we./yr
    QbF = - trapz(glacier.SMB.*glacier.W)*.917*86400*365*dx; %SMB converted from m ie. per second
        %to m we. /yr integrated over the glacier to give m^3 we./yr 
    runoffF = precip.flux + QbF;
   
        
%% Plot the annual profiles & display select variables in the command window
                
        sol.L = [sol.L L];
        [~,aar_i] = min(abs(glacier.h - ELA_year));
        sol.AAR = [sol.AAR glacier.x(aar_i)/L];
        sol.SMB = [sol.SMB trapz(glacier.SMB.*glacier.W)*dx*86400*365];
        sol.V = [sol.V trapz(glacier.H.*glacier.W)*dx];
        sol.t = [sol.t year];
        sol.dx = [sol.dx dx];
        precip.F = [precip.F precip.flux];
        Qb.F = [Qb.F QbF];
        runoff.F = [runoff.F runoffF];
      
        sol.H_basin = [sol.H_basin interp1(glacier.x, glacier.H, x_basin, 'cubic')];
        
        sol.SMB_basin = [sol.SMB_basin interp1(glacier.x, glacier.SMB, x_basin, 'cubic')];
        
        sol.dhdx_basin = [sol.dhdx_basin interp1(glacier.x, glacier.dhdx, x_basin, 'cubic')];
        
        balance_flux = cumtrapz(glacier.SMB)*dx*86400*365;
        balance_velocity = balance_flux./glacier.H;
        
       if year - year_previous >= 1                   
        
         disp(['Year = ' num2str(year) '']);
         disp(['Length (km) = ' num2str(L/1000) '']);
         disp([' ']);
        
        yearlabel = num2str(round(year), '%03.f');
        
        if strcmp(model,'gradual')            
            glacier.x_basin = x_basin;
            
           if SMB.max >= 3 
               mkdir(['glacier_model/output/ELA1500/RCP8.5/slope' num2str(theta) '/maritime/s0=' num2str(s0) ]) 
               save(['glacier_model/output/ELA1500/RCP8.5/slope' num2str(theta) '/maritime/s0=' num2str(s0) '/output.' yearlabel '.mat'],'domain', 'glacier','sol');
           end
           
           if SMB.max <= 3
               mkdir(['glacier_model/output/ELA1500/RCP8.5/slope' num2str(theta) '/continental/s0=' num2str(s0) ])
               save(['glacier_model/output/ELA1500/RCP8.5/slope' num2str(theta) '/continental/s0=' num2str(s0) '/output.' yearlabel '.mat'],'domain', 'glacier','sol');
           end
            
        end
       
       if strcmp(plots,'y') 
       f1 = figure(1); clf
       
        subplot(2,1,1)
            %plotyy(glacier.x/1000,glacier.U*86400*365,glacier.x/1000,rho_i*g*glacier.H.*gradient(glacier.h,glacier.x)*1e-5)
            plot(glacier.x/1000,glacier.U*86400*365,'k','linewidth',1)
            hold on
            plot(x_basin/1000,0,'ro','MarkerSize',4,'MarkerFaceColor','red')
            %plot([0 100],[0 0],'k','linewidth',0.5)
            %hold off
            ylim([-10 500]);
            xlim([0 100])
            xlabel('Longitudinal coordinate, km','fontsize',10)
            ylabel('Velocity, m/yr','fontsize',10)
            text(105,5000*5/6,['year: ' yearlabel ''],'fontsize',10)
            set(gca,'linewidth',1,'fontsize',10)
            box on
    
    
         subplot(2,1,2)
            h=patch([glacier.x/1000 glacier.x(end:-1:1)/1000],[glacier.h glacier.h(end:-1:1)-glacier.H(end:-1:1)],'w');
%             set(h,'linewidth',2);
            h=patch([0 domain.x/1000 50],[-500 domain.hb -500],[0.2 0.2 0]);    
            hold on
            plot(x_basin/1000,interp1(domain.x,domain.hb,x_basin),'ro','MarkerSize',4,'MarkerFaceColor','red')
            
            xlabel('Longitudinal coordinate (km)','fontsize',10)
            ylabel('Elevation (m)','fontsize',10)
             ylim([0 3000]);xlim([0 100]);
            text(70,2300,['year: ' yearlabel ''],'fontsize',10)
            text(30,1900,['rate of advance or retreat: ' num2str(round((sol.L(end)-sol.L(end-1))/(dt/(365*86400)))) ' m/yr'],'fontsize',10)
            box on
            set(gca,'linewidth',1,'fontsize',10)
            
            drawnow
           % saveas(gcf, ['./images/year' yearlabel '.png'],'png')
       end
      end
     
  
    
    %stop the model at a set time & save
    % if (year > year_end && all(abs(sol.L(end-62:end)-sol.L(end-63:end-1))./(dt/(365*86400)) < 2)) | sol.L(end)< 10; %abs(glacier.U(end)-uc)*86400*365 < 1 && year > 50 || glacier.x(end)/1000 > 127.5; %if dL/dt < 1 m/yr, then stop
    % if sol.L(end) < x_basin | glacier.H(1)<10 | (year > year_end && all(abs(sol.L(end-1:end)-sol.L(end-2:end-1))./(dt/(365*86400)) < 1));
  
     if (year > year_end && all(abs(sol.L(end-62:end)-sol.L(end-63:end-1))./(dt/(365*86400)) < 0.5)) | sol.L(end) < x_basin ; % | (year > year_end && all(abs(sol.L(end-1:end)-sol.L(end-2:end-1))./(dt/(365*86400)) < 1));
  
        if strcmp(model,'spinup') && SMB.max >= 3
            save(['spinup.ELA' num2str(ELA_0) '.slope' num2str(theta) '.maritime.v8.mat'],'domain', 'glacier','sol');
        
        elseif strcmp(model,'spinup') && SMB.max < 3
            save(['spinup.ELA' num2str(ELA_0) '.slope' num2str(theta) '.continental.v8.mat'],'domain', 'glacier','sol');
        
        elseif strcmp(model,'advance')
            save(['output.ELA' num2str(ELA) '.alpha' num2str(alpha) '.advance.mat'],'domain','glacier','sol');
                    
        elseif strcmp(model,'gradual') &&  SMB.max >= 3
            save(['output.ELA' num2str(ELA_0) '.RCP' num2str(RCP) '.slope' num2str(theta) '.maritime.v9.mat'],'domain','glacier','sol','precip','Qb','runoff');
        
        elseif strcmp(model,'gradual') &&  SMB.max < 3
            save(['output.ELA' num2str(ELA_0) '.RCP' num2str(RCP) '.slope' num2str(theta) '.continental.v9.mat'],'domain','glacier','sol','precip','Qb','runoff');  
        
        elseif strcmp(model,'stepchange') &&  SMB.max >= 3
            save(['Maritime Outputs 2000 v8/' 'output.ELA' num2str(ELA_0) 'RCP' num2str(RCP) '.slope' num2str(theta) '.stepchange.maritime.v9.mat'],'domain','glacier','sol','precip','Qb','runoff');
        
        elseif strcmp(model,'stepchange') &&  SMB.max < 3
            save(['Continental Outputs 2000 v8/' 'output.ELA' num2str(ELA_0) '-' num2str(ELA) '.slope' num2str(theta) '.stepchange.continental.v9.mat'],'domain','glacier','sol','precip','Qb','runoff');
        
        end
        
        yrs = year_start:1:year_end;
        
        break
    end
   % end
    %refresh the plot counter
    year_previous = floor(year);
    
    %loop through
    j=j+1;

end

end



%% solve the stress balance equations to obtain speed values (U)
function [glacier] = U_convergence(glacier,dx,L,model,year)
global rho_i rho_w g m n A tau_b
    
    % Calculate the thickness on the staggered grid for use in stress balance equations
    % The staggered grid starts between H_1 and H_2 and goes to half a grid
    % point past the terminus. The staggered grid does not need to straddle
    % H_1 because we are using a Dirichlet boundary condition there.
    
    glacier.Hm = (glacier.H(1:end-1)+glacier.H(2:end))/2; % use forward differences
    glacier.Hm = [glacier.Hm,2*glacier.H(end)-glacier.Hm(end)]; % Assume constant slope through the last couple of data points
    %glacier.Hm(glacier.Hm<0) = 0; % cannot have negative values
    
    x = glacier.x;
    U = glacier.U;
    %U(U<0) = 0;
    dUdx = glacier.dUdx;
    Hm = glacier.Hm;
    H = glacier.H;
    h = glacier.h;
    dhdx = glacier.dhdx;
    hb = glacier.hb;
    W = glacier.W;
    
b=1;
while b
    
    
%% CALCULATE THE LINEARIZATION TERMS AND EFFECTIVE VISCOSITY FOR THE INVERSION OF THE STRESS COEFFICIENT MATRIS
    if n == 3;
        gamma = abs(U).^((1-n)/n); % RUNS INTO PROBLEMS WHEN U<0!!!
        gamma(gamma>1e+06) = 1e+06; %set the limit so gamma does not approach infinity (minimum U = 1e-09 m s^-1)
        
        Um = (U(1:end-1) + U(2:end))./2; %speed on the staggered grid from forward differencing
        Um = [Um,2*U(end)-Um(end)]; %linearly interpolate velocity through the terminus to the last staggered grid point
        dUmdx = gradient(Um,x);
        
        vm = ((A).^(-1/n)).*(abs(dUmdx)).^((1-n)/n); %effective viscosity (Pa s)
        vm(1) = vm(2);
        vm(vm>8e+16) = 8e+16; %set a maximum value for very low strain rates
        
        eta = ones(1,length(x)); %if m=1, the basal resistance term does not need to be linearized
    end
    
%% SET UP COEFFICIENT VECTORS FOR THE LINEARIZED STRESS TERMS

    % [C(k)*U(k-1)+E(k)*U(k)+G(k)*U(k+1)=Td]
    G_minus = (2./(dx.^2)).*Hm(1:end-2).*vm(1:end-2); %for U(k-1)
               
    G = (-2/dx^2)*(Hm(1:end-2).*vm(1:end-2)+Hm(2:end-1).*vm(2:end-1))-...
        2*gamma(2:end-1).*H(2:end-1)./W(2:end-1).*((5./(A*W(2:end-1))).^(1/n)) - tau_b/max(U); %for U(k)
    
    %G = (-2/dx^2)*(Hm(1:end-2).*vm(1:end-2)+Hm(2:end-1).*vm(2:end-1))-...
    %    2*gamma(2:end-1).*H(2:end-1)./W(2:end-1).*((5./(A*W(2:end-1))).^(1/n)) - rho_i*g*H(2:end-1).*abs(U(2:end-1)).^(-2/3)*tau_b; %for U(k)
    
    G = [1 G]; % divide coefficient
    
    G_plus = (2./(dx.^2)).*Hm(2:end-1).*vm(2:end-1); %for U(k+1)
    G_plus = [0 G_plus]; % divide coefficient
    
    T = (rho_i*g*H(2:end-1).*dhdx(2:end-1));%+tau_b; % gravitational driving stress

    T = [0 T]; % driving stress at divide

    
%% TERMINUS BOUNDARY CONDITION

% 1st-order approx. of terminus boundary condition --- linearized stress
% terms
%    G_minus = [G_minus -1];
%    G = [G 1];


% 2nd-order approx. of terminus boundary condition --- linearized stress
% terms
    G_minus = [G_minus (2/dx^2)*(Hm(end-1).*vm(end-1)+Hm(end).*vm(end))];%not sure why element multiplication
 
    G = [G (-2/dx^2)*(Hm(end-1)*vm(end-1)+Hm(end)*vm(end))-...
        2*gamma(end)*H(end)/W(end)*(5/(A*W(end)))^(1/n)-tau_b/max(U)];
    
    %G = [G (-2/dx^2)*(Hm(end-1)*vm(end-1)+Hm(end)*vm(end))-...
    %    2*gamma(end)*H(end)/W(end)*(5/(A*W(end)))^(1/n)-rho_i*g*H(end).*abs(U(end)).^(-2/3)*tau_b];
    
    eps_xx = A*(rho_i*g*glacier.H(end)/4)^3; 
    
%     eps_xx = 0; % just testing if this does anything... 
% 1st-order approx. of longitudinal strain rate at the terminus    
%     T = [T eps_xx*dx];

% 2nd-order approx. of longitudinal strain rate at the terminus
    T = [T rho_i*g*H(end)*dhdx(end) - 4/dx*Hm(end)*vm(end)*eps_xx]; % should 2 be replaced with 4???

%%
    %create a sparse tri-diagonal matrix for the velocity coefficient vectors
    M = diag(G_minus,-1) + diag(G) + diag(G_plus,1);
    M = sparse(M);
    
    %make sure Td is a column vector for the inversion 
    T=T(:);
    
    %use the backslash operator to perform the matrix inversion to solve for ice velocities
    Un = M\T; %velocity (m s^-1)
    
    %remove NaNs and apply the ice divide bounday condition
    Un(isnan(Un)) = 0;
    %Un(Un<0) = 0;
    
    %calculate new strain rates (equal to the velocity gradient)
    dUndx = gradient(Un,x); %strain rate (s^-1)
    
    %make sure Un is a row vector so it can be compared with U
    Un = Un(:)';
    
    if isreal(Un)==0
        error('U has become imaginary')
    end
    
    %make sure dUndx is a row vector
    dUndx = dUndx(:)';
    
    %check if the difference in speed between iteratons (U vs. Un) meets a set tolerance
    if abs(sum(U)-sum(Un))<0.01*abs(sum(U)); %determine if U has converged sufficiently
        %use sufficiently converged values for speeds & strain rates
            glacier.U = Un; 
            glacier.dUdx = dUndx;
            
           
            return %break the U iterations
    end
    
    %if not sufficiently converged, set Un to U and solve the stress balance matrix again
    U = Un;
    dUdx = dUndx;
    
    %U(U<0)=0;
    
    
    %terminate the U iteration loop if convergence takes too long
    if str2double(num2str(b)) > str2double(num2str(100));
        b = b;
        return
    end
    
    %loop through
    b = b+1;
    
end

end


%%
function [glacier,dx] = track_glacier_terminus(domain,glacier,L)
% Adjust grid back so that the final grid point lies on the terminus,
% regardless of whether or not the terminus is floating.

%adjust the grid spacing to track the terminus
%xl = floor(L/dx0); % number of ideal grid spaces to reach the terminus

xl = 201;
dx = L/xl; % new grid spacing (should be ~dx0)

xn = 0:dx:L; % new distance vector

% adjust geometry to the new distance vector
    glacier.hb = interp1(domain.x,domain.hb,xn,'linear','extrap');
    glacier.W = interp1(domain.x,domain.W,xn,'linear','extrap');
    
    %adjust the space-dependent variables to the new distance vector
    glacier.H = interp1(glacier.x,glacier.H,xn,'linear','extrap'); %ice thickness (m)
    glacier.U = interp1(glacier.x,glacier.U,xn,'linear','extrap'); %speed (m s^-1)
    glacier.dUdx = interp1(glacier.x,glacier.dUdx,xn,'linear','extrap'); %strain rate (s^-1)
    
    % rename the distance vector
    glacier.x = xn; % distance from the divide (m)
    
    % calculate the new surface elevation and slope
    glacier.h = glacier.hb+glacier.H; % grounded ice surface elevation (m a.s.l.)
    glacier.dhdx = gradient(glacier.h,glacier.x);
end

%% Initialize Model
function [domain,glacier] = initialization_data
% for two overdeepenings

global slope phi mu dx %%%%angle
dx = 100;
max_elev = 2000;

domain.x = 0:dx:120000;

domain.hb = max_elev+slope*domain.x;

% Glacier width
domain.W = 4000*ones(size(domain.x))+tan(phi)*domain.x;


% Initial surface elevation and slope variation to prevent numerical
% instability in first few time steps
if mu == 2
    glacier.x = 0:dx:27923;
    glacier.h = -6.2e-7*(glacier.x).^2-0.026*glacier.x+2.3e+3;
elseif mu == 3
    i_thickness = 350;
    glacier.x = 0:dx:5000;
    angle = -atand((tand(mu)*glacier.x(end)+i_thickness)/glacier.x(end))*pi/180; %set angle of glacier so intital gemoemtry is linear
    glacier.h = (max_elev+i_thickness)+angle*glacier.x;
elseif any(mu == [4,5,6])
    i_thickness = 200;
    glacier.x = 0:dx:2500;
    angle = -atand((tand(mu)*glacier.x(end)+i_thickness)/glacier.x(end))*pi/180; %set angle of glacier so intital gemoemtry is linear
    glacier.h = (max_elev+i_thickness)+angle*glacier.x;
elseif mu ==7
    i_thickness = 150;
    glacier.x = 0:dx:1850;
    angle = -atand((tand(mu)*glacier.x(end)+i_thickness)/glacier.x(end))*pi/180; %set angle of glacier so intital gemoemtry is linear
    glacier.h = (max_elev+i_thickness)+angle*glacier.x;
elseif any(mu == [8,9])
    i_thickness = 100;
    glacier.x = 0:dx:1250;
    angle = -atand((tand(mu)*glacier.x(end)+i_thickness)/glacier.x(end))*pi/180; %set angle of glacier so intital gemoemtry is linear
    glacier.h = (max_elev+i_thickness)+angle*glacier.x;
else 
    i_thickness = 75;
    glacier.x = 0:dx:1000;
    angle = -atand((tand(mu)*glacier.x(end)+i_thickness)/glacier.x(end))*pi/180; %set angle of glacier so intital gemoemtry is linear
    glacier.h = (max_elev+i_thickness)+angle*glacier.x;
end

%find slope of glacier to make initialization shape linear
%  angle = -atand((tand(mu)*glacier.x(end)+i_thickness)/glacier.x(end))*pi/180 %set angle of glacier so intital gemoemtry is linear
%  glacier.h = (max_elev+i_thickness)+angle*glacier.x;


glacier.dhdx = gradient(glacier.h,glacier.x);

% Initial velocity; 0.1 m/day at the terminus
glacier.U = (0.1/86400)*glacier.x/glacier.x(end);
glacier.dUdx = gradient(glacier.U,glacier.x);

% Initial bed geometry
glacier.hb = domain.hb(1:length(glacier.x));
glacier.W = domain.W(1:length(glacier.x));

% Initial glacier thickness
glacier.H = glacier.h-glacier.hb;

%plot(domain.x,domain.hb,glacier.x,glacier.h)
end