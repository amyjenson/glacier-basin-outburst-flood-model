SMB.max = 2; SMB.grad = .005; %continental
%SMB.max = 4; SMB.grad = .01; %maritime
precip.grad = .001;
slope = 6; %in degrees
ELA_0 = 1500;
x_basin = 0;

%'n' means no plot and 'y' means plot

%% Spinup g=lacier geometry
[precip,Qb,runoff,glacier,domain,sol] = glacier_flow_model('spinup',ELA_0,0,slope,SMB,precip,x_basin,'y');

%% Force glacier with climate warming scenario
RCP = 8.5; % can be 8.5 or 2.6; 8.5 is more drastic, 2.6 allows the glacier to restabilize
s0 = 0.75; % basin is located %% percent of the distance from the divide

[precip,Qb,runoff,glacier,domain,sol] = glacier_flow_model('gradual',ELA_0,RCP,slope,SMB,precip,s0,'y');