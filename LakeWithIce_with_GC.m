function basin_geometry = LakeWithIce_with_GC(sol,a,pL,basin)

global rho_i rho_w

H_basin = sol.H_basin(end);         % ice thickness of glacier at the basin [m] 
year0 = 1;                          % year that the basin first formed

if strcmp(basin,'noice')
    %V_ice0 = 0; %no ice in basin
    %V_ice = 0; %no ice in basin
    basin_geometry.h_ice = 0; 
    basin_geometry.H_basin = H_basin;
    basin_geometry.hLi = rho_i/rho_w*(H_basin);
    basin_geometry.VLi = a/pL*basin_geometry.hLi^pL;  % calculate the volume of water stored in the basin

elseif strcmp(basin,'ice')
        dt = int16(1/sol.t(1));
       
        V_ice0 = a/pL*sol.H_basin(year0*dt).^pL; %the amount of ice once basin forms
        V_ice = V_ice0 + trapz(sol.SMB_basin(sol.t>year0)*86400*365*a.*sol.H_basin(sol.t>year0).^(pL-1))*sol.t(1);

        if V_ice < 0
          V_ice = 0;
        end

     % determine floating ice thickness so that effective pressure at basin
     % entrance equals 0
    
%    x = linspace(0,400);                                 % Interval To Evaluate Over
%    f = @(h_ice) volume_ice(a,pL,H_basin,V_ice,h_ice);    % Function
%    fx = f(x); % Function Evaluated Over ?x?
%    fz = [];
%    cs = fx.*circshift(fx,0,400);                         % Product Negative At Zero-Crossings
%    h0 = x(cs <= 0);                                      % Values Of ?x? Near Zero Crossings
%for k1 = 1:length(h0)
%    fz(k1) = fzero(f, h0(k1));                            % Use ?h0? As Initial Zero Estimate
%end
     
   % basin_geometry.h_ice = fz(k1);
     
     h0 = 0; % initial guess for the ice thickness
     basin_geometry.h_ice = fzero(@(h_ice) volume_ice(a,pL,H_basin,V_ice,h_ice), h0);

    % now that we know the floating ice thickness, calculate the water level
    [~,basin_geometry.hLi] = volume_ice(a,pL,H_basin,V_ice,basin_geometry.h_ice);

  

    %basin_geometry.hLi = 117.5102;
    %basin_geometry.h_ice = (rho_i/rho_w*H_basin-basin_geometry.hLi)*rho_w/rho_i;

    
 elseif strcmp(basin,'ice_thickness_defined')
    
     h_ice = 240; 
     basin_geometry.h_ice = h_ice;
     V0 = 0;
     V_ice = fzero(@(V_ice) volume_ice(a,pL,H_basin,V_ice,h_ice), V0);
     
    % now that we know the floating ice thickness, calculate the water level
    [~,basin_geometry.hLi] = volume_ice(a,pL,H_basin,V_ice,basin_geometry.h_ice);

end 


    if rho_i / rho_w * H_basin < abs(basin_geometry.hLi)
        disp('Error in LakeWithIce_GC')
      return
    end
    
    if rho_i / rho_w * H_basin < abs(basin_geometry.h_ice)
             disp('Error in LakeWithIce_GC')
      return
    end
    
    if basin_geometry.hLi  < 0
             disp('Error in LakeWithIce_GC')
      return
    end
     
  % calculate the volume of water stored in the basin
  basin_geometry.VLi = a/pL*basin_geometry.hLi^pL;
  basin_geometry.year0 = year0; 
  basin_geometry.H_basin = H_basin;

 end

function [val,hLi] = volume_ice(a,pL,H_basin,V_ice,h_ice)
global rho_i rho_w

% This subfunction is used for finding the thickness of ice floating in the
% basin, given its volume and the basin geometry

hLi = rho_i/rho_w*(H_basin-h_ice);
val = V_ice - (a/pL*((hLi+h_ice)^pL-hLi^pL));
end