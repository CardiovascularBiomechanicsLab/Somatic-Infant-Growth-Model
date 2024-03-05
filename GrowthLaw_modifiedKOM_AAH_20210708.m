function [r0_new, h0_new,  Fg_new, Fg_R_new, st, sl, sl_ex1, st_ex1, st_ex2, setpoints]=GrowthLaw_modifiedKOM_AAH_20210708( r, h, Fg_old, Fg_old_R, growthparams, HLHS, r0_s, h0_s, iteration_growth, somaticSetpoints, generateSetpoints) 
%GrowthLaw_modifiedKOM: determines the new grown unloaded dimensions of the ventricles

% Inputs:
%   r: loaded radius
%   h: loaded thickness
%   Fg_old: previous step's LV growth tensor
%   Fg_old_R: previous step's RV growth tensor
%   growthparams: input parameters for growth (defined at top of main function)
%   HLHS: flag for hypoplastic left heart syndrome
%   r0_s: somatic unloaded radius
%   h0_s: somatic unloaded thickness
%   iteration_growth: what "day" of growth the simulation is currently on
%   somaticSetpoints: homeostatic strain setpoints
%   generateSetpoints: flag that determines if new homeostatic strain setpoints need to be saved

% Outputs:
%   r0_new: new grown unloaded radius
%   h0_new: new grown unloaded thickness
%   Fg_new: new value of LV growth tensor
%   Fg_R_new: new value of RV growth tensor
%   st: thickening stimulus
%   sl: lengthening stimulus
%   sl_ex1: max fiber strain throughout cardiac cycle
%   st_ex1: max radial strain throughout cardiac cycle
%   st_ex2: baseline radial strain from beginning of cardiac cycle
%   setpoints: homeostatic strain setpoints

% Growth Law Constants
%Fit to PO and VO fitting simulations by CMW
f_ff_max=0.1; f_cc_max=0.1;
f_f=growthparams(1); sl_50=growthparams(2); c_fpos=growthparams(3);  c_fneg=growthparams(4); st_50pos=growthparams(5); st_50neg=growthparams(6); 

% Current Loading
    % (Feg) Calculate the stretch within the ventricle, note that the reference configuration is the original unloaded radius and thickness 
    r_stretch(:,1)=(r(:,1)/r0_s(iteration_growth,1)); % LV
    r_stretch(:,2)=(r(:,2)/r0_s(iteration_growth,2)); % RV
    h_stretch(:,1)=(h(:,1)/h0_s(iteration_growth,1)); % LV
    h_stretch(:,2)=(h(:,2)/h0_s(iteration_growth,2)); % RV
    
    
    %Calculate the elastic fiber and radial stretches and strains where Fe = Feg * 1/Fg 
    fiber_stretch(:,1)=r_stretch(:,1) * 1/Fg_old(1,1,iteration_growth); % LV
    fiber_stretch(:,2)=r_stretch(:,2) * 1/Fg_old_R(1,1,iteration_growth); % RV
    radial_stretch(:,1)=h_stretch(:,1) * 1/Fg_old(3,3,iteration_growth); % LV
    radial_stretch(:,2)=h_stretch(:,2) * 1/Fg_old_R(3,3,iteration_growth); % RV
    fiber_strain(:,1)= 0.5*(fiber_stretch(:,1).^2 - 1); % LV
    fiber_strain(:,2)= 0.5*(fiber_stretch(:,2).^2 - 1); % RV
    radial_strain(:,1)=0.5*(radial_stretch(:,1).^2 - 1); % LV
    radial_strain(:,2)=0.5*(radial_stretch(:,2).^2 - 1); % RV

    % Determining Growth Set Points    
    r_stretch_baseline(:,1)=(r(:,1)/r0_s(iteration_growth,1)); % LV
    r_stretch_baseline(:,2)=(r(:,2)/r0_s(iteration_growth,2)); % RV
    h_stretch_baseline(:,1)=(h(:,1)/h0_s(iteration_growth,1)); % LV
    h_stretch_baseline(:,2)=(h(:,2)/h0_s(iteration_growth,2)); % RV
    
    
    fiber_strain_baseline(:,1)= 0.5*(r_stretch_baseline(:,1).^2 - 1); % LV
    fiber_strain_baseline(:,2)=0.5*(r_stretch_baseline(:,2).^2 - 1); % RV
    radial_strain_baseline(:,1)=0.5*(h_stretch_baseline(:,1).^2 - 1); % LV
    radial_strain_baseline(:,2)=0.5*(h_stretch_baseline(:,2).^2 - 1); % RV
    % LV
    setpoints(1)=max(fiber_strain_baseline(:,1));
    setpoints(2)=max(radial_strain_baseline(:,1));
    % RV
    setpoints(3)=max(fiber_strain_baseline(:,2));
    setpoints(4)=max(radial_strain_baseline(:,2));
    
    
    % to be passed into plotting_stimuli function
    st_ex2 = [setpoints(2),setpoints(4)];
    sl_ex1 = [max(fiber_strain(:,1)), max(fiber_strain(:,2))];
    st_ex1 = [max(radial_strain(:,1)), max(radial_strain(:,2))];
    
      
% Calculating Growth Stimuli 
if generateSetpoints == 1 % if we need to generate a new set of somatic setpoints
    sl=[max(fiber_strain(:,1)) - setpoints(1), max(fiber_strain(:,2)) - setpoints(3)];
    st=[-(max(radial_strain(:,1))-setpoints(2)), -(max(radial_strain(:,2))-setpoints(4))];
else % use established setpoints
    sl=[max(fiber_strain(:,1)) - somaticSetpoints(1), max(fiber_strain(:,2)) - somaticSetpoints(3)];
    st=[-(max(radial_strain(:,1))-somaticSetpoints(2)), -(max(radial_strain(:,2))-somaticSetpoints(4))];
end


%% MODIFIED GROWTH LAW FOR SPHERE 
%fiber direction

% quiescent zone for sl and st if within 1e-4

if (sl(1) < 1e-4 && sl(1) > -1e-4)
    sl(1) = 0;
end

if (sl(2) < 1e-4 && sl(2) > -1e-4)
    sl(2) = 0;
end

if (st(1) < 1e-4 && st(1) > -1e-4) 
    st(1) = 0;
end

if (st(2) < 1e-4 && st(2) > -1e-4)
    st(2) = 0;
end

%% LEFT VENTRICLE
% if sl = 0, set Fgi to eye(3)
if (st(1) == 0 && sl(1) == 0) || HLHS == 1
    Fgi = eye(3);
elseif st(1) == 0 && sl(1) ~=0 % if only a lengthening stimulus occurring
    Fgi(1,1)= logical(sl(1)>0)*sqrt(  (f_ff_max) / (1+exp(-f_f*(sl(1)-sl_50)) )  +1 ) + ...
             logical(sl(1)<0)*sqrt(  -f_ff_max/ (1+exp(f_f*(sl(1)+sl_50)) )+1 ) ;  %eqn 8        
   Fgi(2,2)=Fgi(1,1);
   Fgi(3,3)= 1;
elseif st(1) ~= 0 && sl(1) == 0 % if only a thickening stimulus occurring
    Fgi(1,1) = 1;
    Fgi(2,2)=Fgi(1,1);
    Fgi(3,3)= logical(st(1)>0)*( (f_cc_max) /(1+exp(-c_fpos*(st(1)-st_50pos)) )+1 ) + ...
              logical(st(1)<0) *( -f_cc_max/(1+exp(c_fneg*(st(1)+st_50neg)) )+1 );
else
   Fgi(1,1)= logical(sl(1)>0)*sqrt(  (f_ff_max) / (1+exp(-f_f*(sl(1)-sl_50)) )  +1 ) + ...
             logical(sl(1)<0)*sqrt(  -f_ff_max/ (1+exp(f_f*(sl(1)+sl_50)) )+1 ) ;  %eqn 8        
   Fgi(2,2)=Fgi(1,1); %eqn 10 modified - cross-fiber and fiber growth stretches are the same and radial growth stretch is different

%radial direction        
    Fgi(3,3)= logical(st(1)>0)*( (f_cc_max) /(1+exp(-c_fpos*(st(1)-st_50pos)) )+1 ) + ...
              logical(st(1)<0) *( -f_cc_max/(1+exp(c_fneg*(st(1)+st_50neg)) )+1 );  %eqn 9   
end
    Fg_new = Fgi * Fg_old(:,:,iteration_growth); 

  %Calculate the new unloaded geometry  - note that the reference is the original unloaded state  
  r0_new=zeros(1,2); h0_new=zeros(1,2);
  r0_new(1,1)=r0_s(iteration_growth,1)*Fg_new(1,1);
  h0_new(1,1)=h0_s(iteration_growth,1)*Fg_new(3,3);
  
%% RIGHT VENTRICLE 
% if sl = 0, set Fgi to eye(3)
if st(2) == 0 && sl(2) == 0
    Fgi_R = eye(3);
elseif st(2) == 0 && sl(2) ~=0 % if only a lengthening stimulus occurring
    Fgi_R(1,1)= logical(sl(2)>0)*sqrt(  (f_ff_max) / (1+exp(-f_f*(sl(2)-sl_50)) )  +1 ) + ...
             logical(sl(2)<0)*sqrt(  -f_ff_max/ (1+exp(f_f*(sl(2)+sl_50)) )+1 ) ;  %eqn 8        
   Fgi_R(2,2)=Fgi_R(1,1);
   Fgi_R(3,3)= 1;
elseif st(2) ~= 0 && sl(2) == 0 % if only a thickening stimulus occurring
    Fgi_R(1,1) = 1;
    Fgi_R(2,2)=Fgi_R(1,1);
    Fgi_R(3,3)= logical(st(2)>0)*( (f_cc_max) /(1+exp(-c_fpos*(st(2)-st_50pos)) )+1 ) + ...
              logical(st(2)<0) *( -f_cc_max/(1+exp(c_fneg*(st(2)+st_50neg)) )+1 );
else
   Fgi_R(1,1)= logical(sl(2)>0)*sqrt(  (f_ff_max) / (1+exp(-f_f*(sl(2)-sl_50)) )  +1 ) + ...
             logical(sl(2)<0)*sqrt(  -f_ff_max/ (1+exp(f_f*(sl(2)+sl_50)) )+1 ) ;  %eqn 8        
   Fgi_R(2,2)=Fgi_R(1,1); %eqn 10 modified - cross-fiber and fiber growth stretches are the same and radial growth stretch is different

%radial direction        
    Fgi_R(3,3)= logical(st(2)>0)*( (f_cc_max) /(1+exp(-c_fpos*(st(2)-st_50pos)) )+1 ) + ...
              logical(st(2)<0) *( -f_cc_max/(1+exp(c_fneg*(st(2)+st_50neg)) )+1 );  %eqn 9   
end
    Fg_R_new = Fgi_R * Fg_old_R(:,:,iteration_growth); 

  %Calculate the new unloaded geometry  - note that the reference is the original unloaded state  
  r0_new(1,2)=r0_s(iteration_growth,2)*Fg_R_new(1,1);
  h0_new(1,2)=h0_s(iteration_growth,2)*Fg_R_new(3,3);
  
end
