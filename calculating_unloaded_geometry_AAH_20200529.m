function [Fg, Fg_R, st, sl, r0, h0]= calculating_unloaded_geometry_AAH_20200529(LV_param, RV_param, h0_s, GrowthTime)
% calculating_unloaded_geometry: determine the size of the unloaded LV and RV geometry

% Inputs:
%   LV_param: determines end-systolic and end-diastolic pressure-volume relationships for the LV
%   RV_param: determines end-systolic and end-diastolic pressure-volume relationships for the RV
%   h0_s: somatic unloaded thickness
%   GrowthTime: how many days to run the simulation

% Outputs:
%   Fg: LV growth tensor
%   Fg_R: RV growth tensor
%   st: thickening stimulus
%   sl: lengthening stimulus
%   r0: unloaded radius
%   h0: unloaded thickness

%Initializing Growth 
iteration_growth_max=GrowthTime;
Fg=zeros(3,3,iteration_growth_max); Fg(:,:,1)=eye(3,3); 
Fg_R=zeros(3,3,iteration_growth_max); Fg_R(:,:,1)=eye(3,3); 

% will eventually need separate stimuli for LV and RV
st=zeros(iteration_growth_max, 2); sl=zeros(iteration_growth_max,2);

%Determining r0
V0 = [LV_param(1,4), RV_param(1,4)];
r0=zeros(iteration_growth_max,2);
r0(1,1)= (3/(4*pi())*V0(1)).^(1/3) ; %determined such that r0 =cuberoot(V0*0.75*(1/pi))
r0(1,2)= (3/(4*pi())*V0(2)).^(1/3) ; %determined such that r0 =cuberoot(V0*0.75*(1/pi))

%Determining h0
h0=zeros(iteration_growth_max,2); %normal compartment, total LV, infarct compartment

% Wall volume is currently constant between LV and RV, but that will likely change
h0(1,1) = h0_s(1,1);
h0(1,2) = h0_s(1,2);

end