function [LV_params_new,RV_params_new]=Recomputing_LVRV_params(r0, h0,  a, b, ees, HLHS)
% Recomputing_LVRV_params: computing the altered LV pressure-volume behavior assuming no change 
% (or prescribed changes) in the LV hoop stress - circumferential strain
% behavior

% Inputs:
%   r0: unloaded radius
%   h0: unloaded radius
%   a, b, ees: material parameters that govern the nonlinear relationships between circumferential and radial stretches and hoop stress at ED and ES
%   HLHS: flag for hypoplastic left heart syndrome

% Outputs:
%   LV_params_new: new end-diastolic and end-systolic parameters for the left ventricle 
%   RV_params_new: new end-diastolic and end-systolic parameters for the right ventricle 

% Since we treated the ventricle as a thin-walled spherical pressure vessel the relationship between hoop stress and LV volume at end-systole and end-diastole are 
% sigma_hoop,ED=P_ED*r_ED/(2h_ED ) and sigma_hoop,ES=P_ES*r_ES/(2h_ES )   or
% sigma_hoop,ED=r_ED/(2h_ED )*B*(exp[A*(V_ED-V_0 )]-1) and sigma_hoop,ES=r_ES/(2h_ES )*E*(V_ES-V_0 )

% If we reframe this equation using the unloaded radius and thickness then
% sigma_hoop,ED=r_ED/r0*h0/h_ED*b*(exp[a*(r_ED/r0)^3-a]-1) and sigma_hoop,ES=r_ES/r0*h0/h_ES*e*((r_ES/r0)^3 -1)

% If r0 and h0 are the current unloaded dimensions then this is hoop stress as a function of elastic stretch since Feg = Fe * 1/Fg. 
% In the spherical context Feg = r/r0_initial, thus Fe = r/r0_initial *1/Fg_current.
% But, Fg_current = r0_current/r0_initial so Fe = r/r0_initial * r0_initial/r0_current = r/r0_current
% Thus, in order to calculate the elastic stretch we just divid the loaded radius by the current unloaded radius
% Then the material constants a, b, and e can be used to determine the new
% values of A, B, and E when the unloaded LV grows


%Unloaded LV Volume (ml)
V0=4/3*pi()*r0(1)^3;
LV_params_new(4)=V0;

% Unloaded RV Volume (ml)
V0=4/3*pi()*r0(2)^3;
RV_params_new(4)=V0;

%LV exponential constant in EDPVR  (1/ml)
A=a(1)*3/(4*pi())*1/r0(1)^3;
LV_params_new(1)=A;

%RV exponential constant in EDPVR  (1/ml)
A=a(2)*3/(4*pi())*1/r0(2)^3;
RV_params_new(1)=A;

%LV End-systolic elastance (mmHg/ml)
EES=3/(2*pi())*ees(1)*h0(1)*(1/r0(1)^4);
LV_params_new(3)=EES;

%RV End-systolic elastance (mmHg/ml)
EES=3/(2*pi())*ees(2)*h0(2)*(1/r0(2)^4);
RV_params_new(3)=EES;

%LV linear constant in EDPVR (mmHg)
B=b(1)*2*h0(1)/r0(1);
LV_params_new(2)=B;

%RV linear constant in EDPVR (mmHg)
B=b(2)*2*h0(2)/r0(2);
RV_params_new(2)=B;

if HLHS == 1
    LV_params_new(2) = 0;
    LV_params_new(3) = 0;
end

end

