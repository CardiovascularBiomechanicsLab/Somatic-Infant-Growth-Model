function [sigma_ED, sigma_ES, strainc]= calculate_stress_strain_relationships_AAH_20200529(LV_param, RV_param, r0, h0)
% calculate_stress_strain_relationships: determinines the ED and ES stress-strain relationships of the myocardium

% Inputs:
%   LV_param: determines end-systolic and end-diastolic pressure-volume relationships for the LV
%   RV_param: determines end-systolic and end-diastolic pressure-volume relationships for the RV
%   r0: unloaded radius
%   h0: unloaded thickness

% Outputs:
%   sigma_ED: stress at end diastole
%   sigma_ES: stress at end systole
%   strainc: circumferential strain
 
%Use the unloaded dimensions for the remote healthy ventricle
h0_val=[h0(1,1),h0(1,2)]; % [LV, RV]
r0_val=[r0(1,1),r0(1,2)]; % [LV, RV]
  
  
%Inflate the Ventricle and Calculate ES and ED pressure
V_test = linspace(0.01,200,150000);
r_test(1,:) = (0.75*V_test*(1/pi())).^(1/3);  %LV
r_test(2,:) = (0.75*V_test*(1/pi())).^(1/3); %RV
h_test(1,:)= (  h0_val(1)^3+3* h0_val(1)^2* r0_val(1)+3* h0_val(1)* r0_val(1)^2 + r_test(1,:).^3).^(1/3) - r_test(1,:); %ischoric, page 26 CMW notebook 5  LV
h_test(2,:)= (  h0_val(2)^3+3* h0_val(2)^2* r0_val(2)+3* h0_val(2)* r0_val(2)^2 + r_test(2,:).^3).^(1/3) - r_test(2,:); % RV
PES(1,:) = (  LV_param(1,3)*(V_test-LV_param(1,4))  ); %ESPVR, LV
PES(2,:) = (  RV_param(1,3)*(V_test-RV_param(1,4))  ); % RV
PED(1,:) = (  LV_param(1,2)* ( exp(  LV_param(1,1)*(V_test-LV_param(1,4))  )  -1)  ); %EDPVR, LV
PED(2,:) = (  RV_param(1,2)* ( exp(  RV_param(1,1)*(V_test-RV_param(1,4))  )  -1)  ); % RV
        
%Determine the hoop stress to circumferential strain relationship assuming a thin walled sphere
sigma_ES(1,:,1) = PES(1,:).*r_test(1,:)./(2*h_test(1,:)); % LV
sigma_ES(1,:,2) = PES(2,:).*r_test(2,:)./(2*h_test(2,:)); % RV
sigma_ED(1,:,1) = PED(1,:).*r_test(1,:)./(2*h_test(1,:)); % LV
sigma_ED(1,:,2) = PED(2,:).*r_test(2,:)./(2*h_test(2,:)); % rV
stretch(1,:)=r_test(1,:)/r0_val(1); % LV
stretch(2,:)=r_test(2,:)/r0_val(2); % RV
strainc(1,:,1)=(0.5*(stretch(1,:).^2-1))'; % LV
strainc(1,:,2)=(0.5*(stretch(2,:).^2-1))';
sigma_ED=permute(sigma_ED,[2 1 3]); % to 150000x1x2
sigma_ES=permute(sigma_ES,[2 1 3]);
strainc=permute(strainc,[2,1,3]);
         

end