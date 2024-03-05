function [r, h]=dimensionsofsphere_AAH_20200529(LV_Volume, RV_Volume, r0, h0)

% Inputs:
%   LV_Volume: volume array for the LV
%   RV_Volume: volume array for the RV
%   r0: unloaded radius
%   h0: unloaded thickness

% Outputs:
%   r: loaded radius
%   h: loaded thickness


%Calculating the Dimensions of the LV and RV
r=zeros(size(LV_Volume,1),2); h=zeros(size(LV_Volume,1),2); %intializing radius and thickness 

%Radius of the LV
r(:,1)=(0.75*LV_Volume*(1/pi())).^(1/3); 
h(:,1)=( h0(1,1)^3+3*h0(1,1)^2*r0(1,1)+3*h0(1,1)*r0(1,1)^2 + r(:,1).^3).^(1/3)-r(:,1); %ischoric

% Radius of the RV
r(:,2)=(0.75*RV_Volume*(1/pi())).^(1/3); 
h(:,2)=( h0(1,2)^3+3*h0(1,2)^2*r0(1,2)+3*h0(1,2)*r0(1,2)^2 + r(:,2).^3).^(1/3)-r(:,2); %ischoric



   
end
