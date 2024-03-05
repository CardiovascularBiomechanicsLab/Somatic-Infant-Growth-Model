function [Volumes, TimeVector] =set_initialconditions_noventricles_AAH_20200603(  times, BV, Vint)
%set_initialconditions_noventricles: set circulation parameters and initial circulation volumes

% Inputs
%   times: array containing time to end systole for ventricles and atria, time delay between ventricles and atria, and heart rate
%   BV: current value of stressed blood vlume
%   Vint: proportion of volumes throughout the circulatory model

% Outputs:
%   Volumes: initial compartment volumes
%   TimeVector: cardiac cycle times

%Initialize Timing
HR = times(end);
nRows = 5000; %number of time steps throughout a single cardiac cycle
TimeVector=linspace(0, 60/HR, nRows)'; %setting time vector to length of a single beat (s)
    
%Initialize Volume
Vtotal = BV;  %total effective or stressed blood volume (ml)
Volumes = zeros(nRows,12);   
    

%Estimating initial Volumes
%Weighting compartments by average literature blood volume
Volumes(1,:)=Vint*Vtotal; 
  

end