function [Pressures, V_epsilon, A_epsilon]=calculate_pressures_AAH_20200529(Volumes, capacitances, LVparameters, RVparameters, LAparameters, RAparameters, currentTime, times)
%%calculate_pressures: calculates the current pressure in each vessel bed

   %Inputs: 
   %       capacitances - capacitances for systemic and pulmonary veins and arteries 
   %       LV_parameters - end-diastolic and end-systolic parameters for the left ventricle 
   %       RV_parameters - end-diastolic and end-systolic parameters for the right ventricle 
   %       LA_parameters - end-diastolic and end-systolic parameters for the left atria
   %       RA_parameters - end-diastolic and end-systolic parameters for the right atria
   %       times - vector containing the timing of chamber activation, AV delay, and heart rate
   %       currentTime - current time in the cardiac cycle
   %       Volumes - volumes in each compartment at the current time
   
   %Outputs:
   %       Pressures - pressures in each compartment at the current time
   %       V_epsilon - time-varying elastance function for ventricles
   %       A_epsilon - time-varying elastance function for atria
   
   %Volume and Pressure Compartment Designations
    %   Column  Compartment    
    %   1       pulmonary veins
    %   2       left atrium
    %   3       left ventricle
    %   4       aorta
    %   5       systemic arteries - lower body
    %   6       systemic arteries - upper body
    %   7       systemic veins - lower body
    %   8       systemic veins - upper body
    %   9       right atrium
    %   10      right ventricle
    %   11      main pulmonary artery                  
    %   12      pulmonary arteries      

%%  Calculate Pressures in Circulation
Pressures(1) = Volumes(1)/capacitances(1); 
Pressures(4) = Volumes(4)/capacitances(4); 
Pressures(5) = Volumes(5)/capacitances(5); 
Pressures(6) = Volumes(6)/capacitances(6); 
Pressures(7) = Volumes(7)/capacitances(7);      
Pressures(8) = Volumes(8)/capacitances(8); 
Pressures(11) = Volumes(11)/capacitances(11);       
Pressures(12) = Volumes(12)/capacitances(12); 

%% Calculate Heart Chamber Timing
Tes=times(1:2);  %time to end-systole for atria and ventricles (ms)
AV_delay=times(3); %delay between atrial and ventricular contraction (ms)
HR=times(end); %heart rate (beats/min)

Tes=Tes/1000; %tau=tau/1000; %converting to seconds 

% Calculate the time-varying weight factor (epsilon) for ventricles
if currentTime < 2*Tes(2) 
    V_epsilon = 1 / 2 *  sin(  pi()/Tes(2)* currentTime - pi()/2 )   + 1/2;
else
    V_epsilon = 0; 
end

% Calculate the time-varying weight factor (epsilon) for atria
t_atria = currentTime - (Tes(2)-Tes(1)) + AV_delay/1000; % atria timing offset from ventricles
t_atria=t_atria-(t_atria>( 60/HR))*(60/HR);

if t_atria < 2*Tes(1) 
    A_epsilon = 1 / 2 *  sin( pi()/Tes(1)*t_atria - pi()/2 )   + 1/2;
else
    A_epsilon = 0; 
end
     
%% Calculate Heart Chamber Pressures

%Calculate RV Pressure
RV_ESP = RVparameters(1,3) * (Volumes(10) - RVparameters(1,4));
RV_EDP = RVparameters(1,2) * ( exp(RVparameters(1,1)*(Volumes(10)-RVparameters(1,4))) - 1 );
Pressures(10) = V_epsilon * (RV_ESP - RV_EDP) + RV_EDP;

% Calculate LV Pressure 
LV_ESP = LVparameters(1,3) * (Volumes(3) - LVparameters(1,4));
LV_EDP = LVparameters(1,2) * ( exp(LVparameters(1,1)*(Volumes(3)-LVparameters(1,4))) - 1 );
Pressures(3) = V_epsilon * (LV_ESP - LV_EDP) + LV_EDP;

%Calculate RA Pressure
RA_ESP = RAparameters(1,3) * (Volumes(9) - RAparameters(1,4));
RA_EDP = RAparameters(1,2) * ( exp(RAparameters(1,1)*(Volumes(9)-RAparameters(1,4))) - 1 );
Pressures(9) = A_epsilon * (RA_ESP - RA_EDP) + RA_EDP;

% Calculate LA Pressure 
LA_ESP = LAparameters(1,3) * (Volumes(2) - LAparameters(1,4));
LA_EDP = LAparameters(1,2) * ( exp(LAparameters(1,1)*(Volumes(2)-LAparameters(1,4))) - 1 );
Pressures(2) = A_epsilon * (LA_ESP - LA_EDP) + LA_EDP;


end