function [Volumes, Pressures, TimeVector]=initializecycle_AAH_20200529(times, BV, HLHS, Stage)
%intializecycle: initializes ventors to track the time, volume, and pressure throughout the cardiac cycle

    %Inputs: 
    %       times - vector containing the timing of chamber activation, AV delay, and heart rate
    %       BV - stressed blood volume (ml)
    %       HLHS - HLHS (single ventricle) circulation
    %       Stage - which stage of HLHS surgery has been completed

    %Outputs:
    %       TimeVector - vector of size nRows giving the time in the cardiac cycle in seconds
    %       Volumes - the volume in each compartment at a time step during the cardiac cycle
    %       Pressures - the pressure in each compartment at a time step during the cardiac cycle
   
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

    
nRows = 5000; %number of time steps throughout a single cardiac cycle   
   
HeartRate=times(end);
TimeVector=linspace(0, 60/HeartRate, nRows)'; %setting time vector to length of a single beat (s)

Volumes=zeros(nRows, 12); 
Volumes(1,:)=BV/12; %there are 12 compartments

 if HLHS==1
     Volumes(1,:)=BV/11; %the left ventricle does not exist
     Volumes(1,3)=0; 
 elseif sum(Stage)>0
     error('Input for surgical status is incorrect! A patient without HLHS should have zeros for all stages'); 
 end

% initialize pressure array
Pressures=zeros(nRows, 12); 

end
