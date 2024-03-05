function [Volumes, Pressures, Valves, iES]=RunCirculationModel_AAH_20200529(Volumes, Pressures, TimeVector, capacitances, resistances, LVparam, RVparam, LAparam, RAparam, times, HLHS,  Stage)

%runcirculation: computes the pressure and volume of each compartment throughout the entire cardiac cycle   

   %Inputs: 
   %       Volumes - volumes in each compartment initially
   %       Pressures - pressures in each compartment initially
   %       TimeVector - vector of size nRows giving the time in the cardiac cycle in seconds
   %       capacitances - capacitances for systemic and pulmonary veins and arteries
   %       resistances - systemic and pulmonary venous and pulmonary resistances as well as characteristic
   %                     resistances and shunt resistances 
   %       LV_param - end-diastolic and end-systolic parameters for the left ventricle 
   %       RV_param - end-diastolic and end-systolic parameters for the right ventricle 
   %       LA_param - end-diastolic and end-systolic parameters for the left atria
   %       RA_param - end-diastolic and end-systolic parameters for the right atria
   %       times - vector containing the timing of chamber activation, AV delay, and heart rate
   %       HLHS - flag for hypoplastic left heart syndrome
   %       Stage - HLHS surgical stage flag
   
   %Outputs:
   %       Volumes - volumes in each compartment throughout the cardiac cycle
   %       Pressures - pressures in each compartment throughout the cardiac cycle
   %       Valves - logical array showing whether valves are open or closed
   %       iES - index of end systole in TimeVector

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
    
%% Initialize Cycle
isTransientState = 1; cutoff=0.0001; nIterations=0; absoluteError=100;
Tes = times(2);
Pressures(1,:) = 0;
Volumes(end,:)=Volumes(1,:);
%% Loop Through Cardiac Cycle

while isTransientState
            nIterations = nIterations + 1;           
            Volumes(1,:)=Volumes(end,:);
            
            
            Pressures(1,:)=calculate_pressures_AAH_20200529(Volumes(1,:), capacitances, LVparam, RVparam, LAparam, RAparam, 0, times);
            
            %% ITERATIVELY SOLVE FOR VOLUME AND PRESSURE THROUGHOUT THE CARDIAC CYCLE
            for i=2:length(TimeVector)  
                currentVolumes=Volumes(i-1,:);
                currentTime=TimeVector(i,:);
                currentPressures=Pressures(i-1,:);
                [ Volumes(i,:)] = rk4_AAH_20220405(currentVolumes, currentPressures, resistances, TimeVector, HLHS, Stage) ;
                [Pressures(i,:), Ve(i), Ae(i)]=calculate_pressures_AAH_20200529(Volumes(i,:), capacitances, LVparam, RVparam, LAparam, RAparam, currentTime, times);
            end
            clear i
 
            %% Calculate if Steady-State has Occured
            absoluteError = abs(Volumes(end,:) - Volumes(1,:));
            outOfRangeLogical = absoluteError > cutoff;
            outOfRangeIndices = find(outOfRangeLogical);
            nOutOfRange = numel(outOfRangeIndices);
            isTransientState = nOutOfRange >= 1;
            str = sprintf(['SS iteration #', num2str(nIterations), ...
            '. \t\tAbsolute error: ', num2str(absoluteError)]);
            disp(str);
            if nIterations>99
                break
                fprintf('ERROR: SS iteration # > 100, check cutoff')
            end

         
            
           
end

%Vector showing whether the valves are open or closed
    
    Valves=zeros(size(Pressures,1),4);
    %Valves (MV, AoV , TV, PV)
    Valves(:,1) = logical(Pressures(:,2) > Pressures(:,3)); %PLA > P_LV
    Valves(:,2) = logical(Pressures(:,3) > Pressures(:,4)); %P_LV > Paorta - check with Colleen
    Valves(:,3) = logical(Pressures(:,9) > Pressures(:,10)); %P_RA > P_RV 
    Valves(:,4) = logical(Pressures(:,10) > Pressures(:,11)); %P_RV > P_PA 
     
    %Calculating the time and pressure at maximum contraction
    [tmp, iES]=min(abs(TimeVector-Tes/1000));
    PES=Pressures(iES,3);
    PES_R=Pressures(iES,10);

  
end