clear
clc
close all
tic
%% 12 compartment growth model, healthy infant 0-3 yrs
% Stable somatic and coarc model

% CURRENT PROGRESS: 
%      29 May  2020: Converted to 12 compartments up to computing_dimensionsandhemodynamics
%      01 June 2020: Figure out how to get LV/RV parameters working
%      09 June 2020: Regular growth working, get HLHS setting functional
%      22 June 2020: Non-HLHS fetus has zero pressure in LV and RV 
%      25 June 2020: All settings appear to be working for growth and non growth settings
%      08 July 2020: Better method of altering homeostatic setpoints implemented
%      02 Nov  2020: Begin updating inputs to use infant values
%                       - Changed set_initialconditions to reflect new
%                         compartment splitting
%      08 Dec 2020: generateInfantInputs.m created to automatically generate input set based on birth weight
%      05 Jan 2021: iteration_growth changed to start on Day 0
%                       - General clean up and commenting
%      29 Mar 2021: Updated input parameters
%      3 May  2021: Changed plotting functions to plot at end of program 
%      17 May 2021: Updated generateInfantInputs to lower SVR and PVR with new B values
%      19 June 2021: b now using a vector, changes in early timepoints, bDay0 = 0.12
%      6 July 2021: Day0 SVR now determined by Day0 weight
%      30 July 2021: Added coarctation of the aorta
%      12 Dec 2021: Code clean up, FINAL SOMATIC GROWTH+COARC VERSION
%      05 Apr 2022: Removed 8 and 10 compartment options, fetal option, PDA compartment
%      16 Aug 2022: Fixed Cap and Cap_scaled in generateInfantInputs
%      03 Nov 2022: Adjusted SVR, PVR, Cap, and Cas to compensate for fixing 16 Aug error
%      04 Mar 2024: Prepared for upload to Github

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

% Growth Law Constants
%Fit to PO and VO fitting simulations by CMW; see "Predicting the time
%course of ventricular dilation and thickening using a rapid compartmental
%model" in Journal of Cardiovascular Translational Research
f_f=31;  
sl_50=0.215; 
c_f_neg=553.5130000000000*1.04; 
negst_50=0.034712950000000*0.987; 
c_f_pos=36.42; 
posst_50=0.0971; 
growthparams=[f_f sl_50 c_f_pos c_f_neg posst_50 negst_50];

%% RUNNING CONTROL CIRCULATION LOOPS

weight_percentile = 50; % 50th, 10th, 90th percentiles available (also 'C' and 'D' for Single Ventricle Reconstruction Trial HLHS group data)
if weight_percentile == 50
    load('inputData/20221103_somaticSetpoints_50.mat');
elseif weight_percentile == 10
    load('inputData/20221103_somaticSetpoints_10.mat');
elseif weight_percentile == 90
    load('inputData/20221103_somaticSetpoints_90.mat');
elseif weight_percentile == 'C'
    load('inputData/20221104_somaticSetpoints_C.mat');
elseif weight_percentile == 'D'
    load('inputData/20221104_somaticSetpoints_D.mat');
else
    disp('Weight percentile not correctly input!'); return;
end

somaticSetpoints = setpoints;
iteration_growth=0; % growth step, baseline is now 0 (day of birth) as perturbation step is gone
GrowthTime = 1095; % how long to run the simulation - 1095 days = 3 years
Coarc = 0; % 0 -no, 1 -yes % Coarctation of the aorta
HLHS = 0; %0 - no, 1 - yes   % Experimental
Stage = [0 0 0]; % Surgical Stages involved in palliation of HLHS: 0 if has not occurred, 1 if occurred, Stage(1) is Norwood/Sano, Stage(2) is Glenn, Stage(3) is Fontan
if Coarc == 0 & HLHS == 0 % if this is a somatic simulation
    generateSetpoints = 1; % save new somatic setpoint .mat file (**NOTE: MUST be 1 for any somatic simulations!)
else % this is a pathologic simulation
    generateSetpoints = 0; % use previously loaded somatic setpoints
end
saveOutputs = 1; % 0 - no, 1 - yes % save model output .mat file

% pre-allocate cell arrays
PressureArray = cell(1,GrowthTime);
VolumeArray = cell(1,GrowthTime);
rArray = cell(1,GrowthTime);
hArray = cell(1,GrowthTime);
TimeVectorArray = cell(1,GrowthTime);

% generate set of input parameters - includes capacitances, resistances, and ventricular and atrial parameters
inputSet = generateInfantInputs_AAH_20221103(GrowthTime,weight_percentile,HLHS);

% somatic unloaded radii and thicknesses
r0_s = inputSet(:,30);
r0_s(:,2) = inputSet(:,31);
h0_s = inputSet(:,32);
h0_s(:,2) = inputSet(:,33);

% organize input parameters
[capacitances, resistances, LVparameters, RVparameters, LAparameters, RAparameters, times, SBV_s, ratio_LBR] = set_initialconditions_AAH_20211217(HLHS, Stage, inputSet, Coarc);

% calculate initial volumes and pressures
[Volumes, Pressures, TimeVector] = initializecycle_AAH_20200529(times(1,:), SBV_s(1), HLHS, Stage);

% solve for Day 0 volumes and pressures
[Volumes, Pressures, Valves, iES(1)] = RunCirculationModel_AAH_20200529(Volumes, Pressures, TimeVector, capacitances(1,:), resistances(1,:), ...
                                    LVparameters(1,:), RVparameters(1,:), LAparameters(1,:), RAparameters(1,:), times(1,:), HLHS, Stage)   ;
Baseline_Volumes=Volumes; Baseline_Pressures=Pressures; Baseline_TimeVector=TimeVector; 

PressureArray{1} = Pressures;
VolumeArray{1} = Volumes;
TimeVectorArray{1} = TimeVector;

%% INITIALIZING GROWTH

%determine the size of the unloaded LV and RV geometry
[Fg, Fg_R, st, sl, r0, h0]= calculating_unloaded_geometry_AAH_20200529(LVparameters(1,:), RVparameters(1,:), h0_s, GrowthTime);

%compute the radius and thickness of the LV and RV over the cardiac cycle
[r, h]=dimensionsofsphere_AAH_20200529(Baseline_Volumes(:,3), Baseline_Volumes(:,10), r0(1,:), h0(1,:));

rArray{1} = r;
hArray{1} = h;

% pull a b ees from inputSet
a = inputSet(:,34); a(:,2) = a(:,1);
b = inputSet(:,35); b(:,2) = b(:,1);
ees = inputSet(:,36); ees(:,2) = ees(:,1);

% calculate the stress-strain relationships of the LV and RV
[sigma_ED(:,1,:), sigma_ES(:,1,:), strainc(:,1,:)]= calculate_stress_strain_relationships_AAH_20200529(LVparameters, RVparameters, r0(1,:), h0(1,:));
Baseline_r=r; Baseline_h=h;

%  Initial values of:  
%      dimensions: LV (1)EDV, (2)ESV, (3)ED thickness, (4)ES thickness (thicknesses converted to mm)
%      ValuesofInt: (1)EDP (2)EDV (3)MaxP (4)MAP (5)maxV (6)minV (7)RF (8)CO (9)HR (10)ESP (11)ESP
[dimensions, ValuesofInterest, dimensions_R, ValuesofInterest_R]=computing_dimensionsandhemodynamics_AAH_20200609(Baseline_Pressures, Baseline_Volumes, TimeVector, Valves,  r, h,  iteration_growth,  GrowthTime, HLHS);

iteration_growth=1;

%% IMPLEMENTING GROWTH

Vint= Volumes(1,:)./sum(Volumes(1,:)); %new proportion of volumes throughout the circulatory model

     while iteration_growth <= GrowthTime
         %Calculating the grown unloaded LV dimensions
         [r0(iteration_growth+1,:), h0(iteration_growth+1,:), Fg(:,:, iteration_growth+1), Fg_R(:,:, iteration_growth+1), st(iteration_growth+1,:), sl(iteration_growth+1,:), ...
             fiber_strain(iteration_growth+1,:), radial_strain(iteration_growth+1,:), radial_strain_baseline(iteration_growth+1,:),setpoints(iteration_growth,:)]=GrowthLaw_modifiedKOM_AAH_20210708(r, h, ...
             Fg, Fg_R, growthparams, HLHS, r0_s, h0_s, iteration_growth, somaticSetpoints(iteration_growth,:), generateSetpoints);
         
         %Computing the altered LV pressure-volume parameters
        [LVparameters(iteration_growth+1,:),RVparameters(iteration_growth+1,:)]=Recomputing_LVRV_params(r0(iteration_growth+1,:), h0(iteration_growth+1,:),  ...
             a(iteration_growth+1,:),b(iteration_growth+1,:), ees(iteration_growth+1,:), HLHS); 
         
         %Adjusting "Stressed" Blood Volume      
         %The "Stressed" Blood Volume within the compartmental model actually includes some unstressed volume since the unloaded volumes from the LV and RV are incorporated in this number.
         %Thus, when the unloaded LV volume is increased or decreased with growth the total "Stressed" blood volume parameter within the circulation model should change   
         if Coarc == 1
            SBV_coarc = (0.4302*ratio_LBR(iteration_growth+1)+0.5674)*SBV_s(iteration_growth+1); 
         else
            SBV_coarc = SBV_s(iteration_growth+1);
         end
         
         SBV(iteration_growth+1) = SBV_coarc + ((4/3)*pi)*(r0(iteration_growth+1,1)^3 - r0_s(iteration_growth+1,1)^3 + r0(iteration_growth+1,1)^3 - r0_s(iteration_growth+1,1)^3);
                
         %Re-loading the grown ventricle 
         iteration_growth=iteration_growth+1;
         if iteration_growth<GrowthTime+1
             %initialize the circulatory parameters
             [InitialVolumes, TimeVector] = set_initialconditions_noventricles_AAH_20200603( times(iteration_growth,:), SBV(iteration_growth), Vint); 
              TimeVectorArray{iteration_growth} = TimeVector;
             %run the circulatory model    
             [Volumes, Pressures, Valves, iES(iteration_growth)] = RunCirculationModel_AAH_20200529(InitialVolumes, Pressures, TimeVector, capacitances(iteration_growth,:), resistances(iteration_growth,:), ...
                 LVparameters(iteration_growth,:), RVparameters(iteration_growth,:), LAparameters(iteration_growth,:), RAparameters(iteration_growth,:), times(iteration_growth,:), HLHS, Stage)   ;
                 
             PressureArray{iteration_growth} = Pressures;
             VolumeArray{iteration_growth} = Volumes;
                
             %compute the radius and thickness of the LV over the cardiac cycle
             [r, h]=dimensionsofsphere_AAH_20200529(Volumes(:,3), Volumes(:,10), r0(iteration_growth,:), h0(iteration_growth,:));
             [dimensions(:,iteration_growth), ValuesofInterest(:,iteration_growth), dimensions_R(:,iteration_growth), ValuesofInterest_R(:,iteration_growth)]=computing_dimensionsandhemodynamics_AAH_20200609(Pressures, Volumes, TimeVector, Valves,  r, h,  iteration_growth,  GrowthTime, HLHS);
              
             rArray{iteration_growth} = r;
             hArray{iteration_growth} = h;
             
             
             %compute the stress-strain relationships of the LV. Note here we are NOT resetting the material properties.
             [sigma_ED(:,iteration_growth,:), sigma_ES(:,iteration_growth,:), strainc(:,iteration_growth,:)]= calculate_stress_strain_relationships_AAH_20200529(LVparameters(iteration_growth,:), RVparameters(iteration_growth,:), r0(iteration_growth,:), h0(iteration_growth,:));

             Vint= Volumes(1,:)./sum(Volumes(1,:)); 
         end % circulation loop

     end  % growth loop
  
     if generateSetpoints == 1
        time = datestr(now, 'yyyymmdd_HHMM');
        id = string(weight_percentile);
        filename = sprintf('%s_somaticSetpoints_%s.mat',time,id);
        save(filename,'setpoints');
     end
  
  if saveOutputs == 1
    time = datestr(now, 'yyyymmdd_HHMM');
    id = string(weight_percentile);
    filename = sprintf('%s_output_%s.mat',time,id);
    save(filename,'VolumeArray','PressureArray','LVparameters','RVparameters','rArray','hArray','r0','h0','r0_s','h0_s','TimeVectorArray','Fg', 'Fg_R', 'st', 'sl', 'growthparams','ValuesofInterest','ValuesofInterest_R','dimensions','dimensions_R','SBV');
  end
 
toc


