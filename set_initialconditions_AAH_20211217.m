function [capacitances, resistances, LVparameters, RVparameters, LAparameters, RAparameters, times, SBV, ratio_LBR] = set_initialconditions_AAH_20211217(HLHS, Stage, inputSet, Coarc)
%Inputs: 
   %       HLHS - %0 - no, 1 - yes   HLHS status
   %       Stage -  %0 if has not occurred, 1 if occurred; Stage(1) is Norwood/Sano, Stage(2) is Glenn, Stage(3) is Fontan
   %       inputSet - array of model inputs
   %       Coarc - 0 - no, 1 - yes CoA status
   
   %Outputs:
   %       capacitances - capacitances for systemic and pulmonary veins and arteries
   %       resistances - systemic and pulmonary venous and pulmonary resistances as well as characteristic
   %                     resistances and shunt resistances 
   %       LV_parameters - end-diastolic and end-systolic parameters for the left ventricle 
   %       RV_parameters - end-diastolic and end-systolic parameters for the right ventricle 
   %       LA_parameters - end-diastolic and end-systolic parameters for the left atria
   %       RA_parameters - end-diastolic and end-systolic parameters for the right atria
   %       times - vector containing the timing of chamber activation, AV delay, and heart rate
   %       SBV - stressed blood volume (ml)
   %       ratio_LBR - ratio of (normal) scaled lower body resistance to coarc values
   
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
 
%% Import Parameter File

A = inputSet(:,1);
Aatria = inputSet(:,2);
AVdelay = inputSet(:,3);
B = inputSet(:,4);
Batria = inputSet(:,5);
Cap = inputSet(:,6);
Cap_scaled = inputSet(:,7); 
Cas = inputSet(:,8);
Cvp = inputSet(:,9);
Cvs = inputSet(:,10);
Ees = inputSet(:,11); 
HR = inputSet(:,12);
Rap = inputSet(:,13); % pulmonary arterial resistance, also known as PVR
ratio_AotoSA = inputSet(1,14); % ratios do not change throughout time
ratio_LBtoUB_Arteries = inputSet(1,15); 
ratio_LBtoUB_Veins = inputSet(1,16);
ratio_MPAtoPA = inputSet(1,17);
Rcp = inputSet(:,18);
Rcs = inputSet(:,19);
Rmv = inputSet(:,20); 
Rtv = inputSet(:,21); 
Rvp = inputSet(:,22); 
Rvs = inputSet(:,23);
SBV = inputSet(:,24);
SVR = inputSet(:,25); 
Tmax(:,1) = inputSet(:,26);
Tmax(:,2) = inputSet(:,27);
V0 = inputSet(:,28);
Rap_scaled = inputSet(:,29);
V0_atria = inputSet(:,37);


%% Resistances

resistances=zeros(size(inputSet,1),12); %initialize vector

%Not changed by adding aorta and MPA or splitting systemic body into upper and lower
resistances(:,1)  = Rvp; %Rvp %plumonary venous resistance - resists return of blood to heart (mmHg * s/mL) Rvp 
resistances(:,2)  = Rmv; %Rmv %resistance of mitral valve (mmHg * s/mL) 
resistances(:,3)  = Rcs; %Rcs %systemic characteristic resistance (mmHg * s/mL)     
resistances(:,9)  = Rtv; %Rtv %resistance of tricuspid valve (mmHg * s/mL) 
resistances(:,10) = Rcp; %Rcp %pulmonary characteristic resistance (mmHg * s/mL);

% resistances altered by compartment splitting
resistances(:,6)=inf; %systemic arterial resistance (mmHg * s/mL) - will become upper body 
resistances(:,8)=inf; %systemic venous resistance (mmHg * s/mL) - will become upper body
resistances(:,7)=Rvs; %Rvs %systemic venous resistance (mmHg * s/mL) - will become lower body
resistances(:,5)=SVR; %SVR %systemic arterial resistance (mmHg * s/mL) - will become lower body; 
resistances(:,4)=inf;        %will become aortic resistance (mmHg * s/mL) 
resistances(:,11)=inf;      %will become MPA resistance (mmHg * s/mL) 
resistances(:,12)=Rap_scaled; %Rap %pulmonary arterial resistance (mmHg * s/mL) - will exclude MPA
   
%add aorta and MPA
R_total=resistances(:,5);
resistances(:,5)=R_total*ratio_AotoSA^(3/4) / (1 + ratio_AotoSA^(3/4)); %SA
resistances(:,4)=resistances(:,5)/ratio_AotoSA^(3/4); %Aorta  
R_total=resistances(:,12);
resistances(:,12)=R_total*ratio_MPAtoPA^(3/4) / (1 + ratio_MPAtoPA^(3/4));  %PA
resistances(:,11)=resistances(:,12)/ratio_MPAtoPA^(3/4); %MPA 
% To change PVR after split, need to set resistances(12) to current PVR (apply to distal PA only)
resistances(:,12) = Rap; % custom PVR values

%split systemic arteries and veins into upper and lower body 
% splitting systemic arteries
R_total_SA=resistances(:,5);
resistances(:,6)=R_total_SA*(1+ratio_LBtoUB_Arteries^(-3/4)) / (ratio_LBtoUB_Arteries^(-3/4)); % systemic arteries - upper body
resistances(:,5)=resistances(:,6)*ratio_LBtoUB_Arteries^(-3/4); % systemic arteries - lower body
oldLBresistance = resistances(:,5);

% coarctation of the aorta
if Coarc == 1
    % assuming that resistance of lower body arterial compartment increases
    % with time as surrounding tissue grows but coarctation region remains
    % constricted
    percent_diff(1:37)=0;
    x = 37:245;
    increase = 3.58030573279e-8.*x.^3 - 2.26617600279e-5.*x.^2 + 0.005072046879732.*x - 0.159025716599736;
    percent_diff(38:246) = increase;
    percent_diff(247:1096) = 0.25;
    percent_diff = percent_diff';
    resistances(:,5) = resistances(:,5) + (resistances(:,5).*percent_diff);
 end

% to calculate change in SBV
ratio_LBR = resistances(:,5) ./ oldLBresistance;   

% splitting systemic veins
R_total_SV=resistances(:,7);
resistances(:,8)=R_total_SV*(1+ratio_LBtoUB_Veins^(-3/4)) / (ratio_LBtoUB_Veins^(-3/4)); % systemic veins - upper body
resistances(:,7)=resistances(:,8)*ratio_LBtoUB_Veins^(3/4); % systemic veins - lower body
    
%additional resistances
k=2.66; d=8; %diameters of shunt (mm)
shunt_resistance=1/(k*d^2); 
resistances(:,13)=inf; %resistance of patent ductus arteriosus (PDA)
resistances(:,14)=inf; %resistance of atrial septal defect (ASD), also known as patent formamen ovale (PFO); set to shunt resistance if shunt exists
resistances(:,15)=inf; %resistance of Sano shunt 
resistances(:,16)=inf; %resistance of fenestration between Fontan baffle and right atrium (only present in some patients who have undergone 3rd surgical stage)
resistances(:,(17:20))=inf; %backflow resistances for the mitral, aortic, triscupid, and pulmonary valves

if HLHS==1 %use hypoplastic heart syndrome circulation in which the LV is not present
    resistances(:,2)=inf; resistances(:,3)=inf;
    if Stage(1)==0 %the first surgery has not been completed
            R_ASD = 0.051903364;
            R_PDA = 0.932353168;
            resistances(:,14)=R_ASD;%0.1;%0.005;%shunt_resistance; %there is a ASD %0.0012 lowest bound
            resistances(:,13)=R_PDA;%0.05;%1; %5; %there is a PDA
            
            if sum(Stage)>0
                    error('Input for surgical status is incorrect!'); 
            end
    elseif Stage(1)==1 && Stage(2)==0 && Stage(3)==0 %the first surgery has been completed
            resistances(:,14)=shunt_resistance; %there is a PFO
            resistances(:,15)=  0.46; %0.44; %0.4; %0.48; %0.64; %0.32; %0.16;%0.08; %0.04; %there is a Sano shunt (also called an RV to PA shunt), initially set at same value as pulmonary characteristic resistance 20Mar2019
            resistances(:,11)=resistances(:,11)*2; %the MPA was replaced by the sano shunt, depending on the size of the shunt this value will change
            resistances(:,10)=resistances(:,10)*0.95; %the resistance of the aorta may change with surgery since some of the MPA is typically sewn into the Aorta
                                                % note I do not know whether it typically increases, decreases, or stays the same
            resistances(:,4)=resistances(:,4)*0.95; %the resistance of the aorta may change with surgery since some of the MPA is typically sewn into the Aorta
                                                % note I do not know whether it typically increases, decreases, or stays the same
    elseif Stage(1)==1 && Stage(2)==1 && Stage(3)==0 %the first and second surgeries have been completed
            resistances(:,14)=shunt_resistance; %there is a PFO
            resistances(:,11)=resistances(:,11)*2; %the sano shunt was taken down and thus this resistance may change from post-S1
            resistances(:,4)=resistances(:,4)*0.95; resistances(:,10)=resistances(:,10)*0.95; %should be the same value as after stage 1 
    elseif Stage(1)==1 && Stage(2)==1 && Stage(3)==1 %all three surgeries have been completed
            resistances(:,14)=shunt_resistance; %there is a PFO
            resistances(:,11)=resistances(:,11)*2; %a portion of the RA was used to connect the IVC to the pulmonary arteries which may alter this resistance
            resistances(:,4)=resistances(:,4)*0.95; resistances(:,10)=resistances(:,10)*0.95; %should be the same value as after stage 1
            resistances(:,16)=10; %this is the resistance of the fenestration between Fontan baffle and right atrium, a hole is sometimes left to allow blood to bypass
            %the lungs if the pulmonary pressure gets too high
    end
end


%% Capacitances

capacitances=zeros(size(inputSet,1),12); %initialize vector

%Not changed by adding aorta and MPA or splitting systemic body into upper and lower
%note 2, 3, 9, and 10 are 0 since these are time varying capacitances used for the pumping chambers and are defined elsewhere 
%4, 6, 8, and 11 are 0 for 8 chambers
capacitances(:,1)= Cvp; %Cvp %pulmonary venous compliance (ml/mmHg)
capacitances(:,5)= Cas; %Cas %systemic arterial compliance (ml/mmHg) - will become lower body
capacitances(:,7)= Cvs; %Cvs  %systemic venous compliance (ml/mmHg) - will become lower body
capacitances(:,12)= Cap_scaled; %Cap %pulmonary arterial compliance (ml/mmHg) - will exclude MPA
capacitances(:,4) = inf;
capacitances(:,6) = inf;
capacitances(:,8) = inf;
capacitances(:,11) = inf;

%add aorta and MPA  
Ctotal=capacitances(:,5);
capacitances(:,5)=Ctotal/(1+ratio_AotoSA); %SA 
capacitances(:,4)=capacitances(:,5)*ratio_AotoSA; %Aorta
Ctotal=capacitances(:,12);
capacitances(:,12) = Ctotal/( 1 + ratio_MPAtoPA); %PA
capacitances(:,11) = capacitances(:,12)*ratio_MPAtoPA; %MPA

% Reset distal PA capacitance to maturation function
capacitances(:,12) = Cap;

%split systemic arteries and veins into upper and lower body   
Ctotal_SA=capacitances(:,5);
capacitances(:,6)=Ctotal_SA/(1+ratio_LBtoUB_Arteries); % upper body arteries
capacitances(:,5)=capacitances(:,6).*ratio_LBtoUB_Arteries; % lower body arteries

% adjust capacitance for LB coarctation case
oldLBCapacitance = capacitances(:,5);
capacitances(:,5) = oldLBCapacitance.*(resistances(:,5)./oldLBresistance).^(-4/3);

Ctotal_SV=capacitances(:,7);
capacitances(:,8)=Ctotal_SV./( 1 + ratio_LBtoUB_Veins); % upper body veins
capacitances(:,7)=capacitances(:,8).*ratio_LBtoUB_Veins; % lower body veins

if HLHS==1
    if Stage(1)==1 && Stage(2)==0 && Stage(3)==0 %the first surgery has been completed
        capacitances(:,11)=capacitances(:,11)*0.5; %the MPA was replaced by the sano shunt
        capacitances(:,4)=capacitances(:,4)*1.1; %the aorta was altered to include some of the MPA tissue
        capacitances(:,13)=0; %the PDA no longer exists
    elseif Stage(1)==1  && Stage(2)==1 && Stage(3)==0 %the second surgery has been completed
        capacitances(:,11)=capacitances(:,11)*0.5;   %the sano shunt was taken down and thus this capacitance may change from post-S1
        capacitances(:,4)=capacitances(:,4)*1.1; %should be the same value as after stage 1 
        capacitances(:,13)=0; %the PDA no longer exists
    elseif  Stage(1)==1  && Stage(2)==1 && Stage(3)==1 %all three stages of surgery have been completed
        capacitances(:,11)=capacitances(:,11)*0.5; %a portion of the RA was used to connect the IVC to the pulmonary arteries which may alter this capacitance
        capacitances(:,4)=capacitances(:,4)*1.1; %should be the same value as after stage 1 
        capacitances(:,13)=0; %the PDA no longer exists
    end
end


    
%% Ventricles and Atria
% [A B Ees V0] 
% EDP = B*exp(A*(EDV-V0))-B; ESP = Ees*(ESV-V0)
% A: 1/ml
% B: mmHg
% Ees: mmHg/ml
% V0: ml
% defines end-systolic and end-diastolic pressure volume relationships

x = Ees;
Ees(2:end) = 0;
A(2:end) = 0;
B(2:end) = 0;
V0(2:end) = 0;

A_atria = Aatria;
B_atria = Batria;

RV_Ees = Ees.*(3/7);

LVparameters = [A       B       Ees      V0]; % parameters recalculated during growth
RVparameters = [A       B.*(3/7)       RV_Ees  V0]; % parameters recalculated during growth
LAparameters = [A_atria B_atria x.*0.12 V0_atria]; % parameters not recalculated
RAparameters = [A_atria B_atria x.*0.12 V0_atria]; % parameters not recalculated

if HLHS == 1
    LVparameters(1,1:4) = 0;
end

%% Timing
                                                                                                                                                                          
HeartRate=HR; %beats/min (From Table 1)
Tmax=Tmax;  % time to end-systole for atria and ventricles (ms)
AV_delay=AVdelay; %delay between atrial and ventricular contraction (ms)

times=[Tmax(:,1) Tmax(:,2) AV_delay HeartRate];

%% Volumes
SBV = SBV; %stressed blood volume (ml)
end