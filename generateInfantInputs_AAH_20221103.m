function inputSet = generateInfantInputs_AAH_20221103(GrowthTime,weight_percentile,HLHS)

% Creates set of model inputs based on input birth weight

% Inputs:
%   GrowthTime: how long to run simulation (default is 1095 days)
%   weight_percentile: which weight trajectory to use (10th, 50th, or 90th)
%   HLHS: flag for HLHS anatomy

% Output:
%   inputSet: (GrowthTime+1x37 array) array of model input parameters used throughout growth 


% Read in weight percentiles
if weight_percentile == 50
    weight_input = readmatrix('inputData/weight50th.csv');
elseif weight_percentile == 90
    weight_input = readmatrix('inputData/weight90th.csv');
elseif weight_percentile == 10
    weight_input = readmatrix('inputData/weight10th.csv');
elseif weight_percentile == 'C'
    weight_input = readmatrix('inputData/weight_GrpC');
elseif weight_percentile == 'D'
    weight_input = readmatrix('inputData/weight_GrpD');
end

weight = weight_input(1:GrowthTime+1,2);
age = weight_input(1:GrowthTime+1,1);


%% Generate Inputs
% original
%factors=[0.645903835   1.164913665   0.998228968   1.578667983   ]; %PVR, Cap, SVR, Cas, factors to adjust value we are scaling from
%CMW factors were fitted to match pulmonary and systemic arterial systolic, diastolic, and mean pressures
factors = [0.338597116300287 0.579046871382574 0.848261399241221 1.81511238548566]; % 11/3/22 recalc

% systemic vascular resistance
SVR_scaled=1.75*factors(3)*(weight./20).^-0.75;  %1.75 and 20 are from C. M. Witzenburg and J. W. Holmes, "Predicting the Time Course of Ventricular Dilation and Thickening Using a Rapid Compartmental Model," J. Cardiovasc. Transl. Res., vol. 11, no. 2, pp. 109–122, Apr. 2018.

% pulmonary vascular resistance
PVR_scaled = 0.117238491*factors(1)*(weight./70).^-0.75; % 0.117238491 and 70 are from W. P. Santamore and D. Burkhoff, “Hemodynamic consequences of ventricular interaction as assessed by model analysis.,” Am. J. Physiol., vol. 260, no. 1 Pt 2, pp. H146-57, Jan. 1991.

% maturation function for PVR and SVR
custom_ends=[33 167]; % PVR SVR, day when custom equation ends in days
PVR_x=zeros(GrowthTime+1,1); SVR_x=zeros(GrowthTime+1,1); 
PVR_x(1:custom_ends(1)+1)= -0.03002528 +    1.02993468./(1 + (age(1:custom_ends(1)+1)./0.9589894).^0.9870873);
SVR_x(1:custom_ends(2)+1) = 0.000367464580858626 + (0.9697715 - -0.03227647)./(1 + (age(1:custom_ends(2)+1)./29.9422944623629).^4.11622150939065);
day0values = [3.74649340066676   -0.870922741716052*weight(1)+7.821136893]; % 11/3/2022 recalc
PVR=(1-PVR_x).*PVR_scaled+PVR_x*day0values(1); %5.3341 is the day zero value of PVR and was fitted to normal data
SVR=(1-SVR_x).*SVR_scaled+SVR_x*day0values(2); %5.097332503 is the day zero value of SVR and was fitted to normal data

% systemic arterial capacitance
Cas_ref = 1.4*factors(4)*(weight(501,1)/70)^1; % reference weight at Day 500
% 1.4 and 70 from W. P. Santamore and D. Burkhoff, “Hemodynamic consequences of ventricular interaction as assessed by model analysis.,” Am. J. Physiol., vol. 260, no. 1 Pt 2, pp. H146-57, Jan. 1991.
SVR_ref = SVR(501); % reference SVR from Day 500
Cas = Cas_ref*(SVR./SVR_ref).^(-4/3); % Cas is dependent on SVR
%Eqn. 2.6 in A. Baretta, C. Corsini, W. Yang, I.E. Vignon-Clementel, A.L. Marsden, J.A. Feinstein, T.Y. Hsia, G. Dubini, F. Migliavacca, G. Pennati, Virtual surgeries in patients with congenital heart disease: A multi-scale modelling test case, Philos. Trans. R. Soc. A Math. Phys. Eng. Sci. 369 (2011) 4316–4330.

% pulmonary arterial capacitance
Cap_ref = 7*factors(2)*(weight(181,1)/70)^1; % Cap uses Day 180 values as reference
%7 and 70 are from D. Kaye et al., “Effects of an interatrial shunt on rest and exercise hemodynamics: results of a computer simulation in heart failure.,” J. Card. Fail., vol. 20, no. 3, pp. 212–21, Mar. 2014.
PVR_ref=PVR(181,1);
Cap = Cap_ref*(PVR./PVR_ref).^(-4/3); % Will eventually need to be adjusted for scaled values past Day 90
%Eqn. 2.6 in A. Baretta, C. Corsini, W. Yang, I.E. Vignon-Clementel, A.L. Marsden, J.A. Feinstein, T.Y. Hsia, G. Dubini, F. Migliavacca, G. Pennati, Virtual surgeries in patients with congenital heart disease: A multi-scale modelling test case, Philos. Trans. R. Soc. A Math. Phys. Eng. Sci. 369 (2011) 4316–4330.
Cap_scaled = Cap_ref*(PVR_scaled./PVR_ref).^(-4/3);

% pulmonary venous resistance
Rvp = 0.005861925*(weight./70).^-0.75;%0.005861925 and 70 are from W. P. Santamore and D. Burkhoff, “Hemodynamic consequences of ventricular interaction as assessed by model analysis.,” Am. J. Physiol., vol. 260, no. 1 Pt 2, pp. H146-57, Jan. 1991.

% systemic characteristic resistance
Rcs = 0.05*(weight./20).^-0.75; %0.05 and 20 are from C. M. Witzenburg and J. W. Holmes, “Predicting the Time Course of Ventricular Dilation and Thickening Using a Rapid Compartmental Model,” J. Cardiovasc. Transl. Res., vol. 11, no. 2, pp. 109–122, Apr. 2018.

% systemic venous resistance
Rvs = 0.005861925*(weight./70).^-0.75; %0.005861925 and 70 are from W. P. Santamore and D. Burkhoff, “Hemodynamic consequences of ventricular interaction as assessed by model analysis.,” Am. J. Physiol., vol. 260, no. 1 Pt 2, pp. H146-57, Jan. 1991.

% pulmonary characteristic resistance
Rcp = 0.0075*(weight./20).^-0.75;%0.0075 and 20 are from C. M. Witzenburg and J. W. Holmes, “Predicting the Time Course of Ventricular Dilation and Thickening Using a Rapid Compartmental Model,” J. Cardiovasc. Transl. Res., vol. 11, no. 2, pp. 109–122, Apr. 2018.

% mitral valve resistance
Rmv = 0.0025*(weight./70).^-0.75;%0.0025 and 70 are from D. Kaye et al., “Effects of an interatrial shunt on rest and exercise hemodynamics: results of a computer simulation in heart failure.,” J. Card. Fail., vol. 20, no. 3, pp. 212–21, Mar. 2014.

% tricuspid valve resistance
Rtv = 0.0025*(weight./70).^-0.75; %0.0025 and 70 are from D. Kaye et al., “Effects of an interatrial shunt on rest and exercise hemodynamics: results of a computer simulation in heart failure.,” J. Card. Fail., vol. 20, no. 3, pp. 212–21, Mar. 2014.

% pulmonary venous capacitance
Cvp = 10.5*(weight./70); %10.5 and 70 are from W. P. Santamore and D. Burkhoff, “Hemodynamic consequences of ventricular interaction as assessed by model analysis.,” Am. J. Physiol., vol. 260, no. 1 Pt 2, pp. H146-57, Jan. 1991.

% systemic venous capacitance
Cvs = 59.5*(weight./70);%59.5 and 70 are from W. P. Santamore and D. Burkhoff, “Hemodynamic consequences of ventricular interaction as assessed by model analysis.,” Am. J. Physiol., vol. 260, no. 1 Pt 2, pp. H146-57, Jan. 1991.

% RatioMPA
% set to 8 such that the mean pressure in the main pulmonary artery was ~1.5 mmHg higher than the distal vessels
ratio_MPAtoPA = 8*(weight./70).^0;

% RatioAo
% set to 16 such that the mean pressure in the ascending aorta was ~4.4 mmHg higher than the distal vessels
ratio_AotoSA = 16*(weight./70).^0;

% RatioUL
% for infants, upper to lower body proportions considered approx. equal
ratio_LBtoUB_Arteries = 1*(weight./70).^0;
ratio_LBtoUB_Veins = 1*(weight./70).^0;

% unloaded volume of atria
V0_atria = 17.5.*(weight/70);

% Aatria: exponent for EDPVR (mL)
% 0.048 and 70 from Hay 2005: "Role of impaired myocardial relaxation in the production of elevated left ventricular filling pressure." American Journal of Physiology-Heart and Circulatory Physiology 288, H1203–H1208
Aatria = 0.048*(weight./70).^-1;

% Batria: scaling factor for EDPVR (mmHg)
% 0.6 and 70 from Hay 2005: "Role of impaired myocardial relaxation in the production of elevated left ventricular filling pressure." American Journal of Physiology-Heart and Circulatory Physiology 288, H1203–H1208
Batria = 0.6*(weight./70).^0;

% ventricular unloaded volume
V0 = 15*(weight./20).^1; %15 and 20 are from C. M. Witzenburg and J. W. Holmes, “Predicting the Time Course of Ventricular Dilation and Thickening Using a Rapid Compartmental Model,” J. Cardiovasc. Transl. Res., vol. 11, no. 2, pp. 109–122, Apr. 2018.

% somatic unloaded radius - from scaled values of V0
r0_s = ((3/(4*pi)).*V0).^(1/3); % LV
r0_s(:,2) = ((3/(4*pi)).*V0).^(1/3); % RV

% A (ventricle) 
A = 0.028571429*(1+0.04)*(weight./70).^-1; %0.028571429 and 70 are from W. P. Santamore and D. Burkhoff, “Hemodynamic consequences of ventricular interaction as assessed by model analysis.,” Am. J. Physiol., vol. 260, no. 1 Pt 2, pp. H146-57, Jan. 1991.

% a, b, and e are ventricular material properties
a = 4/3*pi.*A.*r0_s.^3; %Eqn 12 from C. M. Witzenburg and J. W. Holmes, “Predicting the Time Course of Ventricular Dilation and Thickening Using a Rapid Compartmental Model,” J. Cardiovasc. Transl. Res., vol. 11, no. 2, pp. 109–122, Apr. 2018.

% B (ventricle)
B = 0.175*(weight./70).^0; %0.1750 and 70 are from W. P. Santamore and D. Burkhoff, “Hemodynamic consequences of ventricular interaction as assessed by model analysis.,” Am. J. Physiol., vol. 260, no. 1 Pt 2, pp. H146-57, Jan. 1991.

% b/h0 ratio
b_h0 = B.*r0_s/2; % LV only,  %Eqn 13 from C. M. Witzenburg and J. W. Holmes, “Predicting the Time Course of Ventricular Dilation and Thickening Using a Rapid Compartmental Model,” J. Cardiovasc. Transl. Res., vol. 11, no. 2, pp. 109–122, Apr. 2018.

% end-systolic elastance 
Ees = 20*(weight./20).^-1; %20 and 20kg from C. M. Witzenburg and J. W. Holmes, “Predicting the Time Course of Ventricular Dilation and Thickening Using a Rapid Compartmental Model,” J. Cardiovasc. Transl. Res., vol. 11, no. 2, pp. 109–122, Apr. 2018.

% e/h0 ratio
e_h0 = (2*Ees*pi.*r0_s.^4) / 3; % LV only,  %Eqn 14 from C. M. Witzenburg and J. W. Holmes, “Predicting the Time Course of Ventricular Dilation and Thickening Using a Rapid Compartmental Model,” J. Cardiovasc. Transl. Res., vol. 11, no. 2, pp. 109–122, Apr. 2018.

ratio_btoe = b_h0./e_h0; % valid for LV or RV

% h0 = B*r0 /(2*b), we know B and r0 but do not know b
% therefore, we will set b to a value and modify it to match the loaded
% thickness of the ventricle
b=ones(GrowthTime+1,1)*0.14; %*300; %4; %b=0.06;
e=1./ratio_btoe.*b;

% unloaded thickness
h0_s = ((B/2).*r0_s(:,1))./b; % LV, decimal value depends on b,  %Eqn 13 from C. M. Witzenburg and J. W. Holmes, “Predicting the Time Course of Ventricular Dilation and Thickening Using a Rapid Compartmental Model,” J. Cardiovasc. Transl. Res., vol. 11, no. 2, pp. 109–122, Apr. 2018.
h0_s(:,2) = (((B*(3/7)/2).*r0_s(:,2))./b); % RV

% heart rate entered manually
HR = readmatrix('inputData/HR.csv'); % from 50th percentile in Fleming et al., "Normal ranges of heart rate and respiratory rate in children from birth to 18 years of age: a systematic review of observational studies.," Lancet., vol. 377, no. 9770, pp. 1011-1018, Mar. 2011.
weight50=readmatrix('inputData/weight50th.csv');
weight50 = weight50(:,2);
HR = HR(1:GrowthTime+1,2);
% adjusting for different weight percentiles
HRnew=HR.*(weight./weight50(:,1)).^(-1/4);
HR=HRnew;

% AV delay
AVdelay = (160*62)./HR; %160 ms and 62 beats/min are from D. Kaye et al., “Effects of an interatrial shunt on rest and exercise hemodynamics: results of a computer simulation in heart failure.,” J. Card. Fail., vol. 20, no. 3, pp. 212–21, Mar. 2014.

% Tmax(1)% - used to calculate Tmax(1)
percentval(1)=125/ (1/62*60*1000);%atria 125 ms and 62 beats/min are from D. Kaye et al., “Effects of an interatrial shunt on rest and exercise hemodynamics: results of a computer simulation in heart failure.,” J. Card. Fail., vol. 20, no. 3, pp. 212–21, Mar. 2014.
percentval(2)=200/ (1/62*60*1000);%ventricle 125 ms and 62 beats/min are from D. Kaye et al., “Effects of an interatrial shunt on rest and exercise hemodynamics: results of a computer simulation in heart failure.,” J. Card. Fail., vol. 20, no. 3, pp. 212–21, Mar. 2014.
%60 converts from min to s and 1000 converts from ms to s

Tmax1_percent = percentval(1)*(weight./70).^-0.0687814; %-0.0687814 fitted to rat and dog data
%C. Witzenburg, J.W. Holmes, The Impact of Hemodynamic Reflex Compensation Following Myocardial Infarction on Subsequent Ventricular Remodeling., J. Biomech. Eng. 
%J.W. Holmes, Candidate mechanical stimuli for hypertrophy during volume overload., J. Appl. Physiol. 97 (2004) 1453–1460.

% Tmax(1)
Tmax = (Tmax1_percent.*60*1000)./HR; %60 converts from min to s and 1000 converts from ms to s

% Tmax(2)% - used to calculate Tmax(2)
Tmax2_percent =percentval(2)*(weight./70).^-0.0687814; %-0.0687814 fitted to rat and dog data
%C. Witzenburg, J.W. Holmes, The Impact of Hemodynamic Reflex Compensation Following Myocardial Infarction on Subsequent Ventricular Remodeling., J. Biomech. Eng. 
%J.W. Holmes, Candidate mechanical stimuli for hypertrophy during volume overload., J. Appl. Physiol. 97 (2004) 1453–1460.

%Tmax(2)
Tmax(:,2) = (Tmax2_percent.*60*1000)./HR; %60 converts from min to s and 1000 converts from ms to s

% stressed blood volume
SBV_s = 875*(weight./70); %875 and 70 are from W. P. Santamore and D. Burkhoff, “Hemodynamic consequences of ventricular interaction as assessed by model analysis.,” Am. J. Physiol., vol. 260, no. 1 Pt 2, pp. H146-57, Jan. 1991.

inputSet = [A, Aatria, AVdelay, B, Batria, Cap, Cap_scaled, Cas, Cvp, Cvs, Ees, HR, PVR, ratio_AotoSA, ...
            ratio_LBtoUB_Arteries, ratio_LBtoUB_Veins, ratio_MPAtoPA, Rcp, Rcs, Rmv, Rtv, Rvp, Rvs, ...
            SBV_s, SVR, Tmax, V0, PVR_scaled, r0_s, h0_s, a(:,1), b, e(:,1), V0_atria];
end
