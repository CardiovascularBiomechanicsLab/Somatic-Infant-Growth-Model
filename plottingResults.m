clear; close all; clc;

% Generates Figures 5-8 and 10 from Hiebing et al 2023 ("A Computational Model of Ventricular Dimensions and Hemodynamics in Growing Infants")

load('sampleOutput/coarcOutput.mat')

age = 0:1094;
ageTicksLong = 0:6:36; % months
fontSize = 16;

% Plot colors
modelColor = '#4e67c8'; % main model color, dark blue
modelColor_coarc = '#f14124'; % coarc, red

% Figure 5: 
%   (A) LV unloaded radius for somatic and CoA
%   (B) RV unloaded radius for somatic and CoA
%   (C) LV unloaded thickness for somatic and CoA
%   (D) RV unloaded thickness for somatic and CoA
figure(5)
set(gcf,'color','w');
    subplot(2,2,1)
        plot(age,r0_s(1:end-1,1),'LineWidth',2,'Color',modelColor); hold on;
        plot(age,r0(1:end-1,1),'LineWidth',2,'Color',modelColor_coarc);
        title('Left Ventricle');
        xlabel('Age (Months)');
        xticks([0 6*30.4167 12*30.4167 18*30.4167 24*30.4167 30*30.4167 36*30.4167]);
        xticklabels(ageTicksLong);
        xlim([0 1095]); ylim([0.8 1.45]);
        ylabel('Unloaded radius (cm)');
        legend({'r_0_s (somatic)','r_0 (coarctation)'},'Location','southeast')
        set(gca,'FontSize',fontSize)
        legend boxoff
    subplot(2,2,2)
        plot(age,r0_s(1:end-1,2),'LineWidth',2,'Color',modelColor); hold on;
        plot(age,r0(1:end-1,2),'LineWidth',2,'Color',modelColor_coarc);
        title('Right Ventricle');
        xlabel('Age (Months)');
        xticks([0 6*30.4167 12*30.4167 18*30.4167 24*30.4167 30*30.4167 36*30.4167]);
        xticklabels(ageTicksLong);
        xlim([0 1095]); ylim([0.8 1.45]);
        ylabel('Unloaded radius (cm)');
        legend({'r_0_s (somatic)','r_0 (coarctation)'},'Location','southeast')
        set(gca,'FontSize',fontSize)
        legend boxoff
    subplot(2,2,3)
        plot(age,h0_s(1:end-1,1),'LineWidth',2,'Color',modelColor); hold on;
        plot(age,h0(1:end-1,1),'LineWidth',2,'Color',modelColor_coarc);
        xlabel('Age (Months)');
        xticks([0 6*30.4167 12*30.4167 18*30.4167 24*30.4167 30*30.4167 36*30.4167]);
        xticklabels(ageTicksLong);
        xlim([0 1095]); ylim([0.2 1]);
        ylabel('Unloaded thickness (cm)');
        legend({'h_0_s (somatic)','h_0 (coarctation)'},'Location','southeast')
        set(gca,'FontSize',fontSize)
        legend boxoff
    subplot(2,2,4)
        plot(age,h0_s(1:end-1,2),'LineWidth',2,'Color',modelColor); hold on;
        plot(age,h0(1:end-1,2),'LineWidth',2,'Color',modelColor_coarc);
        xlabel('Age (Months)');
        xticks([0 6*30.4167 12*30.4167 18*30.4167 24*30.4167 30*30.4167 36*30.4167]);
        xticklabels(ageTicksLong);
        xlim([0 1095]); ylim([0.2 1]);
        ylabel('Unloaded thickness (cm)');
        legend({'h_0_s (somatic)','h_0 (coarctation)'},'Location','northeast')
        set(gca,'FontSize',fontSize)
        legend boxoff

%% Figure 6:
%   (A) Mean arterial pressure (MAP)
%   (B) Mean pulmonary arterial pressure (MPAP)

load('sampleOutput/somaticOutput.mat')

% Park and Menard 1989: "Normative Oscillometric Blood Pressure Values in the First 5 Years in an Office Setting"
PM_Time = [2, 2.5*7, 3*30.4167, 8.5*30.4167, 12*30.4167, 24*30.4167, 35.5*30.4167]; % 2 days, 2.5 weeks, 3 months, 8.5 months, 12 months, 24 months, 36 months
PM_MBP_50 = [50 58 72 70 71 70 71];
PM_MBP_10 = [40 48 60 61 62 62 63];
PM_MBP_90 = [60 71 86 82 81 79 80];
PM_yneg_MBP = PM_MBP_50-PM_MBP_10;
PM_ypos_MBP = PM_MBP_90-PM_MBP_50;

% Kent et al 2007: “Normative blood pressure data in the early neonatal period”
% and Kent et al 2007: "Blood pressure in the first year of life in healthy infants born at term"
Ke_Time = [0, 1, 26, 52];
Ke_Time = Ke_Time.*7; % convert to days
Ke_MBP_50 = [48, 51, 75, 75];
Ke_MBP_95 = [57, 61, 88, 89];
Ke_MBP_5 = [39, 42, 60, 61];
Ke_yneg_MBP = Ke_MBP_50-Ke_MBP_5;
Ke_ypos_MBP = Ke_MBP_95-Ke_MBP_50;

% Kang et al 2016: "“Dynamic Changes of Pulmonary Arterial Pressure and Ductus Arteriosus in Human Newborns From Birth to 72 Hours of Age”
Ka_Time = [12, 24, 48, 72]; % hrs
% convert to days
Ka_Time = Ka_Time./24;
Ka_PAMP = [40.94, 34.39, 26.23, 25.25];
Ka_PAMP_SD = [9.32, 9.89, 7.49, 8.29];

% Qi et al 2014: “Anatomical and hemodynamic evaluations of the heart and pulmonary arterial pressure in healthy children residing at high altitude in China”
% Subjects sorted into age bins; assuming halfway point
Qi_Time = [0.5, 3.5, 9.5, 24.5]; % in months
% convert to days
Qi_Time = Qi_Time.*30.4368;
Qi_MPAP_SL = [14.5, 15.4, 14.8, 17.7];
Qi_MPAP_SL_SD = [8.6, 7.5, 5.2, 6.9];

% Emmanouilides et al 1964 "Pulmonary arterial pressure changes in human newborn infants from birth to 3 days of age"
Em_Time = [0, 1, 54/24]; % days
Em_PAMP = [40, 31, 20];
Em_PAMP_SD = [8, 9, 9];

sys_S = zeros(1,numel(VolumeArray));
sys_D = zeros(1,numel(VolumeArray));

for i = 1:numel(VolumeArray)
        % Calculate systolic, diastolic, and mean BP
        sys_S(i) = max(PressureArray{i}(:,6));
        sys_D(i) = min(PressureArray{i}(:,6));
end

sys_M = (sys_S + 2.*sys_D)./3;

pulm_S = zeros(1,numel(VolumeArray));
pulm_D = zeros(1,numel(VolumeArray));

for i = 1:numel(VolumeArray)
    pulm_S(i) = max(PressureArray{i}(:,11));
    pulm_D(i) = min(PressureArray{i}(:,11));
end

pulm_M = (pulm_S + 2.*pulm_D)./3;

figure(6)
set(gcf,'color','w');
MAP_subp = subplot(2,1,1); % MAP
    errorbar(PM_Time, PM_MBP_50,PM_ypos_MBP,PM_yneg_MBP, '.', 'color', [0.5 0.5 0.5], 'MarkerSize', 30); hold on;
    errorbar(Ke_Time, Ke_MBP_50,Ke_ypos_MBP,Ke_yneg_MBP, '.', 'color', [0.25 0.25 0.25], 'MarkerSize', 30);
    plot(age,sys_M,'LineWidth',2,'Color',modelColor);
    xlabel('Age (months)');
    xlim([0 1095]);
    xticks([0 6*30.4167 12*30.4167 18*30.4167 24*30.4167 30*30.4167 36*30.4167]);
    xticklabels(ageTicksLong);
    ylabel('MAP (mmHg)');
    legend({'Park 1989','Kent 2007','Model'});
    set(gca,'FontSize',fontSize)
    set(MAP_subp,'Position',[0.1300    0.5901    0.7750    0.3349]);
    legend boxoff
MPAP_subp = subplot(2,1,2); % MPAP
    errorbar(Ka_Time,Ka_PAMP,Ka_PAMP_SD,'.','Color',[0.5,0.5,0.5],'MarkerSize',30); hold on;
    errorbar(Qi_Time,Qi_MPAP_SL,Qi_MPAP_SL_SD,'.','Color',[0.75,0.75,0.75],'MarkerSize',30);
    errorbar(Em_Time, Em_PAMP,Em_PAMP_SD,'.','Color',[0.25,0.25,0.25],'MarkerSize',30);
    plot(age,pulm_M,'LineWidth',2,'Color',modelColor);
    xlabel('Age (months)');
    xlim([0 1095]);
    xticks([0 6*30.4167 12*30.4167 18*30.4167 24*30.4167 30*30.4167 36*30.4167]);
    xticklabels(ageTicksLong);
    ylabel('MPAP (mmHg)');
    legend({'Kang 2016','Qi 2014','Emmanouilides 1964','Model'});
    set(gca,'FontSize',fontSize)
    set(MPAP_subp,'Position',[ 0.1300    0.1162    0.7750    0.3349]);
    legend boxoff
%% Figure 7
%   (A) LV ED volume
%   (B) LV ES volume
%   (C) RV ED volume
%   (D) RV ES volume

% BSA (m^2) to age (months) conversion factors B1, B2, B3, B4
B = [396.6342996	-247.7948178	78.1	-9.339726027];

% Akiba et al 1986: “Echocardiographic Measurements of Left Ventricle in Normal Infants and Children”
% Note: Equation given: LVEDV = 74.1*BSA^1.28. Input BSA range from birth (0.2) to ~3 years (0.7)
% No standard deviations given
Ak_BSA = 0.2:0.01:0.7;
Ak_LVEDV = 74.1.*Ak_BSA.^1.28;
Ak_Age = (B(1).*Ak_BSA.^3 + B(2).*Ak_BSA.^2 + B(3).*Ak_BSA+B(4)).*30.4167; % convert to days

% Colan 2016: “Normal Echocardiographic Values for Cardiovascular Structures”
% LVEDV
    Co_LVEDV_BSA = readmatrix('plottingData/colan_LVEDV_all.csv');
    Co_BSA = Co_LVEDV_BSA(:,1); % m^2
    Co_Age_LVEDV = (B(1).*Co_BSA.^3 + B(2).*Co_BSA.^2 + B(3).*Co_BSA+B(4)).*30.4167; % convert to days
    Co_LVEDV_50 = Co_LVEDV_BSA(:,2); % mL
    Co_LVEDV_2SDpos = Co_LVEDV_BSA(:,3); % mL
    Co_LVEDV_2SDneg = Co_LVEDV_BSA(:,4); % mL
% LVESV
    Co_LVESV_BSA_50 = readmatrix('plottingData/colan_LVESV_50th.csv');
    Co_LVESV_BSA_2SDneg = readmatrix('plottingData/colan_LVESV_2SDneg.csv');
    Co_LVESV_BSA_2SDpos = readmatrix('plottingData/colan_LVESV_2SDpos.csv');

    Co_BSA_50 = Co_LVESV_BSA_50(:,1);
    Co_Age_LVESV = (B(1).*Co_BSA_50.^3 + B(2).*Co_BSA_50.^2 + B(3).*Co_BSA_50+B(4)).*30.4167; % convert to days
    Co_LVESV_50 = Co_LVESV_BSA_50(:,2);
    Co_BSA_2SDneg = Co_LVESV_BSA_2SDneg(:,1);
    Co_Age_2SDneg_LVESV = (B(1).*Co_BSA_2SDneg.^3 + B(2).*Co_BSA_2SDneg.^2 + B(3).*Co_BSA_2SDneg+B(4)).*30.4167; % convert to days
    Co_LVESV_2SDneg = Co_LVESV_BSA_2SDneg(:,2);
    Co_BSA_2SDpos = Co_LVESV_BSA_2SDpos(:,1);
    Co_Age_2SDpos_LVESV = (B(1).*Co_BSA_2SDpos.^3 + B(2).*Co_BSA_2SDpos.^2 + B(3).*Co_BSA_2SDpos+B(4)).*30.4167; % convert to days
    Co_LVESV_2SDpos = Co_LVESV_BSA_2SDpos(:,2);


% Olivieri 2020: “Normal right and left ventricular volumes prospectively obtained from 
% cardiovascular magnetic resonance in awake, healthy, 0-12 year old children”
% LVEDV
    Ol_LVEDV_BSA_50 = readmatrix('plottingData/olivieri_LVEDV_50th.csv');
    Ol_LVEDV_BSA_2SDneg = readmatrix('plottingData/olivieri_LVEDV_2SDneg.csv');
    Ol_LVEDV_BSA_2SDpos = readmatrix('plottingData/olivieri_LVEDV_2SDpos.csv');

    Ol_BSA_50_LVEDV = Ol_LVEDV_BSA_50(:,1);
    Ol_Age_50_LVEDV = (B(1).*Ol_BSA_50_LVEDV.^3 + B(2).*Ol_BSA_50_LVEDV.^2 + B(3).*Ol_BSA_50_LVEDV+B(4)).*30.4167; % convert to days
    Ol_LVEDV_50 = Ol_LVEDV_BSA_50(:,2);
    Ol_BSA_2SDneg_LVEDV = Ol_LVEDV_BSA_2SDneg(:,1);
    Ol_Age_2SDneg_LVEDV = (B(1).*Ol_BSA_2SDneg_LVEDV.^3 + B(2).*Ol_BSA_2SDneg_LVEDV.^2 + B(3).*Ol_BSA_2SDneg_LVEDV+B(4)).*30.4167; % convert to days
    Ol_LVEDV_2SDneg = Ol_LVEDV_BSA_2SDneg(:,2);
    Ol_BSA_2SDpos_LVEDV = Ol_LVEDV_BSA_2SDpos(:,1);
    Ol_Age_2SDpos_LVEDV = (B(1).*Ol_BSA_2SDpos_LVEDV.^3 + B(2).*Ol_BSA_2SDpos_LVEDV.^2 + B(3).*Ol_BSA_2SDpos_LVEDV+B(4)).*30.4167; % convert to days
    Ol_LVEDV_2SDpos = Ol_LVEDV_BSA_2SDpos(:,2);
% LVESV
    Ol_LVESV_BSA_50 = readmatrix('plottingData/olivieri_LVESV_50th.csv');
    Ol_LVESV_BSA_2SDneg = readmatrix('plottingData/olivieri_LVESV_2SDneg.csv');
    Ol_LVESV_BSA_2SDpos = readmatrix('plottingData/olivieri_LVESV_2SDpos.csv');

    Ol_BSA_50_LVESV = Ol_LVESV_BSA_50(:,1);
    Ol_Age_50_LVESV = (B(1).*Ol_BSA_50_LVESV.^3 + B(2).*Ol_BSA_50_LVESV.^2 + B(3).*Ol_BSA_50_LVESV+B(4)).*30.4167; % convert to days
    Ol_LVESV_50 = Ol_LVESV_BSA_50(:,2);
    Ol_BSA_2SDneg_LVESV = Ol_LVESV_BSA_2SDneg(:,1);
    Ol_Age_2SDneg_LVESV = (B(1).*Ol_BSA_2SDneg_LVESV.^3 + B(2).*Ol_BSA_2SDneg_LVESV.^2 + B(3).*Ol_BSA_2SDneg_LVESV+B(4)).*30.4167; % convert to days
    Ol_LVESV_2SDneg = Ol_LVESV_BSA_2SDneg(:,2);
    Ol_BSA_2SDpos_LVESV = Ol_LVESV_BSA_2SDpos(:,1);
    Ol_Age_2SDpos_LVESV = (B(1).*Ol_BSA_2SDpos_LVESV.^3 + B(2).*Ol_BSA_2SDpos_LVESV.^2 + B(3).*Ol_BSA_2SDpos_LVESV+B(4)).*30.4167; % convert to days
    Ol_LVESV_2SDpos = Ol_LVESV_BSA_2SDpos(:,2);
% RVEDV
    Ol_RVEDV_BSA_50 = readmatrix('plottingData/olivieri_RVEDV_50th.csv');
    Ol_RVEDV_BSA_2SDneg = readmatrix('plottingData/olivieri_RVEDV_2SDneg.csv');
    Ol_RVEDV_BSA_2SDpos = readmatrix('plottingData/olivieri_RVEDV_2SDpos.csv');

    Ol_BSA_50_RVEDV = Ol_RVEDV_BSA_50(:,1);
    Ol_Age_50_RVEDV = (B(1).*Ol_BSA_50_RVEDV.^3 + B(2).*Ol_BSA_50_RVEDV.^2 + B(3).*Ol_BSA_50_RVEDV+B(4)).*30.4167; % convert to days
    Ol_RVEDV_50 = Ol_RVEDV_BSA_50(:,2);
    Ol_BSA_2SDneg_RVEDV = Ol_RVEDV_BSA_2SDneg(:,1);
    Ol_Age_2SDneg_RVEDV = (B(1).*Ol_BSA_2SDneg_RVEDV.^3 + B(2).*Ol_BSA_2SDneg_RVEDV.^2 + B(3).*Ol_BSA_2SDneg_RVEDV+B(4)).*30.4167; % convert to days
    Ol_RVEDV_2SDneg = Ol_RVEDV_BSA_2SDneg(:,2);
    Ol_BSA_2SDpos_RVEDV = Ol_RVEDV_BSA_2SDpos(:,1);
    Ol_Age_2SDpos_RVEDV = (B(1).*Ol_BSA_2SDpos_RVEDV.^3 + B(2).*Ol_BSA_2SDpos_RVEDV.^2 + B(3).*Ol_BSA_2SDpos_RVEDV+B(4)).*30.4167; % convert to days
    Ol_RVEDV_2SDpos = Ol_RVEDV_BSA_2SDpos(:,2);


% Lytrivi 2011: “Normal Values for Left Ventricular Volume in Infants and Young Children 
% by the Echocardiographic Subxiphoid Five-Sixth Area by Length (Bullet) Method”
% Note: Equation given: LVEDV = 1/8 (625 BSA^(69/50) - 19). Input BSA range from birth (0.2) to ~3 years (0.6)
% No standard deviations given
Ly_BSA = 0.2:0.01:0.7;
Ly_LVEDV = (1/8)*(625.*Ly_BSA.^(69/50) - 19);
Ly_Age = (B(1).*Ly_BSA.^3 + B(2).*Ly_BSA.^2 + B(3).*Ly_BSA+B(4)).*30.4167; % convert to days

% Graham et al 1972: "Right Ventricular Volume Determinations in Children: 
% Normal Values and Observations with Volume or Pressure Overload”
Gr_RVEDV_BSA = readmatrix('plottingData/graham_RVEDV_all.csv');
Gr_BSA = Gr_RVEDV_BSA(:,1); % m^2
Gr_Age = (B(1).*Gr_BSA.^3 + B(2).*Gr_BSA.^2 + B(3).*Gr_BSA+B(4)).*30.4167; % convert to days
Gr_RVEDV_50 = Gr_RVEDV_BSA(:,2); % mL
Gr_RVEDV_10 = Gr_RVEDV_BSA(:,3); % mL
Gr_RVEDV_90 = Gr_RVEDV_BSA(:,4); % mL

% Thilenius, Arcilla 1974: "Angiographic Right and Left Ventricular Volume
% Determination in Normal Infants and Children"
% Gives RVESV as RVESV = 35.102+0.855*age-0.271*height+1.274*weight-0.125*HR
% age in years, height in cm, weight in kg
% SD +/- 1.28 mL
Th_RVESV_all = readmatrix('plottingData/thilenius_RVESV.csv'); 
Th_Age = Th_RVESV_all(:,1);
Th_RVESV_50 = Th_RVESV_all(:,2);
Th_RVESV_pos = Th_RVESV_all(:,3); % +1 SD
Th_RVESV_neg = Th_RVESV_all(:,4); % -1 SD

% Buechel et al 2009: "Normal right- and left ventricular volumes and myocardial mass in
% children measured by steady state free precession cardiovascular
% magnetic resonance"
% Averaged 50% and 95% CI for males and females
Bu_RVESV_all = readmatrix('plottingData/buechel_RVESV.csv');
Bu_Age = Bu_RVESV_all(:,1);
Bu_RVESV_50 = Bu_RVESV_all(:,2);
Bu_RVESV_95 = Bu_RVESV_all(:,3);
Bu_RVESV_5 = Bu_RVESV_all(:,4);

% Lange et al 1982: "Size and Function of the Human Left and Right Ventricles During Growth
% Normative Angiographic Data"
% Gives following eqn for RVESV: 26.3*BSA^1.18
La_BSA = 0.2:0.01:0.7;
La_RVESV_50 = 26.3.*La_BSA.^1.18;
La_Age = (B(1).*La_BSA.^3 + B(2).*La_BSA.^2 + B(3).*La_BSA+B(4)).*30.4167; % convert to days

figure(7)
set(gcf,'color','w');
LVEDV_subp = subplot(2,2,1); % LVEDV
    Ak = plot(Ak_Age,Ak_LVEDV,'--','LineWidth',2,'Color',[0.5,0.5,0.5]); hold on
    Co_LVEDV = plot(Co_Age_LVEDV,Co_LVEDV_50,'--','LineWidth',2,'Color',[0.75,0.75,0.75]);
    Co_SD_LVEDV = plot(Co_Age_LVEDV,Co_LVEDV_2SDpos,':','LineWidth',2,'Color',[0.75,0.75,0.75]);
    plot(Co_Age_LVEDV,Co_LVEDV_2SDneg,':','LineWidth',2,'Color',[0.75,0.75,0.75]);
    Ol_LVEDV = plot(Ol_Age_50_LVEDV,Ol_LVEDV_50,'--','LineWidth',2,'Color',[0.25,0.25,0.25]);
    Ol_SD_LVEDV = plot(Ol_Age_2SDpos_LVEDV,Ol_LVEDV_2SDpos,':','LineWidth',2,'Color',[0.25,0.25,0.25]);
    plot(Ol_Age_2SDneg_LVEDV,Ol_LVEDV_2SDneg,':','LineWidth',2,'Color',[0.25,0.25,0.25]);
    Ly = plot(Ly_Age,Ly_LVEDV,'--','LineWidth',2,'Color',[0.65,0.65,0.65]);
    model_LVEDV = plot(age,ValuesofInterest(5,:),'LineWidth',2,'Color',modelColor);
    xlim([0 1095]);
    ylim([0 60]);
    %title('Left Ventricular End Diastolic Volume');
    xlabel('Age (Months)');
    xticks([0 6*30.4167 12*30.4167 18*30.4167 24*30.4167 30*30.4167 36*30.4167]);
    xticklabels(ageTicksLong);
    ylabel('LV ED Volume (mL)');
    set(gca,'FontSize',fontSize)
    legend_entries_LVEDV = [Ak, Co_LVEDV, Co_SD_LVEDV, Ol_LVEDV, Ol_SD_LVEDV, Ly, model_LVEDV];
    legend(legend_entries_LVEDV,{'Akiba 1986','Colan 2016','Colan 2016 +/- 2SD,','Olivieri 2020','Olivieri 2020 +/- 2SD','Lytrivi 2011','Model'},'Location','northwest','FontSize',12); 
    volFigSizeTop = get(LVEDV_subp, 'Position');
    legend boxoff
subplot(2,2,2) % LVESV
    Co_LVESV = plot(Co_Age_LVESV,Co_LVESV_50,'--','LineWidth',2,'Color',[0.75,0.75,0.75]); hold on;
    Co_SD_LVESV = plot(Co_Age_2SDneg_LVESV, Co_LVESV_2SDneg,':','LineWidth',2,'Color',[0.75,0.75,0.75]);
    plot(Co_Age_2SDpos_LVESV, Co_LVESV_2SDpos,':','LineWidth',2,'Color',[0.75,0.75,0.75]);
    Ol_LVESV = plot(Ol_Age_50_LVESV,Ol_LVESV_50,'--','LineWidth',2,'Color',[0.25,0.25,0.25]);
    Ol_SD_LVESV = plot(Ol_Age_2SDpos_LVESV,Ol_LVESV_2SDpos,':','LineWidth',2,'Color',[0.25,0.25,0.25]);
    plot(Ol_Age_2SDneg_LVESV,Ol_LVESV_2SDneg,':','LineWidth',2,'Color',[0.25,0.25,0.25]);
    model_LVESV = plot(age,ValuesofInterest(6,:),'LineWidth',2,'Color',modelColor);
    xlim([0 1095]);
    ylim([0 60]);
    % title('Left Ventricular End Systolic Volume');
    xlabel('Age (Months)');
    xticks([0 6*30.4167 12*30.4167 18*30.4167 24*30.4167 30*30.4167 36*30.4167]);
    xticklabels(ageTicksLong);
    ylabel('LV ES Volume (mL)');
    set(gca,'FontSize',fontSize)
    legend_entries_LVESV = [Co_LVESV, Co_SD_LVESV, Ol_LVESV, Ol_SD_LVESV, model_LVESV];
    legend(legend_entries_LVESV,{'Colan 2016','Colan 2016 +/- 2SD,','Olivieri 2020','Olivieri 2020 +/- 2SD','Model'},'Location','northwest','FontSize',12);
    legend boxoff
    
RVEDV_subp = subplot(2,2,3); % RVEDV
    Gr_RVEDV = plot(Gr_Age,Gr_RVEDV_50,'--','LineWidth',2,'Color',[0.75,0.75,0.75]); hold on;
    Gr_SD = plot(Gr_Age,Gr_RVEDV_10,':','LineWidth',2,'Color',[0.75,0.75,0.75]);
    plot(Gr_Age,Gr_RVEDV_90,':','LineWidth',2,'Color',[0.75,0.75,0.75]);
    Ol_RVEDV = plot(Ol_Age_50_RVEDV,Ol_RVEDV_50,'--','LineWidth',2,'Color',[0.25,0.25,0.25]);
    Ol_SD_RVEDV = plot(Ol_Age_2SDpos_RVEDV,Ol_RVEDV_2SDpos,':','LineWidth',2,'Color',[0.25,0.25,0.25]);
    plot(Ol_Age_2SDneg_RVEDV,Ol_RVEDV_2SDneg,':','LineWidth',2,'Color',[0.25,0.25,0.25]);
    model_RVEDV = plot(age,ValuesofInterest_R(5,:),'LineWidth',2,'Color',modelColor);
    xlim([0 1095]);
    ylim([0 60]);
    % title('Right Ventricular End Diastolic Volume');
    xlabel('Age (Months)');
    xticks([0 6*30.4167 12*30.4167 18*30.4167 24*30.4167 30*30.4167 36*30.4167]);
    xticklabels(ageTicksLong);
    ylabel('RV ED Volume (mL)');
    set(gca,'FontSize',fontSize)
    legend_entries_RVEDV = [Gr_RVEDV, Gr_SD, Ol_RVEDV, Ol_SD_RVEDV, model_RVEDV];
    legend(legend_entries_RVEDV,{'Graham 1972','Graham 1972 (P10/90)','Olivieri 2020','Olivieri 2020 +/- 2SD','Model'},'Location','northwest','FontSize',12);   
    volFigSizeBottom = get(RVEDV_subp, 'Position');
    legend boxoff
subplot(2,2,4) % RVESV
    Th_RVESV = plot(Th_Age,Th_RVESV_50,'--','LineWidth',2,'Color',[0.75,0.75,0.75]); hold on;
    Th_SD = plot(Th_Age,Th_RVESV_pos,':','LineWidth',2,'Color',[0.75,0.75,0.75]);
    plot(Th_Age,Th_RVESV_neg,':','LineWidth',2,'Color',[0.75,0.75,0.75]);
    Bu_RVESV = plot(Bu_Age,Bu_RVESV_50,'--','LineWidth',2,'Color',[0.25,0.25,0.25]); hold on;
    Bu_SD = plot(Bu_Age,Bu_RVESV_95,':','LineWidth',2,'Color',[0.25,0.25,0.25]);
    plot(Bu_Age,Bu_RVESV_5,':','LineWidth',2,'Color',[0.25,0.25,0.25]);
    La_RVESV = plot(La_Age,La_RVESV_50,'--','LineWidth',2,'Color',[0.5,0.5,0.5]);
    model_RVESV = plot(age,ValuesofInterest_R(6,:),'LineWidth',2,'Color',modelColor);
    xlim([0 1095]);
    ylim([0 60]);
    % title('Right Ventricular End Systolic Volume');
    xlabel('Age (Months)');
    xticks([0 6*30.4167 12*30.4167 18*30.4167 24*30.4167 30*30.4167 36*30.4167]);
    xticklabels(ageTicksLong);
    ylabel('RV ES Volume (mL)');
    set(gca,'FontSize',fontSize)
    legend_entries_RVESV = [Th_RVESV, Th_SD, Bu_RVESV, Bu_SD, La_RVESV, model_RVESV];
    legend(legend_entries_RVESV,{'Thilenius 1974','Thilenius 1974 (+/- 1SD)','Buechel 2009','Buechel 2009 (95% CI)','Lange 1982','Model'},'Location','northwest','FontSize',12);
    legend boxoff
%% Figure 8
%   (A) LV ED thickness
%   (B) LV ES thickness
%   (C) RV ED thickness
%   (D) RV ES thickness

% BSA (m^2) to age (months) conversion factors B1, B2, B3, B4
B = [396.6342996	-247.7948178	78.1	-9.339726027];

% Kampmann et al 2000: "Normal values of M mode echocardiographic
% measurements of more than 2000 healthy infants
% and children in central Europe"
% LVPWED
    Ka_LVPWED_all = readmatrix('plottingData/kampmann_LVPWED.csv');
    Ka_LVPWED_BSA_50 = Ka_LVPWED_all(:,1);
    Ka_LVPWED_Age_50 = (B(1).*Ka_LVPWED_BSA_50.^3 + B(2).*Ka_LVPWED_BSA_50.^2 + B(3).*Ka_LVPWED_BSA_50+B(4)).*30.4167; % convert to days
    Ka_LVPWED_50 = Ka_LVPWED_all(:,2);
    Ka_LVPWED_BSA_90 = Ka_LVPWED_all(:,3);
    Ka_LVPWED_Age_90 = (B(1).*Ka_LVPWED_BSA_90.^3 + B(2).*Ka_LVPWED_BSA_90.^2 + B(3).*Ka_LVPWED_BSA_90+B(4)).*30.4167; % convert to days
    Ka_LVPWED_90 = Ka_LVPWED_all(:,4);
    Ka_LVPWED_BSA_10 = Ka_LVPWED_all(:,5);
    Ka_LVPWED_Age_10 = (B(1).*Ka_LVPWED_BSA_10.^3 + B(2).*Ka_LVPWED_BSA_10.^2 + B(3).*Ka_LVPWED_BSA_10+B(4)).*30.4167; % convert to days
    Ka_LVPWED_10 = Ka_LVPWED_all(:,6);
% LVPWES
     Ka_LVPWES_all = readmatrix('plottingData/kampmann_LVPWES.csv');
     Ka_LVPWES_BSA = Ka_LVPWES_all(:,1);
     Ka_LVPWES_Age = (B(1).*Ka_LVPWES_BSA.^3 + B(2).*Ka_LVPWES_BSA.^2 + B(3).*Ka_LVPWES_BSA+B(4)).*30.4167; % convert to days
     Ka_LVPWES_50 = Ka_LVPWES_all(:,2);
     Ka_LVPWES_90 = Ka_LVPWES_all(:,3);
     Ka_LVPWES_10 = Ka_LVPWES_all(:,4);
% RVAWED
     Ka_RVAWED_all = readmatrix('plottingData/kampmann_RVAWED.csv');
     Ka_RVAWED_BSA = Ka_RVAWED_all(:,1);
     Ka_RVAWED_Age = (B(1).*Ka_RVAWED_BSA.^3 + B(2).*Ka_RVAWED_BSA.^2 + B(3).*Ka_RVAWED_BSA+B(4)).*30.4167; % convert to days
     Ka_RVAWED_50 = Ka_RVAWED_all(:,2);
     Ka_RVAWED_90 = Ka_RVAWED_all(:,3);
     Ka_RVAWED_10 = Ka_RVAWED_all(:,4);

% Akiba et al 1986: "Echocardiographic Measurements of Left Ventricle in Normal Infants and
% Children" 
% LVPWED
    % Gives equation: LVPWED = 4.4*BSA^0.45
    Ak_LVPWED_BSA = 0.2:0.01:0.7;
    Ak_LVPWED_50 = 4.4.*Ak_LVPWED_BSA.^0.45;
    Ak_LVPWED_Age = (B(1).*Ak_LVPWED_BSA.^3 + B(2).*Ak_LVPWED_BSA.^2 + B(3).*Ak_LVPWED_BSA+B(4)).*30.4167; % convert to days
% LVPWES
    % Gives equation: LVPWES = 9.2*BSA^0.44
    Ak_LVPWES_BSA = 0.2:0.01:0.7;
    Ak_LVPWES_50 = 9.2.*Ak_LVPWES_BSA.^0.44;
    Ak_LVPWES_Age = (B(1).*Ak_LVPWES_BSA.^3 + B(2).*Ak_LVPWES_BSA.^2 + B(3).*Ak_LVPWES_BSA+B(4)).*30.4167; % convert to days

% Qi et al 2015: "Anatomical and hemodynamic evaluations of the heart and pulmonary
% arterial pressure in healthy children residing at high altitude
% in China" 
% LVPWED
    Qi_th_timePoints = [0.5, 3, 9.5, 24.5]; % months (approximate)
    Qi_th_timePoints = Qi_th_timePoints*30.4167;
    Qi_LVED_th_targets = [3, 3, 3, 3];
    Qi_LVED_th_SD = [1, 1, 1, 1];
% LVPWES
    Qi_LVES_th_targets = [5, 6, 6, 7];
    Qi_LVES_th_SD = [1, 1, 1, 1];
% RVAWED
    Qi_RVED_th_targets = [2, 2, 2, 2];
    Qi_RVED_th_SD = [1, 1, 1, 1];

% Epstein et al 1974: "Great Vessel, Cardiac Chamber, and Wall Growth Patterns In Normal
% Children"
% LVPWED
    Ep_LVPWED_all = readmatrix('plottingData/epstein_LVPWED.csv');
    Ep_LVPWED_BSA = Ep_LVPWED_all(:,1);
    Ep_LVPWED_Age = (B(1).*Ep_LVPWED_BSA.^3 + B(2).*Ep_LVPWED_BSA.^2 + B(3).*Ep_LVPWED_BSA+B(4)).*30.4167; % convert to days
    Ep_LVPWED_50 = Ep_LVPWED_all(:,2);
    Ep_LVPWED_5 = Ep_LVPWED_all(:,3);
    Ep_LVPWED_95 = Ep_LVPWED_all(:,4);
% RVAWED
    Ep_RVAWED_all = readmatrix('plottingData/epstein_RVAWED.csv');
    Ep_RVAWED_BSA = Ep_RVAWED_all(:,1);
    Ep_RVAWED_Age = (B(1).*Ep_RVAWED_BSA.^3 + B(2).*Ep_RVAWED_BSA.^2 + B(3).*Ep_RVAWED_BSA+B(4)).*30.4167; % convert to days
    Ep_RVAWED_50 = (Ep_RVAWED_all(:,2)).*10; % convert to mm
    Ep_RVAWED_5 = (Ep_RVAWED_all(:,3)).*10;
    Ep_RVAWED_95 = (Ep_RVAWED_all(:,4)).*10;

figure(8);
set(gcf,'color','w');
subplot(2,2,1); % LVPWED
    Ka_LVPWED = plot(Ka_LVPWED_Age_50,Ka_LVPWED_50,'--','LineWidth',2,'Color',[0.75,0.75,0.75]); hold on;
    Ka_SD_LVPWED = plot(Ka_LVPWED_Age_90,Ka_LVPWED_90,':','LineWidth',2,'Color',[0.75,0.75,0.75]);
    plot(Ka_LVPWED_Age_10,Ka_LVPWED_10,':','LineWidth',2,'Color',[0.75,0.75,0.75]);
    Ak_LVPWED = plot(Ak_LVPWED_Age,Ak_LVPWED_50,'--','LineWidth',2,'Color',[0.5,0.5,0.5]);
    Qi_LVPWED = errorbar(Qi_th_timePoints, Qi_LVED_th_targets, Qi_LVED_th_SD, '.', 'color', [0.66 0.66 0.66], 'MarkerSize', 30);
    Ep_LVPWED = plot(Ep_LVPWED_Age,Ep_LVPWED_50,'--','LineWidth',2,'Color',[0.25,0.25,0.25]);
    Ep_SD_LVPWED = plot(Ep_LVPWED_Age,Ep_LVPWED_95,':','LineWidth',2,'Color',[0.25,0.25,0.25]);
    plot(Ep_LVPWED_Age,Ep_LVPWED_5,':','LineWidth',2,'Color',[0.25,0.25,0.25]);
    model_LVPWED = plot(age,dimensions(3,:),'LineWidth',2,'Color',modelColor);
    %title('LV ED Posterior Wall Thickness');
    xlabel('Age (Months)');
    ylabel('LV ED Thickness (mm)');
    xlim([0 1095]);
    ylim([0 10]);
    xticks([0 6*30.4167 12*30.4167 18*30.4167 24*30.4167 30*30.4167 36*30.4167]);
    xticklabels(ageTicksLong);
    set(gca,'FontSize',fontSize)
    legend_entries_LVPWED = [Ka_LVPWED,Ka_SD_LVPWED,Ak_LVPWED,Qi_LVPWED,Ep_LVPWED,Ep_SD_LVPWED,model_LVPWED];
    legend(legend_entries_LVPWED,{'Kampmann 2000','Kampmann 2000 (P10/90)','Akiba 1986','Qi 2015','Epstein 1974','Epstein 1974 (P5/95)','Model'},'Location','northwest','FontSize',12);
    legend boxoff

subplot(2,2,2); % LVPWES
    Ka_LVPWES = plot(Ka_LVPWES_Age,Ka_LVPWES_50,'--','LineWidth',2,'Color',[0.75,0.75,0.75]); hold on;
    Ka_SD_LVPWES = plot(Ka_LVPWES_Age,Ka_LVPWES_90,':','LineWidth',2,'Color',[0.75,0.75,0.75]);
    plot(Ka_LVPWES_Age,Ka_LVPWES_10,':','LineWidth',2,'Color',[0.75,0.75,0.75]);
    Ak_LVPWES = plot(Ak_LVPWES_Age,Ak_LVPWES_50,'--','LineWidth',2,'Color',[0.5,0.5,0.5]);
    Qi_LVPWES = errorbar(Qi_th_timePoints, Qi_LVES_th_targets, Qi_LVES_th_SD, '.', 'color', [0.66 0.66 0.66], 'MarkerSize', 30);
    model_LVPWES = plot(age,dimensions(4,:),'LineWidth',2,'Color',modelColor);
    %title('LV ES Posterior Wall Thickness');
    xlabel('Age (Months)');
    ylabel('LV ES Thickness (mm)');
    xlim([0 1095]);
    ylim([0 10]);
    xticks([0 6*30.4167 12*30.4167 18*30.4167 24*30.4167 30*30.4167 36*30.4167]);
    xticklabels(ageTicksLong);
    set(gca,'FontSize',fontSize)
    legend_entries_LVPWED = [Ka_LVPWES,Ka_SD_LVPWES,Ak_LVPWES,Qi_LVPWES,model_LVPWES];
    legend(legend_entries_LVPWED,{'Kampmann 2000','Kampmann 2000 (P10/90)','Akiba 1986','Qi 2015','Model'},'Location','southeast','FontSize',12);
    legend boxoff
    
subplot(2,2,3); % RVAWED
    Ka_RVAWED = plot(Ka_RVAWED_Age,Ka_RVAWED_50,'--','LineWidth',2,'Color',[0.75,0.75,0.75]); hold on;
    Ka_SD_RVAWED = plot(Ka_RVAWED_Age,Ka_RVAWED_90,':','LineWidth',2,'Color',[0.75,0.75,0.75]);
    plot(Ka_RVAWED_Age,Ka_RVAWED_10,':','LineWidth',2,'Color',[0.75,0.75,0.75]); 
    Qi_RVAWED = errorbar(Qi_th_timePoints, Qi_RVED_th_targets, Qi_RVED_th_SD, '.', 'color', [0.66 0.66 0.66], 'MarkerSize', 30);
    Ep_RVAWED = plot(Ep_RVAWED_Age,Ep_RVAWED_50,'--','LineWidth',2,'Color',[0.25,0.25,0.25]);
    Ep_SD_RVAWED = plot(Ep_RVAWED_Age,Ep_RVAWED_95,':','LineWidth',2,'Color',[0.25,0.25,0.25]);
    plot(Ep_RVAWED_Age,Ep_RVAWED_5,':','LineWidth',2,'Color',[0.25,0.25,0.25]);
    model_RVAWED = plot(age,dimensions_R(3,:),'LineWidth',2,'Color',modelColor);
    % title('RV ED Anterior Wall Thickness');
    xlabel('Age (Months)');
    ylabel('RV ED Thickness (mm)');
    xlim([0 1095]);
    ylim([0 10]);
    xticks([0 6*30.4167 12*30.4167 18*30.4167 24*30.4167 30*30.4167 36*30.4167]);
    xticklabels(ageTicksLong);
    set(gca,'FontSize',fontSize)
    legend_entries_RVAWED = [Ka_RVAWED,Ka_SD_RVAWED,Qi_RVAWED,Ep_RVAWED, Ep_SD_RVAWED,model_RVAWED];
    legend(legend_entries_RVAWED,{'Kampmann 2000','Kampmann 2000 (P10/90)','Qi 2015','Epstein 1974','Epstein 1974 (P5/95)','Model'},'Location','northwest','FontSize',12);
    legend boxoff
    
subplot(2,2,4); % RVAWES
% no literature data
    plot(age,dimensions_R(4,:),'LineWidth',2,'Color',modelColor); hold on;
    % title('RV ES Anterior Wall Thickness');
    xlabel('Age (Months)');
    ylabel('RV ES Thickness (mm)');
    xlim([0 1095]);
    ylim([0 10]);
    xticks([0 6*30.4167 12*30.4167 18*30.4167 24*30.4167 30*30.4167 36*30.4167]);
    xticklabels(ageTicksLong);
    set(gca,'FontSize',fontSize)
    legend('Model','Location','northwest','FontSize',14);
    legend boxoff

%% Figure 10
% (A) LV ED thickness, normal vs CoA
% (B) LV ED volume, normal vs CoA
% (C) SBP, normal vs CoA
% (D) DBP, normal vs CoA



figure(10)
set(gcf,'color','w');

subplot(3,2,1) % LVED_th
    plot(age,dimensions(3,:),'LineWidth',2,'Color',modelColor); hold on;

subplot(3,2,2); % LVED vol
    plot(age,ValuesofInterest(5,:),'LineWidth',2,'Color',modelColor); hold on;

% Systolic
%Data from Vogt et al 2005: Impaired elastic properties of the ascending aorta in newborns before and early after successful coarctation repair: Proof of a systemic vascular disease of the prestenotic arteries?
%control
subplot(3,2,3:4); hold on
    errorbar(13,82, 13,  '.', 'color', [0.5 0.5 0.5], 'linewidth',2, 'MarkerSize', 30);
    
    %coarc
    errorbar(20, 91, 22,  's', 'color', [0.5 0.5 0.5],'linewidth',2, 'MarkerSize', 10,'MarkerFaceColor',[0.5 0.5 0.5]);
   
%Data from Eerola et al 2007: Left ventricular hypertrophy persists after successful treatment for coarctation of the aorta
    %control
    %errorbar(2.94*365,97, 97-76, 120-97,  '.', 'color', [0.25 0.25 0.25], 'linewidth',2,'MarkerSize', 30);
    errorbar(2.94*365,97, 11,  '.', 'color', [0.25 0.25 0.25], 'linewidth',2,'MarkerSize', 30);

    %coarc
    %errorbar(2.00*365,123, 123-91, 151-123,   's', 'color', [0.25 0.25 0.25], 'linewidth',2,'MarkerSize', 10,'MarkerFaceColor',[0.25 0.25 0.25]);
    errorbar(2.00*365,123, 15,   's', 'color', [0.25 0.25 0.25], 'linewidth',2,'MarkerSize', 10,'MarkerFaceColor',[0.25 0.25 0.25]);

% Diastolic
subplot(3,2,5:6); hold on
    %Data from Vogt: Impaired elastic properties of the ascending aorta in newborns before and early after successful coarctation repair: Proof of a systemic vascular disease of the prestenotic arteries?
    %control
    errorbar(13,50, 9,  '.', 'color', [0.5 0.5 0.5], 'linewidth',2,'MarkerSize', 30);
    
    %coarc
    errorbar(20, 52, 14 ,  's', 'color', [0.5 0.5 0.5], 'linewidth',2,'MarkerSize', 10,'MarkerFaceColor',[0.5 0.5 0.5]);
    
    %Data from Eerola: Left ventricular hypertrophy persists after successful treatment for coarctation of the aorta
    %control
    %errorbar(2.94*365,56, 56-42, 69-56,  '.', 'color', [0.25 0.25 0.25], 'linewidth',2,'MarkerSize', 30);
    errorbar(2.94*365,56, 6.75,  '.', 'color', [0.25 0.25 0.25], 'linewidth',2,'MarkerSize', 30);
    
    %coarc
    %errorbar(2.00*365,65, 65-42, 96-65,  's', 'color', [0.25 0.25 0.25], 'linewidth',2,'MarkerSize', 10,'MarkerFaceColor',[0.25 0.25 0.25]);
    errorbar(2.00*365,65, 13.5,  's', 'color', [0.25 0.25 0.25], 'linewidth',2,'MarkerSize', 10,'MarkerFaceColor',[0.25 0.25 0.25]);


sys_S = zeros(1,numel(VolumeArray));
sys_D = zeros(1,numel(VolumeArray));

for i = 1:numel(VolumeArray)
        % Calculate systolic, diastolic, and mean BP
        sys_S(i) = max(PressureArray{i}(:,6));
        sys_D(i) = min(PressureArray{i}(:,6));
end
   
% title('Systemic Blood Pressure: Systolic');
subplot(3,2,3:4)
    plot(age,sys_S,'LineWidth',2,'Color',modelColor);
    
     
subplot(3,2,5:6); % systemic diastolic
    plot(age,sys_D,'LineWidth',2,'Color',modelColor);
    %title('Systemic Blood Pressure: Diastolic');
    

load('sampleOutput/coarcOutput.mat')

subplot(3,2,1) % LVED_th
    plot(age,dimensions(3,:),'LineWidth',2,'Color',modelColor_coarc); hold on;
    xlabel('Age (months)');
    xlim([0 1095]);
    xticks([0 6*30.4167 12*30.4167 18*30.4167 24*30.4167 30*30.4167 36*30.4167]);
    xticklabels(ageTicksLong);
    ylabel('LV ED Thickness (mm)');
    legend({'Normal','Coarctation'},'Location','southeast','FontSize',12);
    set(gca,'FontSize',12)
    legend boxoff

subplot(3,2,2); % LVED vol
    plot(age,ValuesofInterest(5,:),'LineWidth',2,'Color',modelColor_coarc); hold on;
    xlabel('Age (months)');
    xlim([0 1095]);
    xticks([0 6*30.4167 12*30.4167 18*30.4167 24*30.4167 30*30.4167 36*30.4167]);
    xticklabels(ageTicksLong);
    ylabel('LV ED Volume (mL)');
    legend({'Normal','Coarctation'},'Location','southeast','FontSize',12);
    set(gca,'FontSize',12)
    legend boxoff

sys_S = zeros(1,numel(VolumeArray));
sys_D = zeros(1,numel(VolumeArray));

for i = 1:numel(VolumeArray)
        % Calculate systolic, diastolic, and mean BP
        sys_S(i) = max(PressureArray{i}(:,6));
        sys_D(i) = min(PressureArray{i}(:,6));
end

subplot(3,2,3:4)
    plot(age,sys_S,'LineWidth',2,'Color',modelColor_coarc);
    xlabel('Age (months)');
    xlim([0 1095]);
    xticks([0 6*30.4167 12*30.4167 18*30.4167 24*30.4167 30*30.4167 36*30.4167]);
    xticklabels(ageTicksLong);
    ylabel('Systemic Systolic Pressure (mmHg)');
    ylim([30 155]);
    legend({'Vogt 2005: Control','Vogt 2005: Coarc','Eerola 2007: Control','Eerola 2007: Coarc','Normal','Coarctation'},'Location','eastoutside');
    set(gca,'FontSize',12)
    legend boxoff

subplot(3,2,5:6)
    plot(age,sys_D,'LineWidth',2,'Color',modelColor_coarc);
    xlabel('Age (months)');
    xlim([0 1095]);
    xticks([0 6*30.4167 12*30.4167 18*30.4167 24*30.4167 30*30.4167 36*30.4167]);
    xticklabels(ageTicksLong);
    ylabel('Systemic Diastolic Pressure (mmHg)');
    ylim([30 155]);
    legend({'Vogt 2005: Control','Vogt 2005: Coarc','Eerola 2007: Control','Eerola 2007: Coarc','Normal','Coarctation'},'Location','eastoutside');
    set(gca,'FontSize',12)
    legend boxoff
