function  [ dimensions, ValuesofInt, dimensions_R, ValuesofInt_R]= computing_dimensionsandhemodynamics_AAH_20200609(Pressures, Volumes, TimeVector, Valves,  r, h, growthstep,  iteration_growth_max, HLHS)
% computing_dimensionsandhemodynamics: calculates and saves ED and ES pressures/volumes/thicknesses 

% Inputs:
%   Pressures: compartment pressures throughout cardiac cycle
%   Volumes: compartment volumes throughout cardiac cycle
%   Valves: logical array showing whether valves are open or closed
%   r: loaded radius
%   h: loaded thickness
%   growstep: what "day" the simulation is currently on
%   iteration_growth_max: how many days to run the simulation
%   HLHS: flag for hypoplastic left heart syndrome

% Outputs:
%   dimensions: LV (1)EDV, (2)ESV, (3)ED thickness, (4)ES thickness (thicknesses converted to mm)
%   ValuesofInt: (1)EDP (2)EDV (3)MaxP (4)MAP (5)maxV (6)minV (7)RF (8)CO (9)HR (10)ESP (11)ESP
%   dimensions_R: RV (1)EDV, (2)ESV, (3)ED thickness, (4)ES thickness (thicknesses converted to mm)
%   ValuesofInt_R: (1)EDP (2)EDV (3)MaxP (4)MAP (5)maxV (6)minV (7)RF (8)CO (9)HR (10)ESP (11)ESP
        

%ValuesofInt: (1)EDP (2)EDV (3)MaxP (4)MAP (5)maxV (6)minV (7)RF (8)CO (9)HR (10)ESP (11)ESP
ValuesofInt=zeros(11,  1*logical(growthstep>0)+logical(growthstep<1)*iteration_growth_max);
ValuesofInt_R=zeros(11,  1*logical(growthstep>0)+logical(growthstep<1)*iteration_growth_max);
    
% end diastolic pressure
Valves(end+1,:)=Valves(1,:);
valvechange = diff(Valves); valvechange=circshift(valvechange,1);
if HLHS == 1
    row_MV_closes = 0;
    row_AV_closes = 0;
    EDP_model = 0;
    ESP_model = 0;
    ValuesofInt(2,1)=0;
    ValuesofInt(11,1)=0;
else
    row_MV_closes = find(valvechange(:,1)==-1);
    row_AV_closes = find(valvechange(:,2)==-1);
    EDP_model = Pressures(row_MV_closes,3);
    ESP_model = Pressures(row_AV_closes,3);
    ValuesofInt(2,1)=Volumes(row_MV_closes,3);
    ValuesofInt(11,1)=Volumes(row_AV_closes,3);
end
row_TV_closes = find(valvechange(:,3)==-1); % RV
EDP_model_R = Pressures(row_TV_closes,10);
ValuesofInt(1,1)=EDP_model;
ValuesofInt_R(1,1)=EDP_model_R;

% end systolic pressure
row_PV_closes = find(valvechange(:,4)==-1);
ESP_model_R = Pressures(row_PV_closes,10);
ValuesofInt(10,1)=ESP_model;
ValuesofInt_R(10,1)=ESP_model_R;
    
        
%EDV and Maximum volume    
ValuesofInt_R(2,1)=Volumes(row_TV_closes,10);
[ValuesofInt(5,1), loc]=max(Volumes(:,3));
[ValuesofInt_R(5,1), loc]=max(Volumes(:,10));
if abs(ValuesofInt(2,1)-ValuesofInt(5,1))>0.01
  disp('ERROR: left ventricular end-diastolic volume and maximum left ventricular volume are NOT the same');
end
EDV_model=ValuesofInt(2,1);
EDV_model_R=ValuesofInt_R(2,1);  

%Minimum volume
ValuesofInt(6,1)=min(Volumes(:,3));
ValuesofInt_R(6,1)=min(Volumes(:,10));
ESV_model=ValuesofInt(6,1);
ESV_model_R=ValuesofInt_R(6,1);

%ES Wall Thickness   
[minr, locminr]=min(r(:,1));
[minr_R, locminr_R]=min(r(:,2));
h_es_model=h(locminr,1);
h_es_model_R=h(locminr_R,2);
          
%ED Wall Thickness   
[maxr, locmaxr]=max(r(:,1));
[maxr_R, locmaxr_R]=max(r(:,2));
h_ed_model_R=h(locmaxr_R,2);
h_ed_model=h(locmaxr,1);
    
%MaxP
ValuesofInt(3,1)=max(Pressures(:,3));
ValuesofInt_R(3,1)=max(Pressures(:,10));

%LV Mass     
Vwall=(r(:,1)+h(:,1)).^3 - r(:,1).^3;
% RV Mass
Vwall_R=(r(:,2)+h(:,2)).^3 - r(:,2).^3;
    
    
%dimensions: EDV, ESV, ED thickness, ES thickness (thicknesses converted to mm)   
dimensions=zeros(5,   1*logical(growthstep>-1)+logical(growthstep<0)*iteration_growth_max);
dimensions(:,1)=[EDV_model, ESV_model, h_ed_model.*10, h_es_model.*10, Vwall(1)];
dimensions_R=zeros(5,   1*logical(growthstep>-1)+logical(growthstep<0)*iteration_growth_max);
dimensions_R(:,1)=[EDV_model_R, ESV_model_R, h_ed_model_R.*10, h_es_model_R.*10, Vwall_R(1)];

%ValuesofInterest are EDP(mmHg); EDV(ml); MaxP(mmHg); MAP(mmHg); MaxV(ml); MinV(ml); RF(%); CO(L/min); 
ValuesofInt(4,1)=mean(Pressures(:,5)+Pressures(:,6));
[RF, ForwardSV]=regurgfraction_AAH_20200601(Volumes, Pressures, TimeVector, Valves);
ValuesofInt(7,1)=RF; %regurgitant fraction, not needed for RV
ValuesofInt(8,1)=60/max(TimeVector)*ForwardSV; %cardiac output, not needed for RV   
ValuesofInt(9,1)=60/max(TimeVector);
ValuesofInt_R(9,1)=60/max(TimeVector);

% Actual ESV
ValuesofInt_R(11,1) = Volumes(row_PV_closes,10);
    

end

 