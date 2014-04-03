function mask=autoMask(j_date,b_ref,yr,T_times)
%% Multitepmoral cloud, cloud shadow, & snow masks (global version)
% read in data with 3 more consecutive clear obs & correct data
% Inputs:
% j_date: Julian date
% b_ref: Band 2 & 4 reflectance
% yr: number of years to mask 
%
% Outputs:
% y: corrected Band reflectance

yr=ceil(yr);
w=2*pi/365.25; % anual cycle

TOA_B1_cft=autoRobustFit(j_date,b_ref(:,1),yr); % Band 2
TOA_B2_cft=autoRobustFit(j_date,b_ref(:,2),yr); % Band 5

% predict Band 1 ref
pred_B1=TOA_B1_cft(1)+TOA_B1_cft(2)*cos(j_date*w)+TOA_B1_cft(3)*sin(j_date*w)...
    +TOA_B1_cft(4)*cos(j_date*w/yr)+TOA_B1_cft(5)*sin(j_date*w/yr);

% predict Band 2 ref
pred_B2=TOA_B2_cft(1)+TOA_B2_cft(2)*cos(j_date*w)+TOA_B2_cft(3)*sin(j_date*w)...
    +TOA_B2_cft(4)*cos(j_date*w/yr)+TOA_B2_cft(5)*sin(j_date*w/yr);

% true or false ids
DeltaB1=(b_ref(:,1)-pred_B1);
DeltaB2=(b_ref(:,2)-pred_B2);

mask=zeros(length(j_date),1);
mask(DeltaB1 > T_times|DeltaB2 < -T_times)=1;
end
    


