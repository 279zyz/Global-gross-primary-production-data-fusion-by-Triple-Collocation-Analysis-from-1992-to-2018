%% 1\ data preprcess
% build the mask of data of nan data
clc;clear
load etcgpp.mat gpp_ec_year gpp_nirv_year gpp_tl_year
mask = mean(gpp_ec_year,3)>20 &mean(gpp_nirv_year,3)>20 &mean(gpp_tl_year,3)>20;

save gpp005.mat mask -append

%% 2\ triple collocation analysis
% 1\ the cal of etc: errvar and rho
% 2\ cal the percentage of effective pixels
clc;clear

load gpp005.mat gpp_ec005 gpp_nirv005 gpp_tl005 mask

% Stretch the matrix into 2D
gpp_ec_2d = (ftx_mask_file(gpp_ec005,mask));
gpp_nirv_2d = (ftx_mask_file(gpp_nirv005,mask));
gpp_tl_2d = (ftx_mask_file(gpp_tl005,mask));

save gpp005.mat gpp_ec_2d gpp_nirv_2d gpp_tl_2d -append

% TC analysis
n_grid = sum(mask(:));% 4100523
errVar_ETC = nan(n_grid,3);
rho2_ETC = nan(n_grid,3);
for i_grid = 1:n_grid
    y = [gpp_ec_2d(i_grid,:)',gpp_nirv_2d(i_grid,:)',gpp_tl_2d(i_grid,:)'];
    [errVar_ETC(i_grid,:), rho2_ETC(i_grid,:)] = ETC(y);
    disp(i_grid)
end

errVar = ftx_remask_file(errVar_ETC,mask);
rho2 = ftx_remask_file(rho2_ETC,mask);

save gpp005.mat errVar rho2 -append

mask_effective_etc = errVar(:,:,1)>0 & errVar(:,:,2)>0 & errVar(:,:,3)>0; 
mask_invalid = mask~=mask_effective_etc;

percent_effective = sum(mask_effective_etc(:))/sum(mask(:)); % 90.54%
save gpp005.mat mask_effective_etc mask_invalid -append

%% 3 data fusion
% Least squares based data fusion
clc;clear

load gpp005.mat errVar gpp_ec005 gpp_nirv005 gpp_tl005 mask_effective_etc mask_invalid

e1 = errVar(:,:,1);e2 = errVar(:,:,2);e3 = errVar(:,:,3);

w1 = (e2.*e3)./((e1.*e2)+(e1.*e3)+(e2.*e3));
w2 = (e1.*e3)./((e1.*e2)+(e1.*e3)+(e2.*e3));
w3 = (e1.*e2)./((e1.*e2)+(e1.*e3)+(e2.*e3));

w1_effective = w1.*mask_effective_etc;
w2_effective = w2.*mask_effective_etc;
w3_effective = w3.*mask_effective_etc;

w1_invalid = (nansum(w1_effective(:))/sum(mask_effective_etc(:))).*mask_invalid;
w2_invalid = (nansum(w2_effective(:))/sum(mask_effective_etc(:))).*mask_invalid;
w3_invalid = (nansum(w3_effective(:))/sum(mask_effective_etc(:))).*mask_invalid;

w1_lsm = w1_effective+w1_invalid;
w2_lsm = w2_effective+w2_invalid;
w3_lsm = w3_effective+w3_invalid;

for i = 1:324
    gpp_tcf_month(:,:,i) = w1_lsm.*single(gpp_ec(:,:,i))+w2_lsm.*single(gpp_nirv(:,:,i))+w3_lsm.*single(gpp_tl(:,:,i));
    disp(i)
end
save gpp005.mat gpp_tcf_month w1_lsm w2_lsm w3_lsm -append

