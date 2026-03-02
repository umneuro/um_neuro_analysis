%% MPPCA-ICA code,
function [res_MPPCA] = fxn_MPPCA_ICAv11(result_ca_filt_data_ct, result_data_cell, prms_MPPCA)
%% update info
% 211104, changed xticks angle. 
% 220629, Added parameter decomposing section
% 220721, customization for multiple combination cal, speed-up cal, and decreased RAM usage by Kareem request 
% 220805, customization for more frexible analysis method for USM request.
%% Read grobal
global shift
%% parameter decomposing section
FPS_Ca                   = prms_MPPCA.sample_fps;       % RGECO: 10FPS, G-CaMP:20hz. Dependent on your data fps.
bin_frame_num            = prms_MPPCA.bin_frame_num; 
prms_forced_IC           = prms_MPPCA.prms_forced_IC;    % def=0;    0:MP-PCA, 1-x:forced_IC num input;  
prms_ICA_mode            = prms_MPPCA.prms_ICA_mode;     % def=3;    1:fastica, 2:fastICA, 3:Reconstruction ICA (RICA), 2 or 3 is better.
prms_SD_thr              = prms_MPPCA.prms_SD_thr;       % def=2.5   (SD) ensemble weight threshold
prms_reactivation_SD_thr = prms_MPPCA.prms_reactivation_SD_thr; % def=2;    (SD) filled coloring threshold for reactivation strength 
prms_target_region1      = prms_MPPCA.prms_target_region1 ;     % yellow, input second scale. like as [a:b]. disable: [0:0] 
prms_target_region2      = prms_MPPCA.prms_target_region2 ;     % green, input second scale. like as [c:d]. disable: [0:0]
prms_ticklabel_mode      = prms_MPPCA.prms_ticklabel_mode ;     % def=0;    0:off mode. 1:session name on mode. when you select whole data as target reference

reference_session_num = prms_MPPCA.prms_reference_session_num ; %
extra_session_mode    = prms_MPPCA.prms_extra_session_mode    ; % 0:backpreoject to all session, 1:enable backprojection prior to extra session.
extra_session_num     = prms_MPPCA.prms_extra_session_num     ; % input 1st session num to exclude for backprojection. Latter session also will be removed. 
prms_figure_all_mode  = prms_MPPCA.prms_figure_all_mode       ; % 0:disable figure all, 1:enable figure all.
%% for debug mode
% % clc; clear; close all; addpath('FastICA_25\','-end'); addpath('pca_ica\pca_ica\','-end');
% % load('demo_ca_data.mat'); 
% % data_temp = demo_ca_data; % for debug;
% % 
% % bin_frame_num = 40;
% % FPS_Ca = 20;
% % [~,~,ca_mod_round_bin] = fxn_mod_round_binning_time(data_temp, bin_frame_num);
% % data_bin_z = zscore(ca_mod_round_bin(1:60,:)); % 
% % data_bin_z_full = zscore(ca_mod_round_bin(:,:));
% % 
% % % Parameters
% % prms_forced_IC              = 5;        % def=0;    0:MP-PCA, 1-x:forced_IC num input;  
% % prms_ICA_mode               = 3;        % def=3;    1:fastica, 2:fastICA, 3:Reconstruction ICA (RICA), 2 or 3 is better.
% % prms_SD_thr                 = 2.5;      % def=2.5   (SD) ensemble weight threshold
% % prms_reactivation_SD_thr    = 2;        % def=2;    (SD) filled coloring threshold for reactivation strength 
% % prms_interest_region1       = [60:120] ;   % yellow, input second scale. like as [a:b]. disable: [0:0] 
% % prms_interest_region2       = [400:480] ;   % green, input second scale. like as [c:d]. disable: [0:0]
%% # Assign reference and target section. Select either this section or below, (Disable eiteher here or below)
% # reference vs. all sesion analysis mode. I recomend this first to check data aspect.

data_bin_reference = result_data_cell{reference_session_num+shift,7}; % input data for MPPCA-ICA
data_bin_all_ct = result_ca_filt_data_ct;

% ### assign target session ### 
if extra_session_mode == 0
% # If you backprojection to all session, enable thie line.
disp('Backprojected to all sessions')
data_bin_target = result_ca_filt_data_ct; 
elseif extra_session_mode == 1
disp('Restricted backprojection')
% # If you backprojection to selected session, enable thie line.
data_bin_target = result_ca_filt_data_ct(1:result_data_cell{extra_session_num , 11},:);
end
%% for addpath and ticks
shift =1;
addpath('FastICA_25\','-end'); addpath('pca_ica\pca_ica\','-end');
%% PCA
[pca_coeff, pca_score, pca_latent, pca_tsquared, pca_explained_prop] = pca(data_bin_reference);
pca_eigenvalue = pca_latent; % eigenvalue
%% MP-PCA
thrcov_PC_percnet_indicator = 70; % blue line in fig will visualize the boarder of PCs 
[above_MP_PCA_num, res_MPPCA] = fxn_Marchenko2PCA(data_bin_reference, thrcov_PC_percnet_indicator);
disp([' ', num2str(above_MP_PCA_num), ' PCs are above Marcenko-Pastur threshold.']);
%% Decompose of Correlation matrix
if prms_forced_IC == 0
    i_range = 1:above_MP_PCA_num; disp('MP-PCA mode.');
else
    i_range = 1:prms_forced_IC; disp('Forced IC')
end

% cell_ith_decompose = {};
cell_ith_decompose = {i_range(end),1}; % initialize to speed-up!

% figure;
for i = 1:i_range(end)
ith_Lambda = pca_eigenvalue(i);
ith_eigenvector = pca_coeff(:,i);
ith_decompose =  ith_Lambda .* ith_eigenvector .* ith_eigenvector';
cell_ith_decompose{i} = ith_decompose;
end
%%
Psign = pca_coeff(:,i_range); % correct
Z_proj = Psign' * data_bin_reference';
%% ICA
rng default
r = i_range(end);
tic;
if prms_ICA_mode == 1
    disp('You selected fastica.')
    [A, W] = fastica(Z_proj, 'numOfIC', r); % fastica
elseif prms_ICA_mode == 2
    disp('You selected fastICA.')
    [Zica, W, T, mu] = fastICA(Z_proj,r); % fastICA
elseif prms_ICA_mode == 3
    disp('You selected Reconstruct ICA.')
    Mdl = rica(Z_proj, r, 'IterationLimit',10000,'Standardize',false ,'VerbosityLevel',1); % RICA 
    W = transform(Mdl,Z_proj);
end
toc; disp('ICA cal finished!');
%%
V = Psign * W;
V_sqr = V.^2;
V_sqr_scale = (V_sqr./sum(V_sqr,1));
V_sqr_scale_sum = sum(V_sqr_scale,1); % sum check
V_sqr_scale_z = zscore(V_sqr_scale);

negative2 = find(V_sqr_scale_z < prms_SD_thr);
positive2 = find(V_sqr_scale_z >= prms_SD_thr);
V_sqr_scale_z_cutoff = V_sqr_scale_z;
V_sqr_scale_z_cutoff(negative2) = 0; 
V_sqr_scale_z_cutoff(positive2) = 1;
V_sqr_scale_z_cutoff_sum2 = sum(V_sqr_scale_z_cutoff,2);
%% Tracking the activation-strength of assembly patterns over time
disp([' ', num2str(r), ' ICs are going to be calculated.']);
% cell_p_k_diag_zero = {};
cell_p_k_diag_zero = cell(r,1); % initialize, speed-up
r_kt      = zeros(size(data_bin_reference,1),r); % initialize, speed-up
r_kt_full = zeros(size(data_bin_reference,1),r); % initialize, speed-up

for i_r = 1:r % pattern_vector_id = 1; % 1:r 
    for i_t = 1:size(data_bin_reference,1)% time
    z_t = data_bin_reference(i_t,:);
    pattern_vector = V_sqr_scale(:,i_r); % 1 -> r
    p_k = pattern_vector * pattern_vector';
    p_k_diag_zero = p_k;
        for i = 1:size(p_k,1)
        p_k_diag_zero(i,i) = 0;
        end
    cell_p_k_diag_zero{i_r} = p_k_diag_zero;
    r_kt(i_t,i_r) = z_t * p_k_diag_zero * z_t'; %
    end
    for i_t = 1:size(data_bin_target,1)% time
    z_t_full = data_bin_target(i_t,:);
    r_kt_full(i_t,i_r) = z_t_full * p_k_diag_zero * z_t_full'; %
    end    
%     figure; imagesc(p_k_diag_zero); % for debug;
    disp(['Pattern# ', num2str(i_r),'  out of total ',num2str(r),' Rx(t)s is calculated.']);
end

for i = 1:r
    V3_sort_column = V_sqr_scale_z_cutoff(:,i);
    V3_sort_id = find(V3_sort_column == 1);
    neuron_sig_IDs{i} =  V3_sort_id;
end

%%
% ### start. for figuring off
%%
mkdir result_MPPCA_mat_figs
%% reference data fig, ensemble sorting for check
if prms_figure_all_mode == 0
r = 1; disp('Reduce and speed-up calculation. Off to show all figures.')
elseif prms_figure_all_mode == 1
    disp('Figure all')
% for i = 1:1 % for debug
end

% ### figure cal here ###
for i = 1:r

H1 = figure('Position',[100,250,700,700]); %[left bottom width height] 

% h1 = subaxis(6,1,1, 'SpacingVert',0.06, 'MR',0.15); 
h1 = subplot(611);
stem(V_sqr_scale_z(:,i),'k','MarkerFaceColor','k','MarkerSize', 2); hold on
title(['Reference MPPCA-ICA ensemble #', num2str(i), ' of ' ,num2str(r)])

V3_sqr_scale_z_color = V_sqr_scale_z;
V3_sqr_scale_z_color(negative2) = NaN; 
stem(V3_sqr_scale_z_color(:,i),'r','MarkerFaceColor','r','MarkerSize', 4); 
V3_sqr_scale_z_color_dash = ones(size(data_bin_reference,2),1)*prms_reactivation_SD_thr;
plot(V3_sqr_scale_z_color_dash, 'r--');
xlabel('Neurons'); ylabel({'Neurons weight in';'in assemble (SD)'}); xlim([0 size(data_bin_reference,2)]); % legend below-thr above-thr SD-thr
hold off

V3_sort_column = V_sqr_scale_z_cutoff(:,i);
V3_sort_id = find(V3_sort_column == 1);
V3_sort_data = data_bin_reference(:,V3_sort_id);
fig_time = [1:size(V3_sort_data,1)]/(FPS_Ca/bin_frame_num);
fig_num_cell = [1:size(V3_sort_data,2)];

% h2 = subaxis(6,1,2,'SpacingVert',0.01) ;
h2 = subplot(612);

plot(fig_time, r_kt(:,i),'b'); grid on;
ylabel({'Reactivation';'strength (AU)'}); xlim([0 fig_time(end)]); %ylim([0 1])
xticklabels([]); hold on

% prms_target_region1 = prms_target_region1 *(FPS_Ca/bin_frame_num);
% interest_shadow1 = (ones(prms_target_region1(end)-prms_target_region1(1)+1,1)*max(r_kt(:,i)))';

% interest_shadow1map = (ones(prms_target_region1(end)-prms_target_region1(1)+1,1)*(size(V3_sort_data,2)+0.5))';
% area(prms_target_region1, interest_shadow1,'FaceColor','y','FaceAlpha',.1,'EdgeAlpha',.1);

% prms_target_region2 = prms_target_region2 *(FPS_Ca/bin_frame_num);
% interest_shadow2 = (ones(prms_target_region2(end)-prms_target_region2(1)+1,1)*max(r_kt(:,i)))';

% interest_shadow2map = (ones(prms_target_region2(end)-prms_target_region2(1)+1,1)*(size(V3_sort_data,2)+0.5))';
% area(prms_target_region2, interest_shadow2,'FaceColor','g','FaceAlpha',.1,'EdgeAlpha',.1);

hold off

r_kt_z = zscore(r_kt);
r_kt_z_above_SD = find(r_kt_z >= prms_reactivation_SD_thr);
r_kt_z_below_SD = find(r_kt_z < prms_reactivation_SD_thr);

r_kt_z_color_above_thr = r_kt_z;  
r_kt_z_color_above_thr(r_kt_z_below_SD) = NaN;
r_kt_z_color_above_thr_dash = ones(size(data_bin_reference,1),1)*prms_reactivation_SD_thr;

interest_shadow1sd = (ones(prms_target_region1(end)-prms_target_region1(1)+1,1)*max(r_kt_z(:,i)))';
interest_shadow2sd = (ones(prms_target_region2(end)-prms_target_region2(1)+1,1)*max(r_kt_z(:,i)))';

% h3 = subaxis(6,1,3); 
h3 = subplot(613);
plot(fig_time, r_kt_z(:,i),'k'); hold on; grid on;
area(fig_time, r_kt_z_color_above_thr(:,i),'FaceColor','r','EdgeColor','r'); xticklabels([]); hold on
plot(fig_time, r_kt_z_color_above_thr_dash, 'r--'); 
ylabel({'Reactivation';'strength (SD)'}); xlim([0 fig_time(end)]); 
% area(prms_target_region1, interest_shadow1sd,'FaceColor','y','FaceAlpha',.1,'EdgeAlpha',.1);
% area(prms_target_region2, interest_shadow2sd,'FaceColor','g','FaceAlpha',.1,'EdgeAlpha',.1);
hold off;

% h4 = subaxis(6,1,4:6);
h4 = subplot(6,1,4:6);
imagesc(fig_time, fig_num_cell, V3_sort_data'); ylabel({'Pattern-related';'cell activities (z-scored)'}); hold on;
caxis([-3 3]); colormap(fxn_redblue); xlabel('Time (s)'); 
% c= colorbar; c.Label.String = 'z-score'; % for fig
% area(prms_target_region1, interest_shadow1map,'FaceColor','y','FaceAlpha',.15,'EdgeAlpha',.5);
% area(prms_target_region2, interest_shadow2map,'FaceColor','g','FaceAlpha',.15,'EdgeAlpha',.5);

set(findobj(gcf, 'Type', 'Axes'), 'FontSize', 9, 'FontName','Arial');
linkaxes([h2,h3,h4], 'x');

filename1 = ('ref_');
filename2 = num2str(i);
save_extention = ('_.fig');
filename_ref = [filename1, filename2, save_extention];

cd result_MPPCA_mat_figs
savefig(H1,filename_ref,'compact')
cd ..
% close all
end
%% target data fig, ensemble sorting for check
if prms_figure_all_mode == 0
r = 1; disp('Reduce and speed-up calculation. Off to show all figures.')
elseif prms_figure_all_mode == 1
    disp('Figure all')
% for i = 1:1 % for debug
end

% ### figure cal here ###
for i = 1:r

H2 = figure('Position',[800,250,700,700]); %[left bottom width height] 

% h1 = subaxis(6,1,1, 'SpacingVert',0.06, 'MR',0.15); 
h5 = subplot(711);
stem(V_sqr_scale_z(:,i),'k','MarkerFaceColor','k','MarkerSize', 2); hold on
title(['Backprojected MPPCA-ICA ensemble #', num2str(i), ' of ' ,num2str(r)])

V3_sqr_scale_z_color = V_sqr_scale_z;
V3_sqr_scale_z_color(negative2) = NaN; 
stem(V3_sqr_scale_z_color(:,i),'r','MarkerFaceColor','r','MarkerSize', 4); 
V3_sqr_scale_z_color_dash = ones(size(data_bin_target,2),1)*prms_reactivation_SD_thr;
plot(V3_sqr_scale_z_color_dash, 'r--');
xlabel('Neurons'); ylabel({'Neurons weight in';'in assemble (SD)'}); xlim([0 size(data_bin_target,2)]); % legend below-thr above-thr SD-thr
hold off

V3_sort_column = V_sqr_scale_z_cutoff(:,i);
V3_sort_id = find(V3_sort_column == 1);
V3_sort_id_data{i} =  V3_sort_id;
V3_sort_data = data_bin_target(:,V3_sort_id);
fig_time = [1:size(V3_sort_data,1)]/(FPS_Ca/bin_frame_num);
fig_num_cell = [1:size(V3_sort_data,2)];

% h2 = subaxis(6,1,2,'SpacingVert',0.01) ;
h6 = subplot(712); 
plot(fig_time, r_kt_full(:,i),'b'); ylabel({'Reactivation';'strength (AU)'}); xlim([0 fig_time(end)]); %ylim([0 1])
xticklabels([]); hold on; grid on;

prms_target_region1 = prms_target_region1*(FPS_Ca/bin_frame_num);
interest_shadow1 = (ones(prms_target_region1(end)-prms_target_region1(1)+1,1)*max(r_kt_full(:,i)))';

interest_shadow1map = (ones(prms_target_region1(end)-prms_target_region1(1)+1,1)*(size(V3_sort_data,2)+0.5))';
% area(prms_target_region1, interest_shadow1,'FaceColor','y','FaceAlpha',.1,'EdgeAlpha',.1);

prms_target_region2 = prms_target_region2*(FPS_Ca/bin_frame_num);
interest_shadow2 = (ones(prms_target_region2(end)-prms_target_region2(1)+1,1)*max(r_kt_full(:,i)))';

interest_shadow2map = (ones(prms_target_region2(end)-prms_target_region2(1)+1,1)*(size(V3_sort_data,2)+0.5))';
% area(prms_target_region2, interest_shadow2,'FaceColor','g','FaceAlpha',.1,'EdgeAlpha',.1);

if prms_ticklabel_mode == 0
%     disp('ticklabel mode OFF')
elseif prms_ticklabel_mode == 1
%     disp('ticklabel mode ON')
    xticks([]); % update 211104 to off display
% xticklabels(result_data_cell(1+shift:size(result_data_cell,1),1) );  % for adding ticklabel
% xticks([1; cell2mat(result_data_cell(1+shift:size(result_data_cell,1),9))]); % for adding ticklabel
% xtickangle(45) % update 211104
end

hold off

r_kt_full_z = zscore(r_kt_full);
r_kt_full_z_above_SD = find(r_kt_full_z >= prms_reactivation_SD_thr);
r_kt_full_z_below_SD = find(r_kt_full_z < prms_reactivation_SD_thr);

r_kt_full_z_color_above_thr = r_kt_full_z;  
r_kt_full_z_color_above_thr(r_kt_full_z_below_SD) = NaN;
r_kt_full_z_color_above_thr_dash = ones(size(data_bin_target,1),1)*prms_reactivation_SD_thr;

interest_shadow1sd = (ones(prms_target_region1(end)-prms_target_region1(1)+1,1)*max(r_kt_full_z(:,i)))';
interest_shadow2sd = (ones(prms_target_region2(end)-prms_target_region2(1)+1,1)*max(r_kt_full_z(:,i)))';

% h3 = subaxis(6,1,3); 
h7 = subplot(7,1,3:4);
plot(fig_time, r_kt_full_z(:,i),'k'); hold on; grid on;
area(fig_time, r_kt_full_z_color_above_thr(:,i),'FaceColor','r','EdgeColor','r'); xticklabels([]); hold on
plot(fig_time, r_kt_full_z_color_above_thr_dash, 'r--'); 
ylabel({'Reactivation';'strength (SD)'}); xlim([0 fig_time(end)]); 
% area(prms_target_region1, interest_shadow1sd,'FaceColor','y','FaceAlpha',.1,'EdgeAlpha',.1);
% area(prms_target_region2, interest_shadow2sd,'FaceColor','g','FaceAlpha',.1,'EdgeAlpha',.1);

if prms_ticklabel_mode == 0
%     disp('ticklabel mode OFF')
elseif prms_ticklabel_mode == 1
%     disp('ticklabel mode ON')
xticklabels(result_data_cell(1+shift:size(result_data_cell,1),1) );  % for adding ticklabel
% xticks([0; cell2mat(result_data_cell(1+shift:size(result_data_cell,1),10))]); % for adding ticklabel
xticks([cell2mat(result_data_cell(1+shift:size(result_data_cell,1),10))]); % for adding ticklabel
xtickangle(45) % update 211104
end

hold off;

% h4 = subaxis(6,1,4:6);
h8 = subplot(7,1,6:7);
fig_time = [1:size(V3_sort_data,1)]/(FPS_Ca/bin_frame_num);
fig_num_cell = [1:size(V3_sort_data,2)];
imagesc(fig_time, fig_num_cell, V3_sort_data'); ylabel({'Pattern-related';'cell activities (z-scored)'}); hold on;
caxis([-3 3]); colormap(fxn_redblue); xlabel('Time (s)'); 
% c= colorbar; c.Label.String = 'z-score'; % for fig
% area(prms_target_region1, interest_shadow1map,'FaceColor','y','FaceAlpha',.15,'EdgeAlpha',.5);
% area(prms_target_region2, interest_shadow2map,'FaceColor','g','FaceAlpha',.15,'EdgeAlpha',.5);

set(findobj(gcf, 'Type', 'Axes'), 'FontSize', 8, 'FontName','Arial');
linkaxes([h6,h7,h8], 'x');

filename1 = ('backprop_');
filename2 = num2str(i);
save_extention = ('_.fig');
filename_ref = [filename1, filename2, save_extention];
% mkdir result_figs
cd result_MPPCA_mat_figs
savefig(H2, filename_ref,'compact')
cd ..
% close all
end
%%
% ### end. for figuring off
%%

%% update 211104, add r_kt_full-spilitted data 
result_data_cell{1,13} = ('r_kt'); result_data_cell{1,14} = ('r_kt z'); 


% ### assign target session ### 
if extra_session_mode == 0
    % # If you backprojection to all session, enable thie line.
disp('Backprojected to all sessions')
    for i = 1:size(result_data_cell,1)-1
    result_data_cell{i+shift,13} =  r_kt_full(result_data_cell{i+shift,12},:);
    result_data_cell{i+shift,14} =  r_kt_full_z(result_data_cell{i+shift,12},:);
    end
elseif extra_session_mode == 1
disp('Restricted backprojection')
% # If you backprojection to selected session, enable thie line.
    for i = 1:extra_session_num-1
    result_data_cell{i+shift,13} =  r_kt_full(result_data_cell{i+shift,12},:);
    result_data_cell{i+shift,14} =  r_kt_full_z(result_data_cell{i+shift,12},:);
    end
end

%%
disp('All calculation are finished!')
%% Remind mode
disp('You selected below modes.')
if prms_forced_IC == 0
    disp('   *MP-PCA mode.');
else
    disp(['   *Forced IC mode, ', num2str(prms_forced_IC), ' ICs'])
end
if prms_ICA_mode == 1
    disp('   *fastica mode.')
elseif prms_ICA_mode == 2
    disp('   *fastICA mode.')
elseif prms_ICA_mode == 3
    disp('   *Reconstruct ICA mode.')
end
toc;
%% save data
res_MPPCA.ith_decompose        = ith_decompose;
res_MPPCA.r_strength_ref       = r_kt;
res_MPPCA.r_strength_refz      = r_kt_z;
res_MPPCA.r_strength_target    = r_kt_full;
res_MPPCA.r_strength_targetz   = r_kt_full_z;
res_MPPCA.neuron_weight        = V_sqr_scale;
res_MPPCA.neuron_weight_sigID  = V_sqr_scale_z_cutoff;
% res_PCAICA.neuron_sig_IDs       = V3_sort_id_data;
res_MPPCA.neuron_sig_IDs       = neuron_sig_IDs;
res_MPPCA.result_data_cell     = result_data_cell;

res_MPPCA.V_sqr_scale_z                  = V_sqr_scale_z;
res_MPPCA.r                              = r;
res_MPPCA.prms_reactivation_SD_thr       = prms_reactivation_SD_thr;
res_MPPCA.result_data_cell               = result_data_cell;
res_MPPCA.data_bin_z_full                = data_bin_target;
res_MPPCA.data_bin_z                     = data_bin_reference;
res_MPPCA.negative2                      = negative2;
res_MPPCA.bin_frame_num                  = bin_frame_num;
res_MPPCA.FPS_Ca                         = FPS_Ca;
res_MPPCA.prms_MPPCA                     = prms_MPPCA;

res_MPPCA.data_bin_all_ct    = data_bin_all_ct;
res_MPPCA.data_bin_reference = data_bin_reference;
res_MPPCA.data_bin_target    = data_bin_target; 
%%
filename_mat1 = ('result_MPPCA_ref');
filename_mat2 = num2str(reference_session_num);
filename_mat3 = ('.mat');
filename_mat  = [filename_mat1,filename_mat2,filename_mat3];
% filename_ref = [filename1, filename2, save_extention];
% mkdir result_figs
cd result_MPPCA_mat_figs
save(filename_mat, 'res_MPPCA','-v7.3');
cd ..
disp('Finish file saving!')
%%
% close all
%%
if extra_session_mode == 0
disp('Backprojected to all sessions')
disp(['Reference frame num is     ',num2str(size(data_bin_reference,1)),' frames.'])
disp(['Entire frame num is        ',num2str(size(data_bin_all_ct,1)),' frames.'])
disp(['Backprojected frame num is ',num2str(size(data_bin_target,1)),' frames.'])

elseif extra_session_mode == 1
disp('Restricted backprojection')
disp(['Reference frame num is     ',num2str(size(data_bin_reference,1)),' frames.'])
disp(['Entire frame num is        ',num2str(size(data_bin_all_ct,1)),' frames.'])
disp(['Backprojected frame num is ',num2str(size(data_bin_target,1)),' frames.'])
end
%%
end