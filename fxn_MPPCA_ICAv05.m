%% MPPCA-ICA code
function [res_PCAICA] = fxn_MPPCA_ICAv05(data_bin_z, data_bin_z_full, prms_forced_IC, prms_ICA_mode, prms_SD_thr,...
                                         prms_reactivation_SD_thr, prms_interest_region1, prms_interest_region2, ...
                                         bin_frame_num, result_data_cell, prms_ticklabel_mode)

%% for addpath and ticks
shift =1;
addpath('FastICA_25\','-end'); addpath('pca_ica\pca_ica\','-end');
%% PCA
[pca_coeff, pca_score, pca_latent, pca_tsquared, pca_explained_prop] = pca(data_bin_z);
pca_eigenvalue = pca_latent; % eigenvalue
%% MP-PCA
thrcov_PC_percnet_indicator = 70; % blue line in fig will visualize the boarder of PCs 
[above_MP_PCA_num, res_MPPCA] = fxn_Marchenko2PCA(data_bin_z, thrcov_PC_percnet_indicator);
disp([' ', num2str(above_MP_PCA_num), ' PCs are above Marcenko-Pastur threshold.']);
%% Decompose of Correlation matrix
if prms_forced_IC == 0
    i_range = 1:above_MP_PCA_num; disp('MP-PCA mode.');
else
    i_range = 1:prms_forced_IC; disp('Forced IC')
end
cell_ith_decompose = {};
% figure;
for i = 1:i_range(end)
ith_Lambda = pca_eigenvalue(i);
ith_eigenvector = pca_coeff(:,i);
ith_decompose =  ith_Lambda .* ith_eigenvector .* ith_eigenvector';
cell_ith_decompose{i} = ith_decompose;
end
%%
Psign = pca_coeff(:,i_range); % correct
Z_proj = Psign' * data_bin_z';
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
cell_p_k_diag_zero = {};
for i_r = 1:r % pattern_vector_id = 1; % 1:r 
    for i_t = 1:size(data_bin_z,1)% time
    z_t = data_bin_z(i_t,:);
    pattern_vector = V_sqr_scale(:,i_r); % 1 -> r
    p_k = pattern_vector * pattern_vector';
    p_k_diag_zero = p_k;
        for i = 1:size(p_k,1)
        p_k_diag_zero(i,i) = 0;
        end
    cell_p_k_diag_zero{i_r} = p_k_diag_zero;
    r_kt(i_t,i_r) = z_t * p_k_diag_zero * z_t'; %
    end
    for i_t = 1:size(data_bin_z_full,1)% time
    z_t_full = data_bin_z_full(i_t,:);
    r_kt_full(i_t,i_r) = z_t_full * p_k_diag_zero * z_t_full'; %
    end    
%     figure; imagesc(p_k_diag_zero); % for debug;
    disp(['Pattern# ', num2str(i_r),'  out of total ',num2str(r),' Rx(t)s is calculated.']);
end

%% reference data fig, ensemble sorting for check
% for i = 1:1 % for debug
for i = 1:r

figure('Position',[100,0,700,700]); %[left bottom width height] 

% h1 = subaxis(6,1,1, 'SpacingVert',0.06, 'MR',0.15); 
h1 = subplot(611);
stem(V_sqr_scale_z(:,i),'k','MarkerFaceColor','k','MarkerSize', 2); hold on
title(['Reference MPPCA-ICA ensemble #', num2str(i), ' of ' ,num2str(r)])

V3_sqr_scale_z_color = V_sqr_scale_z;
V3_sqr_scale_z_color(negative2) = NaN; 
stem(V3_sqr_scale_z_color(:,i),'r','MarkerFaceColor','r','MarkerSize', 4); 
V3_sqr_scale_z_color_dash = ones(size(data_bin_z,2),1)*prms_reactivation_SD_thr;
plot(V3_sqr_scale_z_color_dash, 'r--');
xlabel('Neurons'); ylabel({'Neurons weight in';'in assemble (SD)'}); xlim([0 size(data_bin_z,2)]); % legend below-thr above-thr SD-thr
hold off


V3_sort_column = V_sqr_scale_z_cutoff(:,i);
V3_sort_id = find(V3_sort_column == 1);
V3_sort_data = data_bin_z(:,V3_sort_id);
fig_time = [1:size(V3_sort_data,1)]/(20/bin_frame_num);
fig_num_cell = [1:size(V3_sort_data,2)];

h2 = subplot(612);

plot(fig_time, r_kt(:,i),'b'); grid on;
ylabel({'Reactivation';'strength (AU)'}); xlim([0 fig_time(end)]); %ylim([0 1])
xticklabels([]); hold on

prms_interest_region1 = prms_interest_region1;  %*(20/bin_frame_num);
interest_shadow1 = (ones(prms_interest_region1(end)-prms_interest_region1(1)+1,1)*max(r_kt(:,i)))';

interest_shadow1map = (ones(prms_interest_region1(end)-prms_interest_region1(1)+1,1)*(size(V3_sort_data,2)+0.5))';
area(prms_interest_region1, interest_shadow1,'FaceColor','y','FaceAlpha',.1,'EdgeAlpha',.1);

prms_interest_region2 = prms_interest_region2;  %*(20/bin_frame_num);
interest_shadow2 = (ones(prms_interest_region2(end)-prms_interest_region2(1)+1,1)*max(r_kt(:,i)))';

interest_shadow2map = (ones(prms_interest_region2(end)-prms_interest_region2(1)+1,1)*(size(V3_sort_data,2)+0.5))';
area(prms_interest_region2, interest_shadow2,'FaceColor','g','FaceAlpha',.1,'EdgeAlpha',.1);

hold off

r_kt_z = zscore(r_kt);
r_kt_z_above_SD = find(r_kt_z >= prms_reactivation_SD_thr);
r_kt_z_below_SD = find(r_kt_z < prms_reactivation_SD_thr);

r_kt_z_color_above_thr = r_kt_z;  
r_kt_z_color_above_thr(r_kt_z_below_SD) = NaN;
r_kt_z_color_above_thr_dash = ones(size(data_bin_z,1),1)*prms_reactivation_SD_thr;

interest_shadow1sd = (ones(prms_interest_region1(end)-prms_interest_region1(1)+1,1)*max(r_kt_z(:,i)))';
interest_shadow2sd = (ones(prms_interest_region2(end)-prms_interest_region2(1)+1,1)*max(r_kt_z(:,i)))';


h3 = subplot(613);
plot(fig_time, r_kt_z(:,i),'k'); hold on; grid on;
area(fig_time, r_kt_z_color_above_thr(:,i),'FaceColor','r','EdgeColor','r'); xticklabels([]); hold on
plot(fig_time, r_kt_z_color_above_thr_dash, 'r--'); 
ylabel({'Reactivation';'strength (SD)'}); xlim([0 fig_time(end)]); 
area(prms_interest_region1, interest_shadow1sd,'FaceColor','y','FaceAlpha',.1,'EdgeAlpha',.1);
area(prms_interest_region2, interest_shadow2sd,'FaceColor','g','FaceAlpha',.1,'EdgeAlpha',.1);
hold off;

h4 = subplot(6,1,4:6);
imagesc(fig_time, fig_num_cell, V3_sort_data'); ylabel({'Pattern-related';'cell activities (z-scored)'}); hold on;
caxis([-3 3]); colormap(fxn_redblue); xlabel('Time (s)'); 
% c= colorbar; c.Label.String = 'z-score'; % for fig
area(prms_interest_region1, interest_shadow1map,'FaceColor','y','FaceAlpha',.15,'EdgeAlpha',.5);
area(prms_interest_region2, interest_shadow2map,'FaceColor','g','FaceAlpha',.15,'EdgeAlpha',.5);

set(findobj(gcf, 'Type', 'Axes'), 'FontSize', 9, 'FontName','Arial');
linkaxes([h2,h3,h4], 'x');
end
%% target data fig, ensemble sorting for check

for i = 1:r

figure('Position',[800,0,700,700]);

h5 = subplot(711);
stem(V_sqr_scale_z(:,i),'k','MarkerFaceColor','k','MarkerSize', 2); hold on
title(['Backprojected MPPCA-ICA ensemble #', num2str(i), ' of ' ,num2str(r)])

V3_sqr_scale_z_color = V_sqr_scale_z;
V3_sqr_scale_z_color(negative2) = NaN; 
stem(V3_sqr_scale_z_color(:,i),'r','MarkerFaceColor','r','MarkerSize', 4); 
V3_sqr_scale_z_color_dash = ones(size(data_bin_z_full,2),1)*prms_reactivation_SD_thr;
plot(V3_sqr_scale_z_color_dash, 'r--');
xlabel('Neurons'); ylabel({'Neurons weight in';'in assemble (SD)'}); xlim([0 size(data_bin_z_full,2)]); % legend below-thr above-thr SD-thr
hold off

V3_sort_column = V_sqr_scale_z_cutoff(:,i);
V3_sort_id = find(V3_sort_column == 1);
V3_sort_id_data{i} =  V3_sort_id;
V3_sort_data = data_bin_z_full(:,V3_sort_id);
fig_time = [1:size(V3_sort_data,1)]/(20/bin_frame_num);
fig_num_cell = [1:size(V3_sort_data,2)];


h6 = subplot(712); 
plot(fig_time, r_kt_full(:,i),'b'); ylabel({'Reactivation';'strength (AU)'}); xlim([0 fig_time(end)]); %ylim([0 1])
xticklabels([]); hold on; grid on;

prms_interest_region1 = prms_interest_region1;  %*(20/bin_frame_num);
interest_shadow1 = (ones(prms_interest_region1(end)-prms_interest_region1(1)+1,1)*max(r_kt_full(:,i)))';

interest_shadow1map = (ones(prms_interest_region1(end)-prms_interest_region1(1)+1,1)*(size(V3_sort_data,2)+0.5))';
area(prms_interest_region1, interest_shadow1,'FaceColor','y','FaceAlpha',.1,'EdgeAlpha',.1);

prms_interest_region2 = prms_interest_region2;  %*(20/bin_frame_num);
interest_shadow2 = (ones(prms_interest_region2(end)-prms_interest_region2(1)+1,1)*max(r_kt_full(:,i)))';

interest_shadow2map = (ones(prms_interest_region2(end)-prms_interest_region2(1)+1,1)*(size(V3_sort_data,2)+0.5))';
area(prms_interest_region2, interest_shadow2,'FaceColor','g','FaceAlpha',.1,'EdgeAlpha',.1);

if prms_ticklabel_mode == 0
    disp('ticklabel mode OFF')
elseif prms_ticklabel_mode == 1
    disp('ticklabel mode ON')
    xticks([]); % update 211104 to off display

end

hold off

r_kt_full_z = zscore(r_kt_full);
r_kt_full_z_above_SD = find(r_kt_full_z >= prms_reactivation_SD_thr);
r_kt_full_z_below_SD = find(r_kt_full_z < prms_reactivation_SD_thr);

r_kt_full_z_color_above_thr = r_kt_full_z;  
r_kt_full_z_color_above_thr(r_kt_full_z_below_SD) = NaN;
r_kt_full_z_color_above_thr_dash = ones(size(data_bin_z_full,1),1)*prms_reactivation_SD_thr;

interest_shadow1sd = (ones(prms_interest_region1(end)-prms_interest_region1(1)+1,1)*max(r_kt_full_z(:,i)))';
interest_shadow2sd = (ones(prms_interest_region2(end)-prms_interest_region2(1)+1,1)*max(r_kt_full_z(:,i)))';

% h3 = subaxis(6,1,3); 
h7 = subplot(7,1,3:4);
plot(fig_time, r_kt_full_z(:,i),'k'); hold on; grid on;
area(fig_time, r_kt_full_z_color_above_thr(:,i),'FaceColor','r','EdgeColor','r'); xticklabels([]); hold on
plot(fig_time, r_kt_full_z_color_above_thr_dash, 'r--'); 
ylabel({'Reactivation';'strength (SD)'}); xlim([0 fig_time(end)]); 
area(prms_interest_region1, interest_shadow1sd,'FaceColor','y','FaceAlpha',.1,'EdgeAlpha',.1);
area(prms_interest_region2, interest_shadow2sd,'FaceColor','g','FaceAlpha',.1,'EdgeAlpha',.1);

if prms_ticklabel_mode == 0
    disp('ticklabel mode OFF')
elseif prms_ticklabel_mode == 1
    disp('ticklabel mode ON')
xticklabels(result_data_cell(1+shift:size(result_data_cell,1),1) );  % for adding ticklabel
xticks([1; cell2mat(result_data_cell(1+shift:size(result_data_cell,1),9))]); % for adding ticklabel
xtickangle(45) % update 211104
end

hold off;

% h4 = subaxis(6,1,4:6);
h8 = subplot(7,1,6:7);
fig_time = [1:size(V3_sort_data,1)]/(20/bin_frame_num);
fig_num_cell = [1:size(V3_sort_data,2)];
imagesc(fig_time, fig_num_cell, V3_sort_data'); ylabel({'Pattern-related';'cell activities (z-scored)'}); hold on;
caxis([-3 3]); colormap(fxn_redblue); xlabel('Time (s)'); 
% c= colorbar; c.Label.String = 'z-score'; % for fig
area(prms_interest_region1, interest_shadow1map,'FaceColor','y','FaceAlpha',.15,'EdgeAlpha',.5);
area(prms_interest_region2, interest_shadow2map,'FaceColor','g','FaceAlpha',.15,'EdgeAlpha',.5);

set(findobj(gcf, 'Type', 'Axes'), 'FontSize', 8, 'FontName','Arial');
linkaxes([h6,h7,h8], 'x');
end
%% update 211104, add r_kt_full-spilitted data 
result_data_cell{1,13} = ('r_kt'); result_data_cell{1,14} = ('r_kt z'); 

for i = 1:size(result_data_cell,1)-1
result_data_cell{i+shift,13} =  r_kt_full(result_data_cell{i+shift,12},:);
result_data_cell{i+shift,14} =  r_kt_full_z(result_data_cell{i+shift,12},:);
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
res_PCAICA.ith_decompose        = ith_decompose;
res_PCAICA.r_strength_ref       = r_kt;
res_PCAICA.r_strength_refz      = r_kt_z;
res_PCAICA.r_strength_target    = r_kt_full;
res_PCAICA.r_strength_targetz   = r_kt_full_z;
res_PCAICA.neuron_weight        = V_sqr_scale;
res_PCAICA.neuron_weight_sigID  = V_sqr_scale_z_cutoff;
res_PCAICA.neuron_sig_IDs       = V3_sort_id_data;
res_PCAICA.result_data_cell     = result_data_cell;
%%
end