%% MPPCA cell-cell pairwise interation 
function [result_multi, result_final_cell] = fxn_behaviorcells_coreactivation_cal_connectivity_multi_core(result_sorted_cell, prms_multi, target_session_num)

%% system parameter
shift = 1; % system parameter, don't change.
%% decompose data and parameter
binarize_thr          = prms_multi.binarize_thr; % Calcium trace binalizing threshold.
shuffle_mode          = prms_multi.shuffle_mode; % 1:randperm, 2:circshift

thr_percentile        = prms_multi.thr_percentile; % def: 1%, top-bottom %-thresholding

ensemble_overlap_mode = prms_multi.ensemble_overlap_mode;
iteration             = prms_multi.iteration; 

ids_A = prms_multi.ids_A;
ids_B = prms_multi.ids_B;
ids_C = prms_multi.ids_C;
ids_D = prms_multi.ids_D;
%% cell sorting
temp_ca_data = result_sorted_cell{target_session_num+shift,7};

%% find unique cell ids
P1_unique = setdiff(ids_A, [ids_B; ids_C; ids_D]);
P2_unique = setdiff(ids_B, [ids_A; ids_C; ids_D]);
P3_unique = setdiff(ids_C, [ids_A; ids_B; ids_D]);
P4_unique = setdiff(ids_D, [ids_A; ids_B; ids_C]);
%% Calculate connectivity
if ensemble_overlap_mode == 0
[Cells_sum] = fxn_MPPCA_cal_connect_4patterns(temp_ca_data, P1_unique, P2_unique, P3_unique, P4_unique, binarize_thr); % exclude overlap
elseif ensemble_overlap_mode == 1
[Cells_sum] = fxn_MPPCA_cal_connect_4patterns(temp_ca_data, ids_A, ids_B, ids_C, ids_D, binarize_thr); % include overlap
end
%% Calculate shuffled connectivity

temp_ca_shuffled{1,1} = fxn_ensemble_shuffling(temp_ca_data, iteration, shuffle_mode); % ensemble shuffling
temp_ca_shuffled{2,1} = fxn_ensemble_shuffling(temp_ca_data, iteration, shuffle_mode); % ensemble shuffling
temp_ca_shuffled{3,1} = fxn_ensemble_shuffling(temp_ca_data, iteration, shuffle_mode); % ensemble shuffling
temp_ca_shuffled{4,1} = fxn_ensemble_shuffling(temp_ca_data, iteration, shuffle_mode); % ensemble shuffling

%% shuffle data cal, This section takes time. 
Cells_sum_shuffle = cell(iteration,1); % for speed up cal.

for i_shuffle = 1:iteration
    if ensemble_overlap_mode == 0
[Cells_sum_shuffle{i_shuffle,1}] = fxn_MPPCA_cal_connect_4patterns_shuffle ...
    (temp_ca_shuffled, P1_unique, P2_unique, P3_unique, P4_unique, binarize_thr, i_shuffle);
elseif ensemble_overlap_mode == 1
[Cells_sum_shuffle{i_shuffle,1}] = fxn_MPPCA_cal_connect_4patterns_shuffle ...
    (temp_ca_shuffled, ids_A, ids_B, ids_C, ids_D, binarize_thr, i_shuffle);
    end

disp([' ## Cal. surrogate data. ', num2str(i_shuffle),'/',num2str(iteration), ' iterations. in session# ', num2str(target_session_num), '.' ]);
end

%% Distribution thresholding

[Cells_sum_shuffle_res] = fxn_MPPCA_cal_distribution_thr_realine(Cells_sum_shuffle, iteration, thr_percentile);

%% output results
result_multi.Cells_sum             = Cells_sum;
result_multi.Cells_sum_shuffle_res = Cells_sum_shuffle_res;
result_multi.P1_unique = P1_unique;
result_multi.P2_unique = P2_unique;
result_multi.P3_unique = P3_unique;
result_multi.P4_unique = P4_unique;
result_multi.ids_A = ids_A;
result_multi.ids_B = ids_B;
result_multi.ids_C = ids_C;
result_multi.ids_D = ids_D;

%%
disp(' ## Cal. fin ##')

%% cal ABCD index
result_cell_4D = {};
result_cell_4D_surrogate = {};

for i_res = 1:11
                    result_cell_4D_temp = result_multi.Cells_sum{i_res, 4};
                    result_cell_4D_temp(isnan(result_cell_4D_temp)) = 0;
                    result_cell_4D{i_res,1} = result_cell_4D_temp;
                    
                    result_cell_4D_surrogate_temp = ... 
                       mean(result_multi.Cells_sum_shuffle_res.Cells_sum_shuffle_index{i_res, 1});
                    result_cell_4D_surrogate_temp(isnan(result_cell_4D_surrogate_temp)) = 0;
                    result_cell_4D_surrogate{i_res,1} = result_cell_4D_surrogate_temp;

result_final_cell{i_res,2} = mean(cell2mat(result_cell_4D(i_res)),[1 2 3 4]);
result_final_cell{i_res,3} = mean(cell2mat(result_cell_4D_surrogate(i_res)),[1 2 3 4]);
end
%% Output results
result_final_cell(1:11,1)= {'AB','AC','AD','BC','BD','CD','ABC','ABD','ACD','BCD','ABCD'};
disp('finish calculaiton for multi');

ax1 = figure('Position',[100 100 1200 300]);
subplot(131)
plot(cell2mat(result_final_cell(:,2)),'k'); hold on 
plot(cell2mat(result_final_cell(:,3)),'--r'); hold off
xticks(1:11); xticklabels(result_final_cell(1:11,1))
title('Multiple synchronized activity')
ylabel('Synchronized index')
xlabel('Combination of patterns')
legend('Real', 'Surrogate'); box off
set(findobj(gcf, 'Type', 'Axes'), 'FontSize', 8, 'FontName','Arial');

subplot(132)
plot(cell2mat(result_final_cell(1:6,2)),'k'); hold on 
plot(cell2mat(result_final_cell(1:6,3)),'--r'); hold off
xticks(1:6); xticklabels(result_final_cell(1:6,1))
title('Multiple synchronized activity')
ylabel('Synchronized index')
xlabel('Combination of patterns')
legend('Real', 'Surrogate'); box off
set(findobj(gcf, 'Type', 'Axes'), 'FontSize', 8, 'FontName','Arial');

subplot(133)
plot(cell2mat(result_final_cell(7:11,2)),'k'); hold on 
plot(cell2mat(result_final_cell(7:11,3)),'--r'); hold off 
xticks(1:5); xticklabels(result_final_cell(7:11,1))
title('Multiple synchronized activity')
ylabel('Synchronized index')
xlabel('Combination of patterns')
legend('Real', 'Surrogate'); box off
set(findobj(gcf, 'Type', 'Axes'), 'FontSize', 8, 'FontName','Arial');

%% (Multi mode) Parameter input section
% ### input saved filename. To protect calculation resutls.
dt = datetime('now');
filename_DateString = datestr(dt,'yymmdd_HHMM');
filename_for_save = (['zRes_',filename_DateString,'session',num2str(target_session_num),'.mat']); 

% for mat file
save(filename_for_save, 'result_final_cell', 'result_multi', '-mat');

% for pdf file
extension = ('pdf');
file_name2 = [filename_for_save, 'Figure_.', extension];
exportgraphics(ax1, file_name2, 'ContentType','vector');
%%
end