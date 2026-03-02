function [result, coactivity_result] = fxn_MPPCA_ensemble_connectivity_double_multi_session2(result_data_cell, result_MPPCA1, result_MPPCA2, prms)

%% decompose prms (keeping original structure)
session_num_assess_array = prms.session_num_assess; % Now an array of session numbers
binarize_thr        = prms.binarize_thr;
shuffle_mode        = prms.shuffle_mode;
ensemble_num1       = prms.ensemble_num1;
ensemble_num2       = prms.ensemble_num2;
thr_percentile     = prms.thr_percentile;

rng(1,"twister");

shift = 1;  % Added shift variable

%% Initialize storage for all results
all_results = cell(length(session_num_assess_array), 1);
all_coactivity_results = cell(length(session_num_assess_array) + 1, 5); % +1 for header row, changed to 5 columns

% Set header row
all_coactivity_results{1,1} = 'Name of session';  % Modified header
all_coactivity_results{1,2} = 'session number';
all_coactivity_results{1,3} = 'Total count of coactivity';
all_coactivity_results{1,4} = 'Percentage of coactivity';
all_coactivity_results{1,5} = 'Number of coactivity per SoD ensemble';  % Changed index from 4 to 5

%% Loop through each session number
for session_idx = 1:length(session_num_assess_array)
    current_session_num = session_num_assess_array(session_idx);
    
    %% Main processing
    tic;
    result_connect = {};
    parfor i = 1:ensemble_num1
        for ii = 1:ensemble_num2
            [result_connect_cell{i,ii}] = fxn_MPPCA_ensemble_connectivity2_double ...
                (result_data_cell, result_MPPCA1, result_MPPCA2, current_session_num, i, ii, ...
                shuffle_mode, binarize_thr, thr_percentile);
        end
    end
    toc;
    
    %% Results section
    for i = 1:ensemble_num1
        for ii = 1:ensemble_num2
            index_shuffle{i, ii}    = result_connect_cell{i, ii}.shuffle_index;
            index_connect{i, ii}    = result_connect_cell{i, ii}.connect_index;
            vennIDs_spec_A{i, ii}   = result_connect_cell{i, ii}.vennIDs_spec_A;
            vennIDs_spec_B{i, ii}   = result_connect_cell{i, ii}.vennIDs_spec_B;
            vennIDs_both_AB{i, ii}  = result_connect_cell{i, ii}.vennIDs_both_AB;
            index_thr_bottom{i, ii} = result_connect_cell{i, ii}.thr_bottom;
            index_thr_top{i, ii}    = result_connect_cell{i, ii}.thr_top;
        end
    end
    
    %% First figure (Venn diagrams)
    figure('Position', [100 100 700 500]);
    for i = 1:ensemble_num1
        for ii = 1:ensemble_num2
            if ensemble_num1 == ensemble_num2        
                subplot(max(ensemble_num1), max(ensemble_num2), (max(ensemble_num2).*(i-1) + ii))
            elseif ensemble_num1 > ensemble_num2
                subplot(max(ensemble_num1), max(ensemble_num2), (max(ensemble_num2).*(i-1) + ii))
            elseif ensemble_num1 < ensemble_num2
                subplot(max(ensemble_num1), max(ensemble_num2), (max(ensemble_num1).*(i-1) + ii))
            end
            fxn_venn([numel(vennIDs_spec_A{i, ii}), numel(vennIDs_spec_B{i, ii}), numel(vennIDs_both_AB{i, ii})])
            xticklabels([]); yticklabels([]); 
            title(['Session ', num2str(current_session_num), ': Pattern #', num2str(i), ' vs. #', num2str(ii)])
        end
    end
    
    %% Second figure (Histograms)
    figure('Position', [50 50 800 700]);
    for i = 1:ensemble_num1
        for ii = 1:ensemble_num2
            if ensemble_num1 == ensemble_num2        
                subplot(max(ensemble_num1), max(ensemble_num2), (max(ensemble_num2).*(i-1) + ii))
            elseif ensemble_num1 > ensemble_num2
                subplot(max(ensemble_num1), max(ensemble_num2), (max(ensemble_num2).*(i-1) + ii))
            elseif ensemble_num1 < ensemble_num2
                subplot(max(ensemble_num1), max(ensemble_num2), (max(ensemble_num1).*(i-1) + ii))
            end
            histogram(cell2mat(index_shuffle{i, ii}))
            hold on
            if ~isnan(index_connect{i,ii})
                xline((index_connect{i,ii}),'r','Actual value','LineWidth',2)
                xline((index_thr_bottom{i, ii}),'b--','LineWidth',1)
                xline((index_thr_top{i, ii}),'b--','LineWidth',1)
            end
            title(['Session ', num2str(current_session_num), ': Pattern #', num2str(i), ' vs. #', num2str(ii)])
        end
    end
    
    %% Statistical analysis
    index_real   = cell2mat(index_connect);
    index_bottom = cell2mat(index_thr_bottom);
    index_top    = cell2mat(index_thr_top);

    index_res_bottom = index_real <= index_bottom;
    index_res_top    = index_real >= index_top;

    index_res_matrix = (double(index_res_bottom).* -1) + double(index_res_top);
    figure; imagesc(index_res_matrix); colormap(fxn_redblue);
    title(['Session ', num2str(current_session_num), ': Distribution thresholding']);
    xlabel('Pattern #'); ylabel('Pattern #');
    grid on; clim([-1 1]);

    total_count_of_significant_coactivity = sum(index_res_matrix(:) == 1);
    total_size_of_coactivity_matrix = numel(index_res_matrix);
    significant_coactivity_percent = (total_count_of_significant_coactivity/total_size_of_coactivity_matrix)*100;

    number_of_coactivity_per_SoDensemble = sum(index_res_matrix>0,2);
    % number_of_coactivity_per_SoDensemble = sum(index_res_matrix,2); % -1 included

    %% Store results for this session
    current_result.result_connect_cell = result_connect_cell;
    current_result.index_shuffle = index_shuffle;
    current_result.index_connect = index_connect;
    current_result.vennIDs_spec_A = vennIDs_spec_A;
    current_result.vennIDs_spec_B = vennIDs_spec_B;
    current_result.vennIDs_both_AB = vennIDs_both_AB;
    current_result.index_res_matrix = index_res_matrix;
    current_result.total_count_of_significant_coactivity = total_count_of_significant_coactivity;
    current_result.significant_coactivity_percent = significant_coactivity_percent;
    current_result.number_of_coactivity_per_SoDensemble = number_of_coactivity_per_SoDensemble;
    
    all_results{session_idx} = current_result;
    
    %% Store coactivity results for this session
    all_coactivity_results{session_idx + 1, 1} = result_data_cell(current_session_num + shift, 1);  % Modified to include session name with shift
    all_coactivity_results{session_idx + 1, 2} = current_session_num;
    all_coactivity_results{session_idx + 1, 3} = total_count_of_significant_coactivity;
    all_coactivity_results{session_idx + 1, 4} = significant_coactivity_percent;
    all_coactivity_results{session_idx + 1, 5} = number_of_coactivity_per_SoDensemble;  % Changed index from 4 to 5
end

%% Output final results
result = all_results;  % Cell array containing all individual results
coactivity_result = all_coactivity_results;  % Combined table with header and all sessions

end