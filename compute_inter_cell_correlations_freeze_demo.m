function results = compute_inter_cell_correlations_freeze_demo(SI_ICR_Ca_data, CFC_recall_Ca_data, CFC_postshock_Ca_data, ...
    mouse_sod_SI_inactive_binary, mouse_sod_SI_opposite_corners_combined, mouse_freeze_binary, CFC_postshock_stamp, ...
    correlation_method, p_threshold, behavior_cell_method)

%% Input Validation
if ~ismember(correlation_method, {'pearson', 'kendall'})
    error('Invalid correlation type. Use ''pearson'' or ''kendall''.');
end
if p_threshold <= 0 || p_threshold >= 1
    error('Error: p_threshold must be between 0 and 1.');
end
if isempty(SI_ICR_Ca_data) || isempty(CFC_recall_Ca_data) || isempty(CFC_postshock_Ca_data)
    error('Input data matrices cannot be empty.');
end
if nargin < 10 || isempty(behavior_cell_method)
    behavior_cell_method = 'wilcoxon'; % Default
end
if ~ismember(behavior_cell_method, {'wilcoxon', 'foldchange', 'kendall', 'pearson'})
    error('Invalid behavior_cell_method. Use ''wilcoxon'', ''foldchange'', ''kendall'', or ''pearson''.');
end

% Validate dimensions
[num_times_SI, num_neurons] = size(SI_ICR_Ca_data);
[num_times_CFC_recall, num_n] = size(CFC_recall_Ca_data);
[num_times_CFC_postshock, num_m] = size(CFC_postshock_Ca_data);
if num_n ~= num_neurons || num_m ~= num_neurons
    error('Neuron count mismatch: SI (%d), CFC recall (%d), CFC postshock (%d)', ...
        num_neurons, num_n, num_m);
end

% Align behavioral vectors
mouse_sod_SI_inactive_binary = logical(mouse_sod_SI_inactive_binary(1:min(length(mouse_sod_SI_inactive_binary), num_times_SI)));
mouse_sod_SI_inactive_binary(end+1:num_times_SI) = false;

mouse_sod_SI_opposite_corners_combined = logical(mouse_sod_SI_opposite_corners_combined(1:min(length(mouse_sod_SI_opposite_corners_combined), num_times_SI)));
mouse_sod_SI_opposite_corners_combined(end+1:num_times_SI) = false;

mouse_freeze_binary = logical(mouse_freeze_binary(1:min(length(mouse_freeze_binary), num_times_CFC_recall)));
mouse_freeze_binary(end+1:num_times_CFC_recall) = false;

CFC_postshock_stamp = logical(CFC_postshock_stamp(1:min(length(CFC_postshock_stamp), num_times_CFC_postshock)));
CFC_postshock_stamp(end+1:num_times_CFC_postshock) = false;

% Get time indices
SI_inactive_times = find(mouse_sod_SI_inactive_binary);
SI_opposite_times = find(mouse_sod_SI_opposite_corners_combined);
CFC_freeze_times = find(mouse_freeze_binary);
CFC_postshock_times = find(CFC_postshock_stamp);

% Validate time indices
SI_inactive_times = SI_inactive_times(SI_inactive_times <= num_times_SI);
SI_opposite_times = SI_opposite_times(SI_opposite_times <= num_times_SI);
CFC_freeze_times = CFC_freeze_times(CFC_freeze_times <= num_times_CFC_recall);
CFC_postshock_times = CFC_postshock_times(CFC_postshock_times <= num_times_CFC_postshock);

% Debug: Time points
fprintf('SI Inactive: %d time points\n', length(SI_inactive_times));
fprintf('SI Opposite: %d time points\n', length(SI_opposite_times));
fprintf('CFC Freeze: %d time points\n', length(CFC_freeze_times));
fprintf('CFC Postshock: %d time points\n', length(CFC_postshock_times));

if isempty(CFC_freeze_times)
    warning('No CFC Freeze time points found. Returning empty results.');
    results = struct(...
        'SI_inactive_cells', [], ...
        'SI_opposite_cells', [], ...
        'CFC_freeze_cells', [], ...
        'CFC_postshock_cells', [], ...
        'CFC_combined_cells', [], ...
        'CFC_freeze_times', CFC_freeze_times, ...
        'num_sig_corr', [], ...
        'avg_r_corr', [], ...
        'labels', {{}}, ...
        'SI_inactive_vs_CFC_freeze_sig_neurons', [], ...
        'SI_inactive_vs_CFC_postshock_sig_neurons', [], ...
        'SI_inactive_vs_CFC_combined_sig_neurons', [], ...
        'SI_opposite_vs_CFC_freeze_sig_neurons', [], ...
        'SI_opposite_vs_CFC_postshock_sig_neurons', [], ...
        'SI_opposite_vs_CFC_combined_sig_neurons', []);
    return;
end

%% Identify Neurons
SI_inactive_cells = false(num_neurons, 1);
SI_opposite_cells = false(num_neurons, 1);
CFC_freeze_cells = false(num_neurons, 1);
CFC_postshock_cells = false(num_neurons, 1);
active_behaviors = false(num_neurons, 4);
foldchange_threshold = 1.5; % Threshold for foldchange method

for i = 1:num_neurons
    non_behavior_SI = find(~mouse_sod_SI_inactive_binary & ~mouse_sod_SI_opposite_corners_combined);
    non_behavior_CFC_recall = find(~mouse_freeze_binary);
    non_behavior_CFC_postshock = find(~CFC_postshock_stamp);
    
    % SI Inactive
    if length(SI_inactive_times) >= 2 && length(non_behavior_SI) >= 2
        data_on = SI_ICR_Ca_data(SI_inactive_times, i);
        data_off = SI_ICR_Ca_data(non_behavior_SI, i);
        switch behavior_cell_method
            case 'wilcoxon'
                [p, ~] = ranksum(data_on, data_off);
                active_behaviors(i, 1) = p < p_threshold && mean(data_on, 'omitnan') > mean(data_off, 'omitnan');
            case 'foldchange'
                mean_on = mean(data_on, 'omitnan');
                mean_off = mean(data_off, 'omitnan');
                if mean_off == 0
                    mean_off = eps; % Avoid division by zero
                end
                active_behaviors(i, 1) = (mean_on / mean_off) > foldchange_threshold;
            case {'kendall', 'pearson'}
                [r, p] = corr(SI_ICR_Ca_data(:, i), double(mouse_sod_SI_inactive_binary), ...
                              'Type', behavior_cell_method, 'Rows', 'pairwise');
                active_behaviors(i, 1) = p < p_threshold && r > 0;
        end
    end
    
    % SI Opposite
    if length(SI_opposite_times) >= 2 && length(non_behavior_SI) >= 2
        data_on = SI_ICR_Ca_data(SI_opposite_times, i);
        data_off = SI_ICR_Ca_data(non_behavior_SI, i);
        switch behavior_cell_method
            case 'wilcoxon'
                [p, ~] = ranksum(data_on, data_off);
                active_behaviors(i, 2) = p < p_threshold && mean(data_on, 'omitnan') > mean(data_off, 'omitnan');
            case 'foldchange'
                mean_on = mean(data_on, 'omitnan');
                mean_off = mean(data_off, 'omitnan');
                if mean_off == 0
                    mean_off = eps;
                end
                active_behaviors(i, 2) = (mean_on / mean_off) > foldchange_threshold;
            case {'kendall', 'pearson'}
                [r, p] = corr(SI_ICR_Ca_data(:, i), double(mouse_sod_SI_opposite_corners_combined), ...
                              'Type', behavior_cell_method, 'Rows', 'pairwise');
                active_behaviors(i, 2) = p < p_threshold && r > 0;
        end
    end
    
    % CFC Freeze
    if length(CFC_freeze_times) >= 2 && length(non_behavior_CFC_recall) >= 2
        data_on = CFC_recall_Ca_data(CFC_freeze_times, i);
        data_off = CFC_recall_Ca_data(non_behavior_CFC_recall, i);
        switch behavior_cell_method
            case 'wilcoxon'
                [p, ~] = ranksum(data_on, data_off);
                active_behaviors(i, 3) = p < p_threshold && mean(data_on, 'omitnan') > mean(data_off, 'omitnan');
            case 'foldchange'
                mean_on = mean(data_on, 'omitnan');
                mean_off = mean(data_off, 'omitnan');
                if mean_off == 0
                    mean_off = eps;
                end
                active_behaviors(i, 3) = (mean_on / mean_off) > foldchange_threshold;
            case {'kendall', 'pearson'}
                [r, p] = corr(CFC_recall_Ca_data(:, i), double(mouse_freeze_binary), ...
                              'Type', behavior_cell_method, 'Rows', 'pairwise');
                active_behaviors(i, 3) = p < p_threshold && r > 0;
        end
    end
    
    % CFC Postshock
    if length(CFC_postshock_times) >= 2 && length(non_behavior_CFC_postshock) >= 2
        data_on = CFC_postshock_Ca_data(CFC_postshock_times, i);
        data_off = CFC_postshock_Ca_data(non_behavior_CFC_postshock, i);
        switch behavior_cell_method
            case 'wilcoxon'
                [p, ~] = ranksum(data_on, data_off);
                active_behaviors(i, 4) = p < p_threshold && mean(data_on, 'omitnan') > mean(data_off, 'omitnan');
            case 'foldchange'
                mean_on = mean(data_on, 'omitnan');
                mean_off = mean(data_off, 'omitnan');
                if mean_off == 0
                    mean_off = eps;
                end
                active_behaviors(i, 4) = (mean_on / mean_off) > foldchange_threshold;
            case {'kendall', 'pearson'}
                [r, p] = corr(CFC_postshock_Ca_data(:, i), double(CFC_postshock_stamp), ...
                              'Type', behavior_cell_method, 'Rows', 'pairwise');
                active_behaviors(i, 4) = p < p_threshold && r > 0;
        end
    end
    
    % Single-behavior cells
    if sum(active_behaviors(i, :)) == 1
        if active_behaviors(i, 1)
            SI_inactive_cells(i) = true;
        elseif active_behaviors(i, 2)
            SI_opposite_cells(i) = true;
        elseif active_behaviors(i, 3)
            CFC_freeze_cells(i) = true;
        elseif active_behaviors(i, 4)
            CFC_postshock_cells(i) = true;
        end
    end
end

% Combined CFC group
CFC_combined_cells = (CFC_freeze_cells | CFC_postshock_cells) & ...
    ~(SI_inactive_cells | SI_opposite_cells);

% Debug: Cell counts
fprintf('SI Inactive cells: %d\n', sum(SI_inactive_cells));
fprintf('SI Opposite cells: %d\n', sum(SI_opposite_cells));
fprintf('CFC Freeze cells: %d\n', sum(CFC_freeze_cells));
fprintf('CFC Postshock cells: %d\n', sum(CFC_postshock_cells));
fprintf('CFC Combined cells: %d\n', sum(CFC_combined_cells));

% Check for empty groups
if sum(SI_inactive_cells) == 0 && sum(SI_opposite_cells) == 0
    warning('No SI cells found. Returning empty results.');
    results = struct(...
        'SI_inactive_cells', SI_inactive_cells, ...
        'SI_opposite_cells', SI_opposite_cells, ...
        'CFC_freeze_cells', CFC_freeze_cells, ...
        'CFC_postshock_cells', CFC_postshock_cells, ...
        'CFC_combined_cells', CFC_combined_cells, ...
        'CFC_freeze_times', CFC_freeze_times, ...
        'num_sig_corr', [], ...
        'avg_r_corr', [], ...
        'labels', {{}}, ...
        'SI_inactive_vs_CFC_freeze_sig_neurons', [], ...
        'SI_inactive_vs_CFC_postshock_sig_neurons', [], ...
        'SI_inactive_vs_CFC_combined_sig_neurons', [], ...
        'SI_opposite_vs_CFC_freeze_sig_neurons', [], ...
        'SI_opposite_vs_CFC_postshock_sig_neurons', [], ...
        'SI_opposite_vs_CFC_combined_sig_neurons', []);
    return;
end
if sum(CFC_freeze_cells) == 0 && sum(CFC_postshock_cells) == 0
    warning('No CFC cells found. Returning empty results.');
    results = struct(...
        'SI_inactive_cells', SI_inactive_cells, ...
        'SI_opposite_cells', SI_opposite_cells, ...
        'CFC_freeze_cells', CFC_freeze_cells, ...
        'CFC_postshock_cells', CFC_postshock_cells, ...
        'CFC_combined_cells', CFC_combined_cells, ...
        'CFC_freeze_times', CFC_freeze_times, ...
        'num_sig_corr', [], ...
        'avg_r_corr', [], ...
        'labels', {{}}, ...
        'SI_inactive_vs_CFC_freeze_sig_neurons', [], ...
        'SI_inactive_vs_CFC_postshock_sig_neurons', [], ...
        'SI_inactive_vs_CFC_combined_sig_neurons', [], ...
        'SI_opposite_vs_CFC_freeze_sig_neurons', [], ...
        'SI_opposite_vs_CFC_postshock_sig_neurons', [], ...
        'SI_opposite_vs_CFC_combined_sig_neurons', []);
    return;
end

%% Define Pairs
num_pairs = 6;
labels = {...
    'SI Inactive vs CFC Freeze', ...
    'SI Inactive vs CFC Postshock', ...
    'SI Inactive vs CFC Combined', ...
    'SI Opposite vs CFC Freeze', ...
    'SI Opposite vs CFC Postshock', ...
    'SI Opposite vs CFC Combined'};
corr_matrices = cell(num_pairs, 1);
pval_matrices = cell(num_pairs, 1);
num_sig_corr = zeros(num_pairs, 1);
avg_r_corr = zeros(num_pairs, 1);
sig_neuron_pairs = cell(num_pairs, 1); % Store significant neuron pairs

% Get indices
SI_inactive_idx = find(SI_inactive_cells);
SI_opposite_idx = find(SI_opposite_cells);
CFC_freeze_idx = find(CFC_freeze_cells);
CFC_postshock_idx = find(CFC_postshock_cells);
CFC_combined_idx = find(CFC_combined_cells);

% Pair 1: SI Inactive vs CFC Freeze
if ~isempty(SI_inactive_idx) && ~isempty(CFC_freeze_idx)
    n1 = length(SI_inactive_idx);
    n2 = length(CFC_freeze_idx);
    corr_matrices{1} = nan(n1, n2);
    pval_matrices{1} = ones(n1, n2);
    sig_corrs = [];
    sig_pairs = []; % [neuron1_idx, neuron2_idx]
    for i = 1:n1
        for j = 1:n2
            if ~isempty(CFC_freeze_times) && i <= length(SI_inactive_idx) && j <= length(CFC_freeze_idx)
                data1 = CFC_recall_Ca_data(CFC_freeze_times, SI_inactive_idx(i));
                data2 = CFC_recall_Ca_data(CFC_freeze_times, CFC_freeze_idx(j));
                if ~any(isnan(data1)) && ~any(isnan(data2))
                    [r, p] = corr(data1, data2, 'Type', correlation_method);
                    corr_matrices{1}(i,j) = r;
                    pval_matrices{1}(i,j) = p;
                    if p < 0.05 && r > 0
                        sig_corrs = [sig_corrs; r];
                        sig_pairs = [sig_pairs; [SI_inactive_idx(i), CFC_freeze_idx(j)]];
                    end
                end
            end
        end
    end
    num_sig_corr(1) = (length(sig_corrs) / (n1 * n2)) * 100;
    avg_r_corr(1) = mean(sig_corrs, 'omitnan');
    if isnan(avg_r_corr(1))
        avg_r_corr(1) = 0;
    end
    sig_neuron_pairs{1} = sig_pairs;
    fprintf('%s: %.2f%% significant positive correlations, avg r = %.2f\n', ...
        labels{1}, num_sig_corr(1), avg_r_corr(1));
end

% Pair 2: SI Inactive vs CFC Postshock
if ~isempty(SI_inactive_idx) && ~isempty(CFC_postshock_idx)
    n1 = length(SI_inactive_idx);
    n2 = length(CFC_postshock_idx);
    corr_matrices{2} = nan(n1, n2);
    pval_matrices{2} = ones(n1, n2);
    sig_corrs = [];
    sig_pairs = [];
    for i = 1:n1
        for j = 1:n2
            if ~isempty(CFC_freeze_times) && i <= length(SI_inactive_idx) && j <= length(CFC_postshock_idx)
                data1 = CFC_recall_Ca_data(CFC_freeze_times, SI_inactive_idx(i));
                data2 = CFC_recall_Ca_data(CFC_freeze_times, CFC_postshock_idx(j));
                if ~any(isnan(data1)) && ~any(isnan(data2))
                    [r, p] = corr(data1, data2, 'Type', correlation_method);
                    corr_matrices{2}(i,j) = r;
                    pval_matrices{2}(i,j) = p;
                    if p < 0.05 && r > 0
                        sig_corrs = [sig_corrs; r];
                        sig_pairs = [sig_pairs; [SI_inactive_idx(i), CFC_postshock_idx(j)]];
                    end
                end
            end
        end
    end
    num_sig_corr(2) = (length(sig_corrs) / (n1 * n2)) * 100;
    avg_r_corr(2) = mean(sig_corrs, 'omitnan');
    if isnan(avg_r_corr(2))
        avg_r_corr(2) = 0;
    end
    sig_neuron_pairs{2} = sig_pairs;
    fprintf('%s: %.2f%% significant positive correlations, avg r = %.2f\n', ...
        labels{2}, num_sig_corr(2), avg_r_corr(2));
end

% Pair 3: SI Inactive vs CFC Combined
if ~isempty(SI_inactive_idx) && ~isempty(CFC_combined_idx)
    n1 = length(SI_inactive_idx);
    n2 = length(CFC_combined_idx);
    corr_matrices{3} = nan(n1, n2);
    pval_matrices{3} = ones(n1, n2);
    sig_corrs = [];
    sig_pairs = [];
    for i = 1:n1
        for j = 1:n2
            if ~isempty(CFC_freeze_times) && i <= length(SI_inactive_idx) && j <= length(CFC_combined_idx)
                data1 = CFC_recall_Ca_data(CFC_freeze_times, SI_inactive_idx(i));
                data2 = CFC_recall_Ca_data(CFC_freeze_times, CFC_combined_idx(j));
                if ~any(isnan(data1)) && ~any(isnan(data2))
                    [r, p] = corr(data1, data2, 'Type', correlation_method);
                    corr_matrices{3}(i,j) = r;
                    pval_matrices{3}(i,j) = p;
                    if p < 0.05 && r > 0
                        sig_corrs = [sig_corrs; r];
                        sig_pairs = [sig_pairs; [SI_inactive_idx(i), CFC_combined_idx(j)]];
                    end
                end
            end
        end
    end
    num_sig_corr(3) = (length(sig_corrs) / (n1 * n2)) * 100;
    avg_r_corr(3) = mean(sig_corrs, 'omitnan');
    if isnan(avg_r_corr(3))
        avg_r_corr(3) = 0;
    end
    sig_neuron_pairs{3} = sig_pairs;
    fprintf('%s: %.2f%% significant positive correlations, avg r = %.2f\n', ...
        labels{3}, num_sig_corr(3), avg_r_corr(3));
end

% Pair 4: SI Opposite vs CFC Freeze
if ~isempty(SI_opposite_idx) && ~isempty(CFC_freeze_idx)
    n1 = length(SI_opposite_idx);
    n2 = length(CFC_freeze_idx);
    corr_matrices{4} = nan(n1, n2);
    pval_matrices{4} = ones(n1, n2);
    sig_corrs = [];
    sig_pairs = [];
    for i = 1:n1
        for j = 1:n2
            if ~isempty(CFC_freeze_times) && i <= length(SI_opposite_idx) && j <= length(CFC_freeze_idx)
                data1 = CFC_recall_Ca_data(CFC_freeze_times, SI_opposite_idx(i));
                data2 = CFC_recall_Ca_data(CFC_freeze_times, CFC_freeze_idx(j));
                if ~any(isnan(data1)) && ~any(isnan(data2))
                    [r, p] = corr(data1, data2, 'Type', correlation_method);
                    corr_matrices{4}(i,j) = r;
                    pval_matrices{4}(i,j) = p;
                    if p < 0.05 && r > 0
                        sig_corrs = [sig_corrs; r];
                        sig_pairs = [sig_pairs; [SI_opposite_idx(i), CFC_freeze_idx(j)]];
                    end
                end
            end
        end
    end
    num_sig_corr(4) = (length(sig_corrs) / (n1 * n2)) * 100;
    avg_r_corr(4) = mean(sig_corrs, 'omitnan');
    if isnan(avg_r_corr(4))
        avg_r_corr(4) = 0;
    end
    sig_neuron_pairs{4} = sig_pairs;
    fprintf('%s: %.2f%% significant positive correlations, avg r = %.2f\n', ...
        labels{4}, num_sig_corr(4), avg_r_corr(4));
end

% Pair 5: SI Opposite vs CFC Postshock
if ~isempty(SI_opposite_idx) && ~isempty(CFC_postshock_idx)
    n1 = length(SI_opposite_idx);
    n2 = length(CFC_postshock_idx);
    corr_matrices{5} = nan(n1, n2);
    pval_matrices{5} = ones(n1, n2);
    sig_corrs = [];
    sig_pairs = [];
    for i = 1:n1
        for j = 1:n2
            if ~isempty(CFC_freeze_times) && i <= length(SI_opposite_idx) && j <= length(CFC_postshock_idx)
                data1 = CFC_recall_Ca_data(CFC_freeze_times, SI_opposite_idx(i));
                data2 = CFC_recall_Ca_data(CFC_freeze_times, CFC_postshock_idx(j));
                if ~any(isnan(data1)) && ~any(isnan(data2))
                    [r, p] = corr(data1, data2, 'Type', correlation_method);
                    corr_matrices{5}(i,j) = r;
                    pval_matrices{5}(i,j) = p;
                    if p < 0.05 && r > 0
                        sig_corrs = [sig_corrs; r];
                        sig_pairs = [sig_pairs; [SI_opposite_idx(i), CFC_postshock_idx(j)]];
                    end
                end
            end
        end
    end
    num_sig_corr(5) = (length(sig_corrs) / (n1 * n2)) * 100;
    avg_r_corr(5) = mean(sig_corrs, 'omitnan');
    if isnan(avg_r_corr(5))
        avg_r_corr(5) = 0;
    end
    sig_neuron_pairs{5} = sig_pairs;
    fprintf('%s: %.2f%% significant positive correlations, avg r = %.2f\n', ...
        labels{5}, num_sig_corr(5), avg_r_corr(5));
end

% Pair 6: SI Opposite vs CFC Combined
if ~isempty(SI_opposite_idx) && ~isempty(CFC_combined_idx)
    n1 = length(SI_opposite_idx);
    n2 = length(CFC_combined_idx);
    corr_matrices{6} = nan(n1, n2);
    pval_matrices{6} = ones(n1, n2);
    sig_corrs = [];
    sig_pairs = [];
    for i = 1:n1
        for j = 1:n2
            if ~isempty(CFC_freeze_times) && i <= length(SI_opposite_idx) && j <= length(CFC_combined_idx)
                data1 = CFC_recall_Ca_data(CFC_freeze_times, SI_opposite_idx(i));
                data2 = CFC_recall_Ca_data(CFC_freeze_times, CFC_combined_idx(j));
                if ~any(isnan(data1)) && ~any(isnan(data2))
                    [r, p] = corr(data1, data2, 'Type', correlation_method);
                    corr_matrices{6}(i,j) = r;
                    pval_matrices{6}(i,j) = p;
                    if p < 0.05 && r > 0
                        sig_corrs = [sig_corrs; r];
                        sig_pairs = [sig_pairs; [SI_opposite_idx(i), CFC_combined_idx(j)]];
                    end
                end
            end
        end
    end
    num_sig_corr(6) = (length(sig_corrs) / (n1 * n2)) * 100;
    avg_r_corr(6) = mean(sig_corrs, 'omitnan');
    if isnan(avg_r_corr(6))
        avg_r_corr(6) = 0;
    end
    sig_neuron_pairs{6} = sig_pairs;
    fprintf('%s: %.2f%% significant positive correlations, avg r = %.2f\n', ...
        labels{6}, num_sig_corr(6), avg_r_corr(6));
end

% Store results
results = struct(...
    'SI_inactive_cells', SI_inactive_cells, ...
    'SI_opposite_cells', SI_opposite_cells, ...
    'CFC_freeze_cells', CFC_freeze_cells, ...
    'CFC_postshock_cells', CFC_postshock_cells, ...
    'CFC_combined_cells', CFC_combined_cells, ...
    'CFC_freeze_times', CFC_freeze_times, ...
    'num_sig_corr', num_sig_corr, ...
    'avg_r_corr', avg_r_corr, ...
    'labels', {labels}, ...
    'corr_SI_inactive_CFC_freeze', corr_matrices{1}, ...
    'pval_SI_inactive_CFC_freeze', pval_matrices{1}, ...
    'SI_inactive_vs_CFC_freeze_sig_neurons', sig_neuron_pairs{1}, ...
    'corr_SI_inactive_CFC_postshock', corr_matrices{2}, ...
    'pval_SI_inactive_CFC_postshock', pval_matrices{2}, ...
    'SI_inactive_vs_CFC_postshock_sig_neurons', sig_neuron_pairs{2}, ...
    'corr_SI_inactive_CFC_combined', corr_matrices{3}, ...
    'pval_SI_inactive_CFC_combined', pval_matrices{3}, ...
    'SI_inactive_vs_CFC_combined_sig_neurons', sig_neuron_pairs{3}, ...
    'corr_SI_opposite_CFC_freeze', corr_matrices{4}, ...
    'pval_SI_opposite_CFC_freeze', pval_matrices{4}, ...
    'SI_opposite_vs_CFC_freeze_sig_neurons', sig_neuron_pairs{4}, ...
    'corr_SI_opposite_CFC_postshock', corr_matrices{5}, ...
    'pval_SI_opposite_CFC_postshock', pval_matrices{5}, ...
    'SI_opposite_vs_CFC_postshock_sig_neurons', sig_neuron_pairs{5}, ...
    'corr_SI_opposite_CFC_combined', corr_matrices{6}, ...
    'pval_SI_opposite_CFC_combined', pval_matrices{6}, ...
    'SI_opposite_vs_CFC_combined_sig_neurons', sig_neuron_pairs{6});

%% Plotting
% Get screen size for centering figures
screen_size = get(0, 'ScreenSize');
screen_width = screen_size(3);
screen_height = screen_size(4);
vertical_offset = 50; % Pixels to stack figures upward

% Original Correlation Heatmaps
fig_width = 1200;
fig_height = 600;
x_pos = (screen_width - fig_width) / 2;
y_pos = (screen_height - fig_height) / 2;
figure('Name', 'Correlation Heatmaps', 'Position', [x_pos, y_pos, fig_width, fig_height]);
for p = 1:num_pairs
    subplot(2, 3, p);
    if ~isempty(corr_matrices{p})
        imagesc(corr_matrices{p}); colorbar;
        sig_pos = sum((pval_matrices{p} < 0.05) & (corr_matrices{p} > 0), 'all');
        total_pairs = size(corr_matrices{p}, 1) * size(corr_matrices{p}, 2);
        percent_sig_pos = (sig_pos / total_pairs) * 100;
        pos_corrs = corr_matrices{p}(pval_matrices{p} < 0.05 & corr_matrices{p} > 0);
        avg_pos = mean(pos_corrs, 'omitnan');
        if isnan(avg_pos)
            avg_pos = 0;
        end
        title(sprintf('%s\nSig: %d (%.2f%%), Avg r: %.2f', labels{p}, sig_pos, percent_sig_pos, avg_pos), 'FontSize', 10);
        colormap(jet); caxis([-1 1]);
        xlabel(sprintf('Group 2 (n=%d)', size(corr_matrices{p}, 2)));
        ylabel(sprintf('Group 1 (n=%d)', size(corr_matrices{p}, 1)));
    else
        text(0.5, 0.5, ['No Valid ', labels{p}], 'HorizontalAlignment', 'center');
        axis off;
    end
    set(gca, 'FontSize', 8);
end
sgtitle('Correlation Heatmaps During CFC Freeze', 'FontSize', 12);

% Significant Positive Correlations Bar Plot
fig_width = 800;
fig_height = 400;
x_pos = (screen_width - fig_width) / 2;
y_pos = (screen_height - fig_height) / 2 + vertical_offset; % Stack above first figure
figure('Name', 'Significant Positive Correlations', 'Position', [x_pos, y_pos, fig_width, fig_height]);
bar(num_sig_corr);
set(gca, 'XTick', 1:num_pairs, 'XTickLabel', labels, 'XTickLabelRotation', 45, 'FontSize', 8);
ylabel('Percentage of Significant Positive Correlations (%)', 'FontSize', 10);
title('Significant Positive Correlations (p < 0.05, r > 0)', 'FontSize', 12);

% Average r Bar Plot
x_pos = (screen_width - fig_width) / 2;
y_pos = (screen_height - fig_height) / 2 + 2 * vertical_offset; % Stack above second figure
figure('Name', 'Average r of Significant Positive Correlations', 'Position', [x_pos, y_pos, fig_width, fig_height]);
bar(avg_r_corr);
set(gca, 'XTick', 1:num_pairs, 'XTickLabel', labels, 'XTickLabelRotation', 45, 'FontSize', 8);
ylabel('Average r (Significant Positive Correlations)', 'FontSize', 10);
title('Average r of Significant Positive Correlations', 'FontSize', 12);

% Significance Binary Heatmaps with Different Colormaps
colormap_configs = {
    struct('name', 'Parula-Based', 'map', [0 0.5 0; 1 0.8 0], 'desc', 'Non-significant: Emerald Green, Significant: Dark Yellow'),
    struct('name', 'Jet-Based', 'map', [0 0 0.5; 1 0 0], 'desc', 'Non-significant: Dark Blue, Significant: Bright Red'),
    struct('name', 'Hot-Based', 'map', [0.5 0 0; 1 1 0], 'desc', 'Non-significant: Dark Red, Significant: Bright Yellow')
};

fig_width = 1200;
fig_height = 600;1

for c = 1:length(colormap_configs)
    config = colormap_configs{c};
    x_pos = (screen_width - fig_width) / 2;
    y_pos = (screen_height - fig_height) / 2 + (2 + c) * vertical_offset; % Stack above previous figures
    figure('Name', sprintf('Binary Significance Heatmaps - %s', config.name), 'Position', [x_pos, y_pos, fig_width, fig_height]);
    for p = 1:num_pairs
        subplot(2, 3, p);
        if ~isempty(corr_matrices{p})
            % Create binary significance matrix (1 for significant, 0 for non-significant)
            sig_mat = (pval_matrices{p} < 0.05) & (corr_matrices{p} > 0);
            imagesc(sig_mat); colorbar;
            sig_pos = sum(sig_mat, 'all');
            total_pairs = size(corr_matrices{p}, 1) * size(corr_matrices{p}, 2);
            percent_sig_pos = (sig_pos / total_pairs) * 100;
            title(sprintf('%s\nSig: %d (%.2f%%)\n%s', labels{p}, sig_pos, percent_sig_pos, config.desc), 'FontSize', 10);
            colormap(config.map); caxis([0 1]);
            xlabel(sprintf('Group 2 (n=%d)', size(corr_matrices{p}, 2)));
            ylabel(sprintf('Group 1 (n=%d)', size(corr_matrices{p}, 1)));
        else
            text(0.5, 0.5, ['No Valid ', labels{p}], 'HorizontalAlignment', 'center');
            axis off;
        end
        set(gca, 'FontSize', 8);
    end
    sgtitle(sprintf('Binary Significance Heatmaps (%s) During CFC Freeze', config.name), 'FontSize', 12);
end

end