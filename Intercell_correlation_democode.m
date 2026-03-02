
%% First perform PCA/ICA to get the ensemble and neural activity and upload behavioral data timestamp in workspace

load ca1sod6_etho_behavior_timestamps;


%% Create postshock behavior timestamps

% % manual time stamping
CFC_postshock_stamp = zeros(92,1);
fear_stamp = [1:5 61:65];
CFC_postshock_stamp(fear_stamp) = 1;


%% Extract Ca data of SI ICR and CFC freezing

SI_ICR_session_number = 14;
CFC_recall_session_number = 23;
CFC_postshock_session_number = 19;

Sod_day6_homecage_session_number = 7;
Sod_day6_session_number = 8;

SI_ICR_Ca_data = result_data_cell{shift+SI_ICR_session_number,7};
CFC_recall_Ca_data = result_data_cell{shift+CFC_recall_session_number,7};

CFC_postshock_Ca_data = result_data_cell{shift+CFC_postshock_session_number,7};


%% correlation inter cell types during freezing  (using all behavior binarized data tiempoints)

p_threshold = 0.05;
correlation_method = 'kendall'; % or 'pearson', 'kendall'
behavior_cell_method = 'kendall' % 'wilcoxon', 'foldchange', 'kendall', or 'pearson'

results = compute_inter_cell_correlations_freeze_demo(SI_ICR_Ca_data, CFC_recall_Ca_data, CFC_postshock_Ca_data, ...
    mouse_sod_SI_inactive_binary, mouse_sod_SI_opposite_corners_combined, mouse_freeze_binary, CFC_postshock_stamp, ...
    correlation_method, p_threshold, behavior_cell_method)

%% intercell correlation during non-freezing

results = compute_inter_cell_correlations_nonfreeze_demo(SI_ICR_Ca_data, CFC_recall_Ca_data, CFC_postshock_Ca_data, ...
    mouse_sod_SI_inactive_binary, mouse_sod_SI_opposite_corners_combined, mouse_freeze_binary, CFC_postshock_stamp, ...
    correlation_method, p_threshold, behavior_cell_method)

