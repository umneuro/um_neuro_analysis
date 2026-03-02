function [freeze_data_binned,  freeze_data_binned_mean] = funcHF_temporal_binning_freezingdata(freeze_data_for_binning, binning_num)
%% Temporal binning (1 sec)

%         Ca_data     = HFcal01kikgr75ca1korepfinalv;
%         binning_num = 20;
    
    [n_of_frame, n_of_cell]  = size(freeze_data_for_binning);
    freeze_data_binned = zeros(n_of_frame./binning_num, n_of_cell); % m x n initialize
    
    for i = 1:n_of_cell
     freeze_data_for_binning = reshape(freeze_data_for_binning(:,i), binning_num, []);
     freeze_data_binned(:,i) = mean(freeze_data_for_binning,1);
    end
     
    freeze_data_binned_mean   = mean(freeze_data_binned,2);
    
%%
end