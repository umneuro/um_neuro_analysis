%% data cell preparation
function [result_ca_filt_data_ct, result_data_cell] = fxn_MPPCA_ca_data_process(data_cell, ca_raw_data, bin_frame_num);
%% for debug
% bin_frame_num = 20;
%% Parameters
sample_fps   = 20; % 20hz
highpass_frq = 0.01; % 30 sec highpass filter
factor       = 1; % factor
cutoff_val   = 0; % cutoff value
%% Ca data binning
% [ca_bin_time, ca_cell_num, ca_bin_raster] = fxn_mod_round_binning_time(ca_filt_data, bin_frame_num)
shift = 1;
data_cell{1,1} = ('session'); data_cell{1,2} = ('frame range'); 
data_cell{1,3} = ('butter res'); data_cell{1,4} = ('Ca_bin_time');
data_cell{1,5} = ('Ca_bin_time'); data_cell{1,6} = ('Ca_cell_num'); data_cell{1,7} = ('Ca_data_binned');
data_cell{1,8} = ('Cumulative frames'); data_cell{1,9} = ('Cumulative time');

total_session_num = size(data_cell(2:end,1),1);
ca_filt_data_ct = [];
cumulative_frames = 0;
reply = input('Would you like to do filtering? If photo-bleaching obserbed, select yes. (y/n): ','s');

if strcmp(reply,'y') % 
disp('Butter-filt done, and z-score normalization and mod_rounded_binning done!') 
ca_filt_data =  fxn_butter_filt(ca_raw_data, sample_fps, highpass_frq, factor);
    for i = 1:total_session_num
    data_cell{i+shift,3} = ca_filt_data(data_cell{i+shift,2},:); % butter filt
    data_cell{i+shift,4} = fxn_MPPCA_cutoff(zscore(data_cell{i+shift,3}),cutoff_val); % z-score 
    [data_cell{i+shift,5}, data_cell{i+shift,6}, data_cell{i+shift,7}] = ...
        fxn_mod_round_binning_time(data_cell{i+shift,4}, bin_frame_num);
    cumulative_frames = cumulative_frames + length(data_cell{i+shift,5});
    data_cell{i+shift,8} =  cumulative_frames;
    data_cell{i+shift,9} =  data_cell{i+shift,8}/(20/bin_frame_num);
    ca_filt_data_ct = cat(1, ca_filt_data_ct, data_cell{i+shift,7});
    end   
display('Finish filt, z-score normalization, and mod_rounded_binning!')
    
    elseif  strcmp(reply,'n') % 
    disp('Not filtered, then z-score normalization and mod_rounded_binning done!')
    ca_filt_data =  ca_raw_data;
    for i = 1:total_session_num
    data_cell{i+shift,3} = ca_filt_data(data_cell{i+shift,2},:); % butter filt
    data_cell{i+shift,4} = fxn_MPPCA_cutoff(zscore(data_cell{i+shift,3}),cutoff_val); % z-score 
    [data_cell{i+shift,5}, data_cell{i+shift,6}, data_cell{i+shift,7}] = ...
        fxn_mod_round_binning_time(data_cell{i+shift,4}, bin_frame_num);
    cumulative_frames = cumulative_frames + length(data_cell{i+shift,5});
    data_cell{i+shift,8} =  cumulative_frames;
    data_cell{i+shift,9} =  data_cell{i+shift,8}/(20/bin_frame_num);    
    ca_filt_data_ct = cat(1, ca_filt_data_ct, data_cell{i+shift,7});
    end 
disp('   Finish mod_rounded_binning!')

    else
       error('Unexpected situation')
end
%% save data
result_data_cell = data_cell;
result_ca_filt_data_ct = ca_filt_data_ct;
%%
end
