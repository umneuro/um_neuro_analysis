%% data cell preparation
function [result_ca_filt_data_ct, result_data_cell] = fxn_MPPCA_ca_data_process_v4(data_cell, ca_raw_data, prms_MPPCA);
%% for debug
% bin_frame_num = 20;
% sample_fps   = 20; % 20hz

% ### fixed parameter
global shift;
highpass_frq = 0.01; % 30 sec highpass filter
factor       = 1; % factor
cutoff_val   = 0; % cutoff value
%% Decompose parameters
bin_frame_num = prms_MPPCA.bin_frame_num; % 20   -> 1s, 1s binning, For PCA-ICA, I reccomend to not change here, first.
sample_fps    = prms_MPPCA.sample_fps   ; % RGECO: 10FPS, G-CaMP:20hz. Dependent on your data fps.
filter_mode   = prms_MPPCA.filter_mode  ; % 0: disable, 1: enable 0.01hz-highpass butter worth filter.
z_mode   = prms_MPPCA.zscore_mode  ; % 0: disable, 1: enable
cutoff_mode  = prms_MPPCA.cut_off_mode ; % 0: disable, 1: enable (If zscored, cut_off can be active.)
mu_mode       = prms_MPPCA.mu_mode ; % 0: disable, 1: enable
%% Ca data binning
% [ca_bin_time, ca_cell_num, ca_bin_raster] = fxn_mod_round_binning_time(ca_filt_data, bin_frame_num)

data_cell{1,1} = ('session'); data_cell{1,2} = ('frame range'); 
data_cell{1,3} = ('butter res'); data_cell{1,4} = ('Ca_bin_time');
data_cell{1,5} = ('Ca_bin_time'); data_cell{1,6} = ('Ca_cell_num'); data_cell{1,7} = ('Ca_data_binned');
data_cell{1,8} = ('Cumulative frames'); data_cell{1,9} = ('Cumulative time'); 
data_cell{1,10} = ('frame-t onset'); data_cell{1,11} = ('frame-t offset'); data_cell{1,12} = ('frame-t range'); 

total_session_num = size(data_cell(2:end,1),1);
ca_filt_data_ct = [];
cumulative_frames = 0;
% reply = input('Would you like to do filtering? If photo-bleaching obserbed, select yes. (y/n): ','s');

if filter_mode == 1 % strcmp(reply,'y') % 
disp('Butter-filt done, and z-score normalization and mod_rounded_binning done!') 
ca_filt_data =  fxn_butter_filt(ca_raw_data, sample_fps, highpass_frq, factor);
    for i = 1:total_session_num
    data_cell{i+shift,3} = ca_filt_data(data_cell{i+shift,2},:); % butter filt

    if z_mode == 0 && cutoff_mode == 0 && mu_mode == 0
        data_cell{i+shift,4} = data_cell{i+shift,3}; % disable z-score cutoff 
        disp('Simple binning')
    elseif z_mode == 1 && cutoff_mode == 0 && mu_mode == 0
            data_cell{i+shift,4} = zscore(data_cell{i+shift,3}); % disable z-score cutoff 
            disp('Simple Z-socring')
    elseif z_mode ==1 && cutoff_mode == 1 && mu_mode == 0
            data_cell{i+shift,4} = fxn_MPPCA_cutoff(zscore(data_cell{i+shift,3}),cutoff_val); % z-score 
            disp('Z-socring and cutoff')
    elseif z_mode == 0 && cutoff_mode == 1 && mu_mode == 0
            warning('Without zscoring, dF/F value cannot be cut off!!!')
    elseif z_mode == 0 && mu_mode == 1 && cutoff_mode == 0 
            data_cell{i+shift,4} = fxn_NMF_mu_normalization(data_cell{i+shift,3}); 
            disp('Mu-mode ON')
    elseif z_mode == 1 && mu_mode == 1
            warning('Both z-scoreing and mu-normalization cannot be done parallelly!!!')
    end

    [data_cell{i+shift,5}, data_cell{i+shift,6}, data_cell{i+shift,7}] = ...
        fxn_mod_round_binning_time_v2(data_cell{i+shift,4}, bin_frame_num, sample_fps);
    cumulative_frames = cumulative_frames + length(data_cell{i+shift,5});
    data_cell{i+shift,8} =  cumulative_frames;
    data_cell{i+shift,9} =  data_cell{i+shift,8}/(20/bin_frame_num);
    ca_filt_data_ct = cat(1, ca_filt_data_ct, data_cell{i+shift,7});
    end   
    
    % update 211104. give frame info for statistics sorting
    for i = 1:total_session_num
        if i == 1
        data_cell{i+shift,10} = 1;                 data_cell{i+shift,11} = data_cell{i+shift,8};
        data_cell{i+shift,12} = [data_cell{i+shift,10}: data_cell{i+shift,11}];
        else 
        data_cell{i+shift,10} = data_cell{i,8}+1;  data_cell{i+shift,11} = data_cell{i+shift,8};
        data_cell{i+shift,12} = [data_cell{i+shift,10}: data_cell{i+shift,11}];
        end
    end
    
display('Finish filt, z-score normalization, and mod_rounded_binning!')
    
elseif  filter_mode == 0 % strcmp(reply,'n') % 
    disp('Not filtered, then z-score normalization and mod_rounded_binning done!')
    ca_filt_data =  ca_raw_data;
    for i = 1:total_session_num
    data_cell{i+shift,3} = ca_filt_data(data_cell{i+shift,2},:); % butter filt

    if z_mode == 0 && cutoff_mode == 0 && mu_mode == 0
        data_cell{i+shift,4} = data_cell{i+shift,3}; % disable z-score cutoff 
        disp('Simple binning')
    elseif z_mode == 1 && cutoff_mode == 0 && mu_mode == 0
            data_cell{i+shift,4} = zscore(data_cell{i+shift,3}); % disable z-score cutoff 
            disp('Simple Z-socring')
    elseif z_mode ==1 && cutoff_mode == 1 && mu_mode == 0
            data_cell{i+shift,4} = fxn_MPPCA_cutoff(zscore(data_cell{i+shift,3}),cutoff_val); % z-score 
            disp('Z-socring and cutoff')
    elseif z_mode == 0 && cutoff_mode == 1 && mu_mode == 0
            warning('Without zscoring, dF/F value cannot be cut off!!!')
    elseif z_mode == 0 && mu_mode == 1 && cutoff_mode == 0 
            data_cell{i+shift,4} = fxn_NMF_mu_normalization(data_cell{i+shift,3}); 
            disp('Mu-mode ON')
    elseif z_mode == 1 && mu_mode == 1
            warning('Both z-scoreing and mu-normalization cannot be done parallelly!!!')
    end
    
        [data_cell{i+shift,5}, data_cell{i+shift,6}, data_cell{i+shift,7}] = ...
        fxn_mod_round_binning_time_v2(data_cell{i+shift,4}, bin_frame_num, sample_fps);
    cumulative_frames = cumulative_frames + length(data_cell{i+shift,5});
    data_cell{i+shift,8} =  cumulative_frames;
    data_cell{i+shift,9} =  data_cell{i+shift,8}/(20/bin_frame_num);    
    ca_filt_data_ct = cat(1, ca_filt_data_ct, data_cell{i+shift,7});    
    
    end 
    
% update 211104. give frame info for statistics sorting
    for i = 1:total_session_num
        if i == 1
        data_cell{i+shift,10} = 1;                 data_cell{i+shift,11} = data_cell{i+shift,8};
        data_cell{i+shift,12} = [data_cell{i+shift,10}: data_cell{i+shift,11}];
        else 
        data_cell{i+shift,10} = data_cell{i,8}+1;  data_cell{i+shift,11} = data_cell{i+shift,8};
        data_cell{i+shift,12} = [data_cell{i+shift,10}: data_cell{i+shift,11}];
        end
    end
    
disp('   Finish mod_rounded_binning!')

    else
       error('Unexpected situation')
end
%% save data
result_data_cell = data_cell;
result_ca_filt_data_ct = ca_filt_data_ct;
%%

%%
end
