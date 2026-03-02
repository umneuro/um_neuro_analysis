%% MPPCA-ICA calculation code
clc; clear; close all; disp('Previous data is cleared'); tic;
data_cell = {}; global shift; shift=1; % Dont modify this line
% addpath("fxn_box_homemade\"); 
load('ca_data_mousesod6.mat') % Ca data

%% ### load ca data ###
ca_raw_data      = mousesod6 ; % input time x neuron matrix

%% %%  mouse sod 6 session frame information

data_cell{1+shift,1} = ('Control day 1 awake');   data_cell{1+shift,2} = [1:3000] ;   
data_cell{2+shift,1} = ('Control day 1 sleep');    data_cell{2+shift,2} = [3001:14200] ;  

data_cell{3+shift,1} = ('Sod day 1 homecage');    data_cell{3+shift,2} = [14201:17200] ; 
data_cell{4+shift,1} = ('Sod day 1 sensory');    data_cell{4+shift,2} = [17201:20200] ; 
data_cell{5+shift,1} = ('Sod day 1 awake rest');    data_cell{5+shift,2} = [20201:23200] ; 
data_cell{6+shift,1} = ('Sod day 1 sleep');    data_cell{6+shift,2} = [23201:35200] ; 
data_cell{7+shift,1} = ('Sod day 6 homecage');    data_cell{7+shift,2} = [35201:38200] ; 
data_cell{8+shift,1}= ('Sod day 6 sensory');  data_cell{8+shift,2}= [38201:41200] ; 
data_cell{9+shift,1}= ('Sod day 6 awake rest'); data_cell{9+shift,2}= [41201:44200] ;
data_cell{10+shift,1}= ('Sod day 6 sleep'); data_cell{10+shift,2}= [44201:71200] ;
data_cell{11+shift,1}= ('SI homecage'); data_cell{11+shift,2}= [71201:74200] ;
data_cell{12+shift,1}= ('SI empty'); data_cell{12+shift,2}= [74201:80200] ;
data_cell{13+shift,1}= ('SI b6'); data_cell{13+shift,2}= [80201:86200];
data_cell{14+shift,1}= ('SI ICR'); data_cell{14+shift,2}= [86201:95200] ;
data_cell{15+shift,1}= ('SI awake rest'); data_cell{15+shift,2}= [95201:98200];
data_cell{16+shift,1}= ('SI sleep'); data_cell{16+shift,2}= [98201:125000] ;
data_cell{17+shift,1}= ('CFC homecage'); data_cell{17+shift,2}= [125001:128000] ;
data_cell{18+shift,1}= ('CFC preshock'); data_cell{18+shift,2}= [128001:130360] ;
data_cell{19+shift,1}= ('CFC postshock'); data_cell{19+shift,2}= [130361:132200] ;
data_cell{20+shift,1}= ('CFC awake rest'); data_cell{20+shift,2}= [132201:135200] ;
data_cell{21+shift,1}= ('CFC sleep'); data_cell{21+shift,2}= [135201:162000] ;
data_cell{22+shift,1}= ('CFC ret homecage'); data_cell{22+shift,2}= [162001:165000] ;
data_cell{23+shift,1}= ('CFC ret'); data_cell{23+shift,2}= [165001:174600] ;
data_cell{24+shift,1}= ('CFC ret awake rest'); data_cell{24+shift,2}= [174601:177600] ;
data_cell{25+shift,1}= ('CFC ret sleep'); data_cell{25+shift,2}= [177601:195400] ;


%% data processing

% zscore: enable, cut-off: disable (Recommend zscore_mode:1, cut_off_mode:0)

prms_MPPCA.bin_frame_num  = 20  ; 
prms_MPPCA.sample_fps     = 20 ; 
prms_MPPCA.filter_mode    = 0  ; 
prms_MPPCA.zscore_mode    = 1  ; 
prms_MPPCA.cut_off_mode   = 1  ; 
prms_MPPCA.mu_mode        = 0  ; 

[result_ca_filt_data_ct, result_data_cell] = fxn_MPPCA_ca_data_process_v4(data_cell, ca_raw_data, prms_MPPCA);

%% Parameters

prms_MPPCA.prms_forced_IC              = 0;        %  
prms_MPPCA.prms_ICA_mode               = 2;        % 1:fastica, 2:fastICA, 3:Reconstruction ICA (RICA)
prms_MPPCA.prms_SD_thr                 = 2;      
prms_MPPCA.prms_reactivation_SD_thr    = 2;        
prms_MPPCA.prms_target_region1         = [0:0] ;  
prms_MPPCA.prms_target_region2         = [0:0] ;   
prms_MPPCA.prms_ticklabel_mode         = 1 ;      


% Session info
prms_MPPCA.prms_reference_session_num  = 23; %
prms_MPPCA.prms_extra_session_mode     = 1 ; % 
prms_MPPCA.prms_extra_session_num      = 26 ; % input 1st session num to exclude for backprojection. Latter session also will be removed. 
prms_MPPCA.prms_figure_all_mode        = 1 ; % 0:disable figure all, 1:enable figure all.

% ### single mode for showing all figures ### Calculate MPPCA-ICA with backprojection onto all data

[result_MPPCA] = fxn_MPPCA_ICAv11(result_ca_filt_data_ct, result_data_cell, prms_MPPCA); % calculation


%% Statistical analysis
thr_freq_search = 2; 
stat_ref = 23; % reference section
stat_tar = 18; % target section
thr_stat_fold = 0; % def=0, I recommend zero.

result_MPPCA_stat = fxn_MPPCA_ICA_stat(result_MPPCA, thr_freq_search, stat_ref, stat_tar, thr_stat_fold);


 

