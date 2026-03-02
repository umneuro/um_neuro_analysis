%% MPPCA-ICA calculation code
clc; clear; close all; disp('Previous data is cleared')
data_cell = {}; shift=1; % Don't modify this line
% addpath('homemade_functions\') 
%% load data set

load('ca_mouse6.mat') % Ca data

%% load ca data 

ca_raw_data      = mousesod6; % input time x neuron matrix

bin_frame_num    = 20;   % 20   -> 1s, 1s binning, For PCA-ICA, I recommend to not change here, first.

%% %%  mouse sod 6 session frame information 

data_cell{1+shift,1} = ('Control day 1 awake');   data_cell{1+shift,2} = [1:3000] ;   
data_cell{2+shift,1} = ('Control day 1 NREM 1');    data_cell{2+shift,2} = [3001:7800] ;  
data_cell{3+shift,1} = ('Control day 1 REM 1');    data_cell{3+shift,2} = [7801:11800] ;  
data_cell{4+shift,1} = ('Control day 1 NREM 2');    data_cell{4+shift,2} = [11801:13000] ; 
data_cell{5+shift,1} = ('Control day 1 REM 2');    data_cell{5+shift,2} = [13001:14200] ; 

data_cell{6+shift,1} = ('Sod day 1 homecage');    data_cell{6+shift,2} = [14201:17200] ; 
data_cell{7+shift,1} = ('Sod day 1 sensory');    data_cell{7+shift,2} = [17201:20200] ; 
data_cell{8+shift,1} = ('Sod day 1 awake rest');    data_cell{8+shift,2} = [20201:23200] ; 
data_cell{9+shift,1} = ('Sod day 1 NREM 1');    data_cell{9+shift,2} = [23201:33400] ; 
data_cell{10+shift,1} = ('Sod day 1 REM 1');    data_cell{10+shift,2} = [33401:35200] ; 

data_cell{11+shift,1} = ('Sod day 6 homecage');    data_cell{11+shift,2} = [35201:38200] ; 
data_cell{12+shift,1}= ('Sod day 6 sensory');  data_cell{12+shift,2}= [38201:41200] ; 
data_cell{13+shift,1}= ('Sod day 6 awake rest'); data_cell{13+shift,2}= [41201:44200] ;
data_cell{14+shift,1}= ('Sod day 6 NREM 1'); data_cell{14+shift,2}= [44201:59200] ;
data_cell{15+shift,1}= ('Sod day 6 REM 1'); data_cell{15+shift,2}= [59201:59800] ;
data_cell{16+shift,1}= ('Sod day 6 NREM 2'); data_cell{16+shift,2}= [59801:61200] ;
data_cell{17+shift,1}= ('Sod day 6 REM 2'); data_cell{17+shift,2}= [61201:66200] ;
data_cell{18+shift,1}= ('Sod day 6 NREM 3'); data_cell{18+shift,2}= [66201:67600] ;
data_cell{19+shift,1}= ('Sod day 6 REM 3'); data_cell{19+shift,2}= [67601:71200] ;
data_cell{20+shift,1}= ('SI homecage'); data_cell{20+shift,2}= [71201:74200] ;
data_cell{21+shift,1}= ('SI empty'); data_cell{21+shift,2}= [74201:80200] ;
data_cell{22+shift,1}= ('SI b6'); data_cell{22+shift,2}= [80201:86200];
data_cell{23+shift,1}= ('SI ICR'); data_cell{23+shift,2}= [86201:95200] ;
data_cell{24+shift,1}= ('SI awake rest'); data_cell{24+shift,2}= [95201:98200];
data_cell{25+shift,1}= ('SI NREM 1'); data_cell{25+shift,2}= [98201:106200] ;
data_cell{26+shift,1}= ('SI REM 1'); data_cell{26+shift,2}= [106201:109200] ;
data_cell{27+shift,1}= ('SI NREM 2'); data_cell{27+shift,2}= [109201:116600] ;
data_cell{28+shift,1}= ('SI REM 2'); data_cell{28+shift,2}= [116601:119600] ;
data_cell{29+shift,1}= ('SI NREM 3'); data_cell{29+shift,2}= [119601:122200] ;
data_cell{30+shift,1}= ('SI REM 3'); data_cell{30+shift,2}= [122201:122800] ;

data_cell{31+shift,1}= ('CFC homecage'); data_cell{31+shift,2}= [125001:128000] ;
data_cell{32+shift,1}= ('CFC '); data_cell{32+shift,2}= [128001:132200] ;
data_cell{33+shift,1}= ('CFC awake rest'); data_cell{33+shift,2}= [132201:135200] ;
data_cell{34+shift,1}= ('CFC NREM 1'); data_cell{34+shift,2}= [135201:145400] ;
data_cell{35+shift,1}= ('CFC REM 1'); data_cell{35+shift,2}= [145401:147600] ;
data_cell{36+shift,1}= ('CFC NREM 2'); data_cell{36+shift,2}= [147601:156200] ;
data_cell{37+shift,1}= ('CFC REM 2'); data_cell{37+shift,2}= [156201:158800] ;
data_cell{38+shift,1}= ('CFC NREM 3'); data_cell{38+shift,2}= [158801:159000] ;
data_cell{39+shift,1}= ('CFC REM 3'); data_cell{39+shift,2}= [159001:162000] ;
data_cell{40+shift,1}= ('CFC ret homecage'); data_cell{40+shift,2}= [162001:165000] ;
data_cell{41+shift,1}= ('CFC ret'); data_cell{41+shift,2}= [165001:174600] ;
data_cell{42+shift,1}= ('CFC ret awake rest'); data_cell{42+shift,2}= [174601:177600] ;
data_cell{43+shift,1}= ('CFC ret NREM 1'); data_cell{43+shift,2}= [177601:181800] ;

data_cell{44+shift,1}= ('CFC ret REM 1'); data_cell{44+shift,2}= [181801:184800] ;
data_cell{45+shift,1}= ('CFC ret NREM 2'); data_cell{45+shift,2}= [184801:188000] ;
data_cell{46+shift,1}= ('CFC ret REM 2'); data_cell{46+shift,2}= [188001:191000] ;
data_cell{47+shift,1}= ('CFC ret NREM 3'); data_cell{47+shift,2}= [191001:193200] ;
data_cell{48+shift,1}= ('CFC ret REM3 '); data_cell{48+shift,2}= [193201:195400] ;

%% data processing

[result_ca_filt_data_ct, result_data_cell] = fxn_MPPCA_ca_data_process(data_cell, ca_raw_data, bin_frame_num);
 
%% Parameters
prms_forced_IC              = 0;       
prms_ICA_mode               = 2;       
prms_SD_thr                 = 2;      
prms_reactivation_SD_thr    = 2;        
prms_target_region1         = [0:0] ;  
prms_target_region2         = [0:0] ;   
prms_ticklabel_mode         = 1 ;      
%% (Double mode. Input two session values) Assign sessions

reference_session_num1 = 12; % select type 1 session
reference_session_num2 = 32; % select type 2 session

%% PCA/ICA per each session loop

data_bin_z_reference = result_data_cell{reference_session_num1+shift,7}; % input data for MPPCA-ICA
data_bin_z_target    = result_ca_filt_data_ct; % Don't change. for all-session backprojection.

[result_MPPCA1] = fxn_MPPCA_ICAv05(data_bin_z_reference, data_bin_z_target, prms_forced_IC, prms_ICA_mode, prms_SD_thr,...
                                   prms_reactivation_SD_thr, prms_target_region1, prms_target_region2, bin_frame_num, ...
                                   result_data_cell, prms_ticklabel_mode); close all;% calculation

data_bin_z_reference = result_data_cell{reference_session_num2+shift,7}; % input data for MPPCA-ICA
data_bin_z_target    = result_ca_filt_data_ct; % Don't change. for all-session backprojection.

[result_MPPCA2] = fxn_MPPCA_ICAv05(data_bin_z_reference, data_bin_z_target, prms_forced_IC, prms_ICA_mode, prms_SD_thr,...
                                   prms_reactivation_SD_thr, prms_target_region1, prms_target_region2, bin_frame_num, ...
                                   result_data_cell, prms_ticklabel_mode); close all;% calculation%%                               

%% (Double mode) Assesing cell-cell pairwise interaction 

% Don't change
prms.binarize_thr        = 3 ; % Calcium trace binalizing threshold.
prms.shuffle_mode        = 2 ;      % 1:randperm, 2:circshift
prms.ensemble_num1        = numel(result_MPPCA1.neuron_sig_IDs);
prms.ensemble_num2        = numel(result_MPPCA2.neuron_sig_IDs);
prms.thr_percentile      = 1; % def: 1%, top-bottom %-thresholding


%% Looping these process for multiple session

% % Multiple session
prms.session_num_assess  = [1:6 33:39];

[result, coactivity_result] = fxn_MPPCA_ensemble_connectivity_double_multi_session2(result_data_cell, result_MPPCA1, result_MPPCA2, prms);
