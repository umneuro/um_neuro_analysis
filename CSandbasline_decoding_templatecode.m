
%% Preparing binning

clc; clear; close all; disp('Previous data is cleared')

%% 

data_cell = {}; shift=1; % Don't modify this line
% addpath('homemade_functions\')

%% load data set

load('ca1_tracesod2_2023.mat') % Ca data

% % loading defeat pure neuron ID of sleep
load('ca1tracesod2_2023_sleep_defeatneurons_pure.mat');

%% load ca data 

% % 3.Extracting defeat pure neurons of sleep

ca_raw_data = ca1_tracesod2_2023(:,reactivated_defeatneurons_pure);

bin_frame_num    = 20;   % 20   -> 1s, 1s binning, For PCA-ICA, I reccomend to not change here, first.


%% %%  session information

data_cell{1+shift,1} = ('session 1');   data_cell{1+shift,2} = [1:3000] ;   
data_cell{2+shift,1} = ('session 2');    data_cell{2+shift,2} = [3001:23000] ;  
data_cell{3+shift,1} = ('session 3');   data_cell{3+shift,2} = [23001:26000]; 
data_cell{4+shift,1} = ('session 4');   data_cell{4+shift,2} = [26001:47800] ; 

data_cell{5+shift,1} = ('session 5');    data_cell{5+shift,2} = [47801:50800] ; 
data_cell{6+shift,1} = ('session 6');    data_cell{6+shift,2} = [50801:70800] ; 
data_cell{7+shift,1} = ('session 7');    data_cell{7+shift,2} = [70801:73800] ; 
data_cell{8+shift,1} = ('session 8');    data_cell{8+shift,2} = [73801:101000] ; 

data_cell{9+shift,1} = ('session 9');    data_cell{9+shift,2} = [101001:104000] ; 
data_cell{10+shift,1}= ('session 10');  data_cell{10+shift,2}= [104001:140000] ; 
data_cell{11+shift,1}= ('session 11'); data_cell{11+shift,2}= [140001:143000] ;
data_cell{12+shift,1}= ('session 12'); data_cell{12+shift,2}= [143001:173200] ;

data_cell{13+shift,1}= ('session 13'); data_cell{13+shift,2}= [173201:176200] ;
data_cell{14+shift,1}= ('session 14'); data_cell{14+shift,2}= [176201:196200] ;

data_cell{15+shift,1}= ('session 15'); data_cell{15+shift,2}= [196201:199200] ;
data_cell{16+shift,1}= ('session 16'); data_cell{16+shift,2}= [199201:204000] ;
data_cell{17+shift,1}= ('session 17'); data_cell{17+shift,2}= [204001:208800] ;
data_cell{18+shift,1}= ('session 18'); data_cell{18+shift,2}= [208801:211800] ;
data_cell{19+shift,1}= ('session 19'); data_cell{19+shift,2}= [211801:244400] ;

data_cell{20+shift,1}= ('session 20'); data_cell{20+shift,2}= [244401:249200] ;
data_cell{21+shift,1}= ('session 21'); data_cell{21+shift,2}= [249201:254000] ;

data_cell{22+shift,1}= ('session 22'); data_cell{22+shift,2}= [254001:258800] ;
data_cell{23+shift,1}= ('session 23'); data_cell{23+shift,2}= [258801:263600] ;

data_cell{24+shift,1}= ('session 24'); data_cell{24+shift,2}= [263601:268400] ;
data_cell{25+shift,1}= ('session 25'); data_cell{25+shift,2}= [268401:273200] ;

data_cell{26+shift,1}= ('session 26'); data_cell{26+shift,2}= [273201:278000] ;
data_cell{27+shift,1}= ('session 27'); data_cell{27+shift,2}= [278001:282800] ;

data_cell{28+shift,1}= ('session 28'); data_cell{28+shift,2}= [282801:285800] ;
data_cell{29+shift,1}= ('session 29'); data_cell{29+shift,2}= [285801:290600] ;
data_cell{30+shift,1}= ('session 30'); data_cell{30+shift,2}= [290601:295400] ;
data_cell{31+shift,1}= ('session 31'); data_cell{31+shift,2}= [295401:298400] ;
data_cell{32+shift,1}= ('session 32'); data_cell{32+shift,2}= [298401:345000] ;

data_cell{33+shift,1}= ('session 33'); data_cell{33+shift,2}= [345001:348000] ;
data_cell{34+shift,1}= ('session 34'); data_cell{34+shift,2}= [348001:393600] ;
data_cell{35+shift,1}= ('session 35'); data_cell{35+shift,2}= [393601:396600] ;
data_cell{36+shift,1}= ('session 36'); data_cell{36+shift,2}= [396601:420600] ;

data_cell{37+shift,1}= ('session 37'); data_cell{37+shift,2}= [420601:423600] ;
data_cell{38+shift,1}= ('session 38'); data_cell{38+shift,2}= [423601:459600] ;

%% data processing
% This code includes z-score processing, even you select either "yes" or "no" to filter fluctuation.
[result_ca_filt_data_ct, result_data_cell] = fxn_MPPCA_ca_data_process(data_cell, ca_raw_data, bin_frame_num);


%% Extract the target session (trace recall)


target_initial_recall = 38; % Session number of ca imaging binnned session

recall_t = transpose(result_data_cell{target_initial_recall+1 , 7});

recall_t_unbinned = transpose(result_data_cell{target_initial_recall+1 , 4});


% % Only CS+ trials 

CS_trial_cell_recall{1,1} = ('CS');    CS_trial_cell_recall{1,2} = recall_t(:,(416:460)) ;  
CS_trial_cell_recall{2,1} = ('CS');    CS_trial_cell_recall{2,2} = recall_t(:,(496:540)) ;  
CS_trial_cell_recall{3,1} = ('CS');    CS_trial_cell_recall{3,2} = recall_t(:,(576:620)) ;  
CS_trial_cell_recall{4,1} = ('CS');    CS_trial_cell_recall{4,2} = recall_t(:,(656:700)) ;  
CS_trial_cell_recall{5,1} = ('CS');    CS_trial_cell_recall{5,2} = recall_t(:,(736:780)) ;  
CS_trial_cell_recall{6,1} = ('CS');    CS_trial_cell_recall{6,2} = recall_t(:,(816:860)) ;  
CS_trial_cell_recall{7,1} = ('CS');    CS_trial_cell_recall{7,2} = recall_t(:,(896:940)) ; 


% % Only CS+ labels
     
CS_trial_cell_recall{1,3} = ('CS+');   
CS_trial_cell_recall{2,3} = ('CS+');    
CS_trial_cell_recall{3,3} = ('CS+');    
CS_trial_cell_recall{4,3} = ('CS+');    
CS_trial_cell_recall{5,3} = ('CS+');     
CS_trial_cell_recall{6,3} = ('CS+');     
CS_trial_cell_recall{7,3} = ('CS+'); 



number_of_trials_recall = numel(CS_trial_cell_recall(:,1));



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 



%% 1-2 Apply breaking bins for each trials (vertically)

k = 2;    % 2 for untouched binning (1 sec) or 4 for (5 or 10 sec etc bin)

CS_trialdata_at_each_binning_recall = {};

for i = 1:number_of_trials_recall ;   
CS_trialdata_at_each_binning_recall{i,1} = num2cell(CS_trial_cell_recall{i,k},1); % results{i,2} ca data


end   


%% 1-3 Spread all bins of trials in one big spreasheet (cell array sheet)

data_at_each_binning2_recall = vertcat(CS_trialdata_at_each_binning_recall{:});

[a, number_of_cstimebins_recall] = size(data_at_each_binning2_recall);


%% 1-4 Combine all trials per bin (each population activity of trials at each time point)
CS_trials_per_bin_recall= num2cell(data_at_each_binning2_recall,1);


%% 2. Type in each time frames of pre-CS baseline from Target session each time point

  
% % baseline for analysis per bins CS+ only pre-cs -20 to -16 second 
% % 
baseline_cell_recall{1,1} = ('pre-cs');    baseline_cell_recall{1,2} = recall_t(:,401:405) ;  
baseline_cell_recall{2,1} = ('pre-cs');    baseline_cell_recall{2,2} = recall_t(:,481:485) ;  
baseline_cell_recall{3,1} = ('pre-cs');    baseline_cell_recall{3,2} = recall_t(:,561:565) ;  
baseline_cell_recall{4,1} = ('pre-cs');    baseline_cell_recall{4,2} = recall_t(:,641:645) ;  
baseline_cell_recall{5,1} = ('pre-cs');    baseline_cell_recall{5,2} = recall_t(:,721:725) ;  
baseline_cell_recall{6,1} = ('pre-cs');    baseline_cell_recall{6,2} = recall_t(:,801:805) ;  
baseline_cell_recall{7,1} = ('pre-cs');    baseline_cell_recall{7,2} = recall_t(:,881:885) ; 


number_of_baseline_trials_recall = numel(baseline_cell_recall(:,1));


%% 2-1 Replicate the baseline trials to the size of same time bins of trials (if you want to compare baseline to US or Trace, match the timebins of these trials)

% IN case of you use more than 1 second bin, you need mean value of that
% baseline binning then proceed

for i = 1:number_of_baseline_trials_recall ; 
    baseline_cell_recall{i,3} = mean(baseline_cell_recall{i,2},2);
end
 

%% 2-2

baseline_cell_v2_recall = num2cell(baseline_cell_recall);

replicated_baseline_data_per_bin_recall = {};

for i = 1:number_of_baseline_trials_recall;   
replicated_baseline_data_per_bin_recall{i,1} = repmat(baseline_cell_v2_recall{i,3},1,number_of_cstimebins_recall); % results{i,2} ca data


end   


%% 2-3 Spread the replicates of baseline data cell into one spreadsheet 

Spreaded_baseline_replicates_recall = vertcat(replicated_baseline_data_per_bin_recall{:});

%% 2-4 Combine all each trials per bin (each population activity of trials at each time point)
baseline_per_bin_recall = num2cell(Spreaded_baseline_replicates_recall,1);

%% 2-5 Extract trial labels (identity labels)

CS_trial_label_recall = []; baseline_label_recall = []; CS_trial_identity_recall = [];

CS_trial_label_recall = CS_trial_cell_recall(:,1); baseline_label_recall = baseline_cell_recall(:,1); CS_trial_identity_recall = CS_trial_cell_recall(:,3);


%% 3 Combine all trials (baseline trials and cs trials at each time bins)

% 1-1 Combine CS and Baseline

combine_trials_per_bin_recall = vertcat(data_at_each_binning2_recall, Spreaded_baseline_replicates_recall);

combined_labels_recall = vertcat(CS_trial_label_recall, baseline_label_recall);


% 1-2. Combine X and Y and shuffling orders


% all samples (CS+ CS- Baseline with mean population vectors

X_Y_recall_combined_per_binning = cat(2,combine_trials_per_bin_recall,combined_labels_recall);
% 
[row1,column1] = size(X_Y_recall_combined_per_binning);
% 
random_rows1 = transpose(randperm(row1));
% 
X_Y_recall_combined_per_binning_shuffled = X_Y_recall_combined_per_binning(random_rows1,:);

X_recall_per_binning_shuffled_order = X_Y_recall_combined_per_binning_shuffled(:,1:column1-1);

Y_recall_per_binning_shuffled_order = X_Y_recall_combined_per_binning_shuffled(:,column1);


alltrials_per_bin_recall = num2cell(X_recall_per_binning_shuffled_order,1);



%% Loop for SVM each time bin

% prepare data

% 1. prepare X and Y which are sample data and labels

all_trials_per_bin_recall_cat = [];


for ii = 1:number_of_cstimebins_recall ;   
all_trials_per_bin_recall_cat{1,ii} = transpose(horzcat(X_recall_per_binning_shuffled_order{:,ii})); % results{i,2} ca data


end 


%% 2. prepare X and Y which are sample data and label

X4 = all_trials_per_bin_recall_cat;
Y4 = Y_recall_per_binning_shuffled_order;


% basic kernel

SVM_results_table = [];
% 
for ii = 1:number_of_cstimebins_recall ;   
SVM_results_table{1,ii} = fitcsvm(X4{1,ii},Y4,'KernelFunction','linear','Standardize',true,'ClassNames',{'CS','pre-cs'}); % results{i,2} ca data


end 


%% 3. Cross validate each binning

number_of_CV_repeats = 30;

for c = 1:number_of_CV_repeats;
 for ii = 1:number_of_cstimebins_recall ; 

Crossvalidation_results_table{c,ii} = crossval(SVM_results_table{1,ii},KFold=10);

 end
end



%% kfoldloss values
kfoldloss_results_table = [];

for c = 1:number_of_CV_repeats;
 for ii = 1:number_of_cstimebins_recall ; 
    

    kfoldloss_results_table{c,ii} = kfoldLoss(Crossvalidation_results_table{c,ii});
    
 end
end


%% Mean value of kfoldloss of repeated kfold Cross validation


Mean_of_kfoldloss = mean(cell2mat(kfoldloss_results_table));

Prediction_accuracy = 1-Mean_of_kfoldloss;


% 1. Preparing graph (untouched binning)

starting_point_timebins = -5;

last_point_timebins = starting_point_timebins+number_of_cstimebins_recall-1;

format long;xaxis = starting_point_timebins:last_point_timebins;



figure

plot(xaxis,Prediction_accuracy,'Color',[0,0.7,0.9]);


xlim([starting_point_timebins last_point_timebins])
ylim([0 1])

xline(0,'--r',{'CS onset'})

title('Kernel SVM Decoder (CS+/Baseline)')
xlabel('Time (seconds)')
ylabel('Accuracy')

