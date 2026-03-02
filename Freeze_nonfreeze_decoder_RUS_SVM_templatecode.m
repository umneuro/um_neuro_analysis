
%% 1-1 Freezing vector from BehaviorDEPOT and trimout for actual session

load('Behavior_ca1tracesod2recall_filt.mat') % freezing raw data behaviordepot

ca1tracesod2_recall_freezing_behaviordepot = Behavior.Freezing.Vector;

session_start = 60; 

session_duration = 1800; % example duration in seconds
trimrange = getTrimRange(session_start, session_duration);

ca1tracesod2_recall_freezing_behaviordepot_trimmed = ca1tracesod2_recall_freezing_behaviordepot(trimrange);


%% 1-2 Binning into second

% % 1.2 freezing data from behaviordepot
freeze_data_for_binning = ca1tracesod2_recall_freezing_behaviordepot_trimmed';

binning_num = 20;

[freeze_data_binned,  freeze_data_binned_mean] = funcHF_temporal_binning_freezingdata(freeze_data_for_binning, binning_num);


%% 1-3 Rounding decimal values of binned freezing data 


freeze_data_binned = round(freeze_data_binned);

freeze_data_binned = transpose(freeze_data_binned);


%% 1-4 THIS IS FOR CS+ trials only; Extract freezing data during CS+Trace+ (40 second)

accept_frames = [420:459 500:539 580:619 660:699 740:779 820:859 900:939];

% Measure total freezing during all CS+ trials
 freeze_data_binned_during_cstrace  = freeze_data_binned(accept_frames);
 total_duration_of_cstrace = numel(freeze_data_binned_during_cstrace);

 total_freezing_for_whole_trial_period = sum(freeze_data_binned_during_cstrace)/total_duration_of_cstrace*100;


%% 1-5 Extract CS-Trace frames from Ca imaging recall data

% % FIRST YOU NEED recall_t file from original code

 ca_data_binned_during_cstrace  = recall_t(:,accept_frames);

 ca_data_binned_during_cstrace_unbinned  = recall_t_unbinned(:,accept_frames_unbinned);


%% 1-6 Divide freeze/nonfreeze timepoints for Ca response and labels (whole duration)

freezing_timeframes_cstrace = find(freeze_data_binned_during_cstrace);
nonfreezing_timeframes_cstrace = find(freeze_data_binned_during_cstrace == 0);

freezing_labels_cstrace    = cell(1, numel(freezing_timeframes_cstrace));
freezing_labels_cstrace(:) = {'Freeze'};

nonfreezing_labels_cstrace    = cell(1, numel(nonfreezing_timeframes_cstrace));
nonfreezing_labels_cstrace(:) = {'Non-freeze'};

freezing_nonfreezing_all_labels = horzcat(freezing_labels_cstrace,nonfreezing_labels_cstrace);


%% 1-7 Divide Ca response with freeze and nonfreeze then sort (whole duration during CS and trace)

% Extract and divide Ca response of freeze or no freeze state
ca_data_binned_during_cs_trace_freeze = ca_data_binned_during_cstrace(:,freezing_timeframes_cstrace);
ca_data_binned_during_cs_trace_nonfreezing = ca_data_binned_during_cstrace(:,nonfreezing_timeframes_cstrace);

% Combine these response sorted order from freeze to nonfreeze responses

ca_data_binned_freeze_nonfreeze_sorted = horzcat(ca_data_binned_during_cs_trace_freeze,ca_data_binned_during_cs_trace_nonfreezing);

x = num2cell(ca_data_binned_freeze_nonfreeze_sorted');
y = freezing_nonfreezing_all_labels';

%% 1-8 shuffling x predictors and y labels orders and SVM (class imbalance of freeze and non-freeze so use ROC)

X_Y_recall_combined_per_binning = cat(2,x,y);

[row1,column1] = size(X_Y_recall_combined_per_binning);

random_rows1 = transpose(randperm(row1));
% 
X_Y_recall_combined_per_binning_shuffled = X_Y_recall_combined_per_binning(random_rows1,:);

X_shuffled_order = cell2mat(X_Y_recall_combined_per_binning_shuffled(:,1:column1-1));

Y_shuffled_order = X_Y_recall_combined_per_binning_shuffled(:,column1);


% % shuffling y labels only for control
[row3,column3] = size(Y_shuffled_order);

random_rows3 = transpose(randperm(row3));

Y_shuffled_control = Y_shuffled_order(random_rows3,:);


%% 1-9 SVM implementing rnadom undersampling

% Convert Y_shuffled_order to logical (0 and 1)
Y_shuffled_order_ylogical = double(logical(Y_shuffled_order == "Freeze"));

% Combine X and Y into a single dataset for the RUS function
imbalanced_data = [X_shuffled_order, Y_shuffled_order_ylogical];

% Get dimensions
[row1, numFeatures] = size(X_shuffled_order);

% Determine minority and majority classes
classLabels = unique(Y_shuffled_order_ylogical);
classCounts = histcounts(Y_shuffled_order_ylogical, length(classLabels));

disp(['Non-freeze (0): ', num2str(classCounts(1)), ', Freeze (1): ', num2str(classCounts(2))]);

if classCounts(1) < classCounts(2)
    minorityClass = classLabels(1); % Class 0 is minority
    majorityClass = classLabels(2); % Class 1 is majority
else
    minorityClass = classLabels(2); % Class 1 is minority
    majorityClass = classLabels(1); % Class 0 is majority
end

disp(['Minority Class: ', num2str(minorityClass)]);
disp(['Majority Class: ', num2str(majorityClass)]);

minorityClassSize = sum(Y_shuffled_order_ylogical == minorityClass);
Samples_No = floor(minorityClassSize * 0.9); % Use 90% of minority class size

% Parameters
number_of_CV_repeats = 30;
numFolds = 10;

% Create a 10-fold cross-validation partition
cv = cvpartition(Y_shuffled_order_ylogical, 'KFold', numFolds);

% Initialize arrays to store results
lossResults = zeros(numFolds, number_of_CV_repeats);
f1Results = zeros(numFolds, number_of_CV_repeats);
rocAUC_Freeze = zeros(numFolds, number_of_CV_repeats);
rocAUC_Nonfreeze = zeros(numFolds, number_of_CV_repeats);

for c = 1:number_of_CV_repeats
    usedSamples = []; % Track used samples for the current repeat
    for i = 1:numFolds
        % Get training and test indices
        trainIdx = cv.training(i);
        testIdx = cv.test(i);
        
        % Separate training data
        trainData = imbalanced_data(trainIdx, :);
        
        % Apply Random Undersampling to training data
        random_samples = RUS(trainData, Samples_No); % Assumes RUS returns cell array of datasets
        
        % Filter out used samples
        available_samples = setdiff(1:length(random_samples), usedSamples);
        if isempty(available_samples)
            error('No more available samples to use in repeat %d, fold %d.', c, i);
        end
        
        % Randomly select one of the remaining samples
        selectedSampleIndex = available_samples(randi(length(available_samples)));
        newTrainData = random_samples{selectedSampleIndex};
        
        % Update usedSamples
        usedSamples = [usedSamples, selectedSampleIndex];
        
        % Create new training and testing datasets
        X_train = newTrainData(:, 1:end-1);
        Y_train = newTrainData(:, end);
        X_test = imbalanced_data(testIdx, 1:end-1);
        Y_test = imbalanced_data(testIdx, end);
        
        % Train the SVM model
        SVM_model = fitcsvm(X_train, Y_train, 'KernelFunction', 'linear', ...
                            'Standardize', true, 'ClassNames', [0, 1]);
        
        % Predict on test set
        [predictedLabels, validationScores] = predict(SVM_model, X_test);
        
        % Calculate confusion matrix
        confMat = confusionmat(Y_test, predictedLabels, 'Order', [1, 0]);
        
        % Extract TP, FP, FN for F1 score (class 1 = Freeze)
        TP = confMat(1, 1); % True Positives
        FP = confMat(2, 1); % False Positives
        FN = confMat(1, 2); % False Negatives
        
        % Calculate precision and recall with NaN handling
        if (TP + FP) == 0
            precision = 0;
        else
            precision = TP / (TP + FP);
        end
        if (TP + FN) == 0
            recall = 0;
        else
            recall = TP / (TP + FN);
        end
        
        % Calculate F1 Score
        if (precision + recall) == 0
            f1Results(i, c) = 0;
        else
            f1Results(i, c) = 2 * (precision * recall) / (precision + recall);
        end
        
        % Calculate ROC-AUC for Freeze (class 1) and Non-freeze (class 0)
        [~, ~, ~, auc] = perfcurve(Y_test, validationScores(:, 2), 1);
        rocAUC_Freeze(i, c) = auc;
        [~, ~, ~, auc] = perfcurve(Y_test, validationScores(:, 1), 0);
        rocAUC_Nonfreeze(i, c) = auc;
        
        % Calculate loss (error rate)
        lossResults(i, c) = mean(predictedLabels ~= Y_test);
    end
end

% Calculate averages
avg_kfold_loss = mean(mean(lossResults));
avg_f1_score = mean(mean(f1Results));
avg_roc_auc_freeze = mean(mean(rocAUC_Freeze));
avg_roc_auc_nonfreeze = mean(mean(rocAUC_Nonfreeze));
prediction_accuracy = 1 - avg_kfold_loss;

% Display results
disp('--- RUS Original Data (True Labels) ---');
disp(['Average Prediction Accuracy: ', num2str(prediction_accuracy)]);
disp(['Average F1 Score (Freeze): ', num2str(avg_f1_score)]);
disp(['Average ROC-AUC (Freeze): ', num2str(avg_roc_auc_freeze)]);
disp(['Average ROC-AUC (Non-freeze): ', num2str(avg_roc_auc_nonfreeze)]);








