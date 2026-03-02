%% Statistical comparison among MPPCA ensembles across behavioral sesisons
function [result_MPPCA_stat] = fxn_MPPCA_ICA_stat(result_MPPCA, thr_freq_search, stat_ref, stat_tar, thr_stat_fold)
%% update note
% 211104: 1st code
%% debug mode
% thr_freq_search = 2; % def=2 SD, SD threshold for frequency cal. McHugh uses Raw RS 5.
% stat_ref = 16; % 
% stat_tar = 15; % 
% thr_stat_fold = 4; % 
%% data processing
result_data_cell_stat = result_MPPCA.result_data_cell(:,[1,7,12,13,14]);
num_of_session = size(result_data_cell_stat,1)-1;
num_of_ensemble = numel(result_MPPCA.neuron_sig_IDs);
shift = 1;

result_data_cell_stat{1,6} = ('Mean'); result_data_cell_stat{1,7} = ('Frq search');
result_data_cell_stat{1,8} = ('Cut off RS'); result_data_cell_stat{1,9} = ('Cut off zRS');
result_data_cell_stat{1,9} = ('Cut off Mean'); result_data_cell_stat{1,10} = ('Mean fold change');

result_data_cell_stat{1,11} = ('Max'); result_data_cell_stat{1,12} = ('zero pad'); 
result_data_cell_stat{1,13} = ('Local maxima'); result_data_cell_stat{1,14} = ('Pad off local maxima');
result_data_cell_stat{1,15} = ('Maxima address'); result_data_cell_stat{1,16} = ('Maxima val');
result_data_cell_stat{1,17} = ('Num of maxima'); result_data_cell_stat{1,18} = ('Nomarized maxima');
result_data_cell_stat{1,19} = ('Binalized above thr'); result_data_cell_stat{1,20} = ('Num of raster');
result_data_cell_stat{1,21} = ('Frequency'); result_data_cell_stat{1,21} = ('Freq mean fold');
%% Mean cal
% 4 = raw_r_kt, 5 = zscored r_kt.

for i = 1:num_of_session
    result_data_cell_stat{shift+i,6} = mean(result_data_cell_stat{shift+i,5},1); % mean
    
    result_data_cell_stat{shift+i,7} = find(result_data_cell_stat{shift+i,5} < thr_freq_search );
    result_data_cell_stat{shift+i,8} = result_data_cell_stat{shift+i,5};
    result_data_cell_stat{shift+i,8}(result_data_cell_stat{shift+i,7}) = 0; % replace to zero

    result_data_cell_stat{shift+i,9} = mean(result_data_cell_stat{shift+i,8},1); % mean

    result_data_cell_stat{shift+i,11} = max(result_data_cell_stat{shift+i,5},[],1); % Max
    
    for ii = 1:num_of_ensemble
    result_data_cell_stat{shift+i,12}(:,ii) = [0; (result_data_cell_stat{shift+i,8}(:,ii));0] ;
    result_data_cell_stat{shift+i,13}(:,ii) = double(islocalmax(result_data_cell_stat{shift+i,12}(:,ii)));
    result_data_cell_stat{shift+i,14}(:,ii) = result_data_cell_stat{shift+i,13}(2:end-1,ii);
    result_data_cell_stat{shift+i,15}{1,ii} = find(result_data_cell_stat{shift+i,14}(:,ii) == 1);
    result_data_cell_stat{shift+i,16}{1,ii} = result_data_cell_stat{shift+i,8}( (result_data_cell_stat{shift+i,15}{1,ii}),ii);
    result_data_cell_stat{shift+i,17}(:,ii) = numel(result_data_cell_stat{shift+i,16}{1,ii}); 
    result_data_cell_stat{shift+i,18}(:,ii) = sum(result_data_cell_stat{shift+i,16}{1,ii})./ size(result_data_cell_stat{shift+i,8},1);
    
    % binarize for frequency cal
    result_data_cell_stat{shift+i,19} = (result_data_cell_stat{shift+i,8} > 0);
    result_data_cell_stat{shift+i,20} = sum(result_data_cell_stat{shift+i,19});
    result_data_cell_stat{shift+i,21} = (result_data_cell_stat{shift+i,20}) ./ size(result_data_cell_stat{shift+i,8},1);
    end
end
%% stat for reactivation strength
for i = 1:num_of_session   
    result_data_cell_stat{shift+i,10} = result_data_cell_stat{shift+i,9} ./ result_data_cell_stat{shift+stat_ref,9}; % mean   
end

temp_fold_change_stat_RS = {};

for i = 1:num_of_ensemble
temp_fold_change_stat_RS{1,i} = result_data_cell_stat{shift+stat_ref,9}(i) ./ result_data_cell_stat{shift+stat_tar,9}(i);
[temp_fold_change_stat_RS{2,i}, temp_fold_change_stat_RS{3,i}, temp_fold_change_stat_RS{4,i}] = ...
    ranksum(result_data_cell_stat{shift+stat_ref,8}(:,i), result_data_cell_stat{shift+stat_tar,8}(:,i));  
temp_fold_change_stat_RS{5,i} =  temp_fold_change_stat_RS{1,i} >= thr_stat_fold && temp_fold_change_stat_RS{3,i} == 1;
end
%% stat for frequency
for i = 1:num_of_session   
    result_data_cell_stat{shift+i,22} = result_data_cell_stat{shift+i,21} ./ result_data_cell_stat{shift+stat_ref,21}; % mean   
end

temp_fold_change_stat_Frq = {};

for i = 1:num_of_ensemble
temp_fold_change_stat_Frq{1,i} = result_data_cell_stat{shift+stat_ref,21}(i) ./ result_data_cell_stat{shift+stat_tar,21}(i);
[temp_fold_change_stat_Frq{2,i}, temp_fold_change_stat_Frq{3,i}, temp_fold_change_stat_Frq{4,i}] = ...
    ranksum(result_data_cell_stat{shift+stat_ref,19}(:,i), result_data_cell_stat{shift+stat_tar,19}(:,i));  
temp_fold_change_stat_Frq{5,i} =  temp_fold_change_stat_Frq{1,i} >= thr_stat_fold && temp_fold_change_stat_Frq{3,i} == 1;
end
%% Sorting by stat results
% reactivation mean
% posi_ID = find(cell2mat(temp_fold_change_stat_RS(5,:)) == 1);
% nega_ID = 1:num_of_ensemble;
% nega_ID(posi_ID) = [];

% reactivaton frequency
posi_ID = find(cell2mat(temp_fold_change_stat_Frq(5,:)) == 1);
nega_ID = 1:num_of_ensemble;
nega_ID(posi_ID) = [];



%% result stat mean
result_stat_mean = cell2mat(result_data_cell_stat(2:end,9));
result_stat_mean_posi = mean(result_stat_mean(:,posi_ID),2);
result_stat_mean_nega = mean(result_stat_mean(:,nega_ID),2);

figure('Position',[400 100 1400 800]); % imagesc(result_stat_mean);
subplot(221)
plot(result_stat_mean(:,posi_ID),'m','LineWidth',1)
hold on
plot(result_stat_mean(:,nega_ID),'c','LineWidth',1)
hold on
plot(result_stat_mean_posi,'r','LineWidth',2)
hold on
plot(result_stat_mean_nega,'b','LineWidth',2)

xticklabels(result_data_cell_stat(1+shift:size(result_data_cell_stat,1),1) );  % for adding ticklabel
xticks([1:size(result_data_cell_stat,1)-1]); % for adding ticklabel
xtickangle(45) % update 211104
ylabel('Reactivation strength (SD)');
% ylim([0 12])

set(findobj(gcf, 'Type', 'Axes'), 'FontSize', 9, 'FontName','Arial');
hold off; box off
%% result stat maxima
result_stat_maxima = cell2mat(result_data_cell_stat(2:end,18));
result_stat_maxima_posi = mean(result_stat_maxima(:,posi_ID),2);
result_stat_maxima_nega = mean(result_stat_maxima(:,nega_ID),2);

% figure('Position',[400 100 1000 600]); % imagesc(result_stat_maxima);
subplot(222)
plot(result_stat_maxima(:,posi_ID),'m','LineWidth',1)
hold on
plot(result_stat_maxima(:,nega_ID),'c','LineWidth',1)
hold on
plot(result_stat_maxima_posi,'r','LineWidth',2)
hold on
plot(result_stat_maxima_nega,'b','LineWidth',2)

xticklabels(result_data_cell_stat(1+shift:size(result_data_cell_stat,1),1) );  % for adding ticklabel
xticks([1:size(result_data_cell_stat,1)-1]); % for adding ticklabel
xtickangle(45) % update 211104
ylabel('Mean of reactivation maxima (SD)');
% ylim([0 0.5])

set(findobj(gcf, 'Type', 'Axes'), 'FontSize', 8, 'FontName','Arial');
hold off; box off
%% Max
result_stat_max_only = cell2mat(result_data_cell_stat(2:end,11));
result_stat_max_only_posi = mean(result_stat_max_only(:,posi_ID),2);
result_stat_max_only_nega = mean(result_stat_max_only(:,nega_ID),2);

% figure('Position',[400 100 1000 600]); % imagesc(result_stat_max_only);
subplot(223)
plot(result_stat_max_only(:,posi_ID),'m','LineWidth',1)
hold on
plot(result_stat_max_only(:,nega_ID),'c','LineWidth',1)
hold on
plot(result_stat_max_only_posi,'r','LineWidth',2)
hold on
plot(result_stat_max_only_nega,'b','LineWidth',2)

xticklabels(result_data_cell_stat(1+shift:size(result_data_cell_stat,1),1) );  % for adding ticklabel
xticks([1:size(result_data_cell_stat,1)-1]); % for adding ticklabel
xtickangle(45) % update 211104
ylabel('Reactivation max only (SD)');
% ylim([0 0.5])

set(findobj(gcf, 'Type', 'Axes'), 'FontSize', 8, 'FontName','Arial');
hold off; box off
%% Sorting by stat results
posi_ID = find(cell2mat(temp_fold_change_stat_Frq(5,:)) == 1);
nega_ID = 1:num_of_ensemble;
nega_ID(posi_ID) = [];
%% Frequency above threshold
result_stat_frq = cell2mat(result_data_cell_stat(2:end,21));
result_stat_frq_posi = mean(result_stat_frq(:,posi_ID),2);
result_stat_frq_nega = mean(result_stat_frq(:,nega_ID),2);

% figure('Position',[400 100 1000 600]); % imagesc(result_stat_frq);
subplot(224)
plot(result_stat_frq(:,posi_ID),'m','LineWidth',1)
hold on
plot(result_stat_frq(:,nega_ID),'c','LineWidth',1)
hold on
plot(result_stat_frq_posi,'r','LineWidth',2)
hold on
plot(result_stat_frq_nega,'b','LineWidth',2)

xticklabels(result_data_cell_stat(1+shift:size(result_data_cell_stat,1),1) );  % for adding ticklabel
xticks([1:size(result_data_cell_stat,1)-1]); % for adding ticklabel
xtickangle(45) % update 211104
ylabel('Reactivation frequency (Hz)');
% ylim([0 0.5])

set(findobj(gcf, 'Type', 'Axes'), 'FontSize', 8, 'FontName','Arial');
hold off; box off
%% result part
result_MPPCA_stat.fold_change_stat_table       = temp_fold_change_stat_RS;
result_MPPCA_stat.result_stat_mean_nega        = result_stat_mean_nega;
result_MPPCA_stat.result_stat_mean_posi        = result_stat_mean_posi;
result_MPPCA_stat.result_stat_mean_nega_indivi = result_stat_mean(:,nega_ID);
result_MPPCA_stat.result_stat_mean_posi_indivi = result_stat_mean(:,posi_ID);
result_MPPCA_stat.result_stat_maxima           = result_stat_maxima;
result_MPPCA_stat.result_stat_max_only         = result_stat_max_only;
result_MPPCA_stat.result_stat_mean             = result_stat_mean; % reactiovation strength, used generally.
result_MPPCA_stat.result_stat_frq              = result_stat_frq;  % reactiovation frequency, used generally.
%%
end