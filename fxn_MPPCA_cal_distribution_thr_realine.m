%% Distribution thresholding
function [Cells_sum_shuffle_res] = fxn_MPPCA_cal_distribution_thr_realine(Cells_sum_shuffle, iteration, thr_percentile);
% iteration = 1000;
% thr_percentile = 1; % 1%-thresholding
Cells_sum_shuffle_index = {};
distribution_realine = {};
distribution_thr_bottom = {};
distribution_thr_top = {};

for i = 1:11 % till ABCD
    for ii = 1:iteration
    Cells_sum_shuffle_index{i,1}(ii,1) = Cells_sum_shuffle{ii, 1}{i, 4};
    Cells_sum_shuffle_index{i,1}(isnan(Cells_sum_shuffle_index{i,1}))=0; % NaN -> zero conversion
    end

%     Cells_sum_temp(isnan(Cells_sum_temp))=0;
distribution_realine{i,1} = sort(Cells_sum_shuffle_index{i,1});
distribution_thr = (thr_percentile/100)*iteration;
if isinteger(distribution_thr) == 0;
distribution_thr = ceil(distribution_thr);
elseif isinteger(distribution_thr) == 1;
end
distribution_thr_bottom{i,1} = distribution_realine{i}(distribution_thr);
distribution_thr_top{i,1} = distribution_realine{i}(end-distribution_thr+1);
end

%% Results
Cells_sum_shuffle_res.Cells_sum_shuffle_index = Cells_sum_shuffle_index;
Cells_sum_shuffle_res.distribution_thr_bottom = distribution_thr_bottom;
Cells_sum_shuffle_res.distribution_thr_top    = distribution_thr_top;
%%
end
%%