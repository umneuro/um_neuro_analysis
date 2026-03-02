function [Cells_sum] = fxn_MPPCA_cal_connect_4patterns(temp_ca_data, P1_unique, P2_unique, P3_unique, P4_unique, binarize_thr)
%%
temp_ca_A = temp_ca_data(:, P1_unique);
temp_ca_B = temp_ca_data(:, P2_unique);
temp_ca_C = temp_ca_data(:, P3_unique);
temp_ca_D = temp_ca_data(:, P4_unique);

temp_ca_A_logi = temp_ca_A > binarize_thr;
temp_ca_B_logi = temp_ca_B > binarize_thr;
temp_ca_C_logi = temp_ca_C > binarize_thr;
temp_ca_D_logi = temp_ca_D > binarize_thr;

Cells_A_sum = sum(temp_ca_A_logi, 2);
Cells_B_sum = sum(temp_ca_B_logi, 2);
Cells_C_sum = sum(temp_ca_C_logi, 2);
Cells_D_sum = sum(temp_ca_D_logi, 2);
%%
Cells_sum{1,1}='AB'; Cells_sum{2,1}='AC'; Cells_sum{3,1}='AD'; Cells_sum{4,1}='BC';
Cells_sum{5,1}='BD'; Cells_sum{6,1}='CD'; Cells_sum{7,1}='ABC'; Cells_sum{8,1}='ABD'; 
Cells_sum{9,1}='ACD'; Cells_sum{10,1}='BCD'; Cells_sum{11,1}='ABCD';

for i = 1:size(Cells_A_sum,1)
Cells_sum{1,2}(i,1) = Cells_A_sum(i)*Cells_B_sum(i);
Cells_sum{2,2}(i,1) = Cells_A_sum(i)*Cells_C_sum(i);
Cells_sum{3,2}(i,1) = Cells_A_sum(i)*Cells_D_sum(i);
Cells_sum{4,2}(i,1) = Cells_B_sum(i)*Cells_C_sum(i);
Cells_sum{5,2}(i,1) = Cells_B_sum(i)*Cells_D_sum(i);
Cells_sum{6,2}(i,1) = Cells_C_sum(i)*Cells_D_sum(i);
Cells_sum{7,2}(i,1) = Cells_A_sum(i)*Cells_B_sum(i)*Cells_C_sum(i);
Cells_sum{8,2}(i,1) = Cells_A_sum(i)*Cells_B_sum(i)*Cells_D_sum(i);
Cells_sum{9,2}(i,1) = Cells_A_sum(i)*Cells_C_sum(i)*Cells_D_sum(i);
Cells_sum{10,2}(i,1) = Cells_B_sum(i)*Cells_C_sum(i)*Cells_D_sum(i);
Cells_sum{11,2}(i,1) = Cells_A_sum(i)*Cells_B_sum(i)*Cells_C_sum(i)*Cells_D_sum(i);
end

% ### cal old Cells_sum_t_index
Cells_sum{1,3} = Cells_sum{1,2}/(numel(P1_unique)*numel(P2_unique));
Cells_sum{2,3} = Cells_sum{2,2}/(numel(P1_unique)*numel(P3_unique));
Cells_sum{3,3} = Cells_sum{3,2}/(numel(P1_unique)*numel(P4_unique));
Cells_sum{4,3} = Cells_sum{4,2}/(numel(P2_unique)*numel(P3_unique));
Cells_sum{5,3} = Cells_sum{5,2}/(numel(P2_unique)*numel(P4_unique));
Cells_sum{6,3} = Cells_sum{6,2}/(numel(P3_unique)*numel(P4_unique));
Cells_sum{7,3} = Cells_sum{7,2}/(numel(P1_unique)*numel(P2_unique)*numel(P3_unique));
Cells_sum{8,3} = Cells_sum{8,2}/(numel(P1_unique)*numel(P2_unique)*numel(P4_unique));
Cells_sum{9,3} = Cells_sum{9,2}/(numel(P1_unique)*numel(P3_unique)*numel(P4_unique));
Cells_sum{10,3} = Cells_sum{10,2}/(numel(P2_unique)*numel(P3_unique)*numel(P4_unique));
Cells_sum{11,3} = Cells_sum{11,2}/(numel(P1_unique)*numel(P2_unique)*numel(P3_unique)*numel(P4_unique));

% ### cal old Cells_AB_sum_index
% Cells_AB_sum_index   = sum(Cells_AB_sum_t_index)/size(temp_ca_A,1);

for i = 1:11
Cells_sum{i,4} = sum(Cells_sum{i,3})/size(temp_ca_A,1);
Cells_sum{i,1}(isnan(Cells_sum{i,1}))=0; % NaN -> zero conversion
end

%% check results
% figure;
% plot(Cells_sum{2,2}); hold on
% plot(Cells_sum{2,3}); hold on
% plot(Cells_sum{2,4}); hold on
%% referrence
% for i = 1:size(Cells_A_sum,1)
% Cells_ABC_sum_res(i,1) = Cells_A_sum(i)*Cells_B_sum(i)*Cells_C_sum(i);
% end
%%
end