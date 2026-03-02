%% Circshift shuffling
function shuffled_data = fxn_ensemble_shuffling(data_input,iteration,shuffle_mode)

%% 
disp([' ## Generating shuffled surrogate data by ' ,num2str(iteration) ' iterations. mode=3 Individual-shuffling reccomended.']);

[shift_max, ~] = size(data_input); % time
%% circshift
%Y = circshift(A,K,dim)
rng = 1; 
% data_shuffled = {};
%% data shuffling


shuffled_data = cell(iteration,1); % for fast cal.

if shuffle_mode == 1
%     disp('randperm shuffling')
        for i = 1:iteration   
%             disp(['Under shuffling. ' ,num2str(i),'/',num2str(iteration),' iterations.']);
                   shift_each = randperm(shift_max);
            shuffled_data{i,1} = data_input(shift_each,:);
        end

% ### circshift shuffling,
elseif shuffle_mode == 2
    disp('Fixed circshift shuffling for ensemble circshift. for each ensemble shuffling. Abdou et al.')
    for i = 1:iteration
%         disp(['Under shuffling. ' ,num2str(i),'/',num2str(iteration),' iterations.']);
        shift_each = randi(shift_max);
            shuffled_data{i,1} = circshift(data_input, shift_each, 1);
    end     

% ### update circshift shuffling, ver230316 % individual circshift
elseif shuffle_mode == 3
    disp('Individual circshift shuffling. ver231126')
    for i = 1:iteration
%         disp(['Under shuffling. ' ,num2str(i),'/',num2str(iteration),' iterations.']);
        for ii = 1:size(data_input,2)
            shift_each = randi(shift_max); % individual circshift
            shuffled_data{i, 1}(:,ii) = circshift(data_input(:,ii), shift_each, 1); % individual circshift
        end
    end     
end
%%
disp(['### Finished the Generating surrogate data. ###']);
%%
end
