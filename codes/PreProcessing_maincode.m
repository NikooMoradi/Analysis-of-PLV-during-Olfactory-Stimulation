%% 1 : table to array
clear;clc;
% Load the .mat file
load('Subject2.mat');
% Check if the loaded data is a table
if istable(subject2)
    % Convert the table to an array
    data_array = table2array(subject2); 
    disp('done');
else
    error('The loaded variable is not a table.');
end
% Compute transpose
transposed_array = transpose(data_array);
% Save transposed array to a new .mat file
save('Subject2_array.mat', 'transposed_array');

%% 2 : remove the 20th channel

% Load the .mat file
load('Subject2_array.mat'); 
if exist('transposed_array', 'var')
    % Check the size of the data array
    [numChannels, ~] = size(transposed_array);
    if numChannels == 20
        % Remove the 20th channel
        transposed_array(20, :) = [];
        % Save the modified data to a new .mat file
        save('subject2_array_selected.mat', ...
            'transposed_array');
        disp('done')
    else
        error(['The loaded variable does not ' ...
            'contain 20 channels of data.']);
    end
else
    error('The loaded variable does not exist in the .mat file.');
end

%% step 3 (epoch)
x = EEG.data;

result = zeros(19, 600, 120);
for i = 1 : 120
    startIndex = (i - 1) * 2000 + 1400;
    endIndex = startIndex + 600 - 1;
    result(:, :, i) = x(:,startIndex:endIndex);
end
% save('Subject2_epoched.mat', 'result');

%% Step 4 (optinal)

Subject1_epoched = load('Subject1_epoched.mat');
Subject2_epoched = load('Subject2_epoched.mat');
Subject1_epoched = Subject1_epoched.result;
Subject2_epoched = Subject2_epoched.result;
[channels, samples, trials] = size(Subject2_epoched);
noisyTrials = [];

% Iterate over each channel and calculate power spectrum of each trial
for ch = 1 : channels
    
    p = nan(trials, samples); 
    for trial = 1:trials
        trialData = squeeze(Subject2_epoched(ch, :, trial)); % Extract
        %data for the current trial
        power = fft(trialData).^2; % power spectrum for the current channel
        p(trial, :) = power;
    end    
    % suggested commands
    variances = sum(nanstd(p, [], 2).^2, 2);
    noisyTrialsCh = find(abs(zscore(variances)) > 3.5);   
    noisyTrials = union(noisyTrials, noisyTrialsCh);
end

cleanEpochData = Subject2_epoched(:, :, setdiff(1:trials, noisyTrials));
disp(['Noisy Trials: ' num2str(length(noisyTrials))]);
disp('Clean Epoch Dimensions:');
disp(size(cleanEpochData));
save('clean_epoch_data_2.mat', 'cleanEpochData');


%% Step 5: Select and Save channels

Subject1_epoched_cleandata2 = load('Subject1_epoched_cleandata2.mat');
Subject2_epoched_cleandata2 = load('Subject2_epoched_cleandata2.mat');
Subject1_epoched_cleandata2 = Subject1_epoched_cleandata2.x;
Subject2_epoched_cleandata2 = Subject2_epoched_cleandata2.x;

% Define the channels to keep 
channelsToKeep = [1, 5, 10, 15];
numChannels = 19;
% Initialize the channels index array with zeros
channelsIdx = zeros(1, numChannels);
% Set the selected channels to 1 in the index array
for i = channelsToKeep
    channelsIdx(i) = 1;
end
% Keep only the desired channels
selectedEpochData = Subject2_epoched_cleandata2(logical(channelsIdx),:,:);
% Save the selected epoch data
% save('Subject2_final.mat', 'selectedEpochData');
disp(['Epoch data has been saved to ', 'Subject2_final.mat']);


%% Step 5 (Struct): Save Selected Epoch Data
Subject1_final = load('Subject1_final.mat');
Subject2_final = load('Subject2_final.mat');
Subject1_final = Subject1_final.selectedEpochData;
Subject2_final = Subject2_final.selectedEpochData;

Normal = load('Normal.mat','normal');
odor = Normal.normal.odor;
noisy = Normal.normal.noisy;
% Create a struct and assign the matrices as fields
myStruct = struct(...
    'Subject1_final', Subject1_final, ...
    'odor', odor, ...
    'noisy',noisy);

% Save the selected epoch data
% save('subject1_struct.mat', 'myStruct');
disp(['Selected epoch data has been saved to ', 'subject1_struct.mat']);


