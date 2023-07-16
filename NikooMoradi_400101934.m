
% Nikoo Moradi _ 400101934
% Signals and Systems Project

%% Question 3 : Pre-processing ==========================================================================================================================================================
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

%% End of Pre-processing ==========================================================================================================================================================
%% Question 4 : Results ==========================================================================================================================================================
%% 4.1
% Load data files
normalData = load('Normal.mat');
adData = load('AD.mat');
% Define frequency range of interest
frequencyRange = [35, 40]; % Hz
% Initialize arrays to store PLV values
normalPLV = zeros(15, 2); % 15 people x 2 odors
adPLV = zeros(13, 2); 
% Loop through each participant in the Normal group
for person = 1:15
    % Access the epoch
    [~ , ~ , numTrials] = size(normalData.normal(person).epoch);  
    % 4 x 600 x NumTrials matrix
    epochData = normalData.normal(person).epoch;  
     % NumTrials x 1 binary array
    odorData = normalData.normal(person).odor; 
    % Loop through each trial and calculate PLV
    for trial = 1 : numTrials
        numLemon = sum(odorData(:) == 0);
        numRose = sum(odorData(:) == 1);
        odor = odorData(trial,1);
        
        signal1_Fz = epochData(2,:,trial);
        signal2_Cz = epochData(3,:,trial);
        
        normalPLV(person,odor+1) = normalPLV(person,odor+1) + ...
        calculatePLV(signal1_Fz, signal2_Cz, 200, frequencyRange);
    end
    % Average PLV of each odor
    normalPLV(person,1) = normalPLV(person,1)/ numLemon;
    normalPLV(person,2) = normalPLV(person,2)/ numRose;
end
% Loop through each participant in the AD group
for person = 1:13
     % Access the epoch
    [~ , ~ , numTrials] = size(adData.AD(person).epoch);  
    epochData = adData.AD(person).epoch;  % 4 x 600 x NumTrials matrix
    odorData = adData.AD(person).odor;  % NumTrials x 1 binary array
    % Loop through each trial and calculate PLV
    for trial = 1 : numTrials
        numLemon = sum(odorData(:) == 0);
        numRose = sum(odorData(:) == 1);
        odor = odorData(trial,1);
        signal1_Fz = epochData(2,:,trial);
        signal2_Cz = epochData(3,:,trial);
        adPLV(person,odor+1) = adPLV(person,odor+1) + calculatePLV(signal1_Fz, signal2_Cz, 200, frequencyRange);
    end
    % Average PLV of each odor
    adPLV(person,1) = adPLV(person,1)/ numLemon;
    adPLV(person,2) = adPLV(person,2)/ numRose;
end
% Display the results
disp('Normal Group PLV:');
disp(normalPLV);
save('normalPLV.mat', 'normalPLV');
disp('AD Group PLV:');
disp(adPLV);
save('adPLV.mat', 'adPLV');

%% 4.2
clear; clc;
adPLV = load('adPLV.mat');
normalPLV = load('normalPLV.mat');
adPLV = adPLV.adPLV;
normalPLV = normalPLV.normalPLV;
% Create group and odor labels
groupLabels = {'AD', 'Noramal'};
odorLabels = {'Lemon', 'Rose'};
figure;
for i = 1:2
    for j = 1:2
        subplot(2, 2, (i-1)*2+j);
        if (i == 1 && j == 1)
            data = adPLV(:,1);
        elseif (i == 1 && j == 2)
            data = adPLV(:,2); 
        elseif (i == 2 && j == 1)
            data = normalPLV(:,1);
        elseif (i == 2 && j == 2) 
            data = normalPLV(:,2);
        end
        boxplot(data);
        title(sprintf('Group: %s, Odor: %s', ...
        groupLabels{i}, odorLabels{j}));

    end
end

% Fit Gaussian distributions
figure;
for i = 1:2
    for j = 1:2
        subplot(2, 2, (i-1)*2+j);
        if (i == 1 && j == 1)
            data = adPLV(:,1);
        elseif (i == 1 && j == 2)
            data = adPLV(:,2); 
        elseif (i == 2 && j == 1)
            data = normalPLV(:,1);
        elseif (i == 2 && j == 2) 
            data = normalPLV(:,2);
        end
        histogram(data,'Normalization', 'pdf');
        hold on;
        x = linspace(min(data)-2, max(data)+2, 100);
        mu = mean(data);
        sigma = std(data);
        y = normpdf(x, mu, sigma);
        plot(x, y, 'r-', 'LineWidth', 2);
        hold off;
        xlabel('PLV');
        ylabel('Probability');
        title(sprintf('Group: %s, Odor: %s',...
        groupLabels{i}, odorLabels{j}));
        
    end
end

% Perform statistical tests (t-test) between groups and odors
alpha1 = abs(mean(normalPLV(:, 1))- mean(adPLV(:, 1))); % Significance level
alpha2 = abs(mean(normalPLV(:, 2))- mean(adPLV(:, 2))); % Significance level
disp(['alpha1 (Significance level for Group Comparison - Lemon Odor) = ', num2str(alpha1)]);fprintf('\n')
disp(['alpha2 (Significance level for Group Comparison - Rose Odor) = ', num2str(alpha2)]);fprintf('\n')

% Comparison between groups for Lemon odor
[h1 , pval_group_lemon] = ttest2(normalPLV(:, 1), adPLV(:, 1), 'Alpha', alpha1);

% Comparison between groups for Rose odor
[h2 ,pval_group_rose] = ttest2(normalPLV(:, 2), adPLV(:, 2), 'Alpha', alpha2);

disp('Statistical Significance (p-values):');fprintf('\n')

disp(['Group Comparison - Lemon Odor: p = ', num2str(pval_group_lemon)]);
if h1 == 1
    disp('Lemon trials are significantly different between AD and normal groups');fprintf('\n')
else
    disp('not significantly different');
end

disp(['Group Comparison - Rose Odor: p = ', num2str(pval_group_rose)]);
if h2 == 1
    disp('Rose trials are significantly different between AD and normal groups');fprintf('\n')
else
    disp('not significantly different');
end



%% 4.4

% Load data files
normalData = load('Normal.mat');
adData = load('AD.mat');
temp = 0;

% Generate a random number between a and b
randomNormalPerson = floor(1 + (15 - 1) * rand);
randomADPerson = floor(1 + (13 - 1) * rand);

%  Normal group =======================================================
% Access the epoch
    [~ ,~ , numTrials] = size(normalData.normal(randomNormalPerson).epoch);  
    % 4 x 600 x NumTrials matrix
    epochData = normalData.normal(randomNormalPerson).epoch;  
    % NumTrials x 1 binary array
    odorData = normalData.normal(randomNormalPerson).odor;  
    for trial = 1 : numTrials
        numLemon = sum(odorData(:) == 0);
        if odorData(trial,1) == 0 % we just want to do it for lemon odor
        signal1_Fz = epochData(2,:,trial);
        signal2_Cz = epochData(3,:,trial);
        % Calculate phase difference
        phase_diff = abs(angle(signal1_Fz) - angle(signal2_Cz));
        % Convert phase difference to degrees
        phase_diff_deg = rad2deg(phase_diff);
        temp = temp + phase_diff_deg;
        end
    end
    avrg_phase_diff_deg = temp / numLemon;
   % Create polar histogram
   figure;
        subplot(1,2,1);
        num_bins = 20;
        polarhistogram(avrg_phase_diff_deg, num_bins);

        % Customize the plot
        title('Polar Histogram of Phase Difference (Normal)');
        rticks([0 10 20]);
        set(gca, 'Color', 'white');
        

% =========================================================================
temp = 0;
%  AD group ===============================================================

% Access the epoch
    [~ , ~ , numTrials] = size(adData.AD(randomADPerson).epoch);  
    epochData = adData.AD(randomADPerson).epoch;  % 4 x 600 x NumTrials matrix
    odorData = adData.AD(randomADPerson).odor;  % NumTrials x 1 binary array

    
    for trial = 1 : numTrials

        numLemon = sum(odorData(:) == 0);

        if odorData(trial,1) == 0 % we just want to do it for lemon(frequent) odor

        signal1_Fz = epochData(2,:,trial);
        signal2_Cz = epochData(3,:,trial);
        
        % Calculate phase difference
        phase_diff = abs(angle(signal1_Fz) - angle(signal2_Cz));

        % Convert phase difference to degrees
        phase_diff_deg = rad2deg(phase_diff);

        temp = temp + phase_diff_deg;
        end

    end

    avrg_phase_diff_deg = temp / numLemon;

   % Create polar histogram
        subplot(1,2,2);
        num_bins = 20;
        polarhistogram(avrg_phase_diff_deg, num_bins);

        % Customize the plot
        title('Polar Histogram of Phase Difference (AD)');
        rticks([0 10 20]);
        set(gca, 'Color', 'white');
        
        
        
   % =========================================================================
   temp = 0;
   avrg_phase_diff_deg =0;
   % Mean Phase diff for Normal group =========================================================================
   for person = 1 : 15   
    [~ , ~ , numTrials] = size(normalData.normal(person).epoch);  
    % 4 x 600 x NumTrials matrix
    epochData = normalData.normal(person).epoch;  
    % NumTrials x 1 binary array
    odorData = normalData.normal(person).odor;  
    for trial = 1 : numTrials
        numLemon = sum(odorData(:) == 0);
        % we just want to do it for lemon(frequent) odor
        if odorData(trial,1) == 0 
        signal1_Fz = epochData(2,:,trial);
        signal2_Cz = epochData(3,:,trial);
        % Calculate phase difference
        phase_diff = abs(angle(signal1_Fz) - angle(signal2_Cz));
        % Convert phase difference to degrees
        phase_diff_deg = rad2deg(phase_diff);
        temp = temp + phase_diff_deg;
        end
    end
    avrg_phase_diff_deg = avrg_phase_diff_deg + temp / numLemon;
    temp = 0;
   end
   % Create polar histogram
   figure;
        subplot(1,2,1);
        num_bins = 20;
        polarhistogram(avrg_phase_diff_deg / 15, num_bins);

        % Customize the plot
        title('Polar Histogram of Mean Value of Phase Difference (Normal)');
        rticks([0 10 20]);
        set(gca, 'Color', 'white');
% =========================================================================
   temp = 0;
   avrg_phase_diff_deg =0;
   
%  Mean Phase diff for AD group ===============================================================
for person = 1 : 13
% Access the epoch
    [~ , ~ , numTrials] = size(adData.AD(person).epoch);  
    epochData = adData.AD(person).epoch;  % 4 x 600 x NumTrials matrix
    odorData = adData.AD(person).odor;  % NumTrials x 1 binary array

    
    for trial = 1 : numTrials

        numLemon = sum(odorData(:) == 0);

        if odorData(trial,1) == 0 % we just want to do it for lemon(frequent) odor

        signal1_Fz = epochData(2,:,trial);
        signal2_Cz = epochData(3,:,trial);
        
        % Calculate phase difference
        phase_diff = abs(angle(signal1_Fz) - angle(signal2_Cz));

        % Convert phase difference to degrees
        phase_diff_deg = rad2deg(phase_diff);

        temp = temp + phase_diff_deg;
        end

    end

    avrg_phase_diff_deg = avrg_phase_diff_deg + temp / numLemon;
    temp = 0;
end
   % Create polar histogram
        subplot(1,2,2);
        num_bins = 20;
        polarhistogram(avrg_phase_diff_deg / 13, num_bins);

        % Customize the plot
        title('Polar Histogram of Mean Value of Phase Difference (AD)');
        rticks([0 10 20]);
        set(gca, 'Color', 'white');
        
        
        
   % =========================================================================
  
   
%% 4.5

% Normal peaple ===========================================================
% Load data files
normalData = load('Normal.mat');

plv_normal_lemon = zeros(4,4);
plv_normal_rose = zeros(4,4);
    
% Define frequency range of interest
frequencyRange = [35, 40]; % Hz

% Initialize arrays to store PLV values
normalPLV = zeros(15, 2 , 6); % 15 people x 2 odors x 6 pair of channels 

% Loop through each participant in the Normal group
for person = 1:15

    % Access the epoch
    [~ , ~ , numTrials] = size(normalData.normal(person).epoch);  
    epochData = normalData.normal(person).epoch;  % 4 x 600 x NumTrials matrix
    odorData = normalData.normal(person).odor;  % NumTrials x 1 binary array

     % Average epochs across trials
    averagedTrials = zeros(4,600,2); %mean(epochData());  % 4 x 600 matrix x 2 odors
    
    numLemon = sum(odorData(:) == 0);
    numRose = sum(odorData(:) == 1);
    % Loop through each trial and sum the trials , Lemons with Lemons and
    % roses with roses
    for trial = 1 : numTrials
        
        odor = odorData(trial,1);
        averagedTrials(:,:,odor+1) = averagedTrials(:,:,odor+1) + epochData(:,:,trial);

    end
    
    % Average of trials
    averagedTrials(:,:,1) = averagedTrials(:,:,1) / numLemon; % 4 x 600 matrix x 2 odors (Lemon)
    averagedTrials(:,:,2) = averagedTrials(:,:,2) / numRose; % 4 x 600 matrix x 2 odors (Rose)
    
    signal1_Fp1 = averagedTrials(1,:,:); % 1 x 600  x 2 odors
    signal2_Fz = averagedTrials(2,:,:);
    signal1_Cz = averagedTrials(3,:,:);
    signal1_Pz = averagedTrials(4,:,:);
    

    normalPLV(person,1,1) =  calculatePLV(signal1_Fp1(:,:,1), signal2_Fz(:,:,1), 200, frequencyRange);
    normalPLV(person,2,1) =  calculatePLV(signal1_Fp1(:,:,2), signal2_Fz(:,:,2), 200, frequencyRange);
    
    normalPLV(person,1,2) =  calculatePLV(signal1_Fp1(:,:,1), signal1_Cz(:,:,1), 200, frequencyRange);
    normalPLV(person,2,2) =  calculatePLV(signal1_Fp1(:,:,2), signal1_Cz(:,:,2), 200, frequencyRange);
   
    
    normalPLV(person,1,3) =  calculatePLV(signal1_Fp1(:,:,1), signal1_Pz(:,:,1), 200, frequencyRange);
    normalPLV(person,2,3) =  calculatePLV(signal1_Fp1(:,:,2), signal1_Pz(:,:,2), 200, frequencyRange);
    
    
    normalPLV(person,1,4) =  calculatePLV(signal2_Fz(:,:,1), signal1_Cz(:,:,1), 200, frequencyRange);
    normalPLV(person,2,4) =  calculatePLV(signal2_Fz(:,:,2), signal1_Cz(:,:,2), 200, frequencyRange);
    
    
    normalPLV(person,1,5) =  calculatePLV(signal2_Fz(:,:,1), signal1_Pz(:,:,1), 200, frequencyRange);
    normalPLV(person,2,5) =  calculatePLV(signal2_Fz(:,:,2), signal1_Pz(:,:,2), 200, frequencyRange);
    
    
    normalPLV(person,1,6) =  calculatePLV(signal1_Cz(:,:,1), signal1_Pz(:,:,1), 200, frequencyRange);
    normalPLV(person,2,6) =  calculatePLV(signal1_Cz(:,:,2), signal1_Pz(:,:,2), 200, frequencyRange);
    
end

    plv_normal_lemon(1, 2) = abs(mean(normalPLV(:,1,1)));
    plv_normal_rose(1, 2) = abs(mean(normalPLV(:,2,1)));
    plv_normal_lemon(2,1) = plv_normal_lemon(1, 2);
    plv_normal_rose(2, 1) = plv_normal_rose(1, 2);
    
    plv_normal_lemon(1, 3) = abs(mean(normalPLV(:,1,2)));
    plv_normal_rose(1, 3) = abs(mean(normalPLV(:,2,2)));
    plv_normal_lemon(3,1) = plv_normal_lemon(1, 3);
    plv_normal_rose(3, 1) = plv_normal_rose(1, 3);
    
    plv_normal_lemon(1, 4) = abs(mean(normalPLV(:,1,3)));
    plv_normal_rose(1, 4) = abs(mean(normalPLV(:,2,3)));
    plv_normal_lemon(4,1) = plv_normal_lemon(1, 4);
    plv_normal_rose(4, 1) = plv_normal_rose(1, 4);

    plv_normal_lemon(2, 3) = abs(mean(normalPLV(:,1,4)));
    plv_normal_rose(2, 3) = abs(mean(normalPLV(:,2,4)));
    plv_normal_lemon(3,2) = plv_normal_lemon(2, 3);
    plv_normal_rose(3, 2) = plv_normal_rose(2, 3);
    
    plv_normal_lemon(2, 4) = abs(mean(normalPLV(:,1,5)));
    plv_normal_rose(2, 4) = abs(mean(normalPLV(:,2,5)));
    plv_normal_lemon(4,2) = plv_normal_lemon(2,4);
    plv_normal_rose(4, 2) = plv_normal_rose(2, 4);
    
    plv_normal_lemon(3, 4) = abs(mean(normalPLV(:,1,6)));
    plv_normal_rose(3, 4) = abs(mean(normalPLV(:,2,6)));
    plv_normal_lemon(4,3) = plv_normal_lemon(3, 4);
    plv_normal_rose(4, 3) = plv_normal_rose(3, 4);
    
    
    figure(1);

subplot(2, 2, 1)
plotHeatmap(plv_normal_lemon);
title('Normal Participants, Lemon odor');

subplot(2, 2, 2)
plotHeatmap(plv_normal_rose);
title('Normal Participants, Rose odor');

% AD peaple ==============================================================

% Load data files
adData = load('AD.mat');

plv_AD_lemon = zeros(4,4);
plv_AD_rose = zeros(4,4);
    
% Define frequency range of interest
frequencyRange = [35, 40]; % Hz

% Initialize arrays to store PLV values
adPLV = zeros(13, 2,6); 

% Loop through each participant in the Normal group
for person = 1:13

    % Access the epoch
    [~ , ~ , numTrials] = size(adData.AD(person).epoch);  
    epochData = adData.AD(person).epoch;  % 4 x 600 x NumTrials matrix
    odorData = adData.AD(person).odor;  % NumTrials x 1 binary array

     % Average epochs across trials
    averagedTrials = zeros(4,600,2); %mean(epochData());  % 4 x 600 matrix x 2 odors
    
    numLemon = sum(odorData(:) == 0);
    numRose = sum(odorData(:) == 1);
    % Loop through each trial and sum the trials , Lemons with Lemons and
    % roses with roses
    for trial = 1 : numTrials
        
        odor = odorData(trial,1);
        averagedTrials(:,:,odor+1) = averagedTrials(:,:,odor+1) + epochData(:,:,trial);

    end
    
    % Average of trials
    averagedTrials(:,:,1) = averagedTrials(:,:,1) / numLemon; % 4 x 600 matrix x 2 odors (Lemon)
    averagedTrials(:,:,2) = averagedTrials(:,:,2) / numRose; % 4 x 600 matrix x 2 odors (Rose)
    
    signal1_Fp1 = averagedTrials(1,:,:); % 1 x 600  x 2 odors
    signal2_Fz = averagedTrials(2,:,:);
    signal1_Cz = averagedTrials(3,:,:);
    signal1_Pz = averagedTrials(4,:,:);
    

    adPLV(person,1,1) =  calculatePLV(signal1_Fp1(:,:,1), signal2_Fz(:,:,1), 200, frequencyRange);
    adPLV(person,2,1) =  calculatePLV(signal1_Fp1(:,:,2), signal2_Fz(:,:,2), 200, frequencyRange);
    
    adPLV(person,1,2) =  calculatePLV(signal1_Fp1(:,:,1), signal1_Cz(:,:,1), 200, frequencyRange);
    adPLV(person,2,2) =  calculatePLV(signal1_Fp1(:,:,2), signal1_Cz(:,:,2), 200, frequencyRange);
   
    
    adPLV(person,1,3) =  calculatePLV(signal1_Fp1(:,:,1), signal1_Pz(:,:,1), 200, frequencyRange);
    adPLV(person,2,3) =  calculatePLV(signal1_Fp1(:,:,2), signal1_Pz(:,:,2), 200, frequencyRange);
    
    
    adPLV(person,1,4) =  calculatePLV(signal2_Fz(:,:,1), signal1_Cz(:,:,1), 200, frequencyRange);
    adPLV(person,2,4) =  calculatePLV(signal2_Fz(:,:,2), signal1_Cz(:,:,2), 200, frequencyRange);
    
    
    adPLV(person,1,5) =  calculatePLV(signal2_Fz(:,:,1), signal1_Pz(:,:,1), 200, frequencyRange);
    adPLV(person,2,5) =  calculatePLV(signal2_Fz(:,:,2), signal1_Pz(:,:,2), 200, frequencyRange);
    
    
    adPLV(person,1,6) =  calculatePLV(signal1_Cz(:,:,1), signal1_Pz(:,:,1), 200, frequencyRange);
    adPLV(person,2,6) =  calculatePLV(signal1_Cz(:,:,2), signal1_Pz(:,:,2), 200, frequencyRange);
    
end

    plv_AD_lemon(1, 2) = abs(mean(adPLV(:,1,1)));
    plv_AD_rose(1, 2) = abs(mean(adPLV(:,2,1)));
    plv_AD_lemon(2,1) = plv_AD_lemon(1, 2);
    plv_AD_rose(2, 1) = plv_AD_rose(1, 2);
    
    plv_AD_lemon(1, 3) = abs(mean(adPLV(:,1,2)));
    plv_AD_rose(1, 3) = abs(mean(adPLV(:,2,2)));
    plv_AD_lemon(3,1) = plv_AD_lemon(1, 3);
    plv_AD_rose(3, 1) = plv_AD_rose(1, 3);
    
    plv_AD_lemon(1, 4) = abs(mean(adPLV(:,1,3)));
    plv_AD_rose(1, 4) = abs(mean(adPLV(:,2,3)));
    plv_AD_lemon(4,1) = plv_AD_lemon(1, 4);
    plv_AD_rose(4, 1) = plv_AD_rose(1, 4);

    plv_AD_lemon(2, 3) = abs(mean(adPLV(:,1,4)));
    plv_AD_rose(2, 3) = abs(mean(adPLV(:,2,4)));
    plv_AD_lemon(3,2) = plv_AD_lemon(2, 3);
    plv_AD_rose(3, 2) = plv_AD_rose(2, 3);
    
    plv_AD_lemon(2, 4) = abs(mean(adPLV(:,1,5)));
    plv_AD_rose(2, 4) = abs(mean(adPLV(:,2,5)));
    plv_AD_lemon(4,2) = plv_AD_lemon(2,4);
    plv_AD_rose(4, 2) = plv_AD_rose(2, 4);
    
    plv_AD_lemon(3, 4) = abs(mean(adPLV(:,1,6)));
    plv_AD_rose(3, 4) = abs(mean(adPLV(:,2,6)));
    plv_AD_lemon(4,3) = plv_AD_lemon(3, 4);
    plv_AD_rose(4, 3) = plv_AD_rose(3, 4);


subplot(2, 2, 3)
plotHeatmap(plv_AD_lemon);
title('AD Patients, Lemon odor');

subplot(2, 2, 4)
plotHeatmap(plv_AD_rose);
title('AD Patients, Rose odor');



 % Significance levels
alpha1 = abs(mean(normalPLV(:,1,1))- mean(adPLV(:,1,1)));
alpha2 = abs(mean(normalPLV(:,1,2))- mean(adPLV(:,1,2)));

alpha3 = abs(mean(normalPLV(:,1,3))- mean(adPLV(:,1,3)));
alpha4 = abs(mean(normalPLV(:,1,4))- mean(adPLV(:,1,4)));

alpha5 = abs(mean(normalPLV(:,1,5))- mean(adPLV(:,1,5)));
alpha6 = abs(mean(normalPLV(:,1,6))- mean(adPLV(:,1,6)));

% Comparison between Normal and AD for Lemon odor
[h1, pval_group_Fp1_Fz] = ttest2(normalPLV(:,:,1), adPLV(:,:,1),alpha1);

% Comparison between Normal and AD for Rose odor
[h2, pval_group_Fp1_Cz] = ttest2(normalPLV(:,:,2), adPLV(:,:,2),alpha2);

% Comparison between Normal and MCI for Lemon odor
[h3 , pval_group_Fp1_Pz] = ttest2(normalPLV(:,:,3), adPLV(:,:,3),alpha3);

% Comparison between Normal and MCI for Rose odor
[h4,pval_group_Fz_Cz] = ttest2(normalPLV(:,:,4), adPLV(:,:,4),alpha4);

% Comparison between AD and MCI for Lemon odor
[h5,pval_group_Fz_Pz] = ttest2(normalPLV(:,:,5), adPLV(:,:,5),alpha5);

% Comparison between AD and MCI for Rose odor
[h6,pval_group_Cz_Pz] = ttest2(normalPLV(:,:,6), adPLV(:,:,6),alpha6);

disp('Statistical Significance (p-values):');fprintf('\n');

disp(['Group Comparison - Fp1 & Fz : p = ', num2str(pval_group_Fp1_Fz)]);
if h1 == 1
    disp('significantly different');fprintf('\n');
else
    disp('not significantly different');fprintf('\n');
end

disp(['Group Comparison - Fp1 & Cz : p = ', num2str(pval_group_Fp1_Cz)]);
if h2 == 1
    disp('significantly different');fprintf('\n');
else
    disp('not significantly different');fprintf('\n');
end

disp(['Group Comparison - Fp1 & Pz : p = ', num2str(pval_group_Fp1_Pz)]);
if h3 == 1
    disp('significantly different');fprintf('\n');
else
    disp('not significantly different');fprintf('\n');
end

disp(['Group Comparison - Fz & Cz : p = ', num2str(pval_group_Fz_Cz)]);
if h4 == 1
    disp('significantly different');fprintf('\n');
else
    disp('not significantly different');fprintf('\n');
end

disp(['Group Comparison - Fz & Pz : p = ', num2str(pval_group_Fz_Pz)]);
if h5 == 1
    disp('significantly different');fprintf('\n');
else
    disp('not significantly different');fprintf('\n');
end

disp(['Group Comparison - Cz & Pz : p =  ', num2str(pval_group_Cz_Pz)]);
if h6 == 1
    disp('significantly different');fprintf('\n');
else
    disp('not significantly different');fprintf('\n');
end

%% End of Results =========================================================================================================================================
%% Question 5 : Bonus 1 =========================================================================================================================================
%% 5.1.2 : PLV of MCI subjects
% Load data files
mciData = load('MCI.mat');

% Define frequency range of interest
frequencyRange = [35, 40]; % Hz

% Initialize arrays to store PLV values
mciPLV = zeros(7, 2); % 7 people x 2 odors

% Loop through each participant in the MCI group
for person = 1:7
    % Access the epoch
    [~ , ~ , numTrials] = size(mciData.MCI(person).epoch);  

    % 4 x 600 x NumTrials matrix
    epochData = mciData.MCI(person).epoch; 

     % NumTrials x 1 binary array
    odorData = mciData.MCI(person).odor; 

    % Loop through each trial and calculate PLV
    for trial = 1 : numTrials
        numLemon = sum(odorData(:) == 0);
        numRose = sum(odorData(:) == 1);
        odor = odorData(trial,1);
        signal1_Fz = epochData(2,:,trial);
        signal2_Cz = epochData(3,:,trial);
        mciPLV(person,odor+1) = mciPLV(person,odor+1) + ...
        calculatePLV(signal1_Fz, signal2_Cz, 200, frequencyRange);
    end
    % Average PLV of each odor
    mciPLV(person,1) = mciPLV(person,1)/ numLemon;
    mciPLV(person,2) = mciPLV(person,2)/ numRose;
end 

% Display the results
disp('MCI Group PLV:');
disp(mciPLV);
save('mciPLV.mat', 'mciPLV');



%% 5.1.2 : box plot and fitting Guassian for MCI
clear; clc;
mciPLV = load('mciPLV.mat');
mciPLV = mciPLV.mciPLV;

adPLV = load('adPLV.mat');
adPLV = adPLV.adPLV;

mciPLV = load('normalPLV.mat');
mciPLV = mciPLV.normalPLV;

% Create group and odor labels
groupLabels = {'AD', 'Noramal'};
odorLabels = {'Lemon', 'Rose'};

% box plot
figure;
subplot(1,2,1);
boxplot(mciPLV(:,1));
title('Group: MCI, Odor: Lemon');

subplot(1,2,2);
boxplot(mciPLV(:,2));
title('Group: MCI, Odor: Rose');


% Fit Gaussian distributions
figure;
subplot(1,2,1);
histogram(mciPLV(:,1),'Normalization', 'pdf');
        hold on;
        x = linspace(min(mciPLV(:,1))-1, max(mciPLV(:,1))+1, 100);
        mu = mean(mciPLV(:,1));
        sigma = std(mciPLV(:,1));
        y = normpdf(x, mu, sigma);
        plot(x, y, 'r-', 'LineWidth', 2);
        hold off;
        xlabel('PLV');
        ylabel('Probability');
        title('Group: MCI, Odor: Lemon');

subplot(1,2,2);
histogram(mciPLV(:,2),'Normalization', 'pdf');
        hold on;
        x = linspace(min(mciPLV(:,2))-2, max(mciPLV(:,2))+2, 100);
        mu = mean(mciPLV(:,2));
        sigma = std(mciPLV(:,2));
        y = normpdf(x, mu, sigma);
        plot(x, y, 'r-', 'LineWidth', 2);
        hold off;
        xlabel('PLV');
        ylabel('Probability');
        title('Group: MCI, Odor: Rose');
   
%% phase difference for MCI

% Load data files
mciData = load('MCI.mat');
temp = 0;

% Generate a random number between a and b
randomMCIPerson = floor(1 + (7 - 1) * rand);


%  MCI group =======================================================

% Access the epoch
    [~ ,~ , numTrials] = size(mciData.MCI(randomMCIPerson).epoch);  
    % 4 x 600 x NumTrials matrix
    epochData = mciData.MCI(randomMCIPerson).epoch;  
    % NumTrials x 1 binary array
    odorData = mciData.MCI(randomMCIPerson).odor;  
    
    for trial = 1 : numTrials
        numLemon = sum(odorData(:) == 0);
        if odorData(trial,1) == 0 % we just want to do it for lemon odor
        signal1_Fz = epochData(2,:,trial);
        signal2_Cz = epochData(3,:,trial);
        % Calculate phase difference
        phase_diff = abs(angle(signal1_Fz) - angle(signal2_Cz));
        % Convert phase difference to degrees
        phase_diff_deg = rad2deg(phase_diff);
        temp = temp + phase_diff_deg;
        end
    end
    
    avrg_phase_diff_deg = temp / numLemon;
    
   % Create polar histogram
   figure;
       
        num_bins = 20;
        polarhistogram(avrg_phase_diff_deg, num_bins);

        % Customize the plot
        title('Polar Histogram of Phase Difference (MCI)');
        rticks([0 10 20]);
        set(gca, 'Color', 'white');
        

% =========================================================================
temp = 0;
avrg_phase_diff_deg =0;
% Mean Phase diff for MCI group =========================================================================
   for person = 1 : 7   
    [~ , ~ , numTrials] = size(mciData.MCI(person).epoch);  
    % 4 x 600 x NumTrials matrix
    epochData = mciData.MCI(person).epoch;  
    % NumTrials x 1 binary array
    odorData = mciData.MCI(person).odor;  
    
    for trial = 1 : numTrials
        numLemon = sum(odorData(:) == 0);
        
        % we just want to do it for lemon(frequent) odor
        if odorData(trial,1) == 0 
        signal1_Fz = epochData(2,:,trial);
        signal2_Cz = epochData(3,:,trial);
        
        % Calculate phase difference
        phase_diff = abs(angle(signal1_Fz) - angle(signal2_Cz));
        
        % Convert phase difference to degrees
        phase_diff_deg = rad2deg(phase_diff);
        temp = temp + phase_diff_deg;
        
        end
    end
    
    avrg_phase_diff_deg = avrg_phase_diff_deg + temp / numLemon;
    temp = 0;
    
   end
   % Create polar histogram
   figure;
       
        num_bins = 20;
        polarhistogram(avrg_phase_diff_deg / 15, num_bins);

        % Customize the plot
        title('Polar Histogram of Mean Value of Phase Difference (MCI)');
        rticks([0 10 20]);
        set(gca, 'Color', 'white');
% =========================================================================

%% heat map for MCI

% MCI peaple ===========================================================
% Load data files
mciData = load('MCI.mat');

plv_mci_lemon = zeros(4,4);
plv_mci_rose = zeros(4,4);
    
% Define frequency range of interest
frequencyRange = [35, 40]; % Hz

% Initialize arrays to store PLV values
mciPLV = zeros(7, 2 , 6); % 7 people x 2 odors x 6 pair of channels 

% Loop through each participant in the MCI group
for person = 1:7

    % Access the epoch
    [~ , ~ , numTrials] = size(mciData.MCI(person).epoch);  
    epochData = mciData.MCI(person).epoch;  % 4 x 600 x NumTrials matrix
    odorData = mciData.MCI(person).odor;  % NumTrials x 1 binary array

     % Average epochs across trials
    averagedTrials = zeros(4,600,2); %mean(epochData());  % 4 x 600 matrix x 2 odors
    
    numLemon = sum(odorData(:) == 0);
    numRose = sum(odorData(:) == 1);
    
    % Loop through each trial and sum the trials , Lemons with Lemons and
    % roses with roses
    for trial = 1 : numTrials
        
        odor = odorData(trial,1);
        averagedTrials(:,:,odor+1) = averagedTrials(:,:,odor+1) + epochData(:,:,trial);

    end
    
    % Average of trials
    averagedTrials(:,:,1) = averagedTrials(:,:,1) / numLemon; % 4 x 600 matrix x 2 odors (Lemon)
    averagedTrials(:,:,2) = averagedTrials(:,:,2) / numRose; % 4 x 600 matrix x 2 odors (Rose)
    
    signal1_Fp1 = averagedTrials(1,:,:); % 1 x 600  x 2 odors
    signal2_Fz = averagedTrials(2,:,:);
    signal1_Cz = averagedTrials(3,:,:);
    signal1_Pz = averagedTrials(4,:,:);
    

    mciPLV(person,1,1) =  calculatePLV(signal1_Fp1(:,:,1), signal2_Fz(:,:,1), 200, frequencyRange);
    mciPLV(person,2,1) =  calculatePLV(signal1_Fp1(:,:,2), signal2_Fz(:,:,2), 200, frequencyRange);
    
    mciPLV(person,1,2) =  calculatePLV(signal1_Fp1(:,:,1), signal1_Cz(:,:,1), 200, frequencyRange);
    mciPLV(person,2,2) =  calculatePLV(signal1_Fp1(:,:,2), signal1_Cz(:,:,2), 200, frequencyRange);
   
    
    mciPLV(person,1,3) =  calculatePLV(signal1_Fp1(:,:,1), signal1_Pz(:,:,1), 200, frequencyRange);
    mciPLV(person,2,3) =  calculatePLV(signal1_Fp1(:,:,2), signal1_Pz(:,:,2), 200, frequencyRange);
    
    
    mciPLV(person,1,4) =  calculatePLV(signal2_Fz(:,:,1), signal1_Cz(:,:,1), 200, frequencyRange);
    mciPLV(person,2,4) =  calculatePLV(signal2_Fz(:,:,2), signal1_Cz(:,:,2), 200, frequencyRange);
    
    
    mciPLV(person,1,5) =  calculatePLV(signal2_Fz(:,:,1), signal1_Pz(:,:,1), 200, frequencyRange);
    mciPLV(person,2,5) =  calculatePLV(signal2_Fz(:,:,2), signal1_Pz(:,:,2), 200, frequencyRange);
    
    
    mciPLV(person,1,6) =  calculatePLV(signal1_Cz(:,:,1), signal1_Pz(:,:,1), 200, frequencyRange);
    mciPLV(person,2,6) =  calculatePLV(signal1_Cz(:,:,2), signal1_Pz(:,:,2), 200, frequencyRange);
    
end

    plv_mci_lemon(1, 2) = abs(mean(mciPLV(:,1,1)));
    plv_mci_rose(1, 2) = abs(mean(mciPLV(:,2,1)));
    plv_mci_lemon(2,1) = plv_mci_lemon(1, 2);
    plv_mci_rose(2, 1) = plv_mci_rose(1, 2);
    
    plv_mci_lemon(1, 3) = abs(mean(mciPLV(:,1,2)));
    plv_mci_rose(1, 3) = abs(mean(mciPLV(:,2,2)));
    plv_mci_lemon(3,1) = plv_mci_lemon(1, 3);
    plv_mci_rose(3, 1) = plv_mci_rose(1, 3);
    
    plv_mci_lemon(1, 4) = abs(mean(mciPLV(:,1,3)));
    plv_mci_rose(1, 4) = abs(mean(mciPLV(:,2,3)));
    plv_mci_lemon(4,1) = plv_mci_lemon(1, 4);
    plv_mci_rose(4, 1) = plv_mci_rose(1, 4);

    plv_mci_lemon(2, 3) = abs(mean(mciPLV(:,1,4)));
    plv_mci_rose(2, 3) = abs(mean(mciPLV(:,2,4)));
    plv_mci_lemon(3,2) = plv_mci_lemon(2, 3);
    plv_mci_rose(3, 2) = plv_mci_rose(2, 3);
    
    plv_mci_lemon(2, 4) = abs(mean(mciPLV(:,1,5)));
    plv_mci_rose(2, 4) = abs(mean(mciPLV(:,2,5)));
    plv_mci_lemon(4,2) = plv_mci_lemon(2,4);
    plv_mci_rose(4, 2) = plv_mci_rose(2, 4);
    
    plv_mci_lemon(3, 4) = abs(mean(mciPLV(:,1,6)));
    plv_mci_rose(3, 4) = abs(mean(mciPLV(:,2,6)));
    plv_mci_lemon(4,3) = plv_mci_lemon(3, 4);
    plv_mci_rose(4, 3) = plv_mci_rose(3, 4);
    
    
    figure(1);

subplot(1, 2, 1)
plotHeatmap(plv_mci_lemon);
title('MCI Participants, Lemon odor');

subplot(1, 2, 2)
plotHeatmap(plv_mci_rose);
title('MCI Participants, Rose odor');


%% 5.1.2 : p-value and significance testing on the 3 states (Normal & AD & MCI)

mciPLV = load('mciPLV.mat');
mciPLV = mciPLV.mciPLV;

adPLV = load('adPLV.mat');
adPLV = adPLV.adPLV;

normalPLV = load('normalPLV.mat');
normalPLV = normalPLV.normalPLV;


% Perform statistical tests (t-test) between groups and odors

 % Significance levels
alpha1 = abs(mean(normalPLV(:, 1))- mean(adPLV(:, 1)));
alpha2 = abs(mean(normalPLV(:, 2))- mean(adPLV(:, 2))); 

alpha3 = abs(mean(normalPLV(:, 1))- mean(mciPLV(:, 1)));
alpha4 = abs(mean(normalPLV(:, 2))- mean(mciPLV(:, 2)));

alpha5 = abs(mean(mciPLV(:, 1))- mean(adPLV(:, 1)));
alpha6 = abs(mean(mciPLV(:, 2))- mean(adPLV(:, 2)));

% Comparison between Normal and AD for Lemon odor
[h1, pval_normal_AD_lemon] = ttest2(normalPLV(:, 1), adPLV(:, 1),alpha1);

% Comparison between Normal and AD for Rose odor
[h2, pval_normal_AD_rose] = ttest2(normalPLV(:, 2), adPLV(:, 2),alpha2);

% Comparison between Normal and MCI for Lemon odor
[h3 , pval_normal_MCI_lemon] = ttest2(normalPLV(:, 1), mciPLV(:, 1),alpha3);

% Comparison between Normal and MCI for Rose odor
[h4,pval_normal_MCI_rose] = ttest2(normalPLV(:, 2), mciPLV(:, 2),alpha4);

% Comparison between AD and MCI for Lemon odor
[h5,pval_MCI_AD_lemon] = ttest2(adPLV(:, 1), mciPLV(:, 1),alpha5);

% Comparison between AD and MCI for Rose odor
[h6,pval_MCI_AD_rose] = ttest2(adPLV(:, 2), mciPLV(:, 2),alpha6);

disp('Statistical Significance (p-values):');fprintf('\n')
disp(['Comparison between Normal & AD - Lemon Odor: p = ', num2str(pval_normal_AD_lemon)]);
if h1 == 1
    disp('Lemon trials are significantly different between AD and Normal groups');fprintf('\n')
else
    disp('not significantly different');
end

disp(['Comparison between Normal & AD - Rose Odor: p = ', num2str(pval_normal_AD_rose)]);
if h2 == 1
    disp('Rose trials are significantly different between AD and Normal groups');fprintf('\n')
else
    disp('not significantly different');
end

disp(['Comparison between Normal & MCI - Lemon Odor: p = ', num2str(pval_normal_MCI_lemon)]);
if h3 == 1
    disp('Lemon trials are significantly different between MCI and Normal groups');
else
    disp('not significantly different');fprintf('\n')
end

disp(['Comparison between Normal & MCI - Rose Odor: p = ', num2str(pval_normal_MCI_rose)]);
if h4 == 1
    disp('Rose trials are significantly different between MCI and Normal groups');
else
    disp('not significantly different');fprintf('\n')
end

disp(['Comparison between MCI & AD - Lemon Odor: p = ', num2str(pval_MCI_AD_lemon)]);
if h5 == 1
    disp('Lemon trials are significantly different between AD and MCI groups');
else
    disp('not significantly different');fprintf('\n')
end

disp(['Comparison between MCI & AD - Rose Odor: p = ', num2str(pval_MCI_AD_rose)]);
if h6 == 1
    disp('Rose trials are significantly different between AD and MCI groups');
else
    disp('not significantly different');fprintf('\n')
end



%% End of Bonus 1 =========================================================================================================================================
%% Question 5 : Bonus 2 =========================================================================================================================================
%% 5.2.2 : calculating MI for Normal & AD
% Load data files
normalData = load('Normal.mat');
adData = load('AD.mat');
nbins = 18;
% helpMI = zeros( 600 , 600); % 15 people x MI param(600x600)
normalMI = zeros(13, 2);

% Loop through each participant in the Normal group
for person = 1:13
    % Access the epoch
    [~ , ~ , numTrials] = size(normalData.normal(person).epoch);  
    % 4 x 600 x NumTrials matrix
    epochData = normalData.normal(person).epoch;  
     % NumTrials x 1 binary array
    odorData = normalData.normal(person).odor; 
    % Loop through each trial and calculate PLV
    for trial = 1 : numTrials
        numLemon = sum(odorData(:) == 0);
        numRose = sum(odorData(:) == 1);
        odor = odorData(trial,1);

        signal1_Fz = epochData(2,:,trial);
        signal2_Cz = epochData(3,:,trial);
        
        temp = get_mi(angle(signal1_Fz), abs(signal2_Cz), nbins);
        normalMI(person,odor + 1) = normalMI(person,odor + 1) + ...
        mean(mean(temp.MI));
        
    end
    % Average PLV of each odor
    normalMI(person, 1) = normalMI(person,1)/ numLemon;
    normalMI(person, 2) = normalMI(person,2)/ numRose;
    
end

adMI = zeros(11, 2);
% Loop through each participant in the Normal group
for person = 1:11
    % Access the epoch
    [~ , ~ , numTrials] = size(adData.AD(person).epoch);  
    % 4 x 600 x NumTrials matrix
    epochData = adData.AD(person).epoch;  
     % NumTrials x 1 binary array
    odorData = adData.AD(person).odor; 
    % Loop through each trial and calculate PLV
    for trial = 1 : numTrials
        numLemon = sum(odorData(:) == 0);
        numRose = sum(odorData(:) == 1);
        odor = odorData(trial,1);
        
        signal1_Fz = epochData(2,:,trial);
        signal2_Cz = epochData(3,:,trial);
        
        temp = get_mi(angle(signal1_Fz), abs(signal2_Cz), nbins);
        adMI(person,odor + 1) = adMI(person,odor + 1) + ...
            mean(mean(temp.MI));
        
    end
    % Average PLV of each odor
    adMI(person, 1) = adMI(person,1)/ numLemon;
    adMI(person, 2) = adMI(person,2)/ numRose;
    
end

disp('Normal Group MI:');
disp(normalMI);
save('normalMI.mat', 'normalMI');
disp('AD Group MI:');
disp(adMI);
save('adMI.mat', 'adMI');

%% 5.2.2 : box plot & guassian and p-value for MI
clear; clc;
adMI = load('adMI.mat');
normalMI = load('normalMI.mat');
adMI = adMI.adMI;
normalMI = normalMI.normalMI;
% Create group and odor labels
groupLabels = {'AD', 'Noramal'};
odorLabels = {'Lemon', 'Rose'};
figure;
for i = 1:2
    for j = 1:2
        subplot(2, 2, (i-1)*2+j);
        if (i == 1 && j == 1)
            data = adMI(:,1);
        elseif (i == 1 && j == 2)
            data = adMI(:,2); 
        elseif (i == 2 && j == 1)
            data = normalMI(:,1);
        elseif (i == 2 && j == 2) 
            data = normalMI(:,2);
        end
        boxplot(data);
        title(sprintf('Group: %s, Odor: %s', ...
        groupLabels{i}, odorLabels{j}));

    end
end

% Fit Gaussian distributions
figure;
for i = 1:2
    for j = 1:2
        subplot(2, 2, (i-1)*2+j);
        if (i == 1 && j == 1)
            data = adMI(:,1);
        elseif (i == 1 && j == 2)
            data = adMI(:,2); 
        elseif (i == 2 && j == 1)
            data = normalMI(:,1);
        elseif (i == 2 && j == 2) 
            data = normalMI(:,2);
        end
        histogram(data,'Normalization', 'pdf');
        hold on;
        x = linspace(min(data)-0.5, max(data)+0.5, 100);
        mu = mean(data);
        sigma = std(data);
        y = normpdf(x, mu, sigma);
        plot(x, y, 'r-', 'LineWidth', 2);
        hold off;
        xlabel('MI');
        ylabel('Probability');
        title(sprintf('Group: %s, Odor: %s',...
        groupLabels{i}, odorLabels{j}));
        
    end
end

% Perform statistical tests (t-test) between groups and odors
alpha1 = abs(mean(normalMI(:, 1))- mean(adMI(:, 1))); % Significance level
alpha2 = abs(mean(normalMI(:, 2))- mean(adMI(:, 2))); % Significance level
disp(['alpha1 (Significance level for Group Comparison - Lemon Odor) = ', num2str(alpha1)]);fprintf('\n')
disp(['alpha2 (Significance level for Group Comparison - Rose Odor) = ', num2str(alpha2)]);fprintf('\n')

% Comparison between groups for Lemon odor
[h1 , pval_group_lemon] = ttest2(normalMI(:, 1), adMI(:, 1), 'Alpha', alpha1);

% Comparison between groups for Rose odor
[h2 ,pval_group_rose] = ttest2(normalMI(:, 2), adMI(:, 2), 'Alpha', alpha2);

disp('Statistical Significance (p-values):');fprintf('\n');

disp(['Group Comparison - Lemon Odor: p = ', num2str(pval_group_lemon)]);
if h1 == 1
    disp('Lemon trials are significantly different between AD and normal groups');fprintf('\n')
else
    disp('not significantly different');fprintf('\n');
end

disp(['Group Comparison - Rose Odor: p = ', num2str(pval_group_rose)]);
if h2 == 1
    disp('Rose trials are significantly different between AD and normal groups');fprintf('\n')
else
    disp('not significantly different');fprintf('\n');
end




%% End of Bonus 2 =========================================================================================================================================

%% Functions

function plotHeatmap(plvData)
    imagesc(plvData, [0, 1]);
    colorbar;
    xticks(1:4);
    yticks(1:4);
    xticklabels({'Fp1', 'Fz', 'Cz', 'Pz'});
    yticklabels({'Fp1', 'Fz', 'Cz', 'Pz'});
    [x, y] = meshgrid(1:size(plvData, 2), 1:size(plvData, 1));
    text(x(:), y(:), num2str(plvData(:), '%.2f'), 'Color', 'w', 'HorizontalAlignment', 'center');
end


function plv = calculatePLV(signal1, signal2, fs, frequencyRange)
    % signal1: Time series data for channel 1
    % signal2: Time series data for channel 2
    % fs: Sampling frequency (in Hz)
    % frequencyRange(Hz):Frequency range of interest[lowerBound,upperBound]

    % Extract frequency range of interest using a bandpass filter
    filteredSignal1 = bandpass(signal1, frequencyRange, fs);
    filteredSignal2 = bandpass(signal2, frequencyRange, fs);
    % Compute the Hilbert transform to obtain instantaneous phases
    phase1 = angle(hilbert(filteredSignal1));
    phase2 = angle(hilbert(filteredSignal2));
    % Compute phase differences
    phaseDiff = phase2 - phase1;
    % Compute PLV
    plv = abs(mean(exp(1i * phaseDiff)));
end


function out = get_mi(lf_phase,hf_env,nbins)
% out = get_mi(lf_phase,hf_env,nbins) OR
% out = get_mi(lf_phase,hf_env,nbins,nsurr,randType)
%
% This function calls "calc_mi" to calculate MI values 
% Inputs: 
% - lf_phase/ hf_env are TxN/ TxM sized arrays with phase and amplitude information, respectivley.
%   The first two dimensions are "time-freq". Concatenating trials should 
%   be done along the first dimesions. Alternatively, the matrices can be extended into the 
%   third dimension (for example, for single trials analysis)
% - nbins: number of phase bins (eg 18)
% - nsurr (optional): number of randomizations 
% - surrType (optional): either 1,2,3; specifes type of randomization
%       - 1 = scramble phase: (added in later release) phase of low frequency
%       signal is scrambled by taking the inverse FFT. This opertion
%       preserves the power spectrum of the signal.
%       - 2 = timesplice: randomly partitions signal into 2 slices, and
%       rearranges the order of these slices
%       - 3 = randtrial: random permutation of trials of the amplitude signal
%
% Outputs:
% - out: structure with five fields
%        MI: array of MI values, of size=[N,M,# of realizations]
%        MIp: associated p-value (=nan if nargs<4)
%        MIsurrmi: average MI across surroagte MI values (=nan if nargs<4)
%        mean_amps: matrix one dimension higher than out.MI, where the final dimension
%                   corresponds to the average amplitude per phase bin. Useful for
%                   calculating preferred phase
%        ninds: number of instantaneous observations that were used to
%               calculate the MI

% Copyright 2014, Benjamin Voloh
% Distributed under a GNU GENERAL PUBLIC LICENSE

%CHECKS
%check is randomization stats should be retrieved
if nargin<4
    doSurr=false; 
else
    doSurr=true; 
    
    %check what kind of randomizations should be made
    if nargin<5
        error('must specify type of randomization for surrogate stats');
    end

end


%check if amp/ phase matrices have proper dimensions
if ndims(lf_phase) ~= ndims(hf_env) && ...
        size(lf_phase,1) ~= size(hf_env,1) && ... 
        size(lf_phase,3) ~= size(hf_env,3)
    error('mismatched dimensions of phase and ampltiude matrices'); 
end

%check how many dimensions data has, will determine what functions
%subsequently called
if ndims(lf_phase)>2; highdim=true; else highdim=false; end %#ok



%make phase bins
[bins, centers] = make_phase_bins(nbins);


%initialize output: dimensions are "phase freq-amp freq-trial dim
dimtr=3; %dimension for trial
out.MI         = nan(size(lf_phase,2),size(hf_env,2),size(hf_env,dimtr));
out.MIp        = nan(size(lf_phase,2),size(hf_env,2),size(hf_env,dimtr));
out.MIsurrmi   = nan(size(lf_phase,2),size(hf_env,2),size(hf_env,dimtr));

% calculate observed MI
disp('calcualting observed MI...')
[mean_amps, ninds] = wrap_get_amps(lf_phase, hf_env, bins, highdim); %size(mean_amps)= #phase freq, #ampfreq, nbins
uniform = 1/length(centers) .* ones(size(mean_amps)); %calculate here instead of in calc_mi for efficiency
out.MI = calc_mi(centers, mean_amps, uniform);

% calculate surrogate stats
if doSurr
    disp('permutation test...')
    surr_mi = zeros([size(out.MI), nsurr]);

    fn=0;
    for s=1:nsurr
        %updatecounter(s,[1 nsurr],'surrogate run: ')
        disp(['surrogate run # ' num2str(s)])
        
        %randomize signal
        rot_hf_env = randomize_signal(hf_env,randType,highdim);

        [surr_amps, ~] = wrap_get_amps(lf_phase, rot_hf_env, bins, highdim);

        %use linear index to flexibly account for different possible
        %dimensions
        st=fn+1;
        fn=fn+numel(surr_amps)./length(centers);
        
        surr_mi(st:fn) = calc_mi(centers, surr_amps, uniform); %update last
    end
    
    %assign output
    dimsurr=ndims(surr_mi);

    out.MIp       = ( sum(bsxfun(@gt, surr_mi, out.MI),dimsurr)+1 )./(nsurr+1);   
    out.MIsurr    = mean(surr_mi,dimsurr);   
end

%assign outputs
out.mean_amps = mean_amps;  %useful halfstep, can calculate preferred phase later                                                   
out.ninds     = ninds;

disp('...finished MI for one channel pair')
end

function MI = calc_mi(centers, mean_amps, Q)
% calcuates MI based on the KL distance (Tort et al 2010)
% [MI] = calc_mi(centers, mean_amps, uniform)
%
% inputs: 
%   - centres: phase bin centres
%   - mean_amps: average amplitude per phase bin
%   - Q: hypothesized (eg uniform) distribution
% outputs:
%   - MI: normalized KL distance

% Copyright 2014, Benjamin Voloh
% Distributed under a GNU GENERAL PUBLIC LICENSE

nd=ndims(mean_amps);

% normalize mean_amps to make a probability-like distribution
P = bsxfun( @rdivide, mean_amps, nansum(mean_amps,nd) ); 

% MI=normalized KL distance
MI = KL_distance2(P,Q) ./ log10(length(centers));
end

function D = KL_distance2(P, Q)

% if size(P,1) ~= size(Q,1)
%     D = 0;
%     return
% end

D = nansum( P .* log10(P./Q) , ndims(P) );
end

function [mean_amps, ninds]=wrap_get_amps(lf_phase, hf_env, bins, highdim)
% extracts mean amplitde per phase bins
% [mean_amps, ninds]=wrap_get_amps(lf_phase, hf_env, bins, highdim)
%
% Inputs:
%       - lf_phase: TxN array of phases, where T=# of time points and 
%                   N= # of frequencies
%       - hf_env: TxM array of amplitude envelopes, where 
%                 T=# of time points and M= # of frequencies
%       - bins: 2xN array of phase bins, where N= # of phase bins
%       - highdim= rtue or false, determines which functions are called
%
% Outputs:
%       - mean_amps: NxM array of average amplitude in each pahse bin
%       - ninds: number of observations used for each entry of mean_amps

% Copyright 2014, Benjamin Voloh
% Distributed under a GNU GENERAL PUBLIC LICENSE

if highdim
    [mean_amps, ninds] = get_amps3(lf_phase, hf_env, bins);
else
    [mean_amps, ninds] = get_amps2(lf_phase, hf_env, bins);
end
end

function [mean_amps, ninds] = get_amps2(lf_phase, hf_env, bins) 

nbins = length(bins);
mean_amps = zeros(size(lf_phase,2), size(hf_env,2), nbins); %dim=n_phase,n_amp,n_bins
ninds = mean_amps;
%
for i=1:nbins
    ind = sparse( lf_phase >= bins(1,i) & lf_phase < bins(2,i) )';
    notnan = ~isnan(hf_env);

    ninds(:,:,i) = ind*double(notnan);
    hf_env(~notnan)=0;
    mean_amps(:,:,i) = ( ind * hf_env )./ ninds(:,:,i);
end

%set mean amp to 0 if we found no instance of phase in a particular bin
mean_amps(ninds==0)=0;
end

function [mean_amps, ninds] = get_amps3(lf_phase, hf_env, bins)
nbins = length(bins);
mean_amps = zeros(size(lf_phase,2), size(hf_env,2), size(hf_env,3), nbins); %dim=n_phase,n_amp,n_bins
ninds = mean_amps;

notnan = ~isnan(hf_env);
hf_env(~notnan)=0;
for ii=1:nbins
    ind = ( lf_phase >= bins(1,ii) & lf_phase < bins(2,ii) );
    for nT=1:size(hf_env,3)
        ninds(:,:,nT,ii) = ind(:,:,nT)'*double(notnan(:,:,nT));
        mean_amps(:,:,nT,ii) = ( ind(:,:,nT)'*hf_env(:,:,nT) )./ninds(:,:,nT,ii);
    end
end

%set mean amp to 0 if we found no instance of phase in a particular bin
mean_amps(ninds==0)=0;
end

function [bins, centers] = make_phase_bins(nbins)
% creates equally sized phase bins
%
% Copyright 2011, Taufik Valiante
% Distributed under a GNU GENERAL PUBLIC LICENSE

bins = zeros(2,nbins+1);

b = pi/nbins + (0:nbins/2-1)*pi/(nbins/2);

bins(:,1) = [-pi -b(end)];
bins(:,nbins+1) = [b(end) pi];
bins(:,nbins/2+1) = [-pi/nbins pi/nbins];
bins(:,nbins/2+2:end-1) = [b(1:end-1)' b(2:end)']';
bins(:,2:nbins/2) = flipud(fliplr(-bins(:,nbins/2+2:end-1)));

centers = mean(bins,1);  
end

