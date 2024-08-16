%% 5.1.2 : PVL of MCI subjects
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



%% 5.1.2 : box plot and fitting Guassian
clear; clc;
mciPLV = load('mciPLV.mat');
mciPLV = mciPLV.mciPLV;
% 
% adPLV = load('adPLV.mat');
% adPLV = adPLV.adPLV;

% mciPLV = load('normalPLV.mat');
% mciPLV = mciPLV.normalPLV;

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
   
%% phase difference

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

%% heat map

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


%% 5.1.2 : p-value and significance testing on the 3 states 

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

%% function

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

