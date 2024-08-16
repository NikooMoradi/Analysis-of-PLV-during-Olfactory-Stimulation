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
    
    epochData = normalData.normal(person).epoch; % 4 x 600 x NumTrials matrix 
     
    odorData = normalData.normal(person).odor; % NumTrials x 1 binary array
    
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


%% 4.3
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
