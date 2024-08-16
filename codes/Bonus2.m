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

%%
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

%%
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

