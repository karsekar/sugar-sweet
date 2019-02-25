%% quantifyGrowth with Fitting

% Goal: quantify and visualize growth rates of WT and mutants for glycogen
%       utilization under a pulsling nutrient environment.
%
%       input data is data structure produced from glycogen_analysis.m

%       this script contains methods of measuring growth rate:
%
%           section 1. instantaneous rate of volumetric doubling
%           section 2. total area occupied, normalized per cell and initial value
%           section 3. dA/dt, where A is normalized by the number of cells
%           - for specifics, see strategy at the start of each section


% last update: Jen, 2019 Jan 25
% commit: add section to calculate and plot dA/dt

% ok let's go!

%% section 1. instantaneous rate of volumetric doubling


% Strategy:
%                0. initialize experiment parameters and data
%                1. assemble data matrix
%                2. isolate volume (Va), drop, and track number for growth rate calculations
%                3. calculate growth rate
%                4. truncate data to non-erroneous timestamps (e.g. bubbles) 
%                5. bin growth rate into time bins based on timestamp
%                6. isolate selected specific growth rate and remove nans from data analysis
%                7. isolate YFP and CFP intensities
%                8. convert intensities to (+) or (-) fluorophore
%                       - for details on threshold determination,
%                         see scrap_pile section "YFP and CFP thresholds"
%                9. test threshold, throw error if any data points ID as both
%               10. separate growth rates by fluorophore
%               11. bin growth rate by time
%               12. calculate stats for each bin, including mean growth rate
%               13. plot growth rate over time



clear
clc
close all;

% 0. initialize data
%xy = 2;
xy_start = 1;
xy_end = 16;
dt_min = 2;
%dt_min = 30; % reduced frequency dataset

date = '2018-11-23';

%cd(strcat('/Users/jen/Documents/StockerLab/Data/glycogen/',date))
%load(strcat('glycogen-',date,'-allXYs-jiggle-0p5.mat'),'D5');
load("glycogen-2018-11-23-allXYs-jiggle-0p5.mat",'D5');

% 0. define growth rates of interest
specificGrowthRate = 'log2';


% 0. define time binning parameters
specificBinning = 60; % in minutes
binsPerHour = 60/specificBinning;





% 1. assemble data matrix
glycogen_data = buildDM_glycogen(D5, xy_start, xy_end, dt_min);
clear xy_start xy_end


% 2. isolate volume (Va), drop, and track number for growth rate calculations
volumes = glycogen_data(:,5);        % col 4 = calculated va_vals (cubic um)
isDrop = glycogen_data(:,3);         % col 2 = isDrop, 1 marks a birth event
trackNum = glycogen_data(:,12);      % col 12 = track number (not ID from particle tracking)

%get the grouping information

cfp = glycogen_data(:,13);         % col 13 = mean CFP intensity
yfp = glycogen_data(:,14); 
frames = glycogen_data(:,9);

% 3. calculate summary info
dt_sec = dt_min * 60;

%% calculate the summary info
tracks = calculateGrowthRateUsingFitting_glycogen(volumes,isDrop,trackNum,dt_sec,cfp,yfp,frames);

%save the file in case we need to run again and don't want to redo this
%section
save('tracks.mat','tracks');


%% summarizing the data

close all;
load('tracks.mat');

%find where growth rates are negative (shouldn't be cells)
growthRatesAcrossTracks=cellfun(@(x) mean(x.growthRates), tracks);
growthRateThreshold=0.0;
growthRatesAboveThreshold=find(growthRatesAcrossTracks>growthRateThreshold);

figure;
histogram(growthRatesAcrossTracks(growthRatesAboveThreshold));

%find where the tracks are long enough
minTime = 5; %units is hr
maxDurationAcrossTracks=cellfun(@(x) x.time(end), tracks);
tracksExceedingMinTime=find(maxDurationAcrossTracks > minTime);

% 0. define fluorescence intensity threshold
cfp_threshold = 103.4; %103.4;
yfp_threshold = 111;

%get the mean fluorescence datas
cfpAcrossTracks=cellfun(@(x) mean(x.cfp), tracks);
yfpAcrossTracks=cellfun(@(x) mean(x.yfp), tracks);

%find the tracks where the threshold is exceeded
cfpCellIndices = find((cfpAcrossTracks > cfp_threshold)&(cfpAcrossTracks < 350)); %based on histogram. looks like an outlier above 350. let's exclude that
yfpCellIndices = find(yfpAcrossTracks > yfp_threshold);

%only consider where growth rates are positive and times exceed the min
cfpCellIndices = intersect(intersect(cfpCellIndices,growthRatesAboveThreshold),tracksExceedingMinTime);
yfpCellIndices = intersect(intersect(yfpCellIndices,growthRatesAboveThreshold),tracksExceedingMinTime);

%how many do we get with each gating?
cfpNum=length(cfpCellIndices);
yfpNum=length(yfpCellIndices);
intersecting=intersect(cfpCellIndices,yfpCellIndices)

%look at the histogram of the resulting values
figure;
subplot(2,1,1);
histogram(cfpAcrossTracks(cfpCellIndices),'BinWidth',5);
title('CFP distribution');
subplot(2,1,2);
histogram(yfpAcrossTracks(yfpCellIndices),'BinWidth',5);
title('YFP distribution');

%now for each group look at the growth rate over time
%let's first do CFP

figure;
for i=1:cfpNum
    currentTrack=cfpCellIndices(i);
    plot(tracks{currentTrack}.time(1:end-1),tracks{currentTrack}.growthRates,'-c');
    hold on;
end
title('CFP Extension Rates');
xlabel('Time (h)');
ylim([-0.05 0.1]);

%then YFP

figure;
for i=1:yfpNum
    currentTrack=yfpCellIndices(i);
    plot(tracks{currentTrack}.time(1:end-1),tracks{currentTrack}.growthRates,'-m');
    hold on;
end
title('YFP Extension Rates');
xlabel('Time (h)');
ylim([-0.05 0.1]);


% clear isDrop volumes trackNum dt_min
% 
% 
% % 4. truncate data to non-erroneous timestamps (e.g. bubbles) 
% maxTime = 8; % in hours
% frame = glycogen_data(:,9);      % col 9 = frame in image sequence
% timeInSeconds = frame * dt_sec;  % frame = is consequetive images in analysis
% timeInHours = timeInSeconds/3600;
% 
% if maxTime > 0
%     glycogenData_bubbleTrimmed = glycogen_data(timeInHours <= maxTime,:);
%     growthRates_bubbleTrimmed = growthRates(timeInHours <= maxTime,:);
% else
%     glycogenData_bubbleTrimmed = glycogen_data;
%     growthRates_bubbleTrimmed = growthRates;
% end
% clear maxTime frame timeInSeconds timeInHours
% 
% 
% % 5. bin growth rate into time bins based on timestamp
% frame = glycogenData_bubbleTrimmed(:,9);      % col 9 = frame in image sequence
% timeInSeconds = frame * dt_sec;    
% timeInHours = timeInSeconds/3600;
% bins = ceil(timeInHours*binsPerHour);
% clear timeInSeconds frame
% 
% 
% % 6. isolate selected specific growth rate and remove nans from data analysis
% specificColumn = 3;
% growthRate_log2 = growthRates_bubbleTrimmed(:,specificColumn);
% 
% growthRt_noNaNs = growthRate_log2(~isnan(growthRate_log2),:);
% bins_noNaNs = bins(~isnan(growthRate_log2),:);
% glycogenData_noNaNs = glycogenData_bubbleTrimmed(~isnan(growthRate_log2),:);
% 
% 
% % 7. isolate YFP and CFP intensities
% cfp = glycogenData_noNaNs(:,13);         % col 13 = mean CFP intensity
% yfp = glycogenData_noNaNs(:,14);         % col 14 = mean YFP intensity
% 
% 
% % 8. convert intensities to (+) or (-) fluorophore
% isCFP = cfp > threshold;
% isYFP = yfp > threshold;
% 
% 
% 
% % 9. test threshold, throw error if any data points ID as both
% isBoth = isCFP+isYFP;
% if sum(isBoth == 2) > 0
%     error('threshold fail! some cells positive for both fluorophores')
% end
% 
% % this step doesn't matter as long as threshold doesn't allow double-positives
% growthRt_final = growthRt_noNaNs(isBoth < 2);
% glycogenData_final = glycogenData_noNaNs(isBoth < 2,:);
% bins_final = bins_noNaNs(isBoth < 2);
% isYFP_final = isYFP(isBoth < 2);
% isCFP_final = isCFP(isBoth < 2);
% 
% 
% 
% % 10. separate growth rates by fluorophore
% 
% % growthRt_yfp = growthRt_final(isYFP_final == 1);
% % growthRt_cfp = growthRt_final(isCFP_final == 1);
% 
% % above method only takes growth rates if the cell is actively labeled.
% % however, all the growth rates should be taken into account! instead, pull
% % out IDs that are identified with YFP or CFP and use these to 
% trackNum = glycogenData_final(:,12); % col 12 = track num
% IDs_yfp = unique(trackNum(isYFP_final == 1));
% IDs_cfp = unique(trackNum(isCFP_final == 1));
% 
% IDs_labeled = [IDs_yfp; IDs_cfp];
% 
% growthRt_yfp = [];
% growthRt_cfp = [];
% bins_yfp = [];
% bins_cfp = [];
% for id = 1:length(IDs_labeled)
%     
%     currentID = IDs_labeled(id);
%     currentGR = growthRt_final(trackNum == currentID);
%     currentBins = bins_final(trackNum == currentID);
%     
%     if ismember(currentID,IDs_yfp) == 1 % if ID belongs to YFP+
%         growthRt_yfp = [growthRt_yfp; currentGR];
%         bins_yfp = [bins_yfp; currentBins];
%     else % else
%         growthRt_cfp = [growthRt_cfp; currentGR];
%         bins_cfp = [bins_cfp; currentBins];
%     end
%     
% end
% clear currentID currentGR currentBins
% 
% 
% % 11. bin growth rate by time
% binned_yfp = accumarray(bins_yfp,growthRt_yfp,[],@(x) {x});
% binned_cfp = accumarray(bins_cfp,growthRt_cfp,[],@(x) {x});
% 
% 
% % 12. calculate mean, standard dev, counts, and standard error
% y_bin_means = cellfun(@mean,binned_yfp);
% y_bin_stds = cellfun(@std,binned_yfp);
% y_bin_counts = cellfun(@length,binned_yfp);
% y_bin_sems = y_bin_stds./sqrt(y_bin_counts);
% 
% c_bin_means = cellfun(@mean,binned_cfp);
% c_bin_stds = cellfun(@std,binned_cfp);
% c_bin_counts = cellfun(@length,binned_cfp);
% c_bin_sems = c_bin_stds./sqrt(c_bin_counts);
% 
% 
% 
% % 13. plot growth rate over time
% palette = {'DodgerBlue','GoldenRod'};
% 
% yfp_color = 'y'; %rgb(palette(2));
% cfp_color = 'c'; %rgb(palette(1));
% xmark = 'o';
% 
% figure(1)
% errorbar((1:length(y_bin_means))/binsPerHour,y_bin_means,y_bin_sems,'Color',yfp_color)
% hold on
% errorbar((1:length(c_bin_means))/binsPerHour,c_bin_means,c_bin_sems,'Color',cfp_color)
% hold on
% grid on
% axis([0,8.5,-1,1])
% xlabel('Time (hr)')
% ylabel('Growth rate')
% title(strcat(date,': (',specificGrowthRate,')'))
% legend('YFP WT', 'CFP mutant')


