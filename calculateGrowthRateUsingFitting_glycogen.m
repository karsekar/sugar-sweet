% calculateGrowthRateUsingFitting_glycogen

%adapted from Jen Nguyen's code by Karthik on 29.1.2019.

%%
function [tracks] = calculateGrowthRateUsingFitting_glycogen(volumes,isDrop,trackNum,dt_sec,cfp,yfp,frames)

% input data:
%        volumes     =  calculated va_vals (cubic um)
%        isDrop      =  1 marks a birth event, 0 is normal growth

%filter parameters for the smoothing line
b = [0.2 0.2 0.2 0.2 0.2];
a = 1;

% 0. define dt
%dt = 120; % timestep in seconds

%separate the different tracks of cells
tracks = cell(1,max(trackNum));


figure;
%populate different elements into each track
for i=1:length(tracks)
    
    hold off;
    
    %get the indices
    tracks{i}.indices=find(trackNum==i);
    
    %get the volumes
    tracks{i}.volumes=volumes(tracks{i}.indices);
    
    %get the division occurrence indices
    tracks{i}.divisionOccurrences=find(isDrop(tracks{i}.indices));
    
    %get the fluorescence information
    tracks{i}.cfp=cfp(tracks{i}.indices);
    tracks{i}.yfp=yfp(tracks{i}.indices);
    
    %if there is a division occurrence, find up until that point. Otherwise,
    %take the entire range
    if(isempty(tracks{i}.divisionOccurrences))
        tracks{i}.timeRange=1:length(tracks{i}.indices);
        tracks{i}.frames=frames(tracks{i}.indices);
    elseif(tracks{i}.divisionOccurrences(1) > 3) %make sure the first division occurs far enough in
        tracks{i}.timeRange=1:(tracks{i}.divisionOccurrences(1)-2);
        tempframes = frames(tracks{i}.indices);
        tracks{i}.frames=tempframes(1:(tracks{i}.divisionOccurrences(1)-2));
    else
        tracks{i}.timeRange=1:length(tracks{i}.indices);
        tracks{i}.frames=frames(tracks{i}.indices);
    end
    
    %take the transpose
    tracks{i}.frames=tracks{i}.frames';
        
    %get the volumes of interest
    tracks{i}.volumesUntilDivision=tracks{i}.volumes(tracks{i}.timeRange);
    
    %transform the time to the correct units (hours)
    tracks{i}.time=tracks{i}.timeRange'*dt_sec/3600;
    
    %perform linear fit on the data
    tracks{i}.P=polyfit(tracks{i}.time,tracks{i}.volumesUntilDivision,1);
    tracks{i}.volumeFits=tracks{i}.P(1)*tracks{i}.time+tracks{i}.P(2);
    
    %filter the data for smoothing
    %tracks{i}.filtered=filter(b,a,tracks{i}.volumesUntilDivision);
    tracks{i}.filtered=smooth(tracks{i}.time,tracks{i}.volumesUntilDivision,0.9,'rloess');
    

    %plot to see what it looks like
    subplot(1,2,1);
    plot(tracks{i}.time,tracks{i}.volumesUntilDivision,'b.','MarkerSize',10);
    hold on;
    %plot(tracks{i}.time,tracks{i}.volumeFits,'-r');
    
    plot(tracks{i}.time,tracks{i}.filtered,'-r');
    xlabel('Time (h)');
    ylabel('Volume');
    title(['Track ' num2str(i)]);
    hold off;
    
    %growth rate and initial cell size is from the fit
    tracks{i}.growthRates=diff(tracks{i}.filtered);
    
    subplot(1,2,2);
    plot(tracks{i}.time(1:end-1),tracks{i}.growthRates);
    xlabel('Time (h)');
    ylabel('Extension rate');
    hold off;
    
end


end







