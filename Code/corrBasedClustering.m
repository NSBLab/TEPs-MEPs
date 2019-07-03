function [avrgSilhScore, sdSilhScore, chanGroupsInd] = corrBasedClustering(pathIn, Data, tepCond, preStim, winLength, Pr)

% 1 - Finds the maximum correlation between group-averge TEPs recorded at each pair of EEG channels
% 2 - Generates the matrix of correlation values before clustering/reordering the channels
% 3 - Generates the matrix of correlation values after clustering/reordering the channels
% 4 - Makes the following  Plots:
% Fig1 - correlation matrix of channel pairs before clustering/reordering the channels
% Fig2 - correlation matrix of channel pairs after clustering/reordering the channels
% Fig3 - Dendrogram Plots Showing the Electrodes Within Each Cluster Using the Optimal Order
% Fig4 - Display the Electrodes of Each cluster on scalp maps

load([pathIn,Data,'.mat']);
cond = strcmp(condition,tepCond);
dataToExamin = meanSubject{cond}(:, (preStim:preStim+winLength-1));
  
for j = 1:length(eeglabChans)
    tepj=[];
    tepj = dataToExamin(j,:);
    
    for jj = 1:length(eeglabChans)
        tepjj=[];
        tepjj = dataToExamin(jj,:);        
        rho(j,jj) = corr(tepj', tepjj','type','Spearman', 'tail', 'right');
      
    end
end

% Create hierarchical cluster tree and reorder the channels by getting optimal dendrogram reordering
[links,order,~] = BF_ClusterReorder(rho,'Euclidean','average');

% Threshold of clustering
Tr = double(Pr*max(links(:,3)));

% Specify the Indices of the Channels in Each Cluster
c = cluster(links, 'Cutoff', Tr, 'criterion', 'distance');
ct = crosstab(c(:,1),eeglabChans)';
clusterN = size(ct,2);

for i = 1:clusterN
    chanGroupsInd{i} = find(ct(:,i));
end

% Test clusters quality
% figure
% [silh,h] = silhouette(rho, c);

[silh] = silhouette(rho,c);
avrgSilhScore = mean(silh);
sdSilhScore = std(silh);

%%
%---------------------------------------------------------------------------------------------------
% Fig1 - Distance matrix of channel pairs before clustering/reordering the channels
%---------------------------------------------------------------------------------------------------

title1 = 'Channel-Channel Correlation';

% close all;
Fig1 = figure('color', 'w');
imagesc(rho);
set(gca, 'YTick', 1:length(eeglabChans)); 
set(gca, 'YTickLabel', eeglabChans,'FontSize', 20); 
set(gca, 'XTick', 1:length(eeglabChans)); 
set(gca, 'XTickLabel', eeglabChans,'FontSize', 20); 
set(gca,'XTickLabelRotation',90)
box off
colormap(jet(10)); 
title(title1, 'FontSize', 20);
colorbar;

%%
%---------------------------------------------------------------------------------------------------
% Fig2 - Distance matrix of channel pairs after clustering/reordering the channels
%---------------------------------------------------------------------------------------------------

title2 = 'Reordered Channel-Channel Correlations';

% close all;
Fig2 = figure;
imagesc(rho(order,order));
set(gca, 'YTick', 1:length(eeglabChans)); 
set(gca, 'YTickLabel', eeglabChans(order),'FontSize', 20); 
set(gca, 'XTick', 1:length(eeglabChans)); 
set(gca, 'XTickLabel', eeglabChans(order),'FontSize', 20); 
set(gca,'XTickLabelRotation',90)
box off
set(gca, 'color', 'w'); 
colormap(jet(10)); 
title(title2, 'FontSize', 20);
colorbar;

%%
%---------------------------------------------------------------------------------------------------
% Fig3 - Dendrogram Plots Showing the Electrodes Within Each Cluster Using the Optimal Order
%---------------------------------------------------------------------------------------------------

title3 = 'Reordered Channels-Clusters';

% close all;
Fig3 = figure ('color', 'w');
dendrogram(links,0,'reorder',order, 'ColorThreshold',Tr,'Orientation','right', 'Labels',eeglabChans);
title(title3, 'FontSize', 20);

%%
%---------------------------------------------------------------------------------------------------
% Fig4 - Display the Electrodes of Each cluster on scalp maps
%---------------------------------------------------------------------------------------------------

% title4 = 'Electrodes within Each Cluster';
% colorr = {[0 0.50 0.36];[0 0.14 0.53];'k';'r'};
% 
% % close all;
% Fig4 = figure;
% set(gca,'color','w')
% grandAverage = [];
% grandAverage.label = eeglabChans;
% grandAverage.time = [1];
% grandAverage.dimord = 'chan_time';
% grandAverage.avg = zeros(length(eeglabChans),1);
% cfg = [];
% cfg.parameter = 'avg';
% cfg.zlim = [-1,0];
% cfg.comment = 'no';
% cfg.layout = 'quickcap64.mat';
% cfg.highlight = 'Labels';
% cfg.highlightchannel = chanGroupsInd;
% cfg.highlightcolor   = colorr;
% cfg.highlightsymbol = {'.';'.';'.';'.';'.'};
% cfg.highlightsize = [{60};{60};{60};{60};{60}];
% cfg.colormap = 'hot';
% title(title4, 'FontSize', 14);
% ft_topoplotER(cfg, grandAverage);
end