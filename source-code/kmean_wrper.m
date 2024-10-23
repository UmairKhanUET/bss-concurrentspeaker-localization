%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2024, Muhammad Umair Khan (umair.khan.uet59@gmail.com)
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [plotHandle, spotHandle, textBoxHandle, globalSpeakerEstimate] = ...
    kmean_wrper(cluster_size, angle_estimates, num_speakers)

% Create variables to be returned to the caller.
plotHandle = [];
spotHandle = [];
textBoxHandle = [];

hold on;

all_data = [];
for i = 1:1:length(angle_estimates)
    all_data = [all_data; ...
                angle_estimates(i).ch1_locus; angle_estimates(i).ch1p_locus; ...
                angle_estimates(i).ch2_locus; angle_estimates(i).ch2p_locus];
end

IDX = kmeans(all_data, num_speakers);

total_clusters = max(IDX);
clusters_list = zeros(1, total_clusters);

% Cluster the points.
for i = 1:1:total_clusters
    clusters_list(i) = sum(IDX == i);
end

% Sort the clusters in descending order of their sizes.
sorted_clusters_list = sort(clusters_list, 'descend');

if (total_clusters >= num_speakers)
    firstNClusters = num_speakers;
else
    firstNClusters = total_clusters;
end

biggest_clusters = zeros(1, firstNClusters);

% Here we pick first N biggest clusters out of the sorted cluster list. If
% the list contains clusters of equal sizes, we need to pick up the most
% dense ones. If the density of clusters is also same, then all of them
% carry equal probability, so, just pick anyone of them.

% Initialize variables before we start.
i = 1;

while (i <= firstNClusters)
    % Find how many same sized clusters are there and make sure they don't 
    % span more than the length of firstNClusters.
    % Save the indices of the clusters of interest.
    clusterIndices = find(clusters_list == sorted_clusters_list(i));
    numSameSizeClusters = length(clusterIndices);

    if (i + numSameSizeClusters - 1) > firstNClusters
        % Here is where we need to pick clusters based on cluster
        % densitites. 
        % First check how many clusters we need to pick in this
        % iteration. This is equal to the number of empty entries in 
        % biggest_clusters array.
        clustersToPick = length(find(biggest_clusters == 0));
        
        % Create cluster width array.
        clusterWidth = zeros(numSameSizeClusters, 2);

        % Lets find the cluster spans now.
        for j = 1:1:numSameSizeClusters
            cluster = all_data(IDX == clusterIndices(j), :);
            clusterMean = mean(cluster);
            temp = sum(abs(cluster - clusterMean));
            
            % Convert to a scalar.
            clusterWidth(j, 1) = sqrt(temp(1,1)^2 + temp(1,2)^2);

            % Save the index as well.
            clusterWidth(j, 2) = clusterIndices(j);
        end
        
        % Lets sort the clusterWidth array in ascending order and pick
        % first clusterToPick number of clusters.
        [~, idx] = sort(clusterWidth(: ,1), 'ascend');
        clusterWidth = clusterWidth(idx, :);
        
        % Lets pick the cluster indices now, finally.
        biggest_clusters(i : i + clustersToPick - 1) =  ...
                         clusterWidth(1 : clustersToPick, 2);
    else
        clustersToPick = numSameSizeClusters;
        biggest_clusters(i : i + clustersToPick - 1) =  clusterIndices;
    end

    % Increment loop counter
    i = i + numSameSizeClusters;
end

clusterMean = zeros(firstNClusters, 2);
cluster_data = [];
idx = [];

i = 1;

while (i <= firstNClusters)
    indx = biggest_clusters(i);

    % Check if there are more than one indices returned.
    for j = 1:1:length(indx)
        cluster_data = [cluster_data; all_data(IDX == indx(j), :)];
        idx = [idx, i * ones(1, length(find(IDX == indx(j))))];
        clusterMean(i,:) = mean(all_data(IDX == indx(j), :));

        i = i + 1;
        
        if (i > firstNClusters)
            break;
        end
    end
end

if (firstNClusters > 0)
    plotHandle = PlotClusterinResult(cluster_data, idx);
end

% Lets demarkate the estimated speaker locations.
r = 0.01;

% Draw big RED spots at the estimated speaker locations.
for i = 1:1:firstNClusters
    center = clusterMean(i, :);
    pos = [center-r 2*r 2*r];
    spotHandle = [spotHandle, rectangle('Position',pos, ...
                      'Curvature',[1 1], 'FaceColor', 'red', 'Edgecolor',...
                      'none')];
end

globalSpeakerEstimate = zeros(1, firstNClusters);

% Let put the global estimates in a text box.
for i = 1:1:firstNClusters
    globalSpeakerEstimate(i) = atan2(clusterMean(i,2), clusterMean(i,1))*180/pi;
    
    % Transform the estimate as per our reference axis.
    globalSpeakerEstimate(i) = 360 - (globalSpeakerEstimate(i) + 30);
    
    % Lets restrict it between 0 to 360 degress.
    globalSpeakerEstimate(i) = mod(globalSpeakerEstimate(i) + 360, 360);

    textBoxHandle = [textBoxHandle, annotation('textbox', ...
                     [0.47, 0.15 + i*0.05, 0.1, 0.1], 'String', ...
                     num2str(globalSpeakerEstimate(i), '%7.3f'))];
end
