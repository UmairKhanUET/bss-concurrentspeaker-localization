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

% Colors and formats for plotting graphs dynamically.
Colors = 'bkrgycmw';
Styles = '+xo*^v><ph';

% This is separation between two microphones on the diagonal.
L = 0.0926;
ParentFolder = 'respeaker_v2';
Date = '10-03-2019';
Profile = 'dynamic';
audio_name = '270_210_dyn';
RootPath = strcat('../', ParentFolder, '/', Date, '/', Profile, '/',...
                  audio_name);
wavFilesFolder = 'wavResults';
TotalAudioChannels = 6;

ClusterSize = 2;
NumSpeakers = 2;
FramesPerIteration = 5;

% Create and populate the algorithm configuration structure.
config = [];
[config] = SetParams(config, L);
[config] = SetExtraParams(config);

if (config.sig_len == 0)
    audioInfo = audioinfo(strcat(RootPath, '/', audio_name, '.wav'));
    config.sig_len = audioInfo.TotalSamples/audioInfo.SampleRate;
end

NumBuffers = floor(config.sig_len*config.fs/config.B);

% If the individual channels aren't extracted, lets extract them.
if (~exist(strcat(RootPath, '/', 'AudioIn1.wav'), 'file') || ...
    ~exist(strcat(RootPath, '/', 'AudioIn2.wav'), 'file') || ...
    ~exist(strcat(RootPath, '/', 'AudioIn3.wav'), 'file') || ...
    ~exist(strcat(RootPath, '/', 'AudioIn4.wav'), 'file') || ...
    ~exist(strcat(RootPath, '/', 'AudioIn5.wav'), 'file') || ...
    ~exist(strcat(RootPath, '/', 'AudioIn6.wav'), 'file'))
       separateChannels(strcat(RootPath, '/'), audio_name);
end

% List all possible microphone pairs in the geometry.
micPairs = [[4, 1]; [6, 3]; [5, 2]]; %...                         % Diagonals
%micPairs =            [[6, 1]; [1, 2]; [2, 3]; [3, 4]; [5, 4]; [6, 5]];  % Sides
%            micPairs = [[6, 2]; [6, 4]; [4, 2]];                            % Triangular diagonals

% Enter direction in which the angle estimates will be drawn on the 
% hexagon. draw.Outward means the angle will be drawn outward of the closed
% shape. Similarly, draw.Inward means the angle will be drawn inward of the
% surface. draw.OutwardInward means it will be drawn both ways.
vectorDirection = [draw.OutwardInward, draw.OutwardInward, draw.OutwardInward, ...
                   draw.Outward, draw.Outward, draw.Outward, draw.Outward, draw.Inward, draw.Inward, ...
                   draw.Outward, draw.Outward, draw.Outward];

% Geometries to plot.
DrawDiagonals = 1;                  % Estimates on diagonals will be drawn.
DrawFaces = 0;                      % Estimates on hexagon faces will be drawn.
DrawTriangularDiagonals = 1;        
PlotIndividualChannelEstimates = 1; % This would plot the DOA estimates on 
                                    % a 2D graph.

% Maximize the figure.
fig = figure(1);
set(fig, 'units','normalized','outerposition',[0 0 1 1]);

% Clear out any remnants.
clf('reset');

% Setup the subplots.

for i = 1:1:2

    subplot(1, 2, i);

    hold on;

    % Set axis equal to size of the re-speaker
    axis([-L L -L L]);

    % Draw the hexagonal boundary of the respeaker.
    hexagon = drawHexagon(L/2, 30);

    % Make both axis of same length.
    pbaspect([1 1 1]);
    
    hold off;
end

% Select the first subplot.
subplot(1, 2, 1);

hold on;

% Create and initialize the vertices array
vertices = struct('Xdata', zeros(1, TotalAudioChannels), 'Ydata', zeros(1, TotalAudioChannels));

% Save vertices of the hexagon. The index denotes the microphone number
% marked on the ReSpeaker periphery.
for i = 1:1:TotalAudioChannels
    vertices(i).Xdata = hexagon.XData((TotalAudioChannels-i)+2);
    vertices(i).Ydata = hexagon.YData((TotalAudioChannels-i)+2);
end

% Draw inscribed triangle if needed.
if DrawTriangularDiagonals == 1
    
    % Draw the triangle
    % Draw line from vertex 6 to 2
    line([vertices(6).Xdata vertices(2).Xdata], [vertices(6).Ydata vertices(2).Ydata]);
    % Draw line from vertex 6 to 4
    line([vertices(6).Xdata vertices(4).Xdata], [vertices(6).Ydata vertices(4).Ydata]);
    % Draw line from vertex 4 to 2
    line([vertices(4).Xdata vertices(2).Xdata], [vertices(4).Ydata vertices(2).Ydata]);

end

% Calculate angles between each pair of points and horizontal line.
lineAngles = zeros(TotalAudioChannels, TotalAudioChannels);

for i = 1:1:TotalAudioChannels
    for j = 1:1:TotalAudioChannels
        if i == j
            lineAngles(i, j) = 0;
        else
            % Following 'if' checks are needed to get rid of signed zero.
            if (round(vertices(i).Ydata-vertices(j).Ydata, 5)) == 0
                y_component = 0;
            else
                y_component = vertices(i).Ydata-vertices(j).Ydata;
            end
            if (round(vertices(i).Xdata-vertices(j).Xdata, 5)) == 0
                x_component = 0;
            else
                x_component = vertices(i).Xdata-vertices(j).Xdata;
            end
            
            % Find length of hypotenuse.
            hyp = sqrt(x_component^2 + y_component^2);
            
            % Calculate angle.
            lineAngles(i, j) = (atan2(y_component, x_component)*180/pi);
        end
    end
end

% Now lets calculate normals to these line angles.
normals = lineAngles + 90;

% Calculate microphone spacing in units of meters.
% Spacing of the microphones on diagonal. Should be 0.0926m.
diagonalMicSpacing = L;

% Spacing of the adjacent microphones. Should be 0.0463m.
adjacentMicSpacing = sqrt((vertices(1).Xdata-vertices(2).Xdata)^2+ ...
                          (vertices(1).Ydata-vertices(2).Ydata)^2);

% Spacing of the triangular diagonal. Should be 0.0801940m.
triangularDiagonalMicSpacing = sqrt((vertices(1).Xdata-vertices(3).Xdata)^2+ ...
                          (vertices(1).Ydata-vertices(3).Ydata)^2);
% DOA's
doa_ch1 = zeros(TotalAudioChannels, TotalAudioChannels, NumBuffers);
doa_ch2 = zeros(TotalAudioChannels, TotalAudioChannels, NumBuffers);

% Create an array of structures to store all angle estimates.
angle_estimates = repmat(struct(), 1, length(micPairs) * NumBuffers);

% Traverse the microphone pairs list.
for i = 1:1:length(micPairs)

    % Select the microphone pair.
    mic1 = micPairs(i, 1);
    mic2 = micPairs(i, 2);
    
    % Store mic pairs in angle_estimates(i).
    angle_estimates(i).micPair = strcat('[', num2str(mic1), ',', num2str(mic2), ']'); 

    % Find spacing between selected microphone pair with distance formula.
    spacing = sqrt((vertices(mic1).Xdata-vertices(mic2).Xdata)^2+ ...
                   (vertices(mic1).Ydata-vertices(mic2).Ydata)^2);

    % Run BSS algo for the selected microphone pair.
    [doa_ch1(mic1, mic2, :), doa_ch2(mic1, mic2, :)]  ...
      = testBSS(RootPath, 'wavResults', spacing, micPairs(i, :));
end

% Create array for handles returned by the plot function. We shall have
% twice the number of mic pairs sized array since each mic pair returns
% two real estimates and two phantom estimates.
pHandles = zeros(1, length(micPairs) * 4);
globalSpkrEstim = struct();

first_image = 1;
estimateIndex = 1;

for j = 0:1:(NumBuffers - FramesPerIteration)

    pHandleIndex = 1;

    % Select subplot.
    subplot(1, 2, 1);

    hold on;

    startIndex = (j) * length(micPairs) + 1;

    for k = 0:1:(length(micPairs) * FramesPerIteration - 1)

        i = mod(k, length(micPairs)) + 1;
        
        if ( i == 1 )
            j = j + 1;
        end

        % Select the microphone pair.
        mic1 = micPairs(i, 1);
        mic2 = micPairs(i, 2);

        % Find index of angle_estimates structure array.
        ae_index = (j - 1) * length(micPairs) + i;

        % Find spacing between selected microphone pair with distance formula.
        spacing = sqrt((vertices(mic1).Xdata-vertices(mic2).Xdata)^2+ ...
                       (vertices(mic1).Ydata-vertices(mic2).Ydata)^2);

        % Plot DOA for channel 1.
        [angle_estimates(ae_index).ch1, angle_estimates(ae_index).ch1p, ...
         angle_estimates(ae_index).ch1_locus, angle_estimates(ae_index).ch1p_locus, ...
         pHandles(pHandleIndex), pHandles(pHandleIndex + 1)] ...
         = drawGeometryDynamic(normals(mic1, mic2), doa_ch1(mic1, mic2, j), spacing, ...
                     (vertices(mic1).Xdata + vertices(mic2).Xdata)/2, ...
                     (vertices(mic1).Ydata + vertices(mic2).Ydata)/2, ...
                      vectorDirection(i), L);

        % Increment the index.
        pHandleIndex = pHandleIndex + 2;

        % Plot DOA for channel 2.
        [angle_estimates(ae_index).ch2, angle_estimates(ae_index).ch2p, ...
         angle_estimates(ae_index).ch2_locus, angle_estimates(ae_index).ch2p_locus, ...
         pHandles(pHandleIndex), pHandles(pHandleIndex + 1)] ...
         = drawGeometryDynamic(normals(mic1, mic2), doa_ch2(mic1, mic2, j), spacing, ...
                     (vertices(mic1).Xdata + vertices(mic2).Xdata)/2, ...
                     (vertices(mic1).Ydata + vertices(mic2).Ydata)/2, ...
                      vectorDirection(i), L);

        % Increment the index.
        pHandleIndex = pHandleIndex + 2;

    end

    endIndex = ae_index;

    % Make pHandles visible.
    set(pHandles, 'Visible', 'On');

    % Add title to the plot.
    title(strcat(num2str(j), '/', num2str(NumBuffers)));

    hold off;

    % Start clustering
    subplot(1, 2, 2);

    hold on;

    % Apply clustering algorithm.
    [clusterHandle, spotHandle, textBoxHandle, ...
        globalSpeakerEstimatesTmp] = dbscan(ClusterSize, 0.025, ...
                           angle_estimates(startIndex : endIndex), ...
                           NumSpeakers);

    % Save global estimates.
    for i = 1:1:length(globalSpeakerEstimatesTmp)
        globalSpkrEstim.(strcat('e', num2str(i)))(estimateIndex) = ...
            globalSpeakerEstimatesTmp(i);
    end

    estimateIndex = estimateIndex + 1;

    if (~isempty(clusterHandle))
        set(clusterHandle, 'Visible', 'On');
    end
    if (~isempty(spotHandle))
        set(spotHandle, 'Visible', 'On');
    end
    if (~isempty(textBoxHandle))
        set(textBoxHandle, 'Visible', 'On');
    end

    % Now, save the snapshot.
    % Make sure that the directory exists.
    if ~exist(strcat(Profile, '/gif/', Date, '/', audio_name), 'dir')
        mkdir(strcat(Profile, '/gif/', Date, '/', audio_name));
    end

    % Capture the plot as an image 
    frame = getframe(fig);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);

    % Write to the GIF File 
    if first_image == 1 
        imwrite(imind, cm, strcat(Profile, '/gif/', Date, '/', audio_name, '/', ...
                audio_name, '.gif'), 'gif', 'Loopcount', inf); 
        first_image = 0;
    else 
        imwrite(imind, cm, strcat(Profile, '/gif/', Date, '/', audio_name, '/', ...
                audio_name, '.gif'), 'gif', 'WriteMode', 'append'); 
    end

    % Make pHandles invisible.
    if (~isempty(pHandles))
        set(pHandles, 'Visible', 'Off');
    end
    if (~isempty(clusterHandle))
        set(clusterHandle, 'Visible', 'Off');
    end
    if (~isempty(spotHandle))
        set(spotHandle, 'Visible', 'Off');
    end
    if (~isempty(textBoxHandle))
        set(textBoxHandle, 'Visible', 'Off');
    end

    hold off;
end

figure(2)

x_axis = 1:1:length(globalSpkrEstim.e1);
x_axis = x_axis * audioInfo.Duration / length(globalSpkrEstim.e1);

hold on;

% Plot global estimates.
for i = 1:1:length(fieldnames(globalSpkrEstim))
    plot(x_axis, globalSpkrEstim.(strcat('e', num2str(i))), ...
         strcat(Colors(mod(i, length(Colors)) + 1), ...
                Styles(mod(i, length(Styles)) + 1)));
end

axis([0 max(x_axis) 0 360]);

tokens = textscan(audio_name, '%s', 'delimiter', '_');
tokens = tokens{1};

for i = 1:1:NumSpeakers
    y = str2double(tokens(i));
    line([0 length(globalSpkrEstim.e1)], [y y], 'Color', 'blue', ...
        'LineWidth', 1, 'LineStyle','--');
end

xlabel('Time(sec)');
ylabel('Estimated speaker location');

hold off;

if ~exist(strcat(Profile, '/plots/', Date, '/', audio_name), 'dir')
    mkdir(strcat(Profile, '/plots/', Date, '/', audio_name));
end

saveas(gcf, strcat(Profile, '/plots/', Date, '/', audio_name, '/', audio_name, '.png'));

if PlotIndividualChannelEstimates == 1

    % Plot channels for each microphone pair and save as png files.
    figure(9)
    
    for i = 1:1:length(micPairs)
        
        % Clear figure.
        clf
        
        % Save microphone index.
        mic1 = micPairs(i, 1);
        mic2 = micPairs(i, 2);
        
        % DOA estimates.
        ch1 = doa_ch1(mic1, mic2, :);
        ch2 = doa_ch2(mic1, mic2, :);
        
        x_axis = 1:1:length(ch1);
    
        % Plot channel 1.
        subplot(2,1,1);
        plot(x_axis, ch1(:));
        title(strcat('Mic Pair [', num2str(mic1), ',', num2str(mic2), ']', ' ch 1'));
    
        % Plot channel 2.
        subplot(2,1,2);
        plot(x_axis, ch2(:));
        title(strcat('Mic Pair [', num2str(mic1), ',', num2str(mic2), ']', ' ch 2'));
        
        if ~exist(strcat(Profile, '/doa_plots/', Date, '/', audio_name), 'dir')
            mkdir(strcat(Profile, '/doa_plots/', Date, '/', audio_name));
        end
        saveas(gcf, strcat(Profile, '/doa_plots/', Date, '/', audio_name, ...
               '/', num2str(mic1), '_', num2str(mic2), '.png'));
    end
end
