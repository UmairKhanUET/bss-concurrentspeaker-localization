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

function [localEstimates] = FindLocalEstimates( globalEstimates, micPairs )
%FindLocalEstimates - Finds local speaker location estimates.
%   Takes global speaker location estimates in globaEstimates array and
%   prints local estimates on a per microphone pair basis.

L = 0.0926;
TotalAudioChannels = 6;

figure(200);

hold on

% Draw the hexagonal boundary of the respeaker.
hexagon = drawHexagon(L/2, 30);

hold off

% Create and initialize the vertices array
vertices = struct('Xdata', zeros(1, TotalAudioChannels), 'Ydata', zeros(1, TotalAudioChannels));

% Save vertices of the hexagon. The index denotes the microphone number
% marked on the ReSpeaker periphery.
for i = 1:1:TotalAudioChannels
    vertices(i).Xdata = hexagon.XData((TotalAudioChannels-i)+2);
    vertices(i).Ydata = hexagon.YData((TotalAudioChannels-i)+2);
end

close(figure(200));

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

            % Calculate angle.
            lineAngles(i, j) = (atan2(y_component, x_component)*180/pi);
        end
    end
end

% Now lets calculate normals to these line angles.
normals = lineAngles + 90;

[numMicPairs, ~] = size(micPairs);

localEstimates = zeros(numMicPairs, length(globalEstimates));

for j = 1:1:length(globalEstimates)

    for i = 1:1:numMicPairs

        % Select the microphone pair.
        mic1 = micPairs(i, 1);
        mic2 = micPairs(i, 2);

        theta = (normals(mic1, mic2) + (globalEstimates(j) + 30));
        
        % Restrict theta to between 0 to 360.
        theta = mod(theta + 360, 360);
        
        if ((theta > 90) && (theta <= 270))
            theta = 180 - theta;
        elseif (theta > 270)
            theta = theta - 360;
        end
        
        localEstimates(i, j) = theta;
    end
end

end

