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

function [ alpha, alpha_phantom, alpha_locus, alpha_phantom_locus, pHandle1, ...
           pHandle2] ...
         = drawGeometryDynamic( theta_r, theta, length, x1, y1, ...
                                direction, circleRadius )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
L = length;
x_axis = L/2*cosd(theta_r);
y_axis = L/2*sind(theta_r);
alpha = (theta_r - theta);
alpha_phantom = 180 + theta + theta_r;

% Keep alpha and alpha_phantom between 0 - 360 degrees.
alpha = mod(alpha + 360, 360);
alpha_phantom = mod(alpha_phantom + 360, 360);

x2 = x1+(L*cosd(alpha));
y2 = y1+(L*sind(alpha));

% Find intersection points with the circle of Radius L.
[x2, y2] = findIntersection(x1, y1, x2, y2, circleRadius, alpha);
if (direction == draw.Inward || direction == draw.OutwardInward)
    pHandle1 = plot([x1 x2],[y1 y2], '-red');
end
plot([x1 x1+x_axis],[y1 y1+y_axis], '-.blue');
x2_phantom = x1+(L*cosd(alpha_phantom));
y2_phantom = y1+(L*sind(alpha_phantom));

% Find intersection points with the circle of Radius L.
[x2_phantom, y2_phantom] = findIntersection(x1, y1, x2_phantom, ...
                                            y2_phantom, circleRadius, ...
                                            alpha_phantom);

if (direction == draw.Outward || direction == draw.OutwardInward)
    pHandle2 = plot([x1 x2_phantom],[y1 y2_phantom], '-green');
    
    if (direction == draw.Outward)
        pHandle1 = pHandle2;
    end
else
    pHandle2 = pHandle1;
end
plot([x1-x_axis x1],[y1-y_axis y1], '-.blue');

% Populate values to be returned.
alpha_locus = [x2, y2];
alpha_phantom_locus = [x2_phantom, y2_phantom];

end

function [xi, yi] = findIntersection(x1, y1, x2, y2 ,r, alpha)
    % First make sure alpha isn't 90 or -90 degrees. If it is so, linecirc
    % won't return a valid intersection point and we'll have to go about it
    % in an different way.
    if (abs(alpha) ~= 90)
        slope = (y2 - y1) / (x2 - x1);    % Find slope of the line.
        intercept = (y2) - slope * x2;    % Find y intercept. y = mx + c
        [x_i, y_i] = linecirc(slope, intercept, 0, 0, r);

        % Find which of the retuned intersection points matches the
        % direction of the vector. I calculated angle from starting point
        % till the tip and matched it with the alpha. If it matches, this
        % is the desired point. Otherwise, pick the other one.
        if (mod(abs((atan2(y_i(1), x_i(1)) * 180 / pi) - alpha) + 360, 360)...
            < 0.001)
            xi = x_i(1);
            yi = y_i(1);
        else
            xi = x_i(2);
            yi = y_i(2);
        end
    else
        % alpha is either -90 or +90. To find the intersection point, I am
        % using the equation of circle. Since, x^2 + y^2 = r^2, as the
        % circle is cenetered at origin and we already know x (x1 = x2 =
        % x), y can simply be calculated as y = sqrt(r^2 - x^2). Now, if
        % alpha is negative, return -y. Otherwise, return +y.
        xi = x1;
        yi = sqrt(r^2 - x1^2);
        if (alpha == -90)
            yi = yi * -1;
        end
    end
end
