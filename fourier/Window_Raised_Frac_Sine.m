function windowSine = Window_Raised_Frac_Sine(windowLength,taperLength,fractionalOrder,minmaxValues)

% Window_Raised_Frac_Sine: Computes smooth, short-tapered 1D window based on the power
%
% windowSine = Window_Raised_Frac_Sine(windowLength,taperLength,fractionalOrder,minmaxValues) 
%
% Input:        windowLength, taperLength  (int), fractionalOrder (double),
%               minmaxValues (scalar or 2-length array of double)
% Output:       windowSine (array of double)
% Example:
%               windowSine = Window_Raised_Frac_Sine(128,20,2,[0.1 1]);
% Uses:         
% Used in:      
% Comments:     
% Notes:
% Created:      2004/05/06
% Modified:     2006/06/13     
% Modified:     2007/02/13
%               minmaxValues added to avoid vanishing windows
%
%   Author: Laurent C. Duval, laurent.duval@ifpen.fr
%   Institution: IFP Energies nouvelles, Technology Department
%   (c) All right reserved

if nargin < 1
    % Default length is 128
    windowLength = 128;
end
if nargin < 2
    % Default taper length on both sides is approximately 1/16th of the window length
    taperLength = ceil(windowLength/16);
end
if nargin < 3
    % Default exponent for the 'cosine part' is 1, e.g. 'traditional' window
    fractionalOrder = 1;
end
if nargin < 4
    % Default window extend ranges from 0 to 1
    minmaxValues = [0 1]';
end
if (nargin == 4) & (length(minmaxValues) == 1)
    % If a scalar is given, test value with 0.5 to infer min or max value
    if minmaxValues <= 0.5
        minmaxValues = [minmaxValues 1]';
    else
        minmaxValues = [0 minmaxValues]';
    end
end

taper = (cos((0:taperLength-1)'*pi/(taperLength-1))+1)/2;
windowSine = [flipud(taper);ones(windowLength-2*taperLength,1);taper];
windowSine = windowSine.^fractionalOrder;
windowSine = (windowSine - min(windowSine))/(max(windowSine) - min(windowSine)) * (max(minmaxValues) - min(minmaxValues)) + min(minmaxValues);

if nargout < 1
    plot(windowSine);axis tight
end
