function rotation_angles = extract_SSR_angles(particle, nbins, angle_thr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EXTRACT_SSR_ANGLES
% Polar PDF histogram of rotation angles between the FIRST orientation
% values of consecutive valid (non-NaN) segments, pooled over every
% microorganism stored in the "particle" structure (fields p1, p2, ...).
%
% A segment is any maximal run of consecutive non-NaN samples; every
% such segment is used, regardless of its length (no minimum-length
% filtering is applied).
%
% For every pair of consecutive valid segments, the boundary rotation is
%
%   wrapped_diff = wrapToPi( theta_first(next_segment) - theta_first(prev_segment) )
%
% Filtering rule applied to wrapped_diff:
%   1) wrapped_diff > 0  (counter clockwise rotation)   ->  Discarded.
%   2) wrapped_diff <= 0                                ->  if >=  ANGLE_THRESH (20 deg by default), 
%                                                           it is recorded as one rotation angle. 
%                                                           Otherwise it is discarded.
%
% INPUT
%   particle        :   Structure produced by Data_collect, fields p1, p2, ...
%                       each a 3 x Nt matrix (row 3 = force orientation, rad).
%   nbins           :   (optional) Number of angular bins for the polar histogram.
%                       Default = 40.
%
%   angle_thr       :   (optional) If the difference between the angles is,
%                       in magnitude, less than angle_thr (in degrees), we 
%                       discard it. Default = 20.
%
% OUTPUT
%   rotation_angles :   Vector of rotation angles (rad). Every
%                       entry is <= 0 (clockwise convention) 
%                       and >= angle_thr deg (15 degrees by default).
%
% The figure produced is a polarhistogram with Normalization = 'pdf',
% i.e. the area under the circular PDF sums to 1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    nbins = 40;
end

if nargin<3
    angle_thr=20;
end

ANGLE_THRESH = angle_thr * pi/180;      % 15 deg, minimum |rotation| to record standalone

fields = fieldnames(particle);
rotation_angles = [];




for k = 1:numel(fields)
    %% SINCE THE CONVENTION IN THE DATASET REFERS TO CLOCKWISE ANGLES
    %% AS POSITIVE (SEE NOTES.MD) WE CORRECT IT BY TAKING THE OPPOSITE
    %% ANGLES
    theta = -particle.(fields{k})(3, :);  % third row = orientation
    theta = theta(:).';                   % force row vector

    % --- locate contiguous non-NaN segments ---
    valid = ~isnan(theta);
    if ~any(valid)
        continue
    end

    starts = find(valid & [true, ~valid(1:end-1)]);
   

    if numel(starts) < 2
        continue   % need at least two valid segments to form a boundary
    end

    first_vals = theta(starts);   % first sample of each valid segment


    for i = 2:numel(starts)
        raw_diff     = first_vals(i) - first_vals(i-1);
        wrapped_diff = mod(raw_diff + pi, 2*pi) - pi;   % wrapToPi
        if wrapped_diff > 0
            continue   
        end
        
        if abs(wrapped_diff) >= ANGLE_THRESH
            rotation_angles(end+1) = wrapped_diff;
        end
    end
end

if isempty(rotation_angles)
    error('No valid rotation angles could be computed from the input data.');
end

% --- polar PDF histogram ---
figure;
polarhistogram(rotation_angles, nbins, 'Normalization', 'pdf', ...
    'FaceColor', [0.2 0.4 0.8], 'FaceAlpha', 0.7);
h = gca;
h.RTickLabel = [];
end
