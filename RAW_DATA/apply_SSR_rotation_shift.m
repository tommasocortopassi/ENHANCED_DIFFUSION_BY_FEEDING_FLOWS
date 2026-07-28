function particle_out = apply_SSR_rotation_shift(particle, fixed_angle, angle_thr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% APPLY_SSR_ROTATION_SHIFT
% Applies a cumulative fixed angular offset to the orientation trace
% every time an SSR rotation event is detected.
%
% Segments are maximal runs of consecutive non-NaN samples in row 3
% (orientation) of each particle field. For every pair of consecutive
% segments, the boundary is
%
%   wrapped_diff = wrapToPi( theta_first(next_segment) - theta_first(prev_segment) )
%
%   Filtering rule applied to wrapped_diff:
%   1) wrapped_diff > 0  (counter clockwise rotation)   ->  Discarded.
%   2) wrapped_diff <= 0                                ->  if >=  ANGLE_THRESH (20 deg by default), 
%                                                           it is considered a valid SSR. 
%                                      
% Every time an event is detected, FIXED_ANGLE is added to the offset, 
% and from that segment onward the (possibly further incremented) offset is 
% added to every remaining sample until the next event
% increases it again. Segment detection itself always uses the
% ORIGINAL (un-shifted) first-sample values, so the event boundaries
% are found on the raw data regardless of the accumulated correction.
%
% INPUT
%   particle     :  Structure produced by Data_collect, fields p1, p2, ...
%                   each a 3 x Nt matrix (row 3 = force orientation, rad).
%   fixed_angle  :  Angle (degree) added to the running offset every time a
%                   rotation event is detected.
%   angle_thr    :  (optional) Minimum |wrapped_diff| (deg) to qualify as
%                   a rotation event. Default = 20.
%
% OUTPUT
%   particle_out :  Same structure as PARTICLE, with row 3 (orientation)
%                   shifted by the running cumulative offset from the
%                   segment where each event occurs to the end of the
%                   trace, re-wrapped to (-pi, pi]. NaN samples are left
%                   as NaN.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 3
    angle_thr = 20;
end

ANGLE_THRESH = angle_thr * pi/180;
fixed_angle=fixed_angle*pi/180;

particle_out = particle;
fields = fieldnames(particle);

for k = 1:numel(fields)
    %% SINCE THE CONVENTION IN THE DATASET REFERS TO CLOCKWISE ANGLES
    %% AS POSITIVE (SEE NOTES.MD) WE CORRECT IT BY TAKING THE OPPOSITE
    %% ANGLES
    theta = -particle.(fields{k})(3, :);
    theta = theta(:).';

    valid = ~isnan(theta);
    if ~any(valid)
        continue
    end

    starts = find(valid & [true, ~valid(1:end-1)]);
    ends   = find(valid & [~valid(2:end), true]);

    first_vals = theta(starts);   % raw, un-shifted reference values

    offset = 0;
    for i = 1:numel(starts)
        if i > 1
            raw_diff     = first_vals(i) - first_vals(i-1);
            wrapped_diff = mod(raw_diff + pi, 2*pi) - pi;   % wrapToPi

            if abs(wrapped_diff) >= ANGLE_THRESH && wrapped_diff<0 
                offset = offset + fixed_angle;
            end
        end

        idx = starts(i):ends(i);
        theta(idx) = mod(theta(idx) + offset + pi, 2*pi) - pi;
    end

    particle_out.(fields{k})(3, :) = theta;
end
end
