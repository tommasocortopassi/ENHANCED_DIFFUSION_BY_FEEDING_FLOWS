function [msd_xy, msd_z, correlations_angle_xy,correlations_z] = ...
    segmented_trajectories_msd(displacements_xy, displacements_z, ...
    q_norm, side_norm, max_lag)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% SEGMENTED_TRAJECTORIES_MSD
    %
    % The function takes as input data from numerical simulations of advected
    % tracers, and returns MSD (in xy and z) and angular
    % correlations, computed by considering only displacements whose norm is
    % above or below a given threshold quantile, setting other displacements to zero.
    %
    %
    % INPUT
    %
    % - displacements_xy        :   Cell array containing displacements_xy (divided per
    %                               encounter) of the tracer particle across different
    %                               simulated trajectories.
    %
    % - displacements_z         :   Cell array containing z-displacements (divided per
    %                               encounter) of the tracer particle across different
    %                               simulated trajectories.
    %
    % - q_norm                  :   Number in [0,1]. Quantile threshold to consider,
    %                               relative to the norm of the cleaned displacements_xy.
    %
    % - side_norm               :   String. If "left" we consider displacements_xy
    %                               whose norm is below the q_norm threshold. If
    %                               "right", we consider displacements_xy whose norm is
    %                               above.
    %
    % - max_lag                 :   Maximum lag (in terms of number of relevant
    %                               encounters) considered for the angular correlation.
    %                               (Optional) Default = 10
    %
    % OUTPUT
    %
    % - msd_xy                  :   Value of the mean squared displacements at the
    %                               final time with respect to horizontal xy variables,
    %                               computed by considering only relevant displacements
    %                               in each cell of displacement_xy, i.e. only 
    %                               displacements satisfying the threshold condition.
    %
    % - msd_z                   :   Value of the mean squared displacements at the
    %                               final time with respect to vertical z variable,
    %                               computed by considering only relevant displacements
    %                               in each cell of displacement_z, i.e. only 
    %                               displacements satisfying the threshold condition.
    %
    % - correlations_angle_xy   :   Row vector. Correlation between relevant 
    %                               displacement_xy directions, as a function 
    %                               of lag (in units of relevant-encounter 
    %                               index, not physical time).
    %
    % - correlations_z          :   Row vector. Correlation between relevant 
    %                               displacement_z directions (up/down), as a 
    %                               function of lag (in units of 
    %                               relevant-encounter index, not physical time).
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % ==============================
    %% Check arguments
    % ==============================
    
    if nargin < 5 || isempty(max_lag)
        max_lag = 10;
    end
    

    if ~iscell(displacements_xy) || ~iscell(displacements_z) 
        error('displacements_xy, displacements_z, encounter_times and id_actuators must be cell arrays.');
    end
    
    
    if ~(isscalar(q_norm) && q_norm >= 0 && q_norm <= 1)
        error('q_norm must be a number between 0 and 1.');
    end
    
    if ~ismember(side_norm, {'left','right'})
        error('side_norm must be either ''left'' or ''right''.');
    end
    
    if ~isscalar(max_lag) || max_lag < 1 || floor(max_lag) ~= max_lag
        error('max_lag must be a positive integer.');
    end
    
    
    % =========================================================================
    %% STEP 1: Norm threshold, computed on the cleaned (merged) displacements
    % =========================================================================
    n_cells = numel(displacements_xy);
    all_norms = [];
    D_xy_filtered   = cell(1,n_cells);
    D_z_filtered    = cell(1,n_cells);
    for k = 1:n_cells
        
        if isempty(displacements_xy{k}) || isempty(displacements_z{k}) 
            D_xy = zeros(2,0);
            D_z  = zeros(1,0);
            continue;
        else
            D_xy = displacements_xy{k};
            D_z  = displacements_z{k};
        end

        if size(D_xy,1) ~= 2 && size(D_xy,2) == 2
            D_xy = D_xy.';
        end
        
        if size(D_xy,1) ~= 2
            error('Each cell of displacements_xy must contain either a 2xN or a Nx2 matrix.');
        end
    
        if size(D_z,1) ~= 1 && size(D_z,2) == 1
            D_z = D_z.';
        end
        
        if size(D_z,1) ~= 1
            error('Each cell of displacements_z must contain either a 1xN or a Nx1 vector.');
        end
    
        N = size(D_xy,2);
    
        if numel(D_z) ~= N
            error('In cell %d, displacements_xy, displacements_z must have matching length.', k);
        end
    
        norms = vecnorm([D_xy;D_z], 2, 1);
        valid = isfinite(norms) & ~isnan(norms);
        all_norms = [all_norms, norms(valid)]; 
    end
    
    if isempty(all_norms)
        error('No valid encounters found after merging aligned displacements.');
    end
    
    thr_norm = quantile(all_norms, q_norm);
    
    % =========================================================================
    %% STEP 2: Zero-out non-relevant displacements, and
    %%         build a separate compressed sequence for angular correlation
    % =========================================================================
    
    for k = 1:n_cells
        D_xy  = displacements_xy{k};
        D_z = displacements_z{k};
        norms = vecnorm([D_xy;D_z], 2, 1);
        if strcmp(side_norm, 'right')
            mask_norm = norms >= thr_norm;
        else
            mask_norm = norms <= thr_norm;
        end
        mask_keep = mask_norm & isfinite(norms)&  ~isnan(norms);
        
        D_xy_filtered{k}  =  D_xy(:, mask_keep);
        D_z_filtered{k} = D_z(mask_keep);
    end
    
    xy_cum              = zeros(1,n_cells);
    z_cum               = zeros(1,n_cells);
    sum_angle           = zeros(1,max_lag);
    sum_z_directions    = zeros(1, max_lag);
    cnt_angle           = zeros(1,max_lag);
    for k = 1:n_cells
    
        xy_cum(k) = sum(sum(D_xy_filtered{k},2).^2);
        z_cum(k) = sum(D_z_filtered{k})^2;
    
        Dcc  = D_xy_filtered{k};
        Dzcc = D_z_filtered{k};
        
        Ld = size(Dcc, 2);
        if Ld > 1
            maxLagEff = min(max_lag, Ld - 1);
            for lag = 1:maxLagEff
                D1 = Dcc(:, 1:end-lag);
                D2 = Dcc(:, 1+lag:end);
                
                Dz1=Dzcc(1:end-lag);
                Dz2=Dzcc(1+lag:end);
    
                n1 = vecnorm(D1, 2, 1);
                n2 = vecnorm(D2, 2, 1);
                valid = (n1 > 1e-12) & (n2 > 1e-12);
    
                if any(valid)
                    dot_vals = sum(D1(:,valid) .* D2(:,valid), 1);
                    ang_vals = dot_vals ./ (n1(valid) .* n2(valid));
                    sum_z_directions(lag) = sum_z_directions(lag) + sum(sign(Dz1(valid).*Dz2(valid)));
                    sum_angle(lag) = sum_angle(lag) + sum(ang_vals);
                    cnt_angle(lag) = cnt_angle(lag) + numel(ang_vals);
                end
            end
        end
    end
    
    % ============================================================
    %% Assemble outputs
    % ============================================================
    
    msd_xy = sum(xy_cum)/ n_cells;
    msd_z   = sum(z_cum)  / n_cells;
    
    
    valid_lag = cnt_angle > 0;

    correlations_angle_xy = nan(1, max_lag);
    correlations_angle_xy(valid_lag) = sum_angle(valid_lag) ./ cnt_angle(valid_lag);

    correlations_z = nan(1, max_lag);
    correlations_z(valid_lag) = sum_z_directions(valid_lag) ./ cnt_angle(valid_lag);
end