function [times,tracers, indices_NN, min_dis, displacements_xy, displacements_z,encounter_times, MSDxy, MSDz] = ...
    test_displacement_synth_pop(opt)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% TEST_DISPLACEMENT_SYNTH_POP
    %
    % Simulates a tracer particle advected by Np moving actuators confined to
    % a disk of radius R at fixed height h. Returns MSD_xy, MSD_z,
    % nearest-neighbor encounter data, and height distribution PDF. Notice
    % that we use the convention [z;x;y] for 3D Cartesian coordinates.
    %
    % FEATURES
    % ---------------------------------------------------------------------------
    %   • Supports T > dataset duration (concatenates random snippets)
    %   • Two interaction models: nearest-actuator only (only_NN=1) or all
    %     actuators (only_NN=0)
    %   • Reflection boundary: actuators exiting the disk are reflected back at
    %     the diametrically opposed point
    %   • Parallel execution via parfor (Ntrials trials)     
    %   • A maximum of 1000 full tracers trajectories (subsampled according to 
    %     opt.skip value) are saved, to avoid excessive memory consumption
    % INPUT: opt (structure)
    % ---------------------------------------------------------------------------
    % REQUIRED FIELDS:
    %   opt.noise_threshold :   Noise threshold for Data_collect
    %   opt.h               :   Fixed actuator height (z-coordinate)
    %   opt.delta           :   Interaction parameter
    %   opt.f               :   Actuator strength/force
    %   opt.z0              :   Initial tracer z-coordinate
    %   opt.Ntrials         :   Number of independent trials
    %   opt.R               :   Disk radius (actuators confined here)
    %   opt.Np              :   Number of active actuator slots (or NaN for all)
    %
    % OPTIONAL FIELDS (with defaults):
    %   opt.T               :   Simulation time (default: dataset duration)
    %   opt.only_NN         :   Interaction model: 1=nearest only, 0=all (default: 1)
    %   opt.skip            :   Subsampling interval for tracers trajectories (default: 20)
    %   opt.randomise_angle :   Randomise angles flag (default: 0)
    %
    % OUTPUT
    % ---------------------------------------------------------------------------
    %   times               :   Subsampled time vector (Nt_sub x 1)
    %   tracers             :   Structure with tracer.t1...tK (K=min(Ntrials,1000)),
    %                           each 3 x Nt_sub matrix [z; x; y]
    %   indices_NN          :   Cell array, nearest actuator slot indices per encounter
    %   min_dis             :   Cell array, minimum distances per encounter
    %   displacements_xy    :   Cell array, 2 x M [dx; dy] per encounter of
    %                           of horizontal displacements
    %   displacements_z     :   Cell array, 1 x M per encounter
    %                           of vertical displacements
    %   encounter_times     :   Cell array, durations between encounters (seconds)
    %   MSDxy               :   MSD in x-y plane (Nt_sub x 1)
    %   MSDz                :   MSD in z direction (Nt_sub x 1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Check of opt fieldnames and extraction of relevant quantities
    if ~isstruct(opt)
        error('opt must be a structure.');
    end
    
    required_fields = {'noise_threshold','h','delta','f','z0', ...
        'Ntrials','R','Np'};
    
    for k = 1:numel(required_fields)
        if ~isfield(opt, required_fields{k})
            error('Missing required field opt.%s.', required_fields{k});
        end
    end
    
    options         = odeset('AbsTol', 1e-6, 'RelTol', 1e-4);
    noise_threshold = opt.noise_threshold;
    h               = opt.h;
    delta           = opt.delta;
    f               = opt.f;
    z0              = opt.z0;
    Ntrials         = opt.Ntrials;
    R               = opt.R;
    Np              = opt.Np;

    skip=20;
    if isfield(opt, 'skip') && ~isempty(opt.skip)
        skip = opt.skip;
    end
    
    only_NN=0;
    if isfield(opt, 'only_NN')  && ~isempty(opt.only_NN)
        only_NN = opt.only_NN;
    end
    
    randomise_angle = 0;
    if isfield(opt, 'randomise_angle') && ~isempty(opt.randomise_angle)
        randomise_angle = opt.randomise_angle;
    end
    
    if only_NN ~= 0 && only_NN ~= 1
        error('only_NN must be 0 or 1.');
    end
    
    % --- Load data -----------------------------------------------------------
    [times_data, particle] = Data_collect("Dataset1", noise_threshold);
    particle               = select_synth_population(particle, randomise_angle);
    
    if numel(times_data) < 2
        error('times_data must contain at least two samples.');
    end
    
    Nt_data = numel(times_data);
    dt      = times_data(2) - times_data(1);
    t0      = times_data(1);
    
    if ~(isfinite(dt) && dt > 0)
        error('Invalid time step dt.');
    end
    
    if isfield(opt, 'T') && ~isempty(opt.T)
        T = opt.T;
    else
        T = times_data(end) - times_data(1);
    end
    
    N_req         = round(T / dt) + 1;
    times_full    = t0 + (0:N_req-1) * dt;
  
    subsample_idx = 1:skip:N_req;
    times         = times_full(subsample_idx);
    Nt_sub        = numel(times);
    
    % --- Build source library ------------------------------------------------
    all_names = fieldnames(particle);
    Np_names  = all_names(startsWith(string(all_names), "p"));
    Np_total  = numel(Np_names);
    
    if Np_total == 0
        error('No actuator fields found in particle structure.');
    end
    
    pos_lib = cell(1, Np_total);
    ori_lib = cell(1, Np_total);
    
    for i = 1:Np_total
        A          = particle.(Np_names{i});
        %% We use Lc= 200 micrometers as unit of measurement for length. Hence
        %% we divide the coordinates by 200
        pos_lib{i} = [h * ones(1, Nt_data); A(1:2, :) / 200];
        ori_lib{i} = A(3, :);
    end
    
    n_sources = numel(pos_lib);
    
    if isempty(Np) || (isscalar(Np) && isnan(Np))
        Np = n_sources;
    end
    
    % --- Output arrays -------------------------------------------------------
    
    %To avoid excessive use of memory, I save at most 1000 full tracers
    %trajectories
    
    min_dis             =   cell(1, Ntrials);
    indices_NN          =   cell(1, Ntrials);
    displacements_xy    =   cell(1, Ntrials);
    displacements_z     =   cell(1,Ntrials);
    encounter_times     =   cell(1, Ntrials);
    tracers             =   cell(1,min(Ntrials,1000));
    
    MSDxy_sum           = zeros(Nt_sub, 1);
    MSDz_sum            = zeros(Nt_sub, 1);
    
    
    % =========================================================================
    %% MAIN LOOP OVER TRIALS
    % =========================================================================
    parfor trial_idx = 1:Ntrials
        rng(trial_idx, 'twister')
    
        % -----------------------------------------------------------------
        %% STEP 1: build full actuator trajectories over N_req frames
        % -----------------------------------------------------------------
        actuator_pos = zeros(3, N_req, Np);
        actuator_ori = NaN(N_req, Np);
    
        for j = 1:Np
    
            % Initial random placement inside the disk
            r_rand      = R * sqrt(rand);
            theta_rand  = 2 * pi * rand;
            current_x   = r_rand * cos(theta_rand);
            current_y   = r_rand * sin(theta_rand);
    
            % Pick a valid starting snippet: pos_lib{src}(:,frame) and ori_lib{src}(:,frame) have
            % no NaN values 
            src   = randi(n_sources);
            frame = find(~isnan(ori_lib{src}) & all(~isnan(pos_lib{src}(2:3,:)), 1), 1, "first");
    
            t_write = 1;   % next frame to write in the output array
    
            while t_write <= N_req
                frames_left      = N_req - t_write + 1;
                frames_available = Nt_data - frame + 1;
                frames_to_take   = min(frames_left, frames_available);
                frame_end        = frame + frames_to_take - 1;
    
                raw_pos = pos_lib{src}(:, frame:frame_end);   % 3 x frames_to_take
                raw_ori = ori_lib{src}(frame:frame_end);      % 1 x frames_to_take
    
                % Skip snippet if the first point is not finite
                if ~(all(isfinite(raw_pos(2:3, 1))) && isfinite(raw_ori(1)))
                    src   = randi(n_sources);
                    frame = find( ~isnan(sum(ori_lib{src},1)) & ~isnan(ori_lib{src}),1,"first");
                    continue;
                end
    
                % Rigid transform: rotate then translate
                current_ang= 2*pi*rand;
                Rmat = [ cos(current_ang),  sin(current_ang);
                    -sin(current_ang),  cos(current_ang)];
    
                tformed_xy = Rmat * (raw_pos(2:3, :) - raw_pos(2:3, 1)) + [current_x; current_y];
                tformed_ori = mod(raw_ori - current_ang, 2*pi);
    
                write_range = t_write:(t_write + frames_to_take - 1);
                actuator_pos(1, write_range, j)   = h;
                actuator_pos(2:3, write_range, j) = tformed_xy;
                actuator_ori(write_range, j)      = tformed_ori;
    
                % Update to the last finite point written
                finite_cols = all(isfinite(tformed_xy), 1);
                if any(finite_cols)
                    last_col   = find(finite_cols, 1, 'last');
                    current_x  = tformed_xy(1, last_col);
                    current_y  = tformed_xy(2, last_col);
                end
    
                t_write = t_write + frames_to_take;
                frame   = frame_end + 1;
    
                % Snippet exhausted — pick a new one for the next iteration
                if t_write <= N_req
                    src   = randi(n_sources);
                    frame = find( ~isnan(sum(ori_lib{src},1)) & ~isnan(ori_lib{src}),1,"first");
                end
            end
    
            % Compute distance from origin for all time steps
            dist_from_origin = vecnorm(actuator_pos(2:3, :, j), 2, 1);  % 1 x N_req
            
      
            % Find the first time step where actuator exits the circle
            exit_idx = find(dist_from_origin > R, 1);
    
            while ~isempty(exit_idx)
                % Get the position at exit
                xy_exit = actuator_pos(2:3, exit_idx(1), j);     
    
                % Normal vector at the boundary point
                n = xy_exit / norm(xy_exit);
    
                actuator_pos(2:3, exit_idx(1):end, j) = ...
                    actuator_pos(2:3, exit_idx(1):end, j) - 2*R*n;
    
                % Update distances (optional, for verification)
                dist_from_origin(exit_idx(1):end) = ...
                    vecnorm(actuator_pos(2:3, exit_idx(1):end, j), 2, 1);
                
                %Update exit_idx
                exit_idx = find(dist_from_origin > R, 1);
            end
        end
    
    
        % -----------------------------------------------------------------
        %% STEP 3: integrate tracer over full time
        % -----------------------------------------------------------------
        x0 = [z0; 0; 0];
    
        sub_pos_cell = cell(1, Np);
        sub_ori_cell = cell(1, Np);
        for j = 1:Np
            sub_pos_cell{j} = actuator_pos(:, :, j);
            sub_ori_cell{j} = actuator_ori(:, j)';
        end
    
        if only_NN == 0
            [~, X_ode]      =   ode45(@(t, x) total_v(t, x, f, h, delta, sub_ori_cell,...
                                sub_pos_cell, times_full(1), dt),times_full, x0, options);
        else
            pos_matrix_NN   =   actuator_pos(2:3, :, :);
            [~, X_ode]      =   ode45(@(t, x) v_NN(t, x, f, h, delta, sub_ori_cell,...
                                sub_pos_cell, times_full(1), dt, pos_matrix_NN),...
                                        times_full, x0, options);
        end
    
        tracer              =   X_ode.';
        tracer_sub          =   tracer(:, subsample_idx);
        if trial_idx<=1000
            tracers{trial_idx}  =   tracer(:,subsample_idx);
        end
        % -----------------------------------------------------------------
        %% STEP 4: compute nearest-actuator index and distance
        % -----------------------------------------------------------------
        tracer_xy   = tracer(2:3, :);
        tracer_z    = tracer(1,:);
        dist_full   = inf(1, N_req);
        nn_id_full  = zeros(1, N_req);
    
        for t = 1:N_req
            p = tracer_xy(:, t);
            for j = 1:Np
                q = actuator_pos(2:3, t, j);
                if all(isfinite(q))
                    d = norm(q - p);
                    if d < dist_full(t)
                        dist_full(t)  = d;
                        nn_id_full(t) = j;
                    end
                end
            end
        end
    
        % -----------------------------------------------------------------
        %% STEP 5: extract encounter events
        % -----------------------------------------------------------------
        change_flags = [true, diff(nn_id_full) ~= 0];
        event_pos    = find(change_flags);
    
        if event_pos(end) == N_req
            event_pos(end) = [];
        end
    
        if isempty(event_pos)
            indices_NN{trial_idx}       = zeros(1, 0);
            min_dis{trial_idx}          = zeros(1, 0);
            encounter_times{trial_idx}  = zeros(1, 0);
            displacements_xy{trial_idx} = zeros(2, 0);
            displacements_z{trial_idx}  = zeros(1, 0);
    
        else

        event_pos_ext = [event_pos, N_req];
        enc_raw       = diff(times_full(event_pos_ext));
        disp_raw_xy      = diff(tracer_xy(:, event_pos_ext), 1, 2);
        disp_raw_z= diff(tracer_z(event_pos_ext));
        id_raw        = nn_id_full(event_pos);
        mindist_raw   = dist_full(event_pos);

        valid_evt     = isfinite(mindist_raw) & isfinite(enc_raw) & (enc_raw > 0);
        
        indices_NN{trial_idx}      = id_raw(valid_evt);
        min_dis{trial_idx}         = mindist_raw(valid_evt);
        encounter_times{trial_idx} = enc_raw(valid_evt);
        displacements_xy{trial_idx}   = disp_raw_xy(:, valid_evt);
        displacements_z{trial_idx}   =  disp_raw_z(valid_evt); 
            
        end
    
        % -----------------------------------------------------------------
        %% STEP 6: accumulate statistics (MSDxy, MSDz)
        % -----------------------------------------------------------------
        dx = tracer_sub(2, :) - tracer_sub(2, 1);
        dy = tracer_sub(3, :) - tracer_sub(3, 1);
        dz = tracer_sub(1, :) - tracer_sub(1, 1);
    
        msdxy_i = (dx.^2 + dy.^2).';
        msdz_i  = (dz.^2).';
        
    
        MSDxy_sum  = MSDxy_sum + msdxy_i;
        MSDz_sum   = MSDz_sum + msdz_i;
    
    end
    
    % =========================================================================
    %% ASSEMBLE OUTPUTS
    % =========================================================================
    
    MSDxy = MSDxy_sum / Ntrials;
    MSDz  = MSDz_sum  / Ntrials;


end

