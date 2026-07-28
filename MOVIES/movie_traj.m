function v = movie_traj(particles, t, tracer, arrow_length, show_trail,framerate,frames_target)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% MOVIE_TRAJ
    %
    % Creates a 2D (planar) movie of the motion of the tracer under th influence of 
    % the velocity field created by particles. Frame rate is fixed at 24 fps and th total 
    % duration is fixed at 20 s. If the number of time steps exceeds 480 (= 24 fps × 20 s), 
    % trajectories are downscaled to exactly 480 evenly-spaced frames.
    %
    % You can disable particles plotting by passing particles = [].
    %
    % -------------------------------------------------------------------------
    % REQUIRED INPUTS
    %   particles  : Struct of particle data, or [] to disable particle plotting.
    %                Fields: 'p1', 'p2', ..., 'pM'.
    %                Each field is a (3×N) matrix:
    %                  rows 1–2 : x, y positions
    %                  row  3   : orientation of forcing angle theta [rad]   (if 3 rows)
    %
    %   t          : (1×N) or (N×1) time vector.
    %
    % -------------------------------------------------------------------------
    % OPTIONAL INPUTS  (pass [] to use default)
    %   tracer      :   Struct of tracer data. Fields: 't1', 't2', ..., 'tK'.
    %                   Each field is a (3×N) matrix in [z; x; y] row order.
    %                   Pass [] or omit to disable tracers.
    %
    %   arrow_length :  Scalar. Length of the arrow of forcing directions of
    %                   particles. 1 by default.
    %
    %   show_trail   :  Logical (default: false).
    %                   true  → tracers draw a persistent growing path + dot.
    %                   false → tracers show only the current position dot.
    %
    %   framerate    :  Desired framerate. Default = 24 frames/s.
    %
    %   frames_target:  Desired total number of frames. Default = framerate x 20
    %
    %
    % -------------------------------------------------------------------------
    % OUTPUT
    %   v  - VideoWriter object (closed). File written as 'particelle2D_N.avi'
    %        where N is the first integer that does not collide with an existing
    %        file in the current directory.
    %
    % -------------------------------------------------------------------------
    % NOTES
    %   - FrameRate = 24 fps. Duration ≈ 20 s. Frames = min(N, 480).
    %   - Particle markers: blue. Orientation arrows: black. Tracers: red.
    %   - If particles = [], only tracer trajectories are shown.
    %   - Tracer row order [z;x;y] is permuted to [x;y;z] internally.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % ---- Defaults ------------------------------------------------------------
    if nargin < 1, particles = []; end
    if nargin < 3 || isempty(tracer),     tracer     = struct(); end
    if nargin < 4 || isempty(arrow_length),      arrow_length      = 1;        end
    if nargin < 5 || isempty(show_trail), show_trail = false;    end
    if nargin < 6 || isempty(framerate), framerate = 24; end
    
    % =========================================================================
    % 0.  INPUT FLAGS
    % =========================================================================
    has_particles = ~isempty(particles);
    has_tracers   = ~isempty(tracer) && ~isempty(fieldnames(tracer));
    
    if ~has_particles && ~has_tracers
        error('At least one of particles or tracer must contain data.');
    end
    
    
    
    % =========================================================================
    % 1.  PREPROCESSING
    % =========================================================================
    
    % -- 1a. Rescale particles particles (x,y only) -----------------------------------
    if has_particles 
        pnames_r = fieldnames(particles);
        for m = 1:numel(pnames_r)
            nm = pnames_r{m};
            particles.(nm)(1:2,:)=particles.(nm)(1:2,:)/200;
        end
    end
    
    % -- 1b. Permute tracer rows from [z;x;y] → [x;y;z] -----------------------
    if has_tracers
        tnames_r = fieldnames(tracer);
        for m = 1:numel(tnames_r)
            nm = tnames_r{m};
            tracer.(nm) = tracer.(nm)([2, 3, 1], :);
        end
    end
    
    % ---- Dimensions ----------------------------------------------------------
    N = numel(t);
    
    if has_particles
        pnames = fieldnames(particles);
        M      = numel(pnames);
    
        sample_particle = particles.(pnames{1});
        if size(sample_particle,1) >= 4
            ori_row = 4;
            has_ori = true;
        elseif size(sample_particle,1) >= 3
            ori_row = 3;
            has_ori = true;
        else
            ori_row = [];
            has_ori = false;
        end
    else
        pnames   = {};
        M        = 0;
        ori_row  = [];
        has_ori  = false;
    end
    
    if has_tracers
        tnames    = fieldnames(tracer);
        N_tracer  = numel(tnames);
    else
        tnames    = {};
        N_tracer  = 0;
    end
    
    % =========================================================================
    % 2.  DOWNSAMPLE TO 480 FRAMES MAX  (24 fps × 20 s)
    % =========================================================================
    FRAMERATE_TARGET = framerate;
    
    if nargin < 7 || isempty(frames_target)
        FRAMES_TARGET = framerate * 20;
    else
        FRAMES_TARGET    = frames_target;
    end 
    
    if N > FRAMES_TARGET
        idx_ds = round(linspace(1, N, FRAMES_TARGET));
        t      = t(idx_ds);
    
        for m = 1:M
            nm     = pnames{m};
            particles.(nm) = particles.(nm)(:, idx_ds);
        end
    
        for m = 1:N_tracer
            nm          = tnames{m};
            tracer.(nm) = tracer.(nm)(:, idx_ds);
        end
    
        N = FRAMES_TARGET;
    end
    
    FrameRate = FRAMERATE_TARGET;
    fprintf('Frames: %d  |  FrameRate: %d fps  |  Duration: %.1f s\n', ...
        N, FrameRate, N / FrameRate);
    
    % =========================================================================
    % 3.  COMPUTE AXIS LIMITS
    % =========================================================================
    n_cells = M + N_tracer;
    col_x   = cell(n_cells, 1);
    col_y   = cell(n_cells, 1);
    
    idx = 0;
    
    for m = 1:M
        idx        = idx + 1;
        d          = particles.(pnames{m});
        col_x{idx} = d(1, :);
        col_y{idx} = d(2, :);
    end
    
    for m = 1:N_tracer
        idx        = idx + 1;
        d          = tracer.(tnames{m});
        col_x{idx} = d(1, :);
        col_y{idx} = d(2, :);
    end
    
    all_x = [col_x{1:idx}];
    all_y = [col_y{1:idx}];
    
    all_x = all_x(isfinite(all_x));
    all_y = all_y(isfinite(all_y));
    
    if isempty(all_x) || isempty(all_y)
        error('No finite coordinates available to plot.');
    end
    
    margin = 0.05;
    xmin = min(all_x); xmax = max(all_x); dx = xmax - xmin;
    ymin = min(all_y); ymax = max(all_y); dy = ymax - ymin;
    
    if dx == 0, dx = 1; end
    if dy == 0, dy = 1; end
    
    xl = [xmin - margin*dx, xmax + margin*dx];
    yl = [ymin - margin*dy, ymax + margin*dy];
    
    fprintf('Axis limits — particles: [%.3f, %.3f]  Y: [%.3f, %.3f]\n', ...
        xl(1), xl(2), yl(1), yl(2));
    
    % =========================================================================
    % 4.  OPEN VIDEO WRITER
    % =========================================================================
    filename = get_new_filename('particelle2D', 'mp4');
    fprintf('Saving: %s\n', filename);

    v = VideoWriter(filename, 'MPEG-4');
    v.FrameRate = FrameRate;
    v.Quality = 100;
    open(v);
    
    % =========================================================================
    % 5.  FIGURE & AXES SETUP
    % =========================================================================
    fig = figure('Units', 'normalized', 'Position', [0 0 1 1]);
    hold on;
    grid on;
    axis equal;
    xlim(xl);
    ylim(yl);
    
    set(gca, 'YDir', 'reverse');
    
    xlabel('X');
    ylabel('Y');
    
    % =========================================================================
    % 6.  PRE-ALLOCATE GRAPHICS HANDLES
    % =========================================================================
    c_particle = [0, 0, 1];
    c_tracer   = [0.85, 0, 0];
    c_arrow    = [0, 0, 0];
    
    h_pts = gobjects(M, 1);
    h_arr = gobjects(M, 1);
    
    for m = 1:M
        h_pts(m) = plot(NaN, NaN, 'o', ...
            'MarkerSize', 6, ...
            'MarkerFaceColor', c_particle, ...
            'MarkerEdgeColor', 'k');
    
        if has_ori
            h_arr(m) = quiver(NaN, NaN, NaN, NaN, ...
                'AutoScale', 'off', ...
                'Color', c_arrow, ...
                'LineWidth', 1, ...
                'MaxHeadSize', 2);
        end
    end
    
    h_tr     = gobjects(N_tracer, 1);
    h_tr_dot = gobjects(N_tracer, 1);
    
    for m = 1:N_tracer
        if show_trail
            h_tr(m) = plot(NaN, NaN, '-', ...
                'Color', c_tracer, ...
                'LineWidth', 1.5);
        end
    
        h_tr_dot(m) = plot(NaN, NaN, 'o', ...
            'MarkerSize', 5, ...
            'MarkerFaceColor', c_tracer, ...
            'MarkerEdgeColor', c_tracer);
    end
    
    ht = title('');
    
    % Arrow length
    arrow_len = 5 * arrow_length;
    
    % =========================================================================
    % 7.  RENDER LOOP
    % =========================================================================
    for k = 1:N
    
        % -- Particles ---------------------------------------------------------
        for m = 1:M
            pos = particles.(pnames{m})(:, k);
            set(h_pts(m), 'XData', pos(1), 'YData', pos(2));
    
            if has_ori && all(isfinite([pos(1), pos(2), pos(ori_row)]))
                theta = pos(ori_row);
                set(h_arr(m), ...
                    'XData', pos(1), ...
                    'YData', pos(2), ...
                    'UData', arrow_len*cos(theta), ...
                    'VData', arrow_len*sin(theta));
            else
                if has_ori
                    set(h_arr(m), 'XData', NaN, 'YData', NaN, 'UData', NaN, 'VData', NaN);
                end
            end
        end
    
        % -- Tracers -----------------------------------------------------------
        for m = 1:N_tracer
            d = tracer.(tnames{m});
    
            if show_trail
                set(h_tr(m), 'XData', d(1, 1:k), 'YData', d(2, 1:k));
            end
    
            set(h_tr_dot(m), 'XData', d(1, k), 'YData', d(2, k));
        end
    
        set(ht, 'String', sprintf('t = %.3f', t(k)));
        drawnow limitrate;
        writeVideo(v, getframe(fig));
    end
    
    close(v);
    close(fig);
    
    end
    
    
    function filename = get_new_filename(base, ext)
    i = 1;
    filename = sprintf('%s_%d.%s', base, i, ext);
    while isfile(filename)
        i = i + 1;
        filename = sprintf('%s_%d.%s', base, i, ext);
    end
end