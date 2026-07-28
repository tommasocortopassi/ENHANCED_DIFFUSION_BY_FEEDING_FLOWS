function [v_z,D_r,D_z] = Derivation_AD_Coefficients(opt)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% DERIVATION_AD_COEFFICIENTS
    %
    %  Computes the expected vertical net velocity (v_z), the radial diffusion 
    %  coefficient (D_r) and the vertical one (D_z). Convention for coordinates
    %  is [z;x;y].
    %
    %
    %  PHYSICAL MODEL
    %  ---------------------------------------------------------------------------
    %   • Actuator modeled as a regularized Stokeslet at height h with 
    %     regularization radius delta.
    %   • Tracer starts at (z0, r0, 0) and is advected by the actuator's 
    %     velocity field for duration t with probability mass function (PMF)
    %     pmf_t_active.
    %   • Pumping direction θ is uniformly distributed over N_theta discrete 
    %     angles in [0, 2π).
    %   • Radial distance r follows the distribution pmf_r ().
    %
    %  INPUT
    %  ---------------------------------------------------------------------------
    %  opt : structure with fields
    %       - f             : scalar, force / dynamic viscosity
    %       - h             : scalar, Stokeslet height (z-coordinate)
    %       - delta         : scalar, regularization radius
    %       - n             : scalar, number density of actuators
    %       - t             : 1×N_t vector of encounter durations (strictly increasing, positive)
    %       - pmf_t_active  : 1×N_t PMF for active encounter durations
    %       - pmf_t_inactive: 1×N_t PMF for inactive encounter durations
    %       - z             : 1×N_z vector of tracer initial heights
    %       - N_theta       : scalar, number of discrete pumping directions in [0, 2π)
    %       - N_r           : scalar, number of discrete points of distances
    %                         from NN
    %       - R             : Maximal distance from NN (optional, computed 
    %                         autmatically if omitted)
    %
    %  OUTPUT
    %  ---------------------------------------------------------------------------
    %   v_z  : 1×N_z vector. Average vertical speed.
    %   D_r  : 1×N_z vector. Radial diffusion coefficient.
    %   D_z  : 1×N_z vector. Vertical diffusion coefficient.
    %  NOTES
    %  ---------------------------------------------------------------------------
    %   • Average encounter time = dot(t, pmf_t_active + pmf_t_inactive). It is
    %     crucial that pmf_t_active and pmf_t_inactive refer to the SAME vector
    %     t
    %   • Uses parfor parallel loop over N_z heights (requires Parallel Computing Toolbox)
    %   • Velocity field computed via Reg_Stokeslet_velocity_near_wall 
    %     (regularized Stokeslet near wall)
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
    
    % Extract relevant data
    f = opt.f;
    h = opt.h;
    delta = opt.delta;
    N_r = opt.N_r;
    z = opt.z(:).';
    N_theta = opt.N_theta;
    n = opt.n;
    
    % Definition of r 
    if isfield(opt,"R")
        R=opt.R;
    else
        R= sqrt(abs(log(0.001)/(n*pi))); %Only 0.1% of NN aree farther than R
    end
    
    r= linspace(0,R,N_r+1);
    r= r(1:end-1) + diff(r)/2;
    
    %Definition of pmf_r, normalised
    pmf_r=[1-exp(-n*pi*r(1)^2), exp(-n*pi*r(1:end-1).^2)-exp(-n*pi*r(2:end).^2)];
    pmf_r=pmf_r/sum(pmf_r);
    
    t = opt.t(:).';
    
    %Normalise pmf_t_active and pmf_t_inactive for extra safety
    pmf_t_active = opt.pmf_t_active(:).';
    pmf_t_active=pmf_t_active/sum(pmf_t_active);
    
    pmf_t_inactive=opt.pmf_t_inactive(:).';
    pmf_t_inactive=pmf_t_inactive/sum(pmf_t_inactive);
    
    avg_encounter_time=dot(t,pmf_t_active +pmf_t_inactive);
    
    N_z = numel(z);
    N_t = numel(t);
    
    
    if numel(pmf_t_active) ~= N_t
        error('pmf_t_active and t must have the same length.');
    end
    
    if numel(pmf_r) ~= N_r
        error('pmf_r and r must have the same length.');
    end
    
    if any(t <= 0) || any(diff(t) <= 0)
        error('t must be a strictly increasing vector of positive values.');
    end
    
    
    if abs(sum(pmf_t_active) - 1) > 1e-10
        warning('pmf_t_active does not sum to 1.');
    end
    
    if abs(sum(pmf_r) - 1) > 1e-10
        warning('pmf_r does not sum to 1.');
    end
    
    % Definition of theta and pmf_theta
    theta_edges = linspace(0, 2*pi, N_theta+1);
    theta = theta_edges(1:end-1) + diff(theta_edges)/2;
    pmf_theta = 1 / N_theta;
    
    
    v_z=zeros(1,N_z);
    D_z=zeros(1,N_z);
    D_r=zeros(1,N_z);
    
    if isempty(gcp('nocreate'))
        parpool;
    end
    
    % Parllel computation of coefficints at different heights
    parfor i = 1:N_z
        local_MSD_z = 0;
        local_MSD_xy = 0;
        local_drift_z = 0;
        z0 = z(i);
    
        for j = 1:N_r
            r0 = r(j);
    
            for k = 1:N_theta
                x0 = [z0; r0; 0];
                f_dir = [0; cos(theta(k)); sin(theta(k))] * f;
                
                % Compute trajectory of the tracer 
                [~, trajectory] = ode45( ...
                    @(tt,x) Reg_Stokeslet_velocity_near_wall(x, f_dir, h, delta), ...
                    [0, t], x0, options);
                
                % Discard first point --> time t=0
                trajectory = trajectory(2:end, :);
                for l = 1:N_t
                   
                    % Probability of having impact parameters 
                    % (r(j), z(i) , theta(k), t(l))
                    prob_weight = pmf_t_active(l) * pmf_r(j) * pmf_theta;
                    
                    
                    displ_xy_sq = sum((trajectory(l, 2:3) - [r0, 0]).^2);
                    displ_z = trajectory(l,1) - z0;
                    
                    % Update the values of the mean squared displacemnts and
                    % vertical drift
                    local_MSD_xy = local_MSD_xy + displ_xy_sq * prob_weight;
                    local_MSD_z  = local_MSD_z  + displ_z^2  * prob_weight;
                    local_drift_z = local_drift_z + displ_z* prob_weight;
                end
            end
        end
     
        % Definition of the advection-diffusion coefficients at height z(i)
        v_z(i)=local_drift_z/avg_encounter_time;
        D_z(i) = local_MSD_z/(2*avg_encounter_time);
        D_r(i)= local_MSD_xy/(4*avg_encounter_time);
    end
end


