function v = v_NN(t, x, f, h, delta, sub_ori, sub_pos, t0, dt, pos_matrix_NN)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% V_NN 
    % Evaluates the tracer velocity using only the nearest actuator.
    %
    %   INPUT
    %   -t              :Time at which the velocity is computed at
    %   -x              :Position at which the velocity is computed
    %   -f              :Strength of the Stokeslets
    %   -h              :Height of the Stokeslets
    %   -delta          :Radius of the Stokeslets
    %   -sub_ori        :1xNp cell of orientations of the force for each Stokeslet
    %   -sub_pos        :1xNp cell of 2xNt positions for each Stokeslet
    %   -t0             :Initial time
    %   -dt             :Time step
    %   -pos_matrix_NN  :2xNtxNp matrix containing the position of each
    %                    Stokeslet
    %   OUTPUT
    %   -v              :3x1 computed tracer velocity
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    idx_t = round((t - t0) / dt) + 1;
    idx_t = max(1, min(idx_t, size(pos_matrix_NN, 2)));
    dx = pos_matrix_NN(1, idx_t, :) - x(2);
    dy = pos_matrix_NN(2, idx_t, :) - x(3);
    dists_sq = dx.^2 + dy.^2;
    dists_sq = dists_sq(:).';
    [minimum, min_idx] = min(dists_sq);
    
    if  ~isnan(minimum)
        pos_actuator=sub_pos{min_idx};
        pos_actuator_t=pos_actuator(:,idx_t);
        ori_actuator_t=sub_ori{min_idx};
        f=[0;cos(ori_actuator_t(idx_t));sin(ori_actuator_t(idx_t))]*f;
        v= Reg_Stokeslet_velocity_near_wall(x-pos_actuator_t+[h;0;0],f,h,delta);
    else
        v=[0;0;0];
    end
end