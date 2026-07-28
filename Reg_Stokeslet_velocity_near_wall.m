function v = Reg_Stokeslet_velocity_near_wall(x, f, h, delta)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% REG_STOKESLET_VELOCITY_NEAR_WALL
    %  The function computes the velocity field produced by a regularised
    %  Stokeslet of radius delta, placed at [h,0,0] above a no-slip floor. We
    %  remark that we use the convention [z,x,y] for coordinates, following 
    %  Journal of Computational Physics 227 (2008) 4600–4616, where the
    %  formula for the velocity field can be found at equation (21).
    %  Following the convention of the article, we consider the 
    %  strength f of a Stokeselet as F / mu, where F is the force imposed 
    %  on the fluid by the Stokeslet and mu is the dynamic viscosity of the
    %  fluid.
    %
    %  INPUT
    %
    %   - x       :Column vector. Position of the tracer with respect to the Stokeslet.
    %   - f       :Column vector. Strength of the Stokeslet. By strength, 
    %              we consider the force of the Stokeslet divided by the 
    %              dynamic viscosity of the fluid, considered with its direction
    %   - h       :Height of the regularised Stokeslet from the no-slip floor
    %   - delta   :Radius of the regularised Stokeslet
    % 
    %  OUTPUT
    % 
    %   -v         :Column vector. The resulting velocity field.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   

    %  Constants
    inv8pi = 1 / (8 * pi);
    inv4pi = 1 / (4 * pi);
    d2 = delta^2;
    h2 = h^2;

    % Distance vectors (x_k_star and x_k)
    xk_star = x - [h; 0; 0];
    xk = x + [h; 0; 0]; % x - [-h;0;0] = x + [h;0;0]

    % Common terms for xk_star
    r_star_sq = xk_star(1)^2 + xk_star(2)^2 + xk_star(3)^2;
    denom_star_sq = r_star_sq + d2;
    denom_star = sqrt(denom_star_sq);
    denom_star_3 = denom_star_sq * denom_star;

    % H1(||xk_star||) and H2(||xk_star||) 
    H1_star = inv8pi * (1/denom_star + d2/denom_star_3);
    H2_star = inv8pi / denom_star_3;

    % Common terms for xk
    r_xk_sq = xk(1)^2 + xk(2)^2 + xk(3)^2;
    denom_xk_sq = r_xk_sq + d2;
    denom_xk = sqrt(denom_xk_sq);
    denom_xk_3 = denom_xk_sq * denom_xk;
    denom_xk_5 = denom_xk_sq * denom_xk_3;

    % H1(||xk||), H2(||xk||), D1(||xk||), D2(||xk||) 
    H1_xk = inv8pi * (1/denom_xk + d2/denom_xk_3);
    H2_xk = inv8pi / denom_xk_3;
    D1_xk = inv4pi * (1/denom_xk_3 - 3*d2/denom_xk_5);
    D2_xk = -3 * inv4pi / denom_xk_5;

    % Derivatives of H1 and H2, divided by ||xk|| 
    r_xk_norm = sqrt(r_xk_sq);
    if r_xk_norm > 1e-15    %Avoid division by 0
        DH1_over_r = -inv8pi * (3 * d2 / denom_xk_5 + 1 / denom_xk_3);
        DH2_over_r = -3 * inv8pi / denom_xk_5;
    else
        DH1_over_r = 0;
        DH2_over_r = 0;
    end

    % g_k(f) = 2*f1*e1 - f = [f1; -f2; -f3]
    gk = [f(1); -f(2); -f(3)];
    
    % L_k(f) = cross(f, e1) x x_k 
    cross_Lk_xk = [f(2)*xk(2) + f(3)*xk(3); 
                  -f(2)*xk(1); 
                  -f(3)*xk(1)];

    % Scalar products
    dot_f_xk_star = f(1)*xk_star(1) + f(2)*xk_star(2) + f(3)*xk_star(3);
    dot_f_xk = f(1)*xk(1) + f(2)*xk(2) + f(3)*xk(3);
    dot_gk_xk = gk(1)*xk(1) + gk(2)*xk(2) + gk(3)*xk(3);

    % Final formula
    v = f * H1_star + xk_star * (dot_f_xk_star * H2_star) ...
        - (f * H1_xk + xk * (dot_f_xk * H2_xk)) ...
        - h2 * (gk * D1_xk + xk * (dot_gk_xk * D2_xk)) ...
        + (2 * h * (DH1_over_r + H2_xk)) * cross_Lk_xk ...
        + (2 * h) * (gk(1) * xk * H2_xk + ...
                     xk(1) * gk * H2_xk + ...
                     [dot_gk_xk * DH1_over_r; 0; 0] + ...
                     xk(1) * dot_gk_xk * xk * DH2_over_r);
    %safeguard to avoid floor penetration
     if x(1) <=10^-5 && v(1)<0
        v=[0;0;0];
    end
    
end