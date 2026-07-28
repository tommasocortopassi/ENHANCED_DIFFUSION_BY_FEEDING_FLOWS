function P_fitted = fit_trajectories(P,is_tracer)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% FIT_TRAJECTORIES 
    % 
    % Fills missing values in particle or tracer trajectories.
    %
    %
    %   INPUTS
    %   ------
    %   P : struct
    %       Structure containing particle or tracer trajectories. Each field is
    %       expected to contain a 3-by-N matrix, where:
    %           - rows correspond to spatial coordinates,
    %           - columns correspond to successive time steps.
    %
    %       The fields must be named:
    %           p1, p2, ..., pK   if is_tracer == false
    %           t1, t2, ..., tK   if is_tracer == true
    %
    %   is_tracer : logical
    %       Flag indicating whether the structure contains tracer trajectories
    %       (true) or particles trajectories (false).
    %
    %   OUTPUT
    %   ------
    %   P_fitted : struct
    %       Structure with the same fields as the input, where NaN values have
    %       been replaced according to the following rules.
    %
    %   METHOD
    %   ------
    %   For each coordinate independently fill NaN segments in the coordinates:
    %
    %   1. Leading NaNs
    %      Replaced by the first previously available (non-NaN) value.
    %
    %   2. Trailing NaNs
    %      Replaced by the last available value.
    %
    %   3. Interior NaNs
    %      Filled by linear interpolation between the surrounding valid
    %      samples.
    %
    %   4. Entirely missing trajectory
    %      If all entries of a coordinate are NaN, the coordinate is replaced
    %      by zeros.
    %
    %   EXAMPLE OF OUTPUT
    %
    %   [NaN, NaN, 2, 4, NaN, NaN, 10, NaN, NaN] ---> [2, 2, 2, 4, 6, 8, 10, 10, 10]
    %
    %   NOTES
    %   -----
    %   - Consecutive NaN values are treated as a single missing interval.
    %   - Each spatial coordinate is processed independently.
    %   - The output structure preserves the original field names and layout.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Number of trajectories in the structure
    K=numel(fieldnames(P));
    
    % Start from the original structure
    P_fitted=P;
    
    % Determine number of time steps 
    if is_tracer==true
        N=size(P.t1,2);
    elseif is_tracer==false
        N=size(P.p1,2);
    end
    
    % Loop over all trajectories
    for i = 1:K
    
        % Determine the data type of P
        if is_tracer ==false
            name= "p"+num2str(i); 
        elseif is_tracer==true
            name= "t"+num2str(i); 
        end
    
        % Extract the current trajectory matrix
        A=P.(name);
    
        % Logical matrix marking NaN values
        NaNmatrix=isnan(A);
    
        % Number of columns/time steps
        M=size(A,2);
    
        for j=1:1:3
            % Find contiguous NaN blocks in this row
            
            d = diff([0, NaNmatrix(j,:), 0]);
            start_idx = find(d == 1); % beginning of each NaN block
            end_idx = find(d == -1);  % position right after each NaN block
            L = end_idx - start_idx;  % length of each NaN block
            
            % Loop over all NaN blocks in the current row
            for k = 1:length(L)
                idx_start = start_idx(k);
                idx_end = end_idx(k);
    
                % Case 1: the entire row is NaN
                if idx_start==1 && idx_end==M+1
                    P_fitted.(name)(j,:)=zeros(1,N);
    
                % Case 2: NaNs at the beginning of the row
                elseif idx_start==1
                    P_fitted.(name)(j,1:L(k))=ones(1,L(k))*A(j,idx_end);
                
                % Case 3: NaNs at the end of the row    
                elseif idx_end == M+1
                    P_fitted.(name)(j,idx_start:idx_end-1) = A(j,idx_start-1) * ones(1,1,L(k));
                % Case 4: Linear interpolation
                else
                    val_start = A(j,idx_start-1);
                    val_end = A(j,idx_end);
                    v = linspace(val_start, val_end, L(k)+2);
                    P_fitted.(name)(j,idx_start:idx_end-1) = v(2:end-1);  
                end
            end
        end
    end
end