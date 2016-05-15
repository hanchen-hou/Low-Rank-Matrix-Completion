%% Completion Helper Accurate
% The procedure is the same as completionHelperRough, but more accurate
% If the program just use this method, matrix recovery will be slow.
% Use Rough method to approach the correct solution, then use Accurate to
% make the corrections.
function A = completionHelperAccurate(M, M_shadow, step_size)
    % make a copy of M, keep the orgin (M)atrix, change A
    A = M;
    % the result A may get worse, need to save the prev
    prev_result = A; 
    % m = row, n = col
    [m,n] = size(A);

    % initial a huge nuclear norm
    nuclear_norm = 1e100; 
    % S is the diagonal matrix of A, initial value is 0
    [U,S,V] = svd(A);
    % Real loop time, for debug and loop watch
    real_loop_time = 0;
    
    while sum(diag(S))< nuclear_norm
        nuclear_norm = sum(diag(S)); %update nuclear norm
        prev_result = A; % save A before changing A
        
       %% shink singular value
        min_sigma = 0;
        i = min(size(S));
        while i > 0
            % check the round value of S(i,i), but keep the original value 
            % of S(i,i).
            if(round(S(i,i), step_size)>0)
                min_sigma = S(i,i);
                break;
            else
                % Because the matlab may keep the small number each step,
                % and these small entries in the S can never be 0 by
                % calculating, these tiny number should be ignored
                % manually.
                S(i,i) = 0; 
            end
            i = i - 1;
        end
        while i > 0
            S(i,i) = S(i,i) - min_sigma;
            if(S(i,i)<0)
                S(i,i) = 0.0;
            end
            i = i - 1;
        end
        
       %% calculate the new result A
        temp_M = U * S * V';
        if(temp_M == 0)
            %something wrong?
            disp('temp_M is zero!');
            break;
        else
            for i = 1:m
                for j = 1:n
                    if(M_shadow(i,j) == 0)
                        A(i,j) = temp_M(i,j);
                    end
                end
            end
        end 
        [U,S,V] = svd(A);
        real_loop_time = real_loop_time + 1;
    end %loop terminate
    
    A = prev_result; %current result A is worse, go back to prev_result
    disp(['Keep ', num2str(step_size),' significant number for each shink.']);
    disp(['Current best nuclear norm: ', num2str(nuclear_norm,20)]);
	disp(['Iteration time: ',num2str(real_loop_time)]);
end