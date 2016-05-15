%% Completion Helper Method
% Arguments: M = input matrix
%            M_shdow = 1 & 0 matrix=> 1 = known number
%                                     2 = unknown number
%            step_size = round size. e.g. 2 means save 2 significant number
%                        Larger number means more accurate, but each step
%                        will be small ,so need more time. Therefore, we
%                        can use small number at the start, and then use
%                        larger number to make the matrix more accurate.
% Return: A = output matrix
function A = completionHelperRough(M, M_shadow, step_size)
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
    
    %% loop invariant
    % Assume matrix A holds thebest  matrix result and the prev_result
    % holds the previous best result.
    
    % base case: Since no loop is executed, A is the best current result,
    % the loop invarinat holds.
    
    % Inductive step: For each loop, we save the nuclear norm of A to 
    % nuclear_norm and save the A to prev_result before changing A. After 
    % this, decrease the nuclear norm of A to produce a new A. This new A
    % has a smaller nuclear norm compared with previous one. Therefore, A
    % holds the best current result and prev_result holds the previous best 
    % result. Loop invariant holds.
    
    % Termination: Since each step we try to decrease the nuclear norm of A
    % by the step_size (like gradient descent?), the step_size is a rough 
    % number and the new A may get worse. Therefore, we the loop
    % terminates, prev_result holds the best one and A holds a worse one.
    % Then, A = prev_result;
    while sum(diag(S))< nuclear_norm
        nuclear_norm = sum(diag(S)); %update nuclear norm
        prev_result = A; % save A before changing A
        
       %% shink singular value
        min_sigma = 0;
        i = min(size(S));
        while i > 0
            % check the round value of S(i,i), and also round the S(i,i).
            % Thus, the round of S(i,i) ¡Ù S(i,i)
            S(i,i) = round(S(i,i), step_size);
            if(S(i,i)>0)
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