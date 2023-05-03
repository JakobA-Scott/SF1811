%% Laboration 1 SF1811, Alexander Råberg och Jakob Amaya

function f = simplex(A, b, c, Beta, v)

% uppgift 1 SF1811
%max C^T*x
%subject to Ax <= b, x> 0 
%A = [1 1 0 -1 0 0; 1 0 1 0 -1 0; 0 1 1 0 0 -1];
%b = [2;2;2];
%c =[-1 -1 -4 0 0 0]';  
%Beta = [1 2 3];
%v = [4 5 6];

nytt_c = -c;

Loop = true;

iterations = 1;

while Loop == true
    %Rewrite into standard form 
    A_Beta = A(:,Beta);
    C_Beta = nytt_c(Beta);
    A_v = A(:,v);
    c_v = nytt_c(v);    
    b_bar = A_Beta^(-1)*b;
    
    %Check if current solution is a feasible solution
    if min(b_bar) < 0 
        disp('This corresponds to a solution that is not feasible.')
        break
    end
    
    
    y = A_Beta'^-1 * C_Beta;
    rv = c_v - A_v'*y;
    
    %Condition to check if there exists a better solution 
    if rv(1:end) >= 0
        disp('An optimal solution has already been obtained. Current solution is the optimal solution')
        solution = zeros(1, length(Beta)+length(v));
        solution(Beta) = b_bar;
        disp('The optimal solution is')
        cost = solution*nytt_c
        Loop = false;
        break


    else 
        min_rv_index = find(rv == min(rv));

        %Smallest value in the vector v, that will be introduced to beta in next iteration 
        old_v = v(min_rv_index); %Index for v value we get rid off
        
        %Values that we keep in v
        remaning_v = v(v ~= old_v);
        
        %A and a_bar   
        a = A(:, old_v);
        A_bar = A_Beta^(-1)*a;
        comparison = b_bar./A_bar;

        %To find the element in Beta that we are going to take out and
        %replace

        min_value = min(comparison);
        if min_value <= 0     
            dummy = find(comparison ~= min_value);
            min_value = comparison(dummy);
        end
        no_solution = find((comparison ~= min(min_value)));
        
        if length(no_solution) == length(comparison)
            disp('There is unfortunately no solution to this problem')
            break
        end        
        

        index_for_new_v = find(comparison == min_value);
        value_for_leaving_beta = Beta(index_for_new_v);
        

        %Values that we keep in beta
        remaining_beta = Beta(Beta ~= value_for_leaving_beta);

        %New beta vector
        new_BFS = [old_v, remaining_beta];
        Beta = new_BFS;
        Beta = sort(Beta);
        
        %New v vector
        new_v = [remaning_v, value_for_leaving_beta];
        v = new_v;
        v = sort(v);       
    end
end

end