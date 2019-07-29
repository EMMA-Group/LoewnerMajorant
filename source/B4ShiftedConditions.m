function [cineq,ceq] = B4ShiftedConditions(b,Lv,l,l0,p)

if ~iscell(Lv)
    % Null vector n, vecor c and minimum b2 for given b1 ~= Lmax and b4
    n = p./(b(1)-Lv);
    c = b(4)*p - l;
    b2m = -1/(n'*p);

    % Tolerance
    tol = 1e-10; 

    % Equality ceq = 0 and inequality vector cineq <=0 (ineq for b2 and b3)
    ceq = [];
    cineq = [b2m-b(2); 1];

    % Check if on boundary of B2 or in interior of B2
    if abs(cineq(1)) < tol
        % On boundary of B2 (lemma 4)
        % Check b4 for b(1) >= mu_shited > mu_hat
        b4crit = n'*l/(n'*p);
        d = abs(b(4)-b4crit);
        if( d > tol )
            % b4 is inadmissible, set ineq for b3 as false
            fac = 1e6;
            cineq(2) = d * fac;
            % Attention: fac and 1/tol should not differ by more than
            % sqrt(accuracy_fmincon)
        else
            % b4 is admissible
            % evaluate inequality constraint for b3
            ci = pinv(diag(b(1)-Lv) + b(2)*(p*p'))*c;
            cineq(2) = l0 + c'*ci - b(3);
        end
    else
        % Not on boundary of B2 -> in interior of B2 (lemma 5)
        % Evaluate inequality constraint for b3
        ci  = (diag(b(1)-Lv) + b(2)*(p*p')) \ c;
        cineq(2)= l0 + c'*ci - b(3);
    end
else
    % Length of matrix set
    N = length(Lv);

    % For given b1 compute minimum b2 
    [b2m,im] = max(arrayfun(@(i) -1/sum(p{i}.^2./(b(1)-Lv{i})),1:N));
    % and corresponding null vector n and vecor c
    n = p{im}./(b(1)-Lv{im});
    c = b(4)*p{im}-l{im};

    % Equality ceq = 0 and inequality vector cineq <=0 (ineq for b2 and b3)
    ceq = [];
    cineq = [b2m-b(2); 1];

    % Tolerance
    tol = 1e-10;

    % Check if on boundary of B2 or in interior of B2
    if abs(cineq(1)) < tol
        % On boundary of B2 (lemma 4)
        b4c = n'*l/b2m;
        if abs(b(4)-b4c) < tol
            % b(4) accepted
            b3c = -Inf;
            for j=1:N
                c = b4c*p{j}-l{j};
                b3c = max(b3c,l0{j}+c'*pinv(diag(b(1)-Lv{j}) + b2m*p{j}*p{j}')*c);
            end
            % Evaluate condition for b(3)
            cineq(2) = b3c - b(3);
        end
    else
        % In interior of B2 (lemma 5)
        b3c = -Inf;
        for j=1:N
            c = b(4)*p{j}-l{j};
            ci = (diag(b(1)-Lv{j}) + b(2)*p{j}*p{j}')\c;
            b3c = max(b3c,l0{j}+c'*ci);
        end
        % Evaluate condition for b(3)
        cineq(2) = b3c - b(3);
    end
end

end
