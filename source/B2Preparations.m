function [idx,Lv,p,mu_P,mu_Q,c_Q,U] = B2Preparations(A,p0,info)

% Default
if ~exist('info','var')
    info = 0;
end

% A is a matrix (not a cell/set of matrices)
if not(iscell(A))
    % Dimension
    n = length(p0);

    % Check input errors
    if ~iscolumn(p0)
        error('Vector p0 must be a column vector');
    end
    if abs(norm(p0)-1)>1e-14
        error('Vector p0 must be normalized');
    end
    if any(size(A)~=[n,n])
        error('Matrix A must have dimensions [n,n]');
    end

    % Compute spectral decomposition and reorder
    [U,L] = eig(A);
    [Lv,order] = sort(diag(L),'descend');
    L = diag(Lv);
    U = U(:,order);
    L1 = Lv(1);

    % Compute p, mu_P and mu_Q
    p = U' * p0;
    Q = orth(eye(n)-p*p');
    mu_P = p'*L*p;
    mu_Q = max(eig(Q'*L*Q));
    c_Q = p'*pinv(diag(mu_Q-Lv))*p;

    % Tolerance for case identification
    tol = 1e-10;

    % Case identification
    if any(norm(Lv.*p0-mu_P*p0) < tol)
        % Case 1: p is an eigenvector
        idx = 1;
        if info==1
            fprintf('\nCase 1 (lemma 1): p0 is an eigenvector of A\n');
        end
    else
        L2 = Lv(2);
        if or(abs(mu_Q-L1) < tol,abs(mu_Q-L2) < tol)
            % Case 2: p is not an eigenvector and mu_Q \in {L1,L2}
            idx = 2;
            if abs(mu_Q-L1)<tol
                mu_Q = L1;
                if info==1
                    fprintf('\nCase 2 (lemma 2): p0 is NOT an eigenvector of A and mu_Q = L1\n');
                end
            else
                mu_Q = Lv(2);
                if info==1
                    fprintf('\nCase 2 (lemma 2): p0 is NOT an eigenvector of A and mu_Q = L2\n');
                    if c_Q == 0
                        disp('Attention: c_Q = 0, i.e., mu_Q is NOT admissible')
                    end
                end
            end
        else
            % Case 3: p is not an eigenvector and L2 < mu_Q < L1
            idx = 3;
            if info==1
                fprintf('\nCase 3 (lemma 3): p0 is NOT eigenvector of A and L2 < mu_Q < L1\n');
            end
        end
    end
    
    % Further information
    if info==1
        fprintf('mu_P = \t %.4e\nmu_Q = \t %.4e\nL1 = \t %.4e\nL2 = \t %.4e\n',mu_P,mu_Q,Lv(1),Lv(2));
    end
end

% A is a set of matrices/cell
if iscell(A)
    N = length(A);
    idx = cell(N,1);
    Lv = cell(N,1);
    p = cell(N,1);
    mu_P = cell(N,1);
    mu_Q = cell(N,1);
    c_Q = cell(N,1);
    for i=1:N
        [idx{i},Lv{i},p{i},mu_P{i},mu_Q{i},c_Q{i}] = B2Preparations(A{i},p0,info);
    end
end

end
