function [idx,Lv,l,l0,p,mu_P,mu_Q,c_Q] = B4Preparations(A4,p0,info)

if ~exist('info','var')
    info = 0;
end

if ~iscell(A4)
    % Dimension
    n = length(p0);

    % Split the matrix
    A = A4(1:n,1:n);
    a = A4(1:n,n+1);
    a0 = A4(n+1,n+1);

    % Get auxiliary quantities of nxn matrix A
    [idx,Lv,p,mu_P,mu_Q,c_Q,U] = B2Preparations(A,p0,info);

    % Transform a
    l = U'*a;
    l0 = a0;
else
    % Check input errors
    n = length(p0);
    N = length(A4);
    tol = 1e-14;
    if ~iscolumn(p0)
        error('Vector p0 must be a column vector');
    end
    if abs(norm(p0)-1)>tol
        error('Vector p0 must be normalized');
    end
    if any(~cellfun(@(A) size(A,1) == n+1 && size(A,2) == n+1, A4))
        error('All matrices in A must have dimensions [n+1,n+1]');
    end

    % Get auxiliary quantities for each matrix
    idx = cell(N,1);
    Lv = cell(N,1);
    l = cell(N,1);
    l0 = cell(N,1);
    p = cell(N,1);
    mu_P = cell(N,1);
    mu_Q = cell(N,1);
    c_Q = cell(N,1);
    for i=1:N
        [idx{i},Lv{i},l{i},l0{i},p{i},mu_P{i},mu_Q{i},c_Q{i}] = B4Preparations(A4{i},p0,info);
    end
end

end
