function [bmin,phimin] = B4minSharpPath(phi,A4,p0,info)

if ~exist('info','var')
    info = 0;
end

if ~iscell(A4)
    % Preparations
    [idx,Lv,l,l0,p,mu_P,mu_Q,c_Q] = B4Preparations(A4,p0,info);

    % Set b1range such that nullspace special cases are avoided
    if mu_Q > 0
        b1min = 1.01*mu_Q;
    else
        b1min = 0.99*mu_Q;
    end
    b1max = b1min + (Lv(1) - Lv(end))/3;
    b1range = [b1min,b1max];

    % b1-path on boundary of B2
    n_b1 = 1000;
    b = B2Path(idx,Lv,p,mu_P,mu_Q,c_Q,b1range,n_b1);

    % Extend for b1-b4
    b = [b,zeros(n_b1,2)];

    % Compute corresponding b4 and b3
    for i=1:n_b1
        nv = p./(b(i,1)-Lv);
        b(i,4) = nv'*l/(nv'*p);
        c = b(i,4)*p-l;
        b(i,3) = l0+c'*pinv(diag(b(i,1)-Lv)+b(i,2)*p*p')*c;
    end

    % Evaluate function phi along b1-path
    phiv = zeros(n_b1,1);
    for i=1:n_b1
        phiv(i) = phi(b(i,:));
    end
    [phimin,imin] = min(phiv);
    bmin = b(imin,:);
else
    % Length of matrix set
    N = length(A4);

    % Preparations
    [idx,Lv,l,l0,p,mu_P,mu_Q,c_Q] = B4Preparations(A4,p0,info);

    % Set b1range such that nullspace special cases are avoided
    [mu_Q_max,imax] = max(cell2mat(mu_Q));
    b1min = mu_Q_max + 0.01*abs(mu_Q_max);
    b1max = b1min + (Lv{imax}(1) - Lv{imax}(end))/3;
    b1range = [b1min,b1max];

    % b1-path on boundary of B2 of A2 set
    n_b1 = 1000;
    b = B2Path(idx,Lv,p,mu_P,mu_Q,c_Q,b1range,n_b1);
    b = [b,zeros(n_b1,2)];

    % Search along path, shift and stop if new minimum increases
    phimin0 = Inf;
    check = 1;
    delta = 0.01;
    while check
        % Shift b1-path into interior of B2 of A2 set
        b(:,1:2) = b(:,1:2)+delta*ones(n_b1,2);

        % Lemma 5: (b1,b2) in interior of B2, b4 free.
        % Set b4=0, compute sharp b3 for every matrix and choose max
        for i=1:n_b1
            b3temp = zeros(N,1);
            for j=1:N
                c = -l{j};
                b3temp(j) = l0{j} + c'*((diag(b(i,1)-Lv{j}) + b(i,2)*p{j}*p{j}')\c); 
            end
            b(i,3) = max(b3temp);
        end

        % Evaluate function phi along b1-path
        phiv = zeros(n_b1,1);
        for i=1:n_b1
            phiv(i) = phi(b(i,:));
        end
        [phimin,imin] = min(phiv);
        bmin = b(imin,:);

        % Check
        check = phimin < phimin0;
        phimin0 = phimin;
        delta = 2*delta;
    end

end