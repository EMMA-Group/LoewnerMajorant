function b = B2Path(idx,Lv,p,mu_P,mu_Q,c_Q,b1range,n_b1)

if ~iscell(idx)
    % Shifted mu_Q
    mu_Q_s = mu_Q + 0.01*abs(mu_Q);
    
    % Set defaults
    if ~exist('b1range','var')
        b1min = mu_Q;
        if (idx==3)||((idx==2) && (c_Q==0))
            b1min = mu_Q_s;
        end
        b1max = b1min + 1;
    else
        b1min = max(mu_Q,b1range(1));
        if (idx==3)||((idx==2) && (c_Q==0))
            b1min = max(mu_Q_s,b1range(1));
        end
        b1max = max(b1min+0.001,b1range(2));
    end
    if ~exist('n_b1','var')
        n_b1 = 200;
    end

    % Container
    b = zeros(n_b1,2);

    % Discretization
    b(:,1) = linspace(b1min,b1max,n_b1);
    b(:,2) = arrayfun(@(b1) B2b2min(b1,idx,Lv,p,mu_P,mu_Q,c_Q),b(:,1));
    
    % Add point on left border for cases 1 and 2 (with c_Q~=0) if b1min=mu_Q
    if (idx<3)&&(c_Q~=0)&&(b1min==mu_Q)
        b = [[mu_Q,1e5];b];
    end
else
    % Number of matrices in set
    N = length(idx);

    % Determine most restrictive mu_Q
    [mu_Q_max,i_max] = max(cell2mat(mu_Q));
    
    % Shifted mu_Q_s
    mu_Q_s = mu_Q_max+0.01*abs(mu_Q_max);
    
    % Set defaults for optional parameters
    if ~exist('b1range','var')
        b1range = [mu_Q_max,mu_Q_max + 1];
    end
    if ~exist('n_b1','var')
        n_b1 = 100;
    end
    
    % Set b1min and b1max accordingly
    b1min = mu_Q_max;
    if (idx{i_max} == 3)||(idx{i_max}==2 && c_Q{i_max} == 0)
        b1min = mu_Q_s;
    end
    if exist('b1range','var')
        b1min = max(b1min,b1range(1));
    end
    b1max = b1min + 1;
    if exist('b1range','var')
        b1max = max(b1max,b1range(2));
    end
    
    % Container for path
    b = zeros(n_b1,2);

    % Determine path along b1 discretization based on most restrictive b2min,
    % i.e., for each b1 take the max of all b2min
    b(:,1) = linspace(b1min,b1max,n_b1);
    i_max = zeros(n_b1,1);
    b2temp = zeros(N,1);
    for i=1:n_b1
        for j=1:N
            b2temp(j) = B2b2min(b(i,1),idx{j},Lv{j},p{j},mu_P{j},mu_Q{j},c_Q{j});
        end
        [b(i,2),i_max(i)] = max(b2temp);
    end
end