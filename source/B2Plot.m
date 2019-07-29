function B2Plot(A,p0,info,nc,b1range,b2range,n_b1,ls,col)
% B2Plot(A,p0,info,nc,b1range,b2range,n_b1,ls)
% plots the boundary of the B2 set 
%
% Inputs:
% A           nxn symmetric matrix
% p0          n-dimensional normalized column vector
% info        [optional] set to 1 if information is to be displayed
% nc          [optional] set to 1 if necessary conditions are to be plotted
% b1range     [optional] set plot range [b1min,b1max] for b1
% b2range     [optional] set plot range [b2min,b2max] for b2
% n_b1        [optional] number of b1-points
% ls          [optional] line style for plotting
% col         [optional] color
% 
% Output:
% Plot of B2 set.

%% Body

% Set defaults 
if ~exist('info','var')
    % Default: no information displayed
    info = 0;
end
if ~exist('nc','var')
    % Default: necessary conditions not displayed
    nc = 0;
end
if ~exist('ls','var')
    % Default: standard line plot
    ls = '-';
end
if ~exist('col','var')
    % Default: standard plot color string
    col = 'standard';
end

if ~iscell(A)
    % Preparations
    [idx,Lv,p,mu_P,mu_Q,c_Q] = B2Preparations(A,p0,info);

    % Set defaults based on preparations
    if ~exist('b1range','var')
        % Default: set b1range accordingly
        b1min = mu_Q;
        if (idx==3)||(idx==2 && c_Q==0)
            b1min = mu_Q + 0.01*abs(mu_Q);
        end
        b1max = b1min + 1;
        b1range = [b1min,b1max];
    end
    if ~exist('n_b1','var')
        % Default: discretize in at least 200 points
        n_b1 = 200;
    end

    % Lower boundary of B2 set for b1 in [b1min,b1max] 
    b = B2Path(idx,Lv,p,mu_P,mu_Q,c_Q,b1range,n_b1);

    % Set defaults based on boundary of B2
    if ~exist('b2range','var')
        % Default: b2 plot range
        b2range = [min(b(:,2)),1e3];
    end

    % Plot
    if strcmp(col,'standard')
        plot(b(:,1),b(:,2),'linewidth',2,'linestyle',ls);
    else
        plot(b(:,1),b(:,2),'linewidth',2,'linestyle',ls,'color',col);
    end

    % Add necessary conditions, if desired
    if nc == 1
        line([mu_Q,mu_Q],b2range,'color','black','LineStyle','--');
        line(b1range,mu_P-b1range,'color','black','LineStyle','--');
    end

    % Plot limits
    xlim(b1range);
    ylim(b2range);

    % Labels
    xlabel('\beta_1');
    ylabel('\beta_2');
else
    % Preparations for matrix set A
    [idx,Lv,p,mu_P,mu_Q,c_Q] = B2Preparations(A,p0);

    % Determine most restrictive mu_hat
    [mu_Q_max,i_max] = max(cell2mat(mu_Q));
    
    % Set defaults
    if ~exist('b1range','var')
        b1min = mu_Q_max;
        if (idx{i_max}==3)||(idx{i_max}==2 && c_Q{i_max}==0)
            b1min = mu_Q_max + 0.01*abs(mu_Q_max);
        end
        b1max = b1min + 1;
        b1range = [b1min,b1max];
    end
    if ~exist('n_b1','var')
        n_b1 = 100;
    end
    
    % Lower border of B2 set for b1 in [b1min,b1max] 
    bSet = B2Path(idx,Lv,p,mu_P,mu_Q,c_Q,b1range,n_b1);
    
    % Set range for b2, if not specified
    if ~exist('b2range','var')
        b2range = [min(bSet(:,2)),1e3];
    end
    
    % Plot
    hold on;
    % Plot B2 for every single matrix
    for i=1:length(A)
        B2Plot(A{i},p0,info,nc,b1range,b2range,100,'--');
    end
    % Plot B2 for matri set
    plot(bSet(:,1),bSet(:,2),'linewidth',2,'color','black');

    % Plot limits
    xlim(b1range);
    ylim(b2range);

    % Labels
    xlabel('\beta_1');
    ylabel('\beta_2');
    
    hold off;
end