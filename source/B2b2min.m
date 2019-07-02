function out = B2b2min(b1,idx,Lv,p,mu_P,mu_Q,c_Q)

% Eigenvalues and multiplicity of maximum eigenvalue lmax
L1 = Lv(1);
L2 = Lv(2);
m_L1 = length(find(Lv==L1));

% Lemma 1
if idx == 1
    if b1>=mu_Q
        out = mu_P - b1;
    else
        error('b1 value smaller than mu_Q is not admissible.');
    end
end

% Lemma 2
if idx == 2
    % mu_Q = L1
    if mu_Q == L1
        if b1 > L1
            out = -1/sum((p.^2)./(b1 - Lv));
        end
        if b1 == L1
            p_1 = p(1:m_L1);
            if ~any(p_1)
                out = -1/c_Q;
            else
                out = 0;
            end
        end
        if b1 < L1
            error('Case2-L1: b1 < L1 is not admissible');
        end
    end
    % mu_Q = L2 < L1
    if and(L2<L1,mu_Q==L2)
        if (L2<b1)&&(b1~=L1)
            out = -1/sum((p.^2)./(b1 - Lv));
        end
        if b1==L1
            out = 0;
        end
        if b1==L2 && c_Q~=0
            out = -1/c_Q;
        end
        if b1==L2 && c_Q==0
            error('Case2-L2: c_Q=0, i.e., b1=L2 not admissible.');
        end
        if b1<L2
            error('Case2-L2: b1 not admissible');
        end
    end
end

% Lemma 3
if idx == 3
    if b1 > mu_Q && b1 ~= L1
        out = -1/sum((p.^2)./(b1 - Lv));
    else
        if b1 == L1
            out = 0;
        else
            error('b1 <= mu_Q not admissible.');
        end
    end     
end

end
