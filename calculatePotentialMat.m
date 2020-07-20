function [N, N_Y, N_dY] = calculatePotentialMat(g_si,J_si,q_th,q_del,m,K,n)
    N = sym('N',[2*n,6]);
    N(:,:) = zeros(2*n,6);
    N_Y = sym('Ny',[7*n,7*n]);
    N_Y(:,:) = zeros(7*n,7*n);
    N_dY = sym('Ndy',[7*n,7*n]);
    N_dY(:,:) = zeros(7*n,7*n);
    syms V;
    V = sum(m.*permute(g_si(2,4,:),[3,1,2]).*9.81);
    for i = 1:n
        V = V + .5*q_del(i,:)*K(:,:,i)*q_del(i,:).';
    end
    for i = 1:n
        strcat("Calculating N: ", num2str(i))
        N(2*i-1,1) = diff(V,q_th(i));
        %N(2*i,:) = K(:,:,i)*q_del(i,:).';
        for j = 1:6
            N(2*i,j) = N(2*i,j) + diff(V,q_del(i,j));
        end
        %N(:,:) = N(:,:) + squeeze(dmatmul(K(:,:,i)*q_del(i,:).',J_si(:,:,:,2*i-1),true,[2,1],[1,2])); %have to project this onto axis of rotation... J.'? ... diff J.' K J?
    end
    N_ = flattenStateMat(N,n);
    q_ = flattenThDel(q_th,q_del,n);
    for i = 1:7*n
        for j = 1:7*n
            N_Y(i,j) = diff(N_(i),q_(j));
        end
    end
end


%{
del_true = g_si(:,:,2*i)/g_si(:,:,2*i-1);
        theta_rotation = acos(trace(del_true(1:3,1:3)-1)/2);
        if theta_rotation == 0 %might be an issue later
            w = [0;0;0];
        else
            w = 1/sin(theta_rotation).*[del_true(3,2) - del_true(2,3); del_true(1,3) - del_true(3,1); del_true(2,1) - del_true(1,2)];
        end
        g_del = [del_true(1:3,4);w];
        V = V + g_del.'*K(:,:,i)*g_del; % might want dot product
%}