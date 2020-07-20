function [C, C_Y] = calculateCoriolisMat(M,q_th,dq_th,Q,n,tan_del)
    C = sym('C',[2*n,2*n,6,6]);
    C(:,:,:,:) = zeros(2*n,2*n,6,6);
    C_Y = sym('Cy',[14*n,14*n]);
    C_Y(:,:) = zeros(14*n,14*n);
    for i = 1:2*n
        for j = 1:2*n
           strcat("Calculating C: ", num2str(i),num2str(j))
            for k = 1:n
                [i,j,k]
                %C_dY(i,j,:,:,2*k-1,1) = .5*(diff(M(i,j,:,:),q(2*k-1,1))+diff(M(i,2*k-1,:,:),q(j,1))-diff(M(j,2*k-1,:,:),q(i,1)));
                diff_k = diff(M(i,j,:,:),q_th(k));
                if mod(j,2)==1
                    diff_j = diff(M(i,2*k-1,:,:),q_th((j+1)/2));
                else
                    diff_j = sym(0);
                end
                if mod(i,2)==1
                    diff_i = -diff(M(j,2*k-1,:,:),q_th((i+1)/2));
                else
                    diff_i = sym(0);
                end
                C(i,j,:,:) =  C(i,j,:,:) + .5*(diff_k + diff_i + diff_j)*dq_th(k);
                %C_Y(i,j,:,:,2*k-1,1) = diff(C(i,j,:,:),q(k*2-1,1));
                %for s = 1:6
                %    C_Y(i,j,:,:,2*k,s) = diff(C_dY(i,j,:,:,2*k-1,1),q(k*2,s));
                %end
            end
        end
    end
    for i = 1:n
        C(2*i,2*i,:,:) = C(2*i,2*i,:,:) + permute(tan_del*eye(6),[3,4,1,2]);
    end
    C_flat = flattenDynamicMat(permute(C,[3,4,1,2]),n,true);
    dQ_flat = Q(1:7*n);
    C_Y(7*n+1:end,1:7*n) = eye(7*n); %LL
    C_U = sym(zeros(7*n,14*n));
    parfor i = 1:14*n
        i
        C_U(:,i) = diff(C_flat*dQ_flat,Q(i));
    end
    C_Y(1:7*n,:) = C_U;
end