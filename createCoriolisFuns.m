function C_fun = createCoriolisFuns(C, n,q_th,q_del,dq_th,dq_del,save_funs)
    C_fun = cell(2*n,2*n);
    if save_funs
        for i = 1:2*n
            for j = 1:2*n
                C_file = strcat('DynamicMats\Coriolis_fun_',num2str(i),num2str(j))
                C_fun{i,j} = matlabFunction(permute(C(i,j,:,:),[3,4,1,2]), 'Vars', {q_th,q_del,dq_th,dq_del},'File',C_file);
                ['C',num2str(i),num2str(j)]
            end
        end
    else
        for i = 1:2*n
            for j = 1:2*n
                C_fun{i,j} = matlabFunction(permute(C(i,j,:,:),[3,4,1,2]), 'Vars', {q_th,q_del,dq_th,dq_del});
                ['C',num2str(i),num2str(j)]
            end
        end
    end
    %{
for i = 1:2*n
        C_Y_fun{i,1} = matlabFunction(C_Y(:,:,:,:,i,1), 'Vars', {q_th,q_del,dq_th,dq_del});
        C__dY_fun{i,1} = matlabFunction(C_dY(:,:,:,:,i,1), 'Vars', {q_th,q_del,dq_th,dq_del});
        for j = 1:6
             C_Y_fun{i,j} = matlabFunction(C_Y(:,:,:,:,i,j), 'Vars', {q_th,q_del,dq_th,dq_del});
             C_dY_fun{i,j} = matlabFunction(C_dY(:,:,:,:,i,j), 'Vars', {q_th,q_del,dq_th,dq_del});
             ['CY',num2str(i),num2str(j)]
        end
    end
    %}
end