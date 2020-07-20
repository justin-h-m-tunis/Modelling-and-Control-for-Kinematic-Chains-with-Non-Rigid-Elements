function [M_fun,C_fun,N_fun,T_fun] = createDynamicFuns(M,C,N,T,n,q_th,q_del,dq_th,dq_del,K_params,T_th,get_M,get_C,get_N,get_T,save_funs)
    if get_M
    M_fun = {}
    end
    if get_C
    C_fun = {}
    end
    if get_N
    N_fun = {}
    end
    if save_funs
        if get_M
        for i = 1:2*n
            for j = 1:2*n
                M_file = strcat('Mass_fun_',num2str(i),num2str(j))
                M_fun{i,j} = matlabFunction(permute(M(i,j,:,:),[3,4,1,2]), 'Vars', {q_th,q_del},'File',M_file);
            end
        end
        end
        if get_C
        for i = 1:2*n
            for j = 1:2*n
                C_file = strcat('Coriolis_fun_',num2str(i),num2str(j))
                C_fun{i,j} = matlabFunction(permute(C(i,j,:,:),[3,4,1,2]), 'Vars', {q_th,q_del,dq_th,dq_del},'File',C_file);
            end
        end
        end
        if get_N
        for i = 1:2*n
            N_file = strcat('N_fun_',num2str(i))
            N_fun{i} = matlabFunction(N(i,:), 'Vars', {q_th,q_del,K_params},'File',N_file);
        end
        end
        if get_T
        T_file = 'T_fun'
        T_fun = matlabFunction(T, 'Vars', {T_th},'File','Force_mat','File',T_file);
        end
    else
        if get_M
        for i = 1:2*n
            for j = 1:2*n
                M_fun{i,j} = matlabFunction(permute(M(i,j,:,:),[3,4,1,2]), 'Vars', {q_th,q_del});
            end
        end
        end
        if get_C
        for i = 1:2*n
            for j = 1:2*n
                C_fun{i,j} = matlabFunction(permute(C(i,j,:,:),[3,4,1,2]), 'Vars', {q_th,q_del,dq_th,dq_del});
            end
        end
        end
        if get_N
        for i = 1:2*n
            N_fun{i} = matlabFunction(N(i,:), 'Vars', {q_th,q_del,K_params});
        end
        end
        if get_T
        T_fun = matlabFunction(T, 'Vars', {T_th});
        end
    end
end