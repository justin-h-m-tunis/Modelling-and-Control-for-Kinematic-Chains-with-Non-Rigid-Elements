function N_fun = createPotentialFuns(N,n,q_th,q_del,save_funs)
    N_fun = cell(2*n,1)
    %N_Y_fun = matlabFunction(N_Y(:,:), 'Vars',{q_th,q_del,K_params});
    %N_dY_fun = matlabFunction(N_dY(:,:), 'Vars',{q_th,q_del,K_params});
    if save_funs
        for i = 1:2*n
            N_file = strcat('DynamicMats\N_fun_',num2str(i))
            N_fun{i} = matlabFunction(N(i,:), 'Vars', {q_th,q_del},'File',N_file);
            ['N',num2str(i)]
        end
    else
        for i = 1:2*n
            N_fun{i} = matlabFunction(N(i,:), 'Vars', {q_th,q_del});
            ['N',num2str(i)]
        end
    end
end
