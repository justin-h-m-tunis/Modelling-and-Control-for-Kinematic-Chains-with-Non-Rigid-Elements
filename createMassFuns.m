function M_fun = createMassFuns(M,n,q_th,q_del,save_funs)
    M_fun =cell(2*n,2*n)
    if save_funs
        for i = 1:2*n
            for j = 1:2*n
                M_file = strcat('DynamicMats\Mass_fun_',num2str(i),num2str(j))
                M_fun{i,j} = matlabFunction(permute(M(i,j,:,:),[3,4,1,2]), 'Vars', {q_th,q_del},'File',M_file);
                ['M',num2str(i),num2str(j)]
            end
        end
    else
        for i = 1:2*n
            for j = 1:2*n
                M_fun{i,j} = matlabFunction(permute(M(i,j,:,:),[3,4,1,2]), 'Vars', {q_th,q_del});
                ['M',num2str(i),num2str(j)]
            end
        end
    end