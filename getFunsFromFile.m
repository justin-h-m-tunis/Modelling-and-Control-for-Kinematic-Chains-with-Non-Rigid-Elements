function [M_fun,C_fun,N_fun,T_fun] = getFunsFromFile(n)
    M_fun = cell(2*n,2*n)
    C_fun = cell(2*n,2*n)
    N_fun = cell(2*n,1)
    for i = 1:2*n
            for j = 1:2*n
                M_file = strcat('Mass_fun_',num2str(i),num2str(j))
                M_fun{i,j} = str2func(M_file)
                C_file = strcat('Coriolis_fun_',num2str(i),num2str(j))
                C_fun{i,j} = str2func(C_file)
            end
            N_file = strcat('N_fun_',num2str(i))
            N_fun{i} = str2func(N_file)
    end
    T_file = 'T_fun'
    T_fun = str2func(T_file)
    %Y_file = 'DynamicMats\Yacobian_fun'
    %Y_fun = str2func(Y_file)
end
