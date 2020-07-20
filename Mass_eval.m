function M = Mass_eval(t,Q,M_fun,n,ind_th,ind_del)
    [q,dq, q_th, q_del, dq_th, dq_del] = Qtoq(Q,n,ind_th,ind_del);
    fcelleval = @(F,args)feval(F{1},args{:});
    M_args = {q_th,q_del};
    M_eval = permute(cell2mat(arrayfun(@(F)permute(fcelleval(F,M_args),[3,4,1,2]),M_fun, 'UniformOutput',false)),[3,4,1,2]);
    M = zeros(2*7*n,2*7*n);
    %dq part
    M(7*n+1:end,7*n+1:end) = eye(7*n);
    %d2q part
    M(1:7*n,1:7*n) = flattenDynamicMat(M_eval,n);
    %pause()
 end
