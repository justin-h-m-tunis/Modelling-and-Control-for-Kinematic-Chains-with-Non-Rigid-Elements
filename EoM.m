function dQ = EoM(t,Q,C,N,T,n,ind_th,ind_del,K)
    time=t
    %reformat Q:
        %Q = [q_th ; q_del ; dq_th ; dq_del]
    [q, dq, q_th, q_del, dq_th, dq_del] = Qtoq(Q,n,ind_th,ind_del);
    %T = T(Q)
    fcelleval = @(F,args)feval(F{1},args{:});
    C_args = {q_th,q_del,dq_th,dq_del};
    q_del
    C_eval = permute(cell2mat(arrayfun(@(F) permute(fcelleval(F,C_args),[3,4,2,1]), C,'UniformOutput',false)),[3,4,1,2]);
    %C_eval
    N_args = {q_th,q_del};
    N_eval = cell2mat(arrayfun(@(F)fcelleval(F,N_args),N,'UniformOutput',false))
    %N_eval
    T_eval = T(q_th,q_del,dq_th,dq_del);
    %(reshape might cause problems)
    dQ = -squeeze(sum(permute(dmatmul(dq,C_eval,true,[3 2],[1 2]),[3,4,1,2]),2)) - N_eval + T_eval;
    
    %dQ
    dQ = Qfromq(dq_th,dq_del,dQ(ind_th,1),dQ(ind_del,:),n);
    
    
    