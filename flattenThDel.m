function flat = flattenThDel(q_th,q_del,n,symb)
    if nargin==3 || symb
        flat = sym(zeros(7*n,1));
    else 
        flat = zeros(7*n,1);
    end
    flat(mod(1:7*n,7)==1) = q_th;
    flat(mod(1:7*n,7)~=1) = q_del.';
end