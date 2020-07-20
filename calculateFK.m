function [g_si, g_si_th, g_si_del] = calculateFK(e_xi,g_0,n, ind_th, ind_del)
    g_si = sym('gs',[4,4,2*n]);
    for i = 1:2*n
        g_si(:,:,i) = eye(4);
        for j = 1:i
            g_si(:,:,i) = g_si(:,:,i)*e_xi(:,:,j);
        end 
        g_si(:,:,i) = g_si(:,:,i)*g_0(:,:,i);
    end
    g_si_th = g_si(:,:,ind_th);
    g_si_del = g_si(:,:,ind_del);
    