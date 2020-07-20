function J_si = calculateJacobian(e_xi,xi_th,xi_del,g_0,n,simple)
    %Ji
    J_si = sym('j',[6,2*n,6,2*n]);
    J_si(:,:,:,:) = zeros(6,2*n,6,2*n);
     for i = 1:2*n %dagga
        for j = 1:i %col
            g_ji = sym('gji',[4,4]);
            g_ji(:,:) = eye(4);
            for k = j:i
                g_ji = g_ji*e_xi(:,:,k);
            end 
            g_ji = g_ji*g_0(:,:,i);
            R = g_ji(1:3,1:3);
            p = g_ji(1:3,4);
            Adg_inv = [R.' -R.'*auto_skew(p); zeros(3) R.'];
            if mod(j,2)==1
                J_si(:,j,1,i) = Adg_inv*xi_th((j+1)/2,:).';
            else
                J_si(:,j,:,i) = Adg_inv*[zeros(6,1),squeeze(xi_del(j/2,:,2:6))];
            end
        end
     end
     if simple
        J_si = simplify(J_si);
     end
end