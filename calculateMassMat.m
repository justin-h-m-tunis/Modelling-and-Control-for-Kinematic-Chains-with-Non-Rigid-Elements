function M = calculateMassMat(J_si,m,n)
M = sym('M',[2*n,2*n,6,6]);
    M(:,:,:,:) = zeros(2*n,2*n,6,6);
    for i = 1:2*n
        MJ = dmatmul(m(:,:,i),J_si(:,:,:,i),true,[1,2],[1,3]);%6x6xn
        parfor j = 1:2*n %multiply J_6x6,j,i matrices across MJ_6x6,j,i - makes 6x6xn for n rows
            strcat("Calculating M: ", num2str(i),num2str(j))
            M(j,:,:,:) = M(j,:,:,:) + permute(dmatmul(J_si(:,j,:,i),MJ,true,[3,1],[1,2]),[4,3,1,2]);
            "done"
        end
    end    
end