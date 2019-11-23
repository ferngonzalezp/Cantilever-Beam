function [ind] = index(conect,dim, dof_list)
    ind = zeros(1,length(conect)*dim);
    for j = 1:length(conect)
        for k = 2:size(dof_list,2)
                ind(dim*j-(dim+1-k)) = dof_list(conect(j),k);
        end
    end
end