function el_list = find_el_from_node(node_list,connec)
%FIND_EL_FROM_NODE Summary of this function goes here
%   Detailed explanation goes here
NNE = size(connec,1);

el_list = cell(size(node_list,1),1);
for i=1:size(node_list,1)
    node = node_list(i,:);
    
    el1 = my_mod(find(connec == node(1)),NNE);
    el2 = my_mod(find(connec == node(2)),NNE);
    el_list{i} = intersect(el1,el2);
end

    function res  = my_mod(X,Y)
        res = mod(X,Y);
        res(res==0)=Y;
    end

end

