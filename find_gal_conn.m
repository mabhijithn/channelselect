function Galconn_mat = find_gal_conn(chnl_list)
    Galconn_mat = zeros(size(chnl_list,1));
    for i = 1:size(chnl_list, 1)-1
        node_cand = chnl_list(i,:);
        for j = i+1:size(chnl_list,1)
            nxt_cand = chnl_list(j,:);
            uniq_elecs = unique([node_cand,nxt_cand]);
            if(length(uniq_elecs)~=4)
                Galconn_mat(i,j) = 1;
                Galconn_mat(j,i) = 1;
            end
        end
    end
end