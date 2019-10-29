% This function calculates the group-utility of all the groups in the data.Each 
% group is identified by an integer group-ID in a vector grpid of the same
% length as the number of columns in the data. Each column of the data
% should correspond to  group.
% Groups can be of any size. A group need not contain consecutive columns.
%
% Inputs : Rxxinv - (size: M X M) matrix. Inverse of the correlation matrix 
%                   (A'*A) where A is the (T X M) data matrix.
%          
%          W - (size: M X 1) MMSE filter estimate obtained as W = Rxxinv*(A'*b) where b is
%          the desired vector (T X 1) being reconstructed.
%
%          grpid - (size: M X 1) Vector of the same length number of
%          columns of data with each entry containing the group ID that the
%          corresponding column in A belongs to.
function util = grputilcalc (RXXinv, W, grpid)

    uniqgrps = unique(grpid);
    util = zeros(length(uniqgrps),1);
    for k = 1:length(uniqgrps)
        indx = (grpid==uniqgrps(k));
        
        % Select decoder weights of group k 
        Wkq = W(indx);
        
        % utility computation of group k
        S = RXXinv(indx,indx);
        
        util(k) = Wkq'*(S\Wkq);
    end
end