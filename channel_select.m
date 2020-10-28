function ch_selected = channel_select(A, b, N, varargin)

% channel_select A function which selects the best N channels of A for 
%                the least-squared problem min_{x} 1/2 * ||Ax - b ||² 
%                which best estimates b.
%
%  ch_selected = channel_select(A, b, N, 'method', 'utility', 'group_size', group_size)
%                Selects the best N channels of A using greedy 
%                utility-based channel selection. Here the columns of A are
%                assumed to contain groups of group_size channels in consecutive columns.
%
%  ch_selected = channel_select(A, b, N, 'method', 'utility', 'grpid',  grpid)
%                Selects the best N group-IDs of A using greedy 
%                utility-based channel selection. Here the columns of A are
%                assumed to contain groups of channels in non-consecutive columns.
%                Each column of A belongs to a group corresponding to its index 
%                in grpid.
% Required Inputs:
%   A - Data/measurements used in estimation where P is the number of 
%       channels (including delayed versions) and T is the number of
%       time samples. 
%       If A consists of delayed samples of a channel, all delays of 
%       a channel are expected to be consecutive columns
%       Size: (T X P)
%   b - The desired data to be estimated using A. Size: (T X 1)
%   N - Number of channels to be selected.
% Optional Inputs (Name-value pairs: Each optional input string is to be followed by an argument):
%   'method' - A string for the method to be used to select best N channels
%            Currently supported:
%            'utility' - Utility based channel reduction. Here, if A is 
%            singular, utility definition with minimum norm increase 
%            interpretation is used.
%            (ref. A. Bertrand, 2018)
%
%   'group_size' - Number of columns of A to be considered as a group.
%          For eg: Channel-k and its time-lagged copies can be together considered as
%          a group. The channels in the group should be in consecutive
%          columns.
%          The utility of this group is considered as the utility of channel-k.
%          If every channel and their 2 time-lagged versions (a channel itself + 2 lagged versions) are present, then:
%          (group size in the above case would be 3.)
%          default: 0 (no time-lagged versions included. group size is 1)
%  'grpid' - A vector of length equal to number of columns. The vector contains integer entries.
%	     Each integer entry corresponds to the ID of a column. All columns belonging to a group should correspond to same ID.
%	     Every column should have an ID.  This argument allows custom arrangement of columns belonging to same group.
%	     NOTE: This argument is ignored if 'group_size' input is present. 
%               
%   
%
% Outputs: 
%   ch_selected - The best N columns in A or the best N groups of A based 
%                 on group-IDs. 
%                 There is NO ORDER OF SIGNIFICANCE for this list.
%
% author: Abhijith Mundanad Narayanan
% e-mail: mailme@abhimundanad.com
%
%

    % Default conditions
    method = 'utility';
    noflags = 1;
    lagsflag = 1;
    if(nargin<3)
        N = size(A,2);
    end
    
    
    if nargin > 3        
        for i = 1:2:length(varargin)
            Param = varargin{i};
            Value = varargin{i+1};
            if ~isstr(Param)
                error('Flag arguments must be strings')
            end
            Param = lower(Param);
            switch Param
                case 'method'
                    method = Value;
		% The following case is for backward compatibility for an earlier version of this
		% code which had 'lags' argument instead of 'group_size'
                case 'lags'
                    noflags = Value+1;
                    lagsflag = 1;
                case 'group_size'
                    noflags = Value;
                    lagsflag = 1;
                case 'groupid'
                    grpid = Value;
                    nofuniq = unique(grpid); % Number of unique groups
                    lagsflag = 0;
                    if(size(grpid,1)~=size(A,2))
                        error('Length of group-ids != No.of columns in A');
                    end
            end
        end
    end
    
    if(strcmp(method,'utility'))
        % Calculate the auto and crosscovariances
        RXX = (A'*A)/size(A,1);
        RXY = (A'*b)/size(A,1);
        no_of_channels = size(RXX,2)/noflags;

        % Initialise a list of original indices/channels
        chnl_list = (1:no_of_channels)';

        % Indicator matrices to book-keep removed channels

        % selector of columns (channels and lags)
        col_sel = ones(1,size(RXX,2));

        % indicator vector for channels
        node_ids = ones(1,no_of_channels);
        node_ids = logical(node_ids);

        % Vector to store removed channels
        deleted_channels = zeros(no_of_channels,1);
        del_ch_count = 0;
        indep_vars = 0;
        recursive_comp = 0;
        csflag = 1;

        while(csflag)

            % Populate a list of remaining channel numbers
            temp_chnl_list = chnl_list(node_ids,:);
            X_sel = RXX(logical(col_sel), logical(col_sel));
            RXY_sel = RXY(logical(col_sel),:);              

            % Check if normal utility or minimum norm utility definition is
            % to be used based on the covariance being singular
            if(rcond(X_sel)<1e-12)
                min_norm_flag = 1;
            else                
                min_norm_flag = 0;
                % Keeps track of independent columns/variables in the data
                indep_vars = indep_vars+1;
            end
            
            util = zeros(size(temp_chnl_list,1),1);
            
            if(min_norm_flag)                         
                 % Recursive computation of inverse - reduce complexity
                 % Only after the first inverse has been computed
                if(recursive_comp)
                    k = idx;                                       
                    S = Xinv((k-1)*noflags+1:k*noflags,(k-1)*noflags+1:k*noflags);
                    
                    V1 =  Xinv((k-1)*noflags+1:k*noflags,1:(k-1)*noflags);
                    V2 =  Xinv((k-1)*noflags+1:k*noflags,(k*noflags)+1:end);
                    V = [V1,V2];
                    
                    C = [];
                    C1 = Xinv(1:(k-1)*noflags,1:(k-1)*noflags);
                    C = [C1,Xinv(1:(k-1)*noflags,(k*noflags)+1:end)];
                    C2 = Xinv((k*noflags)+1:end,1:(k-1)*noflags);
                    C = [C;C2,Xinv((k*noflags)+1:end,(k*noflags)+1:end)];
                    Xinvnew = C- V'*(S\V);
                    
                    Xinv = Xinvnew;
                else
                    % Find regularized inverse for minimum norm utility definition
                    % Select the minimum positive eigenvalue as lambda
                    eigvals = eig(X_sel);
                    lambda_scaling = min(eigvals(eigvals>1e-4));
                    lambda_I = (lambda_scaling*1.0e-5)*eye(size(X_sel,1)); 
                    Xinv = (X_sel + min_norm_flag*lambda_I)\eye(size(X_sel,1));
                    recursive_comp = 1;
                end
            else
                Xinv = X_sel\eye(size(X_sel,1));
            end

            % Compute the new decoder with remaining channels
            W = Xinv * RXY_sel;
            
            % Create group-IDs if lagged versions in data
            if(lagsflag)
                nofgrps = size(temp_chnl_list,1); 
                grpid = repmat(1:nofgrps, noflags, 1);
                grpid = reshape(grpid(:,1:end),[],1);
            end
            % Compute utility of all remaining channels
%             util = grputilcalc(Xinv, W, grpid);
            
            for k = 1:size(temp_chnl_list,1)
                % Select decoder weights of channel k and its lags
                Wkq = W((k-1)*noflags+1:k*noflags);
                % utility computation of channel k ( if block utility: channel and its lags)
                S = Xinv((k-1)*noflags+1:k*noflags,(k-1)*noflags+1:k*noflags);
                util(k) = Wkq'*(S\Wkq);
            end

            % Find the index of channel with least utility
            [~, idx] = min(util);

            % Pick the actual channel number of least utility
            % from list of remaining channel numbers
            temp_ch_sel = temp_chnl_list(idx, :);

            % Delete that channel from the set
            row_id = find(chnl_list==temp_ch_sel);
            col_sel((row_id-1)*noflags+1:row_id*noflags) = 0;       
            node_ids(row_id) = false;

            % store the deleted/removed channel
            del_ch_count = del_ch_count+1;            
            deleted_channels(del_ch_count,1) = temp_ch_sel;                
            
            if(no_of_channels-del_ch_count==N || del_ch_count==N)
                csflag = 0;
            end
        end

        % Best N channels/group-IDs in ascending order of significance
        % i.e. the last deleted channel comes first        
        % Storing deleted channels to give an order of significance if
        % selecting all channels
        
        if(no_of_channels==N)
            ch_selected = flipud(deleted_channels);
            ch_selected = ch_selected(1:N);        
        % Best N channels
        else
            ch_selected = chnl_list(node_ids,:);
        end
    else
        error('Supported Methods: utility');
    end
        
end
