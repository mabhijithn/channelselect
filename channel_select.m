% A function which selects the best N channels of A for the least-squared 
% problem min_{x} 1/2 * ||Ax - b ||² which best estimates b.
%
% Inputs:
%   A - Data/measurements used in estimation where P is the number of 
%       channels (including delayed versions) and T is the number of
%       time samples. 
%       If A consists of delayed samples of a channel, all delays of 
%       a channel are expected to be consecutive columns
%       Size: (T X P)
%   b - The desired data to be estimated using A. Size: (T X 1)
%   N - Number of channels to be selected.
%   method - A string for the method to be used to select best N channels
%            Currently supported:
%            'utility' - Utility based channel reduction. Here, if A is 
%            singular, utility definition with minimum norm increase 
%            interpretation is used.
%            (ref. A. Bertrand, 2018)
%
%   lags - Number of delays (in number of samples) applied to channels in A
%          default: 0 (assumes no delays)
%
% Outputs: 
%   ch_selected - The indices of the best N channels in A. They are 
%                 ordered in descending order of significance.
function ch_selected = channel_select(A, b, N, method, lags)
    
    
    if(nargin<3)
        method = 'utility';
    elseif(nargin<4)
        method = 'utility';
        lags = 0;
    end
    
    if(isempty(method) || ~strcmp(method,'utility'))
        warning('Method is not supported. Proceeding with the only currently supported method: utility');
    end
    
    
    % Calculate the auto and crosscovariances
    RXX = (A'*A)/size(A,1);
    RXY = (A'*b)/size(A,1);
    noflags = lags+1;
    no_of_channels = size(RXX,2)/noflags;

    % Initialise a list of original indices/channels
    chnl_list = (1:no_of_channels)';
    
    % Check if normal utility or minimum norm utility definition
    % to be used based on rank of covariance
    if(rank(RXX)<size(RXX,2))
        min_norm_flag = 1;
    else
        min_norm_flag = 0;
    end
    
    % Indicator matrices to book-keep removed channels
    
    % selector of columns (channels and lags)
    col_sel = ones(1,size(RXX,2));
    
    % indicator vector for channels
    node_ids = ones(1,no_of_channels);
    node_ids = logical(node_ids);
    
    % Vector to store removed channels
    deleted_channels = zeros(no_of_channels,1);
    del_ch_count = 0;
    
    while(del_ch_count<no_of_channels)
        
        % Populate a list of remaining channel numbers
        temp_chnl_list = chnl_list(node_ids,:);
        X_sel = RXX(logical(col_sel), logical(col_sel));
        RXY_sel = RXY(logical(col_sel));              
        
        util = zeros(size(temp_chnl_list,1),1);

        eigvals = diag(X_sel);
        lambda_scaling = min(eigvals(eigavals>0));
        
        lambda_I = (lambda_scaling*1.0e-5)*eye(size(X_sel,1)); 
        Xinv = inv(X_sel + min_norm_flag*lambda_I);
        
        % Compute the new decoder with remaining channels
        W = Xinv * RXY_sel;

        % Compute utility of all remaining channels
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
    end
    
    % Best N channels in ascending order of significance
    % i.e. the last deleted channel comes first
    ch_selected = flipud(deleted_channels);
    ch_selected = ch_selected(1:N);
end
