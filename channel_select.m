% A function which selects the best N channels of A for the least-squared 
% problem min_{x} 1/2 * ||Ax - b ||² which best estimates b.
%
% Required Inputs:
%   A - Data/measurements used in estimation where P is the number of 
%       channels (including delayed versions) and T is the number of
%       time samples. 
%       If A consists of delayed samples of a channel, all delays of 
%       a channel are expected to be consecutive columns
%       Size: (T X P)
%   b - The desired data to be estimated using A. Size: (T X 1)
%   N - Number of channels to be selected.
% Optional Inputs:
%   'method' - A string for the method to be used to select best N channels
%            Currently supported:
%            'utility' - Utility based channel reduction. Here, if A is 
%            singular, utility definition with minimum norm increase 
%            interpretation is used.
%            (ref. A. Bertrand, 2018)
%            'lasso' - LASSO based channel selection
%
%   'lags' - Number of delays (in number of samples) applied to channels in A
%          default: 0 (assumes no delays)
%   
%   'galiso' - No galvanic isolation (default) : An empty matrix if no galvanic (WESN context data)
%              Galvanic isolation : A matrix which shows galvanic connection
%              between channels (ref: find_gal_conn.m to generate this
%              matrx)
%               
%   
%
% Outputs: 
%   ch_selected - The indices of the best N channels in A. They are 
%                 ordered in descending order of significance.
function [ch_selected, x_mu] = channel_select(A, b, N, varargin)
    
    % Default conditions
    method = 'utility';
    noflags = 1;
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
                case 'lags'
                    noflags = Value+1;
                case 'galiso'
                    Galconn_mat = Value;
                case 'grpnorm'
                    group_norm = Value;
                case 'covar'
                    RXX = Value;
                case 'crossvar'
                    RXY = Value;
            end
        end
    end
    
    if(strcmp(method,'utility'))
        % Calculate the auto and crosscovariances
        if(~exist('RXX'))
            RXX = (A'*A)/size(A,1);
            RXY = (A'*b)/size(A,1);
        end
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
            lambda_scaling = min(eigvals(eigvals>0));

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
    elseif(strcmp(method,'lasso'))
        % Solving LASSO based node selection using YALMIP
%         load('matfiles/nearest_neighbours_50mm.mat','chnl_list');
        n = size(A,2)/noflags;
        m = size(A,1);
        lassochnum = N;
        
        Gmat = zeros(n,size(A,2));
        strt = 1;
        stp = strt+noflags-1;
        for i = 1:n            
            Gmat(i,strt:stp) = 1;
            strt = strt+noflags;
            stp = strt+noflags-1;
        end
        for i = 1:size(Galconn_mat, 1)
            for j = 1:i-1
               if(Galconn_mat(i,j)==1)
                    Galconn_mat(i,j) = 0;
               end               
            end      
        end
        no_of_constraints = length(find(Galconn_mat==1));
        [r,c] = find(Galconn_mat==1);
        Galconn_mat_new = zeros(no_of_constraints,n);
        for i = 1:size(r,1)
            Galconn_mat_new(i,r(i)) = 1;
            Galconn_mat_new(i,c(i)) = 1;
        end
        % Following code adapted from https://osqp.org/docs/examples/lasso.html
        x = sdpvar(size(A,2), 1);
        Rxx = (A'*A);
        Rxy = (A'*b);
        options = sdpsettings('solver', 'gurobi','verbose', 2);
        lambda = lassochnum;
        M = 100;               
        z = binvar(n, 1);
        
        % Group lasso objective        
        k = sdpvar;
        % Contraints for group-LASSO (does not work with lags yet)
%         Constraints = [sum(z)<=k,-M*z<=Gmat*x<=M*z];
        
        % Constraints for galvanic isolation
%         Constraints = [sum(z)<=k, z'*Galconn_mat*z<=1,-M*z<=Gmat*x<=M*z];
%         Constraints = [sum(z)<=k, z'*Galconn_mat*z<=1,-M*z<=Gmat*x<=M*z];
        grpM = zeros(size(A,2), n);
        col_mat = [M*ones(noflags,1);zeros(size(A,2)-noflags,1)];
        for i = 1:n
            grpM(:,i) = circshift(col_mat, noflags*(i-1));
        end
        
        %         objective = 0.5*norm(A*x - b)^2 + norm(x,1);
%         objective_quad = (0.5*x'*Rxx*x - x'*Rxy) + norm(x,1);
        objective_quad = (0.5*x'*Rxx*x - x'*Rxy);
%         Constraints = [sum(z)<=k, z'*Galconn_mat*z<=1, -grpM*z<=x<=grpM*z];
%         Constraints = [sum(z)<=k, Galconn_mat_new*z<=1,-M*z<=Gmat*x<=M*z];
        Constraints = [sum(z)==k,Galconn_mat_new*z<=1,Gmat*abs(x)<=M*z];
%         Constraints = [sum(z)==k,-grpM*z<=x<=grpM*z];
        x_opt = optimizer(Constraints, objective_quad, options, k, {x,z});
        x_mu = x_opt(lambda);
        tmp = Gmat*x_mu{1,1};
        ch_selected = find(abs(tmp)>1e-3);
%         tmp = x_mu{1,1};        
%         ch_selected{1} = {chnl_list(abs(tmp)>1e-3,:)};
%         ch_selected{2} = x_mu;
    else
        error('Supported Methods: utility, lasso');
    end
        
end
