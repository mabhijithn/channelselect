import numpy as np
from itertools import compress
def channel_select(**argV):
''' 
A function which selects the best columns of A in the least-squares problem 
  min_x 1/2 ||Ax-b||_2^2 

  Usage: ch_selected = channel_select(A=A, b=b, method='utility')
  
  Input: A - matrix with channels/features along columns and time/length of features along rows
            size: (T X N) - N channels; T duration/length
        b - Vector of size (T X 1)
         method - a string; currently supported : 'utility' [1,2]
 
 [1] Narayanan, A. M., & Bertrand, A. (2019). Analysis of miniaturization effects and channel selection strategies for EEG sensor networks with application to auditory attention detection. IEEE Transactions on Biomedical Engineering.
 [2] Bertrand, A. (2018). Utility Metrics for Assessment and Subset Selection of Input Variables for Linear Estimation [Tips & Tricks]. IEEE Signal Processing Magazine, 35(6), 93–99.


'''
   if 'A' in argV:
        A = argV['A']
    if 'b' in argV:
        b = argV['b']
    if 'N' in argV:
        N = argV['N']        
    else:
        print('Number of channels to be selected required!')
    if 'method' in argV:
        method = argV['method']
    else:
        # By default channel selection method is utility
        method = 'utility'
    if 'lags' in argV:
        lags = argV['lags'] + 1
    else:
        lags = 1
    
    if method is 'utility':
        # Get the dimensions of the data
        m = A.shape[0]
        n = A.shape[1]
        
        
        no_of_channels = int(n/lags)
        ch_list = list(range(1,no_of_channels+1))
        
        # Compute covariances
        At = np.matrix.transpose(A)
        AtA = At @ A
        Atb = At @ b
        
        # Initialise logical indexing masks to select/remove channels/columns 
        utility_vals = np.zeros((no_of_channels,1), dtype=float)
        del_ch_count = 0
        col_sel_mask = np.ones((n,), dtype=bool)
        ch_sel_mask = np.ones((no_of_channels,), dtype=bool)
        ch_sel = list(compress(ch_list,ch_sel_mask))
        ch_selected = []

        # Start iterative channel deletion based on utility
        while del_ch_count<no_of_channels:
            
            AtA_sel = AtA[:,np.array(col_sel_mask)]
            AtA_sel = AtA_sel[np.array(col_sel_mask),:] 
            
            
            
            # Estimate the current decoder using (A'A)^(-1)(A'b)
            
            # This decoder is used for the fast computation of utlity of a column or a group of columns
            AtA_rank = np.linalg.matrix_rank(AtA_sel)
            if AtA_rank==AtA_sel.shape[1]:
                # If AtA is full-rank, normal definition of utility
                AtA_sel_inv = np.linalg.inv(AtA_sel)
            else:                
                # l2-regularization based definition of utility when AtA is not full-rank
                id_mat = np.identity(AtA_sel.shape[1])
                diag_mat = np.diag(AtA_sel)
                nz_diag = diag_mat>0
                lmb_val = np.amin(diag_mat[nz_diag])
                AtA_sel = AtA_sel + (lmb_val * 1e-5 * id_mat)
                AtA_sel_inv = np.linalg.inv(AtA_sel)
            Atb_sel = Atb[col_sel_mask]
            W_sel = AtA_sel_inv @ Atb_sel            
            indx = np.array(np.where(ch_sel_mask==True))
            util = np.zeros((indx.shape[1],1),dtype='float')
            k = 1      
        # Compute (group)-utility of columns of current A
            while k<=len(util):
                if lags==1:
                    S = AtA_sel_inv[k-1][k-1]
                    util[k-1] = np.matrix.transpose(W_sel[k-1]) @ ((1/S) * W_sel[k-1])
                else:
                    strt = (k-1)*lags
                    stp = strt+lags
                    S = AtA_sel_inv[strt:stp,:]
                    S = S[:,strt:stp]
                    W_block = W_sel[strt:stp]
                    util[k-1] = np.matrix.transpose(W_block) @ np.linalg.inv(S) @ W_block
                k = k+1
                
        # Eliminate the column or group of columns corresponding to lowest utility
            min_util = np.argmin(util)
            ch_sel_mask[ch_sel[min_util]-1] = False
            strt = (ch_sel[min_util]-1)*lags
            stp = strt+lags
            col_sel_mask[strt:stp] = False
       # Append the deleted column number to a list
            ch_selected.append(ch_sel[min_util])
            ch_sel = list(compress(ch_list,ch_sel_mask))
            del_ch_count = del_ch_count+1
      # The order of channels in ch_selected list in the reverse of their importance
      # That is, the first channel in the list is the least important (one that was deleted first)
        return [ch_selected]
            
    return
