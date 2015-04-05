function [I_rs, I_rel, I_max, H_tmp] = MI(H)
%CALCULATES MUTUAL INFORMATION FROM A CONFUSION MATRIX
%  [H_rs, H_rel, H_max] = MI(H)
%
%ARGS:
%  H    - confusion matrix
%
%RETURNS
%  I_rs  - mutual information
%  I_rel - relative mutual information
%  I_max - theoretical value of the max info

%created 07/10/16 Jan


p_rs = H./sum(H(:));%normalize
p_r = sum(p_rs,1);p_r = p_r./sum(p_r);%marginals
p_s = sum(p_rs,2);p_s = p_s./sum(p_s);%marginals

p_rxs = (p_r(:)*p_s(:)')';%independent joint pdf
H_tmp = p_rs.*log2(p_rs./p_rxs);%term to sum over

I_rs = nansum(H_tmp(:));%sum missing nan's
I_max = log2(size(H,1));%maximally achievable Information
I_rel = I_rs/I_max;%relative information

