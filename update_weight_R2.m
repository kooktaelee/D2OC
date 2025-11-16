function [new_weight, weight_distribution, PosWgt_idx, dist_sq] = update_weight_R2(x,z,new_weight,alpha_k,N)
% input arguments
%   x : current location of the agent
%   z : location arrays of sample-points
%   new_weight : weight arrays of sample-points
%   tend : total number of agent-points
%   N : the number of sample-points
% Onput arguments
%   new_weight : updated weight arrays of sample-points
%   weight_distribution : delta new_weight
%   PosWgt_index : index of sample-points with positive weight
%   dist_sq : squared distance between sample-points and the agent
remainder = alpha_k;
weight_distribution = zeros(N,1);

PosWgt_idx = find(new_weight>0); %index of sample-points with positive weight
dist_sq = sum(([z(PosWgt_idx,:)-[x(1,1)*ones(length(PosWgt_idx),1),...
    x(1,2)*ones(length(PosWgt_idx),1)]]').^2); 
[out,idx] = sort(dist_sq);
sp_idx_sorted = PosWgt_idx(idx); %sorted index of sample-points wrt distance

for i=1:length(sp_idx_sorted)
    prev_weight = new_weight;
    new_weight(sp_idx_sorted(i)) = new_weight(sp_idx_sorted(i)) - remainder;
    if new_weight(sp_idx_sorted(i)) >=0 
        weight_distribution(sp_idx_sorted(i)) = remainder;
        break; %after every remainder is distributed, finish the update.
    else
        remainder = -new_weight(sp_idx_sorted(i));
        new_weight(sp_idx_sorted(i)) = 0;
        weight_distribution(sp_idx_sorted(i)) = prev_weight(sp_idx_sorted(i));
    end
end

new_weight(new_weight<0) = 0;
