function [xg, lamb, y_index] = hamilton_optimal_control_Lin_Quad_R3(x,alpha_k,z,new_weight,PosWgt_idx,dist_sq)
%% This is a function to select the local sample-points and their centroid
% It selects the local sample-points until their total weight matches the weight of the agent-point.
% input parameters are:
% x: the location of the current agent-point
% alpha_k : the weight of the agent-point that should be transported to local sample points
% z: the coordinate set of SPs
% new_weight: weight for each point in z
n = size(z,1); % n is the number of point in z
rem_weight = alpha_k;

d_wnE = dist_sq./(transpose(new_weight(PosWgt_idx)).^2);


[out,idx] = sort(d_wnE);
y_index_sorted = PosWgt_idx(idx);

cnt_sample = 1;

%% Select the local sample-points 
while (rem_weight > new_weight(y_index_sorted(cnt_sample)))&&(cnt_sample<length(y_index_sorted))
    y_loc_pos(cnt_sample,:) = z(y_index_sorted(cnt_sample),:);
    y_loc_wgt(cnt_sample) = new_weight(y_index_sorted(cnt_sample));
    rem_weight = rem_weight-y_loc_wgt(cnt_sample);
    cnt_sample = cnt_sample+1;
end
if cnt_sample<=length(y_index_sorted)
y_loc_pos(cnt_sample,:) = z(y_index_sorted(cnt_sample),:); % y_loc_pos: the position of the local sample-points
y_loc_wgt(cnt_sample) = rem_weight; % y_loc_wgt: the weight of the local sample-points   
end

%% Calculate the centroid of the local sample-points
lamb = zeros(1,2);
for i=1:size(y_loc_pos,1) 
    lamb = lamb + (y_loc_wgt(i) * [x(1) - y_loc_pos(i,1), x(2) - y_loc_pos(i,2)]);% / norm([x(1) - y_pos(i,1), x(2) - y_pos(i,2)]) );
end
xg = transpose(-lamb/alpha_k+x);
y_index = y_loc_pos;

