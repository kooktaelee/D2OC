% Generate random gaussian distribution parameters
DF_data = DF.DF16;
Map_info = DF_data.Map_info;
Map_size = [Map_info.x_max-Map_info.x_min,Map_info.y_max-Map_info.y_min];
Margin = 0.05;
Map_margin = Map_size*Margin;
no_gauss = DF_data.no_gauss;
no_sample = 20;

for cnt_gauss = 1:no_gauss
    Mu{cnt_gauss} = ...
        [Map_info.x_min+Map_margin(1) + Map_size(1)*(1-2*Margin)*rand(1),...
        Map_info.y_min+Map_margin(2) +  Map_size(2)*(1-2*Margin)*rand(1)];
    Sig{cnt_gauss} = ...
        [(Map_size(1)*0.008*rand(1)+Map_size(1)*0.004)^2 0;...
        0 (Map_size(2)*0.008*rand(1)+Map_size(2)*0.004)^2];
    N_sample(cnt_gauss) = no_sample;
end

DF_data.Mu = Mu;
DF_data.Sig = Sig;
DF_data.N_sample = N_sample;
