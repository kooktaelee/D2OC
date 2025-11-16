% Sample point generator
%   When the distributions are multimodal distributions consisting of gaussian distributions
clear all;
close all;
fulldir_mfile = mfilename("fullpath");
mfile_name = mfilename();
folder_dir = fulldir_mfile(1:end-length(mfile_name));
DF = load(strcat(folder_dir,'DF.mat')).DF;
DF_id = 'DF42';
Function = 0; % 0: create DF_id, 1: delete DF_id, 2: show DF_id
% Parameters when creating new density function
name = DF_id;
Map_info.x_min = 0; % minimum value of x range in domain
Map_info.y_min = 0; % minimum value of y range in domain
Map_info.x_max = 50; % maximum value of x range in domain
Map_info.y_max = 50; % maximum value of y range in domain
Map_info.no_of_division = 1001; % no of division in x and y range

% Mixture of the gaussian distributions
% no_gauss: the number of gaussian (0 indicates no gaussian distribution is used)
% Mu: (cell array) Mean of each gaussian distribution. {[mean_x1, mean_y1],[mean_x2, mean_y2] ...}
% Sig: (cell array) Covariance of each gaussian distribution. {[Var_xx1, Var_xy1;Var_yx1, Var_yy1],[Var_xx2, Var_xy2;Var_yx2, Var_yy2],...}
% N_sp_gauss: (array) The number of the sample of each gaussian distribution. [# of sp for gauss1, # of sp for gauss2, ...
no_gauss = 20; 
offset_1 = [0,0];
offset_2 = [0,0];
Mu_gauss = {[42.9, 32.82], [41.34, 31.84], [37.57, 29.46], [34.71, 27.58], [31.67, 26.02], [28.56, 24.71], [25.36, 23.48], [21.51, 24.63], [17.42, 26.76], [14.39, 26.51], ...
            [9.5, 27.35], [6.1, 29.15], [ 8.24, 21.1], [33.72, 20.94], [36.1, 21.02], [38.39, 21.19], [33.15, 8.9], [35.52, 9.31] ,[38.72, 10.21], [41.18, 9.47]}; 
Sig_const = 2.5*eye(2);
Sig_gauss = {Sig_const Sig_const Sig_const Sig_const Sig_const Sig_const Sig_const Sig_const Sig_const Sig_const ...
    Sig_const Sig_const Sig_const Sig_const Sig_const Sig_const Sig_const Sig_const Sig_const Sig_const ...
    Sig_const Sig_const Sig_const Sig_const};
Tot_no_sp = 200;
N_sp_gauss = [10 10 10 20 10 20 40 30 30 10 ...
    30 20 10 10 10 10 10 10 8 6]; 
for gauss_cnt = 1:no_gauss
    Sig_gauss{gauss_cnt} = 0.1*N_sp_gauss(gauss_cnt)*eye(2);
end

% Mixture of the lectangular uniform distributions (0 indicates no lect uniform distribution is used)
% no_uniform_rect: the number of lectangular uniform dist'n
% leftdown_rect: (cell array) left-down vertex of the each lectangle {[x1,y1],[x2,y2],...}
% rightup_rect: (cell array) right-up vertex of the each lectangle {[x1,y1],[x2,y2],...}
% N_sp_unif_rect: (array) The number of the sample in each lectangle. [# of sp for lectangle1, # of sp for lectangle2, ...
no_uniform_rect = 0;
leftdown_rect = {[30 30],[60 20]};
rightup_rect = {[40 60],[70,40]};
N_sp_unif_rect = [500,500];

% Mixture of the circular uniform distributions (0 indicates no circular uniform distribution is used)
% no_uniform_cir: the number of circular uniform dist'n
% center_cir: (cell array) the center of the each circle {[x1,y1],[x2,y2],...}
% radius_cir: (array) radius of the each lectangle [radius1, radius2, ...]
% N_sp_unif_cir: (array) The number of the sample of each circle. [# of sp for circle1, # of sp for circle2, ...
no_uniform_cir = 0;
center_cir = {[40 50],[70 60]};
radius_cir = [20,10];
N_sp_unif_cir = [1000,1000];

sp_act_ratio = 0.3/1; % The ratio between sample distn and actual distn

%% create DF_item (when Function == 0)
if Function == 0
    sp = [];
    ActDistn = [];
    % generating sample-points and the actual dist'n for gaussian distributions
    for cnt = 1:no_gauss
        sp_gauss{cnt} = mvnrnd(Mu_gauss{cnt},Sig_gauss{cnt},N_sp_gauss(cnt)); % sample points by Gaussian distribution
        ActDistn_gauss{cnt} = mvnrnd(Mu_gauss{cnt},Sig_gauss{cnt},ceil(sp_act_ratio*N_sp_gauss(cnt)));
        sp = [sp;sp_gauss{cnt}]; 
        ActDistn = [ActDistn;ActDistn_gauss{cnt}];
    end
    
    % generating sample-points and the actual dist'n for rectancular distributions
    for cnt = 1:no_uniform_rect
        Wid_Leng{cnt} = rightup_rect{cnt}-leftdown_rect{cnt};
        if (Wid_Leng{cnt}(1) <= 0)|(Wid_Leng{cnt}(2) <= 0)
            error(strcat("Check the rightup_rect and leftdown_rect of ",num2str(cnt), "th index"));
        end
        sp_rect{cnt} = leftdown_rect{cnt}+...
            [rand(N_sp_unif_rect(cnt),1)*Wid_Leng{cnt}(1) rand(N_sp_unif_rect(cnt),1)*Wid_Leng{cnt}(2)];
        ActDistn_rect{cnt} = leftdown_rect{cnt}+...
            [rand(ceil(sp_act_ratio*N_sp_unif_rect(cnt)),1)*Wid_Leng{cnt}(1) ...
            rand(ceil(sp_act_ratio*N_sp_unif_rect(cnt)),1)*Wid_Leng{cnt}(2)];
        sp = [sp;sp_rect{cnt}]; 
        ActDistn = [ActDistn;ActDistn_rect{cnt}];
    end

    % generating sample-points and the actual dist'n for circular distributions
    for cnt = 1:no_uniform_cir
        random_radius = rand(N_sp_unif_cir(cnt),1);
        random_theta = 2*pi*rand(N_sp_unif_cir(cnt),1);
        sp_cir{cnt} = center_cir{cnt}+radius_cir(cnt)*...
            [sqrt(random_radius).*cos(random_theta) ...
            sqrt(random_radius).*sin(random_theta)];
        random_radius = rand(ceil(sp_act_ratio*N_sp_unif_cir(cnt)),1);
        random_theta = 2*pi*rand(ceil(sp_act_ratio*N_sp_unif_cir(cnt)),1);
        ActDistn_cir{cnt} = center_cir{cnt}+radius_cir(cnt)*...
            [sqrt(random_radius).*cos(random_theta) ...
            sqrt(random_radius).*sin(random_theta)];
        sp = [sp;sp_cir{cnt}]; 
        ActDistn = [ActDistn;ActDistn_cir{cnt}];
    end
    
    % Delete the sample-points outside of the domain
    InDomain_IDX = find(Map_info.x_min<sp(:,1) & Map_info.x_max>sp(:,1) ...
        & Map_info.y_min<sp(:,2) & Map_info.y_max>sp(:,2));
    
    sp = sp(InDomain_IDX,:);

    % generating PDF (gaussian distribution)
    PDF = zeros(Map_info.no_of_division-1);
    x_MS = (Map_info.x_max-Map_info.x_min)/(Map_info.no_of_division-1); %size of the mesh side in x direction
    y_MS = (Map_info.y_max-Map_info.y_min)/(Map_info.no_of_division-1); %size of the mesh side in y direction
    if (no_uniform_rect == 0) & (no_uniform_cir == 0)
        for cnt = 1:no_gauss
            [PDF_gauss{cnt} array1_x1 array1_x2] = gauss2d(linspace(Map_info.x_min+x_MS/2,Map_info.x_max-x_MS/2,Map_info.no_of_division-1),linspace(Map_info.y_min+y_MS/2,Map_info.y_max-y_MS/2,Map_info.no_of_division-1),Mu_gauss{cnt},Sig_gauss{cnt}); % PDF in discretized grid space
            PDF = PDF + N_sp_gauss(cnt)*PDF_gauss{cnt}; % total density distribution
        end
    end
    max_PDF = max(PDF,[],'all'); % find maximum value to normalize the PDF
    norm_PDF = PDF/max_PDF*100; % Density distribution with maximum density of 1.

    % Parameters when creating new density function
    DF_data.name = name;
    DF_data.Map_info=Map_info;
    
    DF_data.no_gauss = no_gauss; 
    DF_data.Mu_gauss = Mu_gauss; 
    DF_data.Sig_gauss = Sig_gauss; 
    DF_data.N_sp_gauss = N_sp_gauss; 
    
    DF_data.no_uniform_rect = no_uniform_rect;
    DF_data.leftdown_rect = leftdown_rect;
    DF_data.rightup_rect = rightup_rect;
    DF_data.N_sp_unif_rect = N_sp_unif_rect;
    
    DF_data.no_uniform_cir = no_uniform_cir;
    DF_data.center_cir = center_cir;
    DF_data.radius_cir = radius_cir;
    DF_data.N_sp_unif_cir = N_sp_unif_cir;
    
    DF_data.sp_act_ratio = sp_act_ratio; % The ratio between sample distn and actual distn
    DF_data.sp = sp;
    DF_data.ActDistn = ActDistn;
    DF_data.norm_PDF = norm_PDF;
    DF_data.PDF_XRange = array1_x1;
    DF_data.PDF_YRange = array1_x2;

elseif Function == 1
    if ~isfield(DF,DF_id)
        error(strcat("There is no field named ",DF_id));
    end
    f=input('\nDo you really want to remove a density function? (y/n)','s');
    if ~strcmp(lower(f),'y')
        error("User choose not to remove a density function");
    end
    DF = rmfield(DF,DF_id);
    save(strcat(folder_dir,'DF.mat'),"DF",'-append');
    return;
else
    DF_data = getfield(DF,DF_id);
    Map_info = DF_data.Map_info;
    sp = DF_data.sp;
    ActDistn = DF_data.ActDistn; % no. of sample points for each Gaussian distribution
end
%% visualization
fig1 = figure(1);
set(gcf,'color', 'none');
set(gca, 'color', 'none');
% clim([-max(new_weight{1})/5 max(new_weight{1})]);
scatter(sp(:,1),sp(:,2),16,[0.8500 0.3250 0.0980],MarkerFaceColor='flat',MarkerEdgeColor=[0.8500 0.3250 0.0980]);
alpha(.2);
axis("equal");
hold on;
% contour(array1_x1, array1_x2,normalized_PDF);
xlim([Map_info.x_min Map_info.x_max]);
ylim([Map_info.y_min Map_info.y_max]);
% xlabel("x [m]");
% ylabel("y [m]");
grid on;
% fig1.Children(1).XTickLabel = [linspace(Map_info.x_min,Map_info.x_max,5)];
% fig1.Children(1).YTickLabel = [linspace(Map_info.y_min,Map_info.y_max,5)];

fig1.Children(1).XTick = [linspace(Map_info.x_min,Map_info.x_max,5)];
fig1.Children(1).YTick = [linspace(Map_info.y_min,Map_info.y_max,5)];
fig1.Children(1).XTickLabel = ["" "" "" "" ""];
fig1.Children(1).YTickLabel = ["" "" "" "" ""];
fig1.Position = [30 30 600 600];

fig2 = figure(2);
scatter(ActDistn(:,1),ActDistn(:,2),16,[0.8500 0.3250 0.0980],MarkerFaceColor='flat',MarkerEdgeColor=[0.8500 0.3250 0.0980]);
alpha(.2);
axis("equal");
hold on;
% contour(array1_x1, array1_x2,normalized_PDF);
xlim([Map_info.x_min Map_info.x_max]);
ylim([Map_info.y_min Map_info.y_max]);
% xlabel("x [m]");
% ylabel("y [m]");
grid on;
fig2.Children(1).XTick = [linspace(Map_info.x_min,Map_info.x_max,5)];
fig2.Children(1).YTick = [linspace(Map_info.y_min,Map_info.y_max,5)];
fig2.Children(1).XTickLabel = ["" "" "" "" ""];
fig2.Children(1).YTickLabel = ["" "" "" "" ""];
fig2.Position = [30 30 600 600];


if isfield(DF_data,"norm_PDF")
fig3 = figure(3);
s2=surf(DF_data.PDF_XRange,DF_data.PDF_YRange,DF_data.norm_PDF);
s2.EdgeColor = "none";
cm2 = colormap("jet");
axis("equal");
hold on;
xlim([Map_info.x_min Map_info.x_max]);
ylim([Map_info.y_min Map_info.y_max]);
grid on;
% s2.Parent.XAxis.Label.String = "x [m]";
% s2.Parent.XAxis.FontSize = 18;
% s2.Parent.XAxis.Label.FontSize = 18;
% s2.Parent.YAxis.Label.String = "y [m]";
% s2.Parent.YAxis.FontSize = 18;
% s2.Parent.YAxis.Label.FontSize = 18;
view(0,90);
fig3.Children(1).XTick = [linspace(Map_info.x_min,Map_info.x_max,5)];
fig3.Children(1).YTick = [linspace(Map_info.y_min,Map_info.y_max,5)];
fig3.Children(1).XTickLabel = ["" "" "" "" ""];
fig3.Children(1).YTickLabel = ["" "" "" "" ""];
fig3.Position = [30 30 600 600];
end

if isfield(DF_data,"norm_PDF")
fig4 = figure(4);
s3=contour(DF_data.PDF_XRange,DF_data.PDF_YRange,DF_data.norm_PDF,LineWidth=1.5);
cm2 = colormap("jet");
axis("equal");
hold on;
xlim([Map_info.x_min Map_info.x_max]);
ylim([Map_info.y_min Map_info.y_max]);
% grid on;
% s3.Children(1).XAxis.Label.String = "x [m]";
% s3.Children(1).XAxis.FontSize = 18;
% s3.Children(1).XAxis.Label.FontSize = 18;
% s3.Children(1).YAxis.Label.String = "y [m]";
% s3.Children(1).YAxis.FontSize = 18;
% s3.Children(1).YAxis.Label.FontSize = 18;
fig4.Children(1).XTick = [linspace(Map_info.x_min,Map_info.x_max,5)];
fig4.Children(1).YTick = [linspace(Map_info.y_min,Map_info.y_max,5)];
fig4.Children(1).XTickLabel = ["" "" "" "" ""];
fig4.Children(1).YTickLabel = ["" "" "" "" ""];
fig4.Position = [30 30 600 600];
end

if Function == 0
    if isfield(DF,DF_id)
        f=input('The same field name exists.\nDo you really want to overwrite an existing density function? (y/n)','s');
        if ~strcmp(lower(f),'y')
            error("Choose another name for the density function");
        end
    end
    DF = setfield(DF,DF_id,DF_data);
    save(strcat(folder_dir,'DF.mat'),"DF",'-append');
end