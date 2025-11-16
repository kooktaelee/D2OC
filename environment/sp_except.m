DF_data = DF.DF16;
y_raw = DF_data.sp;
obs_offset = DF_data.obs.offset;
rec_obs = DF_data.obs.Rectan;
circ_obs = DF_data.obs.Circle;

%% rectacgular obstacle generation
div_ob = 1;
barrier_thickness = 3;
rec_obs_total=[];

for cnt = 1:length(rec_obs)
    x_obs_l = rec_obs{cnt}(1);
    y_obs_t = rec_obs{cnt}(2);

    x_obs_r = x_obs_l + rec_obs{cnt}(3);
    y_obs_b = y_obs_t + rec_obs{cnt}(4);
    
    x_npts = floor(rec_obs{cnt}(3)/div_ob);
    y_npts = floor(rec_obs{cnt}(4)/div_ob);

    x_pt = linspace(x_obs_l,x_obs_r,x_npts)';
    y_pt = linspace(y_obs_t,y_obs_b,y_npts)';
    
    CXY = [];
    for cnt_thickness = 1:barrier_thickness
        CXY = [CXY;
            x_pt(cnt_thickness)*ones(y_npts,1), y_pt;
            x_pt(end-cnt_thickness+1)*ones(y_npts,1), y_pt;
            x_pt(barrier_thickness+1:end-barrier_thickness), y_pt(cnt_thickness)*ones(x_npts-2*barrier_thickness,1);
            x_pt(barrier_thickness+1:end-barrier_thickness), y_pt(end-cnt_thickness+1)*ones(x_npts-2*barrier_thickness,1)];
    end
    
        rec_obs_total = [rec_obs_total ; CXY];
end
%% circular obstacle generation
circ_obs_total=[];
for cnt1 = 1:length(circ_obs)
    num_div_r = floor(circ_obs{cnt1}(3)/div_ob);
    div_r = linspace(0,circ_obs{cnt1}(3),num_div_r);
    circ_obs_center = [circ_obs{cnt1}(1), circ_obs{cnt1}(2)];
    circ_obs_total = [circ_obs_total;circ_obs_center];

    for cnt2 = length(div_r)-barrier_thickness+1:length(div_r)
            circumf = div_r(cnt2)*2*pi; %circumference
            num_div_th = floor(circumf/div_ob);
            div_th = linspace(0,2*pi,num_div_th);
            x_circ = div_r(cnt2)*cos(div_th(1:end-1))'+circ_obs_center(1,1);
            y_circ = div_r(cnt2)*sin(div_th(1:end-1))'+circ_obs_center(1,2);
            circ_obs_total = [circ_obs_total;[x_circ, y_circ]];
    end
end

%%
op = [rec_obs_total;circ_obs_total];

for cnt1=1:length(rec_obs)
    if exist("y_index")
    y_index=[y_index;find((y_raw(:,1)>rec_obs{cnt1}(1)-obs_offset)&(y_raw(:,1)<rec_obs{cnt1}(1)+rec_obs{cnt1}(3)+obs_offset)&(y_raw(:,2)>rec_obs{cnt1}(2)-obs_offset)&(y_raw(:,2)<rec_obs{cnt1}(2)+rec_obs{cnt1}(4)+obs_offset))];
    else
        y_index=find((y_raw(:,1)>rec_obs{cnt1}(1)-obs_offset)&(y_raw(:,1)<rec_obs{cnt1}(1)+rec_obs{cnt1}(3)+obs_offset)&(y_raw(:,2)>rec_obs{cnt1}(2)-obs_offset)&(y_raw(:,2)<rec_obs{cnt1}(2)+rec_obs{cnt1}(4)+obs_offset));
    end
end

for cnt1=1:length(circ_obs)
    if exist("y_index")
    y_index=[y_index;find(vecnorm(y_raw-[circ_obs{cnt1}(1) circ_obs{cnt1}(2)],2,2)<circ_obs{cnt1}(3)+obs_offset)];
    else
        y_index=find(vecnorm(y_raw-[circ_obs{cnt1}(1) circ_obs{cnt1}(2)],2,2)<circ_obs{cnt1}(3)+obs_offset);
    end
end
y = y_raw(setdiff(1:size(y_raw,1), y_index), :);

subplot(1,2,1);
scatter(DF_data.sp(:,1),DF_data.sp(:,2),2*ones(length(y_raw),1));
hold on;
scatter(DF_data.obs.op(:,1),DF_data.obs.op(:,2),2,'r');

subplot(1,2,2);
scatter(y(:,1),y(:,2),2*ones(length(y),1));
hold on;
scatter(DF_data.obs.op(:,1),DF_data.obs.op(:,2),2,'r');

figure(3)
scatter(op(:,1),op(:,2),2,'blue');