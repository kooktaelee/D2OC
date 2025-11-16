%% ========================================================================
%  Density-Driven Optimal Control (D2OC) - Single Simulation Runner
%  This version:
%    • Runs a SINGLE test case from Sim_rev60.mat
%    • Assumes existing sim_result / sim_data structure inside the .mat file
%    • Does NOT create videos
%    • Keeps live plotting active for visualization
%    • Core D2OC algorithm is unchanged
% ========================================================================

clc;
clearvars -except g;   % keep global variable g if it exists
close all; clf;

%% ------------------------------------------------------------------------
%  Paths and data loading
%  - Determine script location
%  - Load simulation database Sim_rev60.mat
% ------------------------------------------------------------------------
mfile_name = mfilename();
mfile_dir  = fileparts(mfilename('fullpath'));

% Path to simulation data folder (contains Sim_rev60.mat)
dir = fullfile(mfile_dir, 'sim_data');
up_dir = fileparts(mfile_dir); % one directory above, if needed

if ~exist(dir,'dir')
    mkdir(dir);   % create folder if missing
end

sim_file = 'Sim_rev60';              % .mat file name (without extension)
sim_struct = load(fullfile(dir, sim_file));  % load struct from mat-file

% sim_result = full simulation run history (if exists), otherwise replicate it
sim_data = sim_struct.sim_data;
if isfield(sim_struct, "sim_result")
    sim_result = sim_struct.sim_result;
else
    sim_result = sim_data;           % fallback structure
    sim_result.date = date();
    sim_result.matlab_script = mfile_name;
end

% Test configuration table (cell array)
% Row 1 = header strings
% Row 2+ = individual test configurations
sim_list = sim_result.tests.inputs;

% Dictionary indexing column names → column numbers
sim_dic  = dictionary([sim_list{1,:}], 1:size(sim_list,2));

% Load density-function structures (Gaussian mixture environment)
DF = load(fullfile(mfile_dir,'environment','DF.mat')).DF;

% ------------------------------------------------------------------------
%  Load parameter set (param07) from folder /param
% ------------------------------------------------------------------------
param_file = "param07";
param_list = load(fullfile(mfile_dir,'param'),'param');
param_list = param_list.param;
param_field = fieldnames(param_list);

% Find parameter index matching 'param07'
for param_cnt = 1:length(param_field)
    if param_field{param_cnt} == param_file
        param_no = param_cnt;
        break;
    end
end
param_data = param_list.(param_field{param_no});

%% ------------------------------------------------------------------------
%  Select ONLY ONE test to run
%  You must choose cnt_sim = 2, 3, 4, ... (since row 1 is header)
% ------------------------------------------------------------------------

    clearvars -except cnt_sim mfile_dir dir up_dir sim_file sim_result sim_list sim_dic ...
        mfile_name sim_data DF param_data Share_law_A Share_law_B TestIndexToRun
    close all;

    warning('off')
    tic_start = tic;   % start timing

    %% --------------------------------------------------------------------
    %  Read simulation configuration for this test index
    % ---------------------------------------------------------------------
    cnt_sim = 2;  % ---> choose which test to run (2 = first test)

    % Extract parameters from sim_list using dictionary indexing
    Test_no      = sim_list{cnt_sim,sim_dic('Test_no')};
    bat_life     = sim_list{cnt_sim,sim_dic('bat_life')};
    Hor_Leng     = sim_list{cnt_sim,sim_dic('Hor_Leng')};
    no_of_agents = sim_list{cnt_sim,sim_dic('num_agent')};
    complete     = sim_list{cnt_sim,sim_dic('complete')};
    Data_save_Int = 20;

    % Optional additional fields
    if sim_dic.isKey('comm_thresh')
        comm_thresh = sim_list{cnt_sim,sim_dic('comm_thresh')};
    else
        comm_thresh = 300;  % default communication threshold
    end
    if sim_dic.isKey('agent_init')
        agent_init = sim_list{cnt_sim,sim_dic('agent_init')};
    end

    % Load density-function ID for this test
    if sim_dic.isKey('DF_id')
        DF_id   = sim_list{cnt_sim,sim_dic('DF_id')};
        DF_data = DF.(DF_id);   % access DF structure by name

        % Extract DF fields
        Ref_map        = DF_data.name;
        N_sample       = DF_data.N_sample;
        Mu             = DF_data.Mu;
        Sig            = DF_data.Sig;
        N_raw          = N_sample * numel(Sig);
        y_ref          = DF_data.sp;
        normalized_PDF = DF_data.normalized_PDF;
        Map_info       = DF_data.Map_info;

        if isfield(DF_data,'obs')
            obs = DF_data.obs;  % optional obstacle set
        end
    else
        error('Error. "DF_id" field not defined')
    end

    fprintf('test_no : %d, op_time : %d, Hor_Leng : %d\n', ...
        cnt_sim-1, bat_life, Hor_Leng);

    if complete
        fprintf('already done \n');
    else

        %% ----------------------------------------------------------------
        %  Control parameters & linearized quadrotor model
        %  - Discrete-time model (Ad, Bd)
        %  - Build KKT matrix for horizon Hor_Leng
        %  Linearize Quadrotor State vector x = [x, x_dot, theta, theta_dot, y, y_dot, phi, phi_dot]'
        %  x, y      : position (m)
        %  theta,phi : pitch/roll angles (rad)
        %  x_dot,y_dot,theta_dot,phi_dot : corresponding velocities
        % -----------------------------------------------------------------
        tend = bat_life * no_of_agents;    % maximum total iterations
        comm_time_thresh = 0;
        iter_weight_thresh = 1e-15;        % (legacy condition)
        alpha_k = 1/tend;                  % per-time-step weight fraction

        T = 0.1;                           % discretization step
        umax = param_data.Agent.umax;
        vmax = param_data.Agent.vmax;
        thresh_ang    = param_data.Agent.thresh_ang;
        thresh_angvel = param_data.Agent.thresh_angvel;

        g  = 9.81;
        Ix = 0.0075;
        Iy = Ix;

        % Linearized dynamics for x (pitch) subsystem
        Ad1 = [1 T 0 0;
               0 1 -T*g 0;
               0 0 1 T;
               0 0 0 1];

        % y (roll) subsystem
        Ad2 = [1 T 0 0;
               0 1 T*g 0;
               0 0 1 T;
               0 0 0 1];

        % Input matrices
        Bd1 = [0 0;
               0 0;
               0 0;
               0 T/Iy];
        Bd2 = [0 0;
               0 0;
               0 0;
               T/Ix 0];

        % Combine into 8-state model
        Ad = [Ad1 zeros(4);
              zeros(4) Ad2];
        Ad_tr = Ad.';
        dim_sv = size(Ad,2);

        Bd = [Bd1; Bd2];
        Bd_tr = Bd.';
        dim_iv = size(Bd,2);

        % LQR weighting matrices
        Q_ = param_data.LQR.Q;
        R_ = param_data.LQR.R;

        % Position-extraction matrix for cost x = [pos_x pos_y ...]
        C_ = [1 0 0 0 0 0 0 0;
              0 0 0 0 1 0 0 0];
        C_tr = C_.';

        % ----------------------------------------------------------------
        % Build block-structured KKT matrix (Hessian + dynamics constraints)
        % ----------------------------------------------------------------
        A_AugE{1,1} = kron(eye(Hor_Leng), alpha_k*C_tr*C_ + Q_);
        A_AugE{1,2} = kron([zeros(Hor_Leng-1,1) eye(Hor_Leng-1,Hor_Leng-1); ...
                             0 zeros(1,Hor_Leng-1)], Ad_tr) ...
                    + kron(eye(Hor_Leng), -eye(dim_sv));
        A_AugE{1,3} = zeros(dim_sv*Hor_Leng, dim_iv*Hor_Leng);

        A_AugE{2,1} = kron([zeros(1,Hor_Leng-1) 0; ...
                             eye(Hor_Leng-1,Hor_Leng-1) zeros(Hor_Leng-1,1)], Ad) ...
                    + kron(eye(Hor_Leng), -eye(dim_sv));
        A_AugE{2,2} = zeros(dim_sv*Hor_Leng);
        A_AugE{2,3} = kron(eye(Hor_Leng), Bd);

        A_AugE{3,1} = zeros(dim_iv*Hor_Leng, dim_sv*Hor_Leng);
        A_AugE{3,2} = A_AugE{2,3}.';
        A_AugE{3,3} = kron(eye(Hor_Leng), R_);

        % Full augmented KKT matrix
        A_aug = [A_AugE{1,1} A_AugE{1,2} A_AugE{1,3};
                 A_AugE{2,1} A_AugE{2,2} A_AugE{2,3};
                 A_AugE{3,1} A_AugE{3,2} A_AugE{3,3}];

        A_aug_inv = inv(A_aug);   % inverse once for reuse

        %% ----------------------------------------------------------------
        %  Domain and reference samples
        % -----------------------------------------------------------------
        domain = Map_info.domain;   % rectangular operational area
        y = y_ref;                  % Gaussian mixture sample points
        N = size(y,1);              % number of samples

        rng(1);                     % reproducible randomness

        %% ----------------------------------------------------------------
        %  Agent initialization
        % -----------------------------------------------------------------
        for cnt1 = 1:no_of_agents
            % Random initial position in domain; zero initial attitude & rates
            x{cnt1} = [ ...
                rand*(domain(2)-domain(1))+domain(1), 0, 0, 0, ...
                rand*(domain(4)-domain(3))+domain(3), 0, 0, 0 ]';
        end

        spatial_dist = ones(N,1)/N;  % uniform initial weight distribution

        for n = 1:no_of_agents
            x_traj{n}     = x{n};
            agent_init{n} = x{n};

            % Weight records for decentralized sharing
            for q = 1:no_of_agents
                agent_wgt{n,q} = zeros(N,1);
            end

            comm_time{n} = zeros(no_of_agents,1);
            new_weight{n} = spatial_dist;

            % Precompute squared distance to all sample points
            dist_sq{n} = sum(([y - [x{n}(1)*ones(N,1), x{n}(5)*ones(N,1)]]').^2);

            PosWgt_idx{n} = 1:N;
            weight_sum_check_index(n) = sum(new_weight{n});
        end

        %% ----------------------------------------------------------------
        %  Live plotting setup (video disabled)
        % -----------------------------------------------------------------
        FrameSampleRate = 10;

        fig1 = figure(1);
        fig1.Position = [100 100 600 600];
        subfig1 = axes(fig1);
        hold(subfig1, 'on');
        axis(subfig1, 'equal');
        colormap(subfig1, flip(gray));
        clim(subfig1, [-max(new_weight{1})/5, max(new_weight{1})*2]);
        set(subfig1,'FontSize',18);
        xlabel(subfig1,'x [m]');
        ylabel(subfig1,'y [m]');

        % Color set for multi-agent visualization
        c_stand = { [0 0.4470 0.7410], [0.8500 0.3250 0.0980], ...
                    [0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560], ...
                    [0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330], ...
                    [0.6350 0.0780 0.1840] };

        %% ----------------------------------------------------------------
        %  Main iteration loop (receding-horizon D2OC)
        %  • Solve local optimal control
        %  • Move agent
        %  • Update weights
        %  • Share weights if within communication range
        % -----------------------------------------------------------------
        iter = 0;
        while (sum(weight_sum_check_index) > 0) && (iter < bat_life)
            iter = iter + 1;

            for n = 1:no_of_agents
                if weight_sum_check_index(n) > 0

                    % ----------------------------------------------------
                    %  Compute optimal target point xg for this iteration
                    %  via Hamiltonian-based D2OC (sampling-based)
                    % ----------------------------------------------------
                    xg = hamilton_optimal_control_Lin_Quad_R3( ...
                        [x{n}(1), x{n}(5)], alpha_k, y, ...
                        new_weight{n}, PosWgt_idx{n}, dist_sq{n});

                    % ----------------------------------------------------
                    %  Compute optimal control input via KKT solution
                    % ----------------------------------------------------
                    B_aug = zeros(size(A_aug,1),1);
                    for cnt_Hor = 1:Hor_Leng
                        B_aug((cnt_Hor-1)*dim_sv+1:cnt_Hor*dim_sv) = alpha_k*C_tr*xg;
                    end
                    B_aug(Hor_Leng*dim_sv+1:Hor_Leng*dim_sv+dim_sv) = -Ad*x{n};

                    % Solve KKT system (finite-horizon optimal control)
                    x_aug = A_aug_inv * B_aug;

                    % Extract control input (first horizon control)
                    u(:,n) = x_aug(dim_sv*(2*Hor_Leng)+1 : dim_sv*(2*Hor_Leng)+dim_iv);
                    u(:,n) = min(umax, max(u(:,n), -umax));   % saturation

                    % ----------------------------------------------------
                    %  Apply dynamics and update state
                    % ----------------------------------------------------
                    x{n} = Ad*x{n} + Bd*u(:,n);

                    % Velocity saturation correction (preserve direction)
                    if ~((abs(x{n}(2))==0)&&(abs(x{n}(6))==0))
                        if abs(x{n}(2)) > abs(x{n}(6))
                            cal_vel = x{n}(2);
                            x{n}(2) = min(vmax, max(x{n}(2), -vmax));
                            x{n}(6) = x{n}(6)/cal_vel*x{n}(2);
                        else
                            cal_vel = x{n}(6);
                            x{n}(6) = min(vmax, max(x{n}(6), -vmax));
                            x{n}(2) = x{n}(2)/cal_vel*x{n}(6);
                        end
                    end

                    % Angle limits
                    x{n}(3) = min(max(x{n}(3), -thresh_ang),    thresh_ang);
                    x{n}(7) = min(max(x{n}(7), -thresh_ang),    thresh_ang);

                    % Angular velocity limits
                    x{n}(4) = min(max(x{n}(4), -thresh_angvel), thresh_angvel);
                    x{n}(8) = min(max(x{n}(8), -thresh_angvel), thresh_angvel);

                    % Log trajectory and target
                    x_traj{n}(:,iter+1) = x{n};
                    x_goal{n}(:,iter+1) = xg;
                    u_traj{n}(:,iter+1) = u(:,n);

                    % ----------------------------------------------------
                    %  Weight update: consume local weight based on sampling
                    % ----------------------------------------------------
                    [new_weight{n}, weight_distribution, PosWgt_idx{n}, dist_sq{n}] = ...
                        update_weight_R2([x{n}(1),x{n}(5)], y, new_weight{n}, alpha_k, N);
                    agent_wgt{n,n} = agent_wgt{n,n} + weight_distribution;

                    % ----------------------------------------------------
                    %  Decentralized weight sharing (if close in distance)
                    %  - Share load only when communication range satisfied
                    % ----------------------------------------------------
                    comm_time{n}(n) = iter+1;

                    for q = 1:no_of_agents
                        if (q ~= n) && (norm(x{n}([1 5],1) - x{q}([1 5],1)) < comm_thresh)
                            comm_time_comp = (comm_time{n}-comm_time{q});
                            comm_req_idx = find(abs(comm_time_comp)>comm_time_thresh);

                            if ~isempty(comm_req_idx)
                                for cnt_ag = 1:length(comm_req_idx)
                                    idx = comm_req_idx(cnt_ag);
                                    if (comm_time{n}(idx)>comm_time{q}(idx))
                                        new_weight{q} = new_weight{q}+agent_wgt{q,idx};
                                        agent_wgt{q,idx} = agent_wgt{n,idx};
                                        comm_time{q}(idx) = comm_time{n}(idx);
                                        new_weight{q} = new_weight{q}-agent_wgt{q,idx};
                                    else
                                        new_weight{n} = new_weight{n}+agent_wgt{n,idx};
                                        agent_wgt{n,idx} = agent_wgt{q,idx};
                                        comm_time{n}(idx) = comm_time{q}(idx);
                                        new_weight{n} = new_weight{n}-agent_wgt{n,idx};
                                    end
                                end

                                % Ensure weight non-negativity
                                new_weight{n}(new_weight{n}<0) = 0;
                                new_weight{q}(new_weight{q}<0) = 0;

                                weight_sum_check_index(q) = sum(new_weight{q});
                            end
                        end
                    end

                    weight_sum_check_index(n) = sum(new_weight{n}); % local exhaustion

                else
                    % If an agent has no remaining weight to process, remove from map
                    x{n}([1 5]) = [inf; inf];
                end
            end  % agent loop

            %% ------------------------------------------------------------
            %  Live plot every FrameSampleRate iterations
            %% ------------------------------------------------------------
            if rem(iter,FrameSampleRate)==0

                % Scatter sample weights
                if exist('subfig1_sp','var')
                    delete(subfig1_sp);
                end
                subfig1_sp = scatter(subfig1, y_ref(:,1),y_ref(:,2),10,new_weight{1},'filled');

                % Delete previous robot plots
                if exist('subfig1_p1','var')
                    for n = 1:no_of_agents
                        delete(subfig1_p1{n});
                        delete(subfig1_p2{n});
                        delete(subfig1_p3{n});
                        delete(subfig1_p4{n});
                    end
                end

                % Draw trajectories, start, end, and goal markers
                for n = 1:no_of_agents
                    subfig1_p1{n} = plot(subfig1, x_traj{n}(1,:), x_traj{n}(5,:), ...
                        'LineWidth',2,'LineStyle','-','Color',c_stand{rem(n,7)+1});
                    subfig1_p2{n} = plot(subfig1, x_traj{n}(1,end), x_traj{n}(5,end), ...
                        'kO','MarkerFaceColor','y');
                    subfig1_p3{n} = plot(subfig1, x_traj{n}(1,1), x_traj{n}(5,1), ...
                        'bX','MarkerFaceColor','y','MarkerSize',15,'LineWidth',2);
                    subfig1_p4{n} = plot(subfig1, x_goal{n}(1,end), x_goal{n}(2,end), ...
                        'X','MarkerFaceColor','y','MarkerSize',12,'LineWidth',2, ...
                        'Color',c_stand{rem(n,7)+1});
                end

                % Adjust axes
                axis(subfig1,"equal");
                margin = (domain(2)-domain(1))/10;
                xlim(subfig1,[domain(1)-margin, domain(2)+margin]);
                ylim(subfig1,[domain(3)-margin, domain(4)+margin]);
                subfig1.XTick = domain(1)-margin : margin : domain(2)+margin;
                subfig1.YTick = domain(3)-margin : margin : domain(4)+margin;
                grid(subfig1,"on");
                drawnow;
            end

        end % while loop

        %% ----------------------------------------------------------------
        %  Wrap-up for this simulation: store results back into sim_result
        % -----------------------------------------------------------------
        m_rem_weight(iter) = mean(weight_sum_check_index);
        tic_end = toc(tic_start);
        ws_dist = 0;    % placeholder (Wasserstein distance disabled)

        fprintf('comp_time : %f, wass_dist : %f\n', tic_end, ws_dist);

        % Save metrics back into sim_result structure
        sim_result.tests.inputs{cnt_sim,sim_dic("comp_time")} = tic_end;
        sim_result.tests.inputs{cnt_sim,sim_dic("Wass_dist")} = ws_dist;

        if sim_dic.isKey("m_rem_weight")
            sim_result.tests.inputs{cnt_sim,sim_dic("m_rem_weight")} = m_rem_weight;
        end
        if sim_dic.isKey("agent_init")
            sim_result.tests.inputs{cnt_sim,sim_dic("agent_init")} = agent_init;
        end
        if sim_dic.isKey("Param")
            sim_result.tests.inputs{cnt_sim,sim_dic("Param")} = param_data;
        end
        if sim_dic.isKey("Term_time")
            for cnt_agent = 1:no_of_agents
                Term_time(cnt_agent) = length(x_traj{cnt_agent});
            end
            sim_result.tests.inputs{cnt_sim,sim_dic("Term_time")} = Term_time;
        end
        sim_result.tests.inputs{cnt_sim,sim_dic("complete")} = true;

        % Periodic save
        if rem(cnt_sim, Data_save_Int) == 0
            save(fullfile(dir, sim_file),"sim_result",'-append');
        end

    end % if ~complete

% final save (disabled by default)
% save(fullfile(dir,sim_file),"sim_result",'-append');
