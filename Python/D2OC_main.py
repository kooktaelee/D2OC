import numpy as np
import matplotlib.pyplot as plt
from PIL import Image  # Required for capturing simulation frames and exporting as GIF

# Close all existing matplotlib windows to start with a clean environment
plt.close('all')

# [Hamiltonian-based Optimal Control Function]
def hamilton_optimal_control(x, alpha_k, z, local_weight, pos_idx, d_sq):
    """
    Computes the optimal centroid (xg) that minimizes the Wasserstein-2 distance.
    This is the core of the 'Density-Driven' aspect, selecting where the agent should move
    to best cover the remaining reference distribution.
    """
    if len(pos_idx) == 0: return None, 0
    
    # Extract weights of remaining density points in the local neighborhood
    w_sub = local_weight[pos_idx]
    
    # Calculate energy-based distance: squared distance divided by weights
    # This guides the OT potential field computation
    d_wnE = d_sq / (w_sub**2 + 1e-12)
    
    # Sort density points to identify the most 'cost-effective' targets for coverage
    idx_sorted = np.argsort(d_wnE)
    sorted_pos_idx = pos_idx[idx_sorted]
    sorted_w = local_weight[sorted_pos_idx]
    
    # Cumulative sum to find the cutoff point that matches the agent's 'consumption capacity' (alpha_k)
    cumsum_w = np.cumsum(sorted_w)
    cutoff = np.searchsorted(cumsum_w, alpha_k)
    
    # Boundary handling for the cutoff index
    if cutoff >= len(sorted_pos_idx): cutoff = len(sorted_pos_idx) - 1
    
    # Select the subset of density points for the current time-step coverage
    final_idx = sorted_pos_idx[:cutoff+1]
    final_w = sorted_w[:cutoff+1].copy()
    
    # Fractional weight adjustment: Ensure the sum of selected weights equals exactly alpha_k
    if cutoff > 0: final_w[-1] = alpha_k - cumsum_w[cutoff-1]
    else: final_w[-1] = alpha_k
    
    # Compute the target position (xg) as the weighted centroid of selected points
    y_loc_pos = z[final_idx]
    y_loc_wgt = final_w.reshape(-1, 1)
    gamma = np.sum(y_loc_wgt)
    
    # The resulting xg is the 'optimal point' in the OT sense
    xg = np.sum(y_loc_wgt * y_loc_pos, axis=0) / (gamma + 1e-12)
    return xg, gamma

def update_weight_local(x, z, local_weight, alpha_k):
    """
    Updates the agent's local knowledge of the density field.
    Simulates 'consumption' or 'coverage' of the area by subtracting probability mass
    from the reference distribution based on the agent's current position.
    """
    pos_idx = np.where(local_weight > 1e-15)[0] # Find points with remaining mass
    if len(pos_idx) == 0: return local_weight, np.zeros_like(local_weight), pos_idx, np.array([])
    
    # Calculate Euclidean distance from the agent to all density points
    dist_sq = np.sum((z[pos_idx] - x)**2, axis=1)
    idx_sort = np.argsort(dist_sq)
    sorted_idx = pos_idx[idx_sort]
    sorted_w = local_weight[sorted_idx]
    cumsum_w = np.cumsum(sorted_w)
    
    # Find points nearest to the agent to 'consume' the density mass
    cutoff = np.searchsorted(cumsum_w, alpha_k)
    weight_dist = np.zeros_like(local_weight)
    
    # Update local density map: subtract alpha_k mass from the nearest points
    if cutoff < len(sorted_idx):
        take_idx = sorted_idx[:cutoff]
        weight_dist[take_idx] = local_weight[take_idx]; local_weight[take_idx] = 0
        remainder = alpha_k - (cumsum_w[cutoff-1] if cutoff > 0 else 0)
        weight_dist[sorted_idx[cutoff]] = remainder; local_weight[sorted_idx[cutoff]] -= remainder
    else:
        weight_dist[sorted_idx] = local_weight[sorted_idx]; local_weight[sorted_idx] = 0
    
    # Numerical thresholding to keep the weight map sparse and clean
    local_weight[local_weight < 1e-12] = 0
    new_pos_idx = np.where(local_weight > 1e-15)[0]
    
    # Return updated weights and distance data for the next control cycle
    new_dist_sq = np.sum((z[new_pos_idx] - x)**2, axis=1) if len(new_pos_idx) > 0 else np.array([])
    return local_weight, weight_dist, new_pos_idx, new_dist_sq

def run_d2oc_decentralized():
    # Simulation Parameters
    no_of_agents = 10     # Total number of UAVs/Agents
    bat_life = 50000     # Maximum simulation steps (battery life)
    hor_leng = 5          # MPC Prediction Horizon length
    T = 0.1               # Sampling time interval (seconds)
    g = 9.81              # Gravity constant for dynamics
    comm_range = 15.0     # Maximum communication radius between agents
    alpha_k = 1.0 / bat_life # Mass consumed per agent per step
    
    frames = [] # Buffer for GIF frames

    # --- Discrete-time System Dynamics (8-state Quadrotor model) ---
    # States: [pos_x, vel_x, pitch, pitch_rate, pos_y, vel_y, roll, roll_rate]
    Ad_sub = np.array([[1, T, 0, 0], [0, 0.7, -T*g, 0], [0, 0, 1, T], [0, 0, 0, 0.7]]) 
    Ad = np.zeros((8, 8)); Ad[:4, :4] = Ad_sub; Ad[4:, 4:] = Ad_sub
    Bd = np.zeros((8, 2)); Bd[3, 0] = T/0.1; Bd[7, 1] = T/0.1 # Control inputs to attitude
    
    # State selection matrix for position tracking (X and Y coordinates)
    CTC = np.zeros((8, 8)); CTC[0,0] = 1.0; CTC[4,4] = 1.0
    
    # MPC Cost Weight Matrices (Q: State error, R: Control effort)
    Q_base, R_block = np.eye(8) * 0.1, np.eye(2) * 30.0 

    # --- KKT Matrix Pre-computation (Finite-Horizon MPC Solver) ---
    # This structure allows for solving the optimal control problem via linear algebra
    E12 = np.zeros((40, 40))
    for i in range(hor_leng):
        E12[i*8:(i+1)*8, i*8:(i+1)*8] = -np.eye(8)
        if i < hor_leng - 1: E12[i*8:(i+1)*8, (i+1)*8:(i+2)*8] = Ad.T
    E12_inv = np.linalg.inv(E12)
    E23 = np.kron(np.eye(hor_leng), Bd)
    E33 = np.kron(np.eye(hor_leng), R_block)
    E12_inv_E23 = E12_inv @ E23
    E23T_E12inv = E23.T @ E12_inv

    # Generate Target Reference Density (Multiple Gaussian Mixtures)
    y_ref = np.vstack([np.random.multivariate_normal(np.random.uniform(10, 90, 2), np.eye(2)*50, 300) for _ in range(10)])
    
    # Decentralized state variables: each agent maintains its own belief (beta) of the global density
    agent_betas = [np.ones(len(y_ref)) / len(y_ref) for _ in range(no_of_agents)]
    agent_wgt_history = np.zeros((no_of_agents, no_of_agents, len(y_ref))) # Track what others have covered
    comm_time = np.zeros((no_of_agents, no_of_agents)) # Timestamp for consensus
    agents = np.zeros((no_of_agents, 8)) # Physical states of all agents
    trajectories = [[] for _ in range(no_of_agents)]

    # Random Initial Deployment of Agents
    initial_positions = []
    for n in range(no_of_agents):
        pos = np.random.uniform(10, 90, 2)
        agents[n, [0, 4]] = pos
        initial_positions.append(pos.copy())

    # --- Plotting and UI Initialization ---
    plt.ion()
    fig, ax = plt.subplots(figsize=(9, 9))
    ax.set_facecolor('#FFFFFF')
    ax.grid(True, color='gray', linestyle=':', linewidth=0.75, alpha=0.8)
    
    # Scatters for visualized density: Gray = Covered, Blue = Remaining
    scat_consumed = ax.scatter([], [], c='#E5E5E5', s=8, alpha=0.8) 
    scat_remain = ax.scatter([], [], c='#007BFF', s=5, alpha=0.6, zorder=5) 
    
    full_paths, tail_paths, agent_plots, agent_circles, comm_lines = [], [], [], [], []
    for n in range(no_of_agents):
        color = plt.cm.tab10(n)
        f_line, = ax.plot([], [], color=color, lw=1, alpha=0.75, zorder=6) # Permanent path
        full_paths.append(f_line)
        t_line, = ax.plot([], [], color=color, lw=5, alpha=0.15, zorder=8) # Recent tail
        tail_paths.append(t_line)
        ax.plot(initial_positions[n][0], initial_positions[n][1], 'x', color='black', ms=8, alpha=0.4)
        point, = ax.plot([], [], 'o', color=color, ms=11, mec='black', mew=0.5, zorder=15) # Agent marker
        agent_plots.append(point)
        circle = plt.Circle((0, 0), comm_range, color=color, fill=True, alpha=0.07, zorder=2) # Comm range
        ax.add_artist(circle)
        agent_circles.append(circle)
        comm_line, = ax.plot([], [], color='orange', lw=2, alpha=0, zorder=12) # Comm link visual
        comm_lines.append(comm_line)

    ax.set_xlim(0, 100); ax.set_ylim(0, 100)
    title_text = ax.set_title("D2OC: Decentralized Coverage (Full Trajectory Mode)", fontsize=14, fontweight='bold')

    max_vel = 0.6 # Velocity limit for stability and realism
    
    # --- Main Time-Stepping Simulation ---
    for t in range(1, bat_life + 1):
        # Reset visual links
        for cl in comm_lines: cl.set_alpha(0)

        for n in range(no_of_agents):
            # 1. Update personal coverage and local weight belief
            trajectories[n].append(agents[n, [0, 4]].copy())
            agent_betas[n], delta_w, pos_idx, d_sq = update_weight_local(agents[n, [0, 4]], y_ref, agent_betas[n], alpha_k)
            agent_wgt_history[n, n] += delta_w
            comm_time[n, n] = t
            
            # 2. Determine Optimal Transport Target Position
            xg, gamma = hamilton_optimal_control(agents[n, [0, 4]], alpha_k, y_ref, agent_betas[n], pos_idx, d_sq)
            
            # 3. Compute MPC Control Input (KKT Solver)
            if xg is not None and gamma > 1e-12:
                # Dynamic Q weighting: gain increases as local density mass (gamma) increases
                Q_track = (gamma * 250.0) * CTC + Q_base 
                xg_full = np.zeros(8); xg_full[0], xg_full[4] = xg[0], xg[1]
                
                # Setup objective function vectors
                F1 = np.tile(Q_track @ xg_full, hor_leng)
                F2 = np.zeros(40); F2[:8] = -Ad @ agents[n]
                
                # Efficient matrix inversion-based solver for the horizon
                M = E12_inv_E23
                E11_M = np.zeros_like(M)
                for i in range(hor_leng): E11_M[i*8:(i+1)*8, :] = Q_track @ M[i*8:(i+1)*8, :]
                
                lhs = E33 + M.T @ E11_M
                temp_F2_trans = E12_inv.T @ F2
                E11_temp_F2 = np.zeros(40)
                for i in range(hor_leng): E11_temp_F2[i*8:(i+1)*8] = Q_track @ temp_F2_trans[i*8:(i+1)*8]
                
                rhs = -(E23T_E12inv @ F1) + E23.T @ (E12_inv @ E11_temp_F2)
                
                # Solve linear system for optimal control sequence u*
                u_star = np.linalg.solve(lhs, rhs)[:2]
            else:
                # If no targets remain, apply active braking
                u_star = -agents[n, [1, 5]] / T 
                agents[n, [1, 2, 3, 5, 6, 7]] = 0

            # 4. State Integration and Safety Constraints
            agents[n] = Ad @ agents[n] + Bd @ np.clip(u_star, -1.0, 1.0) # Apply dynamics
            agents[n, [0, 4]] = np.clip(agents[n, [0, 4]], 0.1, 99.9) # Boundary clipping
            
            # Physical velocity saturation
            vel = agents[n, [1, 5]]; spd = np.linalg.norm(vel)
            if spd > max_vel: agents[n, [1, 5]] = (vel / spd) * max_vel

        # 5. Decentralized Consensus (Local Neighbor Communication)
        for i in range(no_of_agents):
            # Find neighbors within the communication range
            d_sq_mat = np.sum((agents[:, [0, 4]] - agents[i, [0, 4]])**2, axis=1)
            neighbors = np.where((d_sq_mat < comm_range**2) & (np.arange(no_of_agents) != i))[0]
            
            for j in neighbors:
                # Sync logic: update neighbor with the most recent coverage information
                idx_upd = np.where(comm_time[i] > comm_time[j])[0]
                if len(idx_upd) > 0:
                    # Update the neighbor's weight map based on collective coverage data
                    diffs = np.sum(agent_wgt_history[i, idx_upd] - agent_wgt_history[j, idx_upd], axis=0)
                    agent_betas[j] = np.maximum(0, agent_betas[j] - diffs)
                    agent_wgt_history[j, idx_upd] = agent_wgt_history[i, idx_upd].copy()
                    comm_time[j, idx_upd] = comm_time[i, idx_upd].copy()
                    
                    # Visualize communication link
                    if i < j:
                        comm_lines[i].set_data([agents[i,0], agents[j,0]], [agents[i,4], agents[j,4]])
                        comm_lines[i].set_alpha(0.7)

        # 6. Periodic Visualization and Monitoring
        if t % 50 == 0 or t == bat_life:
            # Calculate global coverage consensus
            consensus_rem = np.min(agent_betas, axis=0)
            consumed = consensus_rem <= 1e-8
            scat_consumed.set_offsets(y_ref[consumed])
            scat_remain.set_offsets(y_ref[~consumed])
            
            for n in range(no_of_agents):
                path_data = np.array(trajectories[n])
                full_paths[n].set_data(path_data[:, 0], path_data[:, 1])
                tail_idx = max(0, len(path_data)-300)
                tail_paths[n].set_data(path_data[tail_idx:, 0], path_data[tail_idx:, 1])
                agent_plots[n].set_data([agents[n, 0]], [agents[n, 4]])
                agent_circles[n].center = (agents[n, 0], agents[n, 4])
            
            # Progress calculation based on total probability mass covered
            progress = (1 - np.sum(consensus_rem)) * 100
            title_text.set_text(f"Coverage: {progress:.2f}% | Step: {t}")
            
            # Capture frame for GIF production
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
            image_rgba = renderer.buffer_rgba()
            frames.append(Image.fromarray(np.array(image_rgba)).convert("RGB"))
            
            plt.pause(0.001)
            # Exit if mission objective (full coverage) is achieved
            if progress >= 99.99: break

    # --- Finalize and Save GIF ---
    # if frames:
    #     print(f"Finalizing GIF with {len(frames)} frames...")
    #     durations = [60] * len(frames)
    #     durations[-1] = 2500 # Hold last frame for 2.5 seconds
        
    #     frames[0].save('simulation.gif', 
    #                    save_all=True, 
    #                    append_images=frames[1:], 
    #                    optimize=False, 
    #                    duration=durations, 
    #                    loop=0)
    #     print("Success: simulation.gif saved.")

    plt.ioff(); plt.show()

# Execution Entry Point
if __name__ == "__main__":
    run_d2oc_decentralized()