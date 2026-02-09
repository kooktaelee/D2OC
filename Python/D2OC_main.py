import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from scipy.stats import gaussian_kde 

plt.close('all')

# [Hamiltonian-based Optimal Control]
def hamilton_optimal_control(x, alpha_k, z, local_weight, pos_idx, d_sq):
    """
    Computes the optimal target centroid (xg) using the Optimal Transport (OT) principle.
    The goal is to move towards the region that minimizes the 2-Wasserstein cost.
    
    Parameters:
    - x: Current agent position [x, y]
    - alpha_k: Agent's consumption capacity per time step
    - z: Full set of reference density points (coordinates)
    - local_weight: Current belief of remaining density mass at each point
    - pos_idx: Indices where remaining density weight is non-zero
    - d_sq: Squared distances from agent to points at pos_idx
    """
    if len(pos_idx) == 0: return None, 0
    
    # 1. Extract weights of the remaining density points within the agent's perception
    w_sub = local_weight[pos_idx]
    
    # 2. Compute the 'weighted normalized distance' to identify the most efficient coverage targets
    # This reflects the Hamiltonian-based cost: balancing proximity and density importance
    d_wnE = d_sq / (w_sub**2 + 1e-12)
    
    # 3. Sort points by the energy-distance to prioritize points closest to the agent's 'reach'
    idx_sorted = np.argsort(d_wnE)
    sorted_pos_idx = pos_idx[idx_sorted]
    sorted_w = local_weight[sorted_pos_idx]
    
    # 4. Use cumulative sum to find the specific points that fit within the agent's capacity (alpha_k)
    cumsum_w = np.cumsum(sorted_w)
    cutoff = np.searchsorted(cumsum_w, alpha_k)
    
    # 5. Boundary condition for indexing to prevent overflow
    if cutoff >= len(sorted_pos_idx): cutoff = len(sorted_pos_idx) - 1
    
    # 6. Identify the specific subset of points to be 'targeted' in this control step
    final_idx = sorted_pos_idx[:cutoff+1]
    final_w = sorted_w[:cutoff+1].copy()
    
    # 7. Fractional adjustment: Ensure the sum of targeted weights matches exactly alpha_k
    if cutoff > 0: final_w[-1] = alpha_k - cumsum_w[cutoff-1]
    else: final_w[-1] = alpha_k
    
    # 8. Calculate the centroid (xg) of these points, acting as the 'Hamiltonian Optimal' target
    y_loc_pos = z[final_idx]
    y_loc_wgt = final_w.reshape(-1, 1)
    gamma = np.sum(y_loc_wgt)
    xg = np.sum(y_loc_wgt * y_loc_pos, axis=0) / (gamma + 1e-12)
    return xg, gamma

# [Local Weight Update - Absolute Preservation]
def update_weight_local(x, z, local_weight, alpha_k):
    """
    Simulates the 'coverage' or 'consumption' of the reference distribution.
    Updates the agent's local belief of the remaining density map based on its position.
    """
    # 1. Identify indices where density mass still remains (non-zero weight)
    pos_idx = np.where(local_weight > 1e-15)[0]
    if len(pos_idx) == 0: return local_weight, np.zeros_like(local_weight), pos_idx, np.array([])
    
    # 2. Calculate squared Euclidean distances for nearest-neighbor consumption priority
    dist_sq = np.sum((z[pos_idx] - x)**2, axis=1)
    idx_sort = np.argsort(dist_sq)
    sorted_idx = pos_idx[idx_sort]
    sorted_w = local_weight[sorted_idx]
    cumsum_w = np.cumsum(sorted_w)
    
    # 3. Find points to be 'consumed' by the agent's presence according to alpha_k
    cutoff = np.searchsorted(cumsum_w, alpha_k)
    weight_dist = np.zeros_like(local_weight)
    
    # 4. Update logic: Deduct probability mass from the density map (simulating task completion)
    if cutoff < len(sorted_idx):
        take_idx = sorted_idx[:cutoff]
        weight_dist[take_idx] = local_weight[take_idx]; local_weight[take_idx] = 0
        remainder = alpha_k - (cumsum_w[cutoff-1] if cutoff > 0 else 0)
        weight_dist[sorted_idx[cutoff]] = remainder; local_weight[sorted_idx[cutoff]] -= remainder
    else:
        # If remaining mass is less than alpha_k, consume all available mass
        weight_dist[sorted_idx] = local_weight[sorted_idx]; local_weight[sorted_idx] = 0
    
    # 5. Clean up infinitesimal values for numerical stability and update neighbor info
    local_weight[local_weight < 1e-12] = 0
    new_pos_idx = np.where(local_weight > 1e-15)[0]
    new_dist_sq = np.sum((z[new_pos_idx] - x)**2, axis=1) if len(new_pos_idx) > 0 else np.array([])
    return local_weight, weight_dist, new_pos_idx, new_dist_sq

def run_d2oc_decentralized():
    # --- Simulation Hyperparameters ---
    no_of_agents, bat_life, hor_leng = 10, 100000, 5
    T, g, comm_range = 0.1, 9.81, 15.0 
    alpha_k = 1.0 / bat_life # Individual coverage capacity per step
    frames = []

    # --- System Dynamics: Quadrotor LTI Model (8 States) ---
    # Sub-system: [pos, vel, angle, angle_rate] for X and Y independently
    Ad_sub = np.array([[1, T, 0, 0], [0, 0.7, -T*g, 0], [0, 0, 1, T], [0, 0, 0, 0.7]]) 
    Ad = np.zeros((8, 8)); Ad[:4, :4] = Ad_sub; Ad[4:, 4:] = Ad_sub
    Bd = np.zeros((8, 2)); Bd[3, 0] = T/0.1; Bd[7, 1] = T/0.1
    # CTC: Output matrix to extract position states (X at index 0, Y at index 4)
    CTC = np.zeros((8, 8)); CTC[0,0] = 1.0; CTC[4,4] = 1.0
    # MPC Cost matrices (Q: State error weight, R: Control effort weight)
    Q_base, R_block = np.eye(8) * 0.1, np.eye(2) * 10.0 

    # --- MPC Pre-computation: KKT Matrix Formulation ---
    # Pre-calculating matrix inversions for the prediction horizon to ensure real-time performance
    E12 = np.zeros((40, 40)) # State evolution constraint matrix
    for i in range(hor_leng):
        E12[i*8:(i+1)*8, i*8:(i+1)*8] = -np.eye(8)
        if i < hor_leng - 1: E12[i*8:(i+1)*8, (i+1)*8:(i+2)*8] = Ad.T
    E12_inv = np.linalg.inv(E12)
    E23 = np.kron(np.eye(hor_leng), Bd) # Input influence matrix
    E33 = np.kron(np.eye(hor_leng), R_block) # Control cost matrix
    E12_inv_E23 = E12_inv @ E23
    E23T_E12inv = E23.T @ E12_inv

    # --- Generate Reference Distribution (Target Density) ---
    # Multi-modal Gaussian distribution representing the 'interest map'
    y_ref = np.vstack([np.random.multivariate_normal(np.random.uniform(8, 92, 2), np.eye(2)*10, 300) for _ in range(15)])
    
    # --- Visualization: Contour Generation via KDE ---
    xi, yi = np.mgrid[0:100:100j, 0:100:100j]
    positions = np.vstack([xi.ravel(), yi.ravel()])
    kde = gaussian_kde(y_ref.T)
    zi = np.reshape(kde(positions).T, xi.shape)

    # --- Decentralized State Variables ---
    agent_betas = [np.ones(len(y_ref)) / len(y_ref) for _ in range(no_of_agents)] # Local belief of map
    agent_wgt_history = np.zeros((no_of_agents, no_of_agents, len(y_ref))) # Consensus weight history
    comm_time, agents = np.zeros((no_of_agents, no_of_agents)), np.zeros((no_of_agents, 8))
    trajectories = [[] for _ in range(no_of_agents)]

    # Initial Random Deployment of Agents
    initial_positions = []
    for n in range(no_of_agents):
        pos = np.random.uniform(10, 90, 2)
        agents[n, [0, 4]] = pos
        initial_positions.append(pos.copy())

    # --- Plotting Setup ---
    plt.ion()
    fig, ax = plt.subplots(figsize=(9, 9))
    
    # Static Background: Probability density contour lines
    ax.contour(xi, yi, zi, levels=20, colors='black', linewidths=0.3, alpha=0.35, zorder=0)
    ax.grid(True, color='gray', linestyle=':', linewidth=0.75, alpha=0.8)
    ax.set_axisbelow(True) 

    # Scatter plots for map points (Gray: Consumed, Blue: Remaining)
    scat_consumed = ax.scatter([], [], c='#E5E5E5', s=8, alpha=0.8) 
    scat_remain = ax.scatter([], [], c='#007BFF', s=5, alpha=0.1, zorder=5) 
    
    # Dynamic Plotting Handles for paths, agents, and communication links
    full_paths, tail_paths, agent_plots, agent_circles, comm_lines = [], [], [], [], []
    for n in range(no_of_agents):
        color = plt.cm.tab10(n)
        f_line, = ax.plot([], [], color=color, lw=1, alpha=0.9, zorder=6) # Full path
        full_paths.append(f_line)
        t_line, = ax.plot([], [], color=color, lw=5, alpha=0.5, zorder=8) # Short moving tail
        tail_paths.append(t_line)
        
        # Initial position marker
        ax.plot(initial_positions[n][0], initial_positions[n][1], 'x', color='black', ms=8, alpha=0.4, zorder=3)
        
        point, = ax.plot([], [], 'o', color=color, ms=11, mec='black', mew=0.5, zorder=15) # UAV pos
        agent_plots.append(point)
        circle = plt.Circle((0, 0), comm_range, color=color, fill=True, alpha=0.07, zorder=2) # Communication range
        ax.add_artist(circle)
        agent_circles.append(circle)
        comm_line, = ax.plot([], [], color='orange', lw=2, alpha=0, zorder=12) # Active link indicator
        comm_lines.append(comm_line)

    ax.set_xlim(0, 100); ax.set_ylim(0, 100)
    title_text = ax.set_title("D2OC: Decentralized Coverage with Density Contour", fontsize=14, fontweight='bold')

    
    max_vel = 0.5 # Physical velocity limit for safety
    
    # --- Main Simulation Loop ---
    for t in range(1, bat_life + 1):
        for cl in comm_lines: cl.set_alpha(0) # Reset communication visuals

        for n in range(no_of_agents):
            # 1. State History and Local Task Consumption
            trajectories[n].append(agents[n, [0, 4]].copy())
            agent_betas[n], delta_w, pos_idx, d_sq = update_weight_local(agents[n, [0, 4]], y_ref, agent_betas[n], alpha_k)
            agent_wgt_history[n, n] += delta_w
            comm_time[n, n] = t
            
            # 2. Hamiltonian Optimal Target Generation
            xg, gamma = hamilton_optimal_control(agents[n, [0, 4]], alpha_k, y_ref, agent_betas[n], pos_idx, d_sq)
            
            # 3. Fast MPC Solver (via KKT)
            if xg is not None and gamma > 1e-12:
                # Dynamically weight tracking priority based on the target's mass (gamma)
                Q_track = (gamma * 250.0) * CTC + Q_base 
                xg_full = np.zeros(8); xg_full[0], xg_full[4] = xg[0], xg[1]
                
                # Construct RHS/LHS for the horizon optimization
                F1 = np.tile(Q_track @ xg_full, hor_leng)
                F2 = np.zeros(40); F2[:8] = -Ad @ agents[n]
                M = E12_inv_E23
                E11_M = np.zeros_like(M)
                for i in range(hor_leng): E11_M[i*8:(i+1)*8, :] = Q_track @ M[i*8:(i+1)*8, :]
                
                lhs = E33 + M.T @ E11_M
                temp_F2_trans = E12_inv.T @ F2
                E11_temp_F2 = np.zeros(40)
                for i in range(hor_leng): E11_temp_F2[i*8:(i+1)*8] = Q_track @ temp_F2_trans[i*8:(i+1)*8]
                
                rhs = -(E23T_E12inv @ F1) + E23.T @ (E12_inv @ E11_temp_F2)
                u_star = np.linalg.solve(lhs, rhs)[:2] # Extract current optimal control
            else:
                # Perform braking maneuver if no valid targets are within reach
                u_star = -agents[n, [1, 5]] / T 
                agents[n, [1, 2, 3, 5, 6, 7]] = 0

            # 4. Integrate System Dynamics and Apply Actuator Saturations
            agents[n] = Ad @ agents[n] + Bd @ np.clip(u_star, -1.0, 1.0)
            agents[n, [0, 4]] = np.clip(agents[n, [0, 4]], 0.1, 99.9)
            vel = agents[n, [1, 5]]; spd = np.linalg.norm(vel)
            if spd > max_vel: agents[n, [1, 5]] = (vel / spd) * max_vel

        # 5. Decentralized Communication: Local Belief Synchronization
        for i in range(no_of_agents):
            d_sq_mat = np.sum((agents[:, [0, 4]] - agents[i, [0, 4]])**2, axis=1)
            neighbors = np.where((d_sq_mat < comm_range**2) & (np.arange(no_of_agents) != i))[0]
            for j in neighbors:
                # Exchange and update coverage information based on timestamps
                idx_upd = np.where(comm_time[i] > comm_time[j])[0]
                if len(idx_upd) > 0:
                    diffs = np.sum(agent_wgt_history[i, idx_upd] - agent_wgt_history[j, idx_upd], axis=0)
                    agent_betas[j] = np.maximum(0, agent_betas[j] - diffs)
                    agent_wgt_history[j, idx_upd] = agent_wgt_history[i, idx_upd].copy()
                    comm_time[j, idx_upd] = comm_time[i, idx_upd].copy()
                    # Visual feedback for active communication
                    if i < j:
                        comm_lines[i].set_data([agents[i,0], agents[j,0]], [agents[i,4], agents[j,4]])
                        comm_lines[i].set_alpha(0.7)

        # 6. Periodic Visual Update and Frame Capture
        if t % 50 == 0 or t == bat_life:
            consensus_rem = np.min(agent_betas, axis=0) # Global mission progress belief
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
            
            # Update mission progress in title
            progress = (1 - np.sum(consensus_rem)) * 100
            title_text.set_text(f"Coverage: {progress:.2f}% | Step: {t}")
            
            # Redraw canvas and capture frame for GIF
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
            image_rgba = renderer.buffer_rgba()
            frames.append(Image.fromarray(np.array(image_rgba)).convert("RGB"))
            
            plt.pause(0.001)
            if progress >= 99.99: break # Exit loop if coverage objective achieved

    # --- Save Animation as GIF ---
    # if frames:
    #     # Define frame durations: 60ms for motion, 2500ms for final pause
    #     durations = [60] * len(frames)
    #     durations[-1] = 2500 
    #     frames[0].save('simulation_refined.gif', save_all=True, append_images=frames[1:], duration=durations, loop=0)
    #     print("GIF saved: simulation_refined.gif")

    plt.ioff(); plt.show()

if __name__ == "__main__":

    run_d2oc_decentralized()
