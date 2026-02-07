import numpy as np
import matplotlib.pyplot as plt
from PIL import Image  # GIF 저장용

plt.close('all')

# [기존 연산 함수 - 절대 보존]
def hamilton_optimal_control(x, alpha_k, z, local_weight, pos_idx, d_sq):
    if len(pos_idx) == 0: return None, 0
    w_sub = local_weight[pos_idx]
    d_wnE = d_sq / (w_sub**2 + 1e-12)
    idx_sorted = np.argsort(d_wnE)
    sorted_pos_idx = pos_idx[idx_sorted]
    sorted_w = local_weight[sorted_pos_idx]
    cumsum_w = np.cumsum(sorted_w)
    cutoff = np.searchsorted(cumsum_w, alpha_k)
    if cutoff >= len(sorted_pos_idx): cutoff = len(sorted_pos_idx) - 1
    final_idx = sorted_pos_idx[:cutoff+1]
    final_w = sorted_w[:cutoff+1].copy()
    if cutoff > 0: final_w[-1] = alpha_k - cumsum_w[cutoff-1]
    else: final_w[-1] = alpha_k
    y_loc_pos = z[final_idx]
    y_loc_wgt = final_w.reshape(-1, 1)
    gamma = np.sum(y_loc_wgt)
    xg = np.sum(y_loc_wgt * y_loc_pos, axis=0) / (gamma + 1e-12)
    return xg, gamma

def update_weight_local(x, z, local_weight, alpha_k):
    pos_idx = np.where(local_weight > 1e-15)[0]
    if len(pos_idx) == 0: return local_weight, np.zeros_like(local_weight), pos_idx, np.array([])
    dist_sq = np.sum((z[pos_idx] - x)**2, axis=1)
    idx_sort = np.argsort(dist_sq)
    sorted_idx = pos_idx[idx_sort]
    sorted_w = local_weight[sorted_idx]
    cumsum_w = np.cumsum(sorted_w)
    cutoff = np.searchsorted(cumsum_w, alpha_k)
    weight_dist = np.zeros_like(local_weight)
    if cutoff < len(sorted_idx):
        take_idx = sorted_idx[:cutoff]
        weight_dist[take_idx] = local_weight[take_idx]; local_weight[take_idx] = 0
        remainder = alpha_k - (cumsum_w[cutoff-1] if cutoff > 0 else 0)
        weight_dist[sorted_idx[cutoff]] = remainder; local_weight[sorted_idx[cutoff]] -= remainder
    else:
        weight_dist[sorted_idx] = local_weight[sorted_idx]; local_weight[sorted_idx] = 0
    local_weight[local_weight < 1e-12] = 0
    new_pos_idx = np.where(local_weight > 1e-15)[0]
    new_dist_sq = np.sum((z[new_pos_idx] - x)**2, axis=1) if len(new_pos_idx) > 0 else np.array([])
    return local_weight, weight_dist, new_pos_idx, new_dist_sq

def run_d2oc_decentralized():
    no_of_agents, bat_life, hor_leng = 10, 100000, 5
    T, g, comm_range = 0.1, 9.81, 15.0 
    alpha_k = 1.0 / bat_life
    
    # --- GIF 저장용 리스트 ---
    frames = []

    Ad_sub = np.array([[1, T, 0, 0], [0, 0.7, -T*g, 0], [0, 0, 1, T], [0, 0, 0, 0.7]]) 
    Ad = np.zeros((8, 8)); Ad[:4, :4] = Ad_sub; Ad[4:, 4:] = Ad_sub
    Bd = np.zeros((8, 2)); Bd[3, 0] = T/0.1; Bd[7, 1] = T/0.1
    CTC = np.zeros((8, 8)); CTC[0,0] = 1.0; CTC[4,4] = 1.0
    Q_base, R_block = np.eye(8) * 0.1, np.eye(2) * 30.0 

    E12 = np.zeros((40, 40))
    for i in range(hor_leng):
        E12[i*8:(i+1)*8, i*8:(i+1)*8] = -np.eye(8)
        if i < hor_leng - 1: E12[i*8:(i+1)*8, (i+1)*8:(i+2)*8] = Ad.T
    E12_inv = np.linalg.inv(E12)
    E23 = np.kron(np.eye(hor_leng), Bd)
    E33 = np.kron(np.eye(hor_leng), R_block)
    E12_inv_E23 = E12_inv @ E23
    E23T_E12inv = E23.T @ E12_inv

    y_ref = np.vstack([np.random.multivariate_normal(np.random.uniform(10, 90, 2), np.eye(2)*50, 300) for _ in range(10)])
    agent_betas = [np.ones(len(y_ref)) / len(y_ref) for _ in range(no_of_agents)]
    agent_wgt_history = np.zeros((no_of_agents, no_of_agents, len(y_ref)))
    comm_time, agents = np.zeros((no_of_agents, no_of_agents)), np.zeros((no_of_agents, 8))
    trajectories = [[] for _ in range(no_of_agents)]

    initial_positions = []
    for n in range(no_of_agents):
        pos = np.random.uniform(10, 90, 2)
        agents[n, [0, 4]] = pos
        initial_positions.append(pos.copy())

    # --- 시각화 설정 ---
    plt.ion()
    fig, ax = plt.subplots(figsize=(9, 9))
    ax.set_facecolor('#FFFFFF')
    ax.grid(True, color='gray', linestyle=':', linewidth=0.75, alpha=0.8)
    ax.set_axisbelow(True)
    
    scat_consumed = ax.scatter([], [], c='#E5E5E5', s=8, alpha=0.8) 
    scat_remain = ax.scatter([], [], c='#007BFF', s=5, alpha=0.6, zorder=5) 
    
    full_paths, tail_paths, agent_plots, agent_circles, comm_lines = [], [], [], [], []
    for n in range(no_of_agents):
        color = plt.cm.tab10(n)
        f_line, = ax.plot([], [], color=color, lw=1, alpha=0.75, zorder=6)
        full_paths.append(f_line)
        t_line, = ax.plot([], [], color=color, lw=5, alpha=0.15, zorder=8)
        tail_paths.append(t_line)
        ax.plot(initial_positions[n][0], initial_positions[n][1], 'x', color='black', ms=8, alpha=0.4)
        point, = ax.plot([], [], 'o', color=color, ms=11, mec='black', mew=0.5, zorder=15)
        agent_plots.append(point)
        circle = plt.Circle((0, 0), comm_range, color=color, fill=True, alpha=0.07, zorder=2)
        ax.add_artist(circle)
        agent_circles.append(circle)
        comm_line, = ax.plot([], [], color='orange', lw=2, alpha=0, zorder=12)
        comm_lines.append(comm_line)

    ax.set_xlim(0, 100); ax.set_ylim(0, 100)
    title_text = ax.set_title("D2OC: Decentralized Coverage (Full Trajectory Mode)", fontsize=14, fontweight='bold')

    max_vel = 0.6
    for t in range(1, bat_life + 1):
        for cl in comm_lines: cl.set_alpha(0)

        for n in range(no_of_agents):
            trajectories[n].append(agents[n, [0, 4]].copy())
            agent_betas[n], delta_w, pos_idx, d_sq = update_weight_local(agents[n, [0, 4]], y_ref, agent_betas[n], alpha_k)
            agent_wgt_history[n, n] += delta_w
            comm_time[n, n] = t
            xg, gamma = hamilton_optimal_control(agents[n, [0, 4]], alpha_k, y_ref, agent_betas[n], pos_idx, d_sq)
            
            if xg is not None and gamma > 1e-12:
                Q_track = (gamma * 250.0) * CTC + Q_base 
                xg_full = np.zeros(8); xg_full[0], xg_full[4] = xg[0], xg[1]
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
                u_star = np.linalg.solve(lhs, rhs)[:2]
            else:
                u_star = -agents[n, [1, 5]] / T 
                agents[n, [1, 2, 3, 5, 6, 7]] = 0

            agents[n] = Ad @ agents[n] + Bd @ np.clip(u_star, -1.0, 1.0)
            agents[n, [0, 4]] = np.clip(agents[n, [0, 4]], 0.1, 99.9)
            vel = agents[n, [1, 5]]; spd = np.linalg.norm(vel)
            if spd > max_vel: agents[n, [1, 5]] = (vel / spd) * max_vel

        for i in range(no_of_agents):
            d_sq_mat = np.sum((agents[:, [0, 4]] - agents[i, [0, 4]])**2, axis=1)
            neighbors = np.where((d_sq_mat < comm_range**2) & (np.arange(no_of_agents) != i))[0]
            for j in neighbors:
                idx_upd = np.where(comm_time[i] > comm_time[j])[0]
                if len(idx_upd) > 0:
                    diffs = np.sum(agent_wgt_history[i, idx_upd] - agent_wgt_history[j, idx_upd], axis=0)
                    agent_betas[j] = np.maximum(0, agent_betas[j] - diffs)
                    agent_wgt_history[j, idx_upd] = agent_wgt_history[i, idx_upd].copy()
                    comm_time[j, idx_upd] = comm_time[i, idx_upd].copy()
                    if i < j:
                        comm_lines[i].set_data([agents[i,0], agents[j,0]], [agents[i,4], agents[j,4]])
                        comm_lines[i].set_alpha(0.7)

        if t % 50 == 0 or t == bat_life:
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
            
            progress = (1 - np.sum(consensus_rem)) * 100
            title_text.set_text(f"Coverage: {progress:.2f}% | Step: {t}")
            
            # [수정된 프레임 캡처 로직]
            fig.canvas.draw()
            # 캔버스의 renderer에서 직접 이미지를 가져옵니다 (가장 깨끗한 방법)
            renderer = fig.canvas.get_renderer()
            image_rgba = renderer.buffer_rgba()
            frames.append(Image.fromarray(np.array(image_rgba)).convert("RGB"))
            
            plt.pause(0.001)
            if progress >= 99.99: break

    # --- [GIF 저장: 종료 직후 수행] ---
    if frames:
        print(f"Saving {len(frames)} frames to GIF...")
        
        # 각 프레임의 지속 시간 리스트 생성 (기본 60ms)
        durations = [60] * len(frames)
        
        # 마지막 프레임만 3000ms(3초)로 설정하여 멈춤 효과 부여
        durations[-1] = 2500 
        
        frames[0].save('simulation.gif', 
                       save_all=True, 
                       append_images=frames[1:], 
                       optimize=False, 
                       duration=durations, # 리스트로 전달 가능
                       loop=0)
        print("GIF saved: simulation.gif (with 3s end delay)")

    plt.ioff(); plt.show()

if __name__ == "__main__":
    run_d2oc_decentralized()