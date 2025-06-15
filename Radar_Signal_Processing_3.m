clc
clear
close all

Is_plot = 0;


% % 动画方式逐帧播放目标轨迹
% load('target_traj_data.mat');
% smooth_and_plot_trajectory(target_traj, frame_index_list);
% 

if Is_plot == 1

    load('pointCloudList.mat');

    % DBSCAN 参数设置
    eps = 0.5;    % 距离阈值（单位：m）
    minPts = 3;   % 最小聚类点数
    
    % 运行 DBSCAN 去噪
    denoisedPointCloudList = run_dbscan_denoise(pointCloudList, eps, minPts);

    % === Step 1: 提取轨迹观测位置（从去噪点云中提质心）===
    [obs_positions, valid_frames, cube_info] = extract_observation_from_pointcloud(denoisedPointCloudList, frame_index_list);
    
    % === Step 2: Kalman Filter 平滑轨迹 ===
    tracked_positions = apply_kalman_filter(obs_positions);
    
    % === Step 3: 可视化轨迹 ===
    animate_kalman_tracking(denoisedPointCloudList, valid_frames, tracked_positions, cube_info);

    % 保存可视化对比动画
    save_pointcloud_denoised_video(pointCloudList, denoisedPointCloudList, frame_index_list);

end

if Is_plot == 1
    % 动画方式逐帧播放点云
    load('pointCloudList.mat');
    play_and_save_pointcloud(pointCloudList, frame_index_list)
end


% read .bin radar raw data collected from mmWave Studio CLI Tool
folder_location = './Test/PostProc/';

rawFrameDataFileName = 'radarCube';

data_frames = readDCA1000(folder_location); % (numFrames,numChirps,numRxChan,numSample);

disp('Starting Radar Signal Processing...');

numFrames = size(data_frames,1); % 1820
numChirps = size(data_frames,2); % 288
numRxChan = size(data_frames,3); % 4
numSample = size(data_frames,4); % 128

% parameter setting
params = get_params_value();
% constant parameters
c = params.c; % Speed of light in air (m/s)
fc = params.fc; % Center frequency (Hz)
lambda = params.lambda;
Rx = params.Rx;
Tx = params.Tx;

% configuration parameters
Fs = params.Fs;
sweepSlope = params.sweepSlope;
samples = params.samples;
loop = params.loop;

fft_Rang = params.fft_Rang;
fft_Vel = params.fft_Vel;
fft_Ang_Az = params.fft_Ang_Az; % 方位角 FFT 点数
fft_Ang_El = params.fft_Ang_El; % 俯仰角 FFT 点数
num_crop = params.num_crop;
max_value = params.max_value; % normalization the maximum of data WITH 1843

% Creat grid table
rng_grid = params.rng_grid;
agl_grid_az = params.agl_grid_az; % 方位角网格
agl_grid_el = params.agl_grid_el; % 俯仰角网格
vel_grid = params.vel_grid;

%% Compute Steering Vectors
% Configuration (FOV and step in degrees)
phi_FOV_deg   = 60;    % azimuth field of view (degrees, half-angle)
theta_FOV_deg = 60;    % elevation field of view (degrees, half-angle)
phi_step_deg  = 5;     % azimuth angle step (degrees)
theta_step_deg = 5;    % elevation angle step (degrees)

% Antenna geometry for IWR6843ISK-ODS (virtual antenna indices and phase rotations)
% These arrays define the horizontal (m) and vertical (n) index for each of the 12 virtual antennas,
% and the 180-degree phase rotation coefficients for azimuth (phase_rot = ±1 for each antenna).
m_ind = [0,  0, -1, -1,  -2, -2, -3, -3,  -2, -2, -3, -3].';   % azimuth indices (column vector)
n_ind = [0, -1, -1,  0,   0, -1, -1,  0,  -2, -3, -3, -2].';   % elevation indices (column vector)
phase_rot = [-1, 1, 1, -1,  -1, 1, 1, -1,  -1, 1, 1, -1].';    % phase rotation for each antenna (column vector)

Nr = length(m_ind); % Number of virtual antennas (12)

% Calculate nu_grid (azimuth direction cosine grid)
NA = round(2*phi_FOV_deg/phi_step_deg);
nu_init = -sin(deg2rad(phi_FOV_deg));
nu_step = (2 * abs(nu_init)) / (NA-1);
nu_grid = nu_init : nu_step : -nu_init;

% Calculate mu_grid (elevation direction cosine grid)
NE = round(2*theta_FOV_deg/theta_step_deg);
mu_init = -sin(deg2rad(theta_FOV_deg));
mu_step = (2 * abs(mu_init)) / (NE-1);
mu_grid = mu_init : mu_step : -mu_init;

% Compute Azimuth Steering Vectors (Including phase rotation)
azimuth_steering_vectors = zeros(Nr, NA);
for i = 1:NA
    nu = nu_grid(i);
    azimuth_steering_vectors(:, i) = phase_rot.' .* exp(1j * pi * m_ind.' * nu);
end

% Compute Elevation Steering Vectors
elevation_steering_vectors = zeros(Nr, NE);
for i = 1:NE
    mu = mu_grid(i);
    elevation_steering_vectors(:, i) = exp(1j * pi * n_ind.' * mu);
end

% Results:
% nu_grid and mu_grid: 1 x N_A and 1 x N_E arrays (in radians, as sin of angles).
% azimuth_steering_vectors: N_ant x N_A complex matrix.
% elevation_steering_vectors: N_ant x N_E complex matrix.

%% Compute Azimuth-Elevation Steering Vectors

% Compute Coarse Azimuth-Elevation Steering Vectors
steeringVec_coarse = zeros(Nr, NE, NA);
for az_i = 1:NA
    for el_i = 1:NE
        nu = nu_grid(az_i);
        mu = mu_grid(el_i);
        steeringVec_coarse(:, el_i, az_i) = phase_rot.' .* exp(1j*pi*(m_ind.'*nu + n_ind.'*mu));
    end
end

if Is_plot
    plot_mu_nu_grid(mu_grid, nu_grid);
end

% 在计算steeringVec_coarse后，加入如下代码：
params.agl_grid_az = asind(nu_grid);
params.agl_grid_el = asind(mu_grid);

if Is_plot == 0
    beam_pattern_animation(steeringVec_coarse, params.agl_grid_az, params.agl_grid_el);
end

%% Algorithm parameters
frame_start = 1;
frame_end = 3000;
Is_Windowed = 1;% 1==> Windowing before doing range and angle fft

radarCube = [];

%% MATLAB implementation for MVDR-based range-angle heatmap generation

% 用于记录目标位置轨迹
target_traj = []; % 每行: [range, az, el]
frame_index_list = [];

pointCloudList = {};  % 初始化点云帧列表

for frame_idx = frame_start:frame_end
   
    DataFrame = squeeze(data_frames(frame_idx,:,:,:));
    DataFrame = permute(DataFrame, [3,2,1]); % [samples, Rx, chirp]

    % Separate TDM-MIMO data by Tx antenna
    DataFrame_Tx0 = DataFrame(:,:,1:3:end);
    DataFrame_Tx1 = DataFrame(:,:,2:3:end);
    DataFrame_Tx2 = DataFrame(:,:,3:3:end);

    % --- Add Range-Doppler FFT and CFAR Target Detection ---
    
    % Number of chirps per Tx antenna
    numChirpsTx = size(DataFrame_Tx0,3);
    numRx = size(DataFrame,2);
    
    %% Static Clutter Removal per antenna (Rx) for each Tx
    
    % Initialize mean storage
    clutter_mean = zeros(size(DataFrame,1), numRx, 3);
    
    % For Tx0 (per Rx)
    for rxIdx = 1:numRx
        meanTx0_rx = mean(DataFrame_Tx0(:,rxIdx,:), 3); % Mean over chirps for each Rx
        DataFrame_Tx0(:,rxIdx,:) = DataFrame_Tx0(:,rxIdx,:) - repmat(meanTx0_rx,[1,1,numChirpsTx]);
        clutter_mean(:,rxIdx,1) = meanTx0_rx;
    end
    
    % For Tx1 (per Rx)
    for rxIdx = 1:numRx
        meanTx1_rx = mean(DataFrame_Tx1(:,rxIdx,:), 3);
        DataFrame_Tx1(:,rxIdx,:) = DataFrame_Tx1(:,rxIdx,:) - repmat(meanTx1_rx,[1,1,numChirpsTx]);
        clutter_mean(:,rxIdx,2) = meanTx1_rx;
    end
    
    % For Tx2 (per Rx)
    for rxIdx = 1:numRx
        meanTx2_rx = mean(DataFrame_Tx2(:,rxIdx,:), 3);
        DataFrame_Tx2(:,rxIdx,:) = DataFrame_Tx2(:,rxIdx,:) - repmat(meanTx2_rx,[1,1,numChirpsTx]);
        clutter_mean(:,rxIdx,3) = meanTx2_rx;
    end

    % Range FFT
    [Rangedata_Tx0] = fft_range(DataFrame_Tx0, fft_Rang, Is_Windowed);
    [Rangedata_Tx1] = fft_range(DataFrame_Tx1, fft_Rang, Is_Windowed);
    [Rangedata_Tx2] = fft_range(DataFrame_Tx2, fft_Rang, Is_Windowed);

    % Doppler FFT
    Dopplerdata_Tx0 = fft_doppler(Rangedata_Tx0, fft_Vel, 0);
    Dopplerdata_Tx1 = fft_doppler(Rangedata_Tx1, fft_Vel, 0);
    Dopplerdata_Tx2 = fft_doppler(Rangedata_Tx2, fft_Vel, 0);

    % Dopdata_sum = squeeze(mean(abs(Dopplerdata_Tx0), 2));

    DataFrame_AllTx = cat(2, Dopplerdata_Tx0, Dopplerdata_Tx1, Dopplerdata_Tx2);  % [range, allRx, doppler]
    Dopdata_sum = squeeze(mean(abs(DataFrame_AllTx), 2));

    % Plot range-Doppler heatmap after clutter removal
    if Is_plot
        plot_rangeDop(Dopdata_sum, rng_grid, vel_grid);
    end

    % CFAR detector on Range-Velocity to detect targets 
    % Output format: [doppler index, range index(start from index 1), ...
    % cell power]
    Pfa = 1e-1;  % 设置虚警概率
    Resl_indx = cfar_RV(Dopdata_sum, fft_Rang, 0, Pfa);

    if Is_plot
        plot_rangeDop(Dopdata_sum, rng_grid, vel_grid, Resl_indx);
    end

    %% --- MVDR per CFAR-detected (range, doppler) 点 ---
    current_frame_pointcloud = [];
    
    if ~isempty(Resl_indx) && size(Resl_indx, 2) > 0
        for idx = 1:size(Resl_indx, 2)
            d = Resl_indx(1, idx);  % Doppler bin index
            r = Resl_indx(2, idx);  % Range bin index
    
            % Snapshot: antenna × chirps
            snapshot = squeeze(DataFrame_AllTx(r, :, d)).';  % [allRx x 1]
    
            % MVDR: 单个向量 → 外积作为协方差矩阵
            Ryy = snapshot * snapshot';
            Ryy = Ryy + 1e-5 * trace(Ryy)/size(Ryy,1) * eye(size(Ryy,1));  % diagonal loading
            Ryy_inv = inv(Ryy);
    
            % 波束形成方向谱 (Azimuth x Elevation)
            NAz = length(params.agl_grid_az);
            NEl = length(params.agl_grid_el);
            spectrum = zeros(NEl, NAz);
            for el_idx = 1:NEl
                for az_idx = 1:NAz
                    a = steeringVec_coarse(:, el_idx, az_idx);
                    spectrum(NEl - el_idx + 1, az_idx) = 1 / abs(a' * Ryy_inv * a);
                end
            end
    
            % 取峰值点
            [peak_val, peak_idx] = max(spectrum(:));
            [el_idx, az_idx] = ind2sub(size(spectrum), peak_idx);
    
            % 映射为物理值
            az_angle = params.agl_grid_az(az_idx);
            el_angle = params.agl_grid_el(el_idx);
            rng_val  = params.rng_grid(r);
            vel_val  = params.vel_grid(d);
    
            % 可选绘图
            if Is_plot
                figure;
                surf(params.agl_grid_az, params.agl_grid_el, 10*log10(spectrum));
                shading interp;
                hold on;
                plot3(az_angle, el_angle, 10*log10(peak_val), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
                title(sprintf('MVDR @ Range %.2fm, Velocity %.2fm/s', rng_val, vel_val));
                xlabel('Azimuth Angle (°)'); ylabel('Elevation Angle (°)'); zlabel('Power (dB)');
                colorbar;
                hold off;
            end
    
            % 存入点云: [range, velocity, azimuth, elevation, power]
            current_frame_pointcloud = [current_frame_pointcloud;
                                        rng_val, vel_val, az_angle, el_angle, peak_val];
        end
    else
        fprintf('Frame %d: No CFAR detection -> Skip MVDR\n', frame_idx);
    end
    
    % 存储每帧点云
    pointCloudList{end+1} = current_frame_pointcloud;
    frame_index_list(end+1) = frame_idx;

    fprintf('Frame %d.\n', frame_idx);

end

% %% === 保存target_traj到 .mat 文件 ===
% % save('target_traj_data.mat', 'target_traj', 'frame_index_list');
% 
% load('target_traj_data.mat');
% 
% % 动画方式逐帧播放目标轨迹
% smooth_and_plot_trajectory(target_traj, frame_index_list);


%% === 保存所有帧的点云到 .mat 文件 ===
if Is_plot == 0
    save('pointCloudList.mat', 'pointCloudList', 'frame_index_list');
    play_and_save_pointcloud(pointCloudList, frame_index_list)
end

% Save params and data to mat file
save (rawFrameDataFileName,'radarCube', '-v7.3');

%% Function
%% Function
%% Function

function params = get_params_value()
% Radar parameter settings for IWR6843ISK-ODS radar with MVDR beamforming

% Constants
params.c = physconst('LightSpeed'); % Speed of light (m/s)

% Given parameters
params.samples = 128;        % Number of ADC samples
params.Fs = 4e6;             % Sampling rate (Hz)

% ADC sampling time (单位：秒)
params.T_ADCsampling = params.samples / params.Fs;
fprintf('T_ADCsampling = %.9f s\n', params.T_ADCsampling);

params.sweepSlope = 54.713e12;  % Sweep slope (Hz/s)

% 计算 Bandwidth（Hz）
params.BW = params.sweepSlope * params.T_ADCsampling;
fprintf('Bandwidth = %.2f MHz\n', params.BW / 1e6);

params.f0 = 60.75e9;  % Start frequency in Hz
params.T_ADCstart = 25.00e-6;

% Center frequency
params.fc = params.f0 + params.T_ADCstart * params.sweepSlope + params.BW / 2;

fprintf('Center Frequency = %.3f GHz\n', params.fc / 1e9);

params.lambda = params.c / params.fc;  % Update wavelength

% % 最大距离计算
% R_max = (params.c * params.T_ADCsampling) / (2);
% fprintf('最大可测距离 R_max = %.2f m\n', R_max);

% 距离分辨率
params.range_resolution = params.c / (2 * params.BW);
fprintf('距离分辨率 ΔR = %.4f m\n', params.range_resolution);

% Virtual antenna numbers
params.Rx = 4;
params.Tx = 3;

% 各时间参数（单位秒）
params.T_idle = 30e-6;   % 30 us
params.T_ramp = 59.1e-6; % 59.1 us

% === 计算 Chirp Repetition Time ===
params.T_r = params.Tx * (params.T_idle + params.T_ramp);  % 单位：秒

% 打印结果
fprintf("Chirp repetition time Tr = %.2f us\n", params.T_r * 1e6);

params.T_Excess = params.T_ramp - params.T_ADCstart - params.T_ADCsampling;
fprintf("Excess time = %.2f us\n", params.T_Excess * 1e6);

params.loop = 96;

params.fft_Rang = 128;
params.fft_Vel = 96;

% Angle FFT dimensions for Capon beamforming
params.fft_Ang_Az = 64; % azimuth bins
params.fft_Ang_El = 64; % elevation bins
params.num_crop = 3;
params.max_value = 1e4;

% Create range grid
freq_res = params.Fs / params.fft_Rang;
freq_grid = (0:params.fft_Rang-1).' * freq_res;
params.rng_grid = freq_grid * params.c / params.sweepSlope / 2;

% Create angle grids
az_angles = linspace(-60, 60, params.fft_Ang_Az);
el_angles = linspace(-60, 60, params.fft_Ang_El);
params.agl_grid_az = az_angles; % azimuth angles
params.agl_grid_el = el_angles; % elevation angles

% % Doppler grid
% dop_grid = fftshiftfreqgrid(params.fft_Vel, 1/params.T_r);
% params.vel_grid = dop_grid * params.lambda / 2;

% Construct correct Doppler velocity grid using theory
v_max = params.lambda / (4 * params.T_r);
fprintf("v_max = %.2f m/s\n", v_max);
v_res = params.lambda / (2 * params.loop * params.T_r);
fprintf("v_res = %.2f m/s\n", v_res);
params.vel_grid = (-v_max) : v_res : (v_max - v_res);

end

% function vel_grid = construct_doppler_grid(lambda, Tc, N_fft)
%     % lambda: wavelength (m)
%     % Tc: chirp repetition interval (s)
%     % N_fft: Doppler FFT size
% 
%     % Max unambiguous velocity
%     v_max = lambda / (4 * Tc);
% 
%     % Velocity resolution
%     v_res = lambda / (2 * N_fft * Tc);
% 
%     % Doppler bins (symmetric around 0)
%     vel_grid = (-v_max) : v_res : (v_max - v_res);
% end


function plot_mu_nu_grid(mu_grid, nu_grid)
% plot_mu_nu_grid - 绘制 Steering Vector Grid 的 mu-nu 散点图
%
% 输入:
%   mu_grid : Elevation 方向余弦网格 (vector)
%   nu_grid : Azimuth 方向余弦网格 (vector)
%
% 作用:
%   展示 mu-nu 域中的方向余弦分布（用于 steering vector 可视化）

    [MU, NU] = meshgrid(mu_grid, nu_grid);
    MU = MU(:); NU = NU(:);

    figure;
    scatter(NU, MU, 30, 'b', 'filled');
    % 设置坐标轴刻度字体
    ax = gca;
    set(ax, 'FontSize', 14, 'FontWeight', 'normal');  % 坐标轴刻度字体大小 + 粗细

    xlabel('\nu (nu) (Azimuth \phi )', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('\mu (mu) (Elevation \theta )', 'FontSize', 16, 'FontWeight', 'bold');
    title('Steering Vector in \mu-\nu Domain', 'FontSize', 18, 'FontWeight', 'bold');

    % 设置合理范围
    xlim([-1 1]);
    ylim([-1 1]);

    axis equal;
    grid on;

end

function play_and_save_pointcloud(pointCloudList, frame_index_list)
% PLAY_AND_SAVE_POINTCLOUD
% 可视化并保存 3D Point Cloud 动画（按 velocity 上色）

    % 视频设置
    v = VideoWriter('pointcloud.mp4', 'MPEG-4');
    v.FrameRate = 10;
    v.Quality = 100;
    open(v);

    % 雷达安装角度 & 高度参数
    theta_tilt = deg2rad(15);    % 向下俯视角度
    radar_height = 1.80;         % 雷达高度 (m)

    % 旋转矩阵（绕 X 轴俯视）
    R_x = [1 0 0;
           0 cos(theta_tilt) sin(theta_tilt);
           0 -sin(theta_tilt) cos(theta_tilt)];

    % 创建画布
    fig = figure('Position', [100, 100, 800, 600]);
    set(fig, 'Color', 'w');
    ax = gca;
    set(ax, 'FontSize', 14);
    hold on; grid on; box on;
    view(-88, 10);
    xlabel('X (m) - Azimuth');
    ylabel('Y (m) - Depth');
    zlabel('Z (m) - Height');
    title('Radar Point Cloud Trajectory (Velocity Colored)');
    axis tight manual;
    xlim([-2, 2]); ylim([0.5, 5.5]); zlim([0, 3]);
    % colormap(jet); 
    % colorbar
    % clim([-8 8]);

    % 初始化点云图
    h_points = scatter3(NaN, NaN, NaN, 40, 'filled');

    % 输出目录（自动创建）
    out_dir = './pointcloud_csv/';
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    % === 遍历每帧点云 ===
    for k = 1:length(pointCloudList)
        pc = pointCloudList{k};  % 每帧点云 [range, az, el, velocity, power]

        if isempty(pc)
            continue;
        end

        % 拆分维度
        r   = pc(:,1);
        vel = pc(:,2);
        az  = deg2rad(pc(:,3));
        el  = deg2rad(pc(:,4));
        
        % power = pc(:,5);  % 可选：保留后续做大小映射

        % 极坐标 → 雷达坐标系
        x = r .* cos(el) .* sin(az);
        y = r .* cos(el) .* cos(az);
        z = r .* sin(el);

        % 旋转 + 高度补偿
        pos_cart = R_x * [x.'; y.'; z.'];
        pos_cart(3,:) = pos_cart(3,:) + radar_height;

        % 更新绘图
        h_points.XData = pos_cart(1,:);
        h_points.YData = pos_cart(2,:);
        h_points.ZData = pos_cart(3,:);
        h_points.CData = vel;  % 按 velocity 上色
        h_points.SizeData = 40;

        % 标题
        title(sprintf('Radar Point Cloud - Frame %d', frame_index_list(k)), 'FontSize', 16);

        drawnow;
        frame = getframe(fig);
        writeVideo(v, frame);

        % --- 保存当前帧的 CSV ---
        csv_filename = fullfile(out_dir, sprintf('frame_%04d.csv', frame_index_list(k)));
        csv_data = pc; % 已是 [range, velocity, azimuth, elevation, power]
        csv_header = {'range', 'velocity', 'azimuth', 'elevation', 'power'};
        writecell([csv_header; num2cell(csv_data)], csv_filename);
    end

    close(v);
    disp('视频保存成功：pointcloud.mp4');
    disp(['所有帧 CSV 已导出至文件夹：', out_dir]);
end

function beam_pattern_animation(steeringVec_coarse, agl_grid_az, agl_grid_el)
% 生成波束图扫描动画

[NAz, NEl] = size(steeringVec_coarse(1,:,:));
Nr = size(steeringVec_coarse,1);

azimuth_angles = agl_grid_az;
elevation_angles = agl_grid_el;

% === Step 4: 视频保存设置（高清 1080p） ===
v = VideoWriter('beam_pattern_animation.mp4', 'MPEG-4');
v.FrameRate = 20;  % 增加帧率可选
v.Quality = 100;   % 设置最高画质
open(v);

fig = figure;
set(fig, 'Color', 'w');  % 背景白

% 创建动画
for el_idx = 1:length(elevation_angles)
    for az_idx = 1:length(azimuth_angles)

        % 当前方向波束
        a_steer = steeringVec_coarse(:, el_idx, az_idx);

        % 计算波束图（空间谱）
        beam_pattern = zeros(length(elevation_angles), length(azimuth_angles));
        for el_scan_idx = 1:length(elevation_angles)
            for az_scan_idx = 1:length(azimuth_angles)
                a_scan = steeringVec_coarse(:, el_scan_idx, az_scan_idx);
                beam_pattern(el_scan_idx, az_scan_idx) = abs(a_steer' * a_scan);
            end
        end

        % 绘制波束图
        surf(azimuth_angles, elevation_angles, beam_pattern, 'EdgeColor','none');
        shading interp;
        colormap(jet);
        colorbar;
        xlabel('Azimuth (deg)'); ylabel('Elevation (deg)');
        title(sprintf('Beam Steering at Azimuth: %.1f°, Elevation: %.1f°', ...
            azimuth_angles(az_idx), elevation_angles(el_idx)));

        % 设置坐标轴刻度字体
        ax = gca;
        set(ax, 'FontSize', 14, 'FontWeight', 'normal');  % 坐标轴刻度字体大小 + 粗细

        % 减小 margin 的关键设置
        set(ax, 'LooseInset', [0 0 0 0]);  % 去掉自动边距
        axis tight manual;

        caxis([0 Nr]);
        xlim([-60 60]);
        ylim([-60 60]);
        view(0,90);
        drawnow;
        frame = getframe(fig);   % 获取全高清帧
        writeVideo(v, frame);

        pause(0.01); % 控制动画速度
    end
end

close(v);
disp('高清视频已保存为 beam_pattern_animation.mp4');

end

function clusteredPointCloudList = run_dbscan_on_pointcloud(pointCloudList, eps, minPts)
% 对每一帧点云进行 DBSCAN 聚类，并将聚类标签添加到点云数据中
%
% 输入:
%   pointCloudList - 原始点云帧列表，cell数组，每帧为 [range, vel, az, el, power]
%   eps            - DBSCAN 聚类的邻域半径参数
%   minPts         - DBSCAN 的最小邻居点数
%
% 输出:
%   clusteredPointCloudList - 聚类后的点云，每帧为 [range, vel, az, el, power, cluster_label]

    clusteredPointCloudList = cell(size(pointCloudList));

    for k = 1:length(pointCloudList)
        pc = pointCloudList{k};
        if isempty(pc)
            clusteredPointCloudList{k} = [];
            continue;
        end

        % 提取 az, el, range 转换为 x, y, z 坐标
        r = pc(:,1);
        az = deg2rad(pc(:,3));
        el = deg2rad(pc(:,4));
        x = r .* cos(el) .* sin(az);
        y = r .* cos(el) .* cos(az);
        z = r .* sin(el);

        pos = [x, y, z];

        % 使用 DBSCAN 聚类
        labels = dbscan(pos, eps, minPts);

        % 将标签添加回原点云
        pc_clustered = [pc, labels];
        clusteredPointCloudList{k} = pc_clustered;
    end

end

function plot_pointcloud_with_clusters(originalList, clusteredList, frame_index_list)
% 对比绘制原始点云和聚类结果（左右两个窗口）
%
% 输入:
%   originalList  - 原始点云
%   clusteredList - 聚类后的点云，包含 cluster label
%   frame_index_list - 帧索引列表

    % 雷达旋转角度设置
    theta_tilt = deg2rad(15);  % 俯视角度
    radar_height = 1.8;

    R_x = [1 0 0;
           0 cos(theta_tilt) sin(theta_tilt);
           0 -sin(theta_tilt) cos(theta_tilt)];

    fig = figure('Position', [100, 100, 1200, 600]);
    for k = 1:length(originalList)
        clf;
        
        % === 原始点云绘图 ===
        subplot(1,2,1);
        pc = originalList{k};
        if isempty(pc)
            continue;
        end

        r = pc(:,1);
        az = deg2rad(pc(:,3));
        el = deg2rad(pc(:,4));
        vel = pc(:,2);

        x = r .* cos(el) .* sin(az);
        y = r .* cos(el) .* cos(az);
        z = r .* sin(el);

        pos = R_x * [x.'; y.'; z.'];
        pos(3,:) = pos(3,:) + radar_height;

        scatter3(pos(1,:), pos(2,:), pos(3,:), 40, vel, 'filled');
        title(sprintf('Original - Frame %d', frame_index_list(k)));
        xlabel('X (Azimuth)'); ylabel('Y (Depth)'); zlabel('Z (Height)');
        axis equal; grid on; view(-88, 10);
        xlim([-2 2]); ylim([0.5 5.5]); zlim([0 3]);

        % === 聚类结果绘图 ===
        subplot(1,2,2);
        pc_cl = clusteredList{k};
        if isempty(pc_cl)
            continue;
        end

        r = pc_cl(:,1);
        az = deg2rad(pc_cl(:,3));
        el = deg2rad(pc_cl(:,4));
        labels = pc_cl(:,6);

        x = r .* cos(el) .* sin(az);
        y = r .* cos(el) .* cos(az);
        z = r .* sin(el);

        pos = R_x * [x.'; y.'; z.'];
        pos(3,:) = pos(3,:) + radar_height;

        scatter3(pos(1,:), pos(2,:), pos(3,:), 40, labels, 'filled');
        title(sprintf('DBSCAN Clusters - Frame %d', frame_index_list(k)));
        xlabel('X (Azimuth)'); ylabel('Y (Depth)'); zlabel('Z (Height)');
        axis equal; grid on; view(-88, 10);
        xlim([-2 2]); ylim([0.5 5.5]); zlim([0 3]);

        drawnow;
        pause(0.05);
    end

end


function denoisedPointCloudList = run_dbscan_denoise(pointCloudList, eps, minPts)
% 使用 DBSCAN 去除离群点，返回每帧保留真实目标的点云
%
% 输入:
%   pointCloudList - 原始点云，cell，每帧为 [range, vel, az, el, power]
%   eps            - 邻域半径
%   minPts         - 最小聚类点数
%
% 输出:
%   denoisedPointCloudList - 去除离群点后的点云列表

    denoisedPointCloudList = cell(size(pointCloudList));

    for k = 1:length(pointCloudList)
        pc = pointCloudList{k};
        if isempty(pc)
            denoisedPointCloudList{k} = [];
            continue;
        end

        % 极坐标转换为 3D 坐标
        r = pc(:,1);
        az = deg2rad(pc(:,3));
        el = deg2rad(pc(:,4));

        x = r .* cos(el) .* sin(az);
        y = r .* cos(el) .* cos(az);
        z = r .* sin(el);
        pos = [x, y, z];

        % DBSCAN 聚类
        labels = dbscan(pos, eps, minPts);

        % 只保留非 -1 标签的点（即非噪声点）
        keep_idx = labels ~= -1;
        pc_clean = pc(keep_idx, :);

        denoisedPointCloudList{k} = pc_clean;
    end
end

function plot_pointcloud_denoised(originalList, denoisedList, frame_index_list)
% 左右对比绘制：原始 vs 去除outliers后的点云

    theta_tilt = deg2rad(15);
    radar_height = 1.8;
    R_x = [1 0 0;
           0 cos(theta_tilt) sin(theta_tilt);
           0 -sin(theta_tilt) cos(theta_tilt)];

    fig = figure('Position', [100, 100, 1600, 600]);

    for k = 1:length(originalList)
        clf;

        %% 左侧：原始点云
        subplot(1,2,1);
        pc = originalList{k};
        if isempty(pc)
            continue;
        end

        r = pc(:,1);
        az = deg2rad(pc(:,3));
        el = deg2rad(pc(:,4));
        vel = pc(:,2);

        x = r .* cos(el) .* sin(az);
        y = r .* cos(el) .* cos(az);
        z = r .* sin(el);
        pos = R_x * [x.'; y.'; z.'];
        pos(3,:) = pos(3,:) + radar_height;

        scatter3(pos(1,:), pos(2,:), pos(3,:), 40, vel, 'filled');
        title(sprintf('Original - Frame %d', frame_index_list(k)));
        xlabel('X (Azimuth)'); ylabel('Y (Depth)'); zlabel('Z (Height)');
        axis equal; grid on; view(-88, 10);
        xlim([-2 2]); ylim([0.5 5.5]); zlim([0 3]);

        %% 右侧：去除离群点后的点云
        subplot(1,2,2);
        pc_cl = denoisedList{k};
        if isempty(pc_cl)
            title('No valid clusters');
            continue;
        end

        r = pc_cl(:,1);
        az = deg2rad(pc_cl(:,3));
        el = deg2rad(pc_cl(:,4));
        vel = pc_cl(:,2);

        x = r .* cos(el) .* sin(az);
        y = r .* cos(el) .* cos(az);
        z = r .* sin(el);
        pos = R_x * [x.'; y.'; z.'];
        pos(3,:) = pos(3,:) + radar_height;

        scatter3(pos(1,:), pos(2,:), pos(3,:), 40, vel, 'filled');
        title(sprintf('DBSCAN Filtered - Frame %d', frame_index_list(k)));
        xlabel('X (Azimuth)'); ylabel('Y (Depth)'); zlabel('Z (Height)');
        axis equal; grid on; view(-88, 10);
        xlim([-2 2]); ylim([0.5 5.5]); zlim([0 3]);

        drawnow;
        pause(0.05);
    end
end

function save_pointcloud_denoised_video(originalList, denoisedList, frame_index_list)
% 保存左右对比动画视频：原始 vs DBSCAN 去噪（紧凑高清版）

    v = VideoWriter('pointcloud_denoised_comparison.mp4', 'MPEG-4');
    v.Quality = 100;
    v.FrameRate = 10;
    open(v);

    theta_tilt = deg2rad(15);
    radar_height = 1.8;
    R_x = [1 0 0;
           0 cos(theta_tilt) sin(theta_tilt);
           0 -sin(theta_tilt) cos(theta_tilt)];

    % 创建 figure 和 tiled layout
    fig = figure('Position', [-1920, 10, 1920*0.8, 1080-216-216], 'Color', 'w');
    t = tiledlayout(1,2, 'Padding', 'tight', 'TileSpacing', 'tight');

    for k = 1:length(originalList)
        nexttile(1);
        cla;
        pc = originalList{k};
        if ~isempty(pc)
            r = pc(:,1);
            az = deg2rad(pc(:,3));
            el = deg2rad(pc(:,4));
            vel = pc(:,2);
            x = r .* cos(el) .* sin(az);
            y = r .* cos(el) .* cos(az);
            z = r .* sin(el);
            pos = R_x * [x.'; y.'; z.'];
            pos(3,:) = pos(3,:) + radar_height;

            scatter3(pos(1,:), pos(2,:), pos(3,:), 40, vel, 'filled');
            title(sprintf('Original - Frame %d', frame_index_list(k)), 'FontSize', 20, 'FontWeight', 'bold');
            xlabel('X (Azimuth) [m]', 'FontSize', 16, 'FontWeight', 'normal');
            ylabel('Y (Depth) [m]', 'FontSize', 16, 'FontWeight', 'normal');
            zlabel('Z (Elevation) [m]', 'FontSize', 16, 'FontWeight', 'normal');
            axis equal; grid on; view(-88, 10);
            xlim([-2 2]); ylim([0.5 5.5]); zlim([0 3]);
            set(gca, 'FontSize', 18);
        else
            title('No original data', 'FontSize', 18, 'FontWeight', 'bold');
        end

        nexttile(2);
        cla;
        pc_cl = denoisedList{k};
        if ~isempty(pc_cl)
            r = pc_cl(:,1);
            az = deg2rad(pc_cl(:,3));
            el = deg2rad(pc_cl(:,4));
            vel = pc_cl(:,2);
            x = r .* cos(el) .* sin(az);
            y = r .* cos(el) .* cos(az);
            z = r .* sin(el);
            pos = R_x * [x.'; y.'; z.'];
            pos(3,:) = pos(3,:) + radar_height;

            scatter3(pos(1,:), pos(2,:), pos(3,:), 40, vel, 'filled');
            title(sprintf('DBSCAN Filtered - Frame %d', frame_index_list(k)), 'FontSize', 20, 'FontWeight', 'bold');
            xlabel('X (Azimuth) [m]', 'FontSize', 16, 'FontWeight', 'normal');
            ylabel('Y (Depth) [m]', 'FontSize', 16, 'FontWeight', 'normal');
            zlabel('Z (Elevation) [m]', 'FontSize', 16, 'FontWeight', 'normal');
            axis equal; grid on; view(-88, 10);
            xlim([-2 2]); ylim([0.5 5.5]); zlim([0 3]);
            set(gca, 'FontSize', 18);
        else
            title('No valid clusters', 'FontSize', 18, 'FontWeight', 'bold');
        end

        drawnow;
        frame = getframe(fig);
        writeVideo(v, frame);
    end

    close(v);
    disp('视频保存完成：pointcloud_denoised_comparison.mp4');
end


function [obs_positions, valid_frames, cube_info] = extract_observation_from_pointcloud(denoisedPointCloudList, frame_index_list)
% 使用滑动立方体提取每帧观测值（即使无目标也输出 NaN）
% 输出:
%   obs_positions: N x 3 每行是 [x y z]（可含 NaN）
%   valid_frames: N x 1 对应帧号，始终与 pointCloudList 对应

    theta_tilt = deg2rad(15);
    radar_height = 1.8;
    R_x = [1 0 0;
           0 cos(theta_tilt) sin(theta_tilt);
           0 -sin(theta_tilt) cos(theta_tilt)];

    cube_size = [0.75, 0.75, 1];  % [X Y Z] in meters
    step = 0.2;

    N = length(denoisedPointCloudList);
    obs_positions = nan(N, 3);          % 全部初始化为 NaN
    valid_frames = frame_index_list(:); % 保持维度一致
    cube_info = nan(N, 6);  % [center_x, center_y, center_z, size_x, size_y, size_z]

    for k = 1:N
        pc = denoisedPointCloudList{k};
        if isempty(pc)
            continue;  % 留 NaN
        end

        % 极坐标转笛卡尔
        r = pc(:,1);
        az = deg2rad(pc(:,3));
        el = deg2rad(pc(:,4));

        x = r .* cos(el) .* sin(az);
        y = r .* cos(el) .* cos(az);
        z = r .* sin(el);

        pos = R_x * [x.'; y.'; z.'];
        pos(3,:) = pos(3,:) + radar_height;
        pos = pos.';  % Nx3

        % 滑动 cube 寻找密集点区域
        xmin = -2; xmax = 2;
        ymin = 0.5; ymax = 5.5;
        zmin = 0; zmax = 3;

        max_count = 0;
        best_center = [NaN, NaN, NaN];
        best_box = [NaN NaN NaN];

        for xi = xmin:step:(xmax - cube_size(1))
            for yi = ymin:step:(ymax - cube_size(2))
                for zi = zmin:step:(zmax - cube_size(3))
                    x1 = xi; x2 = xi + cube_size(1);
                    y1 = yi; y2 = yi + cube_size(2);
                    z1 = zi; z2 = zi + cube_size(3);

                    in_cube = ...
                        pos(:,1) >= x1 & pos(:,1) <= x2 & ...
                        pos(:,2) >= y1 & pos(:,2) <= y2 & ...
                        pos(:,3) >= z1 & pos(:,3) <= z2;

                    count = sum(in_cube);

                    if count > max_count
                        max_count = count;
                        best_center = [(x1+x2)/2, (y1+y2)/2, (z1+z2)/2];
                        best_box = [cube_size];
                    end
                end
            end
        end

        obs_positions(k, :) = best_center;
        cube_info(k, :) = [best_center, best_box];
    end
end


function tracked_positions = apply_kalman_filter(obs_positions)
% 支持部分缺失观测（NaN）的 Kalman Filter 实现
% 输入:
%   obs_positions: N x 3，观测点 (可包含 NaN 行)
% 输出:
%   tracked_positions: N x 3，滤波后的目标轨迹

    N = size(obs_positions, 1);
    tracked_positions = nan(N, 3);  % 初始化输出

    dt = 1;  % 时间间隔

    % 状态：[x y z vx vy vz]
    A = [1 0 0 dt 0 0;
         0 1 0 0 dt 0;
         0 0 1 0 0 dt;
         0 0 0 1 0 0;
         0 0 0 0 1 0;
         0 0 0 0 0 1];

    H = [1 0 0 0 0 0;
         0 1 0 0 0 0;
         0 0 1 0 0 0];

    Q = 1e-3 * eye(6);  % 过程噪声
    R = 5e-2  * eye(3); % 观测噪声
    P = eye(6);         % 协方差初值

    % === 初始化状态 ===
    first_valid = find(all(~isnan(obs_positions), 2), 1, 'first');
    if isempty(first_valid)
        warning('No valid observations found.');
        return;
    end
    x = [obs_positions(first_valid,:) 0 0 0]';  % 初始状态

    for i = 1:N
        % === Prediction ===
        x_pred = A * x;
        P_pred = A * P * A' + Q;

        z = obs_positions(i,:)';

        if all(~isnan(z))  % 有观测
            % === Update ===
            K = P_pred * H' / (H * P_pred * H' + R);
            x = x_pred + K * (z - H * x_pred);
            P = (eye(6) - K * H) * P_pred;
        else  % 缺失观测：只预测
            x = x_pred;
            P = P_pred;
        end

        tracked_positions(i,:) = x(1:3)';
    end
end



function animate_kalman_tracking(denoisedList, frame_index_list, tracked_positions, cube_info)
% 显示每帧 Kalman 点（跳过 NaN）+ 红色目标框 + 背景点云，保存为高清动画

    v = VideoWriter('kalman_tracking_animation.mp4', 'MPEG-4');
    v.FrameRate = 10;
    v.Quality = 100;
    open(v);

    theta_tilt = deg2rad(15);
    radar_height = 1.8;
    R_x = [1 0 0;
           0 cos(theta_tilt) sin(theta_tilt);
           0 -sin(theta_tilt) cos(theta_tilt)];

    x_range = [-2, 2];
    y_range = [0.5, 5.5];
    z_range = [0, 3];

    fig = figure('Position', [100, 100, 1920/2, 1080/2], 'Color', 'w');
    ax = axes('Parent', fig);
    hold(ax, 'on');
    set(ax, 'FontSize', 14);
    view(ax, -88, 10);
    xlim(ax, x_range);
    ylim(ax, y_range);
    zlim(ax, z_range);
    grid(ax, 'on');
    box(ax, 'on');
    axis(ax, 'manual');

    xlabel(ax, 'X (Azimuth) [m]', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel(ax, 'Y (Depth) [m]', 'FontSize', 16, 'FontWeight', 'bold');
    zlabel(ax, 'Z (Elevation) [m]', 'FontSize', 16, 'FontWeight', 'bold');

    hXLabel = get(gca,'XLabel');
    set(hXLabel,'rotation',90,'VerticalAlignment','middle')

    for k = 1:length(tracked_positions)
        cla(ax);

        % 背景点云
        pc = denoisedList{k};
        if ~isempty(pc)
            r = pc(:,1);
            az = deg2rad(pc(:,3));
            el = deg2rad(pc(:,4));
            x = r .* cos(el) .* sin(az);
            y = r .* cos(el) .* cos(az);
            z = r .* sin(el);
            pos = R_x * [x.'; y.'; z.'];
            pos(3,:) = pos(3,:) + radar_height;
            scatter3(ax, pos(1,:), pos(2,:), pos(3,:), 10, 'k', 'filled');
        end

        % 绘制红框（目标提取cube）
        box_data = cube_info(k, :);
        if all(~isnan(box_data))
            draw_red_box(ax, box_data(1:3), box_data(4:6));
        end

        % Kalman 当前帧点（若非 NaN）
        p = tracked_positions(k, :);
        if all(~isnan(p))
            scatter3(ax, p(1), p(2), p(3), 100, 'b', 'filled');
        end

        title(ax, sprintf('Kalman Tracking - Frame %d', frame_index_list(k)), ...
              'FontSize', 20, 'FontWeight', 'bold');

        drawnow;
        frame = getframe(fig);
        writeVideo(v, frame);
    end

    close(v);
    disp('Kalman 动画已保存：kalman_tracking_animation.mp4');
end

function draw_red_box(ax, center, size_vec)
% 绘制一个红色立方体框
% center: [x y z]  中心点
% size_vec: [dx dy dz]  宽高深

    x = center(1);
    y = center(2);
    z = center(3);
    dx = size_vec(1)/2;
    dy = size_vec(2)/2;
    dz = size_vec(3)/2;

    % 8个顶点
    corners = [
        x-dx y-dy z-dz;
        x+dx y-dy z-dz;
        x+dx y+dy z-dz;
        x-dx y+dy z-dz;
        x-dx y-dy z+dz;
        x+dx y-dy z+dz;
        x+dx y+dy z+dz;
        x-dx y+dy z+dz;
    ];

    edges = [1 2;2 3;3 4;4 1; % 底面
             5 6;6 7;7 8;8 5; % 顶面
             1 5;2 6;3 7;4 8]; % 垂直边

    for i = 1:size(edges,1)
        pt1 = corners(edges(i,1),:);
        pt2 = corners(edges(i,2),:);
        plot3(ax, [pt1(1), pt2(1)], [pt1(2), pt2(2)], [pt1(3), pt2(3)], ...
              'r-', 'LineWidth', 1.5);
    end
end
