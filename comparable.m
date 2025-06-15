clc
clear
close all

%% === 加载 MVDR 点云数据 ===
load('pointCloudList.mat');  % 包含 pointCloudList, frame_index_list


% 输入起止时间字符串（含微秒）
start_str = '20250410142659140154';
end_str   = '20250410142950336485';

% 转换为 datetime（支持微秒）
start_time = datetime(start_str, 'InputFormat', 'yyyyMMddHHmmssSSSSSS');
end_time   = datetime(end_str,   'InputFormat', 'yyyyMMddHHmmssSSSSSS');

% 获取帧数量
% num_frames = length(frame_index_list);
num_frames = 3077;

% 生成等间隔时间戳
mvdr_timestamps = linspace(start_time, end_time, num_frames).';

% 设置显示格式统一
mvdr_timestamps.Format = 'yyyy-MM-dd HH:mm:ss.SSS';

% 计算实际帧间隔（毫秒）
actual_intervals = milliseconds(diff(mvdr_timestamps));
mean_interval = mean(actual_intervals);
std_interval  = std(actual_intervals);

% 打印结果
fprintf('[MVDR] Total frames: %d\n', num_frames);
fprintf('[MVDR] Actual average frame interval: %.3f ms\n', mean_interval);
% fprintf('[MVDR] Std of frame interval: %.3f ms\n', std_interval);

% % MVDR 点云帧数据
% mvdr_frames = pointCloudList;

% DBSCAN 参数设置
eps = 0.5;    % 距离阈值（单位：m）
minPts = 3;   % 最小聚类点数

% 运行 DBSCAN 去噪
mvdr_frames = run_dbscan_denoise(pointCloudList, eps, minPts);


%%
radar_dir = './Test/radar_F/';

loop_times = 1;

[frames, ts_list, ts_mean_datetime, x_radar_fr, y_radar_fr, z_radar_fr, v_radar_fr, snr, point_ids, target_ids, ...
 posX, posY, posZ, velX, velY, velZ] = load_radarF_csv_frames(radar_dir, loop_times);

ts_mean_datetime.Format = 'yyyy-MM-dd HH:mm:ss.SSS';

%%
sync_table = match_mvdr_to_radar(mvdr_timestamps, ts_mean_datetime);
disp(sync_table(1:20,:));


%%

plot_all_frame_statistics_comparison(sync_table, mvdr_frames, x_radar_fr, y_radar_fr, z_radar_fr, v_radar_fr);


plot_mvdr_vs_radar_comparison(sync_table, mvdr_frames, frames, x_radar_fr, y_radar_fr, z_radar_fr, v_radar_fr);


%%
function sorted = natsortfiles_custom(files)
    expr = '\d+';
    nums = zeros(length(files), 1);
    for i = 1:length(files)
        tokens = regexp(files{i}, expr, 'match');
        if ~isempty(tokens)
            nums(i) = str2double(tokens{end});
        else
            nums(i) = Inf;
        end
    end
    [~, idx] = sort(nums);
    sorted = files(idx);
end

function [radar_frames, radar_ts_list, radar_ts_mean_datetime, ...
          x_list, y_list, z_list, ...
          v_list, snr_list, ...
          point_ids, target_ids, ...
          posX_list, posY_list, posZ_list, ...
          velX_list, velY_list, velZ_list] = load_radarF_csv_frames(radar_dir, loop_times)

    % === 获取文件列表 ===
    file_structs = dir(fullfile(radar_dir, '*.csv'));
    file_names = natsortfiles_custom({file_structs.name});

    df_radar_bag = table();
    
    colNames = {'timestamp', 'point_id', 'elev', 'azim', 'doppler', 'range', 'snr', ...
            'x', 'y', 'z', 'Target_Id', 'posX', 'posY', 'posZ', ...
            'velX', 'velY', 'velZ', 'accX', 'accY', 'accZ'};


    for k = 1:length(file_names)
        file_path = fullfile(radar_dir, file_names{k});
        try
            T = readtable(file_path, ...
                          'Delimiter', ',', ...
                          'FileType', 'text', ...
                          'PreserveVariableNames', true);
        catch
            warning('无法读取文件: %s', file_path);
            continue;
        end

        % T.Properties.VariableNames = matlab.lang.makeValidName(strtrim(T.Properties.VariableNames));
        T.Properties.VariableNames = colNames;

        if isempty(T)
            continue;
        end

        df_radar_bag = [df_radar_bag; T];
    end

    % === 找到 point_id 列 ===
    var_names = df_radar_bag.Properties.VariableNames;
    point_id_col = var_names(contains(lower(var_names), 'point_id'));
    if isempty(point_id_col)
        error('未找到 point_id 列');
    end

    point_ids_all = df_radar_bag.(point_id_col{1});
    if iscell(point_ids_all)
        point_ids_all = str2double(point_ids_all);
    end

    indx_fr = find(point_ids_all == 0);
    indx_fr(end+1) = height(df_radar_bag);  % 结束边界

    % === 初始化输出 ===
    radar_frames = {};
    radar_ts_list = {};
    radar_ts_mean = [];

    x_list = {}; y_list = {}; z_list = {};
    v_list = {}; snr_list = {};

    point_ids = {}; target_ids = {};

    posX_list = {}; posY_list = {}; posZ_list = {};
    velX_list = {}; velY_list = {}; velZ_list = {};

    % === 按帧滑动提取 ===
    for i = 1:(length(indx_fr) - loop_times)
        row_start = indx_fr(i);
        row_end = indx_fr(i + loop_times);

        df_i = df_radar_bag(row_start:row_end-1, :);
        radar_frames{end+1} = df_i;

        radar_ts_list{end+1} = df_i.timestamp;
        radar_ts_mean(end+1) = mean(df_i.timestamp);

        x_list{end+1} = df_i.x;
        y_list{end+1} = df_i.y;
        z_list{end+1} = df_i.z;
        v_list{end+1} = df_i.doppler;
        snr_list{end+1} = df_i.snr;

        point_ids{end+1} = df_i.(point_id_col{1});

        if ismember('Target_Id', df_i.Properties.VariableNames)
            target_ids{end+1} = df_i.Target_Id;
        else
            target_ids{end+1} = zeros(height(df_i), 1);
        end

        posX_list{end+1} = df_i.posX;
        posY_list{end+1} = df_i.posY;
        posZ_list{end+1} = df_i.posZ;

        velX_list{end+1} = df_i.velX;
        velY_list{end+1} = df_i.velY;
        velZ_list{end+1} = df_i.velZ;
    end

    % 转为 datetime
    radar_ts_mean_datetime = convert_radar_ns_timestamp(radar_ts_mean);
end

function datetimes = convert_radar_ns_timestamp(timestamps)
% 将 20 位雷达 timestamp 数字数组转换为 MATLAB datetime 向量
% 输入：timestamps — 数值数组（如 2.025041014270582e+19）
% 输出：datetimes — datetime 类型数组，格式为 yyyy-MM-dd HH:mm:ss.SSS

    % 1. 数字 → 字符串（去掉科学计数法）
    strTimestamps = arrayfun(@(x) num2str(x, '%.0f'), timestamps, 'UniformOutput', false);

    % 2. 取前 17 位用于构造 datetime（忽略最后几位纳秒）
    years   = str2double(cellfun(@(x) x(1:4),   strTimestamps, 'UniformOutput', false));
    months  = str2double(cellfun(@(x) x(5:6),   strTimestamps, 'UniformOutput', false));
    days    = str2double(cellfun(@(x) x(7:8),   strTimestamps, 'UniformOutput', false));
    hours   = str2double(cellfun(@(x) x(9:10),  strTimestamps, 'UniformOutput', false));
    minutes = str2double(cellfun(@(x) x(11:12), strTimestamps, 'UniformOutput', false));
    seconds = str2double(cellfun(@(x) x(13:14), strTimestamps, 'UniformOutput', false));
    msecs   = str2double(cellfun(@(x) x(15:17), strTimestamps, 'UniformOutput', false));

    % 3. 构造 datetime
    datetimes = datetime(years, months, days, hours, minutes, seconds, msecs, ...
                         'Format', 'yyyy-MM-dd HH:mm:ss.SSS');

end


function sync_table = match_mvdr_to_radar(mvdr_timestamps, radar_timestamps)
    % 筛掉 MVDR 比雷达早的帧
    valid_idx = mvdr_timestamps >= min(radar_timestamps);
    mvdr_timestamps = mvdr_timestamps(valid_idx);
    mvdr_frame_id = find(valid_idx);

    % 匹配
    n = length(mvdr_timestamps);
    matched_radar_idx = zeros(n,1);
    matched_radar_ts = NaT(n,1);
    time_diff_ms = NaN(n,1);

    max_diff_ms = 50;  % 设置最大容许误差

    for i = 1:n
        [min_diff, idx] = min(abs(mvdr_timestamps(i) - radar_timestamps));
        % matched_radar_idx(i) = idx;
        % matched_radar_ts(i) = radar_timestamps(idx);
        % time_diff_ms(i) = milliseconds(min_diff);
        diff_ms = milliseconds(min_diff);

        if diff_ms <= max_diff_ms
            matched_radar_idx(i) = idx;
            matched_radar_ts(i) = radar_timestamps(idx);
            time_diff_ms(i) = diff_ms;
        end

    end

    % 设置统一格式
    mvdr_timestamps.Format  = 'yyyy-MM-dd HH:mm:ss.SSS';
    matched_radar_ts.Format = 'yyyy-MM-dd HH:mm:ss.SSS';

    % 构建表格
    sync_table = table(mvdr_frame_id(:), mvdr_timestamps(:), matched_radar_ts, ...
                       matched_radar_idx, time_diff_ms, ...
                       'VariableNames', {'MVDR_Frame', 'MVDR_Timestamp', ...
                                         'Radar_Timestamp', 'Radar_Frame_Index', ...
                                         'Time_Diff_ms'});
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


function plot_mvdr_vs_radar_comparison(sync_table, mvdr_frames, radar_frames, x_radar, y_radar, z_radar, v_radar)
    % 绘制前 N 帧比较
    num_plot_frames = min(3000, height(sync_table)); % 可调

    % 创建视频文件（带时间戳命名）
    timestamp_str = datestr(now, 'yyyymmdd_HHMMSS');
    video_filename = sprintf('Radar_P.C._Generation_vs_IWR6843ISK-ODS.mp4');
    v = VideoWriter(video_filename, 'MPEG-4');
    v.FrameRate = 10;
    v.Quality = 100;
    open(v);

    % 创建图窗
    fig = figure('Position', [10 10 1200 500], 'Color', 'w');

    for i = 1:num_plot_frames
        mvdr_idx = sync_table.MVDR_Frame(i);
        radar_idx = sync_table.Radar_Frame_Index(i);

        if radar_idx == 0 || radar_idx > length(radar_frames)
            continue;  % 没有匹配的雷达帧
        end

        % 获取数据
        mvdr_pc = mvdr_frames{mvdr_idx};
        radar_pc = [x_radar{radar_idx}, ...
            y_radar{radar_idx}, ...
            z_radar{radar_idx}, ...
            v_radar{radar_idx}];  % --> N x 4 matrix


        % 画图
        clf;

        subplot(1,2,1);
        % =====================[ 处理 MVDR 点云 ]===========================
        if ~isempty(mvdr_pc)
            % MVDR 格式: [range, velocity, azimuth, elevation, power]
            r   = mvdr_pc(:,1);
            vel = mvdr_pc(:,2);
            az  = deg2rad(mvdr_pc(:,3));
            el  = deg2rad(mvdr_pc(:,4));
            
            % 极坐标 → 雷达坐标系
            x = r .* cos(el) .* sin(az);
            y = r .* cos(el) .* cos(az);
            z = r .* sin(el);
            
            % 俯视角旋转（绕 X 轴）+ 高度补偿
            theta_tilt = deg2rad(15);  % 俯视角
            radar_height = 1.80;       % 安装高度
            R_x = [1 0 0;
                   0 cos(theta_tilt) sin(theta_tilt);
                   0 -sin(theta_tilt) cos(theta_tilt)];
            
            pos_cart = R_x * [x.'; y.'; z.'];
            pos_cart(3,:) = pos_cart(3,:) + radar_height;
            
            % 更新 MVDR 点云图
            scatter3(pos_cart(1,:), pos_cart(2,:), pos_cart(3,:), ...
                     40, vel, 'filled');  % velocity 上色
            title(sprintf('Enriched Radar Point Clouds', mvdr_idx));
        else
            scatter3(NaN, NaN, NaN, '.');
            title(sprintf('Enriched Radar Point Clouds', mvdr_idx));
        end
        xlabel('X'); ylabel('Y'); zlabel('Z');
        axis equal; grid on;
        xlim([-2, 2]); ylim([0.5, 5.5]); zlim([0, 3]);
        view(-75, 15);

        subplot(1,2,2);
        % =====================[ 处理 Radar 点云 ]===========================
        if ~isempty(radar_pc)
            x_r = radar_pc(:,1);
            y_r = radar_pc(:,2);
            z_r = radar_pc(:,3);
            v_r = radar_pc(:,4);  % Doppler 速度
        
            % === 加入旋转与高度补偿（与 MVDR 一致） ===
            theta_tilt = deg2rad(15);    % 向下俯视角度
            radar_height = 1.80;         % 雷达高度 (m)
        
            R_x = [1 0 0;
                   0 cos(theta_tilt) sin(theta_tilt);
                   0 -sin(theta_tilt) cos(theta_tilt)];
        
            pos_rot = R_x * [x_r.'; y_r.'; z_r.'];
            pos_rot(3,:) = pos_rot(3,:) + radar_height;
        
            scatter3(pos_rot(1,:), pos_rot(2,:), pos_rot(3,:), ...
                     40, v_r, 'filled');
            title(sprintf('Original TI IWR6843ISK-ODS', radar_idx));
        else
            scatter3(NaN, NaN, NaN, '.');
            title(sprintf('Original TI IWR6843ISK-ODS', radar_idx));
        end
        xlabel('X'); ylabel('Y'); zlabel('Z');
        axis equal; grid on;
        xlim([-2, 2]); ylim([0.5, 5.5]); zlim([0, 3]);
        view(-75, 15);


        sgtitle(sprintf('Frame Comparison (∆t = %.2f ms)', sync_table.Time_Diff_ms(i)));

        drawnow;

        % 写入视频
        frame = getframe(fig);
        writeVideo(v, frame);

    end

    close(v);
    disp(['视频保存成功：', video_filename]);

end

function plot_all_frame_statistics_comparison(sync_table, mvdr_frames, ...
    x_radar, y_radar, z_radar, v_radar)

    n = height(sync_table);
    mvdr_count = zeros(n,1);
    radar_count = zeros(n,1);
    mvdr_speed_mean = NaN(n,1);
    radar_speed_mean = NaN(n,1);

    % 汇总用于直方图的数据
    all_mvdr_pos = [];
    all_radar_pos = [];

    for i = 1:n
        mvdr_idx = sync_table.MVDR_Frame(i);
        radar_idx = sync_table.Radar_Frame_Index(i);

        if radar_idx == 0 || radar_idx > length(x_radar)
            continue;
        end

        % ==== MVDR 点云 ====
        mvdr_pc = mvdr_frames{mvdr_idx};
        if ~isempty(mvdr_pc)
            vel = mvdr_pc(:,2);
            mvdr_count(i) = size(mvdr_pc,1);
            mvdr_speed_mean(i) = mean(abs(vel), 'omitnan');

            % 极坐标转笛卡尔坐标
            r   = mvdr_pc(:,1);
            v   = mvdr_pc(:,2);
            az  = deg2rad(mvdr_pc(:,3));
            el  = deg2rad(mvdr_pc(:,4));

            x = r .* cos(el) .* sin(az);
            y = r .* cos(el) .* cos(az);
            z = r .* sin(el);

            theta_tilt = deg2rad(15);
            radar_height = 1.80;
            R_x = [1 0 0; 0 cos(theta_tilt) sin(theta_tilt); 0 -sin(theta_tilt) cos(theta_tilt)];
            pos = R_x * [x.'; y.'; z.'];
            pos(3,:) = pos(3,:) + radar_height;

            all_mvdr_pos = [all_mvdr_pos; pos(1,:)', pos(2,:)', pos(3,:)', v];
        end

        % ==== Radar 点云 ====
        x_r = x_radar{radar_idx};
        y_r = y_radar{radar_idx};
        z_r = z_radar{radar_idx};
        v_r = v_radar{radar_idx};

        radar_count(i) = length(x_r);
        radar_speed_mean(i) = mean(abs(v_r), 'omitnan');

        if ~isempty(x_r)
            theta_tilt = deg2rad(15);
            radar_height = 1.80;
            R_x = [1 0 0; 0 cos(theta_tilt) sin(theta_tilt); 0 -sin(theta_tilt) cos(theta_tilt)];
            pos_r = R_x * [x_r.'; y_r.'; z_r.'];
            pos_r(3,:) = pos_r(3,:) + radar_height;

            all_radar_pos = [all_radar_pos; pos_r(1,:)', pos_r(2,:)', pos_r(3,:)', v_r];
        end
    end

    %% ==== 绘制直方图 ====
    figure('Color', 'w', 'Position', [100, 100, 800, 500]);
    
    % 绘制直方图
    histogram(mvdr_count, 'Normalization', 'pdf', ...
        'DisplayName', 'MVDR', 'FaceAlpha', 0.5, 'FaceColor', [0 0.4470 0.7410]);
    hold on;
    histogram(radar_count, 'Normalization', 'pdf', ...
        'DisplayName', 'Radar', 'FaceAlpha', 0.5, 'FaceColor', [0.8500 0.3250 0.0980]);
    
    % 拟合正态分布
    x = linspace(min([mvdr_count; radar_count]), max([mvdr_count; radar_count]), 500);
    mvdr_fit = normpdf(x, mean(mvdr_count), std(mvdr_count));
    radar_fit = normpdf(x, mean(radar_count), std(radar_count));
    
    plot(x, mvdr_fit, '--', 'Color', [0 0.2 0.6], 'LineWidth', 2.0, 'DisplayName', 'MVDR Fit');
    plot(x, radar_fit, '--', 'Color', [0.8 0.2 0.2], 'LineWidth', 2.0, 'DisplayName', 'Radar Fit');
    
    % 图形美化设置
    xlabel('Number of Points', 'FontSize', 14);
    ylabel('Probability Density', 'FontSize', 14);
    title('PDF of Point Count (MVDR vs Radar)', 'FontSize', 16);
    legend('Location', 'northeast', 'FontSize', 12);
    grid on;
    set(gca, 'FontSize', 13);
    box on;


    %% ==== 分布直方图绘图 ====
    figure('Name','Point Cloud Distribution Histograms','Color','w', 'Position',[1920+100, -216+100, 1200, 800]);

    fields = {'X (Azimuth)', 'Y (Depth)', 'Z (Height)', 'Velocity'};
    for i = 1:4
        subplot(2,2,i);
        h1 = histogram(all_mvdr_pos(:,i), 'BinWidth', 0.1, 'FaceAlpha', 0.6, 'DisplayName', 'MVDR');
        hold on;
        h2 = histogram(all_radar_pos(:,i), 'BinWidth', 0.1, 'FaceAlpha', 0.6, 'DisplayName', 'Radar');

        % 自动设置 X 轴范围
        all_data = [all_mvdr_pos(:,i); all_radar_pos(:,i)];
        xlim([min(all_data)/2, max(all_data)/2]);

        xlabel(fields{i}, 'FontSize', 12);
        ylabel('Count', 'FontSize', 12);
        title(['Histogram of ', fields{i}], 'FontSize', 14);
        legend('FontSize', 11);
        grid on;
    end
end
