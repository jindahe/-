clc;
clear ;
close all;

%% 初始化
tic;
%% 目标图片二值化
image_path = 'zju.jpg';
original_image = imread(image_path);
if size(original_image, 3) == 3
    grayscale_image = rgb2gray(original_image); % 将 RGB 转为灰度
else
    grayscale_image = original_image; % 已是灰度图
end

threshold = 128;
binary_image = imbinarize(grayscale_image, threshold / 255);

%% 调整分辨率
target_resolution = [200, 200];
resized_binary_image = imresize(binary_image, target_resolution);

% 二值化处理
binary_threshold = 0.5; % 自定义阈值
resized_binary_image = resized_binary_image > binary_threshold; % 手动二值化

%{
subplot(1,3,1);
imshow(resized_binary_image);
title('调整分辨率后的二值化图像');
%}
%% 正式过程
target_pattern = resized_binary_image; % 定义目标图案
%% 初始光阑

% 基于目标图案加随机扰动
perturbation_rate = 0.5; % 设置扰动率
random_mask = rand(size(target_pattern)) < perturbation_rate;
initial_aperture = xor(imbinarize(double(target_pattern), 0.5), random_mask);


%smoothed_pattern = imgaussfilt(double(target_pattern), 2); % 高斯模糊，sigma=2
%initial_aperture = imbinarize(smoothed_pattern, 0.5);

%{
subplot(1,3,2);
imshow(initial_aperture);
title('初始光阑');
%}

disp(['初始化完成，运行时间: ', num2str(toc), ' 秒']);
%% 逆向设计光阑
tic;
temp_aperture = monto_carlo_algorithm(target_pattern, initial_aperture);
disp(['monto完成，运行时间: ', num2str(toc), ' 秒']);

tic;
final_aperture = greedy_algorithm(target_pattern, temp_aperture);
disp(['greedy完成，运行时间: ', num2str(toc), ' 秒']);

%subplot(1,3,2);
%imshow(final_aperture);
%title('最终光阑');
%% 模拟结果
tic;
final_pattern = simulate_diffraction(final_aperture);
disp(['衍射模拟完成，运行时间: ', num2str(toc), ' 秒']);

%{
subplot(1,3,3);
imshow(final_pattern);
title('衍射图案');
%}

%% 衍射模拟
function diffraction_pattern = simulate_diffraction(aperture)
    % 参数设置
    lambda = 1; % 波长 (m)
    z = 353; % 传输距离 (m)
    k = 2 * pi / lambda; % 波数
  
    [num_rows, num_cols] = size(aperture);
    dx = 1e-3; % 光阑采样间距 (m)
    dy = dx; % 假设 x 和 y 间距相等

    % 构建频率坐标系
    fx = (-num_cols / 2 : num_cols / 2 - 1) / (num_cols * dx);
    fy = (-num_rows / 2 : num_rows / 2 - 1) / (num_rows * dy);
    [Fx, Fy] = meshgrid(fx, fy);

    % 傅里叶变换求解衍射场
    aperture_fft = fftshift(fft2(ifftshift(aperture))); % 光阑的傅里叶变换
    H = exp(1i * k * z) * exp(-1i * pi * lambda * z * (Fx.^2 + Fy.^2)); % 传输函数
    diffraction_field = aperture_fft .* H; % 传输到观察平面的频域场
    field = fftshift(ifft2(ifftshift(diffraction_field))); % 逆变换得到空间域场

    % 计算光强分布并归一化
    intensity = abs(field).^2;
    diffraction_pattern = intensity / max(intensity(:)); % 归一化强度

    % 二值化处理
    binary_threshold = 0.5; % 自定义阈值
    diffraction_pattern = diffraction_pattern > binary_threshold; % 手动二值化
    
end
%% 蒙特卡洛算法
function optimized_aperture = monto_carlo_algorithm(target_pattern, initial_aperture)

    [rows, cols] = size(target_pattern);
    aperture = initial_aperture;
    T = 1;
    cooling_rate = 0.9;
    iterations_per_temp = 10;

    % 初始化误差
    previous_pattern = simulate_diffraction(aperture);
    previous_error = calculate_error(target_pattern, previous_pattern);
    error_history = previous_error; % 记录误差历史
        
    while T > 0.01
        % 初始化结果存储
        for iter = 1:iterations_per_temp
            % 随机选择并计算初始误差
            r = randi(rows);
            c = randi(cols);
          
            % 改变状态
            aperture(r, c) = 1 -aperture(r, c);

            % 模拟衍射
            current_pattern = simulate_diffraction(aperture);

            % 计算误差
            current_error = calculate_error(target_pattern, current_pattern);

            % 比较
            delta_error = current_error - previous_error;
            if delta_error > 0 && rand() > exp(-delta_error/T)
                aperture(r,c) = 1 -aperture(r, c);
            else
                % 更新误差
                previous_error = current_error;
            end
        end
        T = T * cooling_rate;
         % 更新误差历史
        error_history = [error_history, previous_error]; % 保存当前误差
    end

    optimized_aperture = aperture;

    % 绘制误差曲线
   
    subplot(1,2,1);
    plot(error_history, 'LineWidth', 2);
    xlabel('迭代次数');
    ylabel('误差');
    title('monto误差曲线');
    
    
end
%% 贪心算法
function optimized_aperture = greedy_algorithm(target_pattern, initial_aperture)
    % 输入:
    % target_pattern - 目标衍射图案 (二维矩阵)
    % initial_aperture - 初始光阑状态 (二维矩阵)
    % 输出:
    % optimized_aperture - 经过优化后的光阑状态

    [rows, cols] = size(target_pattern);
    aperture = initial_aperture; % 当前光阑状态

     % 初始化误差
    previous_pattern = simulate_diffraction(aperture);
    previous_error = calculate_error(target_pattern, previous_pattern);
    error_history = previous_error; % 记录误差历史

    for iter = 1:3 % 迭代次数，可根据需要调整
        for r = 1:rows
            for c = 1:cols
                % 改变像素状态（透光或不透光）
                aperture(r, c) = 1 - aperture(r, c); % 状态翻转
                
                % 模拟当前衍射图案
                current_pattern = simulate_diffraction(aperture);
                
                % 计算误差
                current_error = calculate_error(target_pattern, current_pattern);
                
                % 如果误差降低，接受更改；否则回退
                if current_error > previous_error
                    aperture(r, c) = 1 - aperture(r, c); % 恢复状态
                else
                    % 更新误差历史
                    previous_error = current_error;
                end
            end
        end
      
        error_history = [error_history, current_error]; % 保存当前误差
    end

    optimized_aperture = aperture;

    
    
    subplot(1,2,2);
    plot(error_history, 'LineWidth', 2);
    xlabel('迭代次数');
    ylabel('误差');
    title('greedy误差曲线');
    
    
end
%% 误差计算
function error = calculate_error(expected_pattern, actual_pattern)
    
    % 误差计算
    error = sum((expected_pattern(:) - actual_pattern(:)).^2);
end

