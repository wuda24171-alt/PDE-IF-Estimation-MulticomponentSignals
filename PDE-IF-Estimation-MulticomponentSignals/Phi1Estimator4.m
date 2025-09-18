classdef Phi1Estimator4
    properties
        F
        T
        P
        P_u
        P_xi
        dt
    end

    methods
        function obj = Phi1Estimator4(F, T, P, P_u, P_xi)
            % 构造函数：初始化频率、时间和谱图数据
            obj.F = F;
            obj.T = T(:);  % 强制为列向量
            obj.P = P;
            obj.P_u = P_u;
            obj.P_xi = P_xi;
            obj.dt = obj.T(2) - obj.T(1);
        end

        function phi1_2nd = computeSecondDerivative(obj, xi1_1, xi1_2)
            xi1_1 = xi1_1(:);
            xi1_2 = xi1_2(:);
            N = min([length(obj.T), length(xi1_1), length(xi1_2), size(obj.P,2)]);
            phi1_2nd = NaN(1, N);

            for k = 1:N
                if isnan(xi1_1(k)) || isnan(xi1_2(k))
                    continue;
                end
                [~, idx1] = min(abs(obj.F - xi1_1(k)));
                [~, idx2] = min(abs(obj.F - xi1_2(k)));

                if any([idx1, idx2] < 1) || any([idx1, idx2] > size(obj.P,1))
                    continue;
                end

                P1 = obj.P(idx1, k); P2 = obj.P(idx2, k);
                Pu1 = obj.P_u(idx1, k); Pu2 = obj.P_u(idx2, k);
                Pf1 = obj.P_xi(idx1, k); Pf2 = obj.P_xi(idx2, k);

                numerator   = Pu1 * P2 - Pu2 * P1;
                denominator = Pf1 * P2 - Pf2 * P1;

                if abs(denominator) < 1e-50
                    phi1_2nd(k) = 0;
                else
                    phi1_2nd(k) = -numerator / denominator;
                end
            end
        end
function phi1_2nd_updated = smoothSections(obj, phi1_2nd, section1_times, section2_times, s)
    % Savitzky-Golay滤波参数
    windowLength = 2 * round(s / (obj.T(2)-obj.T(1))) + 1;
    if windowLength < 5, windowLength = 5; end
    polyOrder = 3;

    % 中值滤波参数
    median_win = 5;

    % TV滤波参数
    tv_lambda = 0.1;
    tv_iters = 120;
    tv_threshold = 1e-3;

    % 分段多项式阶数
    fit_order = 2; % 二阶多项式，若异常变化剧烈可用3

    T = obj.T;
    phi1_2nd = phi1_2nd(:);
    phi1_2nd_updated = phi1_2nd;

    % 需要平滑的所有区间的点
    all_times = [section1_times(:); section2_times(:)];
    idx_all = [];

    for k = 1:length(all_times)
        u = all_times(k);
        idx_start = find(T >= u - s, 1, 'first');
        idx_end   = find(T <= u + s, 1, 'last');
        if isempty(idx_start) || isempty(idx_end)
            continue;
        end
        idx_all = [idx_all, idx_start:idx_end]; %#ok<AGROW>
    end

    idx_all = unique(idx_all);  % 避免重复

    % 只对这些区间的数据做SG+中值+TV三重滤波
    data_seg = phi1_2nd(idx_all);
    % 1. SG滤波
    data_sg = sgolayfilt(data_seg, polyOrder, min(windowLength, length(data_seg)));
    % 2. 中值滤波
    data_median = medfilt1(data_sg, median_win);
    % 3. TV滤波
    data_tv = tvdenoise(data_median, tv_lambda, tv_iters, tv_threshold);

    % 4. 多项式拟合拉平
    t_seg = T(idx_all);
    p = polyfit(t_seg, data_tv, fit_order);
    data_polyfit = polyval(p, t_seg);

    % 用多项式拟合后的结果替换该区间
    phi1_2nd_updated(idx_all) = data_polyfit;
end



function CR_integrated = tvSmoothAndIntegrate(obj, phi1_2nd)
    % TV滤波参数设置
    lambda = 0.005;   % 平滑强度（可根据实际调整，越小越平滑，越大越保真）
    iters = 300;     % 最大迭代步数
    threshold = 1e-4; % 收敛阈值

    % 先插值到更高分辨率
    T_ext = linspace(min(obj.T), max(obj.T), length(obj.T)*10);
    phi_interp = interp1(obj.T, phi1_2nd, T_ext, 'linear', 'extrap');

    % TV滤波
    phi_tv = tvdenoise(phi_interp(:), lambda, iters, threshold); % 一定要用(:)保证为列向量

    % 数值积分
    dt_ext = T_ext(2) - T_ext(1);
    CR_integrated = cumsum(phi_tv) * dt_ext;
end



function CR_integrated = sgSmoothAndIntegrate(obj, phi1_2nd)
    % SG滤波参数设置
    windowLength = 21; % 必须为奇数，建议结合数据间隔和噪声量调整
    polyOrder = 3;     % 多项式阶数，建议2或3

    % 先插值
    T_ext = linspace(min(obj.T), max(obj.T), length(obj.T)*10);
    phi_interp = interp1(obj.T, phi1_2nd, T_ext, 'linear', 'extrap');
    
    % Savitzky-Golay 滤波平滑
    % （注意窗口不能大于数据长度，必要时可加判别处理）
    if windowLength > length(phi_interp)
        windowLength = length(phi_interp) - mod(length(phi_interp)+1,2); % 保证为奇数
    end
    phi_sg = sgolayfilt(phi_interp, polyOrder, windowLength);

    % 数值积分
    dt_ext = T_ext(2) - T_ext(1);
    %CR_integrated = cumsum(phi_sg) * dt_ext;
    CR_integrated = cumtrapz(T_ext, phi_sg);
end
function CR_integrated = movmeanSmoothAndIntegrate(obj, phi1_2nd)
    % 参数设置
    upsample_rate = 10;           % 插值倍数
    movmean_window = 15;          % 移动平均窗口，建议为奇数，根据分辨率和数据平滑程度调整

    % 1. 插值
    T_ext = linspace(min(obj.T), max(obj.T), length(obj.T) * upsample_rate);
    phi_interp = interp1(obj.T, phi1_2nd, T_ext, 'linear', 'extrap');

    % 2. 移动平均滤波
    phi_ma = movmean(phi_interp, movmean_window);

    % 3. 梯形积分
    CR_integrated = cumtrapz(T_ext, phi_ma);
end



 function plotAll(obj, T, phi1_2nd, phi1_2nd_updated, CR_integrated)
            figure;
            subplot(2,1,1);
            plot(T, phi1_2nd, 'k--'); hold on;
            plot(T, phi1_2nd_updated, 'b');hold on;
            legend('原始', '替换后');
            title('\phi_1''''(u) 替换前后'); xlabel('时间'); grid on;

            subplot(2,1,2);
            T_ext = linspace(min(T), max(T), length(CR_integrated));
            plot(T_ext, CR_integrated, 'LineWidth', 1.5);
            title('\int \phi_1''''(u) du'); xlabel('时间'); ylabel('积分'); grid on;
        end
    end
end
