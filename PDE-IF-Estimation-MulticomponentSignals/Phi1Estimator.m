classdef Phi1Estimator
    properties
        F
        T
        P
        P_u
        P_xi
        dt
    end

    methods
        function obj = Phi1Estimator(F, T, P, P_u, P_xi)
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
            T = obj.T;
             phi1_2nd = phi1_2nd(:);
             phi1_2nd_updated = phi1_2nd;
            
            all_times = [section1_times(:); section2_times(:)];
            for k = 1:length(all_times)
                u = all_times(k);
                idx_start = find(T >= u - s, 1, 'first');
                idx_end   = find(T <= u + s, 1, 'last');

                if isempty(idx_start) || isempty(idx_end)
                    continue;
                end

                idx_range = idx_start:idx_end;
                phi_window = phi1_2nd(idx_range);
                local_avg = (1 / (2 * s)) * sum(phi_window, 'omitnan') * obj.dt;
                phi1_2nd_updated(idx_range) = local_avg;
            end
        end

        function CR_integrated = gaussianSmoothAndIntegrate(obj, phi1_2nd)
            epsilon = 0.05;
            sigma = epsilon / 2;
            u_gauss = linspace(-10, 10, 1000);
            rho = exp(-u_gauss.^2 / (2 * sigma^2));
            rho = rho / sum(rho);

            T_ext = linspace(min(obj.T), max(obj.T), length(obj.T)*10);
            phi_interp = interp1(obj.T, phi1_2nd, T_ext, 'linear', 'extrap');

            CR_normalized = conv(phi_interp, rho, 'same');
            dt_ext = T_ext(2) - T_ext(1);
            CR_integrated = cumsum(CR_normalized) * dt_ext;
        end

        function plotAll(obj, T, phi1_2nd, phi1_2nd_updated, CR_integrated)
            figure;
            subplot(2,1,1);
            plot(T, phi1_2nd, 'k--'); hold on;
            plot(T, phi1_2nd_updated, 'b');
            legend('原始', '替换后');
            title('\phi_1''''(u) 替换前后'); xlabel('时间'); grid on;

            subplot(2,1,2);
            T_ext = linspace(min(T), max(T), length(CR_integrated));
            plot(T_ext, CR_integrated, 'LineWidth', 1.5);
            title('\int \phi_1''''(u) du'); xlabel('时间'); ylabel('积分'); grid on;
        end
    end
end
