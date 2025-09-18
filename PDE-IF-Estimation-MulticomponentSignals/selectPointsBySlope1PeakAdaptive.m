function result =selectPointsBySlope1PeakAdaptive(varargin)


args = struct('P', [], 'F', [], 'PeakIdx', [], 'WindowRatio', 1.0);
args = parseArgs(args, varargin{:});

P = args.P(:);
F = args.F(:);
N = length(P);

% 主峰位置
if isempty(args.PeakIdx)
    [~, peak_idx] = max(P);
else
    peak_idx = args.PeakIdx;
end

% ---- 1. 计算半高宽（FWHM） ----
peak_val = P(peak_idx);
half_val = peak_val / 2;

% 向左找半高点
left_idx = find(P(1:peak_idx) < half_val, 1, 'last');
if isempty(left_idx)
    left_idx = 1;
end

% 向右找半高点
right_idx = peak_idx - 1 + find(P(peak_idx:end) < half_val, 1, 'first');
if isempty(right_idx)
    right_idx = N;
end

% 半高宽
fwhm = abs(F(right_idx) - F(left_idx));
% 若F非等间隔，下面这样更通用
% window_width = round(args.WindowRatio * (right_idx - left_idx) / 2);

% 用半高宽确定自适应窗口左右范围
Lwin = round(args.WindowRatio * (peak_idx - left_idx));
Rwin = round(args.WindowRatio * (right_idx - peak_idx));

Lwin = max(Lwin, 1); % 至少1
Rwin = max(Rwin, 1);

% 限制不要越界
idx_left = max(peak_idx - Lwin, 1):(peak_idx-1);
idx_right = (peak_idx+1):min(peak_idx + Rwin, N);

% ---- 2. 斜率变缓点选择 ----
dL = abs(diff(P(idx_left)));
[~, iL] = min(dL);
idx1 = idx_left(iL);

dR = abs(diff(P(idx_right)));
[~, iR] = min(dR);
idx2 = idx_right(iR);

% ---- 3. 结果输出 ----
result.xi1_1 = F(idx1);
result.xi1_2 = F(idx2);
result.idx1 = idx1;
result.idx2 = idx2;
result.peak_idx = peak_idx;
result.fwhm = fwhm; % 可选，输出半高宽
end

function args = parseArgs(defaults, varargin)
for i = 1:2:length(varargin)
    key = varargin{i};
    val = varargin{i+1};
    if isfield(defaults, key)
        defaults.(key) = val;
    else
        warning('未知参数名: %s', key);
    end
end
args = defaults;
end
