function result = selectPointsBySlope1Peak(varargin)

% 参数解析
args = struct('P', [], 'F', [], 'SearchWidth', 4, 'PeakIdx', []);
args = parseArgs(args, varargin{:});

P = args.P(:);
F = args.F(:);
N = length(P);

if isempty(args.PeakIdx)
    [~, peak_idx] = max(P);
else
    peak_idx = args.PeakIdx;
end

% 限制边界
W = args.SearchWidth;
peak_idx = min(max(peak_idx, W+1), N - W);

% 左边查找斜率变缓点
left = (peak_idx-W):(peak_idx-1);
dL = abs(diff(P(left)));
[~, iL] = min(dL);
idx1 = left(iL);

% 右边查找斜率变缓点
right = (peak_idx+1):(peak_idx+W);
dR = abs(diff(P(right)));
[~, iR] = min(dR);
idx2 = right(iR);

result.xi1_1 = F(idx1);
result.xi1_2 = F(idx2);
result.idx1 = idx1;
result.idx2 = idx2;

end

function args = parseArgs(defaults, varargin)
% 参数解析工具
args = defaults;
for i = 1:2:length(varargin)
    key = varargin{i};
    val = varargin{i+1};
    if isfield(args, key)
        args.(key) = val;
    else
        warning('未知参数名: %s', key);
    end
end
end
