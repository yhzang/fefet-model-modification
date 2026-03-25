%% =========================================================
% extract_Vth_disturb.m
%
% 用途：
%   从扰动仿真结果中提取每个选定累计脉冲数对应的 Vth / DeltaVth
%
% 运行方式：
%   先在工作区中跑完主仿真，至少要有：
%       Stsum, index, epife
%
%   最好同时有：
%       index_pre, index_end, time, volt, Nread
%
% 说明：
%   - 本脚本默认配合 wfdef_acc_fix.m 使用
%   - 若 index_pre 存在，则以“第 1 个 pulse 之前”的状态作为 Vth0
%   - 若 index_pre 不存在，则退化为“第一个有效 Vth 点”为参考
%
% 读出默认参数：
%   这里沿用原仓库注释读出块的思路：
%       VG       = -0.5:0.02:1.7
%       VD       = 0.05
%       Iref     = 1e-7 A
%       VFB_read = -0.5
% =========================================================

%% -------------------------
% 用户可改参数
% -------------------------
trace_id = 1;              % 若 Stsum 是矩阵，取第 trace_id 列
VG = -0.5:0.02:3;        % 读出栅压扫描
VD = 0.05;
VS = 0;
Iref = 1e-7;               % 定义阈值电流
VFB_read = -0.5;           % 沿用原仓库注释中的读出方式

% MOS / stack 参数（与原注释读出块一致）
tfe = 0.8e-6;              % cm
til = 1e-7;                % cm
miu = 50;                  % cm^2 / V / s
Na  = 3e17;                % cm^-3
T   = 300;                 % K
W   = 1;
L   = 1;

show_waveform_debug = true;   % 是否画波形 + 索引标记

%% -------------------------
% 基本检查
% -------------------------
if ~exist('Stsum', 'var')
    error('extract_Vth_disturb:MissingStsum', ...
        '工作区里没有 Stsum，请先跑主仿真。');
end

if ~exist('index', 'var')
    error('extract_Vth_disturb:MissingIndex', ...
        '工作区里没有 index，请先用 wfdef_acc_fix.m 生成索引。');
end

if ~exist('epife', 'var')
    error('extract_Vth_disturb:MissingEpife', ...
        '工作区里没有 epife。');
end

index = index(:);

if exist('index_pre', 'var') && ~isempty(index_pre)
    index_pre = index_pre(:);
end

if exist('index_end', 'var') && ~isempty(index_end)
    index_end = index_end(:);
end

%% -------------------------
% 选取 Stsum 轨迹
% -------------------------
if isvector(Stsum)
    st_trace = Stsum(:);
else
    if trace_id > size(Stsum, 2)
        error('extract_Vth_disturb:BadTraceId', ...
            'trace_id 超出 Stsum 的列数。');
    end
    st_trace = Stsum(:, trace_id);
end

Npulse = numel(index);

%% -------------------------
% 默认读取点
% 若用户没给 Nread，就只取 decade 点 + 最后一点
% -------------------------
if ~exist('Nread', 'var') || isempty(Nread)
    if Npulse <= 0
        error('extract_Vth_disturb:BadPulseCount', 'index 为空。');
    end
    Nread = unique([1; (10.^(0:floor(log10(Npulse)))).'; Npulse]);
end

Nread = unique(round(Nread(:)));
Nread = Nread(Nread >= 1 & Nread <= Npulse);

if isempty(Nread)
    error('extract_Vth_disturb:EmptyNread', 'Nread 过滤后为空。');
end

%% -------------------------
% 先取参考点 Vth0
% 优先使用“第 1 个 pulse 之前”的状态
% -------------------------
if exist('index_pre', 'var') && ~isempty(index_pre)
    idx0 = index_pre(1);
elseif index(1) > 1
    idx0 = index(1) - 1;
else
    idx0 = NaN;
end

valid0 = false;
Vth0 = NaN;
P0   = NaN;

if isfinite(idx0) && idx0 >= 1 && idx0 <= numel(st_trace)
    P0 = st_trace(idx0);
    [Vth0, ID0, ID0min, ID0max, valid0] = local_extract_vth( ...
        P0, epife, tfe, til, miu, Na, T, W, L, VG, VFB_read, VD, VS, Iref);
else
    ID0 = [];
    ID0min = NaN;
    ID0max = NaN;
end

%% -------------------------
% 提取每个 Nread 的 P / Vth
% -------------------------
Np = numel(Nread);

idx_read = NaN(Np,1);
Pread    = NaN(Np,1);
Vth      = NaN(Np,1);
DeltaVth = NaN(Np,1);
IDmin    = NaN(Np,1);
IDmax    = NaN(Np,1);
is_valid = false(Np,1);

for kk = 1:Np
    n = Nread(kk);
    idx_read(kk) = index(n);

    if idx_read(kk) < 1 || idx_read(kk) > numel(st_trace)
        continue;
    end

    Pread(kk) = st_trace(idx_read(kk));

    [Vth(kk), IDtmp, IDmin(kk), IDmax(kk), is_valid(kk)] = local_extract_vth( ...
        Pread(kk), epife, tfe, til, miu, Na, T, W, L, VG, VFB_read, VD, VS, Iref);
end

%% -------------------------
% 参考点处理
% -------------------------
if valid0
    DeltaVth = Vth - Vth0;
    Vth_ref = Vth0;
    ref_mode = 'pre-pulse-1';
else
    first_valid = find(is_valid, 1, 'first');
    if isempty(first_valid)
        error('extract_Vth_disturb:NoValidVth', ...
            '所有读点都没有跨过 Iref，当前 VG / Iref 设置下无法提取 Vth。');
    end
    Vth_ref = Vth(first_valid);
    DeltaVth = Vth - Vth_ref;
    ref_mode = sprintf('first-valid-read-point (N=%d)', Nread(first_valid));
    fprintf(2, ['[warning] 没有拿到 pre-pulse 基准点，已退化为第一个有效读点做参考。\n' ...
                '          参考脉冲数 N = %d\n'], Nread(first_valid));
end

%% -------------------------
% 文本输出
% -------------------------
fprintf('\n================ disturb Vth extraction ================\n');
fprintf('trace_id   = %d\n', trace_id);
fprintf('Iref       = %.3e A\n', Iref);
fprintf('VFB_read   = %.3f V\n', VFB_read);
fprintf('reference  = %s\n', ref_mode);

if valid0
    fprintf('Vth0       = %.6f V   (idx0 = %d, P0 = %.6f)\n', Vth0, idx0, P0);
else
    fprintf('Vth0       = NaN\n');
end

fprintf('\n');
fprintf('%10s %10s %16s %14s %14s %14s %14s %8s\n', ...
    'Nread', 'idx', 'Pread(C/cm^2)', 'Vth(V)', 'DeltaVth(V)', 'IDmin(A)', 'IDmax(A)', 'valid');

for kk = 1:Np
    fprintf('%10d %10d %16.6e %14.6f %14.6f %14.6e %14.6e %8d\n', ...
        Nread(kk), idx_read(kk), Pread(kk), Vth(kk), DeltaVth(kk), ...
        IDmin(kk), IDmax(kk), is_valid(kk));
end

%% -------------------------
% 图 1：波形与索引
% -------------------------
if show_waveform_debug && exist('time', 'var') && exist('volt', 'var')
    figure;
    plot(time, volt, 'b-'); hold on;

    if exist('index_pre', 'var') && ~isempty(index_pre)
        plot(time(index_pre(Nread)), volt(index_pre(Nread)), ...
            'ks', 'MarkerSize', 7, 'LineWidth', 1.2);
    end

    if exist('index_end', 'var') && ~isempty(index_end)
        plot(time(index_end(Nread)), volt(index_end(Nread)), ...
            'md', 'MarkerSize', 7, 'LineWidth', 1.2);
    end

    plot(time(index(Nread)), volt(index(Nread)), ...
        'ro', 'MarkerSize', 7, 'LineWidth', 1.2);

    xlabel('time (s)');
    ylabel('V_G (V)');
    title('Pulse waveform with index\_pre / index\_end / index(read)');
    grid on;
    legend_items = {};
    if exist('index_pre', 'var') && ~isempty(index_pre)
        legend_items{end+1} = 'waveform';
        legend_items{end+1} = 'index\_pre';
        if exist('index_end', 'var') && ~isempty(index_end)
            legend_items{end+1} = 'index\_end';
            legend_items{end+1} = 'index(read)';
            legend('waveform', 'index\_pre', 'index\_end', 'index(read)', 'Location', 'best');
        else
            legend('waveform', 'index\_pre', 'index(read)', 'Location', 'best');
        end
    else
        if exist('index_end', 'var') && ~isempty(index_end)
            legend('waveform', 'index\_end', 'index(read)', 'Location', 'best');
        else
            legend('waveform', 'index(read)', 'Location', 'best');
        end
    end
    hold off;
end

%% -------------------------
% 图 2：DeltaVth vs N
% -------------------------
figure;
semilogx(Nread, DeltaVth, 'o-', 'LineWidth', 1.5, 'MarkerSize', 7);
xlabel('Cumulative pulse count N');
ylabel('\Delta V_{th} (V)', 'Interpreter', 'tex');
title('Disturb-induced threshold-voltage shift');
grid on;

%% -------------------------
% 图 3：等效极化读点
% -------------------------
figure;
semilogx(Nread, Pread, 's-', 'LineWidth', 1.5, 'MarkerSize', 7);
xlabel('Cumulative pulse count N');
ylabel('P_{read} (C/cm^2)', 'Interpreter', 'tex');
title('Read-point ferroelectric state');
grid on;

%% -------------------------
% 保存结果到工作区
% -------------------------
disturb_result = struct();
disturb_result.trace_id = trace_id;
disturb_result.Nread = Nread;
disturb_result.idx_read = idx_read;
disturb_result.idx0 = idx0;
disturb_result.P0 = P0;
disturb_result.Pread = Pread;
disturb_result.Vth0 = Vth0;
disturb_result.Vth = Vth;
disturb_result.DeltaVth = DeltaVth;
disturb_result.IDmin = IDmin;
disturb_result.IDmax = IDmax;
disturb_result.is_valid = is_valid;
disturb_result.Iref = Iref;
disturb_result.VFB_read = VFB_read;
disturb_result.VG = VG;
disturb_result.ref_mode = ref_mode;

assignin('base', 'disturb_result', disturb_result);

%% =========================================================
% Local function
% =========================================================
function [Vth, ID, idmin, idmax, is_valid] = local_extract_vth( ...
    Psum, epi_fe, tfe, til, miu, Na, T, W, L, VG, VFB, VD, VS, Iref)

    ID = get_ID(Psum, epi_fe, tfe, til, miu, Na, T, W, L, VG, VFB, VD, VS);
    ID = real(ID(:)).';

    valid_id = isfinite(ID) & (ID > 0);

    if any(valid_id)
        idmin = min(ID(valid_id));
        idmax = max(ID(valid_id));
    else
        idmin = NaN;
        idmax = NaN;
    end

    % 不做 extrap，只找真正跨过 Iref 的那一段
    cross_idx = find( ...
        isfinite(ID(1:end-1)) & isfinite(ID(2:end)) & ...
        (ID(1:end-1) > 0) & (ID(2:end) > 0) & ...
        ((ID(1:end-1) - Iref) .* (ID(2:end) - Iref) <= 0), ...
        1, 'first');

    if isempty(cross_idx)
        Vth = NaN;
        is_valid = false;
        return;
    end

    x = log10([ID(cross_idx), ID(cross_idx+1)]);
    y = [VG(cross_idx), VG(cross_idx+1)];

    if abs(x(2) - x(1)) < 1e-30
        Vth = mean(y);
    else
        Vth = interp1(x, y, log10(Iref), 'linear');
    end

    is_valid = isfinite(Vth);
end