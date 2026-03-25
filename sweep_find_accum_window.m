%% =========================================================
% sweep_find_accum_window.m
%
% 目的：
%   扫描脉冲幅值 amp，寻找“非饱和累积窗口”：
%   即 DeltaVth(N) 不是全 0，也不是 N=1 后立刻饱和。
%
% 依赖：
%   wfdef_acc_fix.m
%   get_FE_state.m
%   get_ID.m
%
% 输出：
%   1) 每个 amp 一条 DeltaVth(N) 曲线
%   2) summary_table：每个 amp 的统计结果
%   3) 保存每个 case 的 MAT 文件
%
% 建议：
%   先用这个脚本找合适 amp；
%   找到之后，再固定 amp 去扫 delay。
%% =========================================================

clear;
clc;
close all;

rng(0);

%% =========================
% 0) 输出目录
% =========================
outdir = 'sweep_amp_out';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

%% =========================
% 1) 固定参数
% =========================
VDD     = 5;

% 这里只扫 amp
amp_list = [0.2 0.3 0.5 0.65 0.8 1.0 1.2 1.5 2.0 2.5];

pw      = 0.1e-6;
delay   = 1e-6;
transit = 0.02e-6;
cycle   = 100;          % 先看前100个脉冲是否渐进
step    = 0.01e-6;

% FE 参数（沿用你当前版本）
Pr    = 25;
tauo  = 1e-9;
alpha = 2;
bet   = 2;
epife = 28;
Ndom  = 10000;

% 读出参数
VG = -0.5:0.02:1.7;
VD = 0.05;
VS = 0;
Iref = 1e-7;
VFB_read = -0.5;

% MOS / stack 参数
tfe = 0.8e-6;
til = 1e-7;
miu = 50;
Na  = 3e17;
T   = 300;
W   = 1;
L   = 1;

% 早期加密读点
Nread = [1 2 3 4 5 10 20 50 100].';
Nread = unique(Nread);

fprintf('=== Sweep setup ===\n');
fprintf('pw      = %g s\n', pw);
fprintf('delay   = %g s\n', delay);
fprintf('transit = %g s\n', transit);
fprintf('cycle   = %d\n', cycle);
fprintf('step    = %g s\n', step);
fprintf('Ndom    = %d\n', Ndom);
fprintf('Iref    = %g A\n', Iref);
fprintf('VG      = [%g : %g : %g]\n\n', VG(1), VG(2)-VG(1), VG(end));

%% =========================
% 2) 固定随机畴样本
%    说明：
%    为了比较不同 amp，只改激励，不改随机样本
%% =========================
Weight  = ones(Ndom,1) / Ndom;
St_init = -ones(Ndom,1);

Ea_mean  = 1.0;
Ea_sigma = 0.5;
r_Ea = Ea_mean + Ea_sigma * randn(Ndom,1);
r_Ea(r_Ea <= 0) = Ea_mean;

% 先关掉 offset，避免 gap 段额外演化
r_voff = zeros(Ndom,1);

%% =========================
% 3) 预分配结果
% =========================
nAmp = numel(amp_list);
nRead = numel(Nread);

Vth_mat      = nan(nRead, nAmp);
DeltaVth_mat = nan(nRead, nAmp);
Pread_mat    = nan(nRead, nAmp);

Vth0_vec     = nan(nAmp,1);
P0_vec       = nan(nAmp,1);

dV1_vec      = nan(nAmp,1);   % DeltaVth at N=1
dVlast_vec   = nan(nAmp,1);   % DeltaVth at N=100
gradual_vec  = nan(nAmp,1);   % |DeltaVth(100)-DeltaVth(1)|
valid_ratio  = nan(nAmp,1);   % 有效读点比例

%% =========================
% 4) 主循环：扫 amp
% =========================
for ia = 1:nAmp
    amp = amp_list(ia);

    fprintf('\n==================================================\n');
    fprintf('Running case %d / %d : amp = %.3f V\n', ia, nAmp, amp);
    fprintf('==================================================\n');

    %% 4.1 生成波形
    [time, volt, index, index_pre, index_end] = ...
        wfdef_acc_fix(amp, pw, step, delay, transit, cycle);

    index     = index(:);
    index_pre = index_pre(:);
    index_end = index_end(:);

    %% 4.2 FE 状态演化
    tic;
    [vfev, Stsum] = get_FE_state(time, volt, St_init, Weight, r_Ea, r_voff, ...
                                 Pr, tauo, alpha, bet, epife, Ndom);
    telapsed = toc;

    fprintf('get_FE_state done, elapsed = %.3f s\n', telapsed);

    %% 4.3 取单条轨迹
    if isvector(Stsum)
        st_trace = Stsum(:);
    else
        st_trace = Stsum(:,1);
    end

    %% 4.4 参考点 N=0
    idx0 = index_pre(1);
    P0 = st_trace(idx0);

    [Vth0, ~, ~, ~, valid0] = local_extract_vth( ...
        P0, epife, tfe, til, miu, Na, T, W, L, VG, VFB_read, VD, VS, Iref);

    P0_vec(ia) = P0;
    Vth0_vec(ia) = Vth0;

    %% 4.5 各读点
    Vth_case      = nan(nRead,1);
    DeltaVth_case = nan(nRead,1);
    Pread_case    = nan(nRead,1);
    valid_case    = false(nRead,1);

    for k = 1:nRead
        n = Nread(k);
        if n > numel(index)
            continue;
        end

        idx = index(n);
        if idx < 1 || idx > numel(st_trace)
            continue;
        end

        Pread_case(k) = st_trace(idx);

        [Vth_case(k), ~, ~, ~, valid_case(k)] = local_extract_vth( ...
            Pread_case(k), epife, tfe, til, miu, Na, T, W, L, ...
            VG, VFB_read, VD, VS, Iref);
    end

    if valid0
        DeltaVth_case = Vth_case - Vth0;
    else
        % 如果初始点提不出来，就退化到第一个有效点
        first_valid = find(valid_case, 1, 'first');
        if ~isempty(first_valid)
            DeltaVth_case = Vth_case - Vth_case(first_valid);
        end
    end

    %% 4.6 保存到总矩阵
    Vth_mat(:, ia)      = Vth_case;
    DeltaVth_mat(:, ia) = DeltaVth_case;
    Pread_mat(:, ia)    = Pread_case;

    valid_ratio(ia) = mean(valid_case);

    % N=1 与 N=100 的指标
    idx1 = find(Nread == 1, 1);
    idxLast = find(Nread == 100, 1);

    if ~isempty(idx1)
        dV1_vec(ia) = DeltaVth_case(idx1);
    end
    if ~isempty(idxLast)
        dVlast_vec(ia) = DeltaVth_case(idxLast);
    end
    if ~isempty(idx1) && ~isempty(idxLast)
        gradual_vec(ia) = abs(DeltaVth_case(idxLast) - DeltaVth_case(idx1));
    end

    %% 4.7 保存单 case
    casefile = fullfile(outdir, sprintf('case_amp_%0.3fV.mat', amp));
    save(casefile, ...
        'amp', 'pw', 'delay', 'transit', 'cycle', 'step', ...
        'Pr', 'tauo', 'alpha', 'bet', 'epife', 'Ndom', ...
        'VG', 'VD', 'VS', 'Iref', 'VFB_read', ...
        'time', 'volt', 'index', 'index_pre', 'index_end', ...
        'Stsum', 'vfev', 'P0', 'Vth0', ...
        'Nread', 'Pread_case', 'Vth_case', 'DeltaVth_case', 'valid_case');

    %% 4.8 终端输出
    fprintf('P0   = %.6f\n', P0);
    fprintf('Vth0 = %.6f V\n', Vth0);
    for k = 1:nRead
        fprintf('  N=%-4d  P=%.6e  Vth=%.6f  dV=%.6f  valid=%d\n', ...
            Nread(k), Pread_case(k), Vth_case(k), DeltaVth_case(k), valid_case(k));
    end
end

%% =========================
% 5) 汇总表
% =========================
summary_table = table( ...
    amp_list(:), P0_vec, Vth0_vec, dV1_vec, dVlast_vec, gradual_vec, valid_ratio, ...
    'VariableNames', {'amp_V', 'P0', 'Vth0_V', 'dV_N1_V', 'dV_N100_V', 'gradual_metric_V', 'valid_ratio'});

disp(' ');
disp('================ summary_table ================');
disp(summary_table);

save(fullfile(outdir, 'summary_table.mat'), ...
     'summary_table', 'amp_list', 'Nread', 'Vth_mat', 'DeltaVth_mat', 'Pread_mat');

writetable(summary_table, fullfile(outdir, 'summary_table.csv'));

%% =========================
% 6) 图 1：DeltaVth vs N（不同 amp）
% =========================
figure;
hold on;
for ia = 1:nAmp
    semilogx(Nread, DeltaVth_mat(:,ia), 'o-', 'LineWidth', 1.3);
end
xlabel('Cumulative pulse count N');
ylabel('\Delta V_{th} (V)');
title('Sweep of pulse amplitude: \DeltaV_{th}(N)');
grid on;
legend(arrayfun(@(x) sprintf('amp=%.2f V', x), amp_list, 'UniformOutput', false), ...
       'Location', 'best');
hold off;

%% =========================
% 7) 图 2：Pread vs N（不同 amp）
% =========================
figure;
hold on;
for ia = 1:nAmp
    semilogx(Nread, Pread_mat(:,ia), 's-', 'LineWidth', 1.3);
end
xlabel('Cumulative pulse count N');
ylabel('P_{read} (\muC/cm^2)');
title('Sweep of pulse amplitude: P_{read}(N)');
grid on;
legend(arrayfun(@(x) sprintf('amp=%.2f V', x), amp_list, 'UniformOutput', false), ...
       'Location', 'best');
hold off;

%% =========================
% 8) 图 3：找“非饱和窗口”的判据
% =========================
figure;
plot(amp_list, abs(dV1_vec), 'o-', 'LineWidth', 1.5); hold on;
plot(amp_list, abs(dVlast_vec), 's-', 'LineWidth', 1.5);
plot(amp_list, gradual_vec, 'd-', 'LineWidth', 1.5);
xlabel('Pulse amplitude amp (V)');
ylabel('Magnitude (V)');
title('Amplitude sweep metrics');
legend('|dV(N=1)|', '|dV(N=100)|', '|dV(100)-dV(1)|', 'Location', 'best');
grid on;
hold off;

%% =========================
% 9) 自动给出推荐
% =========================
% 经验判据：
%   - 首脉冲不能太小：|dV1| > 5 mV
%   - 首脉冲不能太大：|dV1| < 100 mV
%   - 需要渐进性：gradual_metric > 10 mV
candidate = find( ...
    abs(dV1_vec) > 0.005 & ...
    abs(dV1_vec) < 0.10  & ...
    gradual_vec  > 0.010 & ...
    valid_ratio  > 0.8);

fprintf('\n================ recommendation ================\n');
if isempty(candidate)
    fprintf('没有自动找到理想的“非饱和累积窗口”。\n');
    fprintf('你需要考虑：\n');
    fprintf('1) 降低 amp 扫描下限，或者\n');
    fprintf('2) 改扫 pw（更常见），或者\n');
    fprintf('3) 调整 tauo / alpha / r_Ea 分布。\n');
else
    fprintf('推荐优先检查以下 amp:\n');
    for ii = 1:numel(candidate)
        ia = candidate(ii);
        fprintf('  amp = %.3f V, dV1 = %.4f V, dV100 = %.4f V, gradual = %.4f V\n', ...
            amp_list(ia), dV1_vec(ia), dVlast_vec(ia), gradual_vec(ia));
    end
end

%% =========================================================
% Local function
%% =========================================================
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

    if abs(x(2)-x(1)) < 1e-30
        Vth = mean(y);
    else
        Vth = interp1(x, y, log10(Iref), 'linear');
    end

    is_valid = isfinite(Vth);
end