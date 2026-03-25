%% =========================================================
% sweep_delay_near_transition.m
%
% 目的：
%   在幅值过渡区附近（amp ~ 1.5 V）扫描 delay，
%   检查现有 scalable FeFET 模型对 interval time 是否敏感。
%
% 依赖：
%   wfdef_acc_fix.m
%   get_FE_state.m
%   get_ID.m
%
% 输出：
%   1) 每个 (amp, delay) case 单独保存
%   2) 汇总表 summary_table
%   3) 多张对比图
%
% 核心判断：
%   固定 amp 后，若不同 delay 的 DeltaVth(N) 曲线几乎重合，
%   则说明当前模型对 interval time 不敏感。
%% =========================================================

clear;
clc;
close all;

rng(0);

%% =========================
% 0) 输出目录
% =========================
outdir = 'sweep_delay_out';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

%% =========================
% 1) 固定参数
% =========================
VDD = 5;

% 重点扫过渡区附近
amp_list = [1.4 1.5 1.6];

% 扫 delay
delay_list = [1e-7 3e-7 1e-6 3e-6 1e-5 3e-5 1e-4];

pw      = 0.1e-6;
transit = 0.02e-6;
cycle   = 100;
step    = 0.01e-6;

% FE 参数
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

fprintf('=== Delay sweep near transition ===\n');
fprintf('amp_list   = '); fprintf('%g ', amp_list); fprintf('\n');
fprintf('delay_list = '); fprintf('%g ', delay_list); fprintf('\n');
fprintf('pw         = %g s\n', pw);
fprintf('transit    = %g s\n', transit);
fprintf('cycle      = %d\n', cycle);
fprintf('step       = %g s\n', step);
fprintf('Ndom       = %d\n', Ndom);
fprintf('Iref       = %g A\n\n', Iref);

%% =========================
% 2) 固定随机畴样本
%    为了公平比较 delay，只改激励间隔
%% =========================
Weight  = ones(Ndom,1) / Ndom;
St_init = -ones(Ndom,1);

Ea_mean  = 1.0;
Ea_sigma = 0.5;
r_Ea = Ea_mean + Ea_sigma * randn(Ndom,1);
r_Ea(r_Ea <= 0) = Ea_mean;

% 先关掉 offset，聚焦“原始主干模型”的 gap 效应
r_voff = zeros(Ndom,1);

%% =========================
% 3) 预分配
% =========================
nAmp   = numel(amp_list);
nDelay = numel(delay_list);
nRead  = numel(Nread);

% 三维结果矩阵：Nread x nDelay x nAmp
Vth_cube      = nan(nRead, nDelay, nAmp);
DeltaVth_cube = nan(nRead, nDelay, nAmp);
Pread_cube    = nan(nRead, nDelay, nAmp);

Vth0_mat      = nan(nDelay, nAmp);
P0_mat        = nan(nDelay, nAmp);
valid_ratio   = nan(nDelay, nAmp);

dV1_mat       = nan(nDelay, nAmp);
dV100_mat     = nan(nDelay, nAmp);
gradual_mat   = nan(nDelay, nAmp);

% delay 灵敏度指标：
% 在同一 amp 下，用不同 delay 的 dV100 的 span 来衡量
% 后面再算
delay_sensitivity_dV100 = nan(nAmp,1);
delay_sensitivity_dV1   = nan(nAmp,1);

%% =========================
% 4) 双循环：amp x delay
% =========================
for ia = 1:nAmp
    amp = amp_list(ia);

    for idl = 1:nDelay
        delay = delay_list(idl);

        fprintf('\n==================================================\n');
        fprintf('Running amp = %.3f V, delay = %.3e s\n', amp, delay);
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

        %% 4.3 单条轨迹
        if isvector(Stsum)
            st_trace = Stsum(:);
        else
            st_trace = Stsum(:,1);
        end

        %% 4.4 初始参考点 N=0
        idx0 = index_pre(1);
        P0 = st_trace(idx0);

        [Vth0, ~, ~, ~, valid0] = local_extract_vth( ...
            P0, epife, tfe, til, miu, Na, T, W, L, ...
            VG, VFB_read, VD, VS, Iref);

        P0_mat(idl, ia) = P0;
        Vth0_mat(idl, ia) = Vth0;

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
            first_valid = find(valid_case, 1, 'first');
            if ~isempty(first_valid)
                DeltaVth_case = Vth_case - Vth_case(first_valid);
            end
        end

        %% 4.6 写入总矩阵
        Vth_cube(:, idl, ia)      = Vth_case;
        DeltaVth_cube(:, idl, ia) = DeltaVth_case;
        Pread_cube(:, idl, ia)    = Pread_case;

        valid_ratio(idl, ia) = mean(valid_case);

        idx1 = find(Nread == 1, 1);
        idx100 = find(Nread == 100, 1);

        if ~isempty(idx1)
            dV1_mat(idl, ia) = DeltaVth_case(idx1);
        end
        if ~isempty(idx100)
            dV100_mat(idl, ia) = DeltaVth_case(idx100);
        end
        if ~isempty(idx1) && ~isempty(idx100)
            gradual_mat(idl, ia) = abs(DeltaVth_case(idx100) - DeltaVth_case(idx1));
        end

        %% 4.7 保存单 case
        casefile = fullfile(outdir, sprintf('case_amp_%0.3fV_delay_%0.3es.mat', amp, delay));
        save(casefile, ...
            'amp', 'delay', 'pw', 'transit', 'cycle', 'step', ...
            'Pr', 'tauo', 'alpha', 'bet', 'epife', 'Ndom', ...
            'VG', 'VD', 'VS', 'Iref', 'VFB_read', ...
            'time', 'volt', 'index', 'index_pre', 'index_end', ...
            'Stsum', 'vfev', 'P0', 'Vth0', ...
            'Nread', 'Pread_case', 'Vth_case', 'DeltaVth_case', 'valid_case');

        %% 4.8 简要输出
        fprintf('P0   = %.6f\n', P0);
        fprintf('Vth0 = %.6f V\n', Vth0);

        if ~isempty(idx1)
            fprintf('dV(N=1)   = %.6f V\n', dV1_mat(idl, ia));
        end
        if ~isempty(idx100)
            fprintf('dV(N=100) = %.6f V\n', dV100_mat(idl, ia));
        end
        fprintf('valid_ratio = %.3f\n', valid_ratio(idl, ia));
    end
end

%% =========================
% 5) 汇总表
% =========================
rows = nAmp * nDelay;

amp_col        = nan(rows,1);
delay_col      = nan(rows,1);
P0_col         = nan(rows,1);
Vth0_col       = nan(rows,1);
dV1_col        = nan(rows,1);
dV100_col      = nan(rows,1);
gradual_col    = nan(rows,1);
valid_col      = nan(rows,1);

rr = 0;
for ia = 1:nAmp
    for idl = 1:nDelay
        rr = rr + 1;
        amp_col(rr)     = amp_list(ia);
        delay_col(rr)   = delay_list(idl);
        P0_col(rr)      = P0_mat(idl, ia);
        Vth0_col(rr)    = Vth0_mat(idl, ia);
        dV1_col(rr)     = dV1_mat(idl, ia);
        dV100_col(rr)   = dV100_mat(idl, ia);
        gradual_col(rr) = gradual_mat(idl, ia);
        valid_col(rr)   = valid_ratio(idl, ia);
    end
end

summary_table = table( ...
    amp_col, delay_col, P0_col, Vth0_col, dV1_col, dV100_col, gradual_col, valid_col, ...
    'VariableNames', {'amp_V', 'delay_s', 'P0', 'Vth0_V', 'dV_N1_V', 'dV_N100_V', 'gradual_metric_V', 'valid_ratio'});

disp(' ');
disp('================ summary_table ================');
disp(summary_table);

save(fullfile(outdir, 'summary_table.mat'), ...
     'summary_table', 'amp_list', 'delay_list', 'Nread', ...
     'Vth_cube', 'DeltaVth_cube', 'Pread_cube', ...
     'P0_mat', 'Vth0_mat', 'dV1_mat', 'dV100_mat', 'gradual_mat', 'valid_ratio');

writetable(summary_table, fullfile(outdir, 'summary_table.csv'));

%% =========================
% 6) delay 灵敏度指标
% =========================
for ia = 1:nAmp
    d1 = dV1_mat(:, ia);
    d100 = dV100_mat(:, ia);

    if any(isfinite(d1))
        delay_sensitivity_dV1(ia) = max(d1) - min(d1);
    end
    if any(isfinite(d100))
        delay_sensitivity_dV100(ia) = max(d100) - min(d100);
    end
end

sens_table = table( ...
    amp_list(:), delay_sensitivity_dV1(:), delay_sensitivity_dV100(:), ...
    'VariableNames', {'amp_V', 'delay_span_dV1_V', 'delay_span_dV100_V'});

disp(' ');
disp('================ sens_table ================');
disp(sens_table);

writetable(sens_table, fullfile(outdir, 'delay_sensitivity_table.csv'));

%% =========================
% 7) 图 1：同一 amp 下，不同 delay 的 DeltaVth(N)
% =========================
for ia = 1:nAmp
    figure;
    hold on;
    for idl = 1:nDelay
        semilogx(Nread, DeltaVth_cube(:,idl,ia), 'o-', 'LineWidth', 1.3);
    end
    xlabel('Cumulative pulse count N');
    ylabel('\Delta V_{th} (V)');
    title(sprintf('Delay sweep at amp = %.2f V', amp_list(ia)));
    grid on;
    legend(arrayfun(@(x) sprintf('delay=%.1e s', x), delay_list, 'UniformOutput', false), ...
           'Location', 'best');
    hold off;
end

%% =========================
% 8) 图 2：同一 amp 下，不同 delay 的 Pread(N)
% =========================
for ia = 1:nAmp
    figure;
    hold on;
    for idl = 1:nDelay
        semilogx(Nread, Pread_cube(:,idl,ia), 's-', 'LineWidth', 1.3);
    end
    xlabel('Cumulative pulse count N');
    ylabel('P_{read} (\muC/cm^2)');
    title(sprintf('Delay sweep P_{read}(N) at amp = %.2f V', amp_list(ia)));
    grid on;
    legend(arrayfun(@(x) sprintf('delay=%.1e s', x), delay_list, 'UniformOutput', false), ...
           'Location', 'best');
    hold off;
end

%% =========================
% 9) 图 3：不同 amp 的 delay 灵敏度
% =========================
figure;
plot(amp_list, abs(delay_sensitivity_dV1), 'o-', 'LineWidth', 1.5); hold on;
plot(amp_list, abs(delay_sensitivity_dV100), 's-', 'LineWidth', 1.5);
xlabel('amp (V)');
ylabel('Span across delays (V)');
title('Delay sensitivity metrics');
legend('span of dV(N=1)', 'span of dV(N=100)', 'Location', 'best');
grid on;
hold off;

%% =========================
% 10) 自动解释
% =========================
fprintf('\n================ interpretation helper ================\n');
for ia = 1:nAmp
    fprintf('amp = %.3f V:\n', amp_list(ia));
    fprintf('  span[dV(N=1)]   = %.6f V\n', delay_sensitivity_dV1(ia));
    fprintf('  span[dV(N=100)] = %.6f V\n', delay_sensitivity_dV100(ia));

    if isfinite(delay_sensitivity_dV100(ia))
        if abs(delay_sensitivity_dV100(ia)) < 0.005
            fprintf('  -> 对 delay 基本不敏感（100 次脉冲后差异 < 5 mV）\n');
        elseif abs(delay_sensitivity_dV100(ia)) < 0.02
            fprintf('  -> 对 delay 有弱敏感性（5~20 mV）\n');
        else
            fprintf('  -> 对 delay 有明显敏感性（> 20 mV）\n');
        end
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