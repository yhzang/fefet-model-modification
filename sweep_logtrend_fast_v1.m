%% =========================================================
% sweep_logtrend_fast_v1.m
%
% 目标：
%   快速寻找“小电压脉冲下 Vth 变化趋势”更接近目标的参数：
%   - 前期变化小
%   - 中段逐渐抬升
%   - 后段开始变缓（向饱和趋势靠近）
%
% 已知前提：
%   - 当前模型对 interval/delay 感知弱
%   - 本轮不扫 delay
%   - 不直接展开 1e6 pulses，而用 log 采样趋势重构
%
% 两阶段：
%   Stage-1: Ndom=3000 快筛
%   Stage-2: TopK 用 Ndom=10000 复核
%
% 依赖：
%   wfdef_acc_fix.m
%   get_FE_state.m
%   get_ID.m
%
% 输出：
%   logtrend_fast_v1_out/
%       stage1_all.csv
%       stage1_top.csv
%       stage2_all.csv
%       stage2_top.csv
%       workspace_results.mat
%       figures/
%% =========================================================

clear;
clc;
close all;

rng(0);

%% =========================
% 0) 输出目录
% =========================
outdir = 'logtrend_fast_v1_out';
figdir = fullfile(outdir, 'figures');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
if ~exist(figdir, 'dir')
    mkdir(figdir);
end

diary(fullfile(outdir, 'run_log.txt'));
fprintf('=== sweep_logtrend_fast_v1 ===\n');

%% =========================
% 1) 固定“中等读出映射”
%    不走最强 mapping gain，避免前期过快放大
% =========================
fixed_map = struct();
fixed_map.tfe         = 1.2e-6;
fixed_map.epife_read  = 22;
fixed_map.Iref        = 1e-7;
fixed_map.VFB_read    = -0.5;
fixed_map.VG_name     = 'VGm0p8_to_2p2';

%% =========================
% 2) 本轮主扫参数
%    小范围、快速反馈
% =========================
amp_list      = [1.4 1.6 1.8 2.0];
pw_list       = [0.02e-6 0.05e-6 0.10e-6];
Ea_mean_list  = [0.95 1.00 1.05];
Ea_sigma_list = [0.18 0.20];

%% =========================
% 3) 固定 pulse / FE / MOS 参数
% =========================
delay   = 1e-6;      % 固定，不扫
transit = 0.02e-6;
step    = 0.01e-6;

% FE 参数
Pr    = 25;
tauo  = 1e-9;
alpha = 2;
bet   = 2;

% 写入主干 FE eps
epife_write = 28;

% 暂不引入 voff 演化
r_voff = [];

% MOS / stack 参数
til = 1e-7;
miu = 50;
Na  = 3e17;
T0  = 300;
W   = 1;
L   = 1;
VD  = 0.05;
VS  = 0;

%% =========================
% 4) 两阶段设置
% =========================
% Stage-1：快速粗筛
Ndom_stage1   = 3000;
Nread_stage1  = [1 10 100 1000].';
topK_stage1   = 8;

% Stage-2：精筛复核
Ndom_stage2   = 10000;
Nread_stage2  = [1 3 10 30 100 300 1000 3000].';

%% =========================
% 5) 构造全局固定随机样本（不同 case 共用 z-score）
%    这样比较不同 Ea_mean / Ea_sigma 更公平
% =========================
z_base_stage1 = randn(Ndom_stage1,1);
z_base_stage2 = randn(Ndom_stage2,1);

%% =========================
% 6) 生成 case 列表
% =========================
case_list = {};
icase = 0;

for ia = 1:numel(amp_list)
    for ipw = 1:numel(pw_list)
        for iem = 1:numel(Ea_mean_list)
            for ies = 1:numel(Ea_sigma_list)

                icase = icase + 1;

                cfg = struct();
                cfg.case_idx    = icase;
                cfg.amp         = amp_list(ia);
                cfg.pw          = pw_list(ipw);
                cfg.Ea_mean     = Ea_mean_list(iem);
                cfg.Ea_sigma    = Ea_sigma_list(ies);

                cfg.tfe         = fixed_map.tfe;
                cfg.epife_read  = fixed_map.epife_read;
                cfg.Iref        = fixed_map.Iref;
                cfg.VFB_read    = fixed_map.VFB_read;
                cfg.VG_name     = fixed_map.VG_name;

                cfg.case_id     = make_case_id(cfg);

                case_list{end+1,1} = cfg;
            end
        end
    end
end

nCase = numel(case_list);
fprintf('Total stage-1 cases = %d\n', nCase);

%% =========================
% 7) Stage-1 快筛
% =========================
stage1_results = repmat(init_result_stage1(), nCase, 1);

for ic = 1:nCase
    cfg = case_list{ic};

    fprintf('\n[Stage-1 %d / %d] %s\n', ic, nCase, cfg.case_id);

    try
        % 7.1 生成 r_Ea
        r_Ea = cfg.Ea_mean + cfg.Ea_sigma * z_base_stage1;
        r_Ea(r_Ea <= 0) = cfg.Ea_mean;

        Weight  = ones(Ndom_stage1,1) / Ndom_stage1;
        St_init = -ones(Ndom_stage1,1);
        r_voff  = zeros(Ndom_stage1,1);

        % 7.2 对不同 Nread 逐个从头跑
        [dV_vec, A_vec, valid_vec, P_vec, Vth_vec, elapsed_sum] = ...
            run_case_for_Nread(cfg, Nread_stage1, ...
                               step, delay, transit, ...
                               St_init, Weight, r_Ea, r_voff, ...
                               Pr, tauo, alpha, bet, epife_write, ...
                               til, miu, Na, T0, W, L, VD, VS);

        dV1    = dV_vec(1);
        dV10   = dV_vec(2);
        dV100  = dV_vec(3);
        dV1000 = dV_vec(4);

        A1    = A_vec(1);
        A10   = A_vec(2);
        A100  = A_vec(3);
        A1000 = A_vec(4);

        ratio1   = A1   / max(A1000, 1e-15);
        ratio10  = A10  / max(A1000, 1e-15);
        ratio100 = A100 / max(A1000, 1e-15);

        g1 = A10   - A1;
        g2 = A100  - A10;
        g3 = A1000 - A100;

        is_invalid = any(~isfinite(A_vec));
        is_mono    = (~is_invalid) && all(diff(A_vec) >= -1e-6);

        % 启发式 1e6 外推（保守）
        % 用 100 -> 1000 这一 decade 的增长估计，再乘一个保守系数 1.5
        tail_slope = max(A1000 - A100, 0);
        A1e6_heur  = A1000 + 1.5 * tail_slope;

        trend_score = compute_stage1_score(A1, A10, A100, A1000, ...
                                           ratio1, ratio10, ratio100, ...
                                           g1, g2, g3, is_mono, A1e6_heur);

        stage1_results(ic).case_id      = string(cfg.case_id);
        stage1_results(ic).amp          = cfg.amp;
        stage1_results(ic).pw           = cfg.pw;
        stage1_results(ic).Ea_mean      = cfg.Ea_mean;
        stage1_results(ic).Ea_sigma     = cfg.Ea_sigma;

        stage1_results(ic).tfe          = cfg.tfe;
        stage1_results(ic).epife_read   = cfg.epife_read;
        stage1_results(ic).Iref         = cfg.Iref;
        stage1_results(ic).VFB_read     = cfg.VFB_read;
        stage1_results(ic).VG_name      = string(cfg.VG_name);

        stage1_results(ic).dV1          = dV1;
        stage1_results(ic).dV10         = dV10;
        stage1_results(ic).dV100        = dV100;
        stage1_results(ic).dV1000       = dV1000;

        stage1_results(ic).A1           = A1;
        stage1_results(ic).A10          = A10;
        stage1_results(ic).A100         = A100;
        stage1_results(ic).A1000        = A1000;

        stage1_results(ic).ratio1       = ratio1;
        stage1_results(ic).ratio10      = ratio10;
        stage1_results(ic).ratio100     = ratio100;

        stage1_results(ic).g1_1_10      = g1;
        stage1_results(ic).g2_10_100    = g2;
        stage1_results(ic).g3_100_1000  = g3;

        stage1_results(ic).A1e6_heur    = A1e6_heur;
        stage1_results(ic).is_mono      = is_mono;
        stage1_results(ic).trend_score  = trend_score;
        stage1_results(ic).elapsed_s    = elapsed_sum;

    catch ME
        warning('Stage-1 case failed: %s', cfg.case_id);
        warning('%s', ME.message);

        stage1_results(ic).case_id     = string(cfg.case_id);
        stage1_results(ic).amp         = cfg.amp;
        stage1_results(ic).pw          = cfg.pw;
        stage1_results(ic).Ea_mean     = cfg.Ea_mean;
        stage1_results(ic).Ea_sigma    = cfg.Ea_sigma;
        stage1_results(ic).trend_score = -1e9;
    end
end

T1 = struct2table(stage1_results);
writetable(T1, fullfile(outdir, 'stage1_all.csv'));

T1_valid = T1(isfinite(T1.trend_score), :);
T1_sort  = sortrows(T1_valid, 'trend_score', 'descend');
writetable(T1_sort(1:min(height(T1_sort), 20), :), fullfile(outdir, 'stage1_top.csv'));

disp(' ');
disp('================ Stage-1 Top by trend_score ================');
disp(T1_sort(1:min(height(T1_sort), 20), :));

%% =========================
% 8) Stage-2 只复核 TopK
% =========================
topK = min(topK_stage1, height(T1_sort));
top_case_ids = string(T1_sort.case_id(1:topK));

stage2_cfg_list = {};
for i = 1:numel(case_list)
    if any(string(case_list{i}.case_id) == top_case_ids)
        stage2_cfg_list{end+1,1} = case_list{i};
    end
end

nCase2 = numel(stage2_cfg_list);
fprintf('\nTotal stage-2 cases = %d\n', nCase2);

stage2_results = repmat(init_result_stage2(), nCase2, 1);

for ic = 1:nCase2
    cfg = stage2_cfg_list{ic};

    fprintf('\n[Stage-2 %d / %d] %s\n', ic, nCase2, cfg.case_id);

    try
        r_Ea = cfg.Ea_mean + cfg.Ea_sigma * z_base_stage2;
        r_Ea(r_Ea <= 0) = cfg.Ea_mean;

        Weight  = ones(Ndom_stage2,1) / Ndom_stage2;
        St_init = -ones(Ndom_stage2,1);
        r_voff  = zeros(Ndom_stage2,1);

        [dV_vec, A_vec, valid_vec, P_vec, Vth_vec, elapsed_sum] = ...
            run_case_for_Nread(cfg, Nread_stage2, ...
                               step, delay, transit, ...
                               St_init, Weight, r_Ea, r_voff, ...
                               Pr, tauo, alpha, bet, epife_write, ...
                               til, miu, Na, T0, W, L, VD, VS);

        % 关键点
        A1    = A_vec(1);
        A3    = A_vec(2);
        A10   = A_vec(3);
        A30   = A_vec(4);
        A100  = A_vec(5);
        A300  = A_vec(6);
        A1000 = A_vec(7);
        A3000 = A_vec(8);

        g1 = A10   - A1;
        g2 = A100  - A10;
        g3 = A1000 - A100;
        g4 = A3000 - A1000;

        ratio1   = A1   / max(A3000, 1e-15);
        ratio10  = A10  / max(A3000, 1e-15);
        ratio100 = A100 / max(A3000, 1e-15);

        is_invalid = any(~isfinite(A_vec));
        is_mono    = (~is_invalid) && all(diff(A_vec) >= -1e-6);

        % “你想要的趋势”启发式：
        % 前期小、中段涨、后段开始变缓
        trend_good = false;
        if is_mono
            cond1 = (A10 <= 0.12);           % 前 10 pulse 先小一点
            cond2 = (g2 > g1);               % 10->100 增长比 1->10 更强
            cond3 = (g3 > 0);                % 100->1000 继续涨
            cond4 = (g4 <= g3);              % 1000->3000 开始变缓
            cond5 = (A3000 >= 0.20);         % 到 3000 pulse 至少看得见
            trend_good = cond1 && cond2 && cond3 && cond4 && cond5;
        end

        % 保守 1e6 趋势重构
        % 用最后一段增长继续外推，但乘以 0.8，且不小于 A3000
        tail_per_dec = max(g4, 0) / max(log10(3000) - log10(1000), 1e-12);
        A1e6_pred = A3000 + 0.8 * tail_per_dec * (6 - log10(3000));

        trend_score2 = compute_stage2_score(A1, A10, A100, A1000, A3000, ...
                                            ratio1, ratio10, ratio100, ...
                                            g1, g2, g3, g4, is_mono, trend_good, A1e6_pred);

        stage2_results(ic).case_id      = string(cfg.case_id);
        stage2_results(ic).amp          = cfg.amp;
        stage2_results(ic).pw           = cfg.pw;
        stage2_results(ic).Ea_mean      = cfg.Ea_mean;
        stage2_results(ic).Ea_sigma     = cfg.Ea_sigma;

        stage2_results(ic).tfe          = cfg.tfe;
        stage2_results(ic).epife_read   = cfg.epife_read;
        stage2_results(ic).Iref         = cfg.Iref;
        stage2_results(ic).VFB_read     = cfg.VFB_read;
        stage2_results(ic).VG_name      = string(cfg.VG_name);

        stage2_results(ic).A1           = A1;
        stage2_results(ic).A3           = A3;
        stage2_results(ic).A10          = A10;
        stage2_results(ic).A30          = A30;
        stage2_results(ic).A100         = A100;
        stage2_results(ic).A300         = A300;
        stage2_results(ic).A1000        = A1000;
        stage2_results(ic).A3000        = A3000;

        stage2_results(ic).ratio1       = ratio1;
        stage2_results(ic).ratio10      = ratio10;
        stage2_results(ic).ratio100     = ratio100;

        stage2_results(ic).g1_1_10      = g1;
        stage2_results(ic).g2_10_100    = g2;
        stage2_results(ic).g3_100_1000  = g3;
        stage2_results(ic).g4_1000_3000 = g4;

        stage2_results(ic).trend_good   = trend_good;
        stage2_results(ic).is_mono      = is_mono;
        stage2_results(ic).A1e6_pred    = A1e6_pred;
        stage2_results(ic).trend_score2 = trend_score2;
        stage2_results(ic).elapsed_s    = elapsed_sum;

        % 保存整条曲线，便于后面画图
        stage2_results(ic).Nread_str    = vec2str(Nread_stage2);
        stage2_results(ic).Avec_str     = vec2str(A_vec);

    catch ME
        warning('Stage-2 case failed: %s', cfg.case_id);
        warning('%s', ME.message);

        stage2_results(ic).case_id      = string(cfg.case_id);
        stage2_results(ic).amp          = cfg.amp;
        stage2_results(ic).pw           = cfg.pw;
        stage2_results(ic).Ea_mean      = cfg.Ea_mean;
        stage2_results(ic).Ea_sigma     = cfg.Ea_sigma;
        stage2_results(ic).trend_score2 = -1e9;
    end
end

T2 = struct2table(stage2_results);
writetable(T2, fullfile(outdir, 'stage2_all.csv'));

T2_valid = T2(isfinite(T2.trend_score2), :);
T2_sort  = sortrows(T2_valid, 'trend_score2', 'descend');
writetable(T2_sort(1:min(height(T2_sort), 20), :), fullfile(outdir, 'stage2_top.csv'));

disp(' ');
disp('================ Stage-2 Top by trend_score2 ================');
disp(T2_sort(1:min(height(T2_sort), 20), :));

%% =========================
% 9) 画 Stage-2 Top 曲线
% =========================
nplot = min(5, height(T2_sort));

for i = 1:nplot
    row = T2_sort(i,:);
    Nvec = str2num(row.Nread_str{1}); %#ok<ST2NM>
    Avec = str2num(row.Avec_str{1}); %#ok<ST2NM>

    figure('Visible', 'off');
    semilogx(Nvec, Avec, 'o-', 'LineWidth', 1.5);
    xlabel('Pulse count N');
    ylabel('|DeltaVth| (V)');
    title(sprintf('Top-%d: %s', i, row.case_id{1}), 'Interpreter', 'none');
    grid on;
    saveas(gcf, fullfile(figdir, sprintf('top%d_curve.png', i)));
    close(gcf);
end

%% =========================
% 10) 保存工作区
% =========================
save(fullfile(outdir, 'workspace_results.mat'), ...
    'T1', 'T1_sort', 'T2', 'T2_sort', ...
    'stage1_results', 'stage2_results', 'case_list', ...
    'amp_list', 'pw_list', 'Ea_mean_list', 'Ea_sigma_list', ...
    'fixed_map', 'Nread_stage1', 'Nread_stage2', ...
    'Ndom_stage1', 'Ndom_stage2', ...
    'delay', 'transit', 'step', ...
    'Pr', 'tauo', 'alpha', 'bet', 'epife_write', ...
    'til', 'miu', 'Na', 'T0', 'W', 'L', 'VD', 'VS');

fprintf('\nAll done. Results saved in folder: %s\n', outdir);
diary off;

%% =========================================================
% local functions
%% =========================================================

function s = init_result_stage1()
s = struct( ...
    'case_id', "", ...
    'amp', NaN, ...
    'pw', NaN, ...
    'Ea_mean', NaN, ...
    'Ea_sigma', NaN, ...
    'tfe', NaN, ...
    'epife_read', NaN, ...
    'Iref', NaN, ...
    'VFB_read', NaN, ...
    'VG_name', "", ...
    'dV1', NaN, ...
    'dV10', NaN, ...
    'dV100', NaN, ...
    'dV1000', NaN, ...
    'A1', NaN, ...
    'A10', NaN, ...
    'A100', NaN, ...
    'A1000', NaN, ...
    'ratio1', NaN, ...
    'ratio10', NaN, ...
    'ratio100', NaN, ...
    'g1_1_10', NaN, ...
    'g2_10_100', NaN, ...
    'g3_100_1000', NaN, ...
    'A1e6_heur', NaN, ...
    'is_mono', false, ...
    'trend_score', NaN, ...
    'elapsed_s', NaN ...
    );
end

function s = init_result_stage2()
s = struct( ...
    'case_id', "", ...
    'amp', NaN, ...
    'pw', NaN, ...
    'Ea_mean', NaN, ...
    'Ea_sigma', NaN, ...
    'tfe', NaN, ...
    'epife_read', NaN, ...
    'Iref', NaN, ...
    'VFB_read', NaN, ...
    'VG_name', "", ...
    'A1', NaN, ...
    'A3', NaN, ...
    'A10', NaN, ...
    'A30', NaN, ...
    'A100', NaN, ...
    'A300', NaN, ...
    'A1000', NaN, ...
    'A3000', NaN, ...
    'ratio1', NaN, ...
    'ratio10', NaN, ...
    'ratio100', NaN, ...
    'g1_1_10', NaN, ...
    'g2_10_100', NaN, ...
    'g3_100_1000', NaN, ...
    'g4_1000_3000', NaN, ...
    'trend_good', false, ...
    'is_mono', false, ...
    'A1e6_pred', NaN, ...
    'trend_score2', NaN, ...
    'Nread_str', "", ...
    'Avec_str', "", ...
    'elapsed_s', NaN ...
    );
end

function score = compute_stage1_score(A1, A10, A100, A1000, ...
                                      ratio1, ratio10, ratio100, ...
                                      g1, g2, g3, is_mono, A1e6_heur)
if ~isfinite(A1000)
    score = -1e9;
    return;
end

score = 0;
score = score + 2.0 * A1000;
score = score + 1.5 * max(g2, 0);
score = score + 0.8 * max(g3, 0);
score = score - 2.0 * A10;
score = score - 0.5 * ratio1;
score = score - 0.5 * ratio10;
score = score + 0.5 * min(A1e6_heur, 0.8);

if ~is_mono
    score = score - 2.0;
end

% 想要“先小后涨”
if g2 > g1
    score = score + 0.8;
else
    score = score - 0.8;
end

% 不希望 1 pulse 就过大
if A1 > 0.08
    score = score - 1.0;
end
end

function score = compute_stage2_score(A1, A10, A100, A1000, A3000, ...
                                      ratio1, ratio10, ratio100, ...
                                      g1, g2, g3, g4, is_mono, trend_good, A1e6_pred)
if ~isfinite(A3000)
    score = -1e9;
    return;
end

score = 0;
score = score + 2.5 * A3000;
score = score + 1.0 * min(A1e6_pred, 0.8);
score = score - 2.0 * A10;
score = score - 0.8 * ratio1;

if g2 > g1
    score = score + 1.2;
else
    score = score - 1.2;
end

if g3 > 0
    score = score + 0.8;
else
    score = score - 1.0;
end

if g4 <= g3
    score = score + 0.8;
else
    score = score - 0.8;
end

if ~is_mono
    score = score - 3.0;
end

if trend_good
    score = score + 2.0;
end
end

function [dV_vec, A_vec, valid_vec, P_vec, Vth_vec, elapsed_sum] = ...
    run_case_for_Nread(cfg, Nread, ...
                       step, delay, transit, ...
                       St_init, Weight, r_Ea, r_voff, ...
                       Pr, tauo, alpha, bet, epife_write, ...
                       til, miu, Na, T0, W, L, VD, VS)

nR = numel(Nread);

dV_vec    = nan(nR,1);
A_vec     = nan(nR,1);
valid_vec = false(nR,1);
P_vec     = nan(nR,1);
Vth_vec   = nan(nR,1);

elapsed_sum = 0;

VG = make_VG_from_name(cfg.VG_name);

Vth_ref = NaN;
valid0_ref = false;

for kk = 1:nR
    Npulse = Nread(kk);

    t0 = tic;
    [time, volt, index, index_pre, index_end, info] = ...
        wfdef_acc_fix(cfg.amp, cfg.pw, step, delay, transit, Npulse);
    [vfev, Stsum] = get_FE_state(time, volt, St_init, Weight, r_Ea, r_voff, ...
                                 Pr, tauo, alpha, bet, epife_write, numel(St_init));
    elapsed_sum = elapsed_sum + toc(t0);

    if isvector(Stsum)
        st_trace = Stsum(:);
    else
        st_trace = Stsum(:,1);
    end

    % 参考点：第一次脉冲之前
    idx0 = index_pre(1);
    P0 = st_trace(idx0);

    [Vth0, ~, ~, ~, valid0] = local_extract_vth( ...
        P0, cfg.epife_read, cfg.tfe, til, miu, Na, T0, W, L, ...
        VG, cfg.VFB_read, VD, VS, cfg.Iref);

    idxr = index(end);
    Pread = st_trace(idxr);
    P_vec(kk) = Pread;

    [Vth_read, ~, ~, ~, valid_read] = local_extract_vth( ...
        Pread, cfg.epife_read, cfg.tfe, til, miu, Na, T0, W, L, ...
        VG, cfg.VFB_read, VD, VS, cfg.Iref);

    Vth_vec(kk)   = Vth_read;
    valid_vec(kk) = valid0 && valid_read;

    if valid0 && valid_read
        dV_vec(kk) = Vth_read - Vth0;
        A_vec(kk)  = abs(dV_vec(kk));
    else
        dV_vec(kk) = NaN;
        A_vec(kk)  = NaN;
    end
end
end

function [Vth, ID, idmin, idmax, is_valid] = local_extract_vth( ...
    Psum, epi_fe, tfe, til, miu, Na, T0, W, L, VG, VFB, VD, VS, Iref)

ID = get_ID(Psum, epi_fe, tfe, til, miu, Na, T0, W, L, VG, VFB, VD, VS);
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

if abs(x(2) - x(1)) < 1e-30
    Vth = mean(y);
else
    Vth = interp1(x, y, log10(Iref), 'linear');
end

is_valid = isfinite(Vth);
end

function VG = make_VG_from_name(VG_name)
switch VG_name
    case 'VGm0p8_to_2p2'
        VG = -0.8:0.02:2.2;
    case 'VGm0p5_to_1p7'
        VG = -0.5:0.02:1.7;
    case 'VGm1p0_to_2p5'
        VG = -1.0:0.02:2.5;
    otherwise
        error('Unknown VG_name: %s', VG_name);
end
end

function case_id = make_case_id(cfg)
case_id = sprintf(['A%s_PW%sus_Em%s_Es%s_tfe%s_epi%s_Iref%s_VFB%s_%s'], ...
    fmt_num(cfg.amp), ...
    fmt_num(cfg.pw * 1e6), ...
    fmt_num(cfg.Ea_mean), ...
    fmt_num(cfg.Ea_sigma), ...
    fmt_num(cfg.tfe), ...
    fmt_num(cfg.epife_read), ...
    fmt_num(cfg.Iref), ...
    fmt_num(cfg.VFB_read), ...
    cfg.VG_name);

case_id = strrep(case_id, '.', 'p');
case_id = strrep(case_id, '-', 'm');
end

function s = fmt_num(x)
if abs(x) >= 1e-3 && abs(x) < 1e3
    s = sprintf('%.2f', x);
else
    s = sprintf('%.2e', x);
end
s = strrep(s, '+0', '');
s = strrep(s, '+', '');
s = strrep(s, '-0', '-');
s = strrep(s, '.', 'p');
s = strrep(s, '-', 'm');
end

function s = vec2str(v)
s = mat2str(v(:).', 6);
end