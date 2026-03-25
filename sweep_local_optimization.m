%% =========================================================
% sweep_local_optimization.m
%
% 目的：
%   在上一次结果基础上，进一步细扫并优化，寻找更好的前期不变和后期逐渐增加的Vth趋势。
%
% 当前假设：
%   1) interval / delay 不是主导因素
%   2) 针对甜区附近再进行优化，目标是最终让 A3000 达到 0.5V 或以上。
%   3) 保持中等读出映射，不再走最强 mapping gain
%
% 依赖：
%   wfdef_acc_fix.m
%   get_FE_state.m
%   get_ID.m
%
% 输出：
%   sweep_local_optimization_out/
%       all_cases.csv
%       top_by_score.csv
%       top_by_A3000.csv
%       workspace_results.mat
%       多张图
%% =========================================================

clear;
clc;
close all;

rng(0);

%% =========================
% 0) 输出目录
% =========================
outdir = 'sweep_local_optimization_out';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

diary(fullfile(outdir, 'run_log.txt'));
fprintf('=== sweep_local_optimization ===\n');

%% =========================
% 1) 固定中等读出映射
%    先别再用最强 mapping gain
%% =========================
tfe_fixed        = 1.2e-6;
epife_read_fixed = 22;
Iref_fixed       = 1e-7;
VFB_read_fixed   = -0.5;
VG_name_fixed    = 'VGm0p8_to_2p2';

%% =========================
% 2) 这轮只扫甜区附近
%% =========================
amp_list      = [1.70 1.75 1.80];
pw_list       = [0.03e-6 0.05e-6];
Ea_mean_list  = [0.95 1.00];
Ea_sigma_list = [0.18 0.20];

%% =========================
% 3) 固定 pulse / FE / MOS 参数
%% =========================
delay   = 1e-6;       % 固定，不扫
transit = 0.02e-6;
step    = 0.01e-6;

% FE 参数
Pr    = 25;
tauo  = 1e-9;
alpha = 2;
bet   = 2;
epife_write = 28;

% 快速筛选
Ndom = 3000;

% 读出点：对数采样
Nread = [1 3 10 30 100 300 1000 3000].';

% MOS / stack 参数
til = 1e-7;
miu = 50;
Na  = 3e17;
T0  = 300;
W   = 1;
L   = 1;
VD  = 0.05;
VS  = 0;

% 固定随机样本，保证不同 case 可比
z_base = randn(Ndom,1);

%% =========================
% 4) 评分规则
%    目标趋势：
%      A10 小
%      10->100 增长更明显
%      100->1000 继续涨
%      1000->3000 开始变缓
%% =========================
A10_soft_max  = 0.05;   % 前10脉冲尽量不太大
A100_min      = 0.05;   % 到100脉冲开始抬头
A3000_min     = 0.12;   % 到3000脉冲至少可见

%% =========================
% 5) 生成 case 列表
%% =========================
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

                cfg.tfe         = tfe_fixed;
                cfg.epife_read  = epife_read_fixed;
                cfg.Iref        = Iref_fixed;
                cfg.VFB_read    = VFB_read_fixed;
                cfg.VG_name     = VG_name_fixed;

                cfg.case_id     = make_case_id(cfg);
                case_list{end+1,1} = cfg;
            end
        end
    end
end

nCase = numel(case_list);
fprintf('Total cases = %d\n', nCase);

%% =========================
% 6) 预分配结果
%% =========================
results = repmat(init_result_struct(), nCase, 1);

%% =========================
% 7) 主循环
%% =========================
for ic = 1:nCase
    cfg = case_list{ic};

    fprintf('\n==================================================\n');
    fprintf('[%d / %d] %s\n', ic, nCase, cfg.case_id);
    fprintf('==================================================\n');

    try
        % 7.1 生成当前 case 的 r_Ea
        r_Ea = cfg.Ea_mean + cfg.Ea_sigma * z_base;
        r_Ea(r_Ea <= 0) = cfg.Ea_mean;

        Weight  = ones(Ndom,1) / Ndom;
        St_init = -ones(Ndom,1);
        r_voff  = zeros(Ndom,1);

        % 7.2 跑这一组 Nread
        [dV_vec, A_vec, valid_vec, P_vec, Vth_vec, elapsed_sum] = ...
            run_case_for_Nread(cfg, Nread, ...
                               step, delay, transit, ...
                               St_init, Weight, r_Ea, r_voff, ...
                               Pr, tauo, alpha, bet, epife_write, ...
                               til, miu, Na, T0, W, L, VD, VS);

        % 关键点
        dV1    = dV_vec(1);
        dV3    = dV_vec(2);
        dV10   = dV_vec(3);
        dV30   = dV_vec(4);
        dV100  = dV_vec(5);
        dV300  = dV_vec(6);
        dV1000 = dV_vec(7);
        dV3000 = dV_vec(8);

        A1    = A_vec(1);
        A3    = A_vec(2);
        A10   = A_vec(3);
        A30   = A_vec(4);
        A100  = A_vec(5);
        A300  = A_vec(6);
        A1000 = A_vec(7);
        A3000 = A_vec(8);

        is_invalid = any(~isfinite(A_vec));
        is_mono    = (~is_invalid) && all(diff(A_vec) >= -1e-6);

        g1_1_10      = A10   - A1;
        g2_10_100    = A100  - A10;
        g3_100_1000  = A1000 - A100;
        g4_1000_3000 = A3000 - A1000;

        ratio1   = A1   / max(A3000, 1e-15);
        ratio10  = A10  / max(A3000, 1e-15);
        ratio100 = A100 / max(A3000, 1e-15);

        % 目标趋势判断
        trend_good = false;
        if is_mono
            cond1 = (A10 <= A10_soft_max);
            cond2 = (A100 >= A100_min);
            cond3 = (g2_10_100 > g1_1_10);
            cond4 = (g3_100_1000 > 0);
            cond5 = (g4_1000_3000 <= g3_100_1000 + 1e-6);
            cond6 = (A3000 >= A3000_min);
            trend_good = cond1 && cond2 && cond3 && cond4 && cond5 && cond6;
        end

        % 1e6 脉冲保守外推
        tail_per_dec = max(g4_1000_3000, 0) / max(log10(3000) - log10(1000), 1e-12);
        A1e6_pred = A3000 + 0.8 * tail_per_dec * (6 - log10(3000));

        score = compute_score( ...
            A1, A10, A100, A1000, A3000, ...
            ratio1, ratio10, ratio100, ...
            g1_1_10, g2_10_100, g3_100_1000, g4_1000_3000, ...
            is_mono, trend_good, A1e6_pred);

        % 保存
        results(ic).case_id      = string(cfg.case_id);
        results(ic).amp          = cfg.amp;
        results(ic).pw           = cfg.pw;
        results(ic).Ea_mean      = cfg.Ea_mean;
        results(ic).Ea_sigma     = cfg.Ea_sigma;

        results(ic).tfe          = cfg.tfe;
        results(ic).epife_read   = cfg.epife_read;
        results(ic).Iref         = cfg.Iref;
        results(ic).VFB_read     = cfg.VFB_read;
        results(ic).VG_name      = string(cfg.VG_name);

        results(ic).dV1          = dV1;
        results(ic).dV3          = dV3;
        results(ic).dV10         = dV10;
        results(ic).dV30         = dV30;
        results(ic).dV100        = dV100;
        results(ic).dV300        = dV300;
        results(ic).dV1000       = dV1000;
        results(ic).dV3000       = dV3000;

        results(ic).A1           = A1;
        results(ic).A3           = A3;
        results(ic).A10          = A10;
        results(ic).A30          = A30;
        results(ic).A100         = A100;
        results(ic).A300         = A300;
        results(ic).A1000        = A1000;
        results(ic).A3000        = A3000;

        results(ic).ratio1       = ratio1;
        results(ic).ratio10      = ratio10;
        results(ic).ratio100     = ratio100;

        results(ic).g1_1_10      = g1_1_10;
        results(ic).g2_10_100    = g2_10_100;
        results(ic).g3_100_1000  = g3_100_1000;
        results(ic).g4_1000_3000 = g4_1000_3000;

        results(ic).A1e6_pred    = A1e6_pred;
        results(ic).is_mono      = is_mono;
        results(ic).trend_good   = trend_good;
        results(ic).score        = score;
        results(ic).elapsed_s    = elapsed_sum;

        results(ic).Nread_str    = vec2str(Nread);
        results(ic).Avec_str     = vec2str(A_vec);

    catch ME
        warning('Case failed: %s', cfg.case_id);
        warning('%s', ME.message);

        results(ic).case_id   = string(cfg.case_id);
        results(ic).amp       = cfg.amp;
        results(ic).pw        = cfg.pw;
        results(ic).Ea_mean   = cfg.Ea_mean;
        results(ic).Ea_sigma  = cfg.Ea_sigma;
        results(ic).score     = -1e9;
    end
end

%% =========================
% 8) 生成表
%% =========================
T = struct2table(results);
writetable(T, fullfile(outdir, 'all_cases.csv'));

T_valid = T(isfinite(T.score), :);

T_score = sortrows(T_valid, 'score', 'descend');
writetable(T_score(1:min(20,height(T_score)), :), fullfile(outdir, 'top_by_score.csv'));

T_A3000 = sortrows(T_valid, 'A3000', 'descend');
writetable(T_A3000(1:min(20,height(T_A3000)), :), fullfile(outdir, 'top_by_A3000.csv'));

disp(' ');
disp('================ Top-20 by score ================');
disp(T_score(1:min(20,height(T_score)), :));

disp(' ');
disp('================ Top-20 by A3000 ================');
disp(T_A3000(1:min(20,height(T_A3000)), :));

%% =========================
% 9) 画 top5 曲线
%% =========================
nplot = min(5, height(T_score));

for i = 1:nplot
    row = T_score(i,:);
    Nvec = str2num(row.Nread_str{1}); %#ok<ST2NM>
    Avec = str2num(row.Avec_str{1}); %#ok<ST2NM>

    figure('Visible', 'off');
    semilogx(Nvec, Avec, 'o-', 'LineWidth', 1.5);
    xlabel('Pulse count N');
    ylabel('|DeltaVth| (V)');
    title(sprintf('Top-%d: %s', i, row.case_id{1}), 'Interpreter', 'none');
    grid on;
    saveas(gcf, fullfile(outdir, sprintf('top%d_curve.png', i)));
    close(gcf);
end

%% =========================
% 10) 保存工作区
%% =========================
save(fullfile(outdir, 'workspace_results.mat'), ...
    'T', 'T_score', 'T_A3000', 'results', 'case_list', ...
    'amp_list', 'pw_list', 'Ea_mean_list', 'Ea_sigma_list', ...
    'tfe_fixed', 'epife_read_fixed', 'Iref_fixed', 'VFB_read_fixed', 'VG_name_fixed', ...
    'Nread', 'Ndom', ...
    'delay', 'transit', 'step', ...
    'Pr', 'tauo', 'alpha', 'bet', 'epife_write', ...
    'til', 'miu', 'Na', 'T0', 'W', 'L', 'VD', 'VS');

fprintf('\nAll done. Results saved in folder: %s\n', outdir);
diary off;