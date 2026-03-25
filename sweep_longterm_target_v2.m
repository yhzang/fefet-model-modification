%% =========================================================
% sweep_longterm_target_v2.m
%
% 目标：
%   找到满足下面特征的参数组合：
%   1) 前几十个脉冲（尤其前 30 个）不要快速饱和
%   2) 在 100~1000 脉冲区间仍保持增长
%   3) 依据 log(N) 趋势外推后，在 1e6 脉冲时
%      |DeltaVth| 约 >= 0.5 V
%
% 依赖：
%   wfdef_acc_fix.m
%   get_FE_state.m
%   get_ID.m
%
% 输出目录：
%   longterm_target_v2_out/
%       all_cases.csv
%       top20_by_score.csv
%       progressive_good.csv
%       workspace_results.mat
%       多张诊断图
%
% 说明：
%   1) 本脚本默认实际仿真到 1000 pulses
%   2) 对 1e6 pulses 的结果使用 log10(N) 后段拟合外推
%   3) 这更适合“筛选候选参数”，不是最终定标结论
%% =========================================================

clear;
clc;
close all;

rng(0);

%% =========================
% 0) 输出目录
% =========================
outdir = 'longterm_target_v2_out';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

diary(fullfile(outdir, 'run_log.txt'));
fprintf('=== sweep_longterm_target_v2 ===\n');

%% =========================
% 1) 扫描空间
%    这一轮重点：
%    - 适度提高 amp
%    - 保留较短 pw，避免前期太猛
%    - 提高 Ea_mean / Ea_sigma，拉长累积窗口
%    - 适当扩展 tfe / epife_read
% =========================
amp_list        = [1.8 2.0 2.2 2.4];
pw_list         = [0.01e-6 0.02e-6 0.05e-6 0.1e-6];
Ea_mean_list    = [1.05 1.10 1.15 1.20];
Ea_sigma_list   = [0.10 0.15 0.20 0.25];

tfe_list        = [1.0e-6 1.5e-6 2.0e-6];
epife_read_list = [20 24 28 32];

Iref_list       = [1e-7];
VFB_read_list   = [-0.5];
VG_name_list    = {'VGm0p8_to_2p2'};

%% =========================
% 2) 固定 pulse / FE / MOS 参数
% =========================
VDD = 5;

delay   = 1e-6;
transit = 0.02e-6;
cycle   = 1000;          % 实际仿真到 1000 pulses
step    = 0.01e-6;

% FE 参数
Pr    = 25;
tauo  = 1e-9;
alpha = 2;
bet   = 2;

% 写入主干 FE eps
epife_write = 28;

Ndom = 10000;

Weight  = ones(Ndom,1) / Ndom;
St_init = -ones(Ndom,1);

% 是否全 case 共用同一组 domain sample
use_global_fixed_domain_sample = false;

if use_global_fixed_domain_sample
    Ea_mean0  = 1.10;
    Ea_sigma0 = 0.20;
    r_Ea_base = Ea_mean0 + Ea_sigma0 * randn(Ndom,1);
    r_Ea_base(r_Ea_base <= 0) = Ea_mean0;
end

% 先不引入 voff 演化
r_voff = zeros(Ndom,1);

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
% 3) 读出点
%    更贴近你的目标：
%    看前几十个、几百个、1000 个的趋势
%% =========================
Nread = [1 10 30 100 300 1000].';

%% =========================
% 4) 判据与目标
% =========================
eps_num = 1e-15;
Vmin_eval = 5e-3;     % 小于 5 mV 视作过弱，不做趋势判断

% 你的真实目标
target_dV_1e6 = 0.5;  % 1e6 pulses 时希望有约 0.5 V

% “前几十个不要直接饱和” 的工程判据
ratio30_max          = 0.35;   % |dV30| / |dV1000| 不要太大
ratio100_max         = 0.65;   % |dV100| / |dV1000| 不要太大
late_gain_min        = 0.03;   % |dV1000|-|dV100| 至少有 30mV
late_early_ratio_min = 1.00;   % 后段增长 >= 前段增长

% 进一步偏向“长期可累积”
pred_1e6_min_for_good = 0.40;  % 预测 1e6 先过 0.4V 就认为值得继续看

%% =========================
% 5) 生成 case 列表
% =========================
case_list = {};
icase = 0;

for ia = 1:numel(amp_list)
    for ipw = 1:numel(pw_list)
        for iem = 1:numel(Ea_mean_list)
            for ies = 1:numel(Ea_sigma_list)
                for itfe = 1:numel(tfe_list)
                    for iepi = 1:numel(epife_read_list)
                        for iir = 1:numel(Iref_list)
                            for ivfb = 1:numel(VFB_read_list)
                                for ivg = 1:numel(VG_name_list)
                                    icase = icase + 1;

                                    cfg = struct();
                                    cfg.case_idx    = icase;
                                    cfg.amp         = amp_list(ia);
                                    cfg.pw          = pw_list(ipw);
                                    cfg.Ea_mean     = Ea_mean_list(iem);
                                    cfg.Ea_sigma    = Ea_sigma_list(ies);

                                    cfg.tfe         = tfe_list(itfe);
                                    cfg.epife_read  = epife_read_list(iepi);
                                    cfg.Iref        = Iref_list(iir);
                                    cfg.VFB_read    = VFB_read_list(ivfb);
                                    cfg.VG_name     = VG_name_list{ivg};

                                    cfg.case_id = make_case_id(cfg);

                                    case_list{end+1,1} = cfg;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

nCase = numel(case_list);
fprintf('Total cases = %d\n', nCase);

%% =========================
% 6) 预分配结果
% =========================
results = repmat(init_result_struct(), nCase, 1);

%% =========================
% 7) 主循环
% =========================
for ic = 1:nCase
    cfg = case_list{ic};

    fprintf('\n==================================================\n');
    fprintf('[%d / %d] %s\n', ic, nCase, cfg.case_id);
    fprintf('==================================================\n');

    try
        %% 7.1 生成当前 case 的 r_Ea
        if use_global_fixed_domain_sample
            mu0 = mean(r_Ea_base);
            sg0 = std(r_Ea_base);
            if sg0 < 1e-12
                r_Ea = cfg.Ea_mean * ones(Ndom,1);
            else
                r_Ea = (r_Ea_base - mu0) / sg0 * cfg.Ea_sigma + cfg.Ea_mean;
            end
            r_Ea(r_Ea <= 0) = cfg.Ea_mean;
        else
            r_Ea = cfg.Ea_mean + cfg.Ea_sigma * randn(Ndom,1);
            r_Ea(r_Ea <= 0) = cfg.Ea_mean;
        end

        %% 7.2 波形
        [time, volt, index, index_pre, index_end] = ...
            wfdef_acc_fix(cfg.amp, cfg.pw, step, delay, transit, cycle);

        index     = index(:);
        index_pre = index_pre(:);
        index_end = index_end(:);

        %% 7.3 FE 状态推进
        tic;
        [vfev, Stsum] = get_FE_state(time, volt, St_init, Weight, r_Ea, r_voff, ...
                                     Pr, tauo, alpha, bet, epife_write, Ndom);
        telapsed = toc;
        fprintf('get_FE_state done, elapsed = %.3f s\n', telapsed);

        %% 7.4 单条轨迹
        if isvector(Stsum)
            st_trace = Stsum(:);
        else
            st_trace = Stsum(:,1);
        end

        %% 7.5 读出 VG
        VG = make_VG_from_name(cfg.VG_name);

        %% 7.6 初始参考点
        idx0 = index_pre(1);
        P0 = st_trace(idx0);

        [Vth0, ~, ~, ~, valid0] = local_extract_vth( ...
            P0, cfg.epife_read, cfg.tfe, til, miu, Na, T0, W, L, ...
            VG, cfg.VFB_read, VD, VS, cfg.Iref);

        %% 7.7 提取 Nread 各点
        nR = numel(Nread);

        Vth_vec      = nan(nR,1);
        DeltaVth_vec = nan(nR,1);
        Pread_vec    = nan(nR,1);
        IDmin_vec    = nan(nR,1);
        IDmax_vec    = nan(nR,1);
        valid_vec    = false(nR,1);

        for kk = 1:nR
            n = Nread(kk);
            idxr = index(n);
            Pread_vec(kk) = st_trace(idxr);

            [Vth_vec(kk), ~, IDmin_vec(kk), IDmax_vec(kk), valid_vec(kk)] = ...
                local_extract_vth( ...
                    Pread_vec(kk), cfg.epife_read, cfg.tfe, til, miu, Na, T0, W, L, ...
                    VG, cfg.VFB_read, VD, VS, cfg.Iref);

            if valid0 && valid_vec(kk)
                DeltaVth_vec(kk) = Vth_vec(kk) - Vth0;
            else
                DeltaVth_vec(kk) = NaN;
            end
        end

        % 若 Vth0 无效，则退化到首个有效点为参考
        if ~valid0
            first_valid = find(valid_vec, 1, 'first');
            if ~isempty(first_valid)
                Vth_ref = Vth_vec(first_valid);
                DeltaVth_vec = Vth_vec - Vth_ref;
            end
        end

        %% 7.8 读出关键点
        dV1    = DeltaVth_vec(1);
        dV10   = DeltaVth_vec(2);
        dV30   = DeltaVth_vec(3);
        dV100  = DeltaVth_vec(4);
        dV300  = DeltaVth_vec(5);
        dV1000 = DeltaVth_vec(6);

        A1    = abs(dV1);
        A10   = abs(dV10);
        A30   = abs(dV30);
        A100  = abs(dV100);
        A300  = abs(dV300);
        A1000 = abs(dV1000);

        %% 7.9 增量指标
        inc_1_30      = A30   - A1;
        inc_30_100    = A100  - A30;
        inc_100_300   = A300  - A100;
        inc_300_1000  = A1000 - A300;
        late_gain     = A1000 - A100;

        ratio30  = NaN;
        ratio100 = NaN;
        late_early_ratio = NaN;

        is_invalid = any(~isfinite([dV1 dV10 dV30 dV100 dV300 dV1000])) || ...
                     any(~isfinite([A1 A10 A30 A100 A300 A1000]));

        is_weak_all = (~is_invalid) && (A1000 < Vmin_eval);

        if (~is_invalid) && (~is_weak_all)
            ratio30  = A30  / max(A1000, eps_num);
            ratio100 = A100 / max(A1000, eps_num);

            early_gain = max(A100 - A1, 0);
            late_gain_tmp = max(A1000 - A100, 0);

            if early_gain > eps_num
                late_early_ratio = late_gain_tmp / early_gain;
            else
                if late_gain_tmp > eps_num
                    late_early_ratio = inf;
                else
                    late_early_ratio = 0;
                end
            end
        end

        %% 7.10 “前几十个不要直接饱和”
        is_early_saturate = false;
        if (~is_invalid) && (~is_weak_all)
            cond_sat1 = ratio30  > ratio30_max;
            cond_sat2 = ratio100 > ratio100_max;
            cond_sat3 = late_gain < late_gain_min;
            cond_sat4 = late_early_ratio < late_early_ratio_min;

            % 如果 1->30 已经很大，而后面涨不动，也视为早饱和
            cond_sat5 = (A30 > A10) && (A100 <= A30 + 0.01) && (A1000 <= A100 + 0.02);

            is_early_saturate = cond_sat1 || cond_sat2 || cond_sat3 || cond_sat4 || cond_sat5;
        end

        %% 7.11 拟合 log10(N) 趋势，外推到 1e6
        % 只用后段 [100, 300, 1000] 做拟合，更贴近长期趋势
        dV_1e6_pred = NaN;
        fit_slope = NaN;
        fit_intercept = NaN;
        fit_rmse = NaN;
        fit_npts = 0;

        [dV_1e6_pred, fit_slope, fit_intercept, fit_rmse, fit_npts] = ...
            predict_dv_at_1e6(Nread, DeltaVth_vec);

        A1e6_pred = abs(dV_1e6_pred);

        %% 7.12 是否值得继续
        is_progressive_good = false;
        if (~is_invalid) && (~is_weak_all) && (~is_early_saturate)
            cond_g1 = (A1000 > A100) && (A100 > A30);
            cond_g2 = ratio30 <= ratio30_max;
            cond_g3 = ratio100 <= ratio100_max;
            cond_g4 = late_gain >= late_gain_min;
            cond_g5 = late_early_ratio >= late_early_ratio_min;
            cond_g6 = isfinite(A1e6_pred) && (A1e6_pred >= pred_1e6_min_for_good);

            is_progressive_good = cond_g1 && cond_g2 && cond_g3 && cond_g4 && cond_g5 && cond_g6;
        end

        %% 7.13 评分
        target_score = compute_target_score( ...
            A1, A10, A30, A100, A300, A1000, ...
            ratio30, ratio100, late_gain, late_early_ratio, ...
            A1e6_pred, target_dV_1e6, ...
            is_invalid, is_weak_all, is_early_saturate, fit_rmse);

        %% 7.14 保存结果
        results(ic).case_id           = string(cfg.case_id);
        results(ic).amp               = cfg.amp;
        results(ic).pw                = cfg.pw;
        results(ic).Ea_mean           = cfg.Ea_mean;
        results(ic).Ea_sigma          = cfg.Ea_sigma;
        results(ic).tfe               = cfg.tfe;
        results(ic).epife_read        = cfg.epife_read;
        results(ic).Iref              = cfg.Iref;
        results(ic).VFB_read          = cfg.VFB_read;
        results(ic).VG_name           = string(cfg.VG_name);

        results(ic).dV1               = dV1;
        results(ic).dV10              = dV10;
        results(ic).dV30              = dV30;
        results(ic).dV100             = dV100;
        results(ic).dV300             = dV300;
        results(ic).dV1000            = dV1000;

        results(ic).A1                = A1;
        results(ic).A10               = A10;
        results(ic).A30               = A30;
        results(ic).A100              = A100;
        results(ic).A300              = A300;
        results(ic).A1000             = A1000;

        results(ic).inc_1_30          = inc_1_30;
        results(ic).inc_30_100        = inc_30_100;
        results(ic).inc_100_300       = inc_100_300;
        results(ic).inc_300_1000      = inc_300_1000;
        results(ic).late_gain         = late_gain;

        results(ic).ratio30           = ratio30;
        results(ic).ratio100          = ratio100;
        results(ic).late_early_ratio  = late_early_ratio;

        results(ic).dV_1e6_pred       = dV_1e6_pred;
        results(ic).A1e6_pred         = A1e6_pred;
        results(ic).fit_slope         = fit_slope;
        results(ic).fit_intercept     = fit_intercept;
        results(ic).fit_rmse          = fit_rmse;
        results(ic).fit_npts          = fit_npts;

        results(ic).target_score      = target_score;

        results(ic).is_invalid        = is_invalid;
        results(ic).is_weak_all       = is_weak_all;
        results(ic).is_early_saturate = is_early_saturate;
        results(ic).is_progressive_good = is_progressive_good;

        results(ic).valid0            = valid0;
        results(ic).valid1            = valid_vec(1);
        results(ic).valid10           = valid_vec(2);
        results(ic).valid30           = valid_vec(3);
        results(ic).valid100          = valid_vec(4);
        results(ic).valid300          = valid_vec(5);
        results(ic).valid1000         = valid_vec(6);

        results(ic).P0                = P0;
        results(ic).P1                = Pread_vec(1);
        results(ic).P10               = Pread_vec(2);
        results(ic).P30               = Pread_vec(3);
        results(ic).P100              = Pread_vec(4);
        results(ic).P300              = Pread_vec(5);
        results(ic).P1000             = Pread_vec(6);

        results(ic).Vth0              = Vth0;
        results(ic).Vth1              = Vth_vec(1);
        results(ic).Vth10             = Vth_vec(2);
        results(ic).Vth30             = Vth_vec(3);
        results(ic).Vth100            = Vth_vec(4);
        results(ic).Vth300            = Vth_vec(5);
        results(ic).Vth1000           = Vth_vec(6);

        results(ic).elapsed_s         = telapsed;

    catch ME
        warning('Case failed: %s', cfg.case_id);
        warning('%s', ME.message);

        results(ic).case_id           = string(cfg.case_id);
        results(ic).amp               = cfg.amp;
        results(ic).pw                = cfg.pw;
        results(ic).Ea_mean           = cfg.Ea_mean;
        results(ic).Ea_sigma          = cfg.Ea_sigma;
        results(ic).tfe               = cfg.tfe;
        results(ic).epife_read        = cfg.epife_read;
        results(ic).Iref              = cfg.Iref;
        results(ic).VFB_read          = cfg.VFB_read;
        results(ic).VG_name           = string(cfg.VG_name);

        results(ic).target_score      = -1e9;
        results(ic).is_invalid        = true;
    end
end

%% =========================
% 8) 生成表
% =========================
T = struct2table(results);

writetable(T, fullfile(outdir, 'all_cases.csv'));

T_valid = T(~T.is_invalid, :);
T_score = sortrows(T_valid, 'target_score', 'descend');
writetable(T_score(1:min(20,height(T_score)), :), ...
    fullfile(outdir, 'top20_by_score.csv'));

T_good = T(T.is_progressive_good, :);
if ~isempty(T_good)
    T_good = sortrows(T_good, 'target_score', 'descend');
end
writetable(T_good, fullfile(outdir, 'progressive_good.csv'));

disp(' ');
disp('================ Top-20 by target_score ================');
disp(T_score(1:min(20,height(T_score)), :));

%% =========================
% 9) 画图
% =========================
case_idx = (1:height(T)).';

idx_invalid = T.is_invalid;
idx_weak    = T.is_weak_all & ~idx_invalid;
idx_sat     = T.is_early_saturate & ~idx_invalid & ~idx_weak;
idx_good    = T.is_progressive_good & ~idx_invalid & ~idx_weak;
idx_other   = ~idx_invalid & ~idx_weak & ~idx_sat & ~idx_good;

% 图1：A30 / A100 / A1000 / A1e6_pred
figure('Name', 'Figure 1');
plot(case_idx, T.A30, 'o-', 'DisplayName', '|dV30|'); hold on;
plot(case_idx, T.A100, 's-', 'DisplayName', '|dV100|');
plot(case_idx, T.A1000, 'd-', 'DisplayName', '|dV1000|');
plot(case_idx, T.A1e6_pred, '^-', 'DisplayName', '|dV1e6\_pred|');
yline(target_dV_1e6, 'k--', 'DisplayName', 'target 0.5V');
xlabel('Case index');
ylabel('Magnitude (V)');
title('|dV30|, |dV100|, |dV1000|, |dV1e6\_pred|');
grid on;
legend('Location', 'best');
saveas(gcf, fullfile(outdir, 'fig1_magnitudes_and_pred.png'));

% 图2：ratio30 / ratio100
figure('Name', 'Figure 2');
valid_ratio = ~idx_invalid & ~idx_weak;
plot(case_idx(valid_ratio), T.ratio30(valid_ratio), 'o-', ...
    'DisplayName', 'ratio30 = |dV30|/|dV1000|'); hold on;
plot(case_idx(valid_ratio), T.ratio100(valid_ratio), 's-', ...
    'DisplayName', 'ratio100 = |dV100|/|dV1000|');
yline(ratio30_max, 'k--', 'DisplayName', 'ratio30 max');
yline(ratio100_max, 'r--', 'DisplayName', 'ratio100 max');
xlabel('Case index');
ylabel('Ratio');
title('Early saturation indicators');
grid on;
legend('Location', 'best');
saveas(gcf, fullfile(outdir, 'fig2_ratios.png'));

% 图3：A30 vs A1000
figure('Name', 'Figure 3');
plot(T.A30(idx_weak),  T.A1000(idx_weak),  'o', 'DisplayName', 'weak_all'); hold on;
plot(T.A30(idx_sat),   T.A1000(idx_sat),   's', 'DisplayName', 'early_saturate');
plot(T.A30(idx_good),  T.A1000(idx_good),  'd', 'DisplayName', 'progressive_good');
plot(T.A30(idx_other), T.A1000(idx_other), '^', 'DisplayName', 'other_valid');
xymax = max([T.A30; T.A1000; 1e-3]) * 1.05;
plot([0 xymax], [0 xymax], 'k--', 'DisplayName', 'y=x');
xlabel('|dV30| (V)');
ylabel('|dV1000| (V)');
title('|dV30| vs |dV1000|');
grid on;
legend('Location', 'best');
saveas(gcf, fullfile(outdir, 'fig3_dV30_vs_dV1000.png'));

% 图4：A1000 vs A1e6_pred
figure('Name', 'Figure 4');
plot(T.A1000(idx_weak),  T.A1e6_pred(idx_weak),  'o', 'DisplayName', 'weak_all'); hold on;
plot(T.A1000(idx_sat),   T.A1e6_pred(idx_sat),   's', 'DisplayName', 'early_saturate');
plot(T.A1000(idx_good),  T.A1e6_pred(idx_good),  'd', 'DisplayName', 'progressive_good');
plot(T.A1000(idx_other), T.A1e6_pred(idx_other), '^', 'DisplayName', 'other_valid');
xymax = max([T.A1000; T.A1e6_pred; 1e-3]) * 1.05;
plot([0 xymax], [0 xymax], 'k--', 'DisplayName', 'y=x');
yline(target_dV_1e6, 'r--', 'DisplayName', '0.5V target');
xlabel('|dV1000| (V)');
ylabel('|dV1e6\_pred| (V)');
title('|dV1000| vs |dV1e6\_pred|');
grid on;
legend('Location', 'best');
saveas(gcf, fullfile(outdir, 'fig4_dV1000_vs_dV1e6pred.png'));

% 图5：各段增量
figure('Name', 'Figure 5');
valid_inc = ~idx_invalid;
plot(case_idx(valid_inc), T.inc_1_30(valid_inc), 'o-', ...
    'DisplayName', 'inc_{1->30}'); hold on;
plot(case_idx(valid_inc), T.inc_30_100(valid_inc), 's-', ...
    'DisplayName', 'inc_{30->100}');
plot(case_idx(valid_inc), T.inc_100_300(valid_inc), 'd-', ...
    'DisplayName', 'inc_{100->300}');
plot(case_idx(valid_inc), T.inc_300_1000(valid_inc), '^-', ...
    'DisplayName', 'inc_{300->1000}');
yline(0, 'k--', 'DisplayName', 'zero');
xlabel('Case index');
ylabel('Increment (V)');
title('Incremental growth in different pulse windows');
grid on;
legend('Location', 'best');
saveas(gcf, fullfile(outdir, 'fig5_increments.png'));

% 图6：late_early_ratio
figure('Name', 'Figure 6');
valid_le = ~idx_invalid & ~idx_weak & isfinite(T.late_early_ratio);
plot(case_idx(valid_le), T.late_early_ratio(valid_le), 'o-', ...
    'DisplayName', 'late\_early\_ratio'); hold on;
yline(late_early_ratio_min, 'k--', 'DisplayName', '1.0');
xlabel('Case index');
ylabel('late\_early\_ratio');
title('late\_early\_ratio = (|dV1000|-|dV100|)/(|dV100|-|dV1|)');
grid on;
legend('Location', 'best');
saveas(gcf, fullfile(outdir, 'fig6_late_early_ratio.png'));

% 图7：score
figure('Name', 'Figure 7');
valid_score = ~idx_invalid;
plot(case_idx(valid_score), T.target_score(valid_score), 'o-', ...
    'DisplayName', 'target_score');
xlabel('Case index');
ylabel('Score');
title('Target score across cases');
grid on;
legend('Location', 'best');
saveas(gcf, fullfile(outdir, 'fig7_target_score.png'));

%% =========================
% 10) 保存工作区
% =========================
save(fullfile(outdir, 'workspace_results.mat'), ...
     'T', 'T_score', 'T_good', 'results', 'case_list', ...
     'amp_list', 'pw_list', 'Ea_mean_list', 'Ea_sigma_list', ...
     'tfe_list', 'epife_read_list', 'Iref_list', 'VFB_read_list', 'VG_name_list', ...
     'delay', 'transit', 'cycle', 'step', ...
     'Pr', 'tauo', 'alpha', 'bet', 'epife_write', 'Ndom', ...
     'til', 'miu', 'Na', 'T0', 'W', 'L', 'VD', 'VS', ...
     'Nread', ...
     'eps_num', 'Vmin_eval', ...
     'target_dV_1e6', ...
     'ratio30_max', 'ratio100_max', ...
     'late_gain_min', 'late_early_ratio_min', ...
     'pred_1e6_min_for_good');

fprintf('\nAll done. Results saved in folder: %s\n', outdir);
diary off;

%% =========================================================
% local functions
%% =========================================================

function s = init_result_struct()
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
    'dV30', NaN, ...
    'dV100', NaN, ...
    'dV300', NaN, ...
    'dV1000', NaN, ...
    'A1', NaN, ...
    'A10', NaN, ...
    'A30', NaN, ...
    'A100', NaN, ...
    'A300', NaN, ...
    'A1000', NaN, ...
    'inc_1_30', NaN, ...
    'inc_30_100', NaN, ...
    'inc_100_300', NaN, ...
    'inc_300_1000', NaN, ...
    'late_gain', NaN, ...
    'ratio30', NaN, ...
    'ratio100', NaN, ...
    'late_early_ratio', NaN, ...
    'dV_1e6_pred', NaN, ...
    'A1e6_pred', NaN, ...
    'fit_slope', NaN, ...
    'fit_intercept', NaN, ...
    'fit_rmse', NaN, ...
    'fit_npts', NaN, ...
    'target_score', NaN, ...
    'is_invalid', false, ...
    'is_weak_all', false, ...
    'is_early_saturate', false, ...
    'is_progressive_good', false, ...
    'valid0', false, ...
    'valid1', false, ...
    'valid10', false, ...
    'valid30', false, ...
    'valid100', false, ...
    'valid300', false, ...
    'valid1000', false, ...
    'P0', NaN, ...
    'P1', NaN, ...
    'P10', NaN, ...
    'P30', NaN, ...
    'P100', NaN, ...
    'P300', NaN, ...
    'P1000', NaN, ...
    'Vth0', NaN, ...
    'Vth1', NaN, ...
    'Vth10', NaN, ...
    'Vth30', NaN, ...
    'Vth100', NaN, ...
    'Vth300', NaN, ...
    'Vth1000', NaN, ...
    'elapsed_s', NaN ...
    );
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

function [dV_pred, slope, intercept, rmse, npts] = predict_dv_at_1e6(Nread, DeltaVth_vec)
% 用后段 [100, 300, 1000] 的 DeltaVth 对 log10(N) 做线性拟合
% 再预测 N=1e6 时的 DeltaVth
%
% 若后段不可用，则退化到所有有效点中 N>=30 的部分
% 若还是不够，则返回 NaN

dV_pred   = NaN;
slope     = NaN;
intercept = NaN;
rmse      = NaN;
npts      = 0;

Nread = Nread(:);
DeltaVth_vec = DeltaVth_vec(:);

mask = isfinite(Nread) & isfinite(DeltaVth_vec);

% 优先只用后段
mask_late = mask & (Nread >= 100);

if nnz(mask_late) >= 3
    x = log10(Nread(mask_late));
    y = DeltaVth_vec(mask_late);
else
    mask_mid = mask & (Nread >= 30);
    if nnz(mask_mid) >= 3
        x = log10(Nread(mask_mid));
        y = DeltaVth_vec(mask_mid);
    else
        return;
    end
end

p = polyfit(x, y, 1);
slope = p(1);
intercept = p(2);

yfit = polyval(p, x);
rmse = sqrt(mean((y - yfit).^2));
npts = numel(x);

dV_pred = polyval(p, log10(1e6));
end

function score = compute_target_score( ...
    A1, A10, A30, A100, A300, A1000, ...
    ratio30, ratio100, late_gain, late_early_ratio, ...
    A1e6_pred, target_dV_1e6, ...
    is_invalid, is_weak_all, is_early_saturate, fit_rmse)

if is_invalid
    score = -1e9;
    return;
end

if is_weak_all
    score = -1e6 - A1000;
    return;
end

% 越接近目标越好
term_pred_target = 0;
if isfinite(A1e6_pred)
    term_pred_target = 4.0 * min(A1e6_pred, target_dV_1e6) ...
                     + 2.0 * max(A1e6_pred - target_dV_1e6, 0);
else
    term_pred_target = -1.0;
end

% 当前 1000 脉冲也要有一定量级
term_now = 2.0 * A1000 + 1.0 * A300;

% 不希望前 30 / 100 太大
term_early_penalty = 0;
if isfinite(ratio30)
    term_early_penalty = term_early_penalty - 2.0 * ratio30;
end
if isfinite(ratio100)
    term_early_penalty = term_early_penalty - 2.5 * ratio100;
end

% 希望后段还有增量
term_late_gain = 3.0 * max(late_gain, 0);

% 希望后段增速不差
term_late_ratio = 0;
if isfinite(late_early_ratio)
    term_late_ratio = 0.08 * min(late_early_ratio, 8);
end

% 拟合太差，说明外推不可信，要扣分
term_fit = 0;
if isfinite(fit_rmse)
    term_fit = -2.0 * fit_rmse;
end

term_sat = 0;
if is_early_saturate
    term_sat = -2.0;
end

score = term_pred_target + term_now + term_early_penalty + ...
        term_late_gain + term_late_ratio + term_fit + term_sat;
end