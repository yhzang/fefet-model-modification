%% =========================================================
% sweep_progressive_window_for_longterm.m
%
% 目的：
%   在当前 Scalable-FeFET / 多畴 FE 主干上，扫描一批参数，
%   找到更适合 long-term 渐进累积的 candidate。
%
% 依赖：
%   wfdef_acc_fix.m
%   get_FE_state.m
%   get_ID.m
%
% 输出：
%   progressive_window_out/
%       all_cases.csv
%       top20_by_gradual.csv
%       top20_by_progress_score.csv
%       progressive_only.csv
%       workspace_results.mat
%       多张诊断图
%
% 关键修正：
%   1) weak_all 不再参与 ratio 计算，避免 ratio 爆到 1e12
%   2) 引入 inc_1_10, inc_10_100, late_early_ratio
%   3) progressive / early_saturate 分类更贴近 long-term 目标
%
% 备注：
%   这版是“完整脚本”，不是占位骨架。
%   它直接调用你现有的 wfdef_acc_fix / get_FE_state / get_ID。
%% =========================================================

clear;
clc;
close all;

rng(0);

%% =========================
% 0) 输出目录
% =========================
outdir = 'progressive_window_out';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

diary(fullfile(outdir, 'run_log.txt'));
fprintf('=== sweep_progressive_window_for_longterm ===\n');

%% =========================
% 1) 扫描空间
%    这里按你当前阶段目标，优先找 long-term candidate
%    所以不要再把范围放太激进
% =========================
amp_list        = [1.6 1.8 2.0 2.2];
pw_list         = [0.02e-6 0.05e-6 0.1e-6 0.2e-6 0.5e-6];
Ea_mean_list    = [0.95 1.00 1.05 1.10];
Ea_sigma_list   = [0.10 0.15 0.20 0.25];

% 读出映射参数
tfe_list        = [0.8e-6 1.0e-6 1.5e-6];
epife_read_list = [20 24 28];
Iref_list       = [1e-7];
VFB_read_list   = [-0.5];
VG_name_list    = {'VGm0p8_to_2p2'};

%% =========================
% 2) 固定 pulse / FE / MOS 参数
%    这里保持和你当前调试链一致
% =========================
VDD = 5;

delay   = 1e-6;
transit = 0.02e-6;
cycle   = 100;
step    = 0.01e-6;

% FE 参数
Pr    = 25;
tauo  = 1e-9;
alpha = 2;
bet   = 2;

% 默认 FE eps（写入主干）
epife_write = 28;

Ndom = 10000;

% 固定随机样本：比较不同 case 时，不改 domain realization
Weight  = ones(Ndom,1) / Ndom;
St_init = -ones(Ndom,1);

% 注意：这里每个 case 都重新生成 r_Ea，
% 但同一 case 内不同 Nread 来自同一条轨迹
% 若你想“所有 case 完全共用一套随机样本”，可改成外面固定 r_Ea_base
use_global_fixed_domain_sample = false;

if use_global_fixed_domain_sample
    Ea_mean0  = 1.0;
    Ea_sigma0 = 0.2;
    r_Ea_base = Ea_mean0 + Ea_sigma0 * randn(Ndom,1);
    r_Ea_base(r_Ea_base <= 0) = Ea_mean0;
end

% baseline：先关掉 r_voff，避免人为引入 gap 演化
r_voff = zeros(Ndom,1);

% MOS / stack 参数（与现有读出链保持一致）
til = 1e-7;
miu = 50;
Na  = 3e17;
T   = 300;
W   = 1;
L   = 1;
VD  = 0.05;
VS  = 0;

%% =========================
% 3) 读出点
% =========================
Nread = [1 10 100].';

%% =========================
% 4) 分类阈值
% =========================
eps_num = 1e-15;

% 小于这个 |dV100| 就不讨论 ratio，避免奇异
Vmin_eval = 5e-3;    % 5 mV

% early saturation 判据
ratio10_sat_th   = 0.65;
late_vs_early_th = 0.70;

% progressive 判据
ratio1_prog_max      = 0.20;
ratio10_prog_max     = 0.60;
late_early_ratio_min = 1.00;
gradual_min          = 0.03;   % 30 mV

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
            % 若固定样本，则只做线性重标定，不重新抽样
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

        %% 7.2 生成波形
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

        %% 7.4 取单条轨迹
        if isvector(Stsum)
            st_trace = Stsum(:);
        else
            st_trace = Stsum(:,1);
        end

        %% 7.5 选择读出 VG
        VG = make_VG_from_name(cfg.VG_name);

        %% 7.6 初始参考点
        idx0 = index_pre(1);
        P0 = st_trace(idx0);

        [Vth0, ~, ~, ~, valid0] = local_extract_vth( ...
            P0, cfg.epife_read, cfg.tfe, til, miu, Na, T, W, L, ...
            VG, cfg.VFB_read, VD, VS, cfg.Iref);

        %% 7.7 提取 N=1,10,100
        Vth_vec      = nan(numel(Nread),1);
        DeltaVth_vec = nan(numel(Nread),1);
        Pread_vec    = nan(numel(Nread),1);
        IDmin_vec    = nan(numel(Nread),1);
        IDmax_vec    = nan(numel(Nread),1);
        valid_vec    = false(numel(Nread),1);

        for kk = 1:numel(Nread)
            n = Nread(kk);
            idxr = index(n);
            Pread_vec(kk) = st_trace(idxr);

            [Vth_vec(kk), ~, IDmin_vec(kk), IDmax_vec(kk), valid_vec(kk)] = ...
                local_extract_vth( ...
                    Pread_vec(kk), cfg.epife_read, cfg.tfe, til, miu, Na, T, W, L, ...
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

        dV1   = DeltaVth_vec(1);
        dV10  = DeltaVth_vec(2);
        dV100 = DeltaVth_vec(3);

        A1   = abs(dV1);
        A10  = abs(dV10);
        A100 = abs(dV100);

        gradual    = A100 - A1;
        inc_1_10   = A10 - A1;
        inc_10_100 = A100 - A10;

        %% 7.8 分类
        is_invalid = any(~isfinite([dV1 dV10 dV100])) || ...
                     any(~isfinite([A1 A10 A100]));

        is_weak_all = (~is_invalid) && (A100 < Vmin_eval);

        ratio1 = NaN;
        ratio10 = NaN;
        late_early_ratio = NaN;

        if (~is_invalid) && (~is_weak_all)
            ratio1  = A1  / max(A100, eps_num);
            ratio10 = A10 / max(A100, eps_num);

            if inc_1_10 > eps_num
                late_early_ratio = inc_10_100 / inc_1_10;
            else
                if inc_10_100 > eps_num
                    late_early_ratio = inf;
                else
                    late_early_ratio = 0;
                end
            end
        end

        is_early_saturate = false;
        if (~is_invalid) && (~is_weak_all)
            cond_sat1 = ratio10 > ratio10_sat_th;
            cond_sat2 = inc_10_100 < late_vs_early_th * max(inc_1_10, 0);
            is_early_saturate = cond_sat1 || cond_sat2;
        end

        is_progressive = false;
        if (~is_invalid) && (~is_weak_all) && (~is_early_saturate)
            cond_p1 = ratio1 <= ratio1_prog_max;
            cond_p2 = ratio10 <= ratio10_prog_max;
            cond_p3 = gradual >= gradual_min;
            cond_p4 = late_early_ratio >= late_early_ratio_min;
            cond_p5 = (A100 > A10) && (A10 > A1);

            is_progressive = cond_p1 && cond_p2 && cond_p3 && cond_p4 && cond_p5;
        end

        progress_score = compute_progress_score( ...
            A1, A10, A100, gradual, ...
            ratio1, ratio10, late_early_ratio, ...
            is_invalid, is_weak_all, is_early_saturate);

        %% 7.9 保存
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
        results(ic).dV100             = dV100;

        results(ic).A1                = A1;
        results(ic).A10               = A10;
        results(ic).A100              = A100;

        results(ic).gradual           = gradual;
        results(ic).inc_1_10          = inc_1_10;
        results(ic).inc_10_100        = inc_10_100;
        results(ic).ratio1            = ratio1;
        results(ic).ratio10           = ratio10;
        results(ic).late_early_ratio  = late_early_ratio;
        results(ic).progress_score    = progress_score;

        results(ic).is_invalid        = is_invalid;
        results(ic).is_weak_all       = is_weak_all;
        results(ic).is_early_saturate = is_early_saturate;
        results(ic).is_progressive    = is_progressive;

        results(ic).valid0            = valid0;
        results(ic).valid1            = valid_vec(1);
        results(ic).valid10           = valid_vec(2);
        results(ic).valid100          = valid_vec(3);

        results(ic).P0                = P0;
        results(ic).P1                = Pread_vec(1);
        results(ic).P10               = Pread_vec(2);
        results(ic).P100              = Pread_vec(3);

        results(ic).Vth0              = Vth0;
        results(ic).Vth1              = Vth_vec(1);
        results(ic).Vth10             = Vth_vec(2);
        results(ic).Vth100            = Vth_vec(3);

        results(ic).IDmin1            = IDmin_vec(1);
        results(ic).IDmin10           = IDmin_vec(2);
        results(ic).IDmin100          = IDmin_vec(3);

        results(ic).IDmax1            = IDmax_vec(1);
        results(ic).IDmax10           = IDmax_vec(2);
        results(ic).IDmax100          = IDmax_vec(3);

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

        results(ic).progress_score    = -1e9;
        results(ic).is_invalid        = true;
    end
end

%% =========================
% 8) 生成表
% =========================
T = struct2table(results);

writetable(T, fullfile(outdir, 'all_cases.csv'));

T_gradual = sortrows(T, 'gradual', 'descend');
writetable(T_gradual(1:min(20,height(T_gradual)), :), ...
    fullfile(outdir, 'top20_by_gradual.csv'));

T_valid = T(~T.is_invalid, :);
T_progscore = sortrows(T_valid, 'progress_score', 'descend');
writetable(T_progscore(1:min(20,height(T_progscore)), :), ...
    fullfile(outdir, 'top20_by_progress_score.csv'));

T_prog = T(T.is_progressive, :);
if ~isempty(T_prog)
    T_prog = sortrows(T_prog, 'progress_score', 'descend');
end
writetable(T_prog, fullfile(outdir, 'progressive_only.csv'));

disp(' ');
disp('================ Top-20 by gradual ================');
disp(T_gradual(1:min(20,height(T_gradual)), :));

disp(' ');
disp('================ Top-20 by progress_score ================');
disp(T_progscore(1:min(20,height(T_progscore)), :));

%% =========================
% 9) 画图
% =========================
case_idx = (1:height(T)).';

A1   = T.A1;
A10  = T.A10;
A100 = T.A100;

idx_invalid = T.is_invalid;
idx_weak    = T.is_weak_all & ~idx_invalid;
idx_early   = T.is_early_saturate & ~idx_invalid & ~idx_weak;
idx_prog    = T.is_progressive & ~idx_invalid & ~idx_weak;
idx_other   = ~idx_invalid & ~idx_weak & ~idx_early & ~idx_prog;

% 图1
figure('Name', 'Figure 1');
plot(case_idx, A1,   'o-', 'DisplayName', '|dV1|'); hold on;
plot(case_idx, A10,  's-', 'DisplayName', '|dV10|');
plot(case_idx, A100, 'd-', 'DisplayName', '|dV100|');
xlabel('Case index');
ylabel('Magnitude (V)');
title('|dV1|, |dV10|, |dV100| across cases');
grid on;
legend('Location', 'best');
saveas(gcf, fullfile(outdir, 'fig1_magnitudes.png'));

% 图2：只画有效 ratio
figure('Name', 'Figure 2');
valid_ratio = ~idx_invalid & ~idx_weak;
plot(case_idx(valid_ratio), T.ratio1(valid_ratio), 'o-', ...
    'DisplayName', 'ratio1 = |dV1|/|dV100|'); hold on;
plot(case_idx(valid_ratio), T.ratio10(valid_ratio), 's-', ...
    'DisplayName', 'ratio10 = |dV10|/|dV100|');
xlabel('Case index');
ylabel('Ratio');
title('ratio1, ratio10 (valid cases only)');
grid on;
legend('Location', 'best');
saveas(gcf, fullfile(outdir, 'fig2_ratios_valid_only.png'));

% 图3
figure('Name', 'Figure 3');
plot(A1(idx_weak),  A100(idx_weak),  'o', 'DisplayName', 'weak_all'); hold on;
plot(A1(idx_early), A100(idx_early), 's', 'DisplayName', 'early_saturate');
plot(A1(idx_prog),  A100(idx_prog),  'd', 'DisplayName', 'progressive');
plot(A1(idx_other), A100(idx_other), '^', 'DisplayName', 'other_valid');
xymax = max([A1; A100; 1e-3]) * 1.05;
plot([0 xymax], [0 xymax], 'k--', 'DisplayName', 'y=x');
xlabel('|dV1| (V)');
ylabel('|dV100| (V)');
title('|dV1| vs |dV100|');
grid on;
legend('Location', 'best');
saveas(gcf, fullfile(outdir, 'fig3_dV1_vs_dV100.png'));

% 图4
figure('Name', 'Figure 4');
plot(A10(idx_weak),  A100(idx_weak),  'o', 'DisplayName', 'weak_all'); hold on;
plot(A10(idx_early), A100(idx_early), 's', 'DisplayName', 'early_saturate');
plot(A10(idx_prog),  A100(idx_prog),  'd', 'DisplayName', 'progressive');
plot(A10(idx_other), A100(idx_other), '^', 'DisplayName', 'other_valid');
xymax = max([A10; A100; 1e-3]) * 1.05;
plot([0 xymax], [0 xymax], 'k--', 'DisplayName', 'y=x');
xlabel('|dV10| (V)');
ylabel('|dV100| (V)');
title('|dV10| vs |dV100|');
grid on;
legend('Location', 'best');
saveas(gcf, fullfile(outdir, 'fig4_dV10_vs_dV100.png'));

% 图5
figure('Name', 'Figure 5');
valid_inc = ~idx_invalid;
plot(case_idx(valid_inc), T.inc_1_10(valid_inc), 'o-', ...
    'DisplayName', 'inc_{1->10}'); hold on;
plot(case_idx(valid_inc), T.inc_10_100(valid_inc), 's-', ...
    'DisplayName', 'inc_{10->100}');
yline(0, 'k--', 'DisplayName', 'zero');
xlabel('Case index');
ylabel('Increment (V)');
title('Incremental growth');
grid on;
legend('Location', 'best');
saveas(gcf, fullfile(outdir, 'fig5_increments.png'));

% 图6
figure('Name', 'Figure 6');
valid_le = ~idx_invalid & ~idx_weak & isfinite(T.late_early_ratio);
plot(case_idx(valid_le), T.late_early_ratio(valid_le), 'o-', ...
    'DisplayName', 'late_early_ratio'); hold on;
yline(1.0, 'k--', 'DisplayName', '1.0');
xlabel('Case index');
ylabel('late\_early\_ratio');
title('late\_early\_ratio = (|dV100|-|dV10|)/(|dV10|-|dV1|)');
grid on;
legend('Location', 'best');
saveas(gcf, fullfile(outdir, 'fig6_late_early_ratio.png'));

% 图7
figure('Name', 'Figure 7');
valid_score = ~idx_invalid;
plot(case_idx(valid_score), T.progress_score(valid_score), 'o-', ...
    'DisplayName', 'progress_score');
xlabel('Case index');
ylabel('Score');
title('Progress score across cases');
grid on;
legend('Location', 'best');
saveas(gcf, fullfile(outdir, 'fig7_progress_score.png'));

%% =========================
% 10) 保存工作区
% =========================
save(fullfile(outdir, 'workspace_results.mat'), ...
     'T', 'T_gradual', 'T_progscore', 'T_prog', ...
     'results', 'case_list', ...
     'amp_list', 'pw_list', 'Ea_mean_list', 'Ea_sigma_list', ...
     'tfe_list', 'epife_read_list', 'Iref_list', 'VFB_read_list', 'VG_name_list', ...
     'delay', 'transit', 'cycle', 'step', ...
     'Pr', 'tauo', 'alpha', 'bet', 'epife_write', 'Ndom', ...
     'til', 'miu', 'Na', 'T', 'W', 'L', 'VD', 'VS', ...
     'Nread', ...
     'eps_num', 'Vmin_eval', ...
     'ratio10_sat_th', 'late_vs_early_th', ...
     'ratio1_prog_max', 'ratio10_prog_max', ...
     'late_early_ratio_min', 'gradual_min');

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
    'dV100', NaN, ...
    'A1', NaN, ...
    'A10', NaN, ...
    'A100', NaN, ...
    'gradual', NaN, ...
    'inc_1_10', NaN, ...
    'inc_10_100', NaN, ...
    'ratio1', NaN, ...
    'ratio10', NaN, ...
    'late_early_ratio', NaN, ...
    'progress_score', NaN, ...
    'is_invalid', false, ...
    'is_weak_all', false, ...
    'is_early_saturate', false, ...
    'is_progressive', false, ...
    'valid0', false, ...
    'valid1', false, ...
    'valid10', false, ...
    'valid100', false, ...
    'P0', NaN, ...
    'P1', NaN, ...
    'P10', NaN, ...
    'P100', NaN, ...
    'Vth0', NaN, ...
    'Vth1', NaN, ...
    'Vth10', NaN, ...
    'Vth100', NaN, ...
    'IDmin1', NaN, ...
    'IDmin10', NaN, ...
    'IDmin100', NaN, ...
    'IDmax1', NaN, ...
    'IDmax10', NaN, ...
    'IDmax100', NaN, ...
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

function score = compute_progress_score( ...
    A1, A10, A100, gradual, ...
    ratio1, ratio10, late_early_ratio, ...
    is_invalid, is_weak_all, is_early_saturate)

if is_invalid
    score = -1e9;
    return;
end

if is_weak_all
    score = -1e6 - A100;
    return;
end

term_window   =  2.0 * A100;
term_gradual  =  3.0 * gradual;
term_late     =  0.03 * min(late_early_ratio, 5);

term_ratio1   = -1.0 * ratio1;
term_ratio10  = -1.5 * ratio10;

term_sat = 0;
if is_early_saturate
    term_sat = -0.20;
end

score = term_window + term_gradual + term_late + ...
        term_ratio1 + term_ratio10 + term_sat;
end

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

if abs(x(2) - x(1)) < 1e-30
    Vth = mean(y);
else
    Vth = interp1(x, y, log10(Iref), 'linear');
end

is_valid = isfinite(Vth);
end