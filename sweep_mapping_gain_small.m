%% =========================================================
% sweep_mapping_gain_small.m
%
% 目的：
%   走 B 路线，只围绕 cand3 做读出映射增益的小范围精扫
%   主扫参数：
%       1) tfe
%       2) epife_read
%
%   固定 cand3 写入主干参数：
%       amp     = 2.2 V
%       pw      = 0.2 us
%       Ea_mean = 0.95
%       Ea_sigma= 0.20
%
% 依赖：
%   wfdef_acc_fix.m
%   get_FE_state.m
%   get_ID.m
%
% 输出：
%   mapping_gain_small_out/
%       all_cases.csv
%       top_by_dV100.csv
%       top_by_gradual.csv
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
outdir = 'mapping_gain_small_out';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

diary(fullfile(outdir, 'run_log.txt'));
fprintf('=== sweep_mapping_gain_small ===\n');

%% =========================
% 1) 固定 cand3 写入参数
% =========================
amp      = 2.2;
pw       = 0.2e-6;
Ea_mean  = 0.95;
Ea_sigma = 0.20;

%% =========================
% 2) 只扫 mapping / readout 增益参数
% =========================
tfe_list        = [1.2e-6 1.5e-6 1.8e-6 2.0e-6];
epife_read_list = [16 18 20 22];

Iref_list       = [1e-7];
VFB_read_list   = [-0.5];
VG_name_list    = {'VGm0p8_to_2p2'};

%% =========================
% 3) 固定 pulse / FE / MOS 参数
% =========================
delay   = 1e-6;
transit = 0.02e-6;
cycle   = 100;
step    = 0.01e-6;

% FE 参数
Pr    = 25;
tauo  = 1e-9;
alpha = 2;
bet   = 2;

% 写入主干用 FE eps
epife_write = 28;

Ndom = 10000;

Weight  = ones(Ndom,1) / Ndom;
St_init = -ones(Ndom,1);

% 固定同一组 domain sample，避免 mapping 扫描时随机性污染结论
use_global_fixed_domain_sample = true;

if use_global_fixed_domain_sample
    r_Ea_base = Ea_mean + Ea_sigma * randn(Ndom,1);
    r_Ea_base(r_Ea_base <= 0) = Ea_mean;
end

% 暂不引入 voff 演化
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
% 4) 读出点
% =========================
Nread = [1 10 30 100].';

%% =========================
% 5) 判据
% =========================
eps_num = 1e-15;
Vmin_eval = 5e-3;

% 不希望前期太快吃满，但当前 B 路线允许一定前段抬升
ratio1_soft_max = 0.60;   % |dV1| / |dV100|
ratio30_soft_max = 0.85;  % |dV30| / |dV100|

%% =========================
% 6) case 列表
% =========================
case_list = {};
icase = 0;

for itfe = 1:numel(tfe_list)
    for iepi = 1:numel(epife_read_list)
        for iir = 1:numel(Iref_list)
            for ivfb = 1:numel(VFB_read_list)
                for ivg = 1:numel(VG_name_list)
                    icase = icase + 1;

                    cfg = struct();
                    cfg.case_idx    = icase;

                    cfg.amp         = amp;
                    cfg.pw          = pw;
                    cfg.Ea_mean     = Ea_mean;
                    cfg.Ea_sigma    = Ea_sigma;

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

nCase = numel(case_list);
fprintf('Total cases = %d\n', nCase);

%% =========================
% 7) 预分配结果
% =========================
results = repmat(init_result_struct(), nCase, 1);

%% =========================
% 8) 主循环
% =========================
for ic = 1:nCase
    cfg = case_list{ic};

    fprintf('\n==================================================\n');
    fprintf('[%d / %d] %s\n', ic, nCase, cfg.case_id);
    fprintf('==================================================\n');

    try
        %% 8.1 当前 case 的 r_Ea
        if use_global_fixed_domain_sample
            r_Ea = r_Ea_base;
        else
            r_Ea = cfg.Ea_mean + cfg.Ea_sigma * randn(Ndom,1);
            r_Ea(r_Ea <= 0) = cfg.Ea_mean;
        end

        %% 8.2 波形
        [time, volt, index, index_pre, index_end] = ...
            wfdef_acc_fix(cfg.amp, cfg.pw, step, delay, transit, cycle);

        index     = index(:);
        index_pre = index_pre(:);
        index_end = index_end(:);

        %% 8.3 FE 状态推进
        tic;
        [vfev, Stsum] = get_FE_state(time, volt, St_init, Weight, r_Ea, r_voff, ...
                                     Pr, tauo, alpha, bet, epife_write, Ndom);
        telapsed = toc;
        fprintf('get_FE_state done, elapsed = %.3f s\n', telapsed);

        %% 8.4 单条轨迹
        if isvector(Stsum)
            st_trace = Stsum(:);
        else
            st_trace = Stsum(:,1);
        end

        %% 8.5 VG
        VG = make_VG_from_name(cfg.VG_name);

        %% 8.6 初始参考点
        idx0 = index_pre(1);
        P0 = st_trace(idx0);

        [Vth0, ~, ~, ~, valid0] = local_extract_vth( ...
            P0, cfg.epife_read, cfg.tfe, til, miu, Na, T0, W, L, ...
            VG, cfg.VFB_read, VD, VS, cfg.Iref);

        %% 8.7 提取 N=1,10,30,100
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

        if ~valid0
            first_valid = find(valid_vec, 1, 'first');
            if ~isempty(first_valid)
                Vth_ref = Vth_vec(first_valid);
                DeltaVth_vec = Vth_vec - Vth_ref;
            end
        end

        %% 8.8 关键量
        dV1   = DeltaVth_vec(1);
        dV10  = DeltaVth_vec(2);
        dV30  = DeltaVth_vec(3);
        dV100 = DeltaVth_vec(4);

        A1   = abs(dV1);
        A10  = abs(dV10);
        A30  = abs(dV30);
        A100 = abs(dV100);

        gradual      = A100 - A1;
        inc_1_10     = A10 - A1;
        inc_10_30    = A30 - A10;
        inc_30_100   = A100 - A30;

        %% 8.9 判据
        is_invalid = any(~isfinite([dV1 dV10 dV30 dV100])) || ...
                     any(~isfinite([A1 A10 A30 A100]));

        is_weak_all = (~is_invalid) && (A100 < Vmin_eval);

        ratio1  = NaN;
        ratio10 = NaN;
        ratio30 = NaN;

        if (~is_invalid) && (~is_weak_all)
            ratio1  = A1  / max(A100, eps_num);
            ratio10 = A10 / max(A100, eps_num);
            ratio30 = A30 / max(A100, eps_num);
        end

        % 这里只做“软判据”，不直接扔掉
        is_too_fast = false;
        if (~is_invalid) && (~is_weak_all)
            cond_fast1 = ratio1  > ratio1_soft_max;
            cond_fast2 = ratio30 > ratio30_soft_max;
            is_too_fast = cond_fast1 || cond_fast2;
        end

        score_B = compute_score_B( ...
            A1, A10, A30, A100, gradual, ...
            ratio1, ratio10, ratio30, ...
            inc_1_10, inc_10_30, inc_30_100, ...
            is_invalid, is_weak_all, is_too_fast);

        %% 8.10 保存
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
        results(ic).dV10         = dV10;
        results(ic).dV30         = dV30;
        results(ic).dV100        = dV100;

        results(ic).A1           = A1;
        results(ic).A10          = A10;
        results(ic).A30          = A30;
        results(ic).A100         = A100;

        results(ic).gradual      = gradual;
        results(ic).inc_1_10     = inc_1_10;
        results(ic).inc_10_30    = inc_10_30;
        results(ic).inc_30_100   = inc_30_100;

        results(ic).ratio1       = ratio1;
        results(ic).ratio10      = ratio10;
        results(ic).ratio30      = ratio30;

        results(ic).score_B      = score_B;

        results(ic).is_invalid   = is_invalid;
        results(ic).is_weak_all  = is_weak_all;
        results(ic).is_too_fast  = is_too_fast;

        results(ic).valid0       = valid0;
        results(ic).valid1       = valid_vec(1);
        results(ic).valid10      = valid_vec(2);
        results(ic).valid30      = valid_vec(3);
        results(ic).valid100     = valid_vec(4);

        results(ic).P0           = P0;
        results(ic).P1           = Pread_vec(1);
        results(ic).P10          = Pread_vec(2);
        results(ic).P30          = Pread_vec(3);
        results(ic).P100         = Pread_vec(4);

        results(ic).Vth0         = Vth0;
        results(ic).Vth1         = Vth_vec(1);
        results(ic).Vth10        = Vth_vec(2);
        results(ic).Vth30        = Vth_vec(3);
        results(ic).Vth100       = Vth_vec(4);

        results(ic).IDmin1       = IDmin_vec(1);
        results(ic).IDmin10      = IDmin_vec(2);
        results(ic).IDmin30      = IDmin_vec(3);
        results(ic).IDmin100     = IDmin_vec(4);

        results(ic).IDmax1       = IDmax_vec(1);
        results(ic).IDmax10      = IDmax_vec(2);
        results(ic).IDmax30      = IDmax_vec(3);
        results(ic).IDmax100     = IDmax_vec(4);

        results(ic).elapsed_s    = telapsed;

    catch ME
        warning('Case failed: %s', cfg.case_id);
        warning('%s', ME.message);

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

        results(ic).score_B      = -1e9;
        results(ic).is_invalid   = true;
    end
end

%% =========================
% 9) 生成表
% =========================
T = struct2table(results);

writetable(T, fullfile(outdir, 'all_cases.csv'));

T_valid = T(~T.is_invalid, :);

T_dV100 = sortrows(T_valid, 'A100', 'descend');
writetable(T_dV100(1:min(20,height(T_dV100)), :), ...
    fullfile(outdir, 'top_by_dV100.csv'));

T_gradual = sortrows(T_valid, 'gradual', 'descend');
writetable(T_gradual(1:min(20,height(T_gradual)), :), ...
    fullfile(outdir, 'top_by_gradual.csv'));

T_score = sortrows(T_valid, 'score_B', 'descend');
writetable(T_score(1:min(20,height(T_score)), :), ...
    fullfile(outdir, 'top_by_score_B.csv'));

disp(' ');
disp('================ Top-20 by |dV100| ================');
disp(T_dV100(1:min(20,height(T_dV100)), :));

disp(' ');
disp('================ Top-20 by gradual ================');
disp(T_gradual(1:min(20,height(T_gradual)), :));

disp(' ');
disp('================ Top-20 by score_B ================');
disp(T_score(1:min(20,height(T_score)), :));

%% =========================
% 10) 画图
% =========================
case_idx = (1:height(T)).';

idx_invalid = T.is_invalid;
idx_weak    = T.is_weak_all & ~idx_invalid;
idx_fast    = T.is_too_fast & ~idx_invalid & ~idx_weak;
idx_good    = ~idx_invalid & ~idx_weak & ~idx_fast;
idx_other   = ~idx_invalid & ~idx_weak & ~idx_good & ~idx_fast;

% 图1
figure('Name', 'Figure 1');
plot(case_idx, T.A1,   'o-', 'DisplayName', '|dV1|'); hold on;
plot(case_idx, T.A30,  's-', 'DisplayName', '|dV30|');
plot(case_idx, T.A100, 'd-', 'DisplayName', '|dV100|');
xlabel('Case index');
ylabel('Magnitude (V)');
title('|dV1|, |dV30|, |dV100| across mapping cases');
grid on;
legend('Location', 'best');
saveas(gcf, fullfile(outdir, 'fig1_magnitudes.png'));

% 图2
figure('Name', 'Figure 2');
valid_ratio = ~idx_invalid & ~idx_weak;
plot(case_idx(valid_ratio), T.ratio1(valid_ratio), 'o-', ...
    'DisplayName', 'ratio1 = |dV1|/|dV100|'); hold on;
plot(case_idx(valid_ratio), T.ratio30(valid_ratio), 's-', ...
    'DisplayName', 'ratio30 = |dV30|/|dV100|');
yline(ratio1_soft_max, 'k--', 'DisplayName', 'ratio1 soft max');
yline(ratio30_soft_max, 'r--', 'DisplayName', 'ratio30 soft max');
xlabel('Case index');
ylabel('Ratio');
title('Early-growth ratios');
grid on;
legend('Location', 'best');
saveas(gcf, fullfile(outdir, 'fig2_ratios.png'));

% 图3
figure('Name', 'Figure 3');
plot(T.A1(idx_weak),   T.A100(idx_weak),   'o', 'DisplayName', 'weak_all'); hold on;
plot(T.A1(idx_fast),   T.A100(idx_fast),   's', 'DisplayName', 'too_fast');
plot(T.A1(idx_good),   T.A100(idx_good),   'd', 'DisplayName', 'good');
plot(T.A1(idx_other),  T.A100(idx_other),  '^', 'DisplayName', 'other');
xymax = max([T.A1; T.A100; 1e-3]) * 1.05;
plot([0 xymax], [0 xymax], 'k--', 'DisplayName', 'y=x');
xlabel('|dV1| (V)');
ylabel('|dV100| (V)');
title('|dV1| vs |dV100|');
grid on;
legend('Location', 'best');
saveas(gcf, fullfile(outdir, 'fig3_dV1_vs_dV100.png'));

% 图4
figure('Name', 'Figure 4');
plot(case_idx(~idx_invalid), T.gradual(~idx_invalid), 'o-', 'DisplayName', 'gradual = |dV100|-|dV1|');
xlabel('Case index');
ylabel('Gradual (V)');
title('Gradual window across mapping cases');
grid on;
legend('Location', 'best');
saveas(gcf, fullfile(outdir, 'fig4_gradual.png'));

% 图5
figure('Name', 'Figure 5');
plot(case_idx(~idx_invalid), T.score_B(~idx_invalid), 'o-', 'DisplayName', 'score_B');
xlabel('Case index');
ylabel('Score');
title('score_B across mapping cases');
grid on;
legend('Location', 'best');
saveas(gcf, fullfile(outdir, 'fig5_scoreB.png'));

%% =========================
% 11) 保存工作区
% =========================
save(fullfile(outdir, 'workspace_results.mat'), ...
     'T', 'T_dV100', 'T_gradual', 'T_score', ...
     'results', 'case_list', ...
     'amp', 'pw', 'Ea_mean', 'Ea_sigma', ...
     'tfe_list', 'epife_read_list', ...
     'Iref_list', 'VFB_read_list', 'VG_name_list', ...
     'delay', 'transit', 'cycle', 'step', ...
     'Pr', 'tauo', 'alpha', 'bet', 'epife_write', 'Ndom', ...
     'til', 'miu', 'Na', 'T0', 'W', 'L', 'VD', 'VS', ...
     'Nread', ...
     'ratio1_soft_max', 'ratio30_soft_max', ...
     'Vmin_eval');

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
    'A1', NaN, ...
    'A10', NaN, ...
    'A30', NaN, ...
    'A100', NaN, ...
    'gradual', NaN, ...
    'inc_1_10', NaN, ...
    'inc_10_30', NaN, ...
    'inc_30_100', NaN, ...
    'ratio1', NaN, ...
    'ratio10', NaN, ...
    'ratio30', NaN, ...
    'score_B', NaN, ...
    'is_invalid', false, ...
    'is_weak_all', false, ...
    'is_too_fast', false, ...
    'valid0', false, ...
    'valid1', false, ...
    'valid10', false, ...
    'valid30', false, ...
    'valid100', false, ...
    'P0', NaN, ...
    'P1', NaN, ...
    'P10', NaN, ...
    'P30', NaN, ...
    'P100', NaN, ...
    'Vth0', NaN, ...
    'Vth1', NaN, ...
    'Vth10', NaN, ...
    'Vth30', NaN, ...
    'Vth100', NaN, ...
    'IDmin1', NaN, ...
    'IDmin10', NaN, ...
    'IDmin30', NaN, ...
    'IDmin100', NaN, ...
    'IDmax1', NaN, ...
    'IDmax10', NaN, ...
    'IDmax30', NaN, ...
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

function score = compute_score_B( ...
    A1, A10, A30, A100, gradual, ...
    ratio1, ratio10, ratio30, ...
    inc_1_10, inc_10_30, inc_30_100, ...
    is_invalid, is_weak_all, is_too_fast)

if is_invalid
    score = -1e9;
    return;
end

if is_weak_all
    score = -1e6 - A100;
    return;
end

term_main    = 4.0 * A100;
term_gradual = 2.5 * gradual;

term_ratio1  = 0;
term_ratio30 = 0;
if isfinite(ratio1)
    term_ratio1 = -1.5 * ratio1;
end
if isfinite(ratio30)
    term_ratio30 = -1.5 * ratio30;
end

term_late = 1.5 * max(inc_30_100, 0);

term_fast = 0;
if is_too_fast
    term_fast = -0.8;
end

score = term_main + term_gradual + term_ratio1 + term_ratio30 + ...
        term_late + term_fast;
end