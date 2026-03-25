%% =========================================================
% sweep_mapping_gain_cand3.m
%
% 目的：
%   固定已经找到的 progressive accumulation window（cand3），
%   只扫读出映射参数，检查能否把同样的 FE 累积读成更大的 DeltaVth。
%
% 固定 FE 工作点：
%   cand3:
%     amp = 2.2 V
%     pw = 0.2 us
%     Ea_mean = 0.95
%     Ea_sigma = 0.20
%
% 扫描参数：
%   tfe
%   epife_read
%   Iref
%   VFB_read
%   VG scan range
%
% 依赖：
%   wfdef_acc_fix.m
%   get_FE_state.m
%   get_ID.m
%
% 说明：
%   这版脚本故意采用“固定 FE 轨迹，只扫读出映射”的模式：
%
%   1) get_FE_state 中使用固定 FE 参数：
%      epife_fe = 28
%
%   2) local_extract_vth / get_ID 中使用 sweep 的 epife_read, tfe 等
%
%   这样做的含义是：
%   先回答“同样的 FE 状态变化，能不能被读成更大的 DeltaVth？”
%
%   这不是 fully self-consistent 的全物理 sweep。
%   如果后面你要做“FE 与读出都同步变化”的版本，再单独做一版。
%% =========================================================

clear;
clc;
close all;

rng(0);

%% =========================
% 0) 输出目录
% =========================
outdir = 'mapping_gain_cand3_out';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

%% =========================
% 1) 固定 FE candidate：cand3
% =========================
cand_name = 'cand3_A2p20_PW0p20us_Em0p95_Es0p20';

amp      = 2.2;
pw       = 0.2e-6;
Ea_mean  = 0.95;
Ea_sigma = 0.20;

%% =========================
% 2) 固定时间轴参数
% =========================
delay   = 1e-6;        % 先固定在中间值
transit = 0.02e-6;
cycle   = 100;
step    = 0.01e-6;

%% =========================
% 3) FE 主参数（固定）
% =========================
Pr       = 25;
tauo     = 1e-9;
alpha    = 2;
bet      = 2;
epife_fe = 28;         % 仅用于 get_FE_state
Ndom     = 10000;

%% =========================
% 4) MOS / stack 固定参数
% =========================
til = 1e-7;
miu = 50;
Na  = 3e17;
T   = 300;
W   = 1;
L   = 1;
VD  = 0.05;
VS  = 0;

%% =========================
% 5) 读点
% =========================
Nread = [1 2 3 4 5 10 20 50 100].';
Nread = unique(Nread);

%% =========================
% 6) 扫描参数
% =========================
tfe_list        = [0.8e-6 1.0e-6 1.2e-6 1.5e-6];
epife_read_list = [20 24 28];
Iref_list       = [1e-8 3e-8 1e-7];
VFB_read_list   = [-0.7 -0.5 -0.3];

% 不同 VG 扫描窗口
VG_case = struct([]);

VG_case(1).name = 'VGm0p8_to_2p2';
VG_case(1).VG   = -0.8:0.02:2.2;

VG_case(2).name = 'VGm0p5_to_1p7';
VG_case(2).VG   = -0.5:0.02:1.7;

VG_case(3).name = 'VGm1p0_to_2p5';
VG_case(3).VG   = -1.0:0.02:2.5;

%% =========================
% 7) 固定随机样本
% =========================
Weight  = ones(Ndom,1) / Ndom;
St_init = -ones(Ndom,1);

base_randn_Ea = randn(Ndom,1);
r_Ea = Ea_mean + Ea_sigma * base_randn_Ea;
r_Ea(r_Ea <= 0) = Ea_mean;

% 先保留 baseline：关闭 voff
r_voff = zeros(Ndom,1);

%% =========================
% 8) 先只算一次 FE 轨迹
% =========================
fprintf('=====================================================\n');
fprintf('Precompute FE trajectory for %s\n', cand_name);
fprintf('amp=%.3f V, pw=%g s, Ea_mean=%.3f, Ea_sigma=%.3f\n', ...
    amp, pw, Ea_mean, Ea_sigma);
fprintf('delay=%g s, cycle=%d, step=%g s\n', delay, cycle, step);
fprintf('=====================================================\n');

[time, volt, index, index_pre, index_end] = ...
    wfdef_acc_fix(amp, pw, step, delay, transit, cycle);

index     = index(:);
index_pre = index_pre(:);
index_end = index_end(:);

tic;
[vfev, Stsum] = get_FE_state(time, volt, St_init, Weight, r_Ea, r_voff, ...
    Pr, tauo, alpha, bet, epife_fe, Ndom);
telapsed_fe = toc;

if isvector(Stsum)
    st_trace = Stsum(:);
else
    st_trace = Stsum(:,1);
end

fprintf('FE trajectory ready. elapsed = %.3f s\n\n', telapsed_fe);

%% =========================
% 9) 扫 mapping 参数
% =========================
Result = struct([]);
icase = 0;

nTotal = numel(tfe_list) * numel(epife_read_list) * ...
         numel(Iref_list) * numel(VFB_read_list) * numel(VG_case);

fprintf('Total mapping cases = %d\n\n', nTotal);

for it = 1:numel(tfe_list)
    tfe = tfe_list(it);

    for ie = 1:numel(epife_read_list)
        epife_read = epife_read_list(ie);

        for ii = 1:numel(Iref_list)
            Iref = Iref_list(ii);

            for iv = 1:numel(VFB_read_list)
                VFB_read = VFB_read_list(iv);

                for ig = 1:numel(VG_case)
                    VG_name = VG_case(ig).name;
                    VG = VG_case(ig).VG;

                    icase = icase + 1;
                    fprintf('=====================================================\n');
                    fprintf('Case %d / %d\n', icase, nTotal);
                    fprintf('tfe=%g, epife_read=%g, Iref=%g, VFB_read=%g, %s\n', ...
                        tfe, epife_read, Iref, VFB_read, VG_name);
                    fprintf('=====================================================\n');

                    %% 9.1 初始点 N=0
                    idx0 = index_pre(1);
                    P0 = st_trace(idx0);

                    [Vth0, ID0, idmin0, idmax0, valid0] = local_extract_vth( ...
                        P0, epife_read, tfe, til, miu, Na, T, W, L, ...
                        VG, VFB_read, VD, VS, Iref);

                    %% 9.2 各读点
                    nRead = numel(Nread);

                    Pread_case    = nan(nRead,1);
                    Vth_case      = nan(nRead,1);
                    DeltaVth_case = nan(nRead,1);
                    idmin_case    = nan(nRead,1);
                    idmax_case    = nan(nRead,1);
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

                        [Vth_case(k), ~, idmin_case(k), idmax_case(k), valid_case(k)] = ...
                            local_extract_vth(Pread_case(k), epife_read, tfe, til, miu, Na, T, W, L, ...
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

                    %% 9.3 统计指标
                    idx1   = find(Nread == 1,   1);
                    idx100 = find(Nread == 100, 1);

                    dV1 = nan;
                    dV100 = nan;
                    gradual = nan;
                    dP1 = nan;
                    dP100 = nan;

                    if ~isempty(idx1)
                        dV1 = DeltaVth_case(idx1);
                        dP1 = Pread_case(idx1) - P0;
                    end

                    if ~isempty(idx100)
                        dV100 = DeltaVth_case(idx100);
                        dP100 = Pread_case(idx100) - P0;
                    end

                    if ~isempty(idx1) && ~isempty(idx100)
                        gradual = abs(DeltaVth_case(idx100) - DeltaVth_case(idx1));
                    end

                    valid_ratio = mean(valid_case);

                    % 简单分类
                    is_invalid = valid_ratio < 0.8;
                    is_big_dv100 = abs(dV100) >= 0.30;
                    is_very_big_dv100 = abs(dV100) >= 0.50;

                    %% 9.4 保存 case
                    case_id = sprintf('%s_tfe%0.2e_epi%0.1f_Iref%0.1e_VFB%0.2f_%s', ...
                        cand_name, tfe, epife_read, Iref, VFB_read, VG_name);
                    case_id = strrep(case_id, '.', 'p');
                    case_id = strrep(case_id, '-', 'm');

                    save(fullfile(outdir, [case_id '.mat']), ...
                        'cand_name', 'amp', 'pw', 'Ea_mean', 'Ea_sigma', ...
                        'delay', 'transit', 'cycle', 'step', ...
                        'Pr', 'tauo', 'alpha', 'bet', 'epife_fe', 'Ndom', ...
                        'tfe', 'epife_read', 'til', 'miu', 'Na', 'T', 'W', 'L', ...
                        'VD', 'VS', 'VG', 'VG_name', 'Iref', 'VFB_read', ...
                        'r_Ea', 'r_voff', ...
                        'time', 'volt', 'index', 'index_pre', 'index_end', ...
                        'vfev', 'Stsum', 'st_trace', ...
                        'P0', 'Vth0', 'ID0', 'idmin0', 'idmax0', 'valid0', ...
                        'Nread', 'Pread_case', 'Vth_case', 'DeltaVth_case', ...
                        'idmin_case', 'idmax_case', 'valid_case', ...
                        'dV1', 'dV100', 'gradual', 'dP1', 'dP100', ...
                        'valid_ratio', 'is_invalid', 'is_big_dv100', 'is_very_big_dv100');

                    %% 9.5 汇总
                    Result(icase).case_id          = case_id;
                    Result(icase).tfe              = tfe;
                    Result(icase).epife_read       = epife_read;
                    Result(icase).Iref             = Iref;
                    Result(icase).VFB_read         = VFB_read;
                    Result(icase).VG_name          = VG_name;

                    Result(icase).P0               = P0;
                    Result(icase).Vth0             = Vth0;
                    Result(icase).idmin0           = idmin0;
                    Result(icase).idmax0           = idmax0;
                    Result(icase).valid0           = valid0;

                    Result(icase).dV1              = dV1;
                    Result(icase).dV100            = dV100;
                    Result(icase).gradual          = gradual;
                    Result(icase).dP1              = dP1;
                    Result(icase).dP100            = dP100;
                    Result(icase).valid_ratio      = valid_ratio;

                    Result(icase).is_invalid       = is_invalid;
                    Result(icase).is_big_dv100     = is_big_dv100;
                    Result(icase).is_very_big_dv100= is_very_big_dv100;

                    fprintf('dV1 = %.6f V, dV100 = %.6f V, gradual = %.6f V\n', ...
                        dV1, dV100, gradual);
                    fprintf('dP1 = %.6e, dP100 = %.6e, valid_ratio = %.3f\n', ...
                        dP1, dP100, valid_ratio);
                    fprintf('idmin0 = %.3e, idmax0 = %.3e, valid0 = %d\n\n', ...
                        idmin0, idmax0, valid0);
                end
            end
        end
    end
end

%% =========================
% 10) 汇总表
% =========================
summary_table = struct2table(Result);

save(fullfile(outdir, 'summary_mapping_gain_cand3.mat'), ...
    'summary_table', 'Result', 'Nread');

writetable(summary_table, ...
    fullfile(outdir, 'summary_mapping_gain_cand3.csv'));

disp('================ summary_table ================');
disp(summary_table(:, { ...
    'case_id','tfe','epife_read','Iref','VFB_read','VG_name', ...
    'dV1','dV100','gradual','dP100','valid_ratio', ...
    'is_invalid','is_big_dv100','is_very_big_dv100'}));

%% =========================
% 11) Top-20 by |dV100|
% =========================
[~, idx_sort_dV100] = sort(abs(summary_table.dV100), 'descend');
topk = idx_sort_dV100(1:min(20, height(summary_table)));

top20_by_dV100 = summary_table(topk, { ...
    'case_id','tfe','epife_read','Iref','VFB_read','VG_name', ...
    'dV1','dV100','gradual','dP100','valid_ratio', ...
    'is_invalid','is_big_dv100','is_very_big_dv100'});

disp('================ Top-20 by |dV100| ================');
disp(top20_by_dV100);

writetable(top20_by_dV100, fullfile(outdir, 'top20_by_dV100.csv'));

%% =========================
% 12) Top-20 by gradual
% =========================
[~, idx_sort_grad] = sort(summary_table.gradual, 'descend');
topk_grad = idx_sort_grad(1:min(20, height(summary_table)));

top20_by_gradual = summary_table(topk_grad, { ...
    'case_id','tfe','epife_read','Iref','VFB_read','VG_name', ...
    'dV1','dV100','gradual','dP100','valid_ratio', ...
    'is_invalid','is_big_dv100','is_very_big_dv100'});

disp('================ Top-20 by gradual ================');
disp(top20_by_gradual);

writetable(top20_by_gradual, fullfile(outdir, 'top20_by_gradual.csv'));

%% =========================
% 13) 分类统计
% =========================
n_invalid = sum(summary_table.is_invalid);
n_big     = sum(summary_table.is_big_dv100);
n_verybig = sum(summary_table.is_very_big_dv100);

fprintf('\n================ Class counts ================\n');
fprintf('invalid         = %d / %d\n', n_invalid, height(summary_table));
fprintf('|dV100| >= 0.30 = %d / %d\n', n_big,     height(summary_table));
fprintf('|dV100| >= 0.50 = %d / %d\n', n_verybig, height(summary_table));

%% =========================
% 14) 图 1：|dV100| across cases
% =========================
figure;
plot(abs(summary_table.dV100), 'o-', 'LineWidth', 1.2);
xlabel('Case index');
ylabel('|dV100| (V)');
title('Mapping sweep: |DeltaVth(N=100)| across cases');
grid on;

%% =========================
% 15) 图 2：gradual across cases
% =========================
figure;
plot(summary_table.gradual, 's-', 'LineWidth', 1.2);
xlabel('Case index');
ylabel('gradual = |dV100 - dV1| (V)');
title('Mapping sweep: gradual across cases');
grid on;

%% =========================
% 16) 图 3：|dP100| vs |dV100|
% =========================
figure;
scatter(abs(summary_table.dP100), abs(summary_table.dV100), 35, 'filled');
xlabel('|dP100|');
ylabel('|dV100| (V)');
title('|dP100| vs |dV100|');
grid on;

%% =========================
% 17) 图 4：|dV1| vs |dV100|
% =========================
figure;
scatter(abs(summary_table.dV1), abs(summary_table.dV100), 35, 'filled');
hold on;
mx = max([abs(summary_table.dV1); abs(summary_table.dV100)]) * 1.05;
if isempty(mx) || ~isfinite(mx) || mx <= 0
    mx = 1;
end
plot([0 mx], [0 mx], 'r--', 'LineWidth', 1.2);
xlabel('|dV1| (V)');
ylabel('|dV100| (V)');
title('|dV1| vs |dV100|');
grid on;
axis equal;
xlim([0 mx]);
ylim([0 mx]);
hold off;

%% =========================
% 18) 代表性 case 图
%     选三类：
%       a) |dV100| 最大
%       b) gradual 最大
%       c) 第一个 |dV100| >= 0.50 的 case（如果有）
% =========================
rep_idx = [];

if ~isempty(idx_sort_dV100)
    rep_idx(end+1) = idx_sort_dV100(1);
end
if ~isempty(idx_sort_grad)
    rep_idx(end+1) = idx_sort_grad(1);
end

idx_verybig = find(summary_table.is_very_big_dv100, 1, 'first');
if ~isempty(idx_verybig)
    rep_idx(end+1) = idx_verybig;
end

rep_idx = unique(rep_idx);

for ir = 1:numel(rep_idx)
    ii = rep_idx(ir);
    case_id = summary_table.case_id{ii};
    S = load(fullfile(outdir, [case_id '.mat']));

    figure;
    plot(S.Nread, S.DeltaVth_case, 'o-', 'LineWidth', 1.5);
    xlabel('Cumulative pulse count N');
    ylabel('\DeltaV_{th} (V)');
    title(['Representative case: ' strrep(case_id, '_', '\_')]);
    grid on;

    figure;
    plot(S.Nread, S.Pread_case, 's-', 'LineWidth', 1.5);
    xlabel('Cumulative pulse count N');
    ylabel('P_{read}');
    title(['Representative case P_{read}: ' strrep(case_id, '_', '\_')]);
    grid on;
end

fprintf('\nAll done. Results saved in folder: %s\n', outdir);

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