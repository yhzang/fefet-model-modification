%% =========================================================
% sweep_stage1_distribution_fine.m
%
% 目的：
%   1) 围绕当前最敏感区域做“分布主导细扫”
%   2) 判断当前 case 属于：
%        - small_all  : 整体都小
%        - fast_sat   : 首脉冲基本定型
%        - good_prog  : 有一定后续累积
%   3) 输出 top-20 case，辅助下一步判断
%
% 依赖：
%   wfdef_acc_fix.m
%   get_FE_state.m
%   get_ID.m
%
% 输出：
%   outdir/
%     summary_table_stage1_distribution_fine.csv
%     summary_table_stage1_distribution_fine.mat
%     top20_by_dV100.csv
%     top20_by_gradual.csv
%     每个 case 的 .mat 文件
%
% 作者备注：
%   这版脚本的重点不是直接找最终最优点，
%   而是先把“太弱 / 太快 / 可能可用”分开。
%% =========================================================

clear;
clc;
close all;

rng(0);

%% =========================
% 0) 输出目录
% =========================
outdir = 'stage1_distribution_fine_out';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

%% =========================
% 1) 固定参数
% =========================
VDD = 5;

delay   = 1e-6;
transit = 0.02e-6;
cycle   = 100;
step    = 0.01e-6;

% FE 固定参数
Pr    = 25;
tauo  = 1e-9;
alpha = 2;
bet   = 2;
epife = 28;
Ndom  = 10000;

% 读出固定参数
VG = -0.5:0.02:1.7;
VD = 0.05;
VS = 0;
Iref = 1e-7;
VFB_read = -0.5;

% MOS / stack 固定参数
tfe = 0.8e-6;
til = 1e-7;
miu = 50;
Na  = 3e17;
T   = 300;
W   = 1;
L   = 1;

% 读点
Nread = [1 2 3 4 5 10 20 50 100].';
Nread = unique(Nread);

fprintf('================ Setup ================\n');
fprintf('delay   = %g s\n', delay);
fprintf('transit = %g s\n', transit);
fprintf('cycle   = %d\n', cycle);
fprintf('step    = %g s\n', step);
fprintf('Pr      = %g\n', Pr);
fprintf('tauo    = %g\n', tauo);
fprintf('alpha   = %g\n', alpha);
fprintf('bet     = %g\n', bet);
fprintf('epife   = %g\n', epife);
fprintf('Ndom    = %d\n', Ndom);
fprintf('VG      = [%g : %g : %g]\n', VG(1), VG(2)-VG(1), VG(end));
fprintf('Iref    = %g A\n\n', Iref);

%% =========================
% 2) 扫描参数
%    重点细扫 Ea_mean / Ea_sigma
%% =========================
amp_list      = [1.8 2.2 2.6];
pw_list       = [0.2e-6 0.5e-6 1.0e-6];

Ea_mean_list  = [0.75 0.80 0.85 0.90 0.95 1.00];
Ea_sigma_list = [0.20 0.30 0.40 0.50 0.60 0.70];

%% =========================
% 3) 固定随机样本
%    保证不同 case 可比
%% =========================
Weight  = ones(Ndom,1) / Ndom;
St_init = -ones(Ndom,1);

base_randn_Ea = randn(Ndom,1);

% 先关掉 voff，避免 gap 段附加演化
r_voff = zeros(Ndom,1);

%% =========================
% 4) 预分配
%% =========================
nAmp = numel(amp_list);
nPw  = numel(pw_list);
nEm  = numel(Ea_mean_list);
nEs  = numel(Ea_sigma_list);

nTotal = nAmp * nPw * nEm * nEs;
fprintf('Total cases = %d\n\n', nTotal);

Result = struct([]);
icase = 0;

%% =========================
% 5) 主循环
%% =========================
for ia = 1:nAmp
    amp = amp_list(ia);

    for ipw = 1:nPw
        pw = pw_list(ipw);

        for iem = 1:nEm
            Ea_mean = Ea_mean_list(iem);

            for ies = 1:nEs
                Ea_sigma = Ea_sigma_list(ies);

                icase = icase + 1;
                fprintf('=====================================================\n');
                fprintf('Case %d / %d\n', icase, nTotal);
                fprintf('amp=%.3f V, pw=%g s, Ea_mean=%.3f, Ea_sigma=%.3f\n', ...
                    amp, pw, Ea_mean, Ea_sigma);
                fprintf('=====================================================\n');

                %% 5.1 生成 r_Ea
                r_Ea = Ea_mean + Ea_sigma * base_randn_Ea;
                r_Ea(r_Ea <= 0) = Ea_mean;

                %% 5.2 生成波形
                [time, volt, index, index_pre, index_end] = ...
                    wfdef_acc_fix(amp, pw, step, delay, transit, cycle);

                index     = index(:);
                index_pre = index_pre(:);
                index_end = index_end(:);

                %% 5.3 FE 演化
                tic;
                [vfev, Stsum] = get_FE_state(time, volt, St_init, Weight, r_Ea, r_voff, ...
                    Pr, tauo, alpha, bet, epife, Ndom);
                telapsed = toc;

                if isvector(Stsum)
                    st_trace = Stsum(:);
                else
                    st_trace = Stsum(:,1);
                end

                %% 5.4 初始参考点 N=0
                idx0 = index_pre(1);
                P0 = st_trace(idx0);

                [Vth0, ID0, idmin0, idmax0, valid0] = local_extract_vth( ...
                    P0, epife, tfe, til, miu, Na, T, W, L, ...
                    VG, VFB_read, VD, VS, Iref);

                %% 5.5 各读点
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
                        local_extract_vth(Pread_case(k), epife, tfe, til, miu, Na, T, W, L, ...
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

                %% 5.6 统计指标
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

                %% 5.7 分类
                is_small_all = abs(dV100) < 0.05;

                is_fast_sat = ...
                    abs(dV100) >= 0.05 && ...
                    abs(gradual) < 0.01;

                is_good_prog = ...
                    abs(dV100) >= 0.10 && ...
                    abs(gradual) >= 0.01;

                %% 5.8 case_id
                case_id = sprintf('A%0.2f_PW%0.2eus_Em%0.2f_Es%0.2f', ...
                    amp, pw*1e6, Ea_mean, Ea_sigma);
                case_id = strrep(case_id, '.', 'p');

                %% 5.9 保存单 case
                save(fullfile(outdir, [case_id '.mat']), ...
                    'amp', 'pw', 'delay', 'transit', 'cycle', 'step', ...
                    'Pr', 'tauo', 'alpha', 'bet', 'epife', 'Ndom', ...
                    'VG', 'VD', 'VS', 'Iref', 'VFB_read', ...
                    'tfe', 'til', 'miu', 'Na', 'T', 'W', 'L', ...
                    'Ea_mean', 'Ea_sigma', 'r_Ea', 'r_voff', ...
                    'time', 'volt', 'index', 'index_pre', 'index_end', ...
                    'vfev', 'Stsum', 'P0', 'Vth0', 'ID0', ...
                    'idmin0', 'idmax0', 'valid0', ...
                    'Nread', 'Pread_case', 'Vth_case', 'DeltaVth_case', ...
                    'idmin_case', 'idmax_case', 'valid_case', ...
                    'dV1', 'dV100', 'gradual', 'dP1', 'dP100', ...
                    'valid_ratio', 'is_small_all', 'is_fast_sat', 'is_good_prog', ...
                    'telapsed');

                %% 5.10 汇总
                Result(icase).case_id       = case_id;
                Result(icase).amp           = amp;
                Result(icase).pw            = pw;
                Result(icase).Ea_mean       = Ea_mean;
                Result(icase).Ea_sigma      = Ea_sigma;

                Result(icase).P0            = P0;
                Result(icase).Vth0          = Vth0;
                Result(icase).idmin0        = idmin0;
                Result(icase).idmax0        = idmax0;
                Result(icase).valid0        = valid0;

                Result(icase).dV1           = dV1;
                Result(icase).dV100         = dV100;
                Result(icase).gradual       = gradual;
                Result(icase).dP1           = dP1;
                Result(icase).dP100         = dP100;
                Result(icase).valid_ratio   = valid_ratio;

                Result(icase).is_small_all  = is_small_all;
                Result(icase).is_fast_sat   = is_fast_sat;
                Result(icase).is_good_prog  = is_good_prog;

                Result(icase).elapsed_s     = telapsed;

                fprintf('P0 = %.6e, Vth0 = %.6f V\n', P0, Vth0);
                fprintf('dV1 = %.6f V, dV100 = %.6f V, gradual = %.6f V\n', ...
                    dV1, dV100, gradual);
                fprintf('dP1 = %.6e, dP100 = %.6e, valid_ratio = %.3f\n', ...
                    dP1, dP100, valid_ratio);
                fprintf('class: small_all=%d, fast_sat=%d, good_prog=%d\n', ...
                    is_small_all, is_fast_sat, is_good_prog);
                fprintf('elapsed = %.3f s\n\n', telapsed);
            end
        end
    end
end

%% =========================
% 6) 汇总表
%% =========================
summary_table = struct2table(Result);

save(fullfile(outdir, 'summary_table_stage1_distribution_fine.mat'), ...
    'summary_table', 'Result', 'Nread');

writetable(summary_table, ...
    fullfile(outdir, 'summary_table_stage1_distribution_fine.csv'));

disp('================ summary_table ================');
disp(summary_table(:, { ...
    'case_id','amp','pw','Ea_mean','Ea_sigma', ...
    'dV1','dV100','gradual','dP100','valid_ratio', ...
    'is_small_all','is_fast_sat','is_good_prog'}));

%% =========================
% 7) Top-20: 按 |dV100| 排序
%% =========================
[~, idx_sort_dV100] = sort(abs(summary_table.dV100), 'descend');
topk_dV100 = idx_sort_dV100(1:min(20, height(summary_table)));
top20_by_dV100 = summary_table(topk_dV100, { ...
    'case_id','amp','pw','Ea_mean','Ea_sigma', ...
    'dV1','dV100','gradual','dP100','valid_ratio', ...
    'is_small_all','is_fast_sat','is_good_prog'});

disp('================ Top-20 by |dV100| ================');
disp(top20_by_dV100);

writetable(top20_by_dV100, fullfile(outdir, 'top20_by_dV100.csv'));

%% =========================
% 8) Top-20: 按 gradual 排序
%% =========================
[~, idx_sort_gradual] = sort(summary_table.gradual, 'descend');
topk_gradual = idx_sort_gradual(1:min(20, height(summary_table)));
top20_by_gradual = summary_table(topk_gradual, { ...
    'case_id','amp','pw','Ea_mean','Ea_sigma', ...
    'dV1','dV100','gradual','dP100','valid_ratio', ...
    'is_small_all','is_fast_sat','is_good_prog'});

disp('================ Top-20 by gradual ================');
disp(top20_by_gradual);

writetable(top20_by_gradual, fullfile(outdir, 'top20_by_gradual.csv'));

%% =========================
% 9) 分类统计
%% =========================
n_small_all = sum(summary_table.is_small_all);
n_fast_sat  = sum(summary_table.is_fast_sat);
n_good_prog = sum(summary_table.is_good_prog);

fprintf('\n================ Class counts ================\n');
fprintf('small_all = %d / %d\n', n_small_all, height(summary_table));
fprintf('fast_sat  = %d / %d\n', n_fast_sat,  height(summary_table));
fprintf('good_prog = %d / %d\n', n_good_prog, height(summary_table));

%% =========================
% 10) 图 1: |dV100|
%% =========================
figure;
plot(abs(summary_table.dV100), 'o-', 'LineWidth', 1.2);
xlabel('Case index');
ylabel('|dV100| (V)');
title('Stage 1: |DeltaVth(N=100)| across cases');
grid on;

%% =========================
% 11) 图 2: gradual
%% =========================
figure;
plot(summary_table.gradual, 's-', 'LineWidth', 1.2);
xlabel('Case index');
ylabel('gradual = |dV100 - dV1| (V)');
title('Stage 1: gradual metric across cases');
grid on;

%% =========================
% 12) 图 3: |dV1| vs |dV100|
%% =========================
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
% 13) 图 4: |dP100| vs |dV100|
%% =========================
figure;
scatter(abs(summary_table.dP100), abs(summary_table.dV100), 35, 'filled');
xlabel('|dP100|');
ylabel('|dV100| (V)');
title('|dP100| vs |dV100|');
grid on;

%% =========================
% 14) 图 5: 按类别着色的 |dV1| vs |dV100|
%% =========================
figure;
hold on;

idx_small = find(summary_table.is_small_all);
idx_fast  = find(summary_table.is_fast_sat);
idx_good  = find(summary_table.is_good_prog);

if ~isempty(idx_small)
    scatter(abs(summary_table.dV1(idx_small)), abs(summary_table.dV100(idx_small)), ...
        40, 'o', 'DisplayName', 'small\_all');
end
if ~isempty(idx_fast)
    scatter(abs(summary_table.dV1(idx_fast)), abs(summary_table.dV100(idx_fast)), ...
        40, 's', 'DisplayName', 'fast\_sat');
end
if ~isempty(idx_good)
    scatter(abs(summary_table.dV1(idx_good)), abs(summary_table.dV100(idx_good)), ...
        40, 'd', 'DisplayName', 'good\_prog');
end

mx = max([abs(summary_table.dV1); abs(summary_table.dV100)]) * 1.05;
if isempty(mx) || ~isfinite(mx) || mx <= 0
    mx = 1;
end
plot([0 mx], [0 mx], 'k--', 'LineWidth', 1.2, 'DisplayName', 'y=x');

xlabel('|dV1| (V)');
ylabel('|dV100| (V)');
title('Classification on |dV1| vs |dV100|');
grid on;
axis equal;
xlim([0 mx]);
ylim([0 mx]);
legend('Location', 'best');
hold off;

%% =========================
% 15) 保存代表性 case 图
%     选三个：
%       a) |dV100| 最大
%       b) gradual 最大
%       c) good_prog 中 |dV100| 最大（如果有）
%% =========================
rep_idx = [];

if ~isempty(idx_sort_dV100)
    rep_idx(end+1) = idx_sort_dV100(1);
end
if ~isempty(idx_sort_gradual)
    rep_idx(end+1) = idx_sort_gradual(1);
end

if ~isempty(idx_good)
    [~, ig] = max(abs(summary_table.dV100(idx_good)));
    rep_idx(end+1) = idx_good(ig);
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
    title(['Representative case: ' strrep(case_id,'_','\_')]);
    grid on;

    figure;
    plot(S.Nread, S.Pread_case, 's-', 'LineWidth', 1.5);
    xlabel('Cumulative pulse count N');
    ylabel('P_{read}');
    title(['Representative case P_{read}: ' strrep(case_id,'_','\_')]);
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