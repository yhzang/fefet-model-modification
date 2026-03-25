%% =========================================================
% sweep_delay_3candidates.m
%
% 目的：
%   对 3 组已筛出的“渐进累积” candidate 做 delay sweep，
%   判断原始 Scalable-FeFET 模型对 pulse interval time 是否敏感。
%
% 依赖：
%   wfdef_acc_fix.m
%   get_FE_state.m
%   get_ID.m
%
% 输出：
%   delay_3candidates_out/
%     summary_all_cases.csv
%     summary_by_candidate.csv
%     每个 case 的 .mat
%
% 说明：
%   这版先固定 r_voff = 0，目的是看“原始主干模型”本身的 delay 敏感性。
%% =========================================================

clear;
clc;
close all;

rng(0);

%% =========================
% 0) 输出目录
% =========================
outdir = 'delay_3candidates_out';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

%% =========================
% 1) 固定参数
% =========================
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

% 读点
Nread = [1 2 3 4 5 10 20 50 100].';
Nread = unique(Nread);

%% =========================
% 2) delay 扫描列表
% =========================
delay_list = [1e-7 3e-7 1e-6 3e-6 1e-5 3e-5 1e-4];

%% =========================
% 3) 三组 candidate
% =========================
Cand = struct([]);

% candidate 1
Cand(1).name     = 'cand1_A2p20_PW0p20us_Em1p00_Es0p20';
Cand(1).amp      = 2.2;
Cand(1).pw       = 0.2e-6;
Cand(1).Ea_mean  = 1.00;
Cand(1).Ea_sigma = 0.20;

% candidate 2
Cand(2).name     = 'cand2_A2p20_PW0p50us_Em1p00_Es0p20';
Cand(2).amp      = 2.2;
Cand(2).pw       = 0.5e-6;
Cand(2).Ea_mean  = 1.00;
Cand(2).Ea_sigma = 0.20;

% candidate 3
Cand(3).name     = 'cand3_A2p20_PW0p20us_Em0p95_Es0p20';
Cand(3).amp      = 2.2;
Cand(3).pw       = 0.2e-6;
Cand(3).Ea_mean  = 0.95;
Cand(3).Ea_sigma = 0.20;

%% =========================
% 4) 固定随机样本
%    不同 delay 下只改 delay，不改随机域样本
% =========================
Weight  = ones(Ndom,1) / Ndom;
St_init = -ones(Ndom,1);

base_randn_Ea = randn(Ndom,1);

% 原始主干 baseline：先关掉 voff
r_voff = zeros(Ndom,1);

%% =========================
% 5) 结果容器
% =========================
Result = struct([]);
icase = 0;

nCand  = numel(Cand);
nDelay = numel(delay_list);
nTotal = nCand * nDelay;

fprintf('================ Run setup ================\n');
fprintf('Total cases = %d\n', nTotal);
fprintf('cycle = %d, step = %g s, transit = %g s\n', cycle, step, transit);
fprintf('Iref = %g A, VFB_read = %g V\n\n', Iref, VFB_read);

%% =========================
% 6) 主循环
% =========================
for ic = 1:nCand
    cand_name = Cand(ic).name;
    amp       = Cand(ic).amp;
    pw        = Cand(ic).pw;
    Ea_mean   = Cand(ic).Ea_mean;
    Ea_sigma  = Cand(ic).Ea_sigma;

    r_Ea = Ea_mean + Ea_sigma * base_randn_Ea;
    r_Ea(r_Ea <= 0) = Ea_mean;

    for id = 1:nDelay
        delay = delay_list(id);

        icase = icase + 1;
        fprintf('=====================================================\n');
        fprintf('Case %d / %d\n', icase, nTotal);
        fprintf('%s\n', cand_name);
        fprintf('amp=%.3f V, pw=%g s, Ea_mean=%.3f, Ea_sigma=%.3f, delay=%g s\n', ...
            amp, pw, Ea_mean, Ea_sigma, delay);
        fprintf('=====================================================\n');

        %% 6.1 生成波形
        [time, volt, index, index_pre, index_end] = ...
            wfdef_acc_fix(amp, pw, step, delay, transit, cycle);

        index     = index(:);
        index_pre = index_pre(:);
        index_end = index_end(:);

        %% 6.2 FE 状态演化
        tic;
        [vfev, Stsum] = get_FE_state(time, volt, St_init, Weight, r_Ea, r_voff, ...
            Pr, tauo, alpha, bet, epife, Ndom);
        telapsed = toc;

        if isvector(Stsum)
            st_trace = Stsum(:);
        else
            st_trace = Stsum(:,1);
        end

        %% 6.3 N=0 参考点
        idx0 = index_pre(1);
        P0 = st_trace(idx0);

        [Vth0, ID0, idmin0, idmax0, valid0] = local_extract_vth( ...
            P0, epife, tfe, til, miu, Na, T, W, L, ...
            VG, VFB_read, VD, VS, Iref);

        %% 6.4 各读点
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

        %% 6.5 统计指标
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

        %% 6.6 保存 case
        case_id = sprintf('%s_delay_%0.1e', cand_name, delay);
        case_id = strrep(case_id, '.', 'p');
        case_id = strrep(case_id, '-', 'm');

        save(fullfile(outdir, [case_id '.mat']), ...
            'cand_name', 'amp', 'pw', 'Ea_mean', 'Ea_sigma', ...
            'delay', 'transit', 'cycle', 'step', ...
            'Pr', 'tauo', 'alpha', 'bet', 'epife', 'Ndom', ...
            'VG', 'VD', 'VS', 'Iref', 'VFB_read', ...
            'tfe', 'til', 'miu', 'Na', 'T', 'W', 'L', ...
            'r_Ea', 'r_voff', 'time', 'volt', 'index', 'index_pre', 'index_end', ...
            'vfev', 'Stsum', 'P0', 'Vth0', 'ID0', 'idmin0', 'idmax0', 'valid0', ...
            'Nread', 'Pread_case', 'Vth_case', 'DeltaVth_case', ...
            'idmin_case', 'idmax_case', 'valid_case', ...
            'dV1', 'dV100', 'gradual', 'dP1', 'dP100', 'valid_ratio', 'telapsed');

        %% 6.7 汇总
        Result(icase).case_id      = case_id;
        Result(icase).cand_name    = cand_name;
        Result(icase).amp          = amp;
        Result(icase).pw           = pw;
        Result(icase).Ea_mean      = Ea_mean;
        Result(icase).Ea_sigma     = Ea_sigma;
        Result(icase).delay        = delay;

        Result(icase).P0           = P0;
        Result(icase).Vth0         = Vth0;
        Result(icase).dV1          = dV1;
        Result(icase).dV100        = dV100;
        Result(icase).gradual      = gradual;
        Result(icase).dP1          = dP1;
        Result(icase).dP100        = dP100;
        Result(icase).valid_ratio  = valid_ratio;
        Result(icase).elapsed_s    = telapsed;

        fprintf('dV1 = %.6f V, dV100 = %.6f V, gradual = %.6f V\n', ...
            dV1, dV100, gradual);
        fprintf('dP1 = %.6e, dP100 = %.6e, valid_ratio = %.3f\n', ...
            dP1, dP100, valid_ratio);
        fprintf('elapsed = %.3f s\n\n', telapsed);
    end
end

%% =========================
% 7) 全部结果表
% =========================
summary_all_cases = struct2table(Result);

save(fullfile(outdir, 'summary_all_cases.mat'), ...
    'summary_all_cases', 'Result', 'Nread', 'Cand', 'delay_list');

writetable(summary_all_cases, fullfile(outdir, 'summary_all_cases.csv'));

disp('================ summary_all_cases ================');
disp(summary_all_cases(:, { ...
    'cand_name','delay','dV1','dV100','gradual','dP1','dP100','valid_ratio'}));

%% =========================
% 8) 按 candidate 汇总 span
% =========================
SummaryCand = struct([]);

for ic = 1:nCand
    cand_name = Cand(ic).name;
    idxc = strcmp(summary_all_cases.cand_name, cand_name);

    Tsub = summary_all_cases(idxc, :);

    dV1_vec    = Tsub.dV1;
    dV100_vec  = Tsub.dV100;
    gradual_vec= Tsub.gradual;
    dP100_vec  = Tsub.dP100;

    span_dV1    = max(dV1_vec)   - min(dV1_vec);
    span_dV100  = max(dV100_vec) - min(dV100_vec);
    span_grad   = max(gradual_vec) - min(gradual_vec);
    span_dP100  = max(dP100_vec) - min(dP100_vec);

    [~, i_best_dV100] = max(abs(dV100_vec));
    [~, i_best_grad]  = max(gradual_vec);

    SummaryCand(ic).cand_name          = cand_name;
    SummaryCand(ic).amp                = Cand(ic).amp;
    SummaryCand(ic).pw                 = Cand(ic).pw;
    SummaryCand(ic).Ea_mean            = Cand(ic).Ea_mean;
    SummaryCand(ic).Ea_sigma           = Cand(ic).Ea_sigma;

    SummaryCand(ic).span_dV1           = span_dV1;
    SummaryCand(ic).span_dV100         = span_dV100;
    SummaryCand(ic).span_gradual       = span_grad;
    SummaryCand(ic).span_dP100         = span_dP100;

    SummaryCand(ic).best_abs_dV100     = dV100_vec(i_best_dV100);
    SummaryCand(ic).best_abs_dV100_delay = Tsub.delay(i_best_dV100);

    SummaryCand(ic).best_gradual       = gradual_vec(i_best_grad);
    SummaryCand(ic).best_gradual_delay = Tsub.delay(i_best_grad);
end

summary_by_candidate = struct2table(SummaryCand);

save(fullfile(outdir, 'summary_by_candidate.mat'), 'summary_by_candidate');
writetable(summary_by_candidate, fullfile(outdir, 'summary_by_candidate.csv'));

disp('================ summary_by_candidate ================');
disp(summary_by_candidate);

%% =========================
% 9) 画图：每个 candidate 各两张
%    图1：DeltaVth(N) across delay
%    图2：Pread(N) across delay
% =========================
for ic = 1:nCand
    cand_name = Cand(ic).name;
    idxc = strcmp(summary_all_cases.cand_name, cand_name);
    Tsub = summary_all_cases(idxc, :);

    % 为了按 delay 从小到大画
    [~, order_idx] = sort(Tsub.delay, 'ascend');
    Tsub = Tsub(order_idx, :);

    figure;
    hold on;
    legtxt = cell(height(Tsub),1);

    for ii = 1:height(Tsub)
        S = load(fullfile(outdir, [Tsub.case_id{ii} '.mat']));
        plot(S.Nread, S.DeltaVth_case, 'o-', 'LineWidth', 1.4);
        legtxt{ii} = sprintf('delay=%0.1e s', S.delay);
    end
    xlabel('Cumulative pulse count N');
    ylabel('\DeltaV_{th} (V)');
    title(['Delay sweep: ' strrep(cand_name, '_', '\_')]);
    grid on;
    legend(legtxt, 'Location', 'best');
    hold off;

    figure;
    hold on;
    for ii = 1:height(Tsub)
        S = load(fullfile(outdir, [Tsub.case_id{ii} '.mat']));
        plot(S.Nread, S.Pread_case, 's-', 'LineWidth', 1.4);
    end
    xlabel('Cumulative pulse count N');
    ylabel('P_{read}');
    title(['Delay sweep P_{read}: ' strrep(cand_name, '_', '\_')]);
    grid on;
    legend(legtxt, 'Location', 'best');
    hold off;
end

%% =========================
% 10) 画图：candidate 间对比 span
% =========================
figure;
bar(summary_by_candidate.span_dV1);
set(gca, 'XTickLabel', summary_by_candidate.cand_name, 'XTickLabelRotation', 20);
ylabel('span of dV1 across delay (V)');
title('Delay sensitivity: span(dV1)');
grid on;

figure;
bar(summary_by_candidate.span_dV100);
set(gca, 'XTickLabel', summary_by_candidate.cand_name, 'XTickLabelRotation', 20);
ylabel('span of dV100 across delay (V)');
title('Delay sensitivity: span(dV100)');
grid on;

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