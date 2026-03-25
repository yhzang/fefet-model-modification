%% =========================================================
% sweep_stage1_find_big_dvth_window.m
%
% 目的：
%   第一阶段，先找“能把 DeltaVth 做大”的 FE 工作窗口。
%
% 扫描参数：
%   amp, pw, Ea_mean, Ea_sigma
%
% 固定参数：
%   tfe, epife, Iref, VFB_read 等读出参数先固定
%
% 输出：
%   1) 每个 case 的 MAT
%   2) summary_table.csv
%   3) 自动筛 candidate
%
% 依赖：
%   wfdef_acc_fix.m
%   get_FE_state.m
%   get_ID.m
%% =========================================================

clear;
clc;
close all;

rng(0);

%% =========================
% 0) 输出目录
% =========================
outdir = 'stage1_big_dvth_window';
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

% 早期加密读点
Nread = [1 2 3 4 5 10 20 50 100].';
Nread = unique(Nread);

%% =========================
% 2) 扫描列表
% =========================
amp_list      = [1.8 2.2 2.6 3.0];
pw_list       = [0.2e-6 0.5e-6 1e-6 2e-6];
Ea_mean_list  = [0.5 0.7 0.9];
Ea_sigma_list = [0.4 0.6 0.8];
%% =========================
% 3) 固定随机样本
%    注意：
%    为了不同 case 可比，建议固定 base randn
%% =========================
Weight  = ones(Ndom,1) / Ndom;
St_init = -ones(Ndom,1);

base_randn_Ea = randn(Ndom,1);

% 先关掉 voff，避免 gap 段附加演化
r_voff = zeros(Ndom,1);

%% =========================
% 4) 预分配结果
% =========================
Result = struct([]);

icase = 0;
nTotal = numel(amp_list) * numel(pw_list) * numel(Ea_mean_list) * numel(Ea_sigma_list);

fprintf('Total cases = %d\n', nTotal);

%% =========================
% 5) 主循环
% =========================
for ia = 1:numel(amp_list)
    amp = amp_list(ia);

    for ipw = 1:numel(pw_list)
        pw = pw_list(ipw);

        for im = 1:numel(Ea_mean_list)
            Ea_mean = Ea_mean_list(im);

            for is = 1:numel(Ea_sigma_list)
                Ea_sigma = Ea_sigma_list(is);

                icase = icase + 1;
                fprintf('\n=====================================================\n');
                fprintf('Case %d / %d\n', icase, nTotal);
                fprintf('amp=%.3f V, pw=%g s, Ea_mean=%.3f, Ea_sigma=%.3f\n', ...
                    amp, pw, Ea_mean, Ea_sigma);
                fprintf('=====================================================\n');

                % ---- 生成 r_Ea ----
                r_Ea = Ea_mean + Ea_sigma * base_randn_Ea;
                r_Ea(r_Ea <= 0) = Ea_mean;

                % ---- 生成波形 ----
                [time, volt, index, index_pre, index_end] = ...
                    wfdef_acc_fix(amp, pw, step, delay, transit, cycle);

                index     = index(:);
                index_pre = index_pre(:);
                index_end = index_end(:);

                % ---- FE 演化 ----
                tic;
                [vfev, Stsum] = get_FE_state(time, volt, St_init, Weight, r_Ea, r_voff, ...
                    Pr, tauo, alpha, bet, epife, Ndom);
                telapsed = toc;

                if isvector(Stsum)
                    st_trace = Stsum(:);
                else
                    st_trace = Stsum(:,1);
                end

                % ---- N=0 参考点 ----
                idx0 = index_pre(1);
                P0 = st_trace(idx0);

                [Vth0, ID0, idmin0, idmax0, valid0] = local_extract_vth( ...
                    P0, epife, tfe, til, miu, Na, T, W, L, ...
                    VG, VFB_read, VD, VS, Iref);

                % ---- 各读点 ----
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

                % ---- 统计量 ----
                idx1    = find(Nread == 1,   1);
                idx100  = find(Nread == 100, 1);

                dV1 = nan;
                dV100 = nan;
                gradual = nan;
                dP100 = nan;

                if ~isempty(idx1)
                    dV1 = DeltaVth_case(idx1);
                end
                if ~isempty(idx100)
                    dV100 = DeltaVth_case(idx100);
                    dP100 = Pread_case(idx100) - P0;
                end
                if ~isempty(idx1) && ~isempty(idx100)
                    gradual = abs(DeltaVth_case(idx100) - DeltaVth_case(idx1));
                end

                valid_ratio = mean(valid_case);

                % ---- 保存 case ----
                case_id = sprintf('A%0.2f_PW%0.2eus_Em%0.2f_Es%0.2f', ...
                    amp, pw*1e6, Ea_mean, Ea_sigma);
                case_id = strrep(case_id, '.', 'p');

                save(fullfile(outdir, [case_id '.mat']), ...
                    'amp', 'pw', 'delay', 'transit', 'cycle', 'step', ...
                    'Pr', 'tauo', 'alpha', 'bet', 'epife', 'Ndom', ...
                    'VG', 'VD', 'VS', 'Iref', 'VFB_read', ...
                    'tfe', 'til', 'miu', 'Na', 'T', 'W', 'L', ...
                    'Ea_mean', 'Ea_sigma', 'r_Ea', 'r_voff', ...
                    'time', 'volt', 'index', 'index_pre', 'index_end', ...
                    'Stsum', 'vfev', 'P0', 'Vth0', 'ID0', 'idmin0', 'idmax0', 'valid0', ...
                    'Nread', 'Pread_case', 'Vth_case', 'DeltaVth_case', ...
                    'idmin_case', 'idmax_case', 'valid_case', ...
                    'dV1', 'dV100', 'gradual', 'dP100', 'valid_ratio', 'telapsed');

                % ---- 汇总结构体 ----
                Result(icase).case_id      = case_id;
                Result(icase).amp          = amp;
                Result(icase).pw           = pw;
                Result(icase).Ea_mean      = Ea_mean;
                Result(icase).Ea_sigma     = Ea_sigma;
                Result(icase).P0           = P0;
                Result(icase).Vth0         = Vth0;
                Result(icase).dV1          = dV1;
                Result(icase).dV100        = dV100;
                Result(icase).gradual      = gradual;
                Result(icase).dP100        = dP100;
                Result(icase).idmin0       = idmin0;
                Result(icase).idmax0       = idmax0;
                Result(icase).valid0       = valid0;
                Result(icase).valid_ratio  = valid_ratio;
                Result(icase).elapsed_s    = telapsed;

                fprintf('P0 = %.6e, Vth0 = %.6f, dV1 = %.6f, dV100 = %.6f, gradual = %.6f\n', ...
                    P0, Vth0, dV1, dV100, gradual);
                fprintf('dP100 = %.6e, valid_ratio = %.3f, elapsed = %.2f s\n', ...
                    dP100, valid_ratio, telapsed);
            end
        end
    end
end

%% =========================
% 6) 汇总表
% =========================
summary_table = struct2table(Result);
writetable(summary_table, fullfile(outdir, 'summary_table_stage1.csv'));
save(fullfile(outdir, 'summary_table_stage1.mat'), 'summary_table', 'Result', 'Nread');

disp('================ summary_table_stage1 ================');
disp(summary_table(:, {'case_id','amp','pw','Ea_mean','Ea_sigma','dV1','dV100','gradual','dP100','valid_ratio'}));

%% =========================
% 7) 自动筛 candidate
%    这里的门槛你可以按目标再改
%% =========================
candidate = find( ...
    abs(summary_table.dV100) > 0.3  & ...
    abs(summary_table.dV100) < 1.0  & ...
    summary_table.gradual   > 0.1  & ...
    summary_table.valid_ratio > 0.8);

fprintf('\n================ candidates (stage 1) ================\n');
if isempty(candidate)
    fprintf('没有找到满足条件的候选窗口。\n');
    fprintf('优先考虑：\n');
    fprintf('1) 提高 amp 上限\n');
    fprintf('2) 增加 pw\n');
    fprintf('3) 降低 Ea_mean\n');
    fprintf('4) 增大 Ea_sigma\n');
else
    disp(summary_table(candidate, {'case_id','amp','pw','Ea_mean','Ea_sigma','dV1','dV100','gradual','dP100','valid_ratio'}));
end

%% =========================
% 8) 快速图：dV100 vs case index
%% =========================
figure;
plot(abs(summary_table.dV100), 'o-');
xlabel('Case index');
ylabel('|dV100| (V)');
title('Stage 1: |DeltaVth(N=100)| across cases');
grid on;

figure;
plot(summary_table.gradual, 's-');
xlabel('Case index');
ylabel('gradual = |dV100 - dV1| (V)');
title('Stage 1: gradual metric across cases');
grid on;

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