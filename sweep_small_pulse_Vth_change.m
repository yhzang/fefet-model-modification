%% =========================================================
% sweep_small_pulse_Vth_change.m
%
% 目的：
%   使用小电压脉冲模拟 `Vth` 的变化趋势，
%   在前期变化较小，后期逐渐增加，最终接近饱和。
%
% 依赖：
%   wfdef_acc_fix.m
%   get_FE_state.m
%   get_ID.m
%
% 输出：
%   small_pulse_Vth_out/
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
% 1) 固定脉冲参数
% =========================
total_pulses = 1e6;  % 总脉冲次数

amp      = 2.2;      % 脉冲幅度（V）
pw       = 0.1e-6;   % 脉冲宽度（秒）
Vdisturb = 3.5;      % 电压扰动（V）

% 脉冲间隔可以保持一致（简化模型，避免过多调整）
T_pulse = pw;        % 固定脉冲宽度为 0.1 μs
T_interval = 1e-6;    % 以 1 μs 为脉冲间隔（保持一致）

%% =========================
% 2) 固定 FE 参数
% =========================
tfe_list        = [1.5e-6 1.8e-6 2.0e-6];  % FE 厚度（tfe）调节
epife_read_list = [16 18 20];                % 读出区域的 FE 常数（epife_read）

% 读出点
Nread = [1 10 100].';

% 固定 MOS / stack 参数
til = 1e-7;
miu = 50;
Na  = 3e17;
T0  = 300;    % 温度
W   = 1;      % 宽度（微米）
L   = 1;      % 长度（微米）
VD  = 0.05;   % 漏极电压（V）
VS  = 0;      % 源极电压（V）

% FE 参数
Pr    = 25;   % 极化强度
tauo  = 1e-9; % 退火时间常数
alpha = 2;    % 非线性参数
bet   = 2;    % 非线性参数

% 默认 FE eps（写入主干）
epife_write = 28;

Ndom = 10000; % 域数

Weight  = ones(Ndom,1) / Ndom;
St_init = -ones(Ndom,1);

%% =========================
% 3) 设定参数和输出目录
% =========================
outdir = 'small_pulse_Vth_out';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

diary(fullfile(outdir, 'run_log.txt'));
fprintf('=== sweep_small_pulse_Vth_change ===\n');

%% =========================
% 4) 生成case列表
% =========================
case_list = {};
icase = 0;

for itfe = 1:numel(tfe_list)
    for iepi = 1:numel(epife_read_list)
        icase = icase + 1;

        cfg = struct();
        cfg.case_idx    = icase;
        cfg.amp         = amp;
        cfg.pw          = pw;
        cfg.Vdisturb    = Vdisturb;
        cfg.tfe         = tfe_list(itfe);
        cfg.epife_read  = epife_read_list(iepi);

        cfg.case_id = make_case_id(cfg);

        case_list{end+1,1} = cfg;
    end
end

nCase = numel(case_list);
fprintf('Total cases = %d\n', nCase);

%% =========================
% 5) 预分配结果
% =========================
results = repmat(init_result_struct(), nCase, 1);

%% =========================
% 6) 主循环
% =========================
for ic = 1:nCase
    cfg = case_list{ic};

    fprintf('\n==================================================\n');
    fprintf('[%d / %d] %s\n', ic, nCase, cfg.case_id);
    fprintf('==================================================\n');

    try
        %% 6.1 当前 case 的 r_Ea
        r_Ea = cfg.Vdisturb + 0.2 * randn(Ndom,1);
        r_Ea(r_Ea <= 0) = cfg.Vdisturb;

        %% 6.2 波形生成
        [time, volt, index, index_pre, index_end] = ...
            wfdef_acc_fix(cfg.amp, cfg.pw, T_interval, T_interval, T_interval, total_pulses);

        %% 6.3 FE 状态推进
        tic;
        [vfev, Stsum] = get_FE_state(time, volt, St_init, Weight, r_Ea, r_voff, ...
                                     Pr, tauo, alpha, bet, epife_write, Ndom);
        telapsed = toc;
        fprintf('get_FE_state done, elapsed = %.3f s\n', telapsed);

        %% 6.4 取单条轨迹
        if isvector(Stsum)
            st_trace = Stsum(:);
        else
            st_trace = Stsum(:,1);
        end

        %% 6.5 计算阈值
        VG = make_VG_from_name(cfg.VG_name);
        idx0 = index_pre(1);
        P0 = st_trace(idx0);
        [Vth0, ~, ~, ~, valid0] = local_extract_vth( ...
            P0, cfg.epife_read, cfg.tfe, til, miu, Na, T0, W, L, ...
            VG, cfg.VFB_read, VD, VS, cfg.Iref);

        %% 6.6 提取 N=1,10,30,100
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
                    Pread_vec(kk), cfg.epife_read, cfg.tfe, til, miu, Na, T0, W, L, ...
                    VG, cfg.VFB_read, VD, VS, cfg.Iref);

            if valid0 && valid_vec(kk)
                DeltaVth_vec(kk) = Vth_vec(kk) - Vth0;
            else
                DeltaVth_vec(kk) = NaN;
            end
        end

        %% 6.7 保存结果
        results(ic).case_id      = string(cfg.case_id);
        results(ic).amp          = cfg.amp;
        results(ic).pw           = cfg.pw;
        results(ic).Vdisturb     = cfg.Vdisturb;
        results(ic).tfe          = cfg.tfe;
        results(ic).epife_read   = cfg.epife_read;

        results(ic).dV1          = DeltaVth_vec(1);
        results(ic).dV10         = DeltaVth_vec(2);
        results(ic).dV100        = DeltaVth_vec(3);

        results(ic).A1           = abs(DeltaVth_vec(1));
        results(ic).A10          = abs(DeltaVth_vec(2));
        results(ic).A100         = abs(DeltaVth_vec(3));

        results(ic).gradual      = results(ic).A100 - results(ic).A1;

        results(ic).valid0       = valid0;

    catch ME
        warning('Case failed: %s', cfg.case_id);
        warning('%s', ME.message);

        results(ic).case_id      = string(cfg.case_id);
        results(ic).is_invalid   = true;
    end
end

%% =========================
% 7) 保存结果
% =========================
T = struct2table(results);
writetable(T, fullfile(outdir, 'all_cases.csv'));

T_dV100 = sortrows(T, 'A100', 'descend');
writetable(T_dV100(1:min(20,height(T_dV100)), :), ...
    fullfile(outdir, 'top_by_dV100.csv'));

T_gradual = sortrows(T, 'gradual', 'descend');
writetable(T_gradual(1:min(20,height(T_gradual)), :), ...
    fullfile(outdir, 'top_by_gradual.csv'));

disp('================ Top-20 by |dV100| ================');
disp(T_dV100(1:min(20,height(T_dV100)), :));

disp('================ Top-20 by gradual ================');
disp(T_gradual(1:min(20,height(T_gradual)), :));

%% =========================
% 8) 画图
% =========================
figure;
plot(T.A1, T.A100, 'o-', 'DisplayName', '|dV1| vs |dV100|');
xlabel('|dV1|');
ylabel('|dV100|');
title('Comparison of |dV1| vs |dV100|');
grid on;
legend('show');
saveas(gcf, fullfile(outdir, 'dV1_vs_dV100.png'));

%% =========================
% 9) 保存工作区
% =========================
save(fullfile(outdir, 'workspace_results.mat'), 'T', 'results', 'case_list');

fprintf('All done. Results saved in folder: %s\n', outdir);
diary off;