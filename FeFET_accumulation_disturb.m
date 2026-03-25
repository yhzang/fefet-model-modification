%% =========================================================
%  FeFET_accumulation_disturb.m
%  用于写扰 / 累积脉冲仿真
%  依赖：
%    wfdef_acc_fix.m
%    get_FE_state.m
% =========================================================

clear;
clc;
close all;

rng(0);

%% =========================
% 1) 扰动脉冲参数
% =========================
VDD     = 5;

% 二选一：
amp     = 0.5 * VDD;      % 2.5 V，真正的 V/2 disturb
% amp   = 0.13 * VDD;     % 0.65 V，较弱扰动

pw      = 0.1e-6;
delay   = 1e-6;
transit = 0.02e-6;

% 先不要直接 1e6
cycle   = 1e3;

step    = 0.01e-6;

fprintf('=== 扰动脉冲设置 ===\n');
fprintf('VDD     = %g V\n', VDD);
fprintf('amp     = %g V\n', amp);
fprintf('pw      = %g s\n', pw);
fprintf('delay   = %g s\n', delay);
fprintf('transit = %g s\n', transit);
fprintf('cycle   = %d\n', cycle);
fprintf('step    = %g s\n\n', step);

%% =========================
% 2) 铁电模型参数
% =========================
Pr    = 25;
tauo  = 1e-9;
alpha = 2;
bet   = 2;
epife = 28;

Ndom = 10000;

%% =========================
% 3) 随机畴参数（列向量）
% =========================
Weight  = ones(Ndom,1) / Ndom;
St_init = -ones(Ndom,1);

Ea_mean  = 1.0;
Ea_sigma = 0.5;
r_Ea = Ea_mean + Ea_sigma * randn(Ndom,1);
r_Ea(r_Ea <= 0) = Ea_mean;

% 第一版先关掉 offset，避免零偏压间隔也继续演化
r_voff = zeros(Ndom,1);

%% =========================
% 4) 生成扰动脉冲波形
% =========================
[time, volt, index, index_pre, index_end] = ...
    wfdef_acc_fix(amp, pw, step, delay, transit, cycle);

index     = index(:);
index_pre = index_pre(:);
index_end = index_end(:);

fprintf('波形总点数 = %d\n', length(time));
fprintf('index(read) 点数 = %d\n', length(index));
fprintf('index_pre   点数 = %d\n', length(index_pre));
fprintf('index_end   点数 = %d\n\n', length(index_end));

%% =========================
% 5) 计算铁电状态演化
% =========================
fprintf('开始运行 get_FE_state ...\n');
tic;
[vfev, Stsum] = get_FE_state(time, volt, St_init, Weight, r_Ea, r_voff, ...
                             Pr, tauo, alpha, bet, epife, Ndom);
toc;
fprintf('get_FE_state 完成。\n\n');

%% =========================
% 6) 自动生成累计脉冲读出点 Nread
% =========================
Npulse_total = length(index);

maxPow = floor(log10(Npulse_total));
Nread = 10.^(0:maxPow);

if Nread(end) ~= Npulse_total
    Nread = unique([Nread; Npulse_total]);
end
Nread = Nread(:);

fprintf('总脉冲数 Npulse_total = %d\n', Npulse_total);
fprintf('自动生成读出点 Nread = \n');
disp(Nread.');

%% =========================
% 7) 画脉冲波形与三种索引位置
% =========================
figure;
plot(time, volt, 'b-'); hold on;

idx_pre_plot  = index_pre(index_pre >= 1 & index_pre <= length(time));
idx_end_plot  = index_end(index_end >= 1 & index_end <= length(time));
idx_read_plot = index(index >= 1 & index <= length(time));

plot(time(idx_pre_plot),  volt(idx_pre_plot),  'ks', 'MarkerSize', 6, 'LineWidth', 1.0);
plot(time(idx_end_plot),  volt(idx_end_plot),  'md', 'MarkerSize', 6, 'LineWidth', 1.0);
plot(time(idx_read_plot), volt(idx_read_plot), 'ro', 'MarkerSize', 6, 'LineWidth', 1.0);

xlabel('time (s)');
ylabel('V_G (V)');
title('脉冲波形与 index\_pre / index\_end / index(read)');
legend('waveform', 'index\_pre', 'index\_end', 'index(read)', 'Location', 'best');
grid on;
hold off;

%% =========================
% 8) 保存结果
% =========================
save('disturb_run_workspace.mat', ...
    'time', 'volt', 'index', 'index_pre', 'index_end', ...
    'Stsum', 'vfev', ...
    'Npulse_total', 'Nread', ...
    'VDD', 'amp', 'pw', 'delay', 'transit', 'cycle', 'step', ...
    'Pr', 'tauo', 'alpha', 'bet', 'epife', ...
    'Weight', 'St_init', 'r_Ea', 'r_voff', 'Ndom');

fprintf('\n已完成 FeFET_accumulation_disturb.m\n');
fprintf('下一步请运行：extract_Vth_disturb\n');