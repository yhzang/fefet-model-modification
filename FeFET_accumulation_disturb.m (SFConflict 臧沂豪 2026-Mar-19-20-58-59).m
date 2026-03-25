%% =========================================================
%  FeFET_accumulation_disturb.m
%  用于写扰 / 累积脉冲仿真
%  依赖：
%    wfdef_acc.m
%    get_FE_state.m
% =========================================================

clear;
clc;
close all;

rng(0);

%% =========================
% 1) 扰动脉冲参数
% =========================
VDD     = 5;          % 你的设定
amp     = 0.15 * VDD;  % 1/2 VDD = 2.5 V
pw      = 1e-6;      % 脉宽 10 us
delay   = 1e-6;       % 脉冲间隔 1 us
transit = 0.1e-6;     % 上升/下降时间 0.1 us

% !!! 建议先从 1e3 开始，再试 1e4、1e5，最后再冲 1e6 !!!
cycle   = 1e3;        % 总脉冲数

% 时间步不要太细，否则波形点数会非常大
step    = 0.5e-6;     % 0.5 us

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
% 这些值先沿用你当前工程附近的常见设置；
% 如果你原始 FeFET_accumulation.m 里有更明确的参数，请以后再对齐
Pr    = 25;
tauo  = 1e-9;
alpha = 2;
bet   = 2;
epife = 28;

Ndom = 10000;   % 畴数

%% =========================
% 3) 随机畴参数（关键：全部列向量）
% =========================
Weight  = ones(Ndom,1) / Ndom;

St_init = -ones(Ndom,1);   % 初始全向下

Ea_mean  = 1.0;
Ea_sigma = 0.2;
r_Ea = Ea_mean + Ea_sigma * randn(Ndom,1);
r_Ea(r_Ea <= 0) = Ea_mean;

voff_mean  = 0.0;
voff_sigma = 0.2;
r_voff = voff_mean + voff_sigma * randn(Ndom,1);

%% =========================
% 4) 生成扰动脉冲波形
% =========================
% 默认接口：
% [time, volt, index] = wfdef_acc(amp, pw, step, delay, transit, cycle)
[time, volt, index] = wfdef_acc(amp, pw, step, delay, transit, cycle);

index = index(:);   % 强制列向量

fprintf('波形总点数 = %d\n', length(time));
fprintf('index 点数 = %d\n\n', length(index));

%% =========================
% 5) 计算铁电状态演化
% =========================
% 默认接口：
% [vfev, Stsum] = get_FE_state(time, volt, St_init, Weight, r_Ea, r_voff, ...
%                              Pr, tauo, alpha, bet, epife, Ndom)
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
% 7) 画脉冲波形与 index 位置
% =========================
figure;
plot(time, volt, 'b-'); hold on;
idx_plot = index(index >= 1 & index <= length(time));
plot(time(idx_plot), volt(idx_plot), 'r.', 'MarkerSize', 8);
xlabel('time (s)');
ylabel('V_G (V)');
title('脉冲波形与 index 标记位置');
grid on;
hold off;

%% =========================
% 8) 保存结果（可选）
% =========================
save('disturb_run_workspace.mat', ...
    'time', 'volt', 'index', 'Stsum', 'vfev', ...
    'Npulse_total', 'Nread', ...
    'VDD', 'amp', 'pw', 'delay', 'transit', 'cycle', 'step', ...
    'Pr', 'tauo', 'alpha', 'bet', 'epife', ...
    'Weight', 'St_init', 'r_Ea', 'r_voff', 'Ndom');

fprintf('\n已完成 FeFET_accumulation_disturb.m\n');
fprintf('下一步请运行：extract_Vth_shift_log\n');