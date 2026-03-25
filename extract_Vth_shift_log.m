%% =========================================================
%  extract_Vth_shift_log.m
%  对累计脉冲数 Nread 提取 Vth
%  使用 N=0 初始态作为参考
%  Nread 加密：每个 decade 取 1,2,3,5,7
% =========================================================

clc;

%% =========================
% 1) 检查变量
% =========================
if ~exist('Stsum','var')
    error('未找到 Stsum，请先运行 FeFET_accumulation_disturb.m');
end

if ~exist('index','var')
    error('未找到 index，请先运行 FeFET_accumulation_disturb.m');
end

index = index(:);

if ~exist('Npulse_total','var')
    Npulse_total = length(index);
end

%% =========================
% 2) 更密的读取点 Nread
%    每个 decade 取 1,2,3,5,7
% =========================
basePts = [1 2 3 5 7];
maxPow = floor(log10(Npulse_total));

Nread = [];
for p = 0:maxPow
    Nread = [Nread; basePts(:) * 10^p];
end

Nread = unique(round(Nread));
Nread = Nread(Nread >= 1 & Nread <= Npulse_total);

% 把总脉冲数本身也加进去
if isempty(Nread) || Nread(end) ~= Npulse_total
    Nread = unique([Nread; Npulse_total]);
end

fprintf('总脉冲数 Npulse_total = %d\n', Npulse_total);
fprintf('读取点 Nread = \n');
disp(Nread.');

%% =========================
% 3) 选择样本列 dev_id
% =========================
if isvector(Stsum)
    Stsum_series = Stsum(:);
    dev_id = 1;
else
    dev_id = 1;
    if dev_id > size(Stsum,2)
        error('dev_id 超出 Stsum 列数。');
    end
    Stsum_series = Stsum(:, dev_id);
end

fprintf('当前使用的样本列 dev_id = %d\n', dev_id);

%% =========================
% 4) 读出条件
% =========================
VG   = -0.5:0.005:3.0;
VD   = 0.05;
VS   = 0;
Iref = 1e-7;

epife = 28;
tfe   = 0.8e-6;
til   = 1e-7;
miu   = 50;
Na    = 3e17;
T     = 300;
W     = 1;
L     = 1;
VFB   = 0;

fprintf('\n开始提取 Vth ...\n');
fprintf('Iref = %e A\n', Iref);

%% =========================
% 5) N=0 初始态参考
%    直接取 Stsum 时间序列起点
% =========================
Psum0 = Stsum_series(1);

ID0 = get_ID(Psum0, epife, tfe, til, miu, Na, T, W, L, VG, VFB, VD, VS);
ID0 = ID0(:).';

valid0 = isfinite(ID0) & isfinite(VG);
ID0v = ID0(valid0);
VG0v = VG(valid0);

pos0 = ID0v > 0;
ID0v = ID0v(pos0);
VG0v = VG0v(pos0);

if isempty(ID0v)
    error('N=0 初始态 ID 全部非正，无法提取 Vth0；请增大 VG 扫描范围或降低 Iref');
end

[uid0, uidind0] = unique(ID0v);
VG0_u = VG0v(uidind0);

Imin0 = min(uid0);
Imax0 = max(uid0);

fprintf('N=0, Psum0=%+.6e, ID range=[%e, %e]\n', Psum0, Imin0, Imax0);

if Iref < Imin0 || Iref > Imax0
    error('N=0 初始态下 Iref 不在 Id-Vg 范围内，无法定义参考 Vth0；请增大 VG 扫描范围或降低 Iref');
end

Vth0 = interp1(log10(uid0), VG0_u, log10(Iref), 'linear');
fprintf('   -> Vth0 = %.6f V\n\n', Vth0);

%% =========================
% 6) 预分配
% =========================
Vth        = nan(length(Nread),1);
DeltaVth   = nan(length(Nread),1);
ID_min_all = nan(length(Nread),1);
ID_max_all = nan(length(Nread),1);
Psum_all   = nan(length(Nread),1);
ID_all     = cell(length(Nread),1);

%% =========================
% 7) 在累计脉冲点提取 Vth
% =========================
for jj = 1:length(Nread)
    kk = Nread(jj);

    if kk < 1 || kk > length(index)
        warning('Nread(%d)=%d 超出 index 范围，跳过', jj, kk);
        continue;
    end

    row_idx = index(kk);

    if row_idx < 1 || row_idx > length(Stsum_series)
        warning('累计脉冲 %d: index=%d 超出 Stsum 范围，跳过', kk, row_idx);
        continue;
    end

    Psum = Stsum_series(row_idx);
    Psum_all(jj) = Psum;

    ID = get_ID(Psum, epife, tfe, til, miu, Na, T, W, L, VG, VFB, VD, VS);
    ID = ID(:).';
    ID_all{jj} = ID;

    ID_min_all(jj) = min(ID);
    ID_max_all(jj) = max(ID);

    fprintf('N=%d, Psum=%+.6e, ID range=[%e, %e]\n', ...
        kk, Psum, ID_min_all(jj), ID_max_all(jj));

    valid = isfinite(ID) & isfinite(VG);
    IDv = ID(valid);
    VGv = VG(valid);

    pos = IDv > 0;
    IDv = IDv(pos);
    VGv = VGv(pos);

    if isempty(IDv)
        warning('累计脉冲 %d: ID 没有正值，无法提取 Vth', kk);
        continue;
    end

    [uid, uidind] = unique(IDv);
    VG_u = VGv(uidind);

    Imin = min(uid);
    Imax = max(uid);

    if Iref < Imin || Iref > Imax
        warning(['累计脉冲 %d: Iref 不在当前 Id-Vg 范围内，' ...
                 '不做外推。Iref=%e, ID范围=[%e, %e]'], ...
                 kk, Iref, Imin, Imax);
        continue;
    end

    Vth(jj) = interp1(log10(uid), VG_u, log10(Iref), 'linear');
    fprintf('   -> Vth = %.6f V\n\n', Vth(jj));
end

%% =========================
% 8) 相对 N=0 的阈值漂移
% =========================
DeltaVth = Vth - Vth0;

fprintf('\n提取完成。\n');
fprintf('参考点使用 N=0, Vth0 = %.6f V\n', Vth0);

disp('Vth = ');
disp(Vth);

disp('DeltaVth = ');
disp(DeltaVth);

disp('Psum_all = ');
disp(Psum_all);

%% =========================
% 9) 画 DeltaVth - 累计脉冲数
% =========================
figure;
semilogx(Nread, DeltaVth, 'o-', 'LineWidth', 1.2);
xlabel('累计脉冲数');
ylabel('$\Delta V_{th}$ (V)', 'Interpreter', 'latex');
title(sprintf('累计脉冲作用下的阈值漂移 (dev\\_id=%d, ref=N=0)', dev_id));
grid on;

ax = gca;
xticks_show = 10.^(0:floor(log10(max(Nread))));
ax.XTick = xticks_show;
ax.XTickLabel = arrayfun(@(x) sprintf('10^{%d}', round(log10(x))), ...
                         xticks_show, 'UniformOutput', false);

%% =========================
% 10) 画各累计脉冲点后的 ID 范围
% =========================
figure;

ID_min_plot = ID_min_all;
ID_max_plot = ID_max_all;

ID_min_plot(ID_min_plot <= 0) = NaN;
ID_max_plot(ID_max_plot <= 0) = NaN;

plot(Nread, ID_min_plot, 'o-', 'LineWidth', 1.0); hold on;
plot(Nread, ID_max_plot, 's-', 'LineWidth', 1.0);
yline(Iref, '--', 'Iref');

set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
xlabel('累计脉冲数');
ylabel('电流 (A)');
title(sprintf('各累计脉冲点后的 ID 范围 (dev\\_id=%d)', dev_id));
legend('ID最小值', 'ID最大值', 'Iref', 'Location', 'best');
grid on;

ax = gca;
ax.XTick = xticks_show;
ax.XTickLabel = arrayfun(@(x) sprintf('10^{%d}', round(log10(x))), ...
                         xticks_show, 'UniformOutput', false);

%% =========================
% 11) 画各读取点的 Id-Vg 曲线
% =========================
figure;
hold on;
for jj = 1:length(Nread)
    ID_plot = ID_all{jj};
    if ~isempty(ID_plot)
        valid = isfinite(ID_plot) & (ID_plot > 0);
        if any(valid)
            semilogy(VG(valid), ID_plot(valid), 'LineWidth', 1.0);
        end
    end
end
xlabel('V_G (V)');
ylabel('I_D (A)');
title(sprintf('各累计脉冲点后的 I_D-V_G 曲线 (dev\\_id=%d)', dev_id));
grid on;
hold off;