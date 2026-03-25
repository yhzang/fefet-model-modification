%% ===== 提取每个累积脉冲后的 Vth 漂移 =====
% 运行前默认你已经跑过 FeFET_accumulation
% 并且工作区里已经有:
% Stsum, index, Ncycle
%
% 如果没有 Ncycle，就自动从 index 长度取

clc;

%% 1) 基本检查
if ~exist('Stsum','var')
    error('未找到变量 Stsum，请先运行 FeFET_accumulation.m');
end

if ~exist('index','var')
    error('未找到变量 index，请先确认 FeFET_accumulation.m 已输出 index');
end

if ~exist('Ncycle','var')
    Ncycle = length(index);
end

%% 2) 读出条件与器件参数
VG   = -0.5:0.02:3;   % 栅压扫描范围
VD   = 0.05;            % 漏压
VS   = 0;               % 源压
Iref = 1e-7;            % 定义阈值电流

% ---- 这些参数你可以后续再按模型背景继续调整 ----
epife = 28;             % 铁电层相对介电常数
tfe   = 0.8e-6;         % 铁电层厚度 (cm)
til   = 1e-7;           % 介电层厚度 (cm)
miu   = 50;             % 迁移率 (cm^2/V/s)
Na    = 3e17;           % 基底掺杂浓度 (cm^-3)
T     = 300;            % 温度 (K)
W     = 1;              % 晶体管宽度（单位面积模拟）
L     = 1;              % 晶体管长度（单位面积模拟）
VFB   = 0;              % 平带电压

%% 3) 预分配变量
Vth        = nan(Ncycle,1);
ID_min_all = nan(Ncycle,1);
ID_max_all = nan(Ncycle,1);

% 如需检查每个周期的 Id-Vg 曲线，可保存
ID_all = cell(Ncycle,1);

fprintf('开始提取 Vth ...\n');
fprintf('Iref = %e A\n\n', Iref);

%% 4) 循环提取每个脉冲后的 Vth
for kk = 1:Ncycle
    
    % 防止 index 越界
    if index(kk) > size(Stsum,1)
        warning('第 %d 个脉冲: index 超出 Stsum 行数，跳过', kk);
        continue;
    end

    % 这里默认第 kk 列对应第 kk 次循环
    if kk > size(Stsum,2)
        warning('第 %d 个脉冲: Stsum 列数不足，跳过', kk);
        continue;
    end

    % 第 kk 个循环后的极化值
    Psum = Stsum(index(kk), kk);

    % 计算 Id-Vg
    ID = get_ID(Psum, epife, tfe, til, miu, Na, T, W, L, VG, VFB, VD, VS);
    ID = ID(:).';   % 强制变成行向量，便于后面处理
    ID_all{kk} = ID;

    % 保存原始范围
    ID_min_all(kk) = min(ID);
    ID_max_all(kk) = max(ID);

    fprintf('kk = %d, Psum = %+e, ID range = [%e, %e]\n', ...
        kk, Psum, ID_min_all(kk), ID_max_all(kk));

    % 去除 NaN / Inf
    valid = isfinite(ID) & isfinite(VG);
    IDv = ID(valid);
    VGv = VG(valid);

    if isempty(IDv)
        warning('第 %d 个脉冲: ID 全为空或非法，无法提取 Vth', kk);
        continue;
    end

    % 只保留正电流，避免 log10 出错
    pos = IDv > 0;
    IDv = IDv(pos);
    VGv = VGv(pos);

    if isempty(IDv)
        warning('第 %d 个脉冲: ID 没有正值，无法提取 Vth', kk);
        continue;
    end

    % 去重（interp1 需要自变量单调/至少唯一）
    [uid, uidind] = unique(IDv);
    VG_u = VGv(uidind);

    % 再次检查范围
    Imin = min(uid);
    Imax = max(uid);

    if Iref < Imin || Iref > Imax
        warning(['第 %d 个脉冲: Iref 不在当前 Id-Vg 范围内，' ...
                 '不做外推。Iref=%e, ID范围=[%e, %e]'], ...
                 kk, Iref, Imin, Imax);
        Vth(kk) = NaN;
        continue;
    end

    % 用 log(I)-Vg 插值求 Vth
    Vth(kk) = interp1(log10(uid), VG_u, log10(Iref), 'linear');

    fprintf('   -> Vth(%d) = %f V\n\n', kk, Vth(kk));
end

%% 5) 选择第一个有效 Vth 作为参考点
idx0 = find(~isnan(Vth), 1, 'first');

if isempty(idx0)
    error('所有 Vth 都未成功提取，请检查 Iref 和 Id-Vg 范围。');
end

DeltaVth = Vth - Vth(idx0);

fprintf('\n提取完成。\n');
fprintf('参考点使用第 %d 个有效脉冲，Vth_ref = %f V\n', idx0, Vth(idx0));

%% 6) 画阈值漂移曲线
figure;
plot(1:Ncycle, DeltaVth, 'o-', 'LineWidth', 1.2);
xlabel('脉冲编号');
ylabel('$\Delta V_{\mathrm{th}}$ (V)', 'Interpreter', 'latex');
title('每个累积脉冲后的阈值漂移');
grid on;

%% 7) 画每个脉冲后的 ID 范围，帮助判断 Iref 是否合理
figure;
plot(1:Ncycle, ID_min_all, 'o-', 'LineWidth', 1.0); hold on;
plot(1:Ncycle, ID_max_all, 's-', 'LineWidth', 1.0);
yline(Iref, '--', 'Iref');
set(gca, 'YScale', 'log');
xlabel('脉冲编号');
ylabel('电流 (A)');
title('每个脉冲后的 ID 范围');
legend('ID最小值', 'ID最大值', 'Iref', 'Location', 'best');
grid on;

%% 8) 可选：把所有有效脉冲后的 Id-Vg 曲线画出来
figure;
hold on;
for kk = 1:Ncycle
    if ~isempty(ID_all{kk})
        ID_plot = ID_all{kk};
        valid = isfinite(ID_plot) & (ID_plot > 0);
        if any(valid)
            semilogy(VG(valid), ID_plot(valid), 'LineWidth', 1.0);
        end
    end
end
xlabel('V_G (V)');
ylabel('I_D (A)');
title('各脉冲后的 I_D-V_G 曲线');
grid on;
hold off;
%% DeltaVth
disp('Vth = ');
disp(Vth);

disp('DeltaVth = ');
disp(DeltaVth);