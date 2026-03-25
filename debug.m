%% =========================================
% debug_readout_overlay_fix.m
% 目的：检查 Pread 变化是否真的传到了 ID-VG
%% =========================================

% ---------- 1) 先整理 Stsum ----------
if ~exist('Stsum', 'var')
    error('工作区没有 Stsum，请先跑主仿真。');
end
if ~exist('index', 'var')
    error('工作区没有 index。');
end
if ~exist('index_pre', 'var')
    error('工作区没有 index_pre。');
end

if isvector(Stsum)
    st_trace = Stsum(:);   % 强制列向量
    trace_id = 1;
else
    trace_id = 1;          % 先固定看第1列
    st_trace = Stsum(:, trace_id);
end

fprintf('size(Stsum) = [%d %d]\n', size(Stsum,1), size(Stsum,2));
fprintf('当前使用 trace_id = %d\n', trace_id);

% ---------- 2) 读出参数 ----------
VG = -0.5:0.02:1.7;
VD = 0.05;
VS = 0;
VFB_read = -0.5;
Iref = 1e-7;

% ---------- 3) 器件参数 ----------
tfe = 0.8e-6;
til = 1e-7;
miu = 50;
Na  = 3e17;
T   = 300;
W   = 1;
L   = 1;

% ---------- 4) 选几个累计脉冲点 ----------
selN = [1 10 100 1000];
selN = selN(selN <= numel(index));
selN = selN(:);

% 初始参考点
P0 = st_trace(index_pre(1));

% 读点极化
Psel = [P0; st_trace(index(selN))];

label_txt = cell(numel(Psel),1);
label_txt{1} = sprintf('N=0, P=%.6f', P0);
for k = 1:numel(selN)
    label_txt{k+1} = sprintf('N=%d, P=%.6f', selN(k), Psel(k+1));
end

% ---------- 5) 计算 ID-VG ----------
ID = zeros(numel(Psel), numel(VG));
Vth_repo  = nan(numel(Psel),1);
Vth_cross = nan(numel(Psel),1);

for k = 1:numel(Psel)
    ID(k,:) = get_ID(Psel(k), epife, tfe, til, miu, Na, T, W, L, ...
                     VG, VFB_read, VD, VS);

    % ===== 原仓库风格提取 =====
    [uid, uidind] = unique(ID(k,:));
    uid = real(uid);
    good = isfinite(uid) & (uid > 0);

    if nnz(good) >= 2
        uid2 = uid(good);
        vg2  = VG(uidind(good));
        Vth_repo(k) = interp1(log10(uid2), vg2, log10(Iref), 'linear', NaN);
    end

    % ===== crossing 风格提取 =====
    cross_idx = find( ...
        isfinite(ID(k,1:end-1)) & isfinite(ID(k,2:end)) & ...
        (ID(k,1:end-1) > 0) & (ID(k,2:end) > 0) & ...
        ((ID(k,1:end-1)-Iref) .* (ID(k,2:end)-Iref) <= 0), ...
        1, 'first');

    if ~isempty(cross_idx)
        x = log10([ID(k,cross_idx), ID(k,cross_idx+1)]);
        y = [VG(cross_idx), VG(cross_idx+1)];
        Vth_cross(k) = interp1(x, y, log10(Iref), 'linear');
    end
end

% ---------- 6) 画图 ----------
figure;
semilogy(VG, ID.', 'LineWidth', 1.5);
xlabel('V_G (V)');
ylabel('I_D (A)');
title('ID-VG overlay for selected ferroelectric states');
grid on;
legend(label_txt, 'Location', 'best');

% 画参考电流水平线
hold on;
yline(Iref, '--');
hold off;

% ---------- 7) 文本输出 ----------
fprintf('\n=== Readout debug ===\n');
for k = 1:numel(Psel)
    fprintf('%s : Vth_repo = %.6f sV,   Vth_cross = %.6f V\n', ...
        label_txt{k}, Vth_repo(k), Vth_cross(k));
end

% ---------- 8) 粗略等效电压尺度 ----------
Cfe = 8.85e-14 * epife / tfe;
dP = (max(Psel) - min(Psel)) / 1e6;   % uC/cm^2 -> C/cm^2
dVeq = dP / Cfe;

fprintf('\nCfe  = %.4e F/cm^2\n', Cfe);
fprintf('dP   = %.4e C/cm^2\n', dP);
fprintf('dVeq = %.4f V (rough estimate)\n', dVeq);