function [time, volt, index, index_pre, index_end, info] = wfdef_acc_fix(amp, pw, step, delay, transit, cycle)
% =========================================================
% wfdef_acc_fix.m
%
% 作用：
%   为累积脉冲/写扰分析生成波形，并返回更合理的索引定义
%
% 输入：
%   amp      - 脉冲幅值
%   pw       - 脉宽
%   step     - 时间步
%   delay    - 脉冲间隔（这里按“相邻 pulse 之间只有一次 delay”解释）
%   transit  - 上升/下降时间
%   cycle    - 脉冲个数
%
% 输出：
%   time       - 全局时间轴（行向量）
%   volt       - 全局电压波形（行向量）
%   index      - 推荐读点：每个 pulse 结束并经过 trailing delay 后
%   index_pre  - 每个 pulse 开始前最后一个零电压点
%   index_end  - 每个 pulse 最后一个非零电压点
%   info       - 调试信息结构体
%
% 依赖：
%   wfdef_single.m
%
% 说明：
%   原仓库 wfdef_acc.m 的 index 实际更接近“pulse 开始前”，
%   不适合直接拿来做 disturb 后的 Vth 读出。
% =========================================================

    if nargin ~= 6
        error('wfdef_acc_fix:InputCount', ...
            '需要 6 个输入：amp, pw, step, delay, transit, cycle');
    end

    if cycle < 1 || abs(cycle - round(cycle)) > eps(max(cycle,1))
        error('wfdef_acc_fix:BadCycle', 'cycle 必须是正整数。');
    end
    cycle = round(cycle);

    % -------------------------
    % 初始化 pulse（沿用原仓库思路）
    % -------------------------
    amp_init = -4;
    pw_init  = 1e-6;

    [time, volt] = wfdef_single(amp_init, pw_init, pw_init/30, delay, transit);
    time = reshape(time, 1, []);
    volt = reshape(volt, 1, []);
    tcurr = time(end);

    % -------------------------
    % 预生成一个“标准 disturb pulse 段”
    % 仍用 wfdef_single 生成，但后面会剥掉前导 delay
    % -------------------------
    [pulse_t, pulse_v] = wfdef_single(amp, pw, step, delay, transit);
    pulse_t = reshape(pulse_t, 1, []);
    pulse_v = reshape(pulse_v, 1, []);

    tol = max(1e-15, 1e-12 * max(1, abs(amp)));
    nz = find(abs(pulse_v) > tol);

    if isempty(nz)
        error('wfdef_acc_fix:NoNonZeroPulse', ...
            '生成的 pulse 段里没有非零电压，检查 amp / pw / step / transit。');
    end

    % pulse 开始前最后一个零点
    local_pre = max(nz(1) - 1, 1);

    % pulse 最后一个非零点（通常是下降沿末端）
    local_end = nz(end);

    % 整个 pulse 段结束点（含 trailing delay）
    local_read = numel(pulse_t);

    % -------------------------
    % 输出索引
    % -------------------------
    index     = zeros(cycle, 1);  % 推荐读点
    index_pre = zeros(cycle, 1);  % pulse 前
    index_end = zeros(cycle, 1);  % pulse 末端（最后非零）

    % ------------------------------------------------------
    % 对每个 disturb pulse：
    %   只保留“从 local_pre 开始”的部分
    %   并把 local_pre 对齐到当前全局末尾 tcurr
    %
    % 这样相邻 pulse 之间只保留一次 delay，
    % 不会出现原始串接方式中的“双 delay”
    % ------------------------------------------------------
    for ii = 1:cycle
        L = numel(time);

        % 当前全局末尾点，就是本 pulse 的“pre”读点
        index_pre(ii) = L;

        % 映射到全局：
        % local_pre   -> 现有最后一点 L
        % local_end   -> L + (local_end - local_pre)
        % local_read  -> L + (local_read - local_pre)
        index_end(ii) = L + (local_end  - local_pre);
        index(ii)     = L + (local_read - local_pre);

        % 将 local_pre 对齐到当前 tcurr
        tshift = tcurr - pulse_t(local_pre);

        % 追加 local_pre 之后的点；local_pre 本身不重复加
        time = [time, pulse_t(local_pre+1:end) + tshift];
        volt = [volt, pulse_v(local_pre+1:end)];

        tcurr = time(end);
    end

    % -------------------------
    % 调试信息
    % -------------------------
    info = struct();
    info.amp_init   = amp_init;
    info.pw_init    = pw_init;
    info.local_pre  = local_pre;
    info.local_end  = local_end;
    info.local_read = local_read;
    info.index_read = index;
    info.delay_note = '本 fix 版把相邻 pulse 的静默区解释成一次 delay';
end