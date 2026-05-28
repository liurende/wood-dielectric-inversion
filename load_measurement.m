function out = load_measurement(fp)
% 读取带仪器头的 CSV/TSV（如 Keysight），从真正表头行开始解析。
% 需要的列(优先大小写敏感)： 
%   频率 : 'Freq_Hz_' | 'Freq(Hz)' | 'Freq_Hz' | 'Freq' | 'Frequency' | ...
%   幅度(dB): 'S21(DB)' | 'S21_DB_' | 'S21_DB' | 'S21 dB' | ...
%   相位(deg): 'S21(DEG)' | 'S21_DEG_' | 'S21_DEG' | 'Phase' | 'S21 DEG' | ...
%
% 返回:
%   out.fHz        [N x 1] double
%   out.S_complex  [N x 1] complex

    assert(exist(fp,'file')==2, '文件不存在: %s', fp);

    % ==== 逐行读入，兼容旧版 MATLAB ====
    L = local_read_all_lines(fp);  % cellstr
    if isempty(L), error('空文件: %s', fp); end
    L{1} = local_strip_bom(L{1});

    % ==== 寻找真正表头行 ====
    keyRegex = '(Freq|Frequency).*(Hz)|S21\(DB\)|S21_DB|S21\(DEG\)|S21_DEG|Phase';
    % 先找 BEGIN 行
    idxHeader = local_find_line_contains(L, 'BEGIN');
    if ~isempty(idxHeader)
        % 优先取 BEGIN 下一行，若不含关键字再向下搜
        idxCand = idxHeader + 1;
        if idxCand <= numel(L) && local_contains_regex(L{idxCand}, keyRegex)
            idxHeader = idxCand;
        else
            idx2 = local_find_line_regex(L, keyRegex, idxHeader+1);
            if ~isempty(idx2), idxHeader = idx2; end
        end
    else
        % 没有 BEGIN：直接找第一条含关键字的行
        idxHeader = local_find_line_regex(L, keyRegex, 1);
    end
    assert(~isempty(idxHeader), '未找到表头行（包含 Freq/S21(DB)/S21(DEG) 关键字）: %s', fp);

    headerLine = strtrim(L{idxHeader});

    % ==== 判定分隔符 ====
    delim = local_guess_delim(headerLine);

    % ==== 拆表头（保持原大小写）====
    names = local_split_header(headerLine, delim);  % cellstr, 1xM
    ncol  = numel(names);

    % ==== 收集数据行 ====
    dataLines = L(idxHeader+1:end);
    % 去掉空行和以 '!' 开头的说明行
    mask = true(size(dataLines));
    for i=1:numel(dataLines)
        s = strtrim(dataLines{i});
        if isempty(s) || (~isempty(s) && s(1)=='!')
            mask(i) = false;
        end
    end
    dataLines = dataLines(mask);

    % 对齐到列数
    cols = cell(1,ncol);
    for j=1:ncol, cols{j} = strings(0,1); end
    for i=1:numel(dataLines)
        row = strsplit(dataLines{i}, delim);
        % 裁剪/补齐
        if numel(row) >= ncol
            row = row(1:ncol);
        else
            row(end+1:ncol) = {''};
        end
        for j=1:ncol
            cols{j}(end+1,1) = local_strip_quotes(strtrim(string(row{j})));
        end
    end

    % ==== 锁定三列（先大小写敏感，失败再大小写不敏感）====
    wantF = {'Freq_Hz_','Freq(Hz)','Freq_Hz','Freq','Frequency','Frequency_Hz_','Frequency(Hz)'};
    wantA = {'S21(DB)','S21_DB_','S21_DB','S21 dB','S21dB'};
    wantP = {'S21(DEG)','S21_DEG_','S21_DEG','Phase','S21 DEG','S21deg'};

    iF = pick_by_alias_cs(names, wantF, false);
    if isempty(iF), iF = pick_by_alias_ci(names, wantF, '频率'); end

    iA = pick_by_alias_cs(names, wantA, false);
    if isempty(iA), iA = pick_by_alias_ci(names, wantA, '幅度(S21 dB)'); end

    iP = pick_by_alias_cs(names, wantP, false);
    if isempty(iP), iP = pick_by_alias_ci(names, wantP, '相位(S21 DEG)'); end

    % ==== 转数值 ====
    fHz  = str2double_clean(cols{iF});
    ADB  = str2double_clean(cols{iA});
    Pdeg = str2double_clean(cols{iP});

    % ==== 组 S21 & 清洗 ====
    mag = 10.^(ADB(:)/20);
    ph  = deg2rad(Pdeg(:));
    S   = mag .* exp(1j*ph);

    [fHz, S] = local_clean_series(fHz(:), S(:));

    out = struct('fHz', fHz, 'S_complex', S);
end

%% ====== 本文件内部工具 ======
function L = local_read_all_lines(fp)
    fid = fopen(fp,'r','n','UTF-8');
    if fid<0, error('无法打开文件: %s', fp); end
    C = {};
    i = 0;
    while true
        t = fgetl(fid);
        if ~ischar(t), break; end
        i=i+1; C{i,1} = t; %#ok<AGROW>
    end
    fclose(fid);
    if isempty(C), L = {}; else, L = C; end
end

function s = local_strip_bom(s)
    if isempty(s), return; end
    % UTF-8 BOM (EF BB BF)
    u8 = uint8(s);
    if numel(u8)>=3 && isequal(u8(1:3), uint8([239 187 191]))
        s = char(u8(4:end));
    end
    % Unicode BOM char
    if ~isempty(s) && s(1)==char(65279)
        s = s(2:end);
    end
end

function idx = local_find_line_contains(L, token)
    idx = [];
    for i=1:numel(L)
        if contains(L{i}, token, 'IgnoreCase', true)
            idx = i; return;
        end
    end
end

function tf = local_contains_regex(s, re)
    tf = ~isempty(regexp(s, re, 'once', 'ignorecase'));
end

function idx = local_find_line_regex(L, re, startIdx)
    if nargin<3, startIdx=1; end
    idx = [];
    for i=startIdx:numel(L)
        if local_contains_regex(L{i}, re)
            idx = i; return;
        end
    end
end

function d = local_guess_delim(headerLine)
    % 返回单字符分隔符
    if contains(headerLine, sprintf('\t')), d = sprintf('\t'); return; end
    c_comma = sum(headerLine==',');
    c_tab   = sum(headerLine==sprintf('\t'));
    c_sc    = sum(headerLine==';');
    [~,k]   = max([c_comma,c_tab,c_sc]);
    d = ',';
    if k==2, d=sprintf('\t'); end
    if k==3, d=';'; end
end

function names = local_split_header(headerLine, delim)
    P = strsplit(headerLine, delim);
    names = cellfun(@(s) char(local_strip_quotes(string(strtrim(s)))), P, 'UniformOutput', false);
end

function s = local_strip_quotes(s)
    if strlength(s)>=2
        if s.extractBetween(1,1) == '"' && s.extractBetween(strlength(s),strlength(s)) == '"'
            s = extractBetween(s, 2, strlength(s)-1);
        end
    end
end

function v = str2double_clean(col)
    if iscell(col), col = string(col); end
    col = strrep(col, ",", "");
    col = strrep(col, "，", "");
    col = replace(col, '"', "");
    v = str2double(col);
end

function idx = pick_by_alias_cs(names, cand, must)
    % 大小写敏感优先匹配
    idx = [];
    for i=1:numel(cand)
        j = find(strcmp(names, cand{i}), 1, 'first');
        if ~isempty(j), idx = j; return; end
    end
    if nargin>=3 && must
        error('未找到所需列：%s。当前列名：%s', ...
              strjoin(cand, ', '), strjoin(names, ', '));
    end
end

function idx = pick_by_alias_ci(names, cand, tag)
    % 大小写不敏感兜底
    idx = [];
    for i=1:numel(cand)
        j = find(strcmpi(names, cand{i}), 1, 'first');
        if ~isempty(j)
            fprintf('[提示] %s列使用不区分大小写命中：%s\n', tag, names{j});
            idx = j; return;
        end
    end
    error('未找到所需列（%s）：%s。当前列名：%s', tag, strjoin(cand, ', '), strjoin(names, ', '));
end

function [x, y] = local_clean_series(x, y)
    m = isfinite(x) & isfinite(real(y)) & isfinite(imag(y)) & x>0;
    x = x(m); y = y(m);
    if isempty(x), x = 1; y = 0; return; end
    [x, idx] = sort(x); y = y(idx);
    [ux, ~, ic] = unique(x);
    if numel(ux) < numel(x)
        yr = accumarray(ic, [real(y), imag(y)], [], @(v) mean(v,1,'omitnan'));
        y  = complex(yr(:,1), yr(:,2));
        x  = ux;
    end
    if numel(x)==1
        x = x + [-1;1]*max(1,x*1e-9);
        y = y([1 1]);
    end
end
