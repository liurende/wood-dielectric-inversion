%% ========================================================================
%%  [最终修复版] S11 数据加载函数 (兼容单文件加载 + 批量搜索)
%% ========================================================================
function out = load_meas_data_s11(varargin)
    % ---------------------------------------------------------------------
    % 模式 1: 单文件直接加载 (兼容主程序第 96 行 airRef 的调用)
    % ---------------------------------------------------------------------
    if nargin == 1
        filePath = varargin{1};
        % 直接调用内部读取函数，返回原始数据结构
        out = local_load_measurement_s11_v2(filePath);
        return;
    end

    % ---------------------------------------------------------------------
    % 模式 2: 批量数据搜索与加载 (标准调用)
    % ---------------------------------------------------------------------
    if nargin < 5
        error('load_meas_data_s11: 参数不足。批量模式需要 5 个参数。');
    end

    % 解包参数
    baseName    = varargin{1};
    measDirBase = varargin{2};
    RefData     = varargin{3};
    fGHz        = varargin{4};
    subFolders  = varargin{5};

    measData = struct(); 
    
    % 1. 英文名 -> 中文文件名 映射表
    nameMap = containers.Map();
    nameMap('granite') = '花岗岩';
    nameMap('marble')  = '天然大理石';
    nameMap('pe')      = '板';
    nameMap('board')   = '板';
    nameMap('al')      = '板'; 
    
    lowerName = lower(baseName);
    targetChinese = '';
    
    keys = nameMap.keys;
    for i = 1:length(keys)
        if contains(lowerName, keys{i})
            targetChinese = nameMap(keys{i});
            break;
        end
    end
    
    if isempty(targetChinese)
        searchPatterns = {baseName};
    else
        searchPatterns = {targetChinese, '汉白玉', baseName}; 
    end

    fprintf('\n[Debug] S11 Loading for "%s" (Target: %s):\n', baseName, strjoin(searchPatterns, '/'));

    for k = 1:2
        angStr = subFolders{k}; % '30' or '45'
        measData.(['deg' angStr]).valid = false;
        
        folderPath = fullfile(measDirBase, angStr);
        
        if ~exist(folderPath, 'dir')
            fprintf('   [Error] Folder not found: %s\n', folderPath);
            continue;
        end
        
        % 搜索文件
        foundFile = [];
        foundName = '';
        
        for p = 1:length(searchPatterns)
            pat = searchPatterns{p};
            d = dir(fullfile(folderPath, [pat '*.csv']));
            d = d(~startsWith({d.name}, '.')); % 排除隐藏文件
            
            if ~isempty(d)
                % 优先避开 "板" (除非本身就在找板)
                idx_best = 1;
                if ~contains(pat, '板')
                    for fix = 1:length(d)
                        if ~contains(d(fix).name, '板')
                            idx_best = fix; 
                            break; 
                        end
                    end
                end
                foundName = d(idx_best).name;
                foundFile = fullfile(d(idx_best).folder, foundName);
                break; 
            end
        end
        
        if ~isempty(foundFile)
            fprintf('      -> Angle %s: Found "%s"\n', angStr, foundName);
            try
                % 调用读取函数
                raw = local_load_measurement_s11_v2(foundFile);
                
                % 插值对齐
                measData.(['deg' angStr]).mag = interp1(raw.fHz/1e9, 20*log10(abs(raw.S_complex)), fGHz, 'linear', 'extrap');
                measData.(['deg' angStr]).f = fGHz;
                measData.(['deg' angStr]).valid = true;
            catch ME
                fprintf('      -> [Error] Read failed: %s\n', ME.message);
            end
        else
            fprintf('      -> [Fail] No match in %s\n', folderPath);
        end
    end
    
    out = measData; % 返回处理好的批量数据
end

%% ========================================================================
%%  本地工具: CSV 读取 (支持 S11/LogMag/Refl)
%% ========================================================================
function out = local_load_measurement_s11_v2(fp)
    assert(exist(fp,'file')==2, '文件不存在: %s', fp);
    L = local_read_all_lines(fp);
    if isempty(L), error('空文件: %s', fp); end
    L{1} = local_strip_bom(L{1});
    
    % 正则匹配表头 (兼容 S11, LogMag, S21)
    keyRegex = '(Freq|Frequency).*(Hz)|S11\(DB\)|S11_DB|S11\(DEG\)|S11_DEG|LogMag|LinMag|S21\(DB\)|S21_DB|Phase';
    
    idxHeader = local_find_line_contains(L, 'BEGIN');
    if ~isempty(idxHeader)
        idxCand = idxHeader + 1;
        if idxCand <= numel(L) && local_contains_regex(L{idxCand}, keyRegex)
            idxHeader = idxCand;
        else
            idx2 = local_find_line_regex(L, keyRegex, idxHeader+1);
            if ~isempty(idx2), idxHeader = idx2; end
        end
    else
        idxHeader = local_find_line_regex(L, keyRegex, 1);
    end
    assert(~isempty(idxHeader), '未找到表头行（缺少 Freq/S11/LogMag 关键字）: %s', fp);
    
    headerLine = strtrim(L{idxHeader});
    delim = local_guess_delim(headerLine);
    names = local_split_header(headerLine, delim);
    ncol  = numel(names);
    
    dataLines = L(idxHeader+1:end);
    mask = true(size(dataLines));
    for i=1:numel(dataLines)
        s = strtrim(dataLines{i});
        if isempty(s) || (~isempty(s) && s(1)=='!'), mask(i) = false; end
    end
    dataLines = dataLines(mask);
    
    cols = cell(1,ncol);
    for j=1:ncol, cols{j} = strings(0,1); end
    for i=1:numel(dataLines)
        row = strsplit(dataLines{i}, delim);
        if numel(row) >= ncol, row = row(1:ncol); else, row(end+1:ncol) = {''}; end
        for j=1:ncol, cols{j}(end+1,1) = local_strip_quotes(strtrim(string(row{j}))); end
    end
    
    % 列名匹配
    wantF = {'Freq_Hz_','Freq(Hz)','Freq_Hz','Freq','Frequency','Frequency_Hz_','Frequency(Hz)'};
    wantA = {'S11(DB)','S11_DB','S11 dB','S11dB','LogMag','S11 Log Mag','S21(DB)','S21_DB'}; 
    wantP = {'S11(DEG)','S11_DEG','S11 Phase','Phase','S11 DEG','S21(DEG)','S21_DEG'};
    
    iF = pick_by_alias_cs(names, wantF, false);
    if isempty(iF), iF = pick_by_alias_ci(names, wantF, '频率'); end
    
    iA = pick_by_alias_cs(names, wantA, false);
    if isempty(iA), iA = pick_by_alias_ci(names, wantA, '幅度(dB)'); end
    
    iP = pick_by_alias_cs(names, wantP, false);
    if isempty(iP), iP = pick_by_alias_ci(names, wantP, '相位'); end 
    
    fHz  = str2double_clean(cols{iF});
    ADB  = str2double_clean(cols{iA});
    
    if ~isempty(iP)
        Pdeg = str2double_clean(cols{iP});
    else
        Pdeg = zeros(size(ADB)); 
    end
    
    mag = 10.^(ADB(:)/20);
    ph  = deg2rad(Pdeg(:));
    S   = mag .* exp(1j*ph);
    
    [fHz, S] = local_clean_series(fHz(:), S(:));
    out = struct('fHz', fHz, 'S_complex', S);
end

%% ====== 辅助函数 (压缩版) ======
function L=local_read_all_lines(fp), fid=fopen(fp,'r','n','UTF-8'); C={}; i=0; while 1, t=fgetl(fid); if ~ischar(t), break; end; i=i+1; C{i,1}=t; end; fclose(fid); L=C; end
function s=local_strip_bom(s), if numel(s)>=3 && isequal(uint8(s(1:3)),[239 187 191]), s=s(4:end); end; if ~isempty(s)&&s(1)==65279, s=s(2:end); end; end
function idx=local_find_line_contains(L,t), idx=[]; for i=1:numel(L), if contains(L{i},t,'IgnoreCase',true), idx=i; return; end; end; end
function tf=local_contains_regex(s,r), tf=~isempty(regexp(s,r,'once','ignorecase')); end
function idx=local_find_line_regex(L,r,sI), idx=[]; for i=sI:numel(L), if local_contains_regex(L{i},r), idx=i; return; end; end; end
function d=local_guess_delim(h), if contains(h,sprintf('\t')),d=sprintf('\t');return;end; [~,k]=max([sum(h==','),sum(h==char(9)),sum(h==';')]); d=','; if k==2,d=char(9);elseif k==3,d=';';end; end
function n=local_split_header(h,d), P=strsplit(h,d); n=cellfun(@(s) char(local_strip_quotes(string(strtrim(s)))), P, 'UniformOutput', false); end
function s=local_strip_quotes(s), if strlength(s)>=2 && s{1}(1)=='"' && s{1}(end)=='"', s=extractBetween(s,2,strlength(s)-1); end; end
function v=str2double_clean(c), if iscell(c),c=string(c);end; c=strrep(strrep(replace(c,'"',""),",",""),"，",""); v=str2double(c); end
function idx=pick_by_alias_cs(n,c,m), idx=[]; for i=1:numel(c), j=find(strcmp(n,c{i}),1); if ~isempty(j), idx=j; return; end; end; if m, error('Col not found'); end; end
function idx=pick_by_alias_ci(n,c,t), idx=[]; for i=1:numel(c), j=find(strcmpi(n,c{i}),1); if ~isempty(j), idx=j; return; end; end; end
function [x,y]=local_clean_series(x,y), m=isfinite(x)&isfinite(y)&x>0; x=x(m); y=y(m); [x,I]=sort(x); y=y(I); [ux,~,ic]=unique(x); if numel(ux)<numel(x), yr=accumarray(ic, [real(y),imag(y)],[],@mean); y=complex(yr(:,1),yr(:,2)); x=ux; end; end