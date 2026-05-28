clear; clc;

% =========================================================================
% 1. 指定文件名
% =========================================================================
base_file  = '花岗岩.mat'; % 原始总表
patch_file = 'Patch_188923332778.mat'; % 包含 3 个点的补丁文件

% =========================================================================
% 2. 加载与合并逻辑
% =========================================================================
if ~exist(base_file, 'file') || ~exist(patch_file, 'file')
    error('文件读取失败，请检查文件名是否正确！');
end

S = load(base_file);   % 加载总表
P = load(patch_file);  % 加载补丁 (包含 3 列数据)

fprintf('>>> 补丁文件中包含以下厚度点：\n');
disp(P.d_scan_mm);

% 遍历补丁文件中的每一个点
for j = 1:length(P.d_scan_mm)
    target_d = P.d_scan_mm(j);
    
    % 在总表中寻找最接近该厚度的索引
    [diff_val, idx] = min(abs(S.d_scan_mm - target_d)); 
    
    if diff_val < 0.1  % 容差检查，确保确实找到了对应的点
        fprintf('>>> 正在替换总表中第 %d 列 (厚度: %.2f mm)\n', idx, S.d_scan_mm(idx));
        
        % 执行替换：将补丁的第 j 列 赋给 总表的第 idx 列
        S.RMSE_records(:, idx) = P.RMSE_records(:, j);
        S.CRLB_records(:, idx) = P.CRLB_records(:, j);
    else
        fprintf('⚠️ 警告：补丁中的 %.2f mm 在总表中找不到匹配项，已跳过。\n', target_d);
    end
end

% =========================================================================
% 3. 保存新文件
% =========================================================================
new_filename = ['MultiPoint_Fixed_', base_file];
save_data = S; 
save(new_filename, '-struct', 'save_data');

fprintf('\n>>> ✅ 合并成功！新文件已生成：%s\n', new_filename);
fprintf('>>> 现在你可以用绘图脚本打开这个新文件，那几个尖峰应该都消失了。\n');