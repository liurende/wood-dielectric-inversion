clear; clc;

% 1. 指定文件名（请根据你的实际文件名修改）
old_file = 'Verification_Result_20260514_181649.mat'; % 填入那个有尖峰的文件名
patch_file = 'Patch_3222366741114556500';

% 2. 加载数据
S = load(old_file);   
P = load(patch_file); 

% 3. 寻找 23.33mm 在原数据矩阵中的【列索引】
target_d = 14.44;
[~, idx] = min(abs(S.d_scan_mm - target_d)); 

fprintf('>>> 匹配到原点：%.2f mm (位于第 %d 列)\n', S.d_scan_mm(idx), idx);

% 4. 【核心修复】：按“列”进行数据替换，而不是按“行”！
% 将补丁的第 1 列（4个参数），塞进原矩阵的第 idx 列
S.RMSE_records(:, idx) = P.RMSE_records(:, 1);
S.CRLB_records(:, idx) = P.CRLB_records(:, 1);

% 5. 保存回原来的格式
new_filename = ['Perfect2_Merged_', old_file];
save_data = S; 
save(new_filename, '-struct', 'save_data');

fprintf('>>> 成功！新文件已生成：%s\n', new_filename);
fprintf('>>> 那个平直的直线 BUG 已经修复，尖峰也被削平了，快去画图看看吧！\n');