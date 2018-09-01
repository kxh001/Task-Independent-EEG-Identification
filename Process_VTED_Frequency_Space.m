function [ dataset,datalabel ] = Process_VTED_Frequency_Space( time_length,fs )
%UNTITLED8 此处显示有关此函数的摘要
%   此处显示详细说明
sub_num = 9;
for i = 1:sub_num
    str_i = num2str(i);
    filename1 = strcat('C:\Users\asus\Desktop\论文写作\数据\Class\08竞赛数据4分类\subject',str_i,'_V1.mat');
    load(filename1);
%     X = data.X(1:340000,:);
    X = trials.data(:,1:22);
    X(isnan(X)==1) = 0;
    [ slice ] = STFT(X,time_length);
    dataset{i} = slice;
    datalabel{i} = repmat([i],size(slice,2),1);
%     datalabel{i} = repmat([i],40,1);
end

end

