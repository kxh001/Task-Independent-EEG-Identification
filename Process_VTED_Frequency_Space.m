function [ dataset,datalabel ] = Process_VTED_Frequency_Space( time_length,fs )
%UNTITLED8 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
sub_num = 9;
for i = 1:sub_num
    str_i = num2str(i);
    filename1 = strcat('C:\Users\asus\Desktop\����д��\����\Class\08��������4����\subject',str_i,'_V1.mat');
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

