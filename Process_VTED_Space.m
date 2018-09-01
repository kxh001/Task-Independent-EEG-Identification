function [ dataset,datalabel ] = Process_VTED_Space( time_length,fs )
%UNTITLED8 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

for i = 1:15
    str_i = num2str(i);
    filename1 = strcat('G:\�����޹����ݼ�\signalF',str_i);
    load(filename1);
    for j = 1:(160000/time_length/fs)
         Temp(:,j) = mean(signalF((j-1)*time_length*fs+1:time_length*fs*j,1:62),[62*time_length*fs,1]);
    end
    dataset{i} = Temp;
    datalabel{i} = repmat([i],(160000/time_length/fs),1);
end


end

