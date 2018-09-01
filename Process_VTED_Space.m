function [ dataset,datalabel ] = Process_VTED_Space( time_length,fs )
%UNTITLED8 此处显示有关此函数的摘要
%   此处显示详细说明

for i = 1:15
    str_i = num2str(i);
    filename1 = strcat('G:\任务无关数据集\signalF',str_i);
    load(filename1);
    for j = 1:(160000/time_length/fs)
         Temp(:,j) = mean(signalF((j-1)*time_length*fs+1:time_length*fs*j,1:62),[62*time_length*fs,1]);
    end
    dataset{i} = Temp;
    datalabel{i} = repmat([i],(160000/time_length/fs),1);
end


end

