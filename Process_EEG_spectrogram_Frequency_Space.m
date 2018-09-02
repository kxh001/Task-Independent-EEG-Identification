function [ dataset,datalabel ] = Process_EEG_spectrogram_Frequency_Space( time_length,fs,sub_num,channel_num )

for i = 1:sub_num
    str_i = num2str(i);
    filename1 = 'XXX';
    load(filename1); %named X 
    [ slice ] = STFT(X,time_length,fs,channel_num);
    dataset{i} = slice;
    datalabel{i} = repmat([i],size(slice,2),1);
end

end

