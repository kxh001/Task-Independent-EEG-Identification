function [ slice ] = STFT( X,time_length )

channel_num = 22; 
fs = 250;
windowsize = fs*time_length;
window = rectwin(windowsize); 
nfft = windowsize;
overlap = 0;    % option:  0,windowsize/2

    for k = 1:channel_num
        str_k = num2str(k);
        [S,~,~] = spectrogram(X(:,k),window,overlap,nfft,fs); 
        S = abs(S);
        FM_slice(:,:,k) = S;
    end
    for i = 1:channel_num;
        FM_slice_new(:,i,:) = FM_slice(:,:,i);
    end
    for j = 1:size(FM_slice_new,3)
        slice(:,j) = reshape(FM_slice_new(:,:,j),[1,channel_num*size(S,1)]);
    end
end

