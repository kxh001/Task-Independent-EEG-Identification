%Applying improved GoDec+ for BEEG classification on BCI Graz Dataset A
clear
tic

time_length = 20;
fs = 250;
sub_num = 9;
channel_num = 22;
[ dataset,datalabel ] = Process_EEG_spectrogram_Frequency_Space( time_length,fs,sub_num,channel_num );
r=1:10; 
rep_times = 10;

q = 0;
C_i = 0.0001;

train_num = ceil(0.5*size(datalabel{1},1));
epsilon = 1e-7;

accuracy = zeros(length(r),rep_times);

for r_i = 1:length(r)
    Re = zeros(1,rep_times);
    parfor k = 1:rep_times
%-------------分测试和训练集-------------%
        trainSet = cell(1,length(dataset));
        testSet = cell(1,length(dataset));
        trainLabel = cell(length(dataset),1);
        testLabel = cell(length(dataset),1);
        for i = 1:length(dataset)
            idx = randperm(size(dataset{i},2));
            trainSet{i} = dataset{i}(:,idx(1:train_num));
            testSet{i} = dataset{i}(:,idx(train_num+1:end));
            trainLabel{i} = datalabel{i}(idx(1:train_num));
            testLabel{i} = datalabel{i}(idx(train_num+1:end));
        end
        
        trainSet0 = trainSet;
        trainSet = cell2mat(trainSet);
        testSet = cell2mat(testSet);
        trainLabel = cell2mat(trainLabel);
        testLabel = cell2mat(testLabel);
        
        
%----------------归一化，防止过拟合------------------%
        for i = 1:size(trainSet,2)
            trainSet(:,i) = trainSet(:,i)/norm(trainSet(:,i));
        end
        for i = 1:size(testSet,2)
            testSet(:,i) = testSet(:,i)/norm(testSet(:,i));
        end

%--------------------用训练样本拟合出重构因子h-------------------%
        B = cell(1,length(dataset));
        R = cell(1,length(dataset));
        id = 0;
        for i = 1:length(B)
            trainSet0{i} = trainSet(:,id+1:id+size(trainSet0{i},2)); 
            id = id +size(trainSet0{i},2);
            [B{i},R{i}]=lowrank_corr_RQK(trainSet0{i},r(r_i),C_i,epsilon,q);
            SEEG=trainSet0{i};
            BEEG=B{i};
            REEG=R{i};
            iii=i;
            if iii<6
                 subplot(3,5,iii);
                imagesc(SEEG);axis off;title(strcat('SEEG-',num2str(iii)));
                subplot(3,5,5+iii);
                imagesc(BEEG);axis off;title(strcat('BEEG-',num2str(iii)));
                subplot(3,5,10+iii);
                imagesc(REEG);axis off;title(strcat('REEG-',num2str(iii)));
            end
            [B{i},~]=qr(B{i},0);
            B{i} = B{i}(:,1:r(r_i));
        end
           
                    B = cell2mat(B);
                    H0 = zeros(size(B,2),size(testSet,2));
                    iter = 1;
                    e = zeros(size(testSet));
                    while true
                        testSet_tmp = testSet-e;
                        H = pinv(B)*testSet_tmp;
                        T = testSet - B*H;
                        T_sq = T.*T;
                        e = T - T.*(1-T_sq./(T_sq+C_i));
                        tmp = H-H0;
                        if norm(tmp(:))<1e-7 || iter >100
                            break;
                        end
                        H0 = H;
                        iter = iter + 1;
                    end

                    

%----------------计算测试样本和每个类的交叉熵，最大交叉熵为归属类-----------------%
                    corr = zeros(length(dataset),size(testSet,2));
                    r_ = r(r_i);
                    for i = 1:length(dataset)
                        tmp = testSet - B(:,((i-1)*r_+1):i*r_)*H(((i-1)*r_+1):i*r_,:);
                        tmp_sq = tmp.*tmp;
                        corr(i,:) = sum((1-tmp_sq./(tmp_sq+C_i)),1);
                    end
                    [~,result_label]=max(corr,[],1);
                    result_label = result_label';

%----------------------计算准确率--------------------%
        result = sum(result_label==testLabel)/length(testLabel);
        Re(k) = result;
    end
    toc
    accuracy(r_i,:) = Re;
    time(r_i,:) = toc;
end

