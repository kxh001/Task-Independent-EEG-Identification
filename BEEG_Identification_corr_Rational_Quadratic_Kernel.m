%Applying GoDec+ for BEEG classification
clear
tic

time_length = 20;

fs = 250;
[ dataset,datalabel ] = Process_VTED_Frequency_Space( time_length,fs );
toc
r=1:1; %������ĸ�����Ҳ����rank��1~6����һ��
rep_times = 1;

q = 0;

% sigma = 0.0001;
C_i = 0.0001;

train_num = ceil(0.5*size(datalabel{1},1));
epsilon = 1e-7;

accuracy = zeros(length(r),rep_times);

for r_i = 1:length(r)
    Re = zeros(1,rep_times);
    parfor k = 1:rep_times
%-------------�ֲ��Ժ�ѵ����-------------%
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
        
        
%----------------��һ������ֹ�����------------------%
        for i = 1:size(trainSet,2)
            trainSet(:,i) = trainSet(:,i)/norm(trainSet(:,i));
        end
        for i = 1:size(testSet,2)
            testSet(:,i) = testSet(:,i)/norm(testSet(:,i));
        end

%--------------------��ѵ��������ϳ��ع�����h-------------------%
        L = cell(1,length(dataset));
        id = 0;
        for i = 1:length(L)
            trainSet0{i} = trainSet(:,id+1:id+size(trainSet0{i},2)); %����ʹ�ù�һ��֮������ݣ�֮ǰ��trainSet0���������û����
            id = id +size(trainSet0{i},2);
%             [L{i},RMSE,~,Q]=lowrank_corr(trainSet0{i},r(r_i),sigma,epsilon,q);
            [L{i},RMSE,~,Q]=lowrank_corr_RQK(trainSet0{i},r(r_i),C_i,epsilon,q);
            [L{i},~]=qr(L{i},0);
            L{i} = L{i}(:,1:r(r_i));
        end
           
                    L = cell2mat(L);
                    H0 = zeros(size(L,2),size(testSet,2));
                    iter = 1;
                    e = zeros(size(testSet));
                    while true
                        testSet_tmp = testSet-e;
                        H = pinv(L)*testSet_tmp;
                        T = testSet - L*H;
                        T_sq = T.*T;
%                         e = T - T.*exp(-T_sq/sigma);
                        e = T - T.*(1-T_sq./(T_sq+C_i));
                        tmp = H-H0;
                        if norm(tmp(:))<1e-7 || iter >100
                            break;
                        end
                        H0 = H;
                        iter = iter + 1;
                    end

                    

%----------------�������������ÿ����Ľ����أ���󽻲���Ϊ������-----------------%
                    corr = zeros(length(dataset),size(testSet,2));
                    r_ = r(r_i);
                    for i = 1:length(dataset)
                        tmp = testSet - L(:,((i-1)*r_+1):i*r_)*H(((i-1)*r_+1):i*r_,:);
%                         corr(i,:) = sum(exp(-tmp.*tmp/sigma),1);
                        tmp_sq = tmp.*tmp;
                        corr(i,:) = sum((1-tmp_sq./(tmp_sq+C_i)),1);
                    end
                    [~,result_label]=max(corr,[],1);
                    result_label = result_label';

%----------------------����׼ȷ��--------------------%
        result = sum(result_label==testLabel)/length(testLabel);
        Re(k) = result;
    end
    toc
    accuracy(r_i,:) = Re;
    time(r_i,:) = toc;
end

