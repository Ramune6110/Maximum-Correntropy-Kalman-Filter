%%%% Generate Non-guassian Noise(Shot Noise)
T = 3 ;
num_shot_noise     = 15;
start_of_shotnoise = 60;
index_rand_shot    = [randi([start_of_shotnoise/T N],1,num_shot_noise-1) 21];

%% In process equation
for i=1:num_shot_noise
    temp = randi([1 5],n,1);
    MeasErrX(:,index_rand_shot(i)) =  MeasErrX(:,index_rand_shot(i)) + temp + 0.6*randn(n,1); 
end
 Q = cov(MeasErrX');

%% In Measurment equation
for i=1:num_shot_noise
    temp = randi([1 5],m,1);
    MeasErrZ(:,index_rand_shot(i)) =  MeasErrZ(:,index_rand_shot(i)) + temp + 0.6*randn(m,1); 
end
R = cov(MeasErrZ');