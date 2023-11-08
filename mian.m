%% 基于样本量变化的参数学习效果对比
clc; clear all; close all;
set(0,'RecursionLimit',10000);    % 设置迭代次数

%% 设置实验参数
N_samples = 1000; %      样本量
% % BN 的名字库
bname = 'weather';
% figure(1)
[bnet,nodes,arcs,parameters] = GetBNet(bname);
% draw_graph(bnet.dag);
max_parents = 4;
n = nodes;
ns = bnet.node_sizes;
T = 20;
time = zeros(1,T);
score = zeros(1,T);
%************抽样数据***********************
nsamples=N_samples;
N=n;
order=randperm(n);
samples = cell(N,nsamples);
for m = 1:nsamples                                                    % 抽样数据
    samples(:,m) = sample_bnet(bnet);
end
data=zeros(N,nsamples);
for i=1:N
    for m = 1:nsamples                                                    % 抽样数据
        data(i,m) = samples{i,m};
    end
end
%%***********学习结构***********
for iter=1:T
    odr0 = randperm(n);
    tic
    [dag_gbn_ncos,scor,cache] = yy_HC_order(data,max_parents,ns,odr0);%NCOS
    time(iter) = toc;
    dag{1} = dag_gbn_ncos;
    score(iter) = score_dags(data,ns,dag,'scoring_fn','bic');
end
% figure(2)
% draw_graph(dag)



