function [dag,best_score] = yy_learn_struct_hc(data, nodesizes, seeddag)
%
% LEARN_STRUCT_HC(data,seeddag) learns a structure of Bayesian net by Hill Climbing.
% dag = learn_struct_hc(data, nodesizes, seeddag)
%
% dag: the final structurre matrix
% Data : training data, data(i,m) is the m obsevation of node i
% Nodesizes: the size array of different nodes
% seeddag: given seed Dag for hill climbing, optional
%
% by Gang Li @ Deakin University (gli73@hotmail.com)

[N ncases] = size(data);
if (nargin < 3 ) 
    seeddag = zeros(N,N); % mk_rnd_dag(N); %call BNT function
elseif ~acyclic(seeddag)
    seeddag = mk_rnd_dag(N); %zeros(N,N);
end;
Max_iter = 5000;

done = 0;
best_score = yy_score_dag(data, seeddag, nodesizes);    % Íâ²¿º¯Êý
iter = 0;
while ~done && iter<=Max_iter
    iter = iter+1;
    [dags,op,nodes] = mk_nbrs_of_dag(seeddag);
    nbrs = length(dags);
    scores = yy_score_dags(best_score, seeddag, data, nodesizes, dags,op,nodes);
    max_score = max(scores);
    new = find(scores == max_score );
    if ~isempty(new) & (max_score > best_score)
        p = sample_discrete(normalise(ones(1, length(new))));
        best_score = max_score;
        seeddag = dags{new(p)};
    else
        done = 1;
    end;
end;

dag = seeddag;

%******************************************************************************************************************

function scores = yy_score_dags(best_score, seeddag, data, ns, dags, op, nodes)

N = size(data,1);
orglocscors = zeros(1,N);
for t = 1:N
    orglocscors(t) = yy_fam_scor(data([t; find(seeddag(:,t))],:),ns([t; find(seeddag(:,t))]));
end
scores = zeros(1,length(dags));
for t = 1:length(dags)
    switch op{t}
        case 'del'
            scores(t) = best_score - orglocscors(nodes(t,2));
        case 'rev'
            scores(t) = best_score - sum(orglocscors(nodes(t,:)));
            scores(t) = scores(t) + yy_fam_scor(data([nodes(t,1); find(dags{t}(:,nodes(t,1)))],:),ns([nodes(t,1); find(dags{t}(:,nodes(t,1)))]));
        case 'add'
            scores(t) = best_score - orglocscors(nodes(t,2));
    end
    scores(t) = scores(t) + yy_fam_scor(data([nodes(t,2); find(dags{t}(:,nodes(t,2)))],:),ns([nodes(t,2); find(dags{t}(:,nodes(t,2)))]));
end



