function [dag,best_score] = yy_mmhc_structure_learning(data,nodesizes,alpha)
% 用 mmhc 算法学习 BN 结构

uig = yy_MMPC(data,nodesizes,alpha);

nnds = size(data,1);
seeddag = zeros(nnds);
Max_iter = 5000;

done = 0;
best_score = yy_score_dag(data, seeddag, nodesizes);    % 外部函数
iter = 0;
while ~done && iter<=Max_iter
    iter = iter+1;
    [dags,op,nodes] = yy_cnstrned_mk_nbrs_of_dag(seeddag,uig);
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

%******************************************************************************************************************

function [Gs, op, nodes] = yy_cnstrned_mk_nbrs_of_dag(G0,uig)
% 寻找由 uig 限制的 G0 的邻居
Gs = {};
op = {};
nodes = [];

[I,J] = find(G0);
nnbrs = 1;
% all single edge deletions
for e=1:length(I)
  i = I(e); j = J(e);
  G = G0;
  G(i,j) = 0;
  Gs{nnbrs} = G;
  op{nnbrs} = 'del';
  nodes(nnbrs, :) = [i j];
  nnbrs = nnbrs + 1;
end

% all single edge reversals
for e=1:length(I)
  i = I(e); j = J(e);
  G = G0;
  G(i,j) = 0;
  G(j,i) = 1;
  if acyclic(G)
    Gs{nnbrs} = G;
    op{nnbrs} = 'rev';
    nodes(nnbrs, :) = [i j];
    nnbrs = nnbrs + 1;
  end
end

[I,J] = find(uig-G0);
% all single edge additions
for e=1:length(I)
  i = I(e); j = J(e);
  if i ~= j % don't add self arcs
    G = G0;
    G(i,j) = 1;
    if G(j,i)==0 % don't add i->j if j->i exists already
      if acyclic(G)
	Gs{nnbrs} = G;
	op{nnbrs} = 'add';
	nodes(nnbrs, :) = [i j];
	nnbrs = nnbrs + 1;
      end
    end
  end
end



