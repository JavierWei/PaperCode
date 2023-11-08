function [dag,scor,cache] = yy_HC_order_none(data,max_pares,ns,odr0)
% ���� order ����ɽ��

global yyfieldnames yyindex1 yyindex2
yyfieldnames = '010203040506070809101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960';
yyindex1 = 1:2:59;
yyindex2 = 2:2:60;
[nnds,nsmpls] = size(data);
if nargin<2, max_pares = floor(log(2*nsmpls/log(nsmpls))); end;
if nargin<3, ns = max(data'); end;
if nargin<4, odr0 = randperm(nnds); end;

[dag0,scor_best,cache] = yy_full_order_optimal_dag(data,max_pares,ns,odr0);
%fprintf('��ʼ���ŷ�����%f...\n',scor_best);
done = false;
LongFor = nnds*(nnds-1)/2;
dag_best = dag0;
count = 0;
while ~done
    count = count+1;
 %   fprintf('-----------------------------\n');
 %   fprintf('����������%d...\n',count);
    fam_scors0 = yy_family_scor_dag(data,dag0,ns);
    neigb_index = [1, 1];
    neigb_index_best = neigb_index;
    for t = 1:LongFor
        neigb_index = next_neigb_order(neigb_index,nnds);
        [dag_tmp,scor_tmp,cache] = yy_full_order_optimal_dag(data,max_pares,ns,odr0,neigb_index,dag0,fam_scors0,cache);
        if scor_tmp>scor_best
            dag_best = dag_tmp;
            scor_best = scor_tmp;
            neigb_index_best = neigb_index;
   %         fprintf('�����ڵ�λ�ã�%d<->%d��---->',neigb_index_best);
  %          fprintf('���õķ�����%f...\n',scor_best);
        end
    end
    if isequal(neigb_index_best,[1,1])
        done = true;
    else
        odr0 = exchange_order(odr0,neigb_index_best);
        dag0 = dag_best;
    end
end
dag = dag_best;
scor = scor_best;

%******************************************************************************************************************

function [dag,scor,cache] = yy_full_order_optimal_dag(data,max_pares,ns,odr0,neigb_index,dag0,fam_scors0,cache)
% �� order0 �� neigb_index ��˳���������Žṹ
if nargin<5
    dag = zeros(numel(odr0));
    scor = 0;
    cache = cell(1,numel(odr0));    % �洢��������֣�ÿ���ڵ������ýṹ��洢
    for i = 1:numel(odr0)
        prends = sort(odr0(1:i-1));
        best_fam_scor = yy_fam_scor(data(odr0(i),:),ns(odr0(i)));    % �޸��ڵ�
        cache{odr0(i)}.(assign_fieldname([])) = best_fam_scor;    % �洢����
        best_pares = [];
        for t = 1:min(numel(prends),max_pares)
            ppares = yy_nchoosek(prends,t);
            fam_scors = zeros(size(ppares,1),1);
            for l = 1:size(ppares,1)
                fam_scors(l) = yy_fam_scor(data([odr0(i),ppares(l,:)],:),ns([odr0(i),ppares(l,:)]));
                cache{odr0(i)}.(assign_fieldname(ppares(l,:))) = fam_scors(l);    % �洢����
            end
            [best_scor_tmp,loc] = max(fam_scors);
            if best_scor_tmp>best_fam_scor
                best_fam_scor = best_scor_tmp;
                best_pares = ppares(loc,:);
            end
        end
        dag(best_pares,odr0(i)) = 1;
        scor = scor+best_fam_scor;
        cache{odr0(i)} = orderfields(cache{odr0(i)});
    end
else
    dag = dag0;
    scor = 0;
    odr = exchange_order(odr0,neigb_index);
    for i = 1:numel(odr)
%         if i<neigb_index(1) || i>neigb_index(2)
%             scor = scor+fam_scors0(odr(i));
%         else
            prends = sort(odr(1:i-1));
            [best_fam_scor,cache{odr(i)}] = yy_fam_scor_cache(data(odr(i),:),ns(odr(i)),cache{odr(i)},assign_fieldname([]));
%             best_fam_scor = yy_fam_scor(data(odr(i),:),ns(odr(i)));    % �޸��ڵ�
            best_pares = [];
            for t = 1:min(numel(prends),max_pares)
                ppares = yy_nchoosek(prends,t);
                fam_scors = zeros(size(ppares,1),1);
                for l = 1:size(ppares,1)
                    [fam_scors(l),cache{odr(i)}] = yy_fam_scor_cache(data([odr(i),ppares(l,:)],:),ns([odr(i),ppares(l,:)]),cache{odr(i)},assign_fieldname(ppares(l,:)));
%                     fam_scors(l) = yy_fam_scor(data([odr(i),ppares(l,:)],:),ns([odr(i),ppares(l,:)]));
                end
                [best_scor_tmp,loc] = max(fam_scors);
                if best_scor_tmp>best_fam_scor
                    best_fam_scor = best_scor_tmp;
                    best_pares = ppares(loc,:);
                end
            end
            dag(:,odr(i)) = 0;
            dag(best_pares,odr(i)) = 1;
            scor = scor+best_fam_scor;
%         end
        cache{odr0(i)} = orderfields(cache{odr0(i)});
    end
end


%******************************************************************************************************************

function fam_scors = yy_family_scor_dag(data,dag,ns)
% ���������������
fam_scors = zeros(1,numel(ns));
for i = 1:numel(fam_scors)
    pares_i = find(dag(:,i));
    fam_scors(i) = yy_fam_scor(data([i; pares_i],:),ns([i; pares_i]));
end

%******************************************************************************************************************

function neigb_index = next_neigb_order(neigb_index,nnds)
% ������һ��
neigb_index = yy_next_combination(neigb_index,nnds);


%******************************************************************************************************************

function odr0 = exchange_order(odr0,index)
% ���� index ����λ��
if isempty(index), return; end;
odr0(index) = odr0(fliplr(index));

%******************************************************************************************************************

function [scor,cache] = yy_fam_scor_cache(data_fam,ns_fam,cache,fname)
% �� cache �Ķ�֮���޵���
% try
%     scor = cache.(fname);
% catch
    scor = yy_fam_scor(data_fam,ns_fam);
%     cache.(fname) = scor;
% end
    
%******************************************************************************************************************

function fname = assign_fieldname(A)
% ��һ������ A ת���ɺϷ����ַ����ṹ�� fieldname
% fname = num2str(10*A);
% fname(fname==' ') = '';
% fname = ['f' fname];
global yyfieldnames yyindex1 yyindex2
fname = ['f' yyfieldnames(sort([yyindex1(A),yyindex2(A)]))];


