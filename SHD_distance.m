function [n_red,n_mis,n_rev] = SHD_distance(bnetP, bnetQ)
% ���� �ṹ�������� (SHD)     �趨ǰ��ǰ��һ��Ϊ�ο� graph
% n_red ����ı�, n_mis ȱʧ�ı�, n_rev����ı�


n_red = 0 ; n_mis = 0; n_rev = 0;
% GP = bnetP.dag; GQ = bnetQ.dag;
GP = bnetP; GQ = bnetQ;
G = GP-GQ; [M,N] = find(G);
for i = 1:length(M)
  if G(M(i),N(i))==1
    switch G(N(i),M(i))
        case 0
            n_mis = n_mis+1;
        case -1
            n_rev = n_rev+1;
    end
  else
    switch G(N(i),M(i))
        case 0
            n_red = n_red+1;
        case 1
            n_rev = n_rev+1;
    end
  end
end
n_rev = n_rev/2;   % �ظ��������ʰ�֮