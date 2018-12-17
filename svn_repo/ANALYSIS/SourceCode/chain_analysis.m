function chain_analysis(chains)

top_x = 10; % show the top x

[uchains,m,n] = unique(chains);
occ = hist(n,max(n));
[occ,iocc] = sort(occ,'descend');
top_x = min(top_x,length(occ));
% occ_thresh = mean(occ);
% ind = find(occ > occ_thresh);
iocc = iocc(1:top_x);
rchains = uchains(iocc);

