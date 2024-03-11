function [pval_pns,z_scores, rank] = neujive_test_hipp_caud(res_hipp, res_caud, rank)
% load('aligned_data.mat'); 
% setup;
% all_hipp = [pos_hipp_pdms, neg_hipp_pdms];
% all_caud = [pos_caud_pdms, neg_caud_pdms];
% [normed_hipp, ~, ~] = normalize_data(all_hipp);
% [normed_caud, ~, ~] = normalize_data(all_caud);
% 
% [res_hipp, PNS_hipp] = PNSmain(normed_hipp, 2);
% [res_caud, PNS_caud] = PNSmain(normed_caud, 2);

caBlocks_hipp_caud=[];
caBlocks_hipp_caud{1} = res_hipp;
caBlocks_hipp_caud{2} = res_caud;

caDataSetNames = {'hippocampus', 'caudate'} ;
paramstruct = struct('dataname',{caDataSetNames}, ...
                'iplot', [0 0],...
                           'nresample',1000, ...
                           'ioutput',[1 1 1 1 1 1 1 1 1]) ;
vRank = [rank rank];
outstruct = AJIVEMainMJ(caBlocks_hipp_caud,vRank,paramstruct);
hipp_joint = outstruct.BSSjoint{1, 1};
if isempty(hipp_joint)
    pval_pns = 1;
    z_scores = 0;
    return
end

mat_pos_pns = hipp_joint(:, 1:34); mat_neg_pns = hipp_joint(:, 35:end);
[~,pval_pns,z_scores] = DiProPermSM(mat_pos_pns,...
    mat_neg_pns); 
% save(['../results/neujive_test_hipp_caud_', num2str(rank), '.mat']);
end