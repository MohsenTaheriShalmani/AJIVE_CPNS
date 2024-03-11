clear; close all;
% load('euclideanized_pdms_3_structure.mat');
setup;
load('aligned_data.mat');
num_pos = size(pos_hipp_pdms, 2);
num_neg = size(neg_hipp_pdms, 2);
all_hipp = [pos_hipp_pdms, neg_hipp_pdms];
all_caud = [pos_caud_pdms, neg_caud_pdms];

%% PNS 
[normed_hipp, ~, ~] = normalize_data(all_hipp);
[normed_caud, ~, ~] = normalize_data(all_caud);

[res_hipp, PNS_hipp] = PNSmain(normed_hipp, 9);
[res_caud, PNS_caud] = PNSmain(normed_caud, 9);
pvals = [];
zscores = [];

%% SVD
% [U_hipp, s_hipp, v_hipp] = svds(all_hipp, 61);
% [u_caud, s_caud, v_caud] = svds(all_caud, 61);
% hipp_feat = U_hipp' * all_hipp;
% caud_feat = u_caud' * all_caud;

%% Euclidean + CPCA
% concate_feat = [res_hipp;res_caud];
% [u_cpca, s_cpca, v_cpca] = svds(concate_feat, 61);
% cut_concate_feat = u_cpca' * concate_feat;
% 
% all_pdm = [hipp_feat; caud_feat];
% mat_pos_pns = cut_concate_feat(:, 1:34); mat_neg_pns = cut_concate_feat(:, 35:end);
% [~,pval_pns,z_scores] = DiProPermSM(mat_pos_pns,...
%     mat_neg_pns);

i=42; % from scree plots see AJIVE paper
    [pval, zscore, r] = neujive_test_hipp_caud(res_hipp, res_caud, i);
