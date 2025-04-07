function linear_mixed_models()
load('results.mat','tbl_scatter')
%% Linear Mixed Models
% Using polynomial estimates
delete(['tables' filesep 'LME.txt'])
diary(['tables' filesep 'LME.txt'])
disp('----------------------------------------------------------------------')
disp('**********************************************************************')
disp('----------------------------------------------------------------------')
disp(' ')
disp('Using polynomial estimates:')
disp(' ')

% fit full LMM
tbl_scatter.SI = categorical(tbl_scatter.SI);
lmm             = fitlme(tbl_scatter,'ES ~ SI + (1+SI|obsid) + (1+SI|codenum) + (1+SI|codenum:obsid)');
disp(lmm)
% fit reduced LMM 
tbl_scatter.new_SI = repmat(categorical([1 2 1]'),height(tbl_scatter)/3,1);
reduced_lmm = fitlme(tbl_scatter, 'ES ~ new_SI + (1+new_SI|obsid) + (1+new_SI|codenum) + (1+new_SI|codenum:obsid)');

% to compare the models
bic_full = lmm.ModelCriterion.BIC;
bic_reduced = reduced_lmm.ModelCriterion.BIC;
dBIC = bic_reduced-bic_full;
BF_reduced_over_full = 1/exp((bic_reduced-bic_full)/2);
disp(['dBIC = ' num2str(dBIC) ', and BF_reduced_over_full = ' num2str(BF_reduced_over_full)])

% Using model-free estimates
disp(' ')
disp('----------------------------------------------------------------------')
disp('**********************************************************************')
disp('----------------------------------------------------------------------')
disp(' ')
disp('Using model-free estimates:')
disp(' ')
lmm             = fitlme(tbl_scatter,'bin_scatter ~ SI + (1+SI|obsid) + (1+SI|codenum) + (1+SI|codenum:obsid)');
disp(lmm)
tbl_scatter.new_SI = repmat(categorical([1 2 1]'),height(tbl_scatter)/3,1);
reduced_lmm = fitlme(tbl_scatter, 'bin_scatter ~ new_SI + (1+new_SI|obsid) + (1+new_SI|codenum) + (1+new_SI|codenum:obsid)');

bic_full = lmm.ModelCriterion.BIC;
bic_reduced = reduced_lmm.ModelCriterion.BIC;
dBIC = bic_reduced-bic_full;
BF_reduced_over_full = 1/exp((bic_reduced-bic_full)/2);

disp(['dBIC = ' num2str(dBIC) ', and BF_reduced_over_full = ' num2str(BF_reduced_over_full)])
diary off;

end