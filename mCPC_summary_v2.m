function [survive, pvec, sumstat, h2_snp, h2_ldsc, h2_ldsc_se] = mCPC_summary_v2(info, zmats, pvec, nvec, ds, pthresh, r2thresh) 

% Multistep inference from Aschard et al., 2014
% 
% info_path for the mCPC processing
% globalK for multi-component combination function 
% pthresh for GWAS significance threshold
% r2thresh for LD prunning threshold
%

fprintf('%s -- %s.m: Loading ldsc info \r\n', datestr(now), mfilename);

load('~/UKB_9m_filtered_l2.mat'); % This is a curated ldscores parsed from the reference LD of LDscore regression repo. No need to run this one if you don't need the multivariate regression of LDSC. For file format requirement, see mvn_LDSC.m

fprintf('%s -- %s.m: Pruning  with p in %f r2 in  %f \r\n', datestr(now), mfilename, pthresh, r2thresh);

survive = prune(info.geno_path, pvec, pthresh, r2thresh);

chivec = sum(zmats(survive, :).^2,2);
r2 = (zmats(survive,:).^2)./(zmats(survive,:).^2 + repmat(nvec(survive) - 2, 1, size(zmats,2)));


g = (ds.^2)./sum(ds.^2);
h2_snp = sum(r2*g);

[coef_int, coef_b, est_h2, se_h2, jk_est, jk_se] = mvn_LDSC(zmats(logical(invec),:).^2, nvec(logical(invec)), []);

h2_ldsc = est_h2'*g;
h2_ldsc_se = se_h2*g;

sumstat = table(info.bim{1}(survive), info.bim{2}(survive), info.bim{4}(survive), info.bim{5}(survive), info.bim{6}(survive), nvec(survive), chivec, pvec(survive), -log10(pvec(survive)), 'VariableNames', {'Chr', 'SNP', 'BP', 'A1','A2', 'N', 'Chi2','P', 'nlog10_P'});



