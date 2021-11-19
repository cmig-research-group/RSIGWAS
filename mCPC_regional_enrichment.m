% This script is to testout the idea on getting the regional enrichment analysis

modality = 'ND'


% load ROI probabilities in the atlas space 
% which should contains a v by k matrix, atlas_probvec, v as voxels, k as lables. 

load('~/ukb_atlas_label_prob.mat');

% load summary statistics, including eigen vectors V and z statistics
load(sprintf('~/volmat_%s_k_5000.discovery.summary.mat', modality));

vvec = V*zmats(survive, :)';
expect_v = vvec'*atlas_probvec./sum(atlas_probvec, 1); 

% sampled null SNPs, assuming exchangeable
svec = find(pmat(:,2) > 5e-5);
sarray = 1:50:size(svec,1); % even sampling
svec = svec(sarray);
Titer = length(svec);
null_v = zeros(Titer, size(expect_v,2));
nblock = 1000; % sampling 1000 at a time 
eprob = atlas_probvec./sum(atlas_probvec, 1);
for iiter = 0:nblock:Titer;
  tic;
  if (iiter + nblock) >=Titer;
    iseq = ((iiter+1):Titer);
    vmats = V*zmats(svec(iseq),:)';
  else
    iseq = ((iiter+1):(iiter+nblock));
    vmats = V*zmats(svec(iseq), :)';
  end
  fprintf('%d of %d - ', iiter, Titer);
  null_v(iseq,:) = vmats'*eprob;
  toc
end

mu = mean(null_v,1);sd = std(null_v,1);
enrich_z = (expect_v - mu)./sd;
