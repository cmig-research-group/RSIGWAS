% This is an example script for performing PVS replication on ABCD data
% Assuming the data has been aligned for UKB and ABCD 
% line 20 - 30 is the calculation for PVS

modlist = {'N0','ND','NF'};

% loading mask info
load('~/volinfo.mat','ivec_mask');


for modidx = modlist;

  modality = modidx{:};

  % load results, one modality at a time
  load(sprintf('~/volmat_%s_k_5000.discovery.mat',modality),'V', 'mean_vol', 'ds','uv');
  load(sprintf('~/results/volmat_%s_k_5000.discovery.summary.mat',modality), 'zmats_sig');
  load(sprintf('~/results/volmat_%s_k_5000.ukb.rep.mat',modality), 'sumstat');

  % Generate weight for PVS
  betamat = V*zmats_sig';
  nvar = size(betamat,2);
  betamat_reg = zeros(length(ivec_mask),nvar);

  % load ABCD imaging data
  load(sprintf('~/ABCD/Images/Diffusion/volmat_%s.mat',modality));
  load('~/ABCD/Images/Diffusion/volinfo.mat', 'subjidvec','visitidvec','datevec','timevec','eventvec', 'corrmat_all');

  % calculate PVS as volproj
  volproj = (volmat - mean(volmat,'omitnan'))*betamat_reg;

  info = struct();
  info.geno_path =  '~/ABCD_genotypes/HRC/merged_ABCD'; 
  info.out_path = ['~/results/volmat_' modality '.abcd.rep.mat'];

  % add new censoring info
  info.censor = 0.8;
  invec = mean(corrmat_all,2) > info.censor
  volproj = volproj(invec,:);
  subjidvec = subjidvec(invec,:);
  eventvec = eventvec(invec,:);

  % Get snp info
  fileID = fopen(sprintf('%s.bim', info.geno_path));
  bim_file = textscan(fileID,'%s %s %s %s %s %s');
  fclose(fileID);

  % Get the subject counts
  fileID = fopen(sprintf('%s.fam', info.geno_path));
  fam_file = textscan(fileID,'%s %s %s %s %s %s');
  fclose(fileID);
  fam_nsubj = length(fam_file{1});

  % intersect SNPs to get the indices
  [snplist, ia, ib] = intersect(bim_file{2}, sumstat.SNP, 'stable');

  % read in the data
  geno_int8 = PlinkRead_binary2(fam_nsubj, ia, info.geno_path);
  geno = nan(size(geno_int8), 'single'); for code = int8([0,1,2]), geno(geno_int8==code) = single(code); end;

  flipvec = ~strcmp(bim_file{5}(ia), sumstat.A1(ib));
  geno(:,flipvec) = 2 - geno(:,flipvec);

  [IDs, iia, iib] = intersect(fam_file{1}, strcat('NDAR_', subjidvec), 'stable');

  genoid = fam_file{1};
  subjid = strcat('NDAR_', subjidvec);
  vol = volproj(:,ib);

  save(info.out_path, 'geno', 'vol','genoid','subjid','eventvec');

  % calling R lme4 packages for testing under complex design of ABCD
  status = system(sprintf('/usr/bin/Rscript --vanilla ~/ABCD_ukb_rep.R %s', info.out_path));

  load([info.out_path '.R.results.mat']);
  tab = struct2table(df);

  sumstat.abcd_t = nan(size(sumstat,1),1); sumstat.abcd_t(ib) = table2array(tab(:,4));
  sumstat.abcd_p = nan(size(sumstat,1),1); sumstat.abcd_p(ib) = table2array(tab(:,5));

  save([info.out_path '.abcd.combined.mat'], 'sumstat');
end
















 




