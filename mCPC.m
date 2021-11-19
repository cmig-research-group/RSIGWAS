function [zmats, nvec, freqvec] = mCPC(info); 

dbstop if error

% This is the function for using combined PC approach
% reference:
% Aschard et al., 2014 AJHG
%
example_flag = 0;

if example_flag == 1;

info = struct();
info.vol_path = '~/volmat_N0.mat';
info.subj_info_path = '~/basic_info_rsi.mat';
info.covar_path = '~/ukb_covar_design.mat';
info.geno_path = '~/UKB42k_QCed_160320_filtered';
info.k = 2000;
info.out_path = ['~/results/N0_k_' num2str(info.k) '.mat'];
info.chunk = 5000; % snp chunk to process

end


% print out the parameters
info


% STEP1: Processing phenotypes
%
% NOTE: for generating PC, there is argument not performing it on preresidualized data
%       as it would introduce biases (esp. collider bias).
%       Given the idea about combining PCs, it is better to do covariate control
%       after PC projection rather than before.
%
%

fprintf('%s -- %s.m: Reading imaging data from %s \r\n', datestr(now), mfilename, info.vol_path);

load(info.vol_path);

fprintf('%s -- %s.m: Reading imaging info from %s \r\n', datestr(now), mfilename, info.subj_info_path);

load(info.subj_info_path);

fprintf('%s -- %s.m: Reading the covar from %s \r\n', datestr(now),mfilename, info.covar_path);

covar = load(info.covar_path); 
covar_id = num2str(covar.ID);

fprintf('%s -- %s.m: Intersecting phenotypic data \r\n', datestr(now), mfilename);

[IID, ia, ib] = intersect(idlist, covar_id, 'stable');
R = covar.R(ib,:);
volmat = volmat(ia,:);

fprintf('\t\t\tFound %d individuals in the volmat \r\n \t\t\tFound %d individuals in the covar \r\n \t\t\tFound %d individuals intersect \r\n', length(idlist), length(covar.ID), length(IID));

fprintf('%s -- %s.m: Reading FAM info from %s \r\n', datestr(now), mfilename, sprintf('%s.fam', info.geno_path));
fileID = fopen(sprintf('%s.fam', info.geno_path));
fam_file = textscan(fileID,'%s %s %s %s %s %s');
fclose(fileID);
fam_nsubj = length(fam_file{1});

fprintf('%s -- %s.m: Intersecting genotype data \r\n', datestr(now), mfilename);
[IID2, ip, ig] = intersect(IID, fam_file{1}, 'stable');
fprintf('\t\t\tFound %d individuals in the fam \r\n \t\t\tFound %d individuals intersect \r\n', length(fam_file{1}), length(IID2));
R = R(ip,:);
volmat = volmat(ip,:);

defvec = isfinite(sum(R,2) + sum(volmat, 2));
fprintf('\t\t\tFound %d individuals not missing phenotypes \r\n', sum(defvec));
R = R(defvec,:);
volmat = volmat(defvec, :);
IID2 = IID2(defvec);
ig = ig(defvec);

fprintf('%s -- %s.m: Calculating PCs with num of K in %d \r\n', datestr(now), mfilename, info.k);
mean_vol = mean(volmat,1);
volmat = double(volmat - mean_vol);
[U S V] = svdecon(volmat, info.k);
ds = diag(S);
uv = var(U);
clear volmat; % release the memory load

fprintf('%s -- %s.m: Residualizing the PC scores \r\n', datestr(now), mfilename);
Ri = pinv(R);
beta = U'*Ri';Ures = U - R*beta';


% STEP2: Loop over genome blocks

fprintf('%s -- %s.m: Reading SNP info from %s \r\n', datestr(now), mfilename, sprintf('%s.bim', info.geno_path));
fileID = fopen(sprintf('%s.bim', info.geno_path));
bim_file = textscan(fileID,'%s %s %s %s %s %s');
fclose(fileID);
snps=length(bim_file{1});
fprintf('\t\t\t Found %d SNPs \r\n', snps);


zmats = NaN(snps, info.k);
nvec=zeros(snps, 1, 'single');
freqvec=zeros(snps, 1, 'single');

for i=1:info.chunk:snps
    j=min(i+info.chunk-1, snps);
    fprintf('gwas: loading snps %i to %i... ', i, j);    tic;
    geno_int8 = PlinkRead_binary2(fam_nsubj, i:j, info.geno_path);
    fprintf('processing... ', i, j);
    geno = nan(size(geno_int8), 'single'); for code = int8([0,1,2]), geno(geno_int8==code) = single(code); end;
    geno = geno(ig, :);
    [~, zmat_orig_chunk] = nancorr(Ures, geno);
    zmats(i:j,:) = zmat_orig_chunk';
    nvec(i:j) = sum(isfinite(geno))';
    freqvec(i:j) = (1*sum(geno==1) + 2*sum(geno==2))' ./ (2*nvec(i:j));
    fprintf('done in %.1f sec, %.1f %% completed\n', toc, 100*(j+1)/snps);
end



% STEP3: Save results for testing out different multisteps combination functions
fprintf('%s -- %s.m: saving results \r\n', datestr(now), mfilename, info.out_path);
save(info.out_path, '-v7.3', 'zmats','nvec','freqvec', 'info','V', 'mean_vol', 'ds','uv');
fprintf('%s -- %s.m: Done. \r\n', datestr(now), mfilename);




