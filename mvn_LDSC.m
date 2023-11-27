function [coef_int, coef_b, est_h2, se_h2, jk_est, jk_se] = mvn_LDSC(zz, nvec, ldsc_annot_path);

% This is multiple dependent variable version of LDSC
% 
% The main difference is the weighting factor
% Because this relied on matrix multiplication to get 1000s of outcomes done in one instance
% the reciprocal weights cannot be implemented. 
% However, the difference should be miniscule
%
% Here also assume ldsc annotation files have been pre-intersect with zmats
% therefore, require L2=sumstat.ldsc$L2, invec=invec, chrvec=sumstat.ldsc$CHR, CM=sumstat.ldsc$CM, MAF=sumstat.ldsc$MAF, BP=sumstat.ldsc$BP, in the file. With invec as the subset indices for zmats
%
% while zz is the chisqr stats (z.^2) which already aligned with ldsc annots (invec)

nblocks = 200; % default number of block jackknife 

if isempty(ldsc_annot_path);
  ldsc_annot_path = '~/UKB_9m_filtered_l2.mat'; % now this is fix parameter
end

%
fprintf('%s -- %s.m: Loading LD annotations from %s \r\n', datestr(now), mfilename, ldsc_annot_path);

load(ldsc_annot_path);

% generate the weights 
l = 1./L2;
l(L2 <1) = 1;

% get chisquare
total_SNP = length(l);
total_N = mean(nvec);

fprintf('%s -- %s.m: Calculating Projection Matrix \r\n', datestr(now), mfilename);
W = sparse(1:total_SNP, 1:total_SNP, l, total_SNP, total_SNP);
X = cat(2, ones(total_SNP,1), L2);
XWX = X'*W*X;
XWY = X'*W*zz;

fprintf('%s -- %s.m: Calculating h2 \r\n', datestr(now), mfilename);
coef = inv(XWX)*XWY;

coef_int = coef(1,:)';
coef_b = coef(2:end,:)';
est_h2 = coef_b*total_SNP./total_N;

% TODO: get block JK
% See fast block JK in ldsc, the default is 200 blocks


fprintf('%s -- %s.m: Performing Block Jackknife with default %d blocks \r\n', datestr(now), mfilename, nblocks);

blockvec = int64(linspace(0, total_SNP, nblocks + 1));
XWX_blocks = zeros([size(XWX),nblocks]);
XWY_blocks = zeros([size(XWY),nblocks]);

for iblock = 1:nblocks;
  fprintf('       - %d of %d block \r\b', iblock, nblocks);
  XWX_blocks(:,:,iblock) = X((blockvec(iblock) + 1):blockvec(iblock + 1),:)'*W((blockvec(iblock) + 1):blockvec(iblock + 1),(blockvec(iblock) + 1):blockvec(iblock + 1))*X((blockvec(iblock) + 1):blockvec(iblock + 1),:);  
  XWY_blocks(:,:,iblock) = X((blockvec(iblock) + 1):blockvec(iblock + 1),:)'*W((blockvec(iblock) + 1):blockvec(iblock + 1),(blockvec(iblock) + 1):blockvec(iblock + 1))*zz((blockvec(iblock) + 1):blockvec(iblock + 1),:);
end

XWX_total = sum(XWX_blocks, 3);
XWY_total = sum(XWY_blocks, 3);
delete_values = zeros([size(XWY),nblocks]);
est = inv(XWX_total)*XWY_total;
for iblock = 1:nblocks;
  fprintf('       - %d of %d block \r\b', iblock, nblocks);
  tmpxx = XWX_total - XWX_blocks(:,:,iblock);
  tmpxy = XWY_total - XWY_blocks(:,:,iblock);
  delete_values(:,:,iblock) = inv(tmpxx)*tmpxy;
end
pseudo_values = nblocks*est - (nblocks - 1) * delete_values;
jk_est = mean(pseudo_values,3);
jk_se = sqrt(var(pseudo_values, [], 3)/nblocks);


% TODO: double check the inferential stats
se_h2 = jk_se(2:end,:)*total_SNP./total_N;



