function stde = withinstde(X,unbiased)
% computes the within-subjects normalized standard error (Cousineau,
% 2005)
% X: columns are factors, rows are observations.
% unbiased (default = 1) will apply the Morey (2008) correction to obtain
% unbiased estimates.
% Copyright 2013 Eelke Spaak

% Edited by Claudia Luettke 8-2-2017 

% Changes: 1. comments added
%          2. added "==1" in line 32 

% Error old version: It was not possible to compute within standard error
% without applying Morey's correction.  

% If two inputs were given then the statement in line 23 is not reached.
% That's correct. However, if the input for "unbiased" is not equal to 1,
% the command in line 30 to compute unbiased values is still executed
% although it shouldn't have been.     

if nargin < 2 % If less than 2 inputs are given to function
  unbiased = 1;
end

subjMeans = nanmean(X,2);
grandMean = nanmean(subjMeans);

Y = bsxfun(@minus, X, subjMeans) + grandMean;

yvar = nanvar(Y);
if unbiased==1 % you have to specify 1, otherwise always unbiased method 
  yvar = yvar .* (size(X,2)/(size(X,2)-1));
end

stde = sqrt(yvar)./sqrt(size(X,1));