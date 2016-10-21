function [E, D] = pcamat(vectors, firstEig, lastEig, s_interactive,s_verbose);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default values:
if nargin < 5, s_verbose = 'on'; end
if nargin < 4, s_interactive = 'off'; end
if nargin < 3, lastEig = size(vectors, 1); end
if nargin < 2, firstEig = 1; end
b_verbose = 0;
b_interactive=0;
oldDimension = size (vectors, 1);
if ~(b_interactive)
  if lastEig < 1 | lastEig > oldDimension
    error(sprintf('Illegal value [ %d ] for parameter: ''lastEig''\n', lastEig));
  end
  if firstEig < 1 | firstEig > lastEig
    error(sprintf('Illegal value [ %d ] for parameter: ''firstEig''\n', firstEig));
  end
end
covarianceMatrix = cov(vectors')*(size(vectors',1)-1)/size(vectors',1);
maxLastEig = rank(covarianceMatrix, 1e-9);
[E, D] = eig(covarianceMatrix);
eigenvalues = flipud(sort(diag(D)));
if lastEig > maxLastEig
  lastEig = maxLastEig;
end
if lastEig < oldDimension
  lowerLimitValue = (eigenvalues(lastEig) + eigenvalues(lastEig + 1)) / 2;
else
  lowerLimitValue = eigenvalues(oldDimension) - 1;
end
lowerColumns = diag(D) > lowerLimitValue;
if firstEig > 1
  higherLimitValue = (eigenvalues(firstEig - 1) + eigenvalues(firstEig)) / 2;
else
  higherLimitValue = eigenvalues(1) + 1;
end
higherColumns = diag(D) < higherLimitValue;
selectedColumns = lowerColumns & higherColumns;
E = selcol (E, selectedColumns);
D = selcol (selcol (D, selectedColumns)', selectedColumns);
endfunction