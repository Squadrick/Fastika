function dispsig(signalMatrix, range, titlestr);

if nargin < 3, titlestr = ''; end
if nargin < 2, range = 1:size(signalMatrix, 1); end