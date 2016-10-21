function [Out1, Out2, Out3] = fastica(mixedsig, varargin)

[mixedsig, mixedmean] = remmean(mixedsig);
[Dim, NumOfSampl] = size(mixedsig);

verbose           = 'off';
firstEig          = 1;
lastEig           = Dim;
interactivePCA    = 'off';
approach          = 'defl';
numOfIC           = Dim;
g                 = 'pow3';
finetune          = 'off';
a1                = 1;
a2                = 1;
myy               = 1;
stabilization     = 'off';
epsilon           = 0.0001;
maxNumIterations  = 1000;
maxFinetune       = 5;
initState         = 'rand';
guess             = 0;
sampleSize        = 1;
displayMode       = 'signals';
displayInterval   = 1;
b_verbose = 1;
jumpPCA = 0;
jumpWhitening = 0;
only = 3;
userNumOfIC = 0;
plottype = 'histogram';

[E, D]=pcamat(mixedsig, firstEig, lastEig, interactivePCA, verbose);

if only > 1
    [whitesig, whiteningMatrix, dewhiteningMatrix] = whitenv ...
						     (mixedsig, E, D, verbose);
end

if only > 2
  Dim = size(whitesig, 1);
  if numOfIC > Dim
    numOfIC = Dim;
  end
  
  [A, W] = fpica (whitesig,  whiteningMatrix, dewhiteningMatrix, approach, ...
                  numOfIC, g, finetune, a1, a2, myy, stabilization, epsilon, ...
                  maxNumIterations, maxFinetune, initState, guess, sampleSize, ...
                  displayMode, displayInterval, verbose, plottype);
  
  if ~isempty(W)
    icasig = W * mixedsig;
    if b_verbose & ...
	  (max(abs(W * mixedmean)) > 1e-9) & ...
	  (strcmp(displayMode,'signals') | strcmp(displayMode,'on'))
    end
  else
    icasig = [];
  end
end

if only == 1
  Out1 = E;
  Out2 = D;
elseif only == 2
  if nargout == 2
    Out1 = whiteningMatrix;
    Out2 = dewhiteningMatrix;
  else
    Out1 = whitesig;
    Out2 = whiteningMatrix;
    Out3 = dewhiteningMatrix;
  end
else
  if nargout == 2
    Out1 = A;
    Out2 = W;
  else
    Out1 = icasig;
    Out2 = A;
    Out3 = W;
  end
end

endfunction
