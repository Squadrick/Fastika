function [A, W] = fpica(X, whiteningMatrix, dewhiteningMatrix, approach, ...
			numOfIC, g, finetune, a1, a2, myy, stabilization, ...
			epsilon, maxNumIterations, maxFinetune, initState, ...
			guess, sampleSize, displayMode, displayInterval, ...
			s_verbose, plottype);

if nargin < 3, error('Not enough arguments!'); end
[vectorSize, numSamples] = size(X);
if nargin < 20, s_verbose = 'on'; end
if nargin < 19, displayInterval = 1; end
if nargin < 18, displayMode = 'on'; end
if nargin < 17, sampleSize = 1; end
if nargin < 16, guess = 1; end
if nargin < 15, initState = 'rand'; end
if nargin < 14, maxFinetune = 100; end
if nargin < 13, maxNumIterations = 1000; end
if nargin < 12, epsilon = 0.0001; end
if nargin < 11, stabilization = 'on'; end
if nargin < 10, myy = 1; end
if nargin < 9, a2 = 1; end
if nargin < 8, a1 = 1; end
if nargin < 7, finetune = 'off'; end
if nargin < 6, g = 'pow3'; end
if nargin < 5, numOfIC = vectorSize; end
if nargin < 4, approach = 'defl'; end

b_verbose = 0;
approachMode = 1;
gOrig=10;
stabilizationEnabled = 0;
finetuningEnabled = 0;

if sampleSize > 1
  sampleSize = 1;
elseif sampleSize < 1
  if (sampleSize * numSamples) < 1000
    sampleSize = min(1000/numSamples, 1); 
  end
end

if sampleSize ~= 1
  gOrig = gOrig + 2;
end
if myy ~= 1
  gOrig = gOrig + 1;
end

if (myy ~= 1),
  gFine = gOrig;
else 
  gFine = gOrig + 1;
end

myyOrig = myy;
myyK = 0.01;
failureLimit = 5;

usedNlinearity = gOrig;
stroke = 0;
notFine = 1;
long = 0;
initialStateMode = 0;
usedDisplay = 0;

if displayInterval < 1
  displayInterval = 1;
end
if approachMode == 1
  B = zeros(vectorSize);
  round = 1;
  numFailures = 0;
  while round <= numOfIC,
    myy = myyOrig;
    usedNlinearity = gOrig;
    stroke = 0;
    notFine = 1;
    long = 0;
    endFinetuning = 0;
    if initialStateMode == 0
      w = rand(vectorSize, 1) - .5;
    elseif initialStateMode == 1
      w=whiteningMatrix*guess(:,round);
    end
    w = w - B * B' * w;
    w = w / norm(w);
    
    wOld = zeros(size(w));
    wOld2 = zeros(size(w));
    i = 1;
    gabba = 1;
    while i <= maxNumIterations + gabba
      w = w - B * B' * w;
      w = w / norm(w);
      
      if notFine
	if i == maxNumIterations + 1
	  round = round - 1;
	  numFailures = numFailures + 1;
	  if numFailures > failureLimit
	    if round == 0
	      A=[];
	      W=[];
	    end
	    return;
	  end
	  break;
	end
      else
	if i >= endFinetuning
	  wOld = w;
	end
      end
      if norm(w - wOld) < epsilon | norm(w + wOld) < epsilon
        if finetuningEnabled & notFine
          notFine = 0;
	  gabba = maxFinetune;
          wOld = zeros(size(w));
          wOld2 = zeros(size(w));
          usedNlinearity = gFine;
          myy = myyK * myyOrig;
	  
	  endFinetuning = maxFinetune + i;
	  
        else
          numFailures = 0;
          % Save the vector
          B(:, round) = w;
          A(:,round) = dewhiteningMatrix * w;
          W(round,:) = w' * whiteningMatrix;
	  break;
        end

      elseif stabilizationEnabled
	if (~stroke) & (norm(w - wOld2) < epsilon | norm(w + wOld2) < ...
			epsilon)
	  stroke = myy;
	  myy = .5*myy;
	  if mod(usedNlinearity,2) == 0
	    usedNlinearity = usedNlinearity + 1;
	  end
	elseif stroke
	  myy = stroke;
	  stroke = 0;
	  if (myy == 1) & (mod(usedNlinearity,2) ~= 0)
	    usedNlinearity = usedNlinearity - 1;
	  end
	elseif (notFine) & (~long) & (i > maxNumIterations / 2)
	  long = 1;
	  myy = .5*myy;
	  if mod(usedNlinearity,2) == 0
	    usedNlinearity = usedNlinearity + 1;
	  end
	end
      end
      
      wOld2 = wOld;
      wOld = w;
      
      switch usedNlinearity
	% pow3
       case 10
	w = (X * ((X' * w) .^ 3)) / numSamples - 3 * w;
       case 11
	EXGpow3 = (X * ((X' * w) .^ 3)) / numSamples;
	Beta = w' * EXGpow3;
	w = w - myy * (EXGpow3 - Beta * w) / (3 - Beta);
       case 12
	Xsub=X(:,getSamples(numSamples, sampleSize));
	w = (Xsub * ((Xsub' * w) .^ 3)) / size(Xsub, 2) - 3 * w;
       case 13
	Xsub=X(:,getSamples(numSamples, sampleSize));
	EXGpow3 = (Xsub * ((Xsub' * w) .^ 3)) / size(Xsub, 2);
	Beta = w' * EXGpow3;
	w = w - myy * (EXGpow3 - Beta * w) / (3 - Beta);
	% tanh
       case 20
	hypTan = tanh(a1 * X' * w);
	w = (X * hypTan - a1 * sum(1 - hypTan .^ 2)' * w) / numSamples;
       case 21
	hypTan = tanh(a1 * X' * w);
	Beta = w' * X * hypTan;
	w = w - myy * ((X * hypTan - Beta * w) / (a1 * sum((1-hypTan .^2)') - Beta));
       case 22
	Xsub=X(:,getSamples(numSamples, sampleSize));
	hypTan = tanh(a1 * Xsub' * w);
	w = (Xsub * hypTan - a1 * sum(1 - hypTan .^ 2)' * w) / size(Xsub, 2);
       case 23
	Xsub=X(:,getSamples(numSamples, sampleSize));
	hypTan = tanh(a1 * Xsub' * w);
	Beta = w' * Xsub * hypTan;
	w = w - myy * ((Xsub * hypTan - Beta * w) / (a1 * sum((1-hypTan .^2)') - Beta));
	% gauss
       case 30
	u = X' * w;
	u2=u.^2;
	ex=exp(-a2 * u2/2);
	gauss =  u.*ex;
	dGauss = (1 - a2 * u2) .*ex;
	w = (X * gauss - sum(dGauss)' * w) / numSamples;
       case 31
	u = X' * w;
	u2=u.^2;
	ex=exp(-a2 * u2/2);
	gauss =  u.*ex;
	dGauss = (1 - a2 * u2) .*ex;
	Beta = w' * X * gauss;
	w = w - myy * ((X * gauss - Beta * w) / (sum(dGauss)' - Beta));
       case 32
	Xsub=X(:,getSamples(numSamples, sampleSize));
	u = Xsub' * w;
	u2=u.^2;
	ex=exp(-a2 * u2/2);
	gauss =  u.*ex;
	dGauss = (1 - a2 * u2) .*ex;
	w = (Xsub * gauss - sum(dGauss)' * w) / size(Xsub, 2);
       case 33
	Xsub=X(:,getSamples(numSamples, sampleSize));
	u = Xsub' * w;
	u2=u.^2;
	ex=exp(-a2 * u2/2);
	gauss =  u.*ex;
	dGauss = (1 - a2 * u2) .*ex;
	Beta = w' * Xsub * gauss;
	w = w - myy * ((Xsub * gauss - Beta * w) / (sum(dGauss)' - Beta));
	% skew
       case 40
	w = (X * ((X' * w) .^ 2)) / numSamples;
       case 41
	EXGskew = (X * ((X' * w) .^ 2)) / numSamples;
	Beta = w' * EXGskew;
	w = w - myy * (EXGskew - Beta*w)/(-Beta);
       case 42
	Xsub=X(:,getSamples(numSamples, sampleSize));
	w = (Xub * ((Xub' * w) .^ 2)) / size(Xsub, 2);
       case 43
	Xsub=X(:,getSamples(numSamples, sampleSize));
	EXGskew = (Xsub * ((Xsub' * w) .^ 2)) / size(Xsub, 2);
	Beta = w' * EXGskew;
	w = w - myy * (EXGskew - Beta*w)/(-Beta);
	
       otherwise
	error('Code for desired nonlinearity not found!');
      end
      
      % Normalize the new w.
      w = w / norm(w);
      i = i + 1;
    end
    round = round + 1;
  end
end

if imag(A) ~= 0
  A = real(A);
  W = real(W);
end

endfunction