function [Out1, Out2] = applyW(mixedsig, W, ...)
if nargin==2
  [mixedsig, meanvalues] = remmean(mixedsig);
  icasig = W * mixedsig;
elseif (nargin<2)
  error('Not enough arguments');
  exit(-1);
elseif (nargin>3)
  error('Too many arguments');
  exit(-1);
elseif nargin==3
  param = va_arg();
  if strcmp('mean', param)
    icasig = W * mixedsig;
  else
    error(['Unrecognized parameter: ''' param '''']);
    exit(-1);
  endif
endif

# Determine what to output
if nargout == 1
  Out1 = icasig;
else
  Out1 = icasig;
  Out2 = meanvalues;
endif
endfunction





