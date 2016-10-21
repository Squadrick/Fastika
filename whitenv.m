function [newVectors, whiteningMatrix, dewhiteningMatrix] = whitenv ...
    (vectors, E, D, s_verbose);
whiteningMatrix = inv (sqrt (D)) * E';
dewhiteningMatrix = E * sqrt (D);
newVectors =  whiteningMatrix * vectors;
if max (max (imag (newVectors))) ~= 0,
  error ('Whitened vectors have imaginary values.');
end
