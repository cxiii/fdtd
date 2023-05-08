% Define Custom COLORMAP
CMAP = zeros(256,3);
blue = [0 0 1];
white = [1 1 1];
red = [1 0 0];

for nc = 1:128
  f = (nc-1)/128;
  c = (1-sqrt(f))*blue + sqrt(f)*white;
  CMAP(nc,:) = c;
  c = (1-f^2)*white + f^2*red;
  CMAP(128+nc,:) = c;
endfor

%colormap(CMAP);
