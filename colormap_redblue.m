function cmap = colormap_redblue(n)
if nargin < 1, n = 256; end
half = floor(n/2);
r1 = linspace(0.10, 1.0, half)';
g1 = linspace(0.25, 1.0, half)';
b1 = linspace(0.85, 1.0, half)';
r2 = linspace(1.0, 0.85, n-half)';
g2 = linspace(1.0, 0.15, n-half)';
b2 = linspace(1.0, 0.10, n-half)';
cmap = [[r1;r2], [g1;g2], [b1;b2]];
end