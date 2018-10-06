function alpha = thorp( f )

% Thorp attenuation

alpha = 0.0033 + 0.11 * f.^2 ./ ( 1 + f.^2 ) + 44 * f.^2 ./ ( 4100 + f.^2 ) + 0.0003 * f.^2;

