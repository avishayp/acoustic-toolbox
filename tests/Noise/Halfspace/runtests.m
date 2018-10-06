% tests for halfspace
% If the source is put at lambda/4 then
% the far field should have a level 13.2 dB higher than rho_SL_dB
% which is 10 log10( 6.6 * pi )
%
% See the Kewley paper for a derivation of that based on a sheet of
% monopoles in a halfspace

GreenFile = 'halfspace.grn';

scooter halfspace

sd        = 1.25;  % depth of noise sources
rd        = 0:1:100;
Component = 'P';
rho_SL_dB = 54;   % monopole strength from Kewley for 800 Hz and 40 knots
freq      = 300;

NL = spectral_noise( GreenFile, rho_SL_dB, sd, rd, freq );

figure
plot( rd, NL, 'LineWidth', 3 )
xlabel( 'Receiver depth (m)' )
ylabel( 'NL (dB)' )
title( 'Halfspace, spectral formula' )

