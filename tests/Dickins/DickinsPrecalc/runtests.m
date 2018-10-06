% Dickins seamount test
% mbp

% If you run KRAKEN rather than KRAKENC, the steep angle energy is way too
% low
% KRAKENC has problems if you make the phase speed limit 1e6



global units jkpsflag

units = 'km';


%%
% Coupled mode runs for range-dependent case

% run Dickins case, using the wedge ocean

D = linspace( 3000, 500, 251 );
%D = linspace( 3000, 500, 3 );

preCalcAll(   'DickinsK', D )
%pause( 200 )
fieldLoadAll( 'DickinsK', D )

figure
plotshd( 'DickinsK.shd.mat' )
caxisrev( [ 70 120 ] )
