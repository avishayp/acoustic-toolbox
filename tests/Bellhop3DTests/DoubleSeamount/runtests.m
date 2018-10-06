% run the double seamount test case (created by YT Lin)

global units
units = 'km';

%%
bellhop3d DoubleSeamount3D

%%
figure
plotshdpol( 'DoubleSeamount3D.shd',  0, 0, 400 )
caxis( [ 40 80 ] )
axis equal
axis( [ 0 6 -4 4 ] )

%print -dpng DoubleSeamount3D
