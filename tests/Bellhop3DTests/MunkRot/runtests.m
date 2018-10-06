% run the 3D Munk case

global units
units = 'km';

    
bellhop3d Munk3D

figure
plotshdpol( 'Munk3D.shd', 0, 0, 3000 )
caxisrev( [ 50 100 ] )
axis( [ -1 4 0 100 ] )

bellhop3d Munk3Dz

figure
plotshd( 'Munk3Dz.shd' )
caxisrev( [ 50 100 ] )

