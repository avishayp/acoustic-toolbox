% make a top reflection coefficient for use by Bellhop or other models

for f = 20 : 10 : 40
   ii = 1 : 151;
   U = ( ii - 1 ) / 10;
   
   for jj = ii
      y( jj ) = SL( 10, U( jj ), f );
   end
   
   figure; plot( U, y )
end
