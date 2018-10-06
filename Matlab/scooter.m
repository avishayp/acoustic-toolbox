function scooter( filename )

% run the SCOOTER program
%
% usage: scooter( filename )
% where filename is the environmental file

runscooter = which( 'scooter.exe' );

if ( isempty( runscooter ) )
   error( 'scooter.exe not found in your Matlab path' )
else
   eval( [ '! "' runscooter '" ' filename ] );
end

% Fortran fields routine:
% runfields = which( 'fields.exe' );
% eval( [ '! "' runfields '" ' filename ] );

% Matlab fields routine:
fieldsco( [ filename '.grn' ] );
