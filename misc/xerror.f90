      SUBROUTINE XERROR( MESSG, NMESSG, NERR, LEVEL )

!     Copper penny fuse

!     The original xerror optionally prints a warning message
!     Depending on severity, previous occurence, ...

!     This is a dumb version which always prints the error message
!     Michael B. Porter 4/27/86

      CHARACTER (LEN=72) :: MESSG
      INTEGER NMESSG, NERR, LEVEL

      WRITE( 6, '('' '',A72)' ) MESSG

      RETURN
      END
