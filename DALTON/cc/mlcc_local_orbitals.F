      subroutine mlcc_write_local_orbital_info(cmo)
C
C
C
C
#include "implicit.h"     
#include "ccorb.h"
#include "dummy.h"
#include "priunit.h"
#include "mxcent.h"
#include "center.h"
C
         double precision cmo(nbast,norbt)
C  
         if (nsym .ne. 1) then
            write(lupri,'(/A,A/)')'.MLLOC not compatible with',
     &                      ' calculation using symmetry!'
            call quit('')
         endif
C
C
        mlunit = -1
         call gpopen(mlunit, 'LOCAL_ORB_INFO','UNKNOWN','SEQUENTIAL',
     &        'FORMATTED',idummy,.false.)
         rewind(mlunit)
         write(mlunit,*)nspace
         write(mlunit,*)(nspaoc(1,i), i = 1,nspace)
         write(mlunit,*)(nspavi(1,i), i = 1,nspace)
         write(mlunit,*)ninaoc(1), ninavi(1)
C
         write(mlunit,*)cmo
C
        call gpclose(mlunit, 'KEEP')
C
      end subroutine