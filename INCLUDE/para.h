c para.h
      INTEGER IPROCS, JPROCS
      INTEGER nprocs, myrank, myranki, myrankj
      INTEGER itable
      parameter(IPROCS=1, JPROCS=1)
      common/table/ itable(-1:iprocs, -1:jprocs)
      common/peinfo/ nprocs, myrank, myranki, myrankj

