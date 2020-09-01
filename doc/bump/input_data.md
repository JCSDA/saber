# Input data

Two kinds of input are required for the BUMP executable:
 * A *grid.nc* file containing the coordinates.
 * Ensemble members should have the following name: ens$E_$T_$NNNNNN.nc, where:
   * $E is the ensemble number (1 or 2)
   * $T is the timeslot name (any length)
   * $NNNNNN is the ensemble member index (six digits)
They should all be placed in the directory *datadir*, specified in the YAML input file.
