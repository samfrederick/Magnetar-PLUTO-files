#define  PHYSICS                        HD
#define  DIMENSIONS                     2
#define  COMPONENTS                     2
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     POTENTIAL
#define  FORCED_TURB                    NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  DIMENSIONAL_SPLITTING          NO
#define  NTRACER                        0
#define  USER_DEF_PARAMETERS            0
#define INTERNAL_BOUNDARY               YES
/* -- physics dependent declarations -- */

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */


/* [Beg] user-defined constants (do not change this line) */
// #define  VTK_VECTOR_DUMP                YES
#define  UNIT_DENSITY                   1e10
#define  UNIT_LENGTH                    1e5
#define  UNIT_VELOCITY                  1e5

/* [End] user-defined constants (do not change this line) */
