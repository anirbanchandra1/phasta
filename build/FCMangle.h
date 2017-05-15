#ifndef FortranCInterface_HEADER_INCLUDED
#define FortranCInterface_HEADER_INCLUDED

/* Mangling for Fortran global symbols without underscores. */
#define FortranCInterface_GLOBAL(name,NAME) name##_

/* Mangling for Fortran global symbols with underscores. */
#define FortranCInterface_GLOBAL_(name,NAME) name##_

/* Mangling for Fortran module symbols without underscores. */
#define FortranCInterface_MODULE(mod_name,name, mod_NAME,NAME) __##mod_name##_MOD_##name

/* Mangling for Fortran module symbols with underscores. */
#define FortranCInterface_MODULE_(mod_name,name, mod_NAME,NAME) __##mod_name##_MOD_##name

/*--------------------------------------------------------------------------*/
/* Mangle some symbols automatically.                                       */
#define write_restart FortranCInterface_GLOBAL_(write_restart, WRITE_RESTART)
#define write_error FortranCInterface_GLOBAL_(write_error, WRITE_ERROR)
#define write_displ FortranCInterface_GLOBAL_(write_displ, WRITE_DISPL)
#define write_field FortranCInterface_GLOBAL_(write_field, WRITE_FIELD)
#define write_hessian FortranCInterface_GLOBAL_(write_hessian, WRITE_HESSIAN)

#endif
