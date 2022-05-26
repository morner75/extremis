#include <R.h>
#include <Rinternals.h>
#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Fortran calls */
extern void F77_NAME(ptdensitybetaupmh)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);


static const R_FortranMethodDef FortranEntries[] = {
    {"ptdensitybetaupmh",      (DL_FUNC) &F77_NAME(ptdensitybetaupmh),      35},
    {NULL, NULL, 0}
};

/* .Call calls */
extern SEXP BERN(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"BERN",   (DL_FUNC) &BERN,   11},
    {NULL, NULL, 0}
};

#ifdef _WIN32
# include <fcntl.h>
#endif

void R_init_extremis(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
  
#ifdef _WIN32
  /* gfortran initialization sets these to _O_BINARY */
  setmode(1, _O_TEXT); /* stdout */
  setmode(2, _O_TEXT); /* stderr */
#endif
}

