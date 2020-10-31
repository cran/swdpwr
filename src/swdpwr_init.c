#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(computeparameter)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(computeparametergeelog)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(computeparametergeelogit)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(computeparameterlog)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(computeparameterlogit)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(continuouspowergeenotimewrapper)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(continuouspowergeetimewrapper)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(herzo)(void *, void *, void *);
extern void F77_NAME(legendrehandle)(void *, void *, void *, void *, void *);
extern void F77_NAME(linearpowergeenotimewrapper)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(linearpowergeewrapper)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(linearpowernotimewrapper)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(linearpowertimewrapper)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(logitpowergeenotimewrapper)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(logitpowergeewrapper)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(logitpowernotimewrapper)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(logitpowertimewrapper)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(logpowergeenotimewrapper)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(logpowergeewrapper)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(logpowernotimewrapper)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(logpowertimewrapper)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"computeparameter",                (DL_FUNC) &F77_NAME(computeparameter),                 8},
    {"computeparametergeelog",          (DL_FUNC) &F77_NAME(computeparametergeelog),           8},
    {"computeparametergeelogit",        (DL_FUNC) &F77_NAME(computeparametergeelogit),         8},
    {"computeparameterlog",             (DL_FUNC) &F77_NAME(computeparameterlog),              8},
    {"computeparameterlogit",           (DL_FUNC) &F77_NAME(computeparameterlogit),           12},
    {"continuouspowergeenotimewrapper", (DL_FUNC) &F77_NAME(continuouspowergeenotimewrapper), 13},
    {"continuouspowergeetimewrapper",   (DL_FUNC) &F77_NAME(continuouspowergeetimewrapper),   13},
    {"herzo",                           (DL_FUNC) &F77_NAME(herzo),                            3},
    {"legendrehandle",                  (DL_FUNC) &F77_NAME(legendrehandle),                   5},
    {"linearpowergeenotimewrapper",     (DL_FUNC) &F77_NAME(linearpowergeenotimewrapper),     12},
    {"linearpowergeewrapper",           (DL_FUNC) &F77_NAME(linearpowergeewrapper),           12},
    {"linearpowernotimewrapper",        (DL_FUNC) &F77_NAME(linearpowernotimewrapper),        14},
    {"linearpowertimewrapper",          (DL_FUNC) &F77_NAME(linearpowertimewrapper),          17},
    {"logitpowergeenotimewrapper",      (DL_FUNC) &F77_NAME(logitpowergeenotimewrapper),      12},
    {"logitpowergeewrapper",            (DL_FUNC) &F77_NAME(logitpowergeewrapper),            12},
    {"logitpowernotimewrapper",         (DL_FUNC) &F77_NAME(logitpowernotimewrapper),         12},
    {"logitpowertimewrapper",           (DL_FUNC) &F77_NAME(logitpowertimewrapper),           13},
    {"logpowergeenotimewrapper",        (DL_FUNC) &F77_NAME(logpowergeenotimewrapper),        12},
    {"logpowergeewrapper",              (DL_FUNC) &F77_NAME(logpowergeewrapper),              12},
    {"logpowernotimewrapper",           (DL_FUNC) &F77_NAME(logpowernotimewrapper),           14},
    {"logpowertimewrapper",             (DL_FUNC) &F77_NAME(logpowertimewrapper),             17},
    {NULL, NULL, 0}
};

void R_init_swdpwr(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
