/* DIRect-transp.f -- translated by f2c (version 20050501).

   f2c output hand-cleaned by SGJ (August 2007).
*/

#include <math.h>
#include "direct-internal.h"

/* Common Block Declarations */

/* Table of constant values */

/* +-----------------------------------------------------------------------+ */
/* | Program       : Direct.f                                              | */
/* | Last modified : 07-16-2001                                            | */
/* | Written by    : Joerg Gablonsky (jmgablon@unity.ncsu.edu)             | */
/* |                 North Carolina State University                       | */
/* |                 Dept. of Mathematics                                  | */
/* +-----------------------------------------------------------------------+ */
/* Subroutine */ void direct_direct_(fp fcn, doublereal *x, integer *n, doublereal *eps, doublereal epsabs, integer *maxf, integer *maxt, double starttime, double maxtime, int *force_stop, doublereal *minf, doublereal *l,
	doublereal *u, integer *algmethod, integer *ierror, FILE *logfile,
	doublereal *fglobal, doublereal *fglper, doublereal *volper,
	doublereal *sigmaper, void *fcn_data)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* changed by SGJ to be dynamically allocated ... would be
       even better to use realloc, below, to grow these as needed */
    integer MAXFUNC = *maxf <= 0 ? 101000 : (*maxf + 1000 + *maxf / 2);
    integer MAXDEEP = *maxt <= 0 ? MAXFUNC/5: *maxt + 1000;
    const integer MAXDIV = 5000;

    /* Local variables */
    integer increase;
    doublereal *c__ = 0	/* was [90000][64] */, *f = 0	/*
	    was [90000][2] */;
    integer i__, j, *s = 0	/* was [3000][2] */, t;
    doublereal *w = 0;
    doublereal divfactor;
    integer ifeasiblef, iepschange, actmaxdeep;
    integer actdeep_div__, iinfesiblef;
    integer pos1, newtosample;
    integer ifree, help;
    doublereal *oldl = 0, fmax;
    integer maxi;
    doublereal kmax, *oldu = 0;
    integer oops, *list2 = 0	/* was [64][2] */, cheat;
    doublereal delta;
    integer mdeep, *point = 0, start;
    integer *anchor = 0, *length = 0	/* was [90000][64] */, *arrayi = 0;
    doublereal *levels = 0, *thirds = 0;
    doublereal epsfix;
    integer oldpos, minpos, maxpos, tstart, actdeep, ifreeold, oldmaxf;
    integer numfunc, version;
    integer jones;

#define MY_ALLOC(p, t, n) p = (t *) malloc(sizeof(t) * (n)); \
                          if (!(p)) { *ierror = -100; goto cleanup; }

    MY_ALLOC(c__, doublereal, MAXFUNC * (*n));
    MY_ALLOC(length, integer, MAXFUNC * (*n));
    MY_ALLOC(f, doublereal, MAXFUNC * 2);
    MY_ALLOC(point, integer, MAXFUNC);
    if (*maxf <= 0) *maxf = MAXFUNC - 1000;

    MY_ALLOC(s, integer, MAXDIV * 2);

    MY_ALLOC(anchor, integer, MAXDEEP + 2);
    MY_ALLOC(levels, doublereal, MAXDEEP + 1);
    MY_ALLOC(thirds, doublereal, MAXDEEP + 1);
    if (*maxt <= 0) *maxt = MAXDEEP;

    MY_ALLOC(w, doublereal, (*n));
    MY_ALLOC(oldl, doublereal, (*n));
    MY_ALLOC(oldu, doublereal, (*n));
    MY_ALLOC(list2, integer, (*n) * 2);
    MY_ALLOC(arrayi, integer, (*n));

    /* Parameter adjustments */
    --u;
    --l;
    --x;

    /* Function Body */
    jones = *algmethod;

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	oldu[i__ - 1] = u[i__];
	oldl[i__ - 1] = l[i__];
    }

    version = 204;
    cheat = 0;
    kmax = 1e10;
    mdeep = MAXDEEP;

    direct_dirheader_(logfile, &version, &x[1], n, eps, maxf, maxt, &l[1], &u[1],
	    algmethod, &MAXFUNC, &MAXDEEP, fglobal, fglper, ierror, &epsfix, &
		      iepschange, volper, sigmaper);
    if (*ierror < 0) {
	goto cleanup;
    }
    if (*fglobal == 0.) {
	divfactor = 1.;
    } else {
	divfactor = fabs(*fglobal);
    }
    oldmaxf = *maxf;
    increase = 0;
    direct_dirinitlist_(anchor, &ifree, point, f, &MAXFUNC, &MAXDEEP);
    direct_dirpreprc_(&u[1], &l[1], n, &l[1], &u[1], &oops);
    if (oops > 0) {
	if (logfile)
	     fprintf(logfile,"WARNING: Initialization in DIRpreprc failed.\n");
	*ierror = -3;
	goto cleanup;
    }
    tstart = 2;
    direct_dirinit_(f, fcn, c__, length, &actdeep, point, anchor, &ifree,
	    logfile, arrayi, &maxi, list2, w, &x[1], &l[1], &u[1],
	    minf, &minpos, thirds, levels, &MAXFUNC, &MAXDEEP, n, n, &
	    fmax, &ifeasiblef, &iinfesiblef, ierror, fcn_data, jones,
		    starttime, maxtime, force_stop);
    if (*ierror < 0) {
	if (*ierror == -4) {
	    if (logfile)
		 fprintf(logfile, "WARNING: Error occurred in routine DIRsamplepoints.\n");
	    goto cleanup;
	}
	if (*ierror == -5) {
	    if (logfile)
		 fprintf(logfile, "WARNING: Error occurred in routine DIRsamplef..\n");
	    goto cleanup;
	}
	if (*ierror == -102) goto L100;
    }
    else if (*ierror == DIRECT_MAXTIME_EXCEEDED) goto L100;
    numfunc = maxi + 1 + maxi;
    actmaxdeep = 1;
    oldpos = 0;
    tstart = 2;
    if (ifeasiblef > 0) {
	 if (logfile)
	      fprintf(logfile, "No feasible point found in %d iterations "
		      "and %d function evaluations.\n", tstart-1, numfunc);
    } else {
	 if (logfile)
	      fprintf(logfile, "%d, %d, %g, %g\n",
		      tstart-1, numfunc, *minf, fmax);
    }

/* +-----------------------------------------------------------------------+ */
/* | Main loop!                                                            | */
/* +-----------------------------------------------------------------------+ */
    i__1 = *maxt;
    t = tstart - 1;
    for (t = tstart; t <= i__1; ++t) {

	actdeep = actmaxdeep;
	direct_dirchoose_(anchor, s, &MAXDEEP, f, minf, *eps, epsabs, levels, &maxpos, length,
		&MAXFUNC, &MAXDEEP, &MAXDIV, n, logfile, &cheat, &
		kmax, &ifeasiblef, jones);

	if (*algmethod == 0) {
	     direct_dirdoubleinsert_(anchor, s, &maxpos, point, f, &MAXDEEP, &MAXFUNC,
		     &MAXDIV, ierror);
	    if (*ierror == -6) {
		if (logfile)
		     fprintf(logfile,
"WARNING: Capacity of array S in DIRDoubleInsert reached. Increase maxdiv.\n"
"This means that there are a lot of hyperrectangles with the same function\n"
"value at the center. We suggest to use our modification instead (Jones = 1)\n"
			  );
		goto cleanup;
	    }
	}

	oldpos = minpos;
	newtosample = 0;
	i__2 = maxpos;
	for (j = 1; j <= i__2; ++j) {
	    actdeep = s[j + MAXDIV-1];
	    if (s[j - 1] > 0) {
		actdeep_div__ = direct_dirgetmaxdeep_(&s[j - 1], length, &MAXFUNC, n);
		delta = thirds[actdeep_div__ + 1];
		actdeep = s[j + MAXDIV-1];
		if (actdeep + 1 >= mdeep) {
		    if (logfile)
			 fprintf(logfile, "WARNING: Maximum number of levels reached. Increase maxdeep.\n");
		    *ierror = -6;
		    goto L100;
		}
		actmaxdeep = MAX(actdeep,actmaxdeep);
		help = s[j - 1];
		if (! (anchor[actdeep + 1] == help)) {
		    pos1 = anchor[actdeep + 1];
		    while(! (point[pos1 - 1] == help)) {
			pos1 = point[pos1 - 1];
		    }
		    point[pos1 - 1] = point[help - 1];
		} else {
		    anchor[actdeep + 1] = point[help - 1];
		}
		if (actdeep < 0) {
		    actdeep = (integer) f[(help << 1) - 2];
		}
		direct_dirget_i__(length, &help, arrayi, &maxi, n, &MAXFUNC);
		direct_dirsamplepoints_(c__, arrayi, &delta, &help, &start, length,
			logfile, f, &ifree, &maxi, point, &x[
			1], &l[1], minf, &minpos, &u[1], n, &MAXFUNC, &
			MAXDEEP, &oops);
		if (oops > 0) {
		    if (logfile)
			 fprintf(logfile, "WARNING: Error occurred in routine DIRsamplepoints.\n");
		    *ierror = -4;
		    goto cleanup;
		}
		newtosample += maxi;
		direct_dirsamplef_(c__, arrayi, &delta, &help, &start, length,
			    logfile, f, &ifree, &maxi, point, fcn, &x[
			1], &l[1], minf, &minpos, &u[1], n, &MAXFUNC, &
			MAXDEEP, &oops, &fmax, &ifeasiblef, &iinfesiblef,
				   fcn_data, force_stop);
		if (force_stop && *force_stop) {
		     *ierror = -102;
		     goto L100;
		}
		if (0) {
		     *ierror = DIRECT_MAXTIME_EXCEEDED;
		     goto L100;
		}
		if (oops > 0) {
		    if (logfile)
			 fprintf(logfile, "WARNING: Error occurred in routine DIRsamplef.\n");
		    *ierror = -5;
		    goto cleanup;
		}
		direct_dirdivide_(&start, &actdeep_div__, length, point, arrayi, &
			help, list2, w, &maxi, f, &MAXFUNC, &MAXDEEP, n);
		direct_dirinsertlist_(&start, anchor, point, f, &maxi, length, &
			MAXFUNC, &MAXDEEP, n, &help, jones);
		numfunc = numfunc + maxi + maxi;
	    }
	}
	if (oldpos < minpos) {
	    if (logfile)
		 fprintf(logfile, "%d, %d, %g, %g\n",
			 t, numfunc, *minf, fmax);
	}
	if (ifeasiblef > 0) {
	    if (logfile)
		 fprintf(logfile, "No feasible point found in %d iterations "
			 "and %d function evaluations\n", t, numfunc);
	}

/* +-----------------------------------------------------------------------+ */
/* |                       Termination Checks                              | */
/* +-----------------------------------------------------------------------+ */
    // ВОЗВРАЩЕНА ОРИГИНАЛЬНАЯ ЛОГИКА ГАБЛОНСКОГО
	*ierror = jones;
	jones = 0;
	actdeep_div__ = direct_dirgetlevel_(&minpos, length, &MAXFUNC, n, jones);
	jones = *ierror;
	delta = thirds[actdeep_div__] * 100;
	if (delta <= *volper) {
	    *ierror = 4;
	    if (logfile) {
            if (jones == 1) {
                fprintf(logfile, "DIRECT stopped: Max side length of S_min is "
                         "%g%% <= %g%% of the original side length.\n", delta, *volper);
            } else {
                fprintf(logfile, "DIRECT stopped: Volume of S_min is "
                         "%g%% <= %g%% of the original volume.\n", delta, *volper);
            }
        }
	    goto L100;
	}

	actdeep_div__ = direct_dirgetlevel_(&minpos, length, &MAXFUNC, n, jones);
	delta = levels[actdeep_div__];
	if (delta <= *sigmaper) {
	    *ierror = 5;
	    if (logfile)
		 fprintf(logfile, "DIRECT stopped: Measure of S_min "
			 "= %g < %g.\n", delta, *sigmaper);
	    goto L100;
	}
	if ((*minf - *fglobal) * 100 / divfactor <= *fglper) {
	    *ierror = 3;
	    if (logfile)
		 fprintf(logfile, "DIRECT stopped: minf within fglper of global minimum.\n");
	    goto L100;
	}
	if (iinfesiblef > 0) {
	     direct_dirreplaceinf_(&ifree, &ifreeold, f, c__, thirds, length, anchor,
		    point, &u[1], &l[1], &MAXFUNC, &MAXDEEP, n, n,
		    logfile, &fmax, jones);
	}
	ifreeold = ifree;
	if (iepschange == 1) {
	    d__1 = fabs(*minf) * 1e-4;
	    *eps = MAX(d__1,epsfix);
	}
	if (increase == 1) {
	    *maxf = numfunc + oldmaxf;
	    if (ifeasiblef == 0) {
		if (logfile)
		     fprintf(logfile, "DIRECT found a feasible point.  The "
			     "adjusted budget is now set to %d.\n", *maxf);
		increase = 0;
	    }
	}
	if (numfunc > *maxf) {
	    if (ifeasiblef == 0) {
		*ierror = 1;
		if (logfile)
		     fprintf(logfile, "DIRECT stopped: numfunc >= maxf.\n");
		goto L100;
	    } else {
		increase = 1;
		if (logfile)
                     fprintf(logfile,
"DIRECT could not find a feasible point after %d function evaluations.\n"
"DIRECT continues until a feasible point is found.\n", numfunc);
		*maxf = numfunc + oldmaxf;
	    }
	}
    }
    *ierror = 2;
    if (logfile)
	 fprintf(logfile, "DIRECT stopped: maxT iterations.\n");

L100:
    *maxt = t;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = c__[i__ + minpos * i__1 - i__1-1] * l[i__] + l[i__] * u[i__];
	u[i__] = oldu[i__ - 1];
	l[i__] = oldl[i__ - 1];
    }
    *maxf = numfunc;
    direct_dirsummary_(logfile, &x[1], &l[1], &u[1], n, minf, fglobal, &numfunc,
	    ierror);

 cleanup:
#define MY_FREE(p) if (p) free(p)
    MY_FREE(c__);
    MY_FREE(f);
    MY_FREE(s);
    MY_FREE(w);
    MY_FREE(oldl);
    MY_FREE(oldu);
    MY_FREE(list2);
    MY_FREE(point);
    MY_FREE(anchor);
    MY_FREE(length);
    MY_FREE(arrayi);
    MY_FREE(levels);
    MY_FREE(thirds);
} /* direct_ */
