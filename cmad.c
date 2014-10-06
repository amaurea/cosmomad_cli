#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cosmo_mad.h>

#define MAXFIELD 0x100
#define OMEGA_M 0.3086
#define OMEGA_L 0.6914
#define OMEGA_B 0.0483
#define W_PARAM  -1.0
#define WA_PARAM  0.0
#define HUBBLE_h 0.6777
#define TCMB 2.73

typedef double (*fptr)(Csm_params*,double);
void add_field(char ** fields, int * units, fptr * funcs, int * n, char * name, int unit, fptr fun)
{
	fields[*n] = name;
	funcs[*n]  = fun;
	units[*n]  = unit;
	++*n;
	if(*n >= MAXFIELD) *n = MAXFIELD-1;
}
double afun(Csm_params* p, double a) { return a; }
double rt(Csm_params* p, double a) { return csm_cosmic_time(p, 1.0)-csm_cosmic_time(p,a); }
double zfun(Csm_params* p, double a) { return 1/a-1; }

int main(int argc, char ** argv)
{
	double
		omega_m = OMEGA_M,
		omega_l = OMEGA_L,
		omega_b = OMEGA_B,
		w       = W_PARAM,
		wa      = WA_PARAM,
		h       = HUBBLE_h,
		T       = TCMB;
	double a0 = 0, a1 = 1;
	int na = 1000;
	char * fields[MAXFIELD];
	int  units[MAXFIELD];
	fptr funcs[MAXFIELD];
	int nfield = 0;
	int verbosity = 0;
	int physical = 1;
	char * fmt = " %15.7e";
	char ** i;
	int ia, j;
	for(i = argv+1; *i; i++)
		if(!strcmp(*i, "-n")) na = atoi(*++i);
		else if(!strcmp(*i, "-a0"))  a0 = atof(*++i);
		else if(!strcmp(*i, "-a1"))  a1 = atof(*++i);
		else if(!strcmp(*i, "-T"))   T  = atof(*++i);
		else if(!strcmp(*i, "-h"))   h  = atof(*++i);
		else if(!strcmp(*i, "-w"))   w  = atof(*++i);
		else if(!strcmp(*i, "-wa"))  wa = atof(*++i);
		else if(!strcmp(*i, "-Ob"))  omega_b = atof(*++i);
		else if(!strcmp(*i, "-Ol"))  omega_l = atof(*++i);
		else if(!strcmp(*i, "-Om"))  omega_m = atof(*++i);
		else if(!strcmp(*i, "-v"))   verbosity = 1;
		else if(!strcmp(*i, "-q"))   verbosity = 0;
		else if(!strcmp(*i, "-phys"))physical  = 1;
		else if(!strcmp(*i, "-nophys"))physical= 0;
		else if(!strcmp(*i, "a"))    add_field(fields, units, funcs, &nfield, *i, 0, afun);
		else if(!strcmp(*i, "t"))    add_field(fields, units, funcs, &nfield, *i, 1, csm_cosmic_time);
		else if(!strcmp(*i, "H"))    add_field(fields, units, funcs, &nfield, *i,-1, csm_hubble);
		else if(!strcmp(*i, "dA"))   add_field(fields, units, funcs, &nfield, *i, 1, csm_angular_diameter_distance);
		else if(!strcmp(*i, "dL"))   add_field(fields, units, funcs, &nfield, *i, 1, csm_luminosity_distance);
		else if(!strcmp(*i, "D"))    add_field(fields, units, funcs, &nfield, *i, 0, csm_growth_factor);
		else if(!strcmp(*i, "f"))    add_field(fields, units, funcs, &nfield, *i, 0, csm_f_growth);
		else if(!strcmp(*i, "bao"))  add_field(fields, units, funcs, &nfield, *i, 0, csm_theta_BAO);
		else if(!strcmp(*i, "chi"))  add_field(fields, units, funcs, &nfield, *i, 1, csm_radial_comoving_distance);
		else if(!strcmp(*i, "chi_p"))add_field(fields, units, funcs, &nfield, *i, 1, csm_particle_horizon);
		else if(!strcmp(*i, "r"))    add_field(fields, units, funcs, &nfield, *i, 1, csm_curvature_comoving_distance);
		else if(!strcmp(*i, "rt"))   add_field(fields, units, funcs, &nfield, *i, 1, rt);
		else if(!strcmp(*i, "z"))    add_field(fields, units, funcs, &nfield, *i, 0, zfun);
		else {
			fprintf(stderr, "Unrecognized field '%s'!\n", *i);
			exit(1);
		}
	if(nfield == 0)
	{
		fprintf(stderr, "Specify at least one output field!\n");
		exit(1);
	}
	csm_set_verbosity(verbosity);
	Csm_params * params = csm_params_new();
	csm_background_set(params, omega_m, omega_l, omega_b, w, wa, h, T);
	for(ia=0; ia<na; ia++)
	{
		double a = a0+(a1-a0)*ia/(na-1);
		for(j=0;j<nfield;j++)
		{
			double v = funcs[j](params, a);
			if(physical) v *= pow(h, -units[j]);
			printf(fmt, v);
		}
		printf("\n");
	}
	return 0;
}
