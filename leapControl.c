/* leap_control.c
deciding wether leap year of not

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGCMuSo v7.0.
Copyright 2025, D. Hidy [dori.hidy@gmail.com]
Hungarian Academy of Sciences, Hungary
See the website of Biome-BGCMuSo at http://nimbus.elte.hu/bbgc/ for documentation, model executable and example input files.
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include "ini.h"
#include "bgc_struct.h"
#include "bgc_func.h"
#include "bgc_constants.h"

int leapControl(int year, int* enddays, int* mondays, int* leap)
{
	int enddays_normal[nMONTHS_OF_YEAR]  = {30,58,89,119,150,180,211,242,272,303,333,364};
	int enddays_leap[nMONTHS_OF_YEAR]    = {30,59,90,120,151,181,212,243,273,304,334,364};
	int mondays_normal[nMONTHS_OF_YEAR] = {31,28,31,30, 31, 30, 31, 31, 30, 31, 30, 31};
	int mondays_leap[nMONTHS_OF_YEAR]   = {31,29,31,30, 31, 30, 31, 31, 30, 31, 30, 30};
	
	double leap_control;
	int i;
	int errorCode=0;

	/* leap year or normal year?*/
	leap_control = (double) year /4.;
	if (leap_control - (int) leap_control == 0) 
	{
		for (i=0; i<nMONTHS_OF_YEAR; i++) 
		{
			enddays[i] = enddays_leap[i];
			mondays[i] = mondays_leap[i];
		}
		*leap = 1;
	}				
	else
	{
		for (i=0; i<nMONTHS_OF_YEAR; i++) 
		{
			enddays[i] = enddays_normal[i];
			mondays[i] = mondays_normal[i];
		}
		*leap = 0;
	}

	return (errorCode);

}