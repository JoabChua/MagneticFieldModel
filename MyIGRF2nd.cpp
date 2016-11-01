#include <iostream>
#include <math.h>
#include <iomanip>
#define PI 3.141592654
using namespace std;

void getshc(double gh1[196], double gh2[196]);
double julday(int month, int day, int year);
void extrapsh(double dyear, double gh1[196], double gh2[196], double gha[196], double ghb[196]);
//double interpsh(double dyear, double gh1[195], double gh2[195]);
void shval3(double dlat, double dlon, int igdgc, double alt, double gha[196], double ghb[196], double ans[6], double dyear);
void dihf (double ans[6]);

int main ()
{
double gha[196] = {0}, ghb[196] = {0}, ans[6] = {0};
double yy,mm,dd,alt,dlat,dlon,yy1;
double gh1[196] = {0}, gh2[196] = {0}; double dyear=0,dyearstart=0,dyearend=0,stepsize=0;
int igdgc=0;
int dtype=0,atype=0;

// get info to choose dyear or yymmdd
cout << "Enter date" << endl << "1) Decimal Year" << endl << "2) Year, Month, Day" << endl << "Choice: ";
cin >> dtype;

// get info for computing range or single date
cout << "Would you like output for a single date or for a range of dates?" << endl;
cout << "1) Single date" << endl << "2) Range of date" << endl << "Choice: ";
cin >> atype;

if (atype == 1) // receive single date
{
if (dtype ==2) // get yy, mm, dd
{
  // get info for which year
cout << "Enter Year (2015 - 2020): ";
cin >> yy;
if (yy < 2015 || yy > 2020)
   {
	cout << "Error: unrecognized year" << endl;
   }
	while ((yy < 2015 || yy > 2020))
	   {
  	   cout << "Enter Year (2015 - 2020): ";
   	   cin >> yy;
   	   }

 // get info for which month
cout << "Enter Month (1 - 12): ";
cin >> mm;
if (mm < 1 || mm > 12)
  {
	cout << "Error: unrecognized month" << endl;
  }
   while ((mm < 1 || mm > 12))
   {
   cout << "Enter Month (1 - 12): ";
   cin >> mm;
   }
// get info on which day
cout << "Enter Day (1 - 31): ";
cin >> dd;
if (dd < 1 || dd > 31)
	{
	cout << "Error: unrecognized day" << endl;
	}
	while ((dd < 1 || dd > 31))
 		{
  		cout << "Enter Day (1 - 31): ";
   		cin >> dd;
   		}

//get date in decimal year
dyear = julday(mm, dd, yy);
}
else // get decimal year
{
	cout << "Enter decimal date (2015.00 to 2020): ";
	cin >> dyear;
	if (dyear < 2015.00 || dyear > 2020)
		{
		cout << "Error: unrecognized date here" << endl;
		}
	while ((dyear < 2015 || dyear > 2020))
	   {
  	   cout << "Enter Year (2015 - 2020): ";
   	   cin >> dyear;
   	   }
}
}
else // receive range of dates
{
	if (dtype ==2) // get start yy, mm, dd
	{
	  // get info for the start year
	cout << "Enter Start Date" << endl << "Year (2015 - 2020): ";
	cin >> yy;
	if (yy < 2015 || yy > 2020)
	   {
		cout << "Error: unrecognized year" << endl;
	   }
		while ((yy < 2015 || yy > 2020))
		   {
	  	   cout << "Year (2015 - 2020): ";
	   	   cin >> yy;
	   	   }

	 // get info for which month
	cout << "Month (1 - 12): ";
	cin >> mm;
	if (mm < 1 || mm > 12)
	  {
		cout << "Error: unrecognized month" << endl;
	  }
	   while ((mm < 1 || mm > 12))
	   {
	   cout << "Month (1 - 12): ";
	   cin >> mm;
	   }
	// get info on which day
	cout << "Day (1 - 31): ";
	cin >> dd;
	if (dd < 1 || dd > 31)
		{
		cout << "Error: unrecognized day" << endl;
		}
		while ((dd < 1 || dd > 31))
	 		{
	  		cout << "Day (1 - 31): ";
	   		cin >> dd;
	   		}
	//get start date in decimal year
	dyearstart = julday(mm, dd, yy);

	// get end yy, mm, dd
	cout << "Enter End Date" << endl << "Year (2015 - 2020): ";
	cin >> yy1;
	if (yy1 < yy || yy1 > 2020)
	   {
		cout << "Error: unrecognized year" << endl;
	   }
		while ((yy1 < yy || yy > 2020))
		   {
	  	   cout << "Year (2015 - 2020): ";
	   	   cin >> yy1;
	   	   }

	 // get info for which month
	cout << "Month (1 - 12): ";
	cin >> mm;
	if (mm < 1 || mm > 12)
	  {
		cout << "Error: unrecognized month" << endl;
	  }
	   while ((mm < 1 || mm > 12))
	   {
	   cout << "Month (1 - 12): ";
	   cin >> mm;
	   }
	// get info on which day
	cout << "Day (1 - 31): ";
	cin >> dd;
	if (dd < 1 || dd > 31)
		{
		cout << "Error: unrecognized day" << endl;
		}
		while ((dd < 1 || dd > 31))
	 		{
	  		cout << "Day (1 - 31): ";
	   		cin >> dd;
	}
	//get end date in decimal year
	dyearend = julday(mm, dd, yy1);
	}
	else
	{
		// get start date in decimal year
		cout << "Enter decimal start date (2015.00 to 2020): ";
		cin >> dyearstart;
		if (dyearstart < 2015.00 || dyearstart > 2020)
			{
			cout << "Error: unrecognized date" << endl;
			}
			while ((dyearstart < 2015.00 || dyearstart > 2020))
		 		{
		  		cout << "Enter decimal start date (2015.00 to 2020): ";
		   		cin >> dyearstart;
		   		}
		// get end date in decimal year
		cout << "Enter decimal end date (2015.00 to 2020): ";
		cin >> dyearend;
		if (dyearend < dyearstart || dyearend > 2020)
		{
		cout << "Error: unrecognized date" << endl;
		}
		while ((dyearend < dyearstart || dyearend > 2020))
		{
		cout << "Enter decimal end date (2015.00 to 2020): ";
		cin >> dyearend;
		}
	}
	double nyear = dyearend-dyearstart;
	cout << "Enter Step Size in 2 decimal places (0.01 to " << fixed << setprecision(2) << nyear << "): ";
	cin >> stepsize;
	if (stepsize <=0 || stepsize > (dyearend-dyearstart))
	{
		cout << "Error: Unrecognized Step Size" << endl;
	}
	while (stepsize <=0 || stepsize > (dyearend-dyearstart))
	{
		cout << "Enter Step Size in 2 decimal places (0.01 to " << fixed << setprecision(2) << nyear << "): ";
		cin >> stepsize;
	}
}


//// get info for altitude
//cout << "Enter Coordinate Preferences:" << endl<< "1) Geodetic (WGS84 latitude and altitude above mean sea level)"
//		<< endl << "2) Geocentric (spherical, altitude relative to Earth's center)" << endl << "Choice: ";
//cin >> igdgc;
//if (igdgc < 0 || igdgc > 2)
//{
//	cout << "Error: unrecognized value" << endl;
// }
// while ((igdgc<0)||(igdgc>2))
// 	   {
//	 cout << "Enter" << endl<< "1) Geodetic (WGS84 latitude and altitude above mean sea level)"
//			 << endl << "2) Geocentric (spherical, altitude relative to Earth's center):  " << endl;
//	   cin >> igdgc;
//	   }
igdgc = 1;
// get info for altitude
if (igdgc == 1)
{
cout << "Enter geodetic altitude above mean sea level in km (-1.00 to 600.00): ";
cin >> alt;
if (alt < -1 || alt > 600)
{
	cout << "Error: unrecognized altitude" << endl;
 }
 while ((alt < -1)||(alt > 600))
 	   {
	   cout << "Enter geodetic altitude above mean sea level in km (-1.00 to 600.00): ";
	   cin >> alt;
	   }
}
else // igdgc == 2
{
	cout << "Enter geocentric altitude in km (6370.20 to 6971.20): ";
	cin >> alt;
	if (alt < 6370.20 || alt > 6971.20)
	{
		cout << "Error: unrecognized altitude" << endl;
	 }
	 while ((alt < 6370.20)||(alt > 6971.20))
	{
	cout << "Enter geocentric altitude in km (6370.20 to 6971.20): ";
	cin >> alt;
	}
	alt = alt - 6371.2;
}

 // get info of longitude
cout << "Enter the decimal degree longitude (-180 to 180) (- for Western hemisphere): ";
cin >> dlon;
if (dlon < -180 || dlon > 180)
	{
	cout << "Error: unrecognized longitude" << endl;
	}
   while ((dlon < -180)||(dlon > 180))
 		   {
 		   cout << "Enter the decimal degree longitude (-180 to 180) (- for Western hemisphere): ";
		   cin >> dlon;
 		   }

  // get info for latitude
cout << "Enter the decimal degree latitude (-90 to 90) (- for Southern hemisphere): ";
cin >> dlat;
if (dlat < -90 || dlat > 90)
  {
	cout << "Error: unrecognized latitude" << endl;
   }
   while ((dlat < -90)||(dlat > 90))
   	   {
	   cout << "Enter the decimal degree latitude (-90 to 90) (- for Southern hemisphere): ";
 	   cin >> dlat;
 	   }

// header for the output values
cout << "  Date      X-Coord     Y-Coord     Z-Coord     H-Intensity    T-Intensity    Inclination    Declination" << endl;
cout << "Dec  Year     (nT)        (nT)        (nT)         (nT)           (nT)           (deg)          (deg)   " << endl;

// load spherical harmoinc coefficients into array gh
getshc(gh1, gh2);

if (atype == 2) // calculation for a range of dates
{
	dyear = dyearstart;
while (dyear < dyearend)
	{
	//interpsh(dyear, gh1, gh2);
	//else
	//{
	extrapsh(dyear, gh1, gh2, gha, ghb);
	//}
	// get x y z coordinates
	shval3(dlat, dlon, igdgc, alt, gha, ghb, ans, dyear);
	// get total and horizontal intensities, declination and inclination
	dihf(ans);
	dyear = dyear + stepsize;
	}
	dyear = dyearend;
	//interpsh(dyear, gh1, gh2);
	//else
	//{
	extrapsh(dyear, gh1, gh2, gha, ghb);
	//}
	// get x y z coordinates
	shval3(dlat, dlon, igdgc, alt, gha, ghb, ans,dyear);
	// get total and horizontal intensities, declination and inclination
	dihf(ans);
}
else // for a single date
{
//interpsh(dyear, gh1, gh2);
//else
//{
	extrapsh(dyear, gh1, gh2, gha, ghb);
//}
// get x y z coordinates
shval3(dlat, dlon, igdgc, alt, gha, ghb, ans, dyear);
// get total and horizontal intensities, declination and inclination
dihf(ans);
}
return 0;
}

/****************************************************************************/
/*                                                                          */
/*                           Subroutine getshc                              */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Reads spherical harmonic coefficients from the specified             */
/*     model into an array.                                                 */
/*                                                                          */
/*     Input:                                                               */
/*           stream     - Logical unit number                               */
/*           iflag      - Flag for SV equal to, or not equal to 0           */
/*                        for designated read statements                    */
/*           strec      - Starting record number to read from model         */
/*           nmax_of_gh - Maximum degree and order of model                 */
/*                                                                          */
/*     Output:                                                              */
/*           gh1 or 2   - Schmidt quasi-normal internal spherical           */
/*                        harmonic coefficients                             */
/*                                                                          */
/****************************************************************************/

void getshc(double gh1[196], double gh2[196])
{
double gg[196][2] = {
		{0,0,},
		{-29442.0,10.3,},
		{-1501.0,18.1,},
		{4797.1,-26.6,},
		{-2445.1,-8.7,},
		{3012.9,-3.3,},
		{-2845.6,-27.4,},
		{1676.7,2.1,},
		{-641.9,-14.1,},
		{1350.7,3.4,},
		{-2352.3,-5.5,},
		{-115.3,8.2,},
		{1225.6,-0.7,},
		{244.9,-0.4,},
		{582.0,-10.1,},
		{-538.4,1.8,},
		{907.6,-0.7,},
		{813.7,0.2,},
		{283.3,-1.3,},
		{120.4,-9.1,},
		{-188.7,5.3,},
		{-334.9,4.1,},
		{180.9,2.9,},
		{70.4,-4.3,},
		{-329.5,-5.2,},
		{-232.6,-0.2,},
		{360.1,0.5,},
		{47.3,0.6,},
		{192.4,-1.3,},
		{197.0,1.7,},
		{-140.9,-0.1,},
		{-119.3,-1.2,},
		{-157.5,1.4,},
		{16.0,3.4,},
		{4.1,3.9,},
		{100.2,0.0,},
		{70.0,-0.3,},
		{67.7,-0.1,},
		{-20.8,0.0,},
		{72.7,-0.7,},
		{33.2,-2.1,},
		{-129.9,2.1,},
		{58.9,-0.7,},
		{-28.9,-1.2,},
		{-66.7,0.2,},
		{13.2,0.3,},
		{7.3,0.9,},
		{-70.9,1.6,},
		{62.6,1.0,},
		{81.6,0.3,},
		{-76.1,-0.2,},
		{-54.1,0.8,},
		{-6.8,-0.5,},
		{-19.5,0.4,},
		{51.8,1.3,},
		{5.7,-0.2,},
		{15.0,0.1,},
		{24.4,-0.3,},
		{9.4,-0.6,},
		{3.4,-0.6,},
		{-2.8,-0.8,},
		{-27.4,0.1,},
		{6.8,0.2,},
		{-2.2,-0.2,},
		{24.2,0.2,},
		{8.8,0.0,},
		{10.1,-0.3,},
		{-16.9,-0.6,},
		{-18.3,0.3,},
		{-3.2,0.5,},
		{13.3,0.1,},
		{-20.6,-0.2,},
		{-14.6,0.5,},
		{13.4,0.4,},
		{16.2,-0.2,},
		{11.7,0.1,},
		{5.7,-0.3,},
		{-15.9,-0.4,},
		{-9.1,0.3,},
		{-2.0,0.3,},
		{2.1,0.0,},
		{5.4,0.0,},
		{8.8,0.0,},
		{-21.6,0.0,},
		{3.1,0.0,},
		{10.8,0.0,},
		{-3.3,0.0,},
		{11.8,0.0,},
		{0.7,0.0,},
		{-6.8,0.0,},
		{-13.3,0.0,},
		{-6.9,0.0,},
		{-0.1,0.0,},
		{7.8,0.0,},
		{8.7,0.0,},
		{1.0,0.0,},
		{-9.1,0.0,},
		{-4.0,0.0,},
		{-10.5,0.0,},
		{8.4,0.0,},
		{-1.9,0.0,},
		{-6.3,0.0,},
		{3.2,0.0,},
		{0.1,0.0,},
		{-0.4,0.0,},
		{0.5,0.0,},
		{4.6,0.0,},
		{-0.5,0.0,},
		{4.4,0.0,},
		{1.8,0.0,},
		{-7.9,0.0,},
		{-0.7,0.0,},
		{-0.6,0.0,},
		{2.1,0.0,},
		{-4.2,0.0,},
		{2.4,0.0,},
		{-2.8,0.0,},
		{-1.8,0.0,},
		{-1.2,0.0,},
		{-3.6,0.0,},
		{-8.7,0.0,},
		{3.1,0.0,},
		{-1.5,0.0,},
		{-0.1,0.0,},
		{-2.3,0.0,},
		{2.0,0.0,},
		{2.0,0.0,},
		{-0.7,0.0,},
		{-0.8,0.0,},
		{-1.1,0.0,},
		{0.6,0.0,},
		{0.8,0.0,},
		{-0.7,0.0,},
		{-0.2,0.0,},
		{0.2,0.0,},
		{-2.2,0.0,},
		{1.7,0.0,},
		{-1.4,0.0,},
		{-0.2,0.0,},
		{-2.5,0.0,},
		{0.4,0.0,},
		{-2.0,0.0,},
		{3.5,0.0,},
		{-2.4,0.0,},
		{-1.9,0.0,},
		{-0.2,0.0,},
		{-1.1,0.0,},
		{0.4,0.0,},
		{0.4,0.0,},
		{1.2,0.0,},
		{1.9,0.0,},
		{-0.8,0.0,},
		{-2.2,0.0,},
		{0.9,0.0,},
		{0.3,0.0,},
		{0.1,0.0,},
		{0.7,0.0,},
		{0.5,0.0,},
		{-0.1,0.0,},
		{-0.3,0.0,},
		{0.3,0.0,},
		{-0.4,0.0,},
		{0.2,0.0,},
		{0.2,0.0,},
		{-0.9,0.0,},
		{-0.9,0.0,},
		{-0.1,0.0,},
		{0.0,0.0,},
		{0.7,0.0,},
		{0.0,0.0,},
		{-0.9,0.0,},
		{-0.9,0.0,},
		{0.4,0.0,},
		{0.4,0.0,},
		{0.5,0.0,},
		{1.6,0.0,},
		{-0.5,0.0,},
		{-0.5,0.0,},
		{1.0,0.0,},
		{-1.2,0.0,},
		{-0.2,0.0,},
		{-0.1,0.0,},
		{0.8,0.0,},
		{0.4,0.0,},
		{-0.1,0.0,},
		{-0.1,0.0,},
		{0.3,0.0,},
		{0.4,0.0,},
		{0.1,0.0,},
		{0.5,0.0,},
		{0.5,0.0,},
		{-0.3,0.0,},
		{-0.4,0.0,},
		{-0.4,0.0,},
		{-0.3,0.0,},
		{-0.8,0.0,}};

  int ii=0,jj=0;

      for ( jj = 0; jj < 2; jj++)
        {
          for (ii = 0; ii <= 195; ii++)
          {
        	  if (jj==0)
        	  {
        	  gh1[ii] = gg[ii][jj];
        	  }
        	  else
        	  {
        		  gh2[ii] = gg[ii][jj];
        	  }
          }
        }
  return;
}

/****************************************************************************/
/*                                                                          */
/*                           Subroutine julday                              */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Computes the decimal day of year from month, day, year.              */
/*     								                                        */
/****************************************************************************/

double julday(int month, int day, int year)
{

  int days[12] = { 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};

  int leap_year = (((year % 4) == 0) &&
                   (((year % 100) != 0) || ((year % 400) == 0)));

  double day_in_year = (days[month - 1] + day + (month > 2 ? leap_year : 0));

  return ((double)year + (day_in_year / (365.0 + leap_year)));
}

/****************************************************************************/
/*                                                                          */
/*                           Subroutine extrapsh                            */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Extrapolates linearly a spherical harmonic model with a              */
/*     rate-of-change model.                                                */
/*                                                                          */
/*     Input:                                                               */
/*           date     - date of resulting model (in decimal year)           */
/*           dte1     - date of base model                                  */
/*           nmax1    - maximum degree and order of base model              */
/*           gh1      - Schmidt quasi-normal internal spherical             */
/*                      harmonic coefficients of base model                 */
/*           nmax2    - maximum degree and order of rate-of-change model    */
/*           gh2      - Schmidt quasi-normal internal spherical             */
/*                      harmonic coefficients of rate-of-change model       */
/*                                                                          */
/*     Output:                                                              */
/*           gha or b - Schmidt quasi-normal internal spherical             */
/*                    harmonic coefficients                                 */
/*           nmax   - maximum degree and order of resulting model           */
/****************************************************************************/

void extrapsh(double dyear, double gh1[196], double gh2[196], double gha[196], double ghb[196])
{
double dte1 = 2015;
int   nmax1=13;
int   nmax2=13;
  int   k=0;
  int   i=0;
  double factor=0, factor2=0;

  factor = dyear - dte1;
  factor2 = (dyear+1) - dte1;
  k = nmax1 * (nmax2 + 2);
  for ( i = 1; i <= k; i++)
      {
       gha[i] = gh1[i] + (factor * gh2[i]);
      }
  for ( i = 1; i <= k; i++)
      {
      ghb[i] = gh1[i] + (factor2 * gh2[i]);
      }
  return;
}

/****************************************************************************/
/*                                                                          */
/*                           Subroutine interpsh                            */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Interpolates linearly, in time, between two spherical harmonic       */
/*     models.                                                              */
/*                                                                          */
/*     Input:                                                               */
/*           date     - date of resulting model (in decimal year)           */
/*           dte1     - date of earlier model                               */
/*           nmax1    - maximum degree and order of earlier model           */
/*           gh1      - Schmidt quasi-normal internal spherical             */
/*                      harmonic coefficients of earlier model              */
/*           dte2     - date of later model                                 */
/*           nmax2    - maximum degree and order of later model             */
/*           gh2      - Schmidt quasi-normal internal spherical             */
/*                      harmonic coefficients of internal model             */
/*                                                                          */
/*     Output:                                                              */
/*           gha or b - coefficients of resulting model                     */
/*           nmax     - maximum degree and order of resulting model         */

/****************************************************************************/

//int interpsh(double dyear, dte1, nmax1, dte2, nmax2, gh)
//{
//	 double date;
//     double dte1;
//     int   nmax1;
//     double dte2;
//     int   nmax2;
//     int   gh;
//
//  int   nmax;
//  int   k=0, l=0;
//  int   ii;
//  double factor=0,factor2=0;
//
//  factor = (dyear - dte1) / (dte2 - dte1);
//  factor2 = (dyear+1 - dte1) / (dte2 - dte1);
//  if (nmax1 == nmax2)
//    {
//      k =  nmax1 * (nmax1 + 2);
//      nmax = nmax1;
//    }
//  else
//    {
//      if (nmax1 > nmax2)
//        {
//          k = nmax2 * (nmax2 + 2);
//          l = nmax1 * (nmax1 + 2);
//  for ( ii = k + 1; ii <= l; ++ii)
//     {
//     gha[ii] = gh1[ii] + factor * (-gh1[ii]);
//     }
//  for ( ii = k + 1; ii <= l; ++ii)
//  	  {
//      ghb[ii] = gh1[ii] + factor * (-gh1[ii]);
//      }
//          nmax = nmax1;
//        }
//      else
//        {
//          k = nmax1 * (nmax1 + 2);
//          l = nmax2 * (nmax2 + 2);
//   for ( ii = k + 1; ii <= l; ++ii)
//      {
//      gha[ii] = factor * gh2[ii];
//      }
//   for ( ii = k + 1; ii <= l; ++ii)
//      {
//      ghb[ii] = factor * gh2[ii];
//      }
//      nmax = nmax2;
//        }
//    }
//  for ( ii = 1; ii <= k; ++ii)
//        {
//          gha[ii] = gh1[ii] + factor * (gh2[ii] - gh1[ii]);
//        }
//  for ( ii = 1; ii <= k; ++ii)
//        {
//          ghb[ii] = gh1[ii] + factor * (gh2[ii] - gh1[ii]);
//        }
//  return 0;
//}

/****************************************************************************/
/*                                                                          */
/*                           Subroutine shval3                              */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Calculates field components from spherical harmonic (sh)             */
/*     models.                                                              */
/*                                                                          */
/*     Input:                                                               */
/*           latitude  - north latitude, in degrees                         */
/*           longitude - east longitude, in degrees                         */
/*           elev      - WGS84 altitude above ellipsoid (igdgc=1), or       */
/*                       radial distance from earth's center (igdgc=2)      */
/*           a2,b2     - squares of semi-major and semi-minor axes of       */
/*                       the reference spheroid used for transforming       */
/*                       between geodetic and geocentric coordinates        */
/*                       or components                                      */
/*           nmax      - maximum degree and order of coefficients           */
/*           iext      - external coefficients flag (=0 if none)            */
/*           ext1,2,3  - the three 1st-degree external coefficients         */
/*                       (not used if iext = 0)                             */
/*                                                                          */
/*     Output:                                                              */
/*           x         - northward component                                */
/*           y         - eastward component                                 */
/*           z         - vertically-downward component                      */
/*                                                                          */
/****************************************************************************/


void shval3(double dlat, double dlon, int igdgc, double alt, double gha[196], double ghb[196],double ans[6], double dyear)
{
  int   nmax=13;
// WGS84
  double Re = 6371.2; // radius of earth
  double a2 = 40680631.59;            /* WGS84 */
  double b2 = 40408299.98;            /* WGS84 */
// preallocate vectors for calculations
  double sl[14] = {0};
  double cl[14] = {0};
  double p[119] = {0};
  double q[119] = {0};
//  double gtg = 0.9933;
  double dtr = 0.01745329; // pi / 180 == change from deg to radian
  double slat; // sine of latitude
  double clat; // cosine of latitude
  double ratio=0;
  double aa=0, bb=0, cc=0, dd=0;
  double sd=0.0;
  double cd=1.0;
  double r=0;
  double rr=0;
  double fm=0,fn=0;
  int ii,j,k,l,m,n;
  int npq=0;
  double argument=0;
  double power=0;

  // initialize outputs
  double x = 0, y = 0, z = 0, xtemp = 0, ytemp = 0, ztemp = 0;

  // initialize counters
  l = 1;
  n = 0;
  m = 1;

// convert geodetic latitude to geocentric latitude
//  if (igdgc == 2)
//  {
//	  dlat = atan2(gtg*tan(dlat*dtr),1);
//	  dlat = dlat / dtr;
//  }
  r = alt;
  argument = dlat * dtr;
  slat = sin( argument );
  if ((90.0 - dlat) < 0.001)
    {
      aa = 89.999;            /*  300 ft. from North pole  */
    }
  else
    {
      if ((90.0 + dlat) < 0.001)
        {
          aa = -89.999;        /*  300 ft. from South pole  */
        }
      else
        {
          aa = dlat;
        }
    }
  argument = aa * dtr;
  clat = cos( argument );
  argument = dlon * dtr;
  sl[1] = sin( argument );
  cl[1] = cos( argument );
  npq = (nmax * (nmax + 3)) / 2; // equals to 104.

//  convert geodetic to geocentric
  aa = a2 * clat * clat;
  bb = b2 * slat * slat;
  cc = aa + bb;
  dd = sqrt( cc );
  argument = alt * (alt + 2.0 * dd) + (a2 * aa + b2 * bb) / cc;
  r = sqrt( argument );
  cd = (alt + dd) / r;
  sd = (a2 - b2) / dd * slat * clat / r;
  argument = slat;
  slat = slat * cd - clat * sd;
  clat = clat * cd + argument * sd;

  ratio = Re / r;
  argument = sqrt( 3.0 );

  p[1] = 2.0 * slat;
  p[2] = 2.0 * clat;
  p[3] = 4.5 * slat * slat - 1.5;
  p[4] = 3.0 * argument * clat * slat;
  q[1] = -clat;
  q[2] = slat;
  q[3] = -3.0 * clat * slat;
  q[4] = argument * (slat * slat - clat * clat);

  for ( k = 1; k <= npq; k++)
    {
      if (n < m)
        {
          m = 0;
          n = n + 1;
          power =  n + 2;
          rr = pow(ratio,power);
          fn = n;
        }
      fm = m;
      if (k >= 5)
        {
          if (m == n)
            {
              aa = sqrt(1.0 - 0.5 / fm);
              j = k - n - 1;
              p[k] = (1.0 + 1.0/fm) * aa * clat * p[j];
              q[k] = aa * (clat * q[j] + slat/fm * p[j]);
              sl[m] = sl[m-1] * cl[1] + cl[m-1] * sl[1];
              cl[m] = cl[m-1] * cl[1] - sl[m-1] * sl[1];
            }
          else
            {
              aa = sqrt( fn*fn - fm*fm );
              bb = sqrt( ((fn - 1.0)*(fn-1.0)) - (fm * fm) ) / aa;
              cc = (2.0 * fn - 1.0)/aa;
              ii = k - n;
              j = k - 2 * n + 1;
              p[k] = (fn + 1.0) * (cc * slat/fn * p[ii] - bb/(fn - 1.0) * p[j]);
              q[k] = cc * (slat * q[ii] - clat/fn * p[ii]) - bb * q[j];
            }
        }
      aa = rr * gha[l];
      double aatemp = rr * ghb[l];

      if (m == 0)
        {
  	  	  x = x + aa * q[k];
          z = z - aa * p[k];
          xtemp = xtemp + aatemp * q[k];
          ztemp = ztemp - aatemp * p[k];
          l = l + 1;
        }
      else
        {
    	  bb = rr * gha[l+1];
          cc = aa * cl[m] + bb * sl[m];
          x = x + cc * q[k];
          z = z - cc * p[k];

      double bbtemp = rr * ghb[l+1];
      double cctemp = aatemp * cl[m] + bbtemp * sl[m];
      xtemp = xtemp + cctemp * q[k];
      ztemp = ztemp - cctemp * p[k];

      if (clat > 0)
       {
        y = y + (aa * sl[m] - bb * cl[m]) * fm * p[k]/((fn + 1.0) * clat);
        ytemp = ytemp + (aatemp * sl[m] - bbtemp *cl[m]) * fm * p[k] / ((fn + 1.0) * clat);
       }
       else
        {
    	  y = y + (aa * sl[m] - bb * cl[m]) * q[k] * slat;
          ytemp = ytemp + (aatemp * sl[m] - bbtemp * cl[m]) * q[k] * slat;
        }
              l = l + 2;
        }
      m = m + 1;
    }

  	  double xold = x;
  	  x = x * cd + z * sd;
  	  z = z * cd - xold * sd;

  	  double xtempold = xtemp;
  	  xtemp = xtemp * cd + ztemp * sd;
  	  ztemp = ztemp * cd - xtempold * sd;

  	  ans[0] = x; ans[1] = y; ans[2] = z;
  	  ans[3] = xtemp; ans[4] = ytemp; ans[5] = ztemp;

  	 cout << setw(12) << left << dyear;
  	 cout << setw(12) << left << x;
  	 cout << setw(12) << left << y;
  	 cout << setw(12) << left << z;
//  	  cout << "Xtemp-Coord: " << xtemp << endl;
//  	  cout << "Ytemp-Coord: " << ytemp << endl;
//  	  cout << "Ztemp-Coord: " << ztemp << endl;
  return;
}

/****************************************************************************/
/*                                                                          */
/*                           Subroutine dihf                                */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Computes the geomagnetic d, i, h, and f from x, y, and z.            */
/*                                                                          */
/*     Input:                                                               */
/*           x  - northward component                                       */
/*           y  - eastward component                                        */
/*           z  - vertically-downward component                             */
/*                                                                          */
/*     Output:                                                              */
/*           d  - declination                                               */
/*           i  - inclination                                               */
/*           h  - horizontal intensity                                      */
/*           f  - total intensity                                           */
/****************************************************************************/

void dihf (double ans[6])
{
  double d=0,i=0,dtemp=0,htemp=0,itemp=0,ftemp=0;
  double h2,h,f;
  double argument, argument2;
  double sn = 0.0001;
  double x = ans[0], y = ans[1], z = ans[2], xtemp = ans[3], ytemp = ans[4], ztemp = ans[5];

  h2 = x*x + y*y;
  h = sqrt( h2 );       /* calculate horizontal intensity */
  f = sqrt(h2 + z*z);      /* calculate total intensity */
  if (f < sn)
     {
      d = 0;        /* If d and i cannot be determined, */
      i = 0;        /*       set equal to NaN         */
     }
  else
     {
      i = atan2(z,h);
   if (h < sn)
     {
      d = 0;
     }
   else
     {
   if ((h+x) < sn)
     {
      d = PI;
     }
   else
     {
      d = 2.0 * atan2(y,(h+x));
     }
   }
}

  h2 = xtemp*xtemp + ytemp*ytemp;
  htemp = sqrt(h2);
  ftemp = sqrt(h2 + ztemp*ztemp);
  if (ftemp < sn)
  {
  dtemp = 0;    /* If d and i cannot be determined, */
  itemp = 0;    /*       set equal to 999.0         */
  }
  else
  {
  argument = ztemp;
  argument2 = htemp;
  itemp = atan2(argument,argument2);
  if (htemp < sn)
    {
    dtemp = 0;
    }
   else
    {
   if ((htemp + xtemp) < sn)
    {
    dtemp = PI;
    }
     else
    {
    argument = ytemp;
    argument2 = htemp + xtemp;
    dtemp = 2.0 * atan2(argument,argument2);
    }
    }
  }
	 cout << setw(15) << left << h;
	 cout << setw(15) << left << f;
	 cout << setw(15) << left << (i*180/PI) + itemp - itemp;
	 cout << setw(15) << left << (d*180/PI) + dtemp - dtemp << endl;
//  cout << "Horizontal Intensity 2: " << htemp << endl;
//  cout << "Total Intensity 2: " << ftemp << endl;
//  cout << "Inclination 2: " << itemp << endl;
//  cout << "declination 2: " << dtemp << endl;
  return;
}
