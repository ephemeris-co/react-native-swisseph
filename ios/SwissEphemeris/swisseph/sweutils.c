#include "swephexp.h"
#include "sweutils.h"
#include "SwissEphGlue.h"

double get_search_interval_for_planet(int planet)
{
  switch (planet)
  {
  case SE_SUN:
    return 1.0;
  case SE_MOON:
    return 0.5;
  case SE_MERCURY:
    return 1.0;
  case SE_VENUS:
    return 1.0;
  case SE_MARS:
    return 1.0;
  case SE_JUPITER:
    return 3.0;
  case SE_SATURN:
    return 4.0;
  case SE_URANUS:
    return 10.0;
  case SE_NEPTUNE:
    return 15.0;
  case SE_PLUTO:
    return 30.0;
  default:
    return 1.0;
  }
}

double get_longitude_difference(double jd, int planet, double targetLongitude)
{
  double result[6];
  int flag = SEFLG_SPEED | SEFLG_JPLEPH;
  swe_calc_ut(jd, planet, flag, result, NULL);
  double longitude = result[0];
  return swe_difdeg2n(longitude, targetLongitude);
}

double secant_method(double (*func)(double, int32, double), double jd, int32 ipl, double target_longitude, double a, double b)
{
  double TOLERANCE = 1e-6;
  double MAX_ITERARIONS = 100;
  double iter = 0;
  double c = a;

  while (fabs(b - a) >= TOLERANCE && iter < MAX_ITERARIONS)
  {
    double a_result = func(a, ipl, target_longitude);
    double b_result = func(b, ipl, target_longitude);

    c = b - (b_result * (b - a)) / (b_result - a_result);

    if (func(c, ipl, target_longitude) == 0.0)
    {
      break;
    }
    else
    {
      a = b;
      b = c;
    }
    iter++;
  }

  return c;
}

double get_distance(double jd, int32 ipl, double target_longitude)
{
  return get_longitude_difference(jd, ipl, target_longitude);
}

double find_next_crossing(int ipl, double target_longitude, double start_jd)
{

  double jd = start_jd;
  double prevSign = signbit(get_distance(jd, ipl, target_longitude));
  double searchInterval = get_search_interval_for_planet(ipl);

  while (1)
  {
    jd += searchInterval;
    double currentSign = signbit(get_distance(jd, ipl, target_longitude));

    if (currentSign != prevSign)
    {
      return secant_method(&get_distance, jd, ipl, target_longitude, jd - searchInterval, jd);
    }

    prevSign = currentSign;
  }
}

double find_previous_crossing(int ipl, double target_longitude, double start_jd)
{

  double jd = start_jd;
  double prevSign = signbit(get_distance(jd, ipl, target_longitude));
  double searchInterval = get_search_interval_for_planet(ipl);

  while (1)
  {
    jd -= searchInterval;
    double currentSign = signbit(get_distance(jd, ipl, target_longitude));

    if (currentSign != prevSign)
    {
      return secant_method(&get_distance, jd, ipl, target_longitude, jd + searchInterval, jd);
    }

    prevSign = currentSign;
  }
}