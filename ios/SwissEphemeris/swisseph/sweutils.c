#include <stdio.h>
#include <stdlib.h>
#include "swephexp.h"
#include "sweutils.h"
#include "SwissEphGlue.h"

const int NUMBER_OF_RETROGRADATION_PLANETS = 8;
const int RETROGRADATION_IPLS[NUMBER_OF_RETROGRADATION_PLANETS] = {SE_MERCURY, SE_VENUS, SE_MARS, SE_JUPITER, SE_SATURN, SE_URANUS, SE_NEPTUNE, SE_PLUTO};

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

double secant_method_longitude(double (*func)(double, int32, double), double jd, int32 ipl, double target_longitude, double a, double b)
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
      return secant_method_longitude(&get_distance, jd, ipl, target_longitude, jd - searchInterval, jd);
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
      return secant_method_longitude(&get_distance, jd, ipl, target_longitude, jd, jd + searchInterval);
    }

    prevSign = currentSign;
  }
}

Ephemeris create_ephemeris(double jd, int ipl)
{

  double result[6];
  int flag = SEFLG_SPEED | SEFLG_JPLEPH;

  swe_calc_ut(jd, ipl, flag, result, NULL);

  Ephemeris ephemeris;

  ephemeris.jd = jd;
  ephemeris.longitude = result[0];
  ephemeris.latitude = result[1];
  ephemeris.longitude_speed = result[3];
  ephemeris.latitude_speed = result[4];
  ephemeris.distance = result[2];
  ephemeris.distance_speed = result[5];

  return ephemeris;
}

double secant_method_longitude_speed(double (*func)(double, int32), double jd, int32 ipl, double a, double b)
{
  double TOLERANCE = 1e-6;
  double MAX_ITERATIONS = 100;
  double iter = 0;
  double c = a;

  while (fabs(b - a) >= TOLERANCE && iter < MAX_ITERATIONS)
  {
    double a_result = func(a, ipl);
    double b_result = func(b, ipl);

    c = b - (b_result * (b - a)) / (b_result - a_result);

    if (func(c, ipl) == 0.0)
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

double get_longitude_speed(double jd, int ipl)
{
  Ephemeris ephemeris = create_ephemeris(jd, ipl);
  return ephemeris.longitude_speed;
}

double find_retrogradation_start(int ipl, double start_jd)
{

  double jd = start_jd;
  double prevSign = signbit(get_longitude_speed(jd, ipl));

  while (1)
  {
    jd -= 1;
    double currentSign = signbit(get_longitude_speed(jd, ipl));

    if (currentSign != prevSign)
    {
      return secant_method_longitude_speed(&get_longitude_speed, jd, ipl, jd, jd + 1);
    }

    prevSign = currentSign;
  }
}

double find_retrogradation_end(int ipl, double start_jd)
{

  double jd = start_jd;
  double prevSign = signbit(get_longitude_speed(jd, ipl));

  while (1)
  {
    jd += 1;
    double currentSign = signbit(get_longitude_speed(jd, ipl));

    if (currentSign != prevSign)
    {
      return secant_method_longitude_speed(&get_longitude_speed, jd, ipl, jd - 1, jd);
    }

    prevSign = currentSign;
  }
}

RetrogradesList get_new_retrogrades_in_day(double jd)
{

  int32 current_year;
  int32 current_month;
  int32 current_day;
  int32 current_hour;
  int32 current_minute;
  double current_second;
  swe_jdut1_to_utc(jd, SE_GREG_CAL, &current_year, &current_month, &current_day, &current_hour, &current_minute, &current_second);

  double start_of_day_jd;
  swe_utc_to_jd(current_year, current_month, current_day, 0, 0, 0.0, SE_GREG_CAL, &start_of_day_jd, NULL);

  double end_of_day_jd = start_of_day_jd + 1.0;

  RetrogradesList new_retrogrades;
  new_retrogrades.count = 0;
  new_retrogrades.retrogrades = malloc(sizeof(Retrograde));

  int count = 0;
  int i;

  for (i = 0; i < NUMBER_OF_RETROGRADATION_PLANETS; i++)
  {
    int ipl = RETROGRADATION_IPLS[i];
    Ephemeris start_of_day_ephemeris = create_ephemeris(start_of_day_jd, ipl);
    Ephemeris end_of_day_ephemeris = create_ephemeris(end_of_day_jd, ipl);

    if (start_of_day_ephemeris.longitude_speed > 0 && end_of_day_ephemeris.longitude_speed < 0)
    {

      double retrograde_start_jd = find_retrogradation_start(ipl, end_of_day_jd);
      double retrograde_end_jd = find_retrogradation_end(ipl, end_of_day_jd);
      Ephemeris start_ephemeris = create_ephemeris(retrograde_start_jd, ipl);
      Ephemeris end_ephemeris = create_ephemeris(retrograde_end_jd, ipl);

      Retrograde retrograde;
      retrograde.ipl = ipl;
      retrograde.start_jd = retrograde_start_jd;
      retrograde.end_jd = retrograde_end_jd;
      retrograde.start_ephemeris = start_ephemeris;
      retrograde.end_ephemeris = end_ephemeris;

      int size = count + 1;
      new_retrogrades.retrogrades = realloc(new_retrogrades.retrogrades, sizeof(Retrograde) * size);

      new_retrogrades.retrogrades[count] = retrograde;

      count++;

      new_retrogrades.count = count;
    }
  }

  return new_retrogrades;
}

RetrogradesList get_new_directs_in_day(double jd)
{

  int32 current_year;
  int32 current_month;
  int32 current_day;
  int32 current_hour;
  int32 current_minute;
  double current_second;
  swe_jdut1_to_utc(jd, SE_GREG_CAL, &current_year, &current_month, &current_day, &current_hour, &current_minute, &current_second);

  double start_of_day_jd;
  swe_utc_to_jd(current_year, current_month, current_day, 0, 0, 0.0, SE_GREG_CAL, &start_of_day_jd, NULL);

  double end_of_day_jd = start_of_day_jd + 1.0;

  RetrogradesList new_directs;
  new_directs.count = 0;
  new_directs.retrogrades = malloc(sizeof(Retrograde));

  int count = 0;
  int i = 0;

  for (i = 0; i < NUMBER_OF_RETROGRADATION_PLANETS; i++)
  {

    int ipl = RETROGRADATION_IPLS[i];
    Ephemeris start_of_day_ephemeris = create_ephemeris(start_of_day_jd, ipl);
    Ephemeris end_of_day_ephemeris = create_ephemeris(end_of_day_jd, ipl);

    if (start_of_day_ephemeris.longitude_speed < 0 && end_of_day_ephemeris.longitude_speed > 0)
    {

      double retrograde_start_jd = find_retrogradation_start(ipl, start_of_day_jd);
      double retrograde_end_jd = find_retrogradation_end(ipl, start_of_day_jd);
      Ephemeris start_ephemeris = create_ephemeris(retrograde_start_jd, ipl);
      Ephemeris end_ephemeris = create_ephemeris(retrograde_end_jd, ipl);

      Retrograde retrograde;
      retrograde.ipl = ipl;
      retrograde.start_jd = retrograde_start_jd;
      retrograde.end_jd = retrograde_end_jd;
      retrograde.start_ephemeris = start_ephemeris;
      retrograde.end_ephemeris = end_ephemeris;

      int size = count + 1;
      new_directs.retrogrades = realloc(new_directs.retrogrades, sizeof(Retrograde) * size);

      new_directs.retrogrades[count] = retrograde;
      count++;
      new_directs.count = count;
    }
  }

  return new_directs;
}

void free_retrogrades_list(RetrogradesList retrogrades_list)
{
  free(retrogrades_list.retrogrades);
}

RetrogradesList get_new_retrogrades(double start_jd, double end_jd)
{
  double current_jd;

  RetrogradesList new_retrogrades;
  new_retrogrades.count = 0;
  new_retrogrades.retrogrades = malloc(sizeof(Retrograde));

  int count = 0;
  for (current_jd = start_jd; current_jd < end_jd; current_jd++)
  {

    RetrogradesList new_retrogrades_in_day = get_new_retrogrades_in_day(current_jd);

    if (new_retrogrades_in_day.count == 0)
    {
      continue;
    }

    int i;
    for (i = 0; i < new_retrogrades_in_day.count; i++)
    {
      int size = count + 1;
      new_retrogrades.retrogrades = realloc(new_retrogrades.retrogrades, sizeof(Retrograde) * size);
      new_retrogrades.retrogrades[count] = new_retrogrades_in_day.retrogrades[i];
      count++;
      new_retrogrades.count = count;
    }

    free(new_retrogrades_in_day.retrogrades);
  }

  return new_retrogrades;
}

RetrogradesList get_new_directs(double start_jd, double end_jd)
{
  double current_jd;

  RetrogradesList new_directs;
  new_directs.count = 0;
  new_directs.retrogrades = malloc(sizeof(Retrograde));

  int count = 0;
  for (current_jd = start_jd; current_jd < end_jd; current_jd++)
  {

    RetrogradesList new_directs_in_day = get_new_directs_in_day(current_jd);

    if (new_directs_in_day.count == 0)
    {
      continue;
    }

    int i;
    for (i = 0; i < new_directs_in_day.count; i++)
    {
      int size = count + 1;
      new_directs.retrogrades = realloc(new_directs.retrogrades, sizeof(Retrograde) * size);
      new_directs.retrogrades[count] = new_directs_in_day.retrogrades[i];
      count++;
      new_directs.count = count;

      free(new_directs_in_day.retrogrades);
    }
  }

  return new_directs;
}