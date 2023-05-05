

#ifdef __cplusplus
extern "C"
{
#endif

  typedef struct
  {
    double jd;
    double longitude;
    double latitude;
    double longitude_speed;
    double latitude_speed;
    double distance;
    double distance_speed;
  } Ephemeris;

  typedef struct
  {
    double start_jd;
    double end_jd;
    int ipl;
    Ephemeris start_ephemeris;
    Ephemeris end_ephemeris;
  } Retrograde;

  typedef struct
  {
    Retrograde *retrogrades;
    int count;
  } RetrogradesList;

#ifdef __cplusplus
}
#endif

double find_next_crossing(int ipl, double target_longitude, double start_jd);
double find_previous_crossing(int ipl, double target_longitude, double start_jd);
RetrogradesList get_new_retrogrades(double start_jd, double end_jd);
RetrogradesList get_new_directs(double start_jd, double end_jd);
Retrograde get_retrograde(double jd, int ipl);
Retrograde find_next_retrograde(double since_jd, int ipl);