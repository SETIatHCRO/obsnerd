import datetime

TIME_FORMATS = ['%Y-%m-%dT%H:%M', '%y-%m-%dT%H:%M',
                '%Y-%m-%d %H:%M', '%y-%m-%d %H:%M',
                '%Y/%m/%dT%H:%M', '%y/%m/%dT%H:%M',
                '%Y/%m/%d %H:%M', '%y/%m/%d %H:%M',
                '%d/%m/%YT%H:%M', '%d/%m/%yT%H:%M',
                '%d/%m/%Y %H:%M', '%d/%m/%y %H:%M',
                '%Y%m%dT%H%M', '%y%m%dT%H%M',
                '%Y%m%d %H%M', '%y%m%d %H%M',
                '%Y%m%d_%H%M', '%y%m%d_%H%M',
                '%Y%m%d%H%M', '%y%m%d%H%M'
            ]


def strip_tz(this_datetime):
    if not isinstance(this_datetime, str):
        return this_datetime
    try:
        if this_datetime[-6] in ['+', '-']:
            return this_datetime[:-6]
    except IndexError:
        this_datetime


def create_tz(tz, default='server'):
    if isinstance(default, str):
        if default.lower() == 'server':
            default = datetime.datetime.now().astimezone().tzinfo
        elif default.lower() == 'utc':
            default = datetime.timezone(datetime.timedelta(0), 'UTC+00:00')
    if tz is None:
        return default
    if isinstance(tz, datetime.datetime):
        return tz.tzinfo
    if isinstance(tz, datetime.timezone):
        return tz

    if isinstance(tz, str) and ':' in tz:  # Provided as e.g. -08:00 or +02:00 or isoformat with tz
        if len(tz) == 6:
            sgn = tz[0]            
            hr, mn = [float(x) for x in tz[1:].split(':')]
        else:
            try:
                sgn = tz[-6]
                hr, mn = [float(x) for x in tz[-6:].split(':')]
            except (IndexError, ValueError):
                return default
        vsgn = 1.0 if sgn == '+' else -1.0
        hroffset = vsgn * (hr + mn / 60.0)
        name = f"UTC{tz}"
    else:
        from math import modf
        try:
            tz = float(tz)
        except (ValueError, TypeError):
            return default
        sgn = '-' if tz < 0.0 else '+'
        vsgn = 1.0 if sgn == '+' else -1.0
        hroffset = tz
        fhr, hr = modf(abs(tz))
        name = f"UTC{sgn}{int(hr)}:{int(60.0*fhr)}"
    return datetime.timezone(datetime.timedelta(hours=hroffset), name)


def make_datetime(**kwargs):
    """
    Take various datetime/str/offset/timezone options and return timezone-aware datetimes

    """
    for ptype in ['date', 'time', 'datetime', 'datestamp', 'timestamp', 'offset',
                  'tstart', 'tstop', 'tle', 'expected']:
        if ptype in kwargs:
            datetimes = kwargs[ptype]
            break
    else:
        raise ValueError("No valid datetime term included")
    if isinstance(datetimes, str):
        datetimes = datetimes.split(',')
    elif not isinstance(datetimes, list):
        datetimes = [datetimes]

    timezones = None
    for p in ['timezone', 'tz']:
        if p in kwargs:
            timezones = kwargs[p]
            break
    if isinstance(timezones, str):
        timezones = timezones.split(',')
    elif not isinstance(timezones, list):
        timezones = [timezones] * len(datetimes)

    if len(timezones) != len(datetimes):
        raise ValueError(f"Datetime lengths ({len(datetimes)}) differs from timezone lengths ({len(timezones)})")

    datetime_out = []
    for dt, tz in zip(datetimes, timezones):
        datetime_out.append(proc_datetime(dt, tz, ptype))

    if len(datetime_out) == 1:
        return datetime_out[0]
    else:
        return datetime_out


def proc_datetime(this_datetime, this_timezone, ptype=''):
    """
    Handles one datetime/timezone pair

    """
    this_timezone = create_tz(this_timezone, default=None)
    if this_timezone is None:
        this_timezone = create_tz(this_datetime, default='utc')

    # Process datetime value and timezone
    if this_datetime == 'now' or this_datetime is None:
        return datetime.datetime.now().astimezone(this_timezone)

    if isinstance(this_datetime, (float, int)):
        return datetime.datetime.now().astimezone(this_timezone) + datetime.timedelta(minutes=this_datetime)

    if isinstance(this_datetime, datetime.datetime):
        return this_datetime.replace(tzinfo=this_timezone)

    # ... it is a str, make sure no tz info left
    this_datetime = strip_tz(this_datetime)
    this_dt = None
    for this_tf in TIME_FORMATS:
        try:
            this_dt = datetime.datetime.strptime(this_datetime, this_tf)
            break
        except (TypeError, ValueError):
            try:
                if ':' in this_tf:
                    this_tf += ':%S'
                else:
                    this_tf += '%S'
                this_dt = datetime.datetime.strptime(this_datetime, this_tf)
                break
            except (TypeError, ValueError):
                try:
                    this_dt = datetime.datetime.strptime(this_datetime, this_tf+'.%f')
                    break
                except (TypeError, ValueError):
                    continue
    if isinstance(this_dt, datetime.datetime):
        return this_dt.replace(tzinfo=this_timezone)

    try:
        dt = float(this_datetime)
        return datetime.datetime.now().astimezone(this_timezone) + datetime.timedelta(minutes=dt)
    except (TypeError, ValueError):
        return None