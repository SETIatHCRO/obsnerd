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

def make_datetime(**kwargs):
    """
    Take various datetime/str/offset/timezone options and return timezone-aware datetimes

    """
    for p in ['date', 'time', 'datetime', 'datestamp', 'timestamp', 'offset']:
        if p in kwargs:
            datetimes = kwargs[p]
            break
    else:
        raise InputError("No valid datetime term included")
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
        datetime_out.append(proc_datetime(dt, tz))

    if len(datetime_out) == 1:
        return datetime_out[0]
    else:
        return datetime_out


def proc_datetime(this_datetime, this_timezone):
    """
    Handles one datetime/timezone pair

    """
    tzindt = None
    if isinstance(this_datetime, str):
        try:
            sgn = this_datetime[-6]
        except IndexError:
            sgn = False
        if  sgn in ['+', '-']:
            vsgn = 1.0 if sgn == '+' else -1.0
            a, b = [float(x) for x in this_datetime[-5:].split(':')]
            name = f"UTC{sgn}{abs(a):.0f}"
            tzindt = datetime.timezone(datetime.timedelta(hours=vsgn * (a + b / 60.0)), name)
            this_datetime = this_datetime[:-6]
    elif isinstance(this_datetime, datetime.datetime):
        tzindt = this_datetime.tzinfo

    # Get timezone
    if not isinstance(this_timezone, datetime.timezone):
        try:
            hr = float(this_timezone)
            name = f"UTC{'+' if hr>=0.0 else '-'}{hr:.0f}"
            this_timezone = datetime.timezone(datetime.timedelta(hours=hr), name)
        except (ValueError, KeyError, TypeError):
            this_timezone = None

    # If supplied, the this_timezone overrides any timezone sent in this_datetime
    if this_timezone is None:
        this_timezone = tzindt

    # Process datetime value and timezone
    if this_datetime == 'now' or this_datetime is None:
        return datetime.datetime.now().astimezone(this_timezone)

    if isinstance(this_datetime, (float, int)):
        return datetime.datetime.now().astimezone(this_timezone) + datetime.timedelta(minutes=this_datetime)

    if isinstance(this_datetime, datetime.datetime):
        return this_datetime.replace(tzinfo=this_timezone)

    # ... it is a str
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