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
    Take various datetime/str/offset/timezone options and return timezone-aware datetime

    """
    # Get datetime value
    for p in ['date', 'time', 'datetime', 'datestamp', 'timestamp', 'offset']:
        if p in kwargs:
            this_datetime = kwargs[p]
            break
        else:
            continue
    tzindt = None
    if isinstance(this_datetime, str):
        try:
            sgn = this_datetime[-6]
        except IndexError:
            sgn = False
        if  sgn in ['+', '-']:
            sgn = 1.0 if this_datetime[-6] == '+' else -1.0
            a, b = [float(x) for x in this_datetime[-5:].split(':')]
            name = f"UTC{'+' if sgn > 0.0 else '-'}{a:.0f}"
            tzindt = datetime.timezone(datetime.timedelta(hours=sgn * (a + b / 60.0)), name)
            this_datetime = this_datetime[:-6]

    # Get timezone
    timezone = None
    for p in ['timezone', 'tz']:
        try:
            if isinstance(kwargs[p], datetime.timezone):
                timezone = kwargs[p]
                break
            else:
                hr = float(kwargs[p])
                name = f"UTC{'+' if hr>=0.0 else '-'}{hr:.0f}"
                timezone = datetime.timezone(datetime.timedelta(hours=hr), name)
                break
        except (ValueError, KeyError, TypeError):
            continue
    if timezone is None:
        timezone = tzindt

    # Process datetime value and timezone
    if this_datetime == 'now' or this_datetime is None:
        return datetime.datetime.now().astimezone(timezone)

    try:
        dt = float(this_datetime)
        return datetime.datetime.now().astimezone(timezone) + datetime.timedelta(minutes=dt)
    except (TypeError, ValueError):
        pass

    if isinstance(this_datetime, datetime.datetime):
        return this_datetime.replace(tzinfo=timezone)

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
    if not isinstance(this_dt, datetime.datetime):
        return None
    return this_dt.replace(tzinfo=timezone)
