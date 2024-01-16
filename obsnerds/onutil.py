import datetime

TIME_FORMATS = ['%Y-%m-%dT%H:%M:%S', '%y-%m-%dT%H:%M:%S',
                '%Y-%m-%d %H:%M:%S', '%y-%m-%d %H:%M:%S',
                '%Y/%m/%dT%H:%M:%S', '%y/%m/%dT%H:%M:%S',
                '%Y/%m/%d %H:%M:%S', '%y/%m/%d %H:%M:%S',
                '%d/%m/%YT%H:%M:%S', '%d/%m/%yT%H:%M:%S',
                '%d/%m/%Y %H:%M:%S', '%d/%m/%y %H:%M:%S',
                '%Y%m%dT%H%M%S', '%y%m%dT%H%M%S',
                '%Y%m%d %H%M%S', '%y%m%d %H%M%S',
                '%Y%m%d_%H%M%S', '%y%m%d_%H%M%S',
                '%Y%m%d%H%M%S', '%y%m%d%H%M%S'
            ]

def make_datetime(**kwargs):
    for p in ['date', 'time', 'datetime', 'datestamp', 'timestamp']:
        if p in kwargs:
            this_datetime = kwargs[p]
            break
        else:
            continue
    if not isinstance(this_datetime, (str, datetime.datetime)):
        return None

    timezone = None
    for p in ['timezone', 'tz']:
        try:
            if isinstance(kwargs[p], datetime.timezone):
                break
            else:
                hr = float(kwargs[p])
                name = f"UTC{'+' if hr>=0.0 else '-'}{hr:.0f}"
                timezone = datetime.timezone(datetime.timedelta(hours=hr), name)
                break
        except (ValueError, KeyError):
            continue
    
    if this_datetime == 'now':
        return datetime.datetime.now().astimezone(timezone)

    if isinstance(this_datetime, datetime.datetime):
        if this_datetime.tzinfo == timezone:
            return this_datetime
        return this_datetime.replace(tzinfo=timezone)

    this_dt = None
    for this_tf in TIME_FORMATS:
        try:
            this_dt = datetime.strptime(this_datetime, this_tf)
            break
        except ValueError:
            try:
                this_dt = datetime.strptime(this_datetime, this_tf+'.%f')
                break
            except ValueError:
                continue
    return this_dt.replace(tzinfo=timezone)
