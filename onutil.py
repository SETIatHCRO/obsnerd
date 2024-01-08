from datetime import datetime, timedelta

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
    this_datetime = None
    for p in ['date', 'time', 'datetime', 'datestamp', 'timestamp']:
        if p in kwargs:
            this_datetime = kwargs[p]
            break
        else:
            continue
    print("ONUTIL22")
    print(this_datetime)
    if not isinstance(this_datetime, (str, datetime)):
        return None
    timezone = 0.0
    for p in ['timezone', 'tz']:
        try:
            timezone = float(kwargs[p])
            break
        except (ValueError, KeyError):
            continue
    if isinstance(this_datetime, datetime):
        return this_datetime + timedelta(hours=timezone)

    this_dt = None
    for this_tf in TIME_FORMATS:
        try:
            this_dt = datetime.strptime(this_datetime, this_tf) + timedelta(hours=timezone)
            break
        except ValueError:
            continue
    return this_dt
