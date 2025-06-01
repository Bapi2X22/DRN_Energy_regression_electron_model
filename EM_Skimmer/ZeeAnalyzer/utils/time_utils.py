from  datetime import datetime

def get_utc_time():
    tstmp_ = datetime.utcnow()
    tstr_ = '{}/{}/{} {}:{}:{}'.format(tstmp_.year, tstmp_.month, tstmp_.day,
                                      tstmp_.hour, tstmp_.minute, tstmp_.second)
    return tstr_

def print_creation_timestamp():
    return 'Created on {} UTC'.format(get_utc_time())
