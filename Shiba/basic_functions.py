
def cell(x,y,z):
    return (abs(x)+abs(y)+abs(z))

def logger(s):
    logfile = 'out.log'
    print(s, end='')
    with open(logfile, 'a') as f:
        f.write(s)