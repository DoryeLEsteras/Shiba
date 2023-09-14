import os

def cell(x:int,y:int,z:int) -> int:
    return (abs(x)+abs(y)+abs(z))

def logger(s:str) -> None:
    logfile = 'out.log'
    print(s, end='')
    with open(logfile, 'a') as f:
        f.write(s)

def get_time(value) -> str:
    timestr = ''
    if(value == 0):
        timestr += '00'
    elif(value < 10):
        timestr += '0%s'%value
    else:
        timestr += str(value)
    return timestr

def stop_watch(start, finish) -> str:
    s_in_h = 3600
    s_in_min = 60
    elapsed = int(finish-start)
    mystr = ''
    mystr += 'Time wasted: '
    if(elapsed > 0):
        mystr += get_time(int(elapsed / s_in_h))
        mystr += ':'
        mystr += get_time(int((elapsed % s_in_h) / s_in_min))
        mystr += ':'
        mystr += get_time(int((elapsed % s_in_h) % (s_in_min)))
    else:
        mystr += '00:00:00'
    mystr += '\n'
    return mystr

def estimate_memory(Parameters,Wannier_h): 
   dbl=8
   cdbl=16
   Mmels=Wannier_h.norb*Wannier_h.norb*Wannier_h.nws*dbl
   Mh=(Wannier_h.blocks*4*(Wannier_h.norb//2))*(Wannier_h.blocks*4*(Wannier_h.norb//2))*cdbl
   MId=Wannier_h.blocks*4*(Wannier_h.norb//2)*cdbl
   MIdnorb2=(Wannier_h.norb//2)*cdbl
   Mhuu=(Wannier_h.norb//2)*(Wannier_h.norb//2)*cdbl
   Mgam1=(4*(Wannier_h.norb//2))*(4*(Wannier_h.norb//2))*cdbl
   Mspec=Parameters.vpts*dbl*3+Mh*3
   Mcond=Parameters.vpts*dbl*4+Mh*5
   Mtot=Mmels+Mh*24+MId*2+MIdnorb2+Mhuu*4+Mgam1*2+(Mspec*(Parameters.spec==1))+Mcond
   logger("Estimated memory requirement = %.2f GB"%(Mtot*1e-9) + '\n\n')
   #total_memory = melsup.size * melsup.itemsize * 1e-9
   #print(total_memory)