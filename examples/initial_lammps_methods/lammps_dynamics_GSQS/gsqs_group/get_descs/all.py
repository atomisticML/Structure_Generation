from subprocess import call

for i in range(1,10+1):
    call('python get_target.py sample.%d.dat' % i,shell=True)
    
