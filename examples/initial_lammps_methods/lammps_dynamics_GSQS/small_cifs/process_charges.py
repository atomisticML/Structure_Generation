import numpy as np

def get_charges(nats,f='charges.dat'):
    q_qeq = []
    q_exact = []
    q_approx = []
    chi = []
    with open(f,'r') as readin:
        lines = readin.readlines()
        for line in lines:
            l = line.split()
            llst = [float(k) for k in line.split()]
            qall = llst[1:]
            q_qeqi=qall[0:nats]
            q_qexacti=qall[nats:2*nats]
            #q_ni=qall[2*nats:3*nats]
            #chii=qall[3*nats:4*nats]
            q_qeqi = np.array(q_qeqi)
            q_qexacti = np.array(q_qexacti)
            #q_ni = np.array(q_ni)
            #q_fixedi = np.array(q_fixedi)
            #chii = np.array(chii)
            q_qeq.append(q_qeqi)
            #q_fixed.append(q_fixedi)
            q_exact.append(q_qexacti)
            #q_approx.append(q_ni)
            #chi.append(chii)
    #return np.array(q_qeq).T, np.array(q_exact).T, np.array(q_approx).T, np.array(q_fixed).T ,np.array(chi).T
    return np.array(q_qeq).T, np.array(q_exact).T 
