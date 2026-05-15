import qiskit
import scipy
from qiskit.primitives import Sampler
from qiskit.primitives import Estimator as Estimator
from qiskit.quantum_info import SparsePauliOp
import stim
import openfermion
from copy import deepcopy
from numpy.random import choice
import numpy as np
import math
import sys
import os



#tableau class has a lot of baggage from previous uses but here I just use it to produce the diagonalization circuits

class tableau: 
    def __init__(self, init_stabs):
        self.init_stabs=init_stabs
        self.init_stabslist=self.stabslist_from_strings(self.init_stabs)
        self.curr_stabs=self.init_stabslist
        self.transformed_stabs=deepcopy(self.curr_stabs)
        self.n=len(self.curr_stabs[0])
        self.m=len(self.init_stabslist)
        self.circuit=[]
        self.rank=0
        self.qorder=[i for i in range(self.n)]

    def stabslist_from_strings(self, stabs):
        stabslist=[]
        for stab in stabs:
            stabslist.append(stim.PauliString(str(stab)))
        return stabslist

    def stabslist_from_tab(self, tab):
        stabslist=[]
        for i in range(len(tab)):
            stabslist.append(tab.z_output(i))
        return stabslist

    def h(self, q): 
        for stab in self.curr_stabs:
            if stab[q]==1:
                stab[q]=3
            elif stab[q]==3:
                stab[q]=1
            elif stab[q]==2:
                stab*=-1
        return self.curr_stabs

    def s(self, q):
        for stab in self.curr_stabs:
            if stab[q]==1:
                stab[q]=2
            elif stab[q]==2:
                stab[q]=1
                stab*=-1

        return self.curr_stabs

    def cx(self, c, t): #this is a hard-coded CNOT gate. Hadamard and Phase (S) gates above.
        for stab in self.curr_stabs:
            if stab[c]==0:
                if stab[t]==2:
                    stab[c]=3
                elif stab[t]==3:
                    stab[c]=3
            elif stab[c]==1:
                if stab[t]==0:
                    stab[t]=1
                elif stab[t]==1:
                    stab[t]=0
                elif stab[t]==2:
                    stab[c]=2
                    stab[t]=3
                elif stab[t]==3:
                    stab[c]=2
                    stab[t]=2
                    stab*=-1
            elif stab[c]==2:
                if stab[t]==0:
                    stab[t]=1
                elif stab[t]==1:
                    stab[t]=0
                elif stab[t]==2:
                    stab[c]=1
                    stab[t]=3
                    stab*=-1
                elif stab[t]==3:
                    stab[c]=1
                    stab[t]=2
            elif stab[c]==3:
                if stab[t]==2:
                    stab[c]=0
                elif stab[t]==3:
                    stab[c]=0


        return self.curr_stabs

#the three functions below perform tableau row/column swaps and row multiplications

    def rowswap(self, i, j):

        newj=self.curr_stabs[i]
        newi=self.curr_stabs[j]
        self.curr_stabs[i]=newi
        self.curr_stabs[j]=newj

        return self.curr_stabs

    def columnswap(self, i, j):
        qorderi=self.qorder[j]
        qorderj=self.qorder[i]
        self.qorder[i]=qorderi
        self.qorder[j]=qorderj
        for k in range(len(self.curr_stabs)):
            s=self.curr_stabs[k]
            newj=s[i]
            newi=s[j]
            news=s
            news[j]=newj
            news[i]=newi
            self.curr_stabs[k]=news

        return self.curr_stabs

    def rowsum(self, i, j):
        if i == j:
            print('stabilizer '+str(i)+': you trynna kill me??')
            return
        self.curr_stabs[j]=self.curr_stabs[i]*self.curr_stabs[j]
        return self.curr_stabs


    def h_remove(self):
        swap = False
        for q in range(self.n):
            if swap:
                break
            for i in range(self.m):
                if self.curr_stabs[i][q]>0:
                    self.rowswap(i,0)
                    swap = True
                    break

        for i in range(1, len(self.curr_stabs)):
            self.rowsum(0, i)


        return self.curr_stabs


    def upper_triang(self):
        k=1
        swap = False
        for j in range(self.n):
            for i in range(k, len(self.curr_stabs)):
                swap = False
                if self.curr_stabs[i][j]!=0:
                    self.rowswap(i, k)
                    k+=1
                    swap = True
                    break

            if swap:
                for i in range(k, len(self.curr_stabs)):
                    if self.curr_stabs[i][j]==self.curr_stabs[k-1][j]:
                        self.rowsum(k-1,i)


        for i in reversed(range(len(self.curr_stabs))):
            if self.curr_stabs[i]==stim.PauliString(self.n):
                del self.curr_stabs[i]


        for i in reversed(range(len(self.curr_stabs))):
            for j in range(self.n):
                pauli = -1
                if self.curr_stabs[i][j]!=0:
                    pauli=self.curr_stabs[i][j]
                for k in reversed(range(i)):
                    if self.curr_stabs[k][j]==pauli:
                        self.rowsum(i,k)


        return self.curr_stabs






    def cz_clear(self):
        k = self.rank
        for i in range(1,k):
            for j in range(i):
                if self.curr_stabs[i][j]>1:
                    self.cz(i,j)
                    self.circuit.append('CZ '+str(self.qorder[i])+' '+str(self.qorder[j]))
        for i in range(k):
            if self.curr_stabs[i][i]>1:

                self.s(i)
                self.circuit.append('S '+str(self.qorder[i]))
            self.h(i)
            self.circuit.append('H '+str(self.qorder[i]))

        return self.curr_stabs





    def cx_clear(self):
        k = self.rank
        for i in range(k):
            s=0
            for j in range(i+1):
                if self.curr_stabs[i][j]>1:
                    s+=1
            if s%2==0:
                self.s(i)
                self.circuit.append('S '+str(self.qorder[i]))
            for j in range(i):
                if self.curr_stabs[i][j]>1:
                    self.cx(i,j)
                    self.circuit.append('CX '+str(self.qorder[i])+' '+str(self.qorder[j]))
                    self.rowsum(i,j)
        for i in range(k):
            self.s(i)
            self.circuit.append('S '+str(self.qorder[i]))
            self.h(i)
            self.circuit.append('H '+str(self.qorder[i]))


        return self.curr_stabs

    def local_diagonalize(self):
        for q in range(self.n):
            pauli=0
            for i in range(len(self.curr_stabs)):
                if self.curr_stabs[i][q]!=0 and pauli == 0:
                    pauli = self.curr_stabs[i][q]

                    continue
                elif self.curr_stabs[i][q]!=0 and pauli != self.curr_stabs[i][q]:
                    pauli=0

                    break


            if pauli == 3:
                continue

            elif pauli == 1:
                self.h(q)
                self.circuit.append('H '+str(q))


            elif pauli == 2:
                self.s(q)
                self.circuit.append('S '+str(q))

                self.h(q)
                self.circuit.append('H '+str(q))


        return self.curr_stabs

    def cx_diagonalize(self):
        for i in range(1, len(self.curr_stabs)):
            q1=q2=-1
            for j in range(self.n):
                if 3>self.curr_stabs[i][j]>0 and q1 == -1:
                    q1 = j
                elif 3>self.curr_stabs[i][j]>0 and q1 != -1:
                    q2 = j
                    break
            if q2!=-1:
                if self.curr_stabs[i][q2]==2:
                    self.s(q2)
                    self.circuit.append('S '+str(q2))
                self.cx(q1,q2)
                self.circuit.append('CX '+str(q1)+' '+str(q2))

                if self.curr_stabs[i][q1]==2:
                    self.s(q1)
                    self.circuit.append('S '+str(q1))

                self.h(q1)
                self.circuit.append('H '+str(q1))


        return self.curr_stabs

    def x_diag(self):
        k = 0
        found=True
        while  found == True:
            found = False
            for i in range(k, self.m):
                for j in range(k, self.n):
                    if 3>self.curr_stabs[i][j]>0:

                        found = True
                        self.rowswap(i,k)
                        self.columnswap(j,k)

                        for l in list(range(k))+list(range(k+1, self.m)):
                            if 3>self.curr_stabs[l][k]>0:
                                self.rowsum(k,l)
                        k+=1
                        break
                    if found == True:
                        break

        kx=k
        found=True
        while  found == True:
            found = False
            for i in range(k, self.m):
                for j in range(k, self.n):
                    if self.curr_stabs[i][j]>1:

                        found = True
                        self.rowswap(i,k)

                        self.columnswap(j,k)
                        for l in list(range(k))+list(range(k+1, self.m)):
                            if self.curr_stabs[l][k]>1:
                                self.rowsum(k,l)

                        k+=1
                        break
                    if found == True:
                        break


        for j in range(kx, k):
            self.h(j)
            self.circuit.append('H '+str(self.qorder[j]))
        for i in range(k):
            for j in range(k, self.n):
                if 3>self.curr_stabs[i][j]>0:
                    self.cx(i,j)
                    self.circuit.append('CX '+str(self.qorder[i])+' '+str(self.qorder[j]))

        self.rank=k

        return

    def cx_clear(self):
        k = self.rank
        for i in range(k):
            s=0
            for j in range(i+1):
                if self.curr_stabs[i][j]>1:
                    s+=1
            if s%2==0:
                self.s(i)
                self.circuit.append('S '+str(self.qorder[i]))
            for j in range(i):
                if self.curr_stabs[i][j]>1:
                    self.cx(i,j)
                    self.circuit.append('CX '+str(self.qorder[i])+' '+str(self.qorder[j]))
                    self.rowsum(i,j)
        for i in range(k):
            self.s(i)
            self.circuit.append('S '+str(self.qorder[i]))
            self.h(i)
            self.circuit.append('H '+str(self.qorder[i]))




        return


    def diagonalize(self):
        self.h_remove()
        self.upper_triang()
        self.local_diagonalize()
        self.cx_diagonalize()
        self.local_diagonalize()

        self.curr_stabs = self.transformed_stabs

        self.do_circuit()

        success='SUCCESS'
        for i in self.curr_stabs:
            for j in i:
                if 3>j>0:
                    success='FAILURE'
                    print(success)
                    return


        return self.circuit, self.curr_stabs



    def diagonalize_general(self):
        self.x_diag()
        self.cx_clear()
        self.curr_stabs = self.transformed_stabs
        self.do_circuit()
        success='SUCCESS'
        for i in self.curr_stabs:
            for j in i:
                if 3>j>0:
                    success='FAILURE'
                    print(success)
                    return


        return self.circuit, self.curr_stabs

    def do_circuit(self):
        for gate in self.circuit:
            if gate.split()[0] == "H":
                self.h(int(gate.split()[1]))
            elif gate.split()[0] == "CX":
                self.cx(int(gate.split()[1]), int(gate.split()[2]))
            elif gate.split()[0] == "S":
                self.s(int(gate.split()[1]))
        return



def oftermtoqpauli(ofterm, n): #openfermion term to qiskit-compatible string (zeroth qubit is right-most in output)
    qplist=['I']*n
    for nonI in ofterm:
        qplist[nonI[0]]=nonI[1]

    return '+'+''.join(qplist)

def stimtermtoqpauli(stimterm, n): #stim term to qiskit-compatible string (zeroth qubit is right-most in output)
    qplist=['I']*n
    for i in range(len(stimterm)):
        if stimterm[i]==0:
            continue
        elif stimterm[i]==3:
            qplist[i]='Z'
        else:
            print('NONDIAGONAL PAULI')
            return

    if stimterm.sign.real<0:
        sign='-'
    elif stimterm.sign.real>0:
        sign='+'

    return sign+''.join(qplist)

def ofoptoqop(ofop, n): #openfermion operator to qiskit operator
    paulis=[]
    coeffs=[]
    for term in ofop.terms:
        paulis.append(oftermtoqpauli(term, n))
        coeffs.append(ofop.terms[term])
    return SparsePauliOp(paulis, coeffs)

def evaleval(bitstring, zstring): #return eigenvalue (+-1) given bitstring and zstring
    if len(zstring) != len(bitstring)+1:
        print("Z-string and bitstring lengths don't match up")
        return
    if zstring[0]!='+' and zstring[0]!='-':
        print("Z-string missing sign")
        return
    evalue=1
    for i in range(1, len(zstring)):
        if zstring[i]=='Z' and bitstring[i-1]=='1':
            evalue*=-1
    if zstring[0] == '-':
        evalue*=-1
    return evalue


def zpauliexp(circuit, shots, listofpaulis): #calculates expectation values for a list of diagonal paulis
    expectations=[]
    circ_copy=deepcopy(circuit)
    circ_copy.measure_all()
    job = Sampler().run(circ_copy, shots=shots)
    job_result = job.result()
    for pauli in listofpaulis:
        expectation=0
        for bitstring in list(job_result.quasi_dists[0].binary_probabilities()):
            expectation+=job_result.quasi_dists[0].binary_probabilities()[bitstring]*evaleval(bitstring, pauli)
        expectations.append(expectation)
    return(expectations)



def circ2qcirc(circlist, qcirc): #take qiskit circuit, append gates in circlist to it, and spit out final circuit
    fcirc=deepcopy(qcirc)
    for gate in circlist:
        if gate.split()[0] == "H":
            fcirc.h(int(gate.split()[1]))
        elif gate.split()[0] == "CX":
            fcirc.cx(int(gate.split()[1]), int(gate.split()[2]))
        elif gate.split()[0] == "S":
            fcirc.s(int(gate.split()[1]))
    return fcirc

def pool(n):
    pool_list=[]
    for first in range(n):
        for second in range(first+1, n):
            for third in range(second+1, n):
                for fourth in range(third+1,n):
                    if (first+second+third+fourth)%2 == 0:
                        doublist=[]
                        doublist.append('X%d X%d X%d Y%d ' %(first, second, third, fourth))
                        doublist.append('Y%d Y%d Y%d X%d ' %(first, second, third, fourth))
                        doublist.append('X%d X%d X%d Y%d ' %(second, third, fourth, first))
                        doublist.append('Y%d Y%d Y%d X%d ' %(second, third, fourth, first))
                        doublist.append('X%d X%d X%d Y%d ' %(third, fourth, first, second))
                        doublist.append('Y%d Y%d Y%d X%d ' %(third, fourth, first, second))
                        doublist.append('X%d X%d X%d Y%d ' %(fourth, first, second, third))
                        doublist.append('Y%d Y%d Y%d X%d ' %(fourth, first, second, third))

                        for double_string in doublist:
                            double = openfermion.QubitOperator(double_string, 0 + 1j)
                            if double not in pool_list:
                                pool_list.append(double)
    return pool_list

def poolgroups(n):
    pool_list=[]
    for first in range(n):
        Ygroup=[]
        Xgroup=[]
        for second in range(n):
            for third in range(second+1, n):
                for fourth in range(third+1,n):
                    if (first+second+third+fourth)%2 == 0 and (first != second) and (first != third) and (first != fourth):

                        ydstring=('Y%d X%d X%d X%d ' %(first, second, third, fourth))
                        xdstring=('X%d Y%d Y%d Y%d ' %(first, second, third, fourth))



                        xdouble = openfermion.QubitOperator(xdstring, 0 + 1j)
                        ydouble = openfermion.QubitOperator(ydstring, 0 + 1j)

                        Ygroup.append(ydouble)
                        Xgroup.append(xdouble)

        pool_list.append(Ygroup)
        pool_list.append(Xgroup)

    return pool_list

def poolgroup(xy, anchor, n):


    group=[]
    first=anchor
    for second in range(n):
        for third in range(second+1, n):
            for fourth in range(third+1,n):
                if (first+second+third+fourth)%2 == 0 and (first != second) and (first != third) and (first != fourth):
                    if xy == 'Y':
                        dstring=('Y%d X%d X%d X%d ' %(first, second, third, fourth))
                    elif xy == 'X':
                        dstring=('X%d Y%d Y%d Y%d ' %(first, second, third, fourth))
                    else:
                        print('bad input')
                        return
                    double = openfermion.QubitOperator(dstring, 0 + 1j)

                    group.append(double)

    return group

def diagonalize_single(pstring):
    if pstring[0]!="+" and pstring[0]!="-":
        print("missing sign")
        return
    qplist=['I']*len(pstring)
    if pstring[0] == '+':
        qplist[0] = '+'
    if pstring[0] == '-':
        qplist[0] = '-'
    circlist=[]
    for i in range(1, len(pstring)):
        if (pstring[i] == "I") or (pstring[i] == "Z"):
            qplist[i] = pstring[i]
        elif pstring[i] == "X":
            circlist.append('H '+str(len(pstring)-i-1))

            qplist[i] = "Z"
        elif pstring[i] == "Y":
            circlist.append('S '+str(len(pstring)-i-1))
            circlist.append('H '+str(len(pstring)-i-1))
            qplist[i] = "Z"
            if qplist[0]=="-":
                qplist[0]='+'
            else:
                qplist[0]='-'
    return ''.join(qplist), circlist



def measure_single(circ, pauli, n, shots): #state, op to measure, n qubits
    pstringnd = oftermtoqpauli(list(pauli.terms)[0], n)
    pstring, circlist = diagonalize_single(pstringnd)
    circ_c = deepcopy(circ)
    for gate in circlist:
        if gate.split()[0] == "H":
            circ_c.h(int(gate.split()[1]))
        elif gate.split()[0] == "S":
            circ_c.s(int(gate.split()[1]))
        elif gate.split()[0] == "CX":
            circ_c.cx(int(gate.split()[1]), int(gate.split()[2]))

    exp = zpauliexp(circ_c, shots, [pstring])
    return exp[0]



def measure_multiple(circ, diagcirclist, paulilist, n, shots): #state, diagonalizing_circuit, ops to measure, n qubits, shots
    qpaulilist=[]
    for term in paulilist:
        qpaulilist.append(stimtermtoqpauli(term, n))
    circ_c = deepcopy(circ)
    for gate in diagcirclist:
        if gate.split()[0] == "H":
            circ_c.h(n-1-int(gate.split()[1]))
        elif gate.split()[0] == "S":
            circ_c.s(n-1-int(gate.split()[1]))
        elif gate.split()[0] == "CX":
            circ_c.cx(n-1-int(gate.split()[1]), n-1-int(gate.split()[2]))

    exp = zpauliexp(circ_c, shots, qpaulilist)
    return exp


def keytostim(key, n):
    stkey=stim.PauliString(n)
    for nonI in key:
        stkey[nonI[0]]=nonI[1]
    return stkey


def cov_key(key1, key2):
    n1=int.from_bytes(str(key1).encode())
    n2=int.from_bytes(str(key2).encode())
    if n1<n2:
        return str(key2)+str(key1)
    else:
        return str(key1)+str(key2)




class poolgradients:
    def __init__(self, h, pool, circ, spstate, n, meas_ratio): # h and pool (list) are OpenFermion QubitOperators, circ is a qiskit state, n is an integer (n of qubits)
        self.circ=circ
        self.h=h
        self.spstate=spstate
        self.pool=pool
        self.energy={}
        self.gradients_first={}
        self.hterms_first={}
        self.observables={}
        self.observables[()]=[1, 0, 1]
        self.observables_cov={}
        self.energy_setup = False
        self.gradients_setup = False
        self.hterms_setup = False
        self.sorted_grads_setup = False
        self.sorted_energy_setup = False
        self.n = n
        self.e_ind_shots=0
        self.g_ind_shots=0
        self.g_sim_shots=0
        self.e_sort_shots=0
        self.g_sort_shots=0
        self.guess_var=1
        self.meas_ratio=meas_ratio


    def energy_setup_dict(self):
        self.energy=self.h.terms
        for term in self.energy:
            if term not in self.observables:
                self.observables[term]=[0, self.guess_var, 1]
        self.energy_setup = True




    def measure_energy_independent(self, precision):

        print(' ')
        print(' ')
        print(' ')
        print('Now measuring energy term by term.')
        print(' ')
        print(' ')
        print(' ')
        print('Measurement ratios: '+str(self.meas_ratio))
        print(' ')
        print(' ')
        if self.energy_setup == False:
            self.energy_setup_dict()

        list_of_ops=list(self.energy)
        batch_est=[]
        for term in list_of_ops:
            coeff=self.energy[term]
            obs=self.observables[term]
            batch_est.append(abs(coeff)*np.sqrt(obs[1]))
        rounds=0
        evar=self.energy_variance_independent()
        while evar > precision**2:
            for ratio in self.meas_ratio:
                print('Ratio: '+str(ratio))
                if evar < precision**2:
                    break
                for opind in range(len(list_of_ops)):
                    if evar < precision**2:
                        break
                    draw=list_of_ops[opind]
                    rounds+=1
                    coeff=self.energy[draw]
                    obs=self.observables[draw]
                    batch_size=math.ceil(ratio*((sum(batch_est).real*batch_est[opind]/precision**2))-obs[2])
                    if batch_size < 1:
                        continue

                    exp=measure_single(self.circ, openfermion.QubitOperator(draw), self.n, batch_size)
                    self.e_ind_shots+=batch_size
                    obs[0]=((obs[0]*(obs[2]-1))+(exp*batch_size))/(obs[2]-1+batch_size)
                    obs[2]+=batch_size
                    obs[1]=(1-obs[0]**2)
                    self.observables[draw]=obs
                    batch_est[opind]=(abs(coeff)*np.sqrt(obs[1]))

                    evar=self.energy_variance_independent()


                    print('Iter:', rounds, 'Term:', openfermion.QubitOperator(draw), 'E variance:', evar, 'Tot Shots:', self.e_ind_shots, 'Iter Shots:', batch_size)



        return


    def energy_variance_independent(self):
        self.energy_var=0
        for term in self.energy:
            termobs=self.observables[term]
            self.energy_var += (((self.energy[term])**2).real)*termobs[1]/termobs[2]
        return self.energy_var

    def energy_eval(self):
        self.E=0
        for hterm in self.energy:
            self.E += ((self.energy[hterm]).real)*(self.observables[hterm][0])
        exact=(self.spstate.transpose().conjugate()@openfermion.linalg.get_sparse_operator(self.h)@self.spstate)[0,0].real

        print("Measured Energy: ", self.E, "Error: ", self.E-exact)
        return

    def gradients_setup_dict(self):
        for op in self.pool:
            comm=openfermion.utils.commutator(op, self.h)
            self.gradients_first[str(op)]=comm.terms
            for commterm in self.gradients_first[str(op)]:
                if commterm not in self.observables:
                    self.observables[commterm]=[0, self.guess_var, 1]
        self.gradients_setup = True
        return

    def measure_gradients_first_independent(self, precision):

        print(' ')
        print(' ')
        print(' ')
        print('Now measuring gradients term by term.')
        print(' ')
        print(' ')
        print(' ')
        print('Measurement ratios: '+str(self.meas_ratio))
        print(' ')
        print(' ')

        if self.gradients_setup == False:
            self.gradients_setup_dict()

        for op in self.pool:
            list_of_ops=list(self.gradients_first[str(op)])
            batch_est=[]
            for term in list_of_ops:
                coeff=self.gradients_first[str(op)][term]
                obs=self.observables[term]
                batch_est.append(abs(coeff)*np.sqrt(obs[1]))
            rounds=0
            gvar=self.grad_variance_independent(str(op))
            while gvar > precision**2:
                for ratio in self.meas_ratio:
                    print('Ratio: '+str(ratio))
                    if gvar < precision**2:
                        break
                    for opind in range(len(list_of_ops)):
                        if gvar < precision**2:
                            break
                        rounds+=1
                        draw=list_of_ops[opind]
                        obs=self.observables[draw]
                        coeff=self.gradients_first[str(op)][draw]
                        batch_size=math.ceil(ratio*((sum(batch_est).real*batch_est[opind]/precision**2))-obs[2])
                        if batch_size < 1:
                                continue
                        exp=measure_single(self.circ, openfermion.QubitOperator(draw), self.n, batch_size)
                        self.g_ind_shots+=batch_size
                        obs[0]=((obs[0]*(obs[2]-1))+(exp*batch_size))/((obs[2]-1)+batch_size)
                        obs[2]+=batch_size
                        obs[1]=(1-obs[0]**2)
                        self.observables[draw]=obs
                        batch_est[opind]=(abs(coeff)*np.sqrt(obs[1]))
                        gvar=self.grad_variance_independent(str(op))

                        print('Operator: ', str(op),'Iter: ',  rounds, 'Grad Var: ', self.grad_variance_independent(str(op)), 'Tot Shots: ',  self.g_ind_shots,'Iter Shots: ',  batch_size, flush = True)

        return



    def grad_variance_independent(self, op):
        grad_var=0
        for term in self.gradients_first[op]:
            termobs=self.observables[term]
            grad_var += ((self.gradients_first[op][term]).real**2)*termobs[1]/termobs[2]

        return grad_var

    def max_grad_variance_independent(self, poolgroup):
        max_grad_var=0
        for op in poolgroup:
            grad_var = self.grad_variance_independent(str(op))
            if grad_var > max_grad_var:
                max_grad_var = grad_var

        return max_grad_var

    def gradients_first_eval(self):

        gradients=[]
        for op in self.gradients_first:
            grad=0
            for gterm in self.gradients_first[op]:
                grad += ((self.gradients_first[op][gterm]).real)*(self.observables[gterm][0])
            exact=(self.spstate.transpose().conjugate()@openfermion.linalg.get_sparse_operator(openfermion.utils.commutator(openfermion.QubitOperator(op),self.h))@self.spstate)[0,0].real

            print(op, "Measured grad: ", grad, "Error: ", grad-exact)
        return

    def hterms_setup_dict(self):
        for hterm in self.energy:
            if hterm not in self.hterms_first:
                self.hterms_first[hterm]={}

            for op in self.pool:
                comm=openfermion.utils.commutator(op, openfermion.QubitOperator(hterm))
                if comm.terms:
                    comm*=self.energy[hterm]
                    commterm=list(comm.terms)[0]
                    self.hterms_first[hterm][commterm]=comm.terms[commterm]
                    if commterm not in self.observables:
                        self.observables[commterm]=[0, self.guess_var, 1]
        self.hterms_setup = True
        return

    def measure_hterms_first(self, poolgroup, precision):


        print(' ')
        print(' ')
        print(' ')
        print('Now simultaneously measuring gradients, term by term of the Hamiltonian.')
        print(' ')
        print(' ')
        print(' ')
        print('Measurement ratios: '+str(self.meas_ratio))
        print(' ')
        print(' ')


        for i in range(len(poolgroup)):
            for j in range(i+1, len(poolgroup)):
                if openfermion.utils.commutator(poolgroup[i], poolgroup[j]).terms:
                    print('ops in group do not commute')
                    return

        if self.energy_setup == False:
            self.energy_setup_dict()
        if self.hterms_setup == False:
            self.hterms_setup_dict()
        if self.gradients_setup == False:
            self.gradients_setup_dict()

        rounds=0
        list_of_hterms=list(self.energy)
        opstomeasure=[]
        gradvarsum=0
        for op in poolgroup:
            gradvar=self.grad_variance_independent(str(op))
            gradvarsum+=gradvar
            if gradvar > precision**2:
                opstomeasure.append(op)

        batch_est=[]
        for term in list_of_hterms:
            term_sum=0
            coeff2=(abs(self.energy[term]))**2
            for op in opstomeasure:
                comm=openfermion.utils.commutator(openfermion.QubitOperator(term), openfermion.QubitOperator(str(op)))
                if comm.terms:
                    term_sum+=(coeff2*self.observables[list(comm.terms)[0]][1])
            batch_est.append(math.sqrt(term_sum))



        while len(opstomeasure)>0:
            for ratio in self.meas_ratio:
                print('Ratio: '+str(ratio))
                if not opstomeasure:
                    break

                for termind in range(len(list_of_hterms)):
                    if not opstomeasure:
                        break
                    rounds+=1
                    draw=list_of_hterms[termind]
                    stimobslist=[]
                    obslist=[]
                    batch_sum=sum(batch_est)
                    batch_size=0

                    termkeys=list(self.hterms_first[draw].keys())
                    min_sample=float('inf')
                    for key in termkeys:
                        obslist.append(key)
                        stimobslist.append(keytostim(key, self.n))
                        if self.observables[key][2] < min_sample:
                            min_sample=self.observables[key][2]

                    if min_sample == float('inf'):
                        continue

                    batch_size=math.ceil((4*ratio*batch_sum.real*batch_est[termind]/(len(opstomeasure)*precision**2))-min_sample)

                    if batch_size < 1:
                        continue

                    if len(obslist) == 0:
                        continue


                    tab=tableau(stimobslist)
                    circlist, diag_obs = tab.diagonalize()
                    ncx=0
                    for g in circlist:
                        if g[0:2]=='CX':
                            ncx+=1
                    #coeff=self.energy[draw]
                    exp=measure_multiple(self.circ, circlist, diag_obs, self.n, batch_size)
                    self.g_sim_shots+=batch_size
                    for j in range(len(obslist)):
                        obs=self.observables[obslist[j]]
                        obs[0]=((obs[0]*(obs[2]-1))+(exp[j]*batch_size))/((obs[2]-1)+batch_size)
                        obs[2]+=batch_size
                        obs[1]=(1-obs[0]**2)#/obs[3]
                        self.observables[obslist[j]]=obs
                    print('Iter:', rounds,'H Term:',  openfermion.QubitOperator(draw),'Total Var:', gradvarsum ,'Tot Shots:',  self.g_sim_shots,'Iter Shots:',  batch_size, 'CX:', ncx)




                    newopstomeasure=[]
                    gradvarsum=0
                    for op in opstomeasure:
                        gradvar=self.grad_variance_independent(str(op))
                        gradvarsum+=gradvar
                        if gradvar > precision**2:
                            newopstomeasure.append(op)

                    if opstomeasure == newopstomeasure:
                        term_sum=0
                        coeff=self.energy[draw]
                        for op in opstomeasure:
                            comm=openfermion.utils.commutator(openfermion.QubitOperator(draw), openfermion.QubitOperator(str(op)))
                            if comm.terms:
                                term_sum+=((abs(coeff)**2)*self.observables[list(comm.terms)[0]][1])
                        batch_est[termind]=(math.sqrt(term_sum))

                    else:
                        batch_est=[]
                        for term in list_of_hterms:
                            term_sum=0
                            coeff=self.energy[term]
                            for op in opstomeasure:
                                comm=openfermion.utils.commutator(openfermion.QubitOperator(term), openfermion.QubitOperator(str(op)))
                                if comm.terms:
                                    term_sum+=((abs(coeff)**2)*self.observables[list(comm.terms)[0]][1])
                            batch_est.append(math.sqrt(term_sum))
                            opstomeasure = newopstomeasure





        return







    def sorted_grad_setup_dict(self):
        totgroups=0
        if self.gradients_setup == False:
            self.gradients_setup_dict()
        self.sorted_grads={}
        for op in self.pool:
            self.sorted_grads[str(op)]=[]
            groups=[[]]
            terms=list(self.gradients_first[str(op)])
            abscoeffs=[]
            for gradterm in terms:
                abscoeffs.append(abs(self.gradients_first[str(op)][gradterm]))
            while len(abscoeffs)>=1:
                indmax=np.argmax(abscoeffs)
                if groups[-1] != []:
                    groups.append([])
                for group in groups:
                    dnc = False
                    for term in group:
                        if abs(openfermion.utils.commutator(openfermion.QubitOperator(terms[indmax]), openfermion.QubitOperator(term)).induced_norm()) != 0 :
                            dnc=True
                            break

                    if dnc == False:
                        group.append(terms[indmax])
                        terms.pop(indmax)
                        abscoeffs.pop(indmax)
                        break


            if groups[-1] == []:
                groups.pop(-1)
            for group in groups:
                for i in range(len(group)):
                    term1=group[i]
                    for j in range(i+1, len(group)):
                        term2=group[j]
                        prod=list((openfermion.QubitOperator(term1)*openfermion.QubitOperator(term2)).terms)[0]
                        if prod not in self.observables:
                            self.observables[prod]=[0, self.guess_var, 1]
                        self.observables_cov[cov_key(term1, term2)]={str(term1):0, str(term2):0, str(prod):0, 'samples':0}


            self.sorted_grads[str(op)]=groups
            totgroups+=len(groups)


        self.sorted_grads_setup = True
        print('Total number of groups: ', totgroups)

        return

    def grad_variance_sorted(self, op):
        grad_var=0
        opterms=list(self.gradients_first[op])
        for i in range(len(opterms)):
            term1=opterms[i]
            term1obs=self.observables[term1]
            grad_var += ((self.gradients_first[op][term1]).real**2)*term1obs[1]/term1obs[2]
            for j in range(i+1, len(opterms)):
                term2=opterms[j]
                if cov_key(term1, term2) in self.observables_cov:
                    covdict=self.observables_cov[cov_key(term1, term2)]
                    joint_samples=covdict['samples']
                    term1exp=covdict[str(term1)]
                    term2exp=covdict[str(term2)]
                    term2obs=self.observables[term2]

                    prod=(openfermion.QubitOperator(term1)*openfermion.QubitOperator(term2))
                    prodterm=list(prod.terms)[0]
                    prodsign=prod.terms[prodterm].real
                    prodexp=covdict[str(prodterm)]
                    grad_var += 2*((self.gradients_first[op][term1]*self.gradients_first[op][term2]).real)*(joint_samples/(term1obs[2]*term2obs[2]))*(prodsign*prodexp-(term1exp*term2exp))

        return grad_var



    def measure_gradients_sorted(self, precision):

        print(' ')
        print(' ')
        print(' ')
        print('Now measuring gradients one by one, using sorted insertion.')
        print(' ')
        print(' ')
        print(' ')
        print('Measurement ratios: '+str(self.meas_ratio))
        print(' ')
        print(' ')

        if self.gradients_setup == False:
            self.gradients_setup_dict()

        if self.sorted_grads_setup == False:
            self.sorted_grad_setup_dict()

        for op in self.pool:
            groups=self.sorted_grads[str(op)]
            batch_est=[]
            for grpind in range(len(groups)):
                group=groups[grpind]
                group_batch=0
                for i in range(len(group)):
                    term1=group[i]
                    term1obs=self.observables[term1]
                    group_batch += ((self.gradients_first[str(op)][term1].real)**2)*term1obs[1]

                    for j in range(i+1, len(group)):
                        term2=group[j]
                        term2obs=self.observables[term2]
                        covdict=self.observables_cov[cov_key(term1, term2)]
                        term1exp=covdict[str(term1)]
                        term2exp=covdict[str(term2)]
                        prod=(openfermion.QubitOperator(term1)*openfermion.QubitOperator(term2))
                        prodterm=list(prod.terms)[0]
                        prodsign=prod.terms[prodterm].real
                        prodexp=covdict[str(prodterm)]
                        group_batch += 2*((self.gradients_first[str(op)][term1]*self.gradients_first[str(op)][term2]).real)*((covdict['samples']**2)/(term1obs[2]*term2obs[2]))*(prodsign*prodexp-(term1exp*term2exp))
                if group_batch > 0:
                    batch_est.append(group_batch)
                else:
                    batch_est.append(0)







            rounds=0
            #while self.grad_variance_independent(str(op)) > precision**2:
            gvar=self.grad_variance_sorted(str(op))
            while gvar > precision**2:
                for ratio in self.meas_ratio:
                    print('Ratio: '+str(ratio))
                    if gvar <= precision**2:
                        break
                    for grpind in range(len(groups)):
                        if gvar <= precision**2:
                            break

                        rounds+=1
                        draw=groups[grpind]

                        stimobslist=[]
                        obslist=[]
                        batchsum=0
                        for g in range(len(batch_est)):
                            batchsum+=math.sqrt(batch_est[g])

                        min_sample=float('inf')
                        for i in range(len(draw)):
                            if self.observables[draw[i]][2] < min_sample:
                                min_sample=self.observables[draw[i]][2]







                        batch_size=math.ceil(ratio*((batchsum*math.sqrt(batch_est[grpind]))/precision**2)-min_sample)


                        if batch_size < 1:
                            continue



                        for i in range(len(draw)):
                            term1=draw[i]
                            stimobslist.append(keytostim(term1, self.n))
                            obslist.append(term1)
                            for j in range(i+1, len(draw)):
                                term2=draw[j]
                                prod=list((openfermion.QubitOperator(term1)*openfermion.QubitOperator(term2)).terms)[0]
                                stimobslist.append(keytostim(prod, self.n))
                                obslist.append(prod)



                        if len(obslist) == 0:
                            continue



                        tab=tableau(stimobslist)
                        circlist, diag_obs = tab.diagonalize_general()
                        ncx=0
                        for g in circlist:
                            if g[0:2]=='CX':
                                ncx+=1
                        exp=measure_multiple(self.circ, circlist, diag_obs, self.n, batch_size)
                        self.g_sort_shots+=batch_size
                        for i in range(len(obslist)):
                            obs=self.observables[obslist[i]]
                            obs[0]=((obs[0]*(obs[2]-1))+(exp[i]*batch_size))/((obs[2]-1)+batch_size)
                            obs[2]+=batch_size
                            obs[1]=(1-obs[0]**2)
                            self.observables[obslist[i]]=obs
                            for j in range(i+1, len(obslist)):
                                if obslist[j] not in draw:
                                  continue
                                if obslist[i] not in draw:
                                  continue
                                covdict=self.observables_cov[cov_key(obslist[i], obslist[j])]
                                group_samples=covdict['samples']
                                term1exp=covdict[str(obslist[i])]
                                term2exp=covdict[str(obslist[j])]
                                prod=list((openfermion.QubitOperator(obslist[i])*openfermion.QubitOperator(obslist[j])).terms)[0]
                                prodexp=covdict[str(prod)]
                                term1exp=((term1exp*(group_samples))+(exp[i]*batch_size))/((group_samples)+batch_size)
                                term2exp=((term2exp*(group_samples))+(exp[j]*batch_size))/((group_samples)+batch_size)
                                prodexp=((prodexp*(group_samples))+(exp[obslist.index(prod)]*batch_size))/((group_samples)+batch_size)
                                covdict[str(obslist[i])]=term1exp
                                covdict[str(obslist[j])]=term2exp
                                covdict[str(prod)]=prodexp
                                covdict['samples']=group_samples+batch_size
                                self.observables_cov[cov_key(obslist[i], obslist[j])]=covdict


                        group_batch=0
                        for i in range(len(draw)):
                            term1=draw[i]
                            term1obs=self.observables[term1]
                            group_batch += ((self.gradients_first[str(op)][term1]).real**2)*term1obs[1]
                            for j in range(i+1, len(draw)):
                                term2=draw[j]
                                term2obs=self.observables[term2]
                                covdict=self.observables_cov[cov_key(term1, term2)]
                                term1exp=covdict[str(term1)]
                                term2exp=covdict[str(term2)]
                                prod=openfermion.QubitOperator(term1)*openfermion.QubitOperator(term2)
                                prodterm=list(prod.terms)[0]
                                prodsign=prod.terms[prodterm].real
                                prodexp=covdict[str(prodterm)]
                                group_batch += 2*((self.gradients_first[str(op)][term1]*self.gradients_first[str(op)][term2]).real)*((covdict['samples']**2)/(term1obs[2]*term2obs[2]))*(prodsign*prodexp-(term1exp*term2exp))


                        batch_est[grpind]=group_batch





                        gvar=self.grad_variance_sorted(str(op))


                        print('Grad: ', str(op),'Iter: ',  rounds,'Grad Var: ',  gvar,'Tot Shots: ',  self.g_sort_shots,'Iter Shots: ',  batch_size, 'CX: ', ncx, flush = True)

        return






