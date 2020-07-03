# -*- coding: utf-8 -*-
"""
Created on Tue May  5 17:35:23 2020

@author: Maia

En esta versión intentamos optimizar el algoritmo:
    - No llamamos a la funcion clean()
    
"""

import numpy as np
import json
from copy import deepcopy
from collections import OrderedDict
import math

class tree(OrderedDict):
    """Autovivified dictionary."""
    
    def __missing__(self, key):
        value = self[key] = type(self)()
        return value

pprint = lambda Q: print(json.dumps(Q, indent = 4))

hits = lambda r, v: list(filter(lambda n: v[n]==r,range(len(v)))) 


def p2o_v1(x,y,z,A): 
    """Forms a number in base 8.

    Keyword arguments:
    x -- the x component (integer)
    y -- the y component (integer)
    z -- the z component
    A -- length (integer)
    """
    xmax = 2**A-1  
    px = lambda N_b8 : int(N_b8.replace('1','0').replace('2','0').replace('3','0').replace('4','1').replace('5','1').replace('6','1').replace('7','1'),2)  # proy x (b8 -> b10)
    py = lambda N_b8 : int(N_b8.replace('1','0').replace('2','1').replace('3','1').replace('4','0').replace('5','0').replace('6','1').replace('7','1'),2)  # proy y (b8 -> b10)
    pz = lambda N_b8 : int(N_b8.replace('1','1').replace('2','0').replace('3','1').replace('4','0').replace('5','1').replace('6','0').replace('7','1'),2)  # proy z (b8 -> b10)

    v2s = lambda v : str(v)[1:-1].replace(', ','')                                # vector -> string

    IP = lambda x,y,z,a : int(v2s([x//a,y//a,z//a]),2) 

    assert 0 <= x <= xmax,"it must be the case that 0 <= x <= xmax"
    assert 0 <= y <= xmax,"it must be the case that 0 <= y <= ymax"
    assert 0 <= z <= xmax,"it must be the case that 0 <= z <= zmax"

    ov = [0]*A           
    ov[0] = IP(x,y,z,2**(A-1))

    for n in range(1,A):
        os = v2s(ov) # os es la version string de ov
        ov[n] = IP(x - px(os),y - py(os),z - pz(os),2**(A-1-n)) 
        # en cada etapa del ciclo for obtenemos un digito más del numero de celda

    return v2s(ov)


def p2o(x,y,z,A):
    """Forms a number in base 8.

    Keyword arguments:
    x -- the x component (integer)
    y -- the y component (integer)
    z -- the z component
    A -- length (integer)
    """
    xmax = 2**A-1
    assert 0 <= x <= xmax,"it must be the case that 0 <= x <= xmax"
    assert 0 <= y <= xmax,"it must be the case that 0 <= y <= ymax"
    assert 0 <= z <= xmax,"it must be the case that 0 <= z <= zmax"
    
    b = len(bin(2**A-1))
        
    bx = "{{:0{:d}d}}".format(b).format(int(bin(x)[2:]))
    by = "{{:0{:d}d}}".format(b).format(int(bin(y)[2:]))
    bz = "{{:0{:d}d}}".format(b).format(int(bin(z)[2:]))
            
    bx = bx.replace('1','4')
    by = by.replace('1','2')
       
    q = "{{:0{:d}d}}".format(A).format(int(bx) + int(by) + int(bz))
    return q

    
def MASS(Q):
    """Adds all masses in Q.
    
    Keyword arguments:
    Q -- Tree   
    """
    masa = 0
    for k1,v1 in Q.items():
        for k2,v2 in v1.items():
            masa += v2['M']
    return masa    


def combine(I,*L):
    """Combines identical parameters in the list L.
    
    Keyword arguments:
    I -- Parameter identifier (integer)
    L -- list of dictionaries
    """
    R = {
    'M': L[0]['M'] + L[1]['M'],
    'S': L[0]['S'] + L[1]['S'],
    'X': np.round((L[0]['X'] * L[0]['M'] + L[1]['X'] * L[1]['M'])/(L[0]['M'] + L[1]['M'])),
    'Y': np.round((L[0]['Y'] * L[0]['M'] + L[1]['Y'] * L[1]['M'])/(L[0]['M'] + L[1]['M'])),
    'Z': np.round((L[0]['Z'] * L[0]['M'] + L[1]['Z'] * L[1]['M'])/(L[0]['M'] + L[1]['M'])),
    'ancs': [L[0]['I']] + [L[1]['I']],
    'nodes': L[0]['nodes'] + L[1]['nodes'],
    'I': I}
    for m in L[2:]: # combinacion de 3+ cuerpos
        temp = R['ancs'] + [m['I']]
        R = combine(I,R,m) 
        R['ancs'] = temp
    return R   


def promotion(C,a):
    """Promotes bodies in trees to lower levels.
    
    Keyword arguments:
    C -- Tree  
    a -- level of the tree (integer)
    """
    cindx = lambda n,m: list(map(lambda x: x[:m],C[n].keys()))  # cindex(nivel,nº decimales), el nivel x tiene x+1 decimales
    for k,v in C[a+1].items():  # k es de nivel: a+1
        r = k[:-1]                # r es de nivel: a
        H1 = hits(r,cindx(a+1,a+1)) # examino posibles choques/coincidencias al nivel: a+1 si recorto
        H0 = hits(r,cindx(a,a+1))   # examino posibles choques/coincidencias al nivel: a 
        # los vectores H contienen las posiciones de los cuerpos de cindx que son iguales a r 
        if len(H1) == 1 and len(H0) == 0:
            aux = cindx(a+1,a+2)[H1[0]]
            C[a][r] = deepcopy(C[a+1][aux]) # necesito desligar C[a] de C[a+1]
            C[a+1][aux]['del'] = 1  # marco para limpiar (delete)
    return C


def clean(C):
    """Erases bodies from trees.
    
    Keyword arguments:
    C -- Tree
    """
    SPQR = deepcopy(C)
    for k1,v1 in C.items():
        for k2,v2 in v1.items():
            for k3,v3 in v2.items():
                if k3 == 'del':
                    del SPQR[k1][k2]  # borro
    C = SPQR
    return C


def combination(D,a,P,I):
    """Adds new level and its information to a tree.
    
    Keyword arguments:
    D -- Tree  
    a -- level of the tree (integer)
    P -- Dicctionary
    I -- Parameter identifier (integer)
    """
    assert a >=1,"en Q no hay llaves negativas"
    cindx = lambda n,a: list(map(lambda x: x[:a],D[n].keys()))  # cindex(nivel,nº decimales) devuelve las posiciones de todos los cuerpos del nivel a elegir con el numero de decimales a elegir
    aux = cindx(a,a) # lista de posiciones de los cuerpos que contiene el nivel a
     # asi si el p se hace más pequeño por que se han juntado masas lo tenemos en cuenta
    
    for u in list(set(aux)): # indices unicos del nivel a (inferior) para rellenar
        H = hits(u,aux) # indices en los que las posiciones de aux son iguales a u
        if len(H) > 1:
            I += 1
            D[a-1][u] = combine(I,*map(lambda n : D[a][cindx(a,a+1)[n]],H))
            P.update({I:D[a-1][u]}) # guardo info de nuevos cuerpos (estrellas combinadas)
            
            for h in H:
                D[a][cindx(a,a+1)[h]]['del'] = 1 # marco para limpiar
        
    return D,P,I


def P_builder(j,v,r,m,A,P2,count,cells):
    """Updates the information of the dictionary P2.
    
    Keyword arguments:
    j -- Parameter identifier (integer)
    v -- velocity vector (list: [velocity in x, velocity in y, velocity in z])
    r -- position vector (list: [position in x, position in y, position in z])
    m -- mass (integer)
    A -- length (integer)
    P2 -- Dicctionary
    count -- counter (integer)
    cells -- list of positions in base-8
    """
    
    b = '0' + p2o(r[0],r[1],r[2],A)
    
    if b in cells.values():
        count -= 1
        
        i = list(cells.keys())[list(cells.values()).index(b)] # indice del primer planeta en esa celda

        P2[i] = {'X':r[0],'Y':r[1],'Z':r[2], 'VX':(v[0]*m +P2[i]['VX']*P2[i]['M'])/(m+P2[i]['M']),'VY':(v[1]*m +P2[i]['VY']*P2[i]['M'])/(m+P2[i]['M']),'VZ':(v[2]*m +P2[i]['VZ']*P2[i]['M'])/(m+P2[i]['M']),'M':P2[i]['M']+m,'cell':cells[i],'S':1,'I':i,'ancs':[i],'nodes':[i] }
        
    else:
        P2[j] = {'X':r[0],'Y':r[1],'Z':r[2],'VX':v[0],'VY':v[1],'VZ':v[2],'M':m,'cell':b,'S':1,'I':j,'ancs':[j],'nodes':[j]}
        cells[j] = '0' + p2o(r[0],r[1],r[2],A)
    
    return P2,count,cells


def tree_builder(P):
    """Builds a Tree.
    
    Keyword arguments:
    P -- Dictionary   
    """    
    def tree_builder_aux(T,n,level): # empieza construyendo el nodo root
        T['level'] = level
        T['nodes'] = P[str(n)]['nodes']
        ANCS = P[str(n)]['ancs']
        if len(ANCS) > 1 :
            level += 1
            for a in ANCS:
                tree_builder_aux(T[a],a,level)  # llamada recursiva
        return T
    ROOT = max(map(int,P.keys()))
    T = tree()					    # defino el arbol
    tree_builder_aux(T[ROOT],ROOT,0)	# lo relleno
    return T


def Tree_Renewer(P,N,A):
    """Updates a dictionary P and a tree.
    
    Keyword arguments:
    Q -- Tree   
    """
    T = tree()
    Q = tree()    
    
    Q[A] = {P[n]['cell']:{'X':P[n]['X'],'Y':P[n]['Y'],'Z':P[n]['Z'],'VX':P[n]['VX'],'VY':P[n]['VY'],'VZ':P[n]['VZ'],'M':P[n]['M'],'S':1,'I':n,'ancs':[n],'nodes':[n]} for n in range(N) if n in P.keys()}
    
    I = max(map(int,P.keys())) 
    for a in range(A-1,-1,-1): Q[a] = {}
    
    for a in range(A-1,0,-1):
        Q = promotion(Q,a)
        #Q = clean(Q)
        Q,P,I = combination(Q,a,P,I)
        #Q = clean(Q)
        
#    SPQR = deepcopy(Q)
#    for k,v in Q.items():
#        if len(v) == 0:
#            del SPQR[k]
#    Q = SPQR
    
    # el tree_builder tal y como está escrito funciona a partir del P leído del fichero save_run
    # es decir, que está serializado (pasado por json).

    # Alternativamente lo podemos serializar hic et nunc:
    
    P = eval(json.dumps(P, indent = 10))
    
    T = tree_builder(P)
            
    return T,P

def new_velocity_and_position(a,F,xmax,P):
    """Updates velocity and position using Euler algorithm.
    
    Keyword arguments:
    a -- Parameter identifier (integer)
    F -- acceleration vector (list: [acceeration in x, acceleration in y, acceleration in z])
    xmax -- size (integer)
    P -- Dictionary
    """
    r = [P[str(a)]["X"],P[str(a)]["Y"],P[str(a)]["Z"]]
    v = [P[str(a)]["VX"],P[str(a)]["VY"],P[str(a)]["VZ"]]
    #m = P[str(a)]["M"]
    tau = 0.1
    
    for i in range (3):  
                    
        v[i] = v[i] + F[i]*tau
        
        r[i] = round(r[i] + v[i]*tau)
        
        if (r[i] >= xmax):
            v[i] = -v[i]
            r[i] = xmax
        
        elif (r[i] <= 0):
            v[i] = -v[i]
            r[i] = 0
                       
    # hacer rebotar, cambio signo v  
    return v,r

def IvsSize (D,I,L,A):
    """Creates a list.
    
    Keyword arguments:
    D -- Dicctionary
    I -- Parameter identifier (integer)
    L -- List
    A -- length
    """
    IN = {}
    for el in list(D[I].keys())[2:]:
        L[el] = 2**(A-D[I][el]["level"])
        IN[el] = D[I][el]
        IvsSize (IN,el,L,A)
        IN = {}
    return L


cff=[]

def BH(D,I,j,F,P,L,N,theta):
    """Calculates gravitational acceleration exerted in body j.
    
    Keyword arguments:
    D -- Dicctionary
    I -- Parameter identifier (integer). Labels of bodies in the system.
    j -- Parameter identifier (integer). Calculate the acceleration exerted in this body.
    F -- acceleration vector (list: [acceeration in x, acceleration in y, acceleration in z])
    P -- Dictionary
    L -- list
    N -- number (integer)
    theta -- invariable parameter. BH threashold criterium.
    """
    G = 0.4
    IN = {}
    for el in list(D[I].keys())[2:]:
                
        IN[el] = D[I][el]
        dx = P[str(j)]["X"] - P[str(el)]["X"]
        dy = P[str(j)]["Y"] - P[str(el)]["Y"]
        dz = P[str(j)]["Z"] - P[str(el)]["Z"]
        d = math.sqrt(dx**2 + dy**2 + dz**2)
        s = L[el]
        if j == el: # no tenemos en cuenta la autofuerza. fx,fy = 0
            continue
                
        elif d == 0: # no tenemos en cuenta la autofuerza. Redundante o excluye otros casos? fx,fy = 0
            continue
            
        elif el < N or s/d < theta: # CALCULAMOS ACELERACION. Es un solo cuerpo o el cúmulo está lejos
            F[0] += -G*P[str(el)]["M"]*dx/d**3
            F[1] += -G*P[str(el)]["M"]*dy/d**3
            F[2] += -G*P[str(el)]["M"]*dz/d**3

            cff.append(1) # Contador de numero de veces que se calcula la fuerza            
            
        else:
            BH(IN,el,j,F,P,L,N,theta) # el es un nodo que está demasiado cerca

  
        IN = {}
    
    
    return F,cff
