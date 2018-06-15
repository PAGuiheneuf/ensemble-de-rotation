# -*- coding: utf-8 -*-
"""
Éditeur de Spyder

Ce script temporaire est sauvegardé ici :
/home/pag/.spyder2/.temp.py
"""

import math
import time
import matplotlib.pyplot as plt
import numpy
import sys
#Imports pour le parallélisme
from functools import partial # nécessaire si on veut passer des paramètres à
import multiprocessing


compteur  = 0
compteur2 = 0


###############################################################"

class Node:
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.edges = []
        self.edges2 = []

    def Add_edge(self, edge):
        self.edges.append(edge)

    def Add_edge2(self, edge2):
        self.edges2.append(edge2)

    def Get_name(self):
        return "n%u_%u" % (self.x, self.y)


class Edge:
    def __init__(self, startnode, endnode, edge2, etiquette = [1,0]):
        self.etiquette = etiquette
        self.start = startnode
        self.end = endnode
        self.start.Add_edge(self)
        self.edge2 = edge2

    def Get_etiquette_str(self):
        s = ""
        f = True
        for i in self.etiquette:
            if f:
                f = False
            else:
                s+= ";"
            s+= "%u" % i
        return s


class Edge2:
    def __init__(self, startnode, endnode, etiquette = [0]):
        self.etiquette = etiquette
        self.start = startnode
        self.end = endnode
        self.start.Add_edge2(self)

    def Get_etiquette_str(self):
        s = ""
        f = True
        for i in self.etiquette:
            if f:
                f = False
            else:
                s+= ";"
            s+= "%u" % i
        return s


class Graph:
    def __init__(self):
        self.nodes = []
        self.edges = []
        self.edges2 = []

    def Get_node(self, x, y):
        while len(self.nodes) <= x:
            self.nodes.append([])
        while len(self.nodes[x]) <= y:
            self.nodes[x].append(None)
        node = self.nodes[x][y]
        if node == None:
            self.nodes[x][y] = Node(x, y)
            node = self.nodes[x][y]
        return node

    def Add_edge(self, startx, starty, endx, endy, edgelink, etiquette):
        start = self.Get_node(startx, starty)
        end = self.Get_node(endx, endy)
        edge = Edge(start, end, edgelink, etiquette)
        self.edges.append(edge)

    def Add_edge2(self, startx, starty, endx, endy, etiquette):
        start = self.Get_node(startx, starty)
        end = self.Get_node(endx, endy)
        edge2 = Edge2(start, end, etiquette)
        self.edges2.append(edge2)
        return edge2

    def Export_dot(self):
        s = "digraph {\n"

        for e in self.edges:
            s += "    %s -> %s [label=\"%s\"];\n" % (e.start.Get_name(), e.end.Get_name(), e.Get_etiquette_str())
        s += "}\n"

        for e in self.edges2:
            s += "    %s -> %s [label=\"%s\"];\n" % (e.start.Get_name(), e.end.Get_name(), e.Get_etiquette_str())
        s += "}\n"

        return s


def Max_Etiq(graph, n):
    MaxE = False
    for e in graph.edges2:
        if e.etiquette[2] == n+1:
            if MaxE == False:
                MaxE = e.etiquette[1]
            else:
                MaxE = max(MaxE,e.etiquette[1])
    return MaxE


def Path_Bary(graph, n, t, u):
    current = graph
    maxE=0
    for e in current.edges2:
        if e.etiquette[2] >= n:
            if e.etiquette[2] == n:
                eti = e.etiquette[1]
            elif e.etiquette[2] == n+1:
                eti = e.etiquette[0]
            end = e.end
            for f in end.edges:
                e2 = f.edge2
                if e2.etiquette[2] == n+1:
                    e2.etiquette[1] = max(e2.etiquette[1] , eti + t*f.etiquette[0] + u*f.etiquette[1])
                    maxE=max(maxE,e2.etiquette[1])
                elif e2.etiquette[2] == n:
                    e2.etiquette[0] = e2.etiquette[1]
                    e2.etiquette[1] = eti + t*f.etiquette[0] + u*f.etiquette[1]
                    e2.etiquette[2] = n+1
                    maxE=max(maxE,e2.etiquette[1])

    return [current,maxE]


def map(coord,param): # La fonction à itérer
    x=coord[0] + param*math.sin(2*math.pi*coord[1])
    y=coord[1] + param*math.sin(2*math.pi*x)
    cp=[float(x),float(y)]
    return(cp)


def trans(g,x,y,N,d, param): #Construction du graphe de transition : coordonnées, nb de gros carrés, nb de subdivisions de ces gros carrés
    graphTrans=g
    for i in range(-d-1,d+1):  #Une boucle for pour chaque côté du carré
        j=-d-1
        temp=map([float(x+float(i)/(2*d))/N,float(y+float(j)/(2*d))/N], param)
        temp=[int(round(k*N)) for k in temp]
        temp2=[int(k%int(N)) for k in temp]
        temp[0] = int((temp[0]-temp2[0])/N)
        temp[1] = int((temp[1]-temp2[1])/N)

        start=graphTrans.Get_node(x,y)
        end=graphTrans.Get_node(temp2[0],temp2[1])

        edg2 = True
        for e in start.edges2:
            if e.end == end:
                edg2 = e
        if edg2 == True:
            edg2=graphTrans.Add_edge2(x,y,temp2[0],temp2[1],[0])

        edg = True
        for f in start.edges:
            if f.end == end:
                if f.etiquette == temp:
                    edg = False
                    global compteur
                    compteur+=1
        if edg == True:
            graphTrans.Add_edge(x,y,temp2[0],temp2[1],edg2,temp)

        global compteur2
        compteur2+=1




    for i in range(-d-1,d+1):
        j=d+1
        temp=map([float(x+float(i)/(2*d))/N,float(y+float(j)/(2*d))/N], param)
        temp=[int(round(k*N)) for k in temp]
        temp2=[int(k%int(N)) for k in temp]
        temp[0] = int((temp[0]-temp2[0])/N)
        temp[1] = int((temp[1]-temp2[1])/N)

        start=graphTrans.Get_node(x,y)
        end=graphTrans.Get_node(temp2[0],temp2[1])

        edg2 = True
        for e in start.edges2:
            if e.end == end:
                edg2 = e
        if edg2 == True:
            edg2=graphTrans.Add_edge2(x,y,temp2[0],temp2[1],[0])

        edg = True
        for f in start.edges:
            if f.end == end:
                if f.etiquette == temp:
                    edg = False
                    global compteur
                    compteur+=1
        if edg == True:
            graphTrans.Add_edge(x,y,temp2[0],temp2[1],edg2,temp)

        global compteur2
        compteur2+=1

    for j in range(-d-1,d+1):
        i=-d-1
        temp=map([float(x+float(i)/(2*d))/N,float(y+float(j)/(2*d))/N], param)
        temp=[int(round(k*N)) for k in temp]
        temp2=[int(k%int(N)) for k in temp]
        temp[0] = int((temp[0]-temp2[0])/N)
        temp[1] = int((temp[1]-temp2[1])/N)

        start=graphTrans.Get_node(x,y)
        end=graphTrans.Get_node(temp2[0],temp2[1])

        edg2 = True
        for e in start.edges2:
            if e.end == end:
                edg2 = e
        if edg2 == True:
            edg2=graphTrans.Add_edge2(x,y,temp2[0],temp2[1],[0])

        edg = True
        for f in start.edges:
            if f.end == end:
                if f.etiquette == temp:
                    edg = False
                    global compteur
                    compteur+=1
        if edg == True:
            graphTrans.Add_edge(x,y,temp2[0],temp2[1],edg2,temp)

        global compteur2
        compteur2+=1

    for j in range(-d-1,d+1):
        i=d+1
        temp=map([float(x+float(i)/(2*d))/N,float(y+float(j)/(2*d))/N], param)
        temp=[int(round(k*N)) for k in temp]
        temp2=[int(k%int(N)) for k in temp]
        temp[0] = int((temp[0]-temp2[0])/N)
        temp[1] = int((temp[1]-temp2[1])/N)

        start=graphTrans.Get_node(x,y)
        end=graphTrans.Get_node(temp2[0],temp2[1])
        edg2 = True
        for e in start.edges2:
            if e.end == end:
                edg2 = e
        if edg2 == True:
            edg2=graphTrans.Add_edge2(x,y,temp2[0],temp2[1],[0])
        edg = True
        for f in start.edges:
            if f.end == end:
                if f.etiquette == temp:
                    edg = False
                    global compteur
                    compteur+=1
        if edg == True:
            graphTrans.Add_edge(x,y,temp2[0],temp2[1],edg2,temp)

        global compteur2
        compteur2+=1


    return(graphTrans)



# "Wrappers" autour des fonctions
def compute_the_graph(param):
    g = Graph()
    for x in range(N):
        for y in range(N):
            g = trans(g,x,y,N,d, param)
    print "Graphe a %i noeuds, %i edges1 et %i edges2" % (len(g.nodes), len(g.edges), len(g.edges2) )
    return g
def compute_the_rotation_ensemble(graph):
    rot=[]
    #Pour plusieurs angles
    for k in range(M+1):
        for edg2 in graph.edges2:
            edg2.etiquette=[0,0,1]
        theta = float(math.pi*k)/(4*M)
        t=math.cos(theta)
        u=math.sin(theta)

        rotTemp=[]
        for n in range (1,Tps+1):
            [graph,maxE] = Path_Bary(graph, n, t, u)
            rotTemp.append(maxE)
        print rotTemp
        fit = numpy.polyfit(
            range(len(rotTemp)-30,len(rotTemp)+1),
            rotTemp[len(rotTemp)-31:len(rotTemp)],
            1
        )
        print fit
        rot.append(float(fit[0]))
        print(k)
    return rot
def plot_the_rotation_ensemble(rot, param, save=False):
    plt.close()
    plt.figure(figsize=(12,12))
    for k in range(M+1):
        theta = float(math.pi*k)/(4*M)
        t=math.cos(theta)
        u=math.sin(theta)
        r2=rot[k]

        plt.plot([r2*t+u,r2*t-u], [r2*u-t,r2*u+t], "b")
        plt.plot([r2*u-t,r2*u+t], [r2*t+u,r2*t-u], "b")
        plt.plot([-r2*t-u,-r2*t+u], [r2*u-t,r2*u+t], "b")
        plt.plot([r2*u-t,r2*u+t], [-r2*t-u,-r2*t+u], "b")
        plt.plot([r2*t+u,r2*t-u], [-r2*u+t,-r2*u-t], "b")
        plt.plot([-r2*u+t,-r2*u-t], [r2*t+u,r2*t-u], "b")
        plt.plot([-r2*u+t,-r2*u-t], [-r2*t-u,-r2*t+u], "b")
        plt.plot([-r2*t-u,-r2*t+u], [-r2*u+t,-r2*u-t], "b")

    plt.axes().set_aspect('equal', 'datalim')
    plt.axis([-1.2,1.2,-1.2,1.2])
    s = "N"+repr(N)+"d"+repr(d)+"T"+repr(Tps)+"p"+repr(param).zfill(3)

    plt.title(", ".join([ X[0] + "= " + str(X[1]).zfill(3) for X in [["N",N],["d",d],["T",Tps],["p",param]] ]))

    if save:
        plt.savefig('%s.png' % s)

#Fonction qui sera apellée par multiprocessing en //, avec p variable de 0 à 100
def do_the_whole_work(p):
    param=float(100-p)/100

    ##### Construction du graphe
    debut = time.time()
    g = compute_the_graph(param)
    print "compteurs = " + str(compteur) + ", " + str(compteur2)
    print('construction graphe : ' + repr(round(100*(time.time()-debut))/float(100)) + ' sec')

    ##### Exploration du graphe
    debut2=time.time()
    rot = compute_the_rotation_ensemble(g)
    print rot
    print('calcul rotation : ' + repr(round(100000*(time.time()-debut2))/float(100000)) + ' sec')

    duration = time.time() - debut

    ##### dessin de l'ensemble
    plot_the_rotation_ensemble(rot, p, save=True)

    return duration




#Paramètres d'itérations et de résolution
N   = 10   #Nb de gros carrés
d   = 16   #Nb de petits carrés dedans
Tps = 50  #Nb d'itérations
M   = 4   #Nb d'angles utilisés pour la représentation

if __name__ == "__main__":
    pool      = multiprocessing.Pool(1) # 15 = nb CPUS à utiliser (mon ordinateur a 16 coeurs)
    t0 = time.time()
    do_the_whole_work(0)
    #durations = pool.map(do_the_whole_work, range(1)) # range(100) = liste distribuée à la fonction do_the_whole_work
    print durations
    print "temps total = ", time.time() - t0

# (bash) : pour sortir une vidéo à partir des images
# ffmpeg -framerate 30 -pattern_type glob -i '*.png' -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4
