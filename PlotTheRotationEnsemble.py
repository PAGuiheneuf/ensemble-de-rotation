# -*- coding: utf-8 -*-
"""
Plot des ensembles de rotation issus du code C++
Le script attend en entrée le fichier .txt contenant:
- premiere ligne: les valeurs de N, d, Tps, M
- seconde ligne: l'ensemble des valeurs de paramètres
- toutes les autres lignes: ensembles de rotations pour une valeur de paramètre
"""

import math
import os
import sys
import matplotlib.pyplot as plt



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

    lxTemp=[]
    lyTemp=[]

    for k in range(M):
        theta = float(math.pi*k)/(4*M)
        thetap = float(math.pi*(k+1))/(4*M)
        r=rot[k]
        rp=rot[k+1]
        lxTemp.append((-r*math.sin(thetap)+rp*math.sin(theta))/math.sin(theta-thetap))
        lyTemp.append((r*math.cos(thetap)-rp*math.cos(theta))/math.sin(theta-thetap))

    lx=lxTemp+list(reversed(lyTemp))
    ly=lyTemp+list(reversed(lxTemp))

    lx=lx+list(reversed([ -x for x in lx]))
    ly=ly+list(reversed(ly))

    lx=lx+[-x for x in lx]
    ly=ly+[-y for y in ly]

    aire=0
    for i in range(len(lx)-1):
        aire+=lx[i]*ly[i+1] - ly[i]*lx[i+1]

    aire+=lx[len(lx)-1]*ly[0] - ly[len(lx)-1]*lx[0]
    aire=aire/2


    plt.plot(lx, ly, 'ro-')
    plt.fill(lx, ly, '#fffe71')


    plt.axes().set_aspect('equal', 'datalim')
    plt.axis([-1.2,1.2,-1.2,1.2])
    s = "N"+repr(N)+"d"+repr(d)+"T"+repr(Tps)+"p"+repr(int(param)).zfill(3)

    plt.title(", ".join([ X[0] + "= " + str(X[1]).zfill(3) for X in [["N",N],["d",d],["T",Tps],["p",param/100],["Aire",aire]] ]))

    if save:
        plt.savefig('%s.png' % s)

    plt.close()
    print(rot)
    return aire



if __name__ == "__main__":
    #On vérifie que l'argument passé au code contienne bien le chemin vers le fichier ensembleRotation.txt
    if len(sys.argv) != 2:
        print("Too few or too many arguments to the script, please specify only the path to the C++ generated file")
        sys.exit(0)
    if not os.path.exists(sys.argv[1]):
        print("The file " + sys.argv[1] + " does not exist")
        sys.exit(0)

    #Lecture du fichier de résultats (rotationEnsemble.txt) généré par le code C++
    with open(sys.argv[1]) as f:
        lines = [l.strip() for l in f.readlines()]

        #Première ligne
        try:
            N   = int(lines[0].split()[0])
            d   = int(lines[0].split()[1])
            Tps = int(lines[0].split()[2])
            M   = int(lines[0].split()[3])
        except:
            print("Did not manage to parse the first line correctly!")
            sys.exit()

        #Seconde ligne (valeurs des paramètres)
        try:
            params = [float(x) for x in lines[1].split()]
        except:
            print("Did not manage to parse the second line correctly!")
            sys.exit()

        #Ensembles de rotation pour chaque paramètre:
        try:
            ROTS = [ [float(x) for x in l.split()] for l in lines[2:] ]
        except:
            print("Did not manage to parse the rotation ensembles correctly!")
            sys.exit()

        #Plot des images
        aire=[]
        for p, rot in zip(params, ROTS):
            a=plot_the_rotation_ensemble(rot, p*100, save=True)
            aire.append(a)

        print(aire)
        plt.close()
        plt.plot(aire)
        plt.savefig('aire.png')
        #Création de la vidéo avec ffmpeg
        if True:
            pattern = 'N%id%iT%ip' % (N, d, Tps)
            cmd = "ffmpeg -framerate 6 -i '%*.png' -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4"
            print(cmd)
            os.system(cmd)

# (bash) : pour sortir une vidéo à partir des images
# ffmpeg -framerate 30 -pattern_type glob -i '*.png' -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4









"""


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

    lxTemp=[]
    lyTemp=[]

    for k in range(M):
        theta = float(math.pi*k)/(4*M)
        thetap = float(math.pi*(k+1))/(4*M)
        r=rot[k]
        rp=rot[k+1]
        lxTemp.append((-r*math.sin(thetap)+rp*math.sin(theta))/math.sin(theta-thetap))
        lyTemp.append((r*math.cos(thetap)-rp*math.cos(theta))/math.sin(theta-thetap))

    lx=lxTemp+list(reversed(lyTemp))
    ly=lyTemp+list(reversed(lxTemp))

    lx=lx+list(reversed([ -x for x in lx]))
    ly=ly+list(reversed(ly))

    lx=lx+[-x for x in lx]
    ly=ly+[-y for y in ly]


    plt.plot(lx, ly, 'ro-')
    plt.fill(lx,ly, '#fffe71')


    plt.axes().set_aspect('equal', 'datalim')
    plt.axis([-1.2,1.2,-1.2,1.2])
    s = "N"+repr(N)+"d"+repr(d)+"T"+repr(Tps)+"p"+repr(int(param)).zfill(3)

    plt.title(", ".join([ X[0] + "= " + str(X[1]).zfill(3) for X in [["N",N],["d",d],["T",Tps],["p",param]] ]))

    if save:
        plt.savefig('%s.png' % s)

    plt.close()
    print(rot)



if __name__ == "__main__":
    #On vérifie que l'argument passé au code contienne bien le chemin vers le fichier ensembleRotation.txt
    if len(sys.argv) != 2:
        print "Too few or too many arguments to the script, please specify only the path to the C++ generated file"
        sys.exit(0)
    if not os.path.exists(sys.argv[1]):
        print "The file " + sys.argv[1] + " does not exist"
        sys.exit(0)

    #Lecture du fichier de résultats (rotationEnsemble.txt) généré par le code C++
    with open(sys.argv[1]) as f:
        lines = [l.strip() for l in f.readlines()]

        #Première ligne
        try:
            N   = int(lines[0].split()[0])
            d   = int(lines[0].split()[1])
            Tps = int(lines[0].split()[2])
            M   = int(lines[0].split()[3])
        except:
            print "Did not manage to parse the first line correctly!"
            sys.exit()

        #Seconde ligne (valeurs des paramètres)
        try:
            params = [float(x) for x in lines[1].split()]
        except:
            print "Did not manage to parse the second line correctly!"
            sys.exit()

        #Ensembles de rotation pour chaque paramètre:
        try:
            ROTS = [ [float(x) for x in l.split()] for l in lines[2:] ]
        except:
            print "Did not manage to parse the rotation ensembles correctly!"
            sys.exit()

        #Plot des images
        for p, rot in zip(params, ROTS):
            plot_the_rotation_ensemble(rot, p*100, save=True)

        #Création de la vidéo avec ffmpeg
        if True:
            pattern = 'N%id%iT%ip' % (N, d, Tps)
            cmd = "ffmpeg -framerate 6 -i '%*.png' -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4"
            print cmd
            os.system(cmd)

# (bash) : pour sortir une vidéo à partir des images
# ffmpeg -framerate 30 -pattern_type glob -i '*.png' -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4
"""
