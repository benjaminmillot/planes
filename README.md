# Planes

Ce programme permet de calculer l'angle entre deux domaines
d'une même protéine. Il a été réalisé dans le cadre d'un projet long
pour le Master M2BI de l'Université Paris Diderot.

## Installation

Le programme a été testé fonctionnel avec la version 2.7 de Python.

Le package numpy pour python doit etre installé. Veuillez vous rendre sur
http://www.numpy.org/ pour télécharger la version correspondante a votre systeme
d'exploitation. Sur Ubuntu, il peut etre installé par :

    sudo apt-get install python-numpy

Le package matplotlib de génération de graphes pour python doit etre installé.
Veuillez vous rendre sur http://matplotlib.org/ pour télécharger la version
correspondante a votre distribution. Sur Ubuntu, il peut etre installé en ligne
de commande par :

    sudo apt-get install python-matplotlib

Le package GromacsTrajectory (Konrad Hinsen) doit être installé pour
assurer le bon fonctionnnement du programme.

    cd GromacsTrajectoyReader-0.12/
    python setup.py build
    python setup.py install     (peut necessiter les droits d'administrateurs)

Veuillez vous référer au fichier README present dans le dossier si besoin.

Certaines dépendances python doivent être installés et fonctionnelles, mais
sont usuellement installées de base sur la plupart des distribution: re,
argparse, matplotlib.

## Données d'entrée

Le programme prends en entrée des résultats de dynamique moléculaire:
    - un fichier topologie (.gro)
    - un fichier binaire de trajectoire (.xtc)

Le projet et rapport ont été réalisé dans le cadre de l'étude des deux domaines
S et P de la molécule A de 1IHM, constituant de la capside du Norovirus.

Ces deux fichiers participant à un projet de recherche encore non publié à ce
jour, et le fichier de la dynamique étant très volumineux, ils ne seront pas
présents dans cette archive. Veuillez prendre contact avec M. Thibault Tubiana
(thibault.tubiana@gmail.com) afin d'obtenir son autorisation et ses données.

## Lancement

Exemple de ligne de commande démarrant le programme:
    
    cd bin/
    python module.py -top ../input/VP1.gro -traj ../input/traj.xtc \
    -fd 1058 1230 -sd 1259 1549 -o ../log/log.txt

Les arguments d'entrée sont les suivants:

- -top [chaine de caractères] : chemin d'accès du fichier topologie
- -traj [chaine de caractères] : chemin d'accès du fichier de trajectoire
- -fd [deux entiers] : numéro des résidus du début et de fin du premier
                       domaine
- -sd [deux entiers] : numéro des résidus du début et de fin du second
                       domaine
- -o [optionnel, chaine de caractères] : redirige la sortie standard vers
                                         un fichier log de sortie

Note : les numéros qui délimitent les résidus à considérer doivent être
cohérents avec le fichier de topologie fourni.

## Description des dossier

- bin : scripts impliqués dans le programme.
- doc : contient le sujet, les instructions pour le rapport et la version finale
        du compte-rendu.
- graphs : dossier ou les graphes résultats seront produits. Des résultats sont
           fournis en exemple.
- GromacsTrajectoryReader-0.12 : dépendance à installer pour lire la trajectoire
- input [non-elementaire] : dossier ou l'utilisateur à la possibilité de placer
                            sa topologie et sa trajectoire. Un fichier topologie
                            exemple est fourni.
- log [non-elementaire] : dossierou l'utilisateur peut stocker ses fichiers log.
                          Un log de démonstration est founi.

## Contact

Auteur : Benjamin Millot

Email : benjamin.millot@outlook.fr
