from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP

from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import numpy as np
from scipy.stats import uniform_direction

class Residu:
    def __init__(self, numero, acide_amine, coord,structure_secondaire=None, surface_accessible=None):
        self.numero = numero
        self.acide_amine = acide_amine
        self.coord = coord
        self.structure_secondaire = structure_secondaire
        self.surface_accessible = surface_accessible
        

# Chemin Linux sous WSL
pdb_file = "/mnt/c/Users/nhoun/Desktop/projetcourt/1CRN.pdb"

def surface_accessible(file_path):
    """fonction qui calcule la surface accessible d'un proteine a partir d'un fichier pdb
    parametre:
    file_path: string chemin du fichier pdb
    return: dictionnaire avec le numero de residu et l'acide amine en clé et la surface accessible en valeur"""
    parser=PDBParser(QUIET=True)
    structure=parser.get_structure("protein", file_path)
    model=structure[0]
    dssp=DSSP(model, file_path, dssp="mkdssp")
    dico_AS={}
    for key in list(dssp.keys()):
        acide_amine=dssp[key][1]
        numero=dssp[key][0]
        if numero not in dico_AS:
            dico_AS[numero]={}
        dico_AS[numero][acide_amine]=dssp[key][3]
    return dico_AS

def structure_secondaire(file_path):
    """Fonction qui détermine la structure secondaire pour chaque résidu d'une protéine
    à partir d'un fichier PDB en utilisant DSSP.
    parametre:
    file_path : string chemin du fichier PDB
    Returns: dictionnaire avec le numéro de résidu en clé et la structure secondaire en valeur
    """
    parser=PDBParser(QUIET=True)
    structure=parser.get_structure("protein", file_path)
    model=structure[0]
    dssp=DSSP(model, file_path, dssp="mkdssp")
    dico_struct={}
    for key in list(dssp.keys()):
        structure_secondaire=dssp[key][2]
        numero=dssp[key][0]
        dico_struct[numero]={structure_secondaire}
    return dico_struct

#print(structure_secondaire(pdb_file))
    
def read_pdb(file_path):
    """
    Lit un fichier PDB et extrait les coordonnées des CA (carbone alpha) des résidus.
    parametre:
    file_path : string chemin du fichier PDB
    return: dictionnaire avec le numéro de résidu en clé et les coordonnées (x,y,z) des CA en valeur
    """
    CA_coord={}
    with open(file_path, 'r') as fichier:
        for ligne in fichier:
            if ligne.startswith("ATOM") and ligne[12:16].strip() == "CA":
                res_num =int(ligne[22:26].strip())
                x = float(ligne[30:38].strip())
                y = float(ligne[38:46].strip())
                z = float(ligne[46:54].strip())
                if res_num not in CA_coord:
                    CA_coord[res_num]={}
                CA_coord[res_num]= (x, y, z)
    return CA_coord


def invert_CA_dict(CA_coord):
    """
    Inverse le dictionnaire : coordonnées (x,y,z) comme clé,
    numéro de résidu comme valeur.
    parametre:
    CA_coord : dict dictionnaire avec le numéro de résidu en clé et les coordonnées (x,y,z) en valeur
    return: dictionnaire avec les coordonnées (x,y,z) des CA en clé et le numéro de résidu en valeur
    """
    coord_to_res = {}
    for res_num, coord in CA_coord.items():
        coord_to_res[coord] = res_num
    return coord_to_res

coordonnes_CA=read_pdb(pdb_file)
dico_inverse=invert_CA_dict(coordonnes_CA)

def calc_centremasse(coords):
    """
    Calcule le centre de masse d'un ensemble de coordonnées.
    parametre:
    coords : dict dictionnaire avec le numéro de résidu en clé et les coordonnées (x,y,z) des CA en valeur
    return: tuple (x,y,z) des coordonnées du centre de masse
    """
    coord_list=[]
    for coord in coords.values():
            coord_list.append(coord)
    coord_array=np.array(coord_list)
    #print(coord_array.shape)
    centre_masse=coord_array.mean(axis=0)
    return centre_masse 
   
centre=calc_centremasse(coordonnes_CA)

    
  
#print("centre :", centre) 
def generate_directions(center, n_theta=40):
    """Genere des directions uniformement reparties sur une sphere ici le centre de masse de la proteine
    parametre:
    center : tuple (x,y,z) des coordonnées du centre de masse
    n_theta : int nombre de directions a generer
    

    Returns:
        _type_: _description_
    """
    n_phi = 2 * n_theta
    center = np.array(center)
    lines = []

    costheta = np.linspace(-1, 1, n_theta)
    theta = np.arccos(costheta)
    phi = np.linspace(0, 2*np.pi, n_phi, endpoint=False)

    for t in theta:
        for p in phi:
            x = np.sin(t) * np.cos(p)
            y = np.sin(t) * np.sin(p)
            z = np.cos(t)
            dir_vec = np.array([x, y, z])
            point1 = center - dir_vec
            point2 = center + dir_vec
            lines.append((point1, point2))
    return lines


    
def numeroresidu_tranche(tranche, dico_coord, tol=1):
    """
    Retourne les numéros des résidus dont le CA est dans la tranche.
    """
    point1, point2 = tranche
    dir_vec = np.array(point2, dtype=float) - np.array(point1, dtype=float)
    dir_vec = dir_vec / np.linalg.norm(dir_vec)  # vecteur normalisé

    Liste_residu = []
    for coords, res_num in dico_coord.items():
        P = np.array(coords, dtype=float)
        # distance point-plan
        dist = np.abs(np.dot(P - point1, dir_vec))
        if dist <= tol:
            Liste_residu.append(res_num)
    return Liste_residu

directions=generate_directions(centre)
#print(numeroresidu_tranche(directions[0],dico_inverse))

#for direction in directions:
    #liste_residu=numeroresidu_tranche(direction,dico_inverse)
    #print(liste_residu)
def calc_hydrophobicity(residue_numbers, dico_AS):
    """
    Calcule l'hydrophobicité totale pour une liste de numéros de résidus.
    """
    aa_hydrophobe = ["F", "G", "I", "L", "M", "V", "W", "Y"]
    hydrophobic_surface = 0
    surface_totale = 0
    for res_num in residue_numbers:
        for aa, surface in dico_AS[res_num].items(): 
            surface_totale += surface
            if aa in aa_hydrophobe:
                hydrophobic_surface += surface
    if surface_totale == 0:
        return 0
    else:
        hydrophobic_factor= hydrophobic_surface / surface_totale
    return hydrophobic_factor


def deplacer_droite(droite, distance=1.0):
    """
    Déplace une droite de 'distance' Å dans la direction normale (orthogonale à la droite).
    La droite est définie par 2 points (point1, point2).
    """
    point1, point2 = droite
    point1 = np.array(point1, dtype=float)
    point2 = np.array(point2, dtype=float)

    # Vecteur directeur de la droite
    vec = point2 - point1
    vec = vec / np.linalg.norm(vec)  # normalisation

    # Trouver un vecteur normal à vec
    ref = np.array([1, 0, 0])
    if np.allclose(vec, ref):  # si parallèle, on change de référence
        ref = np.array([0, 1, 0])
    normal = np.cross(vec, ref)
    normal = normal / np.linalg.norm(normal)

    # Translation de la droite le long de la normale
    new_point1 = point1 + distance * normal
    new_point2 = point2 + distance * normal

    return (new_point1, new_point2)
def tranches_deplacees(droites, ntranche, step=1.0):
    """
    Prend une liste de droites et renvoie une liste de listes de tranches.
    Chaque sous-liste contient les ntranche droites décalées successivement
    de step Å le long de la normale.
    """
    toutes_tranches = []
    for droite in droites:
        tranches = [droite]
        current = droite
        for i in range(1, ntranche):
            current = deplacer_droite(current, step)
            tranches.append(current)
        toutes_tranches.append(tranches)
    return toutes_tranches


tranches=tranches_deplacees(directions,5)
#print(tranches[0])

def calc_qvalue(droites,ntranche,dico_coord,dico_AS):
    liste_Qvaleurs = []
    tranches = tranches_deplacees(droites, ntranche)
    for tranche in tranches:
        L=[]
        for droite in tranche:
            liste_residu = numeroresidu_tranche(droite, dico_coord)
            if len(liste_residu) > 0:
                hydrophobicity = calc_hydrophobicity(liste_residu, dico_AS)
                L.append(hydrophobicity)
        liste_Qvaleurs.append(L)
    L_Qvalues=[]
    for Qval in liste_Qvaleurs:
        Qval=np.array(Qval)
        L_Qvalues.append(np.nanmean(Qval))
    return L_Qvalues


    
a=calc_qvalue(directions,5,dico_inverse,surface_accessible(pdb_file))
print(a[0])
#dico_AS=surface_accessible(pdb_file)
liste_residu=[1,2,3,4,5,6,7,8,9,10]
#print(calc_hydrophobicity(liste_residu,dico_AS))
