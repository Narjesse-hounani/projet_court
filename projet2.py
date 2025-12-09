import numpy as np
import math
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP

class Residu:
    """Classe représentant un résidu d'une protéine.
    Attributs:
    numero : int numéro du résidu dans la chaîne
    acide_amine : str code à une lettre de l'acide aminé
    coord : tuple (x, y, z) coordonnées 3D de l'atome CA
    structure_secondaire : str structure secondaire 
    surface_accessible : float surface accessible 
    """
    def __init__(self, numero, acide_amine, coord,
                 structure_secondaire=None, surface_accessible=None):
        self.numero = numero
        self.acide_amine = acide_amine
        self.coord = coord
        self.structure_secondaire = structure_secondaire
        self.surface_accessible = surface_accessible

    def __repr__(self):
        return (f"Residu({self.numero}, {self.acide_amine}, "
                f"coord={self.coord}, "
                f"SS={self.structure_secondaire}, "
                f"ASA={self.surface_accessible})")


def lire_residus(file_path):
    """
    Construit une liste d'objets Residu à partir d'un fichier PDB
    en utilisant DSSP pour les surfaces accessibles et structures secondaires.
    parametre:
    file_path : chemin vers le fichier PDB
    return: liste d'objets Residu
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", file_path)
    model = structure[0]
    # Calcul DSSP
    dssp = DSSP(model, file_path, dssp="mkdssp")
    residus = []
    for key in dssp.keys():
        # chain_id identifiant de chaîne
        chain_id = key[0]
        # res_id identifiant complet du résidu
        res_id = key[1]       
        # aa acide aminé, ss structure secondaire, asa surface accessible
        aa = dssp[key][1]      
        ss = dssp[key][2]      
        asa = dssp[key][3]     
        # Vérification que le résidu existe bien dans la chaîne
        if chain_id not in model:
            continue
        chain = model[chain_id]
        if res_id not in chain:
            continue
        residue = chain[res_id]
        if "CA" in residue:  
            coord = tuple(residue["CA"].coord)
            r = Residu(
                numero=res_id[1],   # numéro simple
                acide_amine=aa,
                coord=coord,
                structure_secondaire=ss,
                surface_accessible=asa)
            residus.append(r)
    return residus


def calc_centremasse(liste_residu):
    """
    Calcule le centre de masse d'un ensemble de coordonnées.
    parametre:
    coords : dict dictionnaire avec le numéro de résidu en clé et les coordonnées (x,y,z) des CA en valeur
    return: tuple (x,y,z) des coordonnées du centre de masse
    """
    coord_list = []
    # extraction des coordonnées
    for residu in liste_residu:
        if residu.coord is not None:
            coord_list.append(residu.coord)
    # calcul du centre de masse
    coord_array = np.array(coord_list)
    centre_masse = coord_array.mean(axis=0)
    return centre_masse 

def generate_directions(center, n_theta=10):
    """
    Génère des directions uniformément réparties sur une sphère autour du centre de masse d'une protéine.

    Paramètres
    ----------
    center : tuple (x, y, z)
        Coordonnées du centre de masse.
    n_theta : int, optionnel (par défaut = 40)
        Nombre de divisions pour l’angle polaire (contrôle le nombre de directions).

    Retour
    ------
    lines : list of tuples
        Chaque élément est (point1, point2), deux extrémités d’une droite passant par 'center'.
    """
    # Nombre de divisions en phi (longitude)
    n_phi = 2 * n_theta
    # Conversion du centre en numpy array
    center = np.array(center)
    # Liste pour stocker les droites
    lines = []
    # Échantillonnage uniforme en cos(theta) pour une bonne répartition des points sur la sphère
    costheta = np.linspace(-1, 1, n_theta)
    theta = np.arccos(costheta)  # angle polaire
    phi = np.linspace(0, 2*np.pi, n_phi, endpoint=False)  # angle azimutal
    # Boucle sur toutes les combinaisons de (theta, phi)
    for t in theta:
        for p in phi:
            # Conversion coordonnées sphériques → cartésiennes (vecteur unitaire)
            x = np.sin(t) * np.cos(p)
            y = np.sin(t) * np.sin(p)
            z = np.cos(t)
            dir_vec = np.array([x, y, z])
            # Construction des deux extrémités de la droite (symétriques autour du centre)
            point1 = center - dir_vec
            point2 = center + dir_vec
            # Ajout de la droite dans la liste
            lines.append((point1, point2))
    return lines

def construire_parallelepipede(pt1, pt2, width, height):
    """
    Construit les sommets d'un parallélépipède à partir de deux points en 3D, largeur et hauteur.
    pt1, pt2 : tuples (x, y, z)
    width : largeur 
    height : hauteur
    Retourne la liste des 8 sommets du parallélépipède.
    """
    pt1 = np.array(pt1)
    pt2 = np.array(pt2)

    # Vecteur de la base (de pt1 à pt2)
    base_vec = pt2 - pt1
    base_length = np.linalg.norm(base_vec)

    # base_dir vecteur unitaire base
    base_dir = base_vec / base_length 

    # Pour trouver une direction perpendiculaire à base_dir,
    # on choisit arbitrairement un vecteur non colinéaire
    arbitrary_vec = np.array([0, 0, 1]) if abs(base_dir[2]) < 0.9 else np.array([0, 1, 0])

    # Vecteur largeur perpendiculaire à base_dir
    width_dir = np.cross(base_dir, arbitrary_vec)
    width_dir = width_dir / np.linalg.norm(width_dir)

    # Vecteur hauteur (perpendiculaire aux deux autres)
    height_dir = np.cross(base_dir, width_dir)
    height_dir = height_dir / np.linalg.norm(height_dir)

    # Construire les 8 sommets
    base_point_1 = pt1
    base_point_2 = pt2
    base_point_3 = pt2 + width_dir * width
    base_point_4 = pt1 + width_dir * width

    # Sommets du haut
    top_point_1 = base_point_1 + height_dir * height
    top_point_2 = base_point_2 + height_dir * height
    top_point_3 = base_point_3 + height_dir * height
    top_point_4 = base_point_4 + height_dir * height

    # Liste des sommets
    sommets = [
        base_point_1, base_point_2, base_point_3, base_point_4,
        top_point_1, top_point_2, top_point_3, top_point_4
    ]
    return [tuple(pt) for pt in sommets]


def point_dans_parallelepipede(sommets, point):
    """
    Détermine si un point est dans un parallélépipède défini par ses 8 sommets.
    
    sommets : liste de 8 tuples (x, y, z)
    point : tuple (x, y, z)
    
    Retourne True si point est dans le parallélépipède, False sinon.
    """    
    pts = []
    for s in sommets:
        pts.append(np.array(s))
    P = np.array(point)
    # Origine O
    O = pts[0]
    # Vecteurs formant le parallélépipède (arêtes)
    u = pts[1] - O  # base longueur
    v = pts[3] - O  # base largeur
    w = pts[4] - O  # hauteur
    # Résoudre le système: P - O = a*u + b*v + c*w
    vec = P - O
    M = np.column_stack((u, v, w))    
    # Vérifie si la matrice est inversible (rang 3)
    if np.linalg.matrix_rank(M) < 3:
        return False
    a, b, c = np.linalg.solve(M, vec)
    # Vérifier si a, b, c sont dans [0, 1]
    if (0 <= a <= 1) and (0 <= b <= 1) and (0 <= c <= 1):
        return True
    else:
        return False

def residu_tranche(tranche, width, height, liste_residu):
    """
    Retourne les résidus dont le CA est à l’intérieur du parallélépipède
    défini par les deux points de 'tranche', une largeur et une hauteur.
    Args:
        tranche (tuple): deux points (point1, point2) définissant la base du parallélépipède
        width (float): largeur du parallélépipède
        height (float): hauteur du parallélépipède
        liste_residu (list[Residu]): liste d'objets Residu avec attributs `coord` et `numero`

    Returns:
        list[int]: numéros des résidus dans le parallélépipède
    """
    sommets = construire_parallelepipede(tranche[0], tranche[1], width, height)
    liste_residu_tranche= []
    for residu in liste_residu:
        if point_dans_parallelepipede(sommets, residu.coord):
            liste_residu_tranche.append(residu)
    return set(liste_residu_tranche)

def bounding_box_coords(liste_residus):
    """
    Retourne les coordonnées min et max de la protéine (boîte englobante cubique).
    parametre: liste_residus : liste d'objets Residu 
    Returns: (cube_min, cube_max) : tuples des coordonnées min et max
    """
    coords = []
    for r in liste_residus:
        if r.coord is not None:
            coords.append(r.coord)
    if len(coords) == 0:
        return None, None
    coords = np.array(coords, dtype=float)
    # min et max sur chaque axe
    min_vals = coords.min(axis=0)
    max_vals = coords.max(axis=0)
    # cube : côté = plus grande dimension
    lengths = max_vals - min_vals
    max_len = lengths.max()
    # centre du cube
    center = (min_vals + max_vals) / 2.0
    cube_min = center - max_len / 2.0
    cube_max = center + max_len / 2.0
    return cube_min, cube_max


def max_tranches_from_axis(liste_residus, droite, epaisseur_min=1.0):
    """
    Calcule le nombre max de tranches possibles dans la direction normale à une droite,
    en utilisant la boîte englobante de la protéine.
    Args:
        liste_residus (list[Residu]): liste d'objets avec attribut coord
        droite (tuple): deux points (p1, p2) définissant l’axe
        epaisseur_min (float): épaisseur minimale en Å
    Returns:
        int: nombre de tranches possibles
    """
    cube_min, cube_max = bounding_box_coords(liste_residus)
    if cube_min is None:
        return 0
    # vecteur directeur de la droite
    p1=np.array(droite[0])
    p2=np.array(droite[1])
    axis_vec = p2 - p1
    axis_vec = axis_vec / np.linalg.norm(axis_vec)
    # vecteur normal (perpendiculaire à axis_vec)
    arbitrary = np.array([1, 0, 0])
    if np.allclose(np.cross(axis_vec, arbitrary), 0):
        arbitrary = np.array([0, 1, 0])
    normal = np.cross(axis_vec, arbitrary)
    normal = normal / np.linalg.norm(normal)
    # générer les 8 coins du cube
    corners = []
    for x in [cube_min[0], cube_max[0]]:
        for y in [cube_min[1], cube_max[1]]:
            for z in [cube_min[2], cube_max[2]]:
                corners.append(np.array([x, y, z]))
    # projeter les coins sur le vecteur normal
    projections = [np.dot(corner, normal) for corner in corners]
    min_proj = min(projections)
    max_proj = max(projections)
    longueur = max_proj - min_proj
    # nombre de tranches
    n_tranches = math.ceil(longueur / float(epaisseur_min))
    return n_tranches


def calc_hydrophobicity(liste_residu):
    """
    Calcule l'hydrophobicité totale pour une liste de numéros de résidus
    à partir de la liste d'objets Residu.
    
    Args:
    liste_residu (list[Residu]): liste d'objets Residu
    
    Returns:
        float: facteur d'hydrophobicité (0 à 1)
    """
    # aa_hydrophobe : liste des acides aminés hydrophobes
    aa_hydrophobe = ["F", "G", "I", "L", "M", "V", "W", "Y"]
    hydrophobic_surface = 0.0
    surface_totale = 0.0
    for residu in liste_residu:
        surface_totale += residu.surface_accessible
        if residu.acide_amine in aa_hydrophobe:
            hydrophobic_surface += residu.surface_accessible
    if surface_totale == 0:
        return 0.0
    else:
        return hydrophobic_surface / surface_totale


def deplacer_droite(droite, droite_reference, distance=1.0):
    """
    Déplace une droite de 'distance'  dans la direction normale à une droite de référence.
    
    Args:
        droite (tuple): droite à déplacer, définie par (point1, point2)
        droite_reference (tuple): droite fixe pour définir la direction normale
        distance (float): distance de déplacement en Angströms
    
    Returns:
        tuple: nouvelle droite déplacée (new_point1, new_point2)
    """
    point1, point2 = droite
    point1 = np.array(point1, dtype=float)
    point2 = np.array(point2, dtype=float)
    # Vecteur directeur de la droite de référence
    ref_vec = np.array(droite_reference[1], dtype=float) - np.array(droite_reference[0], dtype=float)
    ref_vec = ref_vec / np.linalg.norm(ref_vec)
    # Vecteur directeur de la droite à déplacer
    vec = point2 - point1
    vec = vec / np.linalg.norm(vec)
    # Vecteur normal à la droite de référence et à la droite actuelle
    normal = np.cross(vec, ref_vec)
    norm_len = np.linalg.norm(normal)
    #1e-8 tolérance pour les flottants
    if norm_len < 1e-8:
        # Si vec et ref_vec sont colinéaires, on choisit un vecteur arbitraire perpendiculaire
        normal = np.cross(vec, [1, 0, 0])
        if np.linalg.norm(normal) < 1e-8:
            normal = np.cross(vec, [0, 1, 0])
        normal = normal / np.linalg.norm(normal)
    else:
        normal = normal / norm_len
    # Translation de la droite le long de la normale
    new_point1 = point1 + distance * normal
    new_point2 = point2 + distance * normal
    return (new_point1, new_point2)

def tranches_deplacees(droites, ntranche, step=1.0):
    """
    Prend une liste de droites et renvoie une liste de listes de tranches.
    Chaque sous-liste contient les ntranche droites décalées successivement
    de step Å le long de la normale à la droite originale.

    Args:
        droites (list[tuple]): liste de droites à déplacer, chaque droite = (point1, point2)
        ntranche (int): nombre de tranches à générer
        step (float): distance entre chaque tranche (en Å)
    
    Returns:
        list[list[tuple]]: liste de listes de droites déplacées
    """
    toutes_tranches = []
    for droite in droites:
        tranches = [droite]  # première tranche = droite originale
        current = droite
        for i in range(1, ntranche):
            # On passe la droite originale comme référence pour que la normale reste constante
            current = deplacer_droite(current, droite_reference=droite, distance=step)
            tranches.append(current)
        toutes_tranches.append(tranches)
    return toutes_tranches

def filtrage_tranches(tranches, liste_residu, width, height):
    """
    Filtre les tranches pour ne garder que celles contenant au moins un résidu.
    
    Args:
        tranches (list[tuple]): liste de liste de droites (tranches)
        liste_residu (list[Residu]): liste d'objets Residu
        width (float): largeur du parallélépipède
        height (float): hauteur du parallélépipède
    
    Returns:
        list[tuple]: tranches filtrées
    """
    tranches_filtrees = []
    for tranche in tranches:
        for droite in tranche:
            residue_numbers = residu_tranche(tranche, width, height, liste_residu)
        if len(residue_numbers) > 0:
            tranches_filtrees.append(tranche)
    return tranches_filtrees

def calc_qvalue(droites, ntranche, liste_residu):
    """
    Calcule la Q-value pour chaque tranche de droites décalées.
    
    Args:
        droites (list[tuple]): liste de droites initiales
        ntranche (int): nombre de tranches à générer
        liste_residu (list[Residu]): liste d'objets Residu    
    
    Returns:
        tuple: (Q-value maximale, tranche correspondante)
    """
    #Liste de liste des Q-values pour chaque axe
    liste_Qvaleurs = []
    tranches = tranches_deplacees(droites, ntranche)
    for tranche in tranches:
        # Liste des Q-values pour cette tranche
        L = []
        for droite in tranche:
            residue_numbers = residu_tranche(droite, 10, 10, liste_residu)
            hydrophobicity = calc_hydrophobicity(residue_numbers)
            L.append(hydrophobicity)
        liste_Qvaleurs.append(L)
    # Calcul de la moyenne des Q-values pour chaque axe
    L_Qvalues = []
    if len(liste_Qvaleurs) == 0:
        return (0.0, None)
    # liste_Qvaleurs est une liste de listes
    for Qval in liste_Qvaleurs:
        if Qval is not None and len(Qval) > 0:
            mean_val = np.nanmean(np.array(Qval))
            L_Qvalues.append(mean_val)
        else:
            L_Qvalues.append(0.0)
    # Trouver la Q-value maximale et l'axe correspondant
    id_axe = np.argmax(L_Qvalues)
    axe_membrane = tranches[id_axe]
    return (max(L_Qvalues), axe_membrane)



def position_membrane(axe, width_init, height_init, liste_residus, n_tranches=5, max_epaisseur=20):
    """
    Trouve la position de la membrane en cherchant la droite 
    avec l'hydrophobicité maximale, en testant différentes épaisseurs.

    Args:
        axe (list[tuple]): axe de droites (comme renvoyé par calc_qvalue)
        width_init (float): largeur de base des parallélépipèdes
        height_init (float): hauteur de base des parallélépipèdes
        liste_residus (list[Residu]): liste des résidus
        n_tranches (int): nombre de tranches à générer le long de l’axe
        max_epaisseur (int): nombre de pas à tester pour augmenter width et height

    Returns:
        dict: informations sur la membrane optimale :
            {
                "residus": liste des résidus sélectionnés,
                "hydrophobicite": valeur max,
                "tranche_id": index de la tranche optimale,
                "epaisseur": valeur de l’augmentation (i),
                "droite": la meilleure droite trouvée
            }
    """
    tranches = tranches_deplacees(axe, n_tranches)
    meilleure_hydro = 0.0
    meilleurs_residus = []
    meilleure_epaisseur = 0
    meilleure_tranche_id = -1
    meilleure_droite = None
    # Parcours des tranches
    for tranche_id, tranche in enumerate(tranches):
        # tester différentes épaisseurs
        for i in range(max_epaisseur):
            width = width_init + i
            height = height_init + i
            for droite in tranche:
                # résidus capturés par cette droite
                nums = residu_tranche(droite, width, height, liste_residus)
                if len(nums) > 0:
                    hydro = calc_hydrophobicity(nums)
                    # mise à jour globale si cette droite est meilleure
                    if hydro > meilleure_hydro:
                        meilleure_hydro = hydro
                        meilleurs_residus = nums
                        meilleure_epaisseur = i
                        meilleure_tranche_id = tranche_id
                        meilleure_droite = droite
    dico_position = {
        "residus": meilleurs_residus,
        "hydrophobicite": meilleure_hydro,
        "tranche_id": meilleure_tranche_id,
        "epaisseur": meilleure_epaisseur,
        "droite": meilleure_droite
    }
    return dico_position

def classification_prot(droites, ntranche, liste_residu):
    """
    Classifie une protéine en soluble ou membranaire en fonction de la Q-value maximale.
    
    Args:
        droites (list[tuple]): liste de droites initiales
        ntranche (int): nombre de tranches à générer
        liste_residu (list[Residu]): liste d'objets Residu    
    
    Returns:
        tuple: (classification, Q-value maximale, tranche correspondante)
    """
    q_value, axe_membrane = calc_qvalue(droites, ntranche, liste_residu)
    classification="soluble"
    if q_value>=0.45:
        classification="transmembranaire"
        resultat_position_membrane = position_membrane(axe_membrane,10,10,liste_residu,ntranche)
        residu_membrane=[]
        for residu in resultat_position_membrane["residus"]:
            residu_membrane.append(residu)
        structure_secondaire_membrane = []
        for residu in residu_membrane:
            structure_secondaire_membrane.append(residu.structure_secondaire)
        # set pour extraire les structures secondaires uniques
        structure_secondaire= set(structure_secondaire_membrane)
        # Comptage de la fréquence de chaque type de structure
        dico_structure = {}
        for ss in structure_secondaire_membrane:
            dico_structure[ss] = structure_secondaire_membrane.count(ss)
        type_membrane=max(dico_structure, key=dico_structure.get)
        return (classification, resultat_position_membrane,type_membrane)
    # Si la protéine est soluble : retour classification + valeur Q
    return (classification, q_value)

# Exemple d’utilisation
if __name__ == "__main__":
    pdb_file = input("Entrez le chemin du fichier PDB : ").strip()
    # 1. Lecture des résidus depuis un fichier PDB
    liste_residus = lire_residus(pdb_file)
    print("Exemple de résidus lus :", liste_residus[:5])  # Affiche les 5 premiers résidus pour vérifier
    # 2. Calcul du centre de masse de la protéine
    centre = calc_centremasse(liste_residus)
    print("Centre de masse de la protéine :", centre)
    # 3. Génération de directions à partir du centre de masse
    directions = generate_directions(centre)
    # 4. Calcul du nombre maximum de tranches possibles selon chaque direction
    L_nombre_tranches = []
    for droite in directions:
        n_tranches = max_tranches_from_axis(liste_residus, droite, 1.0)
        L_nombre_tranches.append(n_tranches)
    max_tranches = max(L_nombre_tranches)
    print("Nombre max de tranches possibles :", max_tranches)
    # 5. Filtrage des directions pour ne garder que les tranches pertinentes
    droites_filtrees = filtrage_tranches(directions, liste_residus, 10, 10)
    # 6. Calcul de la Q-value sur les droites filtrées
    qvalue, axe = calc_qvalue(droites_filtrees, max_tranches, liste_residus)
    print("Q-value maximale :", qvalue)
    # 7. Classification finale de la protéine
    class_prot = classification_prot(droites_filtrees, max_tranches, liste_residus)
    print("Classification de la protéine :", class_prot)


