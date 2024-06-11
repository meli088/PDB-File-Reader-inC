#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h> // Pour sqrt

// Structure pour stocker les informations d'un atome
typedef struct Atom {
    char recordName[7]; // "ATOM  " ou "HETATM"
    int serial; // Numéro de série
    char name[5]; // Nom de l'atome
    char altLoc[2]; // Indicateur d'emplacement alternatif
    char resName[4]; // Nom de résidu
    char chainID[2]; // Identifiant de chaîne
    int resSeq; // Nombre seq résidu
    char iCode[2]; // Code pour insertion de résidus
    float x, y, z; // Coordonnées
    float occupancy;  // Occupance
    float tempFactor;  // Température
    char element[3]; // Symbole de l'élément
    char charge[3]; // Charge de l'atome
} Atom;

// Structure pour stocker une liste d'atomes
typedef struct Structure {
    Atom *atom; // Pointeur vers un atome
    struct Structure *next; // Pointeur vers le prochain élément de la liste
} Structure;

// Structure pour stocker les sommes et les compteurs
typedef struct ChainSum { 
    float sumX, sumY, sumZ, totalMass;
    int count;
} ChainSum;

// Structure pour stocker des informations sur les ligands
typedef struct Ligands {
    char **names; // Tableau de pointeurs vers les noms des ligands
    int count;    // Nombre de ligands uniques trouvés
} Ligands;

// Fonction pour lire un fichier PDB et construire une liste chaînée d'atomes
Structure *read_pdb_file(char* filename) {
    FILE *file = fopen(filename, "r");  // Ouvre le fichier PDB en lecture
    if (!file) {
        printf("Impossible d'ouvrir le fichier %s\n", filename);
        return NULL;
    }
    
    Structure *head = NULL, *current = NULL;  // Initialise les pointeurs de tête (head) et de l'élément courant (current) de la liste chaînée  
    char line[201];
    while (fgets(line, sizeof(line), file)) {  // Lit le fichier ligne par ligne
        if (strncmp(line, "ATOM", 4) == 0 || strncmp(line, "HETATM", 6) == 0) {  // Vérifie si la ligne correspond à un atome (ATOM ou HETATM)
            Atom *new_atom = (Atom*)malloc(sizeof(Atom));  // Allocation mémoire pour un nouvel atome
            sscanf(line, "%6s%5d %4s%c%3s %c%4d%c %8f%8f%8f%6f%6f %2s%2s",  // Lit les informations de l'atome depuis la ligne lue
                   new_atom->recordName,
                   &new_atom->serial,
                   new_atom->name,
                   new_atom->altLoc,
                   new_atom->resName,
                   new_atom->chainID,
                   &new_atom->resSeq,
                   new_atom->iCode,
                   &new_atom->x,
                   &new_atom->y,
                   &new_atom->z,
                   &new_atom->occupancy,
                   &new_atom->tempFactor,
                   new_atom->element,
                   new_atom->charge);

            Structure *new_noeud = (Structure*)malloc(sizeof(Structure));  // Crée un nouveau nœud de la liste chaînée pour stocker l'atome
            new_noeud ->atom = new_atom;
            new_noeud ->next = NULL;

            if (!head) {  // Si la liste est vide (head est NULL), initialise la liste avec le nouvel atome
                head = new_noeud ;
                current = new_noeud ;
            } else {  // Sinon, ajoute le nouvel atome à la fin de la liste chaînée
                current->next = new_noeud ;
                current = new_noeud ;
            }
        }
    }

    fclose(file);
    return head;
}

double calculate_distance(Atom atom1, Atom atom2) {
    double dx = atom2.x - atom1.x;
    double dy = atom2.y - atom1.y;
    double dz = atom2.z - atom1.z;

    return sqrt(dx * dx + dy * dy + dz * dz);
}

void center_geometry(Structure *protein) {
    ChainSum chainSums[256] = {0}; // Initialiser tous les champs à 0

    // Parcourir tous les atomes et calculer les sommes pour chaque chaîne
    Structure *current = protein;
    while (current != NULL) {
        Atom *atom = current->atom;
            int index = atom->chainID[0]; // Utiliser le premier caractère de chainID comme index
            chainSums[index].sumX += atom->x;
            chainSums[index].sumY += atom->y;
            chainSums[index].sumZ += atom->z;
            chainSums[index].count++;
        current = current->next;
    }

    // Calculer et afficher le centre géométrique pour chaque chaîne
    for (int i = 0; i < 256; i++) {
        if (chainSums[i].count > 0) {
            float centerX = chainSums[i].sumX / chainSums[i].count;
            float centerY = chainSums[i].sumY / chainSums[i].count;
            float centerZ = chainSums[i].sumZ / chainSums[i].count;
            printf("CHAINE %c: Centre géométrique -> x = %.3f, y = %.3f, z = %.3f\n", i, centerX, centerY, centerZ);
        }
    }
}

const float massC = 12;
const float massN = 14;
const float massO = 16;
const float massS = 32;
const float massH = 1;

float get_atomic_mass(char elementSymbol[3]) {
    if (strncmp(elementSymbol, "C", 1) == 0) return massC;
    if (strncmp(elementSymbol, "N", 1) == 0) return massN;
    if (strncmp(elementSymbol, "O", 1) == 0) return massO;
    if (strncmp(elementSymbol, "S", 1) == 0) return massS;
    if (strncmp(elementSymbol, "H", 1) == 0) return massH;
    return 0; // En cas d'élément non géré
}


void center_of_mass_chains(Structure *protein) {
    ChainSum chainSum[256] = {0};

    // Parcourir tous les atomes et calculer les sommes pondérées par la masse pour chaque chaîne
    Structure *current = protein;
    while (current != NULL) {
        Atom *atom = current->atom;
            int index = atom->chainID[0];
            float mass = get_atomic_mass(atom->element);
            chainSum[index].sumX += atom->x * mass;
            chainSum[index].sumY += atom->y * mass;
            chainSum[index].sumZ += atom->z * mass;
            chainSum[index].totalMass += mass;
        
        current = current->next;
    }

    // Calculer et afficher le centre de masse pour chaque chaîne
    for (int i = 0; i < 256; i++) {
        if (chainSum[i].totalMass > 0) {
            float centerX = chainSum[i].sumX / chainSum[i].totalMass;
            float centerY = chainSum[i].sumY / chainSum[i].totalMass;
            float centerZ = chainSum[i].sumZ / chainSum[i].totalMass;
            printf("CHAINE %c: Centre de masse -> x = %.3f, y = %.3f, z = %.3f\n", i, centerX, centerY, centerZ);
        }
    }
}


void detect_contacts(Structure *head, float seuil) {
    Structure *current1, *current2;     // current1 et current2 sont des pointeurs utilisés pour parcourir la liste chaînée de structures (protéines)
    for (current1 = head; current1 != NULL; current1 = current1->next) {  // Cette boucle imbriquée parcourt la liste à partir de l'élément suivant current1 jusqu'au dernier
        for (current2 = current1->next; current2 != NULL; current2 = current2->next) {
            if (current1->atom->chainID[0] != current2->atom->chainID[0]) {
                double distance = calculate_distance(*(current1->atom), *(current2->atom));
                if (distance < seuil) {
                    printf("Contact entre l'atome %d (chaîne %s) et l'atome %d (chaîne %s) à moins de %.2f Å\n",
                           current1->atom->serial, current1->atom->chainID,
                           current2->atom->serial, current2->atom->chainID, seuil);
                }
            }
        }
    }
}


Ligands* extract_ligand_names(char* filename) { // La variable 'file' sera un pointeur vers le fichier ouvert
    FILE *file = fopen(filename, "r");
    if (!file) {
        printf("Impossible d'ouvrir le fichier %s\n", filename);
        return NULL;
    }

    Ligands* ligands = malloc(sizeof(Ligands));  // Allocation initiale de mémoire pour stocker la structure Ligands.
    ligands->names = NULL;  // Pointeur vers le tableau de noms de ligands initialisé à NULL
    ligands->count = 0;

    char line[201];
    char ligandResName[4];  // Déclaration d'un tableau pour stocker temporairement le nom du ligand extrait de la ligne
    while (fgets(line, sizeof(line), file)) {
        if (strncmp(line, "HETATM", 6) == 0) {
            strncpy(ligandResName, line + 17, 3);
            ligandResName[3] = '\0'; // Chaîne correctement terminée
            
            // Boucle pour vérifier si le nom du ligand extrait est déjà stocké dans la structure Ligands
            int i;
            for (i = 0; i < ligands->count; i++) {
                if (strcmp(ligands->names[i], ligandResName) == 0) break;  // Si le nom est déjà présent, sort de la boucle.
            }
            if (i == ligands->count) { // Nouveau ligand trouvé
                ligands->names = realloc(ligands->names, (ligands->count + 1) * sizeof(char*)); // Réalloue de la mémoire pour agrandir le tableau de noms de ligands et tenir compte du contenu du bloc original
                ligands->names[ligands->count] = strdup(ligandResName);  // Duplique le nom du ligand extrait et stocke le pointeur dans le tableau agrandi sinon toutes les entrées pointent vers la même chaîne de caractères ligandResName
                ligands->count++;
            }
        }
    }

    fclose(file);
    return ligands;
}


void detect_protein_ligand_contacts(Structure *protein, Ligands* ligands, float seuil) {  // Parcourt tous les noms de ligands stockés dans la structure Ligands
    for (int l = 0; l < ligands->count; l++) {
        char* ligandResName = ligands->names[l];  // Récupère le nom du ligand actuel à partir du tableau des noms
        Structure *current = protein;  // Parcourt tous les atomes de la protéine pour trouver ceux qui correspondent au ligand actuel
        while (current != NULL) {
            if (strcmp(current->atom->resName, ligandResName) == 0) { // Compare le nom du résidu de l'atome courant avec le nom du ligand
                Structure *other = protein;
                while (other != NULL) {
                    if (strcmp(other->atom->resName, ligandResName) != 0) {  // Vérifie que l'atome "other" n'est pas le même ligand que "current"
                        double distance = calculate_distance(*(current->atom), *(other->atom));
                        if (distance <= seuil) {
                            printf("Contact entre le ligand %s (atome %d) et la protéine (atome %d) à une distance de %.2f Å\n",
                                   ligandResName, current->atom->serial, other->atom->serial, distance);
                        }
                    }
                    other = other->next;  // Passe à l'atome suivant dans la protéine
                }
            }
            current = current->next;
        }
    }
}

void free_structure(Structure *head) {
    Structure *tmp;
    while (head != NULL) {
        tmp = head;
        head = head->next;
        free(tmp->atom); // Libère la mémoire de l'atome
        free(tmp); // Libère la mémoire du noeuf
}
}

int main(int argc, char *argv[]) {
    // Vérifie que l'utilisateur a bien passé un nom de fichier en argument
    if (argc != 2) {
        fprintf(stderr, "Pas de fichier fourni avec %s\n", argv[0]);
        return 1;
    }

    // Affiche un message indiquant le début de la lecture du fichier
    printf("Lecture du fichier PDB '%s'en cours :)\n", argv[1]);

    // Lecture du fichier PDB
    Structure *protein = read_pdb_file(argv[1]);

    if (!protein) {
        printf("Aucune donnée lue du fichier PDB. Assurez-vous que le fichier existe et est au format correct.\n");
        return 1;
    }

    // Affichage du nombre d'atomes lus pour vérification
    Structure *current = protein;
    int count = 0;
    while (current != NULL) {
        count++;
        current = current->next;
    }
    printf("Nombre total d'atomes lus : %d\n", count);

    // Calcul et affichage du centre géométrique
    printf("\nCalcul et affichage du centre géométrique pour chaque chaîne :\n");
    center_geometry(protein);

    // Calcul et affichage du centre de masse
    printf("\nCalcul et affichage du centre de masse pour chaque chaîne :\n");
    center_of_mass_chains(protein);

    // Detection contact
    printf("\nDétection des contacts entre chaînes à un seuil de 3 Å\n");
    detect_contacts(protein, 3.0);


    // Lecture du fichier PDB et extraction du nom du ligand
    Ligands* ligands = extract_ligand_names(argv[1]);
    if (!ligands || ligands->count == 0) {
        printf("Aucun ligand trouvé dans le fichier %s.\n", argv[1]);
        return 1;
    }

    printf("\nLigands trouvés (%d): \n", ligands->count);
    for (int i = 0; i < ligands->count; i++) {
        printf("- %s\n", ligands->names[i]);
    }


    // Détection des contacts entre la protéine et les ligands à un seuil de 3 Å
    printf("\nDétection des contacts entre la protéine et les ligands à un seuil de 3 Å\n");
    detect_protein_ligand_contacts(protein, ligands, 3.0);


    // Libérer la mémoire allouée
    free_structure(protein);
    free(ligands->names);
    free(ligands);

    return 0;
}

