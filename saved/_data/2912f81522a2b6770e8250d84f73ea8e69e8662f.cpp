#include <stdio.h>

// a decommenter si vous voulez utiliser les couts
#include <iostream>
#include <ratio>
#include <thread>
#include <chrono>
#include <omp.h>

typedef std::chrono::high_resolution_clock Clock;

#define VALMAX 48

// Ordre interne
//  aiguille d'une montre :
// j=0 5 4 3 2 1 O
//j=1 5 4 3 2 1 0
// voir tableau des next

#include "minimax_graph.h"

void init_position(Position* pos) {
	#pragma omp parallel for
	for(int i=0;i<6;i++){
		pos->_Cases[0][i]=0;
		pos->_Cases[1][i]=0;
	}
	pos->_PionsPris[0]=pos->_PionsPris[1]=0;
}

void print_position(Position* pos){
	printf("--------------------------\n");
	for(int i = 5 ; i >= 0 ; i--){
		printf("[%d] " , pos->_Cases[0][i]);
	}
	printf(" PP=%d\n",pos->_PionsPris[0]);
	for(int i = 5 ; i >= 0 ; i--){
		printf("[%d] " , pos->_Cases[1][i]);
	}
	printf(" PP=%d\n--------------------------\n",pos->_PionsPris[1]);
}

void print_position_ordi_bas_inv(Position* pos){
	printf("--------------------------\n");
	for(int i = 0 ; i <6 ; i++){
		printf("[%d] " , pos->_Cases[1][i]);
	}
	printf(" PP=%d\n",pos->_PionsPris[1]);
	for(int i = 0 ; i <6 ; i++){
		printf("[%d] " , pos->_Cases[0][i]);
	}
	printf(" PP=%d\n--------------------------\n",pos->_PionsPris[0]);
}

void print_position_ordi_haut_inv(Position* pos){
	printf("--------------------------\n");
	for(int i = 5 ; i >= 0 ; i--){
		printf("[%d] " , pos->_Cases[0][i]);
	}
	printf(" PP=%d\n",pos->_PionsPris[0]);
	for(int i = 5 ; i >=0 ; i--){
		printf("[%d] " , pos->_Cases[1][i]);
	}
	printf(" PP=%d\n--------------------------\n",pos->_PionsPris[1]);
}

// Affichage
// 12 11 10 9 8 7
// 1 2 3 4 5 6
// on tourne dans le sens inverse des aiguilles d'une montre

void print_position_ordi_bas(Position* pos){
	printf("--------------------------\n");
	for(int i = 5 ; i >= 0 ; i--){
		printf("[%d] " , pos->_Cases[1][i]);
	}
	printf(" PP=%d\n",pos->_PionsPris[1]);
	for(int i = 5 ; i >= 0 ; i--){
		printf("[%d] " , pos->_Cases[0][i]);
	}
	printf(" PP=%d\n--------------------------\n",pos->_PionsPris[0]);
}

void print_position_ordi_haut(Position* pos){
	printf("--------------------------\n");
	for(int i = 0 ; i < 6 ; i++){
		printf("[%d] " , pos->_Cases[0][i]);
	}
	printf(" PP=%d\n",pos->_PionsPris[0]);
	for(int i = 0 ; i <6 ; i++){
		printf("[%d] " , pos->_Cases[1][i]);
	}
	printf(" PP=%d\n--------------------------\n",pos->_PionsPris[1]);
}

int est_affame(Position* pos, const int joueur){
	// TODO faire une somme cela ira plus vite
	int somme=0;
	for(int i=0;i<6;i++){
//		if (pos->_Cases[joueur][i] > 0) return 0;
		somme += pos->_Cases[joueur][i];
	}
	return !somme;
//	return 1;
}

void copier(Position* newPos, Position* pos){
	#pragma omp parallel for
	for(int i=0;i<6;i++){
		newPos->_Cases[0][i]=pos->_Cases[0][i];
		newPos->_Cases[1][i]=pos->_Cases[1][i];
	}
	newPos->_PionsPris[0]=pos->_PionsPris[0];
	newPos->_PionsPris[1]=pos->_PionsPris[1];
}

int jouer_coup(CaseSuivante* cs, Position* newPos, Position* pos, const int joueur, const int coup){
	// retourne 1 si le coup est jouable et 0 sinon
	const int nbpions=pos->_Cases[joueur][coup];
//	cout << "joueur: " << joueur << " coup: " << coup << " nbpions : " << nbpions << endl;
	if (nbpions == 0){
//		cout << "pas de pions" << endl;
		return 0; // le cout n'est pas jouable
	}
	// on met a 0 la structure
	copier(newPos,pos);
	newPos->_Cases[joueur][coup]=0;
	int j=joueur;
	int c=coup;
	for(int i=0;i<nbpions;i++){
		const int tj=j;
		j=cs->_Jnext[j][c];
		c=cs->_Cnext[tj][c];
//		cout << "[" << j << "][" << c << "] ";
		newPos->_Cases[j][c]++;
	}
//	cout << endl;
//	print_position(newPos);
	// on regarde si la case de d�part est vide ou pas.
	int nbp=newPos->_Cases[joueur][coup];
	while (nbp != 0){
		// on distribue ces pions
		newPos->_Cases[joueur][coup]=0;
		for(int i=0;i<nbp;i++){
			const int tj=j;
			j=cs->_Jnext[j][c];
			c=cs->_Cnext[tj][c];
			newPos->_Cases[j][c]++;
		}
		nbp=newPos->_Cases[joueur][coup];
	}
	if (j != joueur){// on regarde si on doit prendre des pions
		if (j ==0){
			for(int i=c;i<=5;i++){
				if (newPos->_Cases[j][i] == 2 || newPos->_Cases[j][i] == 3){
					newPos->_PionsPris[joueur] += newPos->_Cases[j][i];
					newPos->_Cases[j][i]=0;
				} else {
					break;
				}
			}
		} else {
			for(int i=c;i>=0;i--){
				if (newPos->_Cases[j][i] == 2 || newPos->_Cases[j][i] == 3){
					newPos->_PionsPris[joueur] += newPos->_Cases[j][i];
					newPos->_Cases[j][i]=0;
				} else {
					break;
				}
			}
		}
	}
	// on regarde si l'adversaire est affame
	if (est_affame(newPos,!joueur)) return 0;
	return 1;
}

inline
int evaluer(Position* pos){
	return pos->_PionsPris[0] - pos->_PionsPris[1];
}

int test_fin(Position* pos){
	if (pos->_PionsPris[0] + pos->_PionsPris[1] > 24){
		if (48-pos->_PionsPris[0]-pos->_PionsPris[1] <= 6){
			return 1;
		}
		if (pos->_PionsPris[0] >= 25) return 1;
		if (pos->_PionsPris[1] >= 25) return 1;
	}
	return 0;
}


struct ECoup {
	int _Val[6];
	int _Coup[6];
};

typedef struct ECoup EvalCoup;

long long NUM_MINIMAX=0;
int VALMM=0;

void print_eval_coup(EvalCoup* ec, int nb){
	for(int i=0;i<nb;i++){
		std::cout << "coup: " << ec->_Coup[i] << " eval: " << ec->_Val[i] << std::endl;
	}
}


int calculer_eval_coup(EvalCoup* ec, CaseSuivante* cs, Position* pos,const int joueur, const int alpha, const int beta, const int pmax){
	int nbv=0;
	Position newPos;
	for(int i=0;i<6;i++){
		if (jouer_coup(cs,&newPos,pos,joueur,i)){
			ec->_Val[nbv]=valeur_minimax(cs,&newPos,!joueur,alpha,beta,pmax-1);
			ec->_Coup[nbv]=i;
			nbv++;
		}
	}
	return nbv;
}

int valeur_minimax(CaseSuivante* cs, Position* pos,const int joueur, int alpha,int beta, const int pmax){
	// Calcule la valeur de e pour le joueur J selon que e EstUnEtatMax ou pas et pmax la
	NUM_MINIMAX++;
	if (pos->_PionsPris[0] + pos->_PionsPris[1] > 24){
		if (48-pos->_PionsPris[0]-pos->_PionsPris[1] <= 6){
			if (pos->_PionsPris[0] > pos->_PionsPris[1]) return VALMAX;
			if (pos->_PionsPris[0] < pos->_PionsPris[1]) return -VALMAX;
			return 0;
		}
		if (pos->_PionsPris[0] >= 25) return VALMAX;
		if (pos->_PionsPris[1] >= 25) return -VALMAX;
	}
	if (pmax == 0) return evaluer(pos);
	EvalCoup ec;
	const int nbv=calculer_eval_coup(&ec,cs,pos,joueur,alpha,beta,pmax);
	//print_eval_coup(&ec,nbv );
	int imin=0;
	if (joueur==0){ // on prend le max de vals
		for(int i=1;i<nbv;i++){
			if (ec._Val[i] > ec._Val[imin]){
				imin=i;
			}
		}
	} else { // on prend le min de vals
		for(int i=1;i<nbv;i++){
			if (ec._Val[i] < ec._Val[imin]){
				imin=i;
			}
		}
	}
//	cout << "CHOIX joueur: " << joueur << " coup: " << ec._Coup[imin] << " minimax: " << ec._Val[imin] << endl;
	return ec._Val[imin];
}

int decision(CaseSuivante* cs,Position* pos,int pmax){
	// on regarde le nombre de case vide et on ajoute de la profondeur eventuellement
	// k case vide = profmax * k/12
	int k=0;
	for(int i=0;i<6;i++){
		if (pos->_Cases[0][i] == 0) k++;
//		if (pos->_Cases[1][i] == 0) k++;
	}
//	k--;
	if (k >0){
		if (pos->_PionsPris[1] + pos->_PionsPris[0] >= 20) pmax++;

		if (k <= 3) pmax++;
		if (k > 3) pmax +=2;
	}
	int alpha=-VALMAX-50;
	int beta=VALMAX+50;
	EvalCoup ec;
	const int nbv=calculer_eval_coup(&ec,cs,pos,0,alpha,beta,pmax);
	//print_eval_coup(&ec,nbv );
	int imin=0;
	for(int i=1;i<nbv;i++){
		if (ec._Val[i] > ec._Val[imin]){
			imin=i;
		}
	}
	VALMM=ec._Val[imin];
	return ec._Coup[imin];
}

int valeur_minimaxAB(CaseSuivante* cs, Position* pos,const int joueur, int alpha,int beta, const int pmax, const bool gagne);

int calculer_coup(CaseSuivante* cs, Position* pos,const int joueur, int alpha, int beta, const int pmax,const bool gagne){
	Position newPos;
	if (joueur==0){// MAX
		//cout << "MAX";
		for(int i=0;i<6;i++){
			if (jouer_coup(cs,&newPos,pos,joueur,i)){
				const int val=valeur_minimaxAB(cs,&newPos,!joueur,alpha,beta,pmax-1,gagne);
			//	cout << " coup: " << i << " val: " << val << " alpha: " << alpha << " beta:" << beta << endl;
				if (val > alpha) {
					alpha=val;
				}
				if (alpha >= beta){
					//cout << "a";
					return alpha;
				}
			}
		}
		return alpha;
	}
	//	cout << "MIN";
	for(int i=0;i<6;i++){
		if (jouer_coup(cs,&newPos,pos,joueur,i)){
			const int val=valeur_minimaxAB(cs,&newPos,!joueur,alpha,beta,pmax-1,gagne);
		//		cout << " coup: " << i << " val: " << val << " alpha: " << alpha << " beta:" << beta << endl;
			if (val < beta){
				beta=val;
			}
			if (beta <= alpha){
				//cout << "b";
				return beta;
			}
		}
	}
	return beta;
}

/// <summary> Fonction qui calcule la valeur minimax avec coupes Alpha-Beta </summary>
int valeur_minimaxAB(CaseSuivante* cs, Position* pos,const int joueur, int alpha,int beta, const int pmax, const bool gagne){
	// Calcule la valeur de e pour le joueur J selon que e EstUnEtatMax ou pas et pmax la
	NUM_MINIMAX++;
	int ajoutProf=(gagne) ? pmax : 0;
	// TODO rajouter un parametre gagne qui est vrai ou faux
	// puis on definit un parametre ajoutProf en fn de cette variable
	if (pos->_PionsPris[0] + pos->_PionsPris[1] > 24){
		if (48-pos->_PionsPris[0]-pos->_PionsPris[1] <= 6){
			if (pos->_PionsPris[0] > pos->_PionsPris[1]){
				//cout << "position finale gagnante joueur: " << joueur << " pmax: " << pmax << endl;
				//print_position(pos);
				return VALMAX + ajoutProf;
			}
			if (pos->_PionsPris[0] < pos->_PionsPris[1]) return -VALMAX-ajoutProf;
			return 0;
		}
		if (pos->_PionsPris[0] >= 25){
			//cout << "position finale gagnante joueur: " << joueur << " pmax: " << pmax << endl;
			return VALMAX+ajoutProf;
		}
		if (pos->_PionsPris[1] >= 25) return -VALMAX-ajoutProf;
	}
	if (pmax == 0) return evaluer(pos);
	return calculer_coup(cs,pos,joueur,alpha,beta,pmax,gagne);
}

/// Fonction qui d�termine la premi�re d�cision avec utilisation d'un algo minimax et de coupes
/// Alpha-Beta.
// Jean Charles, 19/05/2011.

int decisionAB(CaseSuivante* cs,Position* pos,int pmax, bool gagne){
	// on regarde le nombre de case vide et on ajoute de la profondeur eventuellement
	// k case vide = profmax * k/12
	int k=0;
	for(int i=0;i<6;i++){
		if (pos->_Cases[0][i] == 0) k++;
//		if (pos->_Cases[1][i] == 0) k++;
	}
//	k--;
	if (k >0){
		if (pos->_PionsPris[1] + pos->_PionsPris[0] >= 20) pmax++;
		if (pos->_PionsPris[1] + pos->_PionsPris[0] >= 25) pmax++;

		if (k <= 3) pmax++;
		if (k > 3) pmax +=2;
	}

	int alpha=-VALMAX-50; // avant -1
	int beta=VALMAX+50; // avant +1
	Position newPos;
	int coup;
	for(int i=0;i<6;i++){
		if (jouer_coup(cs,&newPos,pos,0,i)){
			const int val=valeur_minimaxAB(cs,&newPos,1,alpha,beta,pmax-1,gagne);
			if (val > alpha) {
				alpha=val;
				coup=i;
			}
		}
	}

	VALMM=alpha;
	return coup;
}

/// <summary> Position a partir de laquelle on commence � jouer </summary>
void position_debut(Position* pos) {
	for(int i=0;i<6;i++){
		pos->_Cases[0][i]=4;
		pos->_Cases[1][i]=4;
	}
	pos->_PionsPris[0]=pos->_PionsPris[1]=0;

}
int main(int argc, char* argv[]){
	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();
	CaseSuivante cs;
	cs._Cnext[0][0]=0;
	cs._Jnext[0][0]=1;
	cs._Cnext[0][1]=0;
	cs._Jnext[0][1]=0;
	cs._Cnext[0][2]=1;
	cs._Jnext[0][2]=0;
	cs._Cnext[0][3]=2;
	cs._Jnext[0][3]=0;
	cs._Cnext[0][4]=3;
	cs._Jnext[0][4]=0;
	cs._Cnext[0][5]=4;
	cs._Jnext[0][5]=0;

	cs._Cnext[1][0]=1;
	cs._Jnext[1][0]=1;
	cs._Cnext[1][1]=2;
	cs._Jnext[1][1]=1;
	cs._Cnext[1][2]=3;
	cs._Jnext[1][2]=1;
	cs._Cnext[1][3]=4;
	cs._Jnext[1][3]=1;
	cs._Cnext[1][4]=5;
	cs._Jnext[1][4]=1;
	cs._Cnext[1][5]=5;
	cs._Jnext[1][5]=0;

	Position pos;
	Position newPos;
	position_debut(&pos);
/**
	printf("DEBUT DU JEU AWALE\n");
	printf("(0) l'ordinateur commence\n");
	printf("(1) le joueur commence\n");
**/
	int joueur;
	if(scanf("%d",&joueur)){}


	int ordiCommence= (joueur==0)? 1 : 0;
	//printf("C'est parti!\n");
	int fin=0;
	bool gagne=false;
	int numeroCoup = 1;
	while(!fin){
		int coup;
		if (joueur == 0){ // l'ordi JOUE
//			coup=decision(&cs,&pos,10); // ATTENTION C'est LA PROFONDEUR MAX
			coup=decisionAB(&cs,&pos,17,gagne); // ATTENTION C'est LA PROFONDEUR MAX 11
			if (!gagne && VALMM==48){gagne=true;}
//			if (!gagne && VALMM >= 48) { gagne = true; }
			//printf("Noeuds traites: %d valeur minimax:%d\n", NUM_MINIMAX, VALMM);
			int cj;
			if (ordiCommence){
				cj=6-coup;
			} else {
				cj=12-coup;
			}
			//printf("COUP JOUE: %d\n",cj);
			NUM_MINIMAX=0;
			jouer_coup(&cs,&newPos,&pos,joueur,coup);
			copier(&pos,&newPos);
/**
			if (ordiCommence){
				print_position_ordi_haut_inv(&pos);
			} else {
				print_position_ordi_bas_inv(&pos);
			}
**/
		} else { // le JOUEUR JOUE
			if (ordiCommence){
				//printf("selectionner une case 12 11 10 9 8 7\n");
				if(scanf("%d",&coup)){}
				coup -=7;
			} else {
				//printf("selectionner une case 1 2 3 4 5 6\n");
				if(scanf("%d",&coup)){}
				//coup=6-coup;
				coup--;
			}
			//printf("COUP interne JOUE: %d\n",coup);

			jouer_coup(&cs,&newPos,&pos,joueur,coup);

			copier(&pos,&newPos);
/**
			if (ordiCommence){
				print_position_ordi_haut_inv(&pos);
			} else {
				print_position_ordi_bas_inv(&pos);
			}
**/
		}
		fin=test_fin(&pos);
		joueur = !joueur;
	}
		//printf("Fin de la partie!\n");
    end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds = end - start;
		std::cout << elapsed_seconds.count() << std::endl;
}

// REMARQUE : Une fois que l'on a trouv� un chemin gagnant a coup sur il faut le prendre et ne pas laisser le systeme en calculer un autre
// car dans ce cas il peut changer d'avis et en fait faire ces cycles !!!
// TODO : il faut implementer cela ici pour eviter ce mauvais comportement

// TODO quand on est gagant il faut simplement obliger le plus court chemin !!!
// On sait qu'on a gagne ALORS ON recherche le chemin le plus court
