/**
 * Programme  de l'attaque TEMPLATE sur le calcul f(x+k)
 *
 *
 * brancher le canal A sur la r�sistance et le canal B sur l'I/O
 *
 * Le programme les algorithmes et 1 et 2 du papier.
 *
 **/

// pour l'attaque sur le xor, activer cette directive
// sinon, l'attaque est sur la sbox
//#define XOR

// pour l'attaque am�lior�e activer cette directive
//#define IMPROVED
// d�calage des mesures

#define DELTAT 525

// pour �crire les fichiers activer cette d�finition
//#define ECRIRE_FICHIERS
// pour analyser directement les r�sultats, activer cette d�finition
#define ANALYSER

// activer cette directive pour afficher les informations
// (APDU et SW)
//#define VERBOSE
// nombre d'it�rations (pour les essais)



#define SEUIL_CORR 0.7

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include <setjmp.h>

#include <time.h>

unsigned int CENTRAL_POINT = 350; // Default value



#define TRUE_KEY 0x00 //0xA5
#define NB_SAMPLES 700// 64000 //32*1024  //nombre d'�chantillon

// chemin o� �crire les donn�es
#define CHEMIN  "./temp/"







///////////////////////////////////////////////////////////////////
// calcul de la transformee de Fourier de la fonction f
void Fourier_float(float*f,int n)
{
    unsigned int nn; // taille du tableau = 2^n
    unsigned int jj; // taille de la sous-fonction
    int i;
    float t;   // variable temporaire
    nn=1<<n;   // 2^n
    jj=nn>>1;  // premi�re sous-fonction de taille 2^{n-1}

    do
    {
        i=0;
        do
        {  // papillon sur sous-fonction
            t=f[i]-f[i+jj];
            f[i]+=f[i+jj];
            f[i+jj]=t;
            i++;
            if ((i&(jj-1))==0) i+=jj; // si sous-fonction finie, passer � la suivante
        }
        while (i<nn); // jusqu'� terminer la table
        jj>>=1;       // recommencer avec des sous-fonctions de taille moiti�
    }
    while (jj!=0);  // jusqu'aux sous-fonctions de taille 1
}

/////////////////////////////////////////////////////////////

float bufferA[NB_SAMPLES];





////////////////////////////////////////////////////////////////
// r�glage du Picoscope



///////////////////////////////////////////////////////

// affecte l'atr au tableau atr et sa longueur � lenAtr
// lenAtr doit �tre initialis� avec la taille acceptable pour le tableau
//void ATR(uint8_t*atr, DWORD*lenAtr)
////////////////////////////////////////////////////////////////////////////


// cla, ins, p1, p2, p3 est le header APDU
// data_buffer contient les donn�es introduites dans la carte pour une commande entrante
// les donn�es en provenance de la carte pour les commandes sortantes
// len_adata vaut 0 pour une commande sortante
// vaut p3 pour une commande entrante

//void SendAPDU(uint8_t*sw1, uint8_t* sw2

//////////////////////////////////////////////////////////////////////

//void atr()





//////////////////////////////////////////////////////////





// nombre de mesures traitees
#define NBMESURES 1 //700//40000// 1023 Max //c'est le D (la fenetre du temps)
#include "Matrice.h"
#define NB_ITE_ATTACKK 200 //7// le nombre maximal de mesure lors de l'attaque
//#define NB_ITE_PROFILING 256//le nombre maximal de mesure lors du profiling

#define Q 51200//2560// profilling Q il ya 51200 pour MSP430 et 50000 pour ASCAD
#define ATTACKQ 11776 //2560//Attack Q il ya 11776 pour MSP430 et  10000  pour ASCAD

#define NBM 256


float ymoyen[NBMESURES][256]; // un tableau de fonctions models de 256 valeurs
float ymoyentrans[256][NBMESURES]; // un tableau de fonctions models de 256 valeurs
float x[NBMESURES][Q]; // Q trace
float xmoyentran[NBMESURES][256];
float segmainverse_fois_ymoyen[NBMESURES][256];
float ymoyen_fois_segmainverse_fois_ymoyen[256][256]; //ycal
int n[256]; // le nombre de  trace par message

//float xmoyen[NBMESURES][256]; // un tableau de fonctions models de 256 valeurs
//int t[Q];//Q messages

// seul � partir duquel on teste la fin du calcul
// lorsque la carte r�pond, elle envoie le premier octet du status word
// Le niveau d'E/S baisse au moment du bit start et permet de synchroniser la mesure
// ce seuil change selon un lecteur 3V ou un lecteur 5V.
// pour un lecteur 3V 19000 marche et pour un lecteur 5V, 32000 marche

////////////////////////////////////////////////////////////////////////////////


char nom[128];

FILE * Profiling_traces_traces_file;

FILE * Profiling_traces_labels_file;

FILE * Attack_traces_traces_file;

FILE * Attack_traces_labels_file;


void open_files(void)
{
///sprintf(nom,  "C:/Users/A/Desktop/Pratique/template_offline/ASCADtxt/Profiling_traces_traces.txt");
sprintf(nom,  "./Profiling_traces_traces.txt");

   			   Profiling_traces_traces_file = fopen(nom,"r");
   			   if (Profiling_traces_traces_file==NULL)
    		         {
    		     	    printf("** Erreur d'ouverture du fichier1 '%s'** \n", nom);
    			   }

sprintf(nom,  "./Profiling_traces_labels.txt");

   			   Profiling_traces_labels_file = fopen(nom,"r");
   			   if (Profiling_traces_labels_file==NULL)
    		         {
    		     	    printf("** Erreur d'ouverture du fichier2 '%s'** \n", nom);
    			   }

sprintf(nom,  "./Attack_traces_traces.txt");

   			   Attack_traces_traces_file = fopen(nom,"r");
   			   if (Attack_traces_traces_file==NULL)
    		         {
    		     	    printf("** Erreur d'ouverture du fichier3 '%s'** \n", nom);
    			   }

sprintf(nom,  "./Attack_traces_labels.txt");

   			   Attack_traces_labels_file = fopen(nom,"r");
   			   if (Attack_traces_labels_file==NULL)
    		         {
    		     	    printf("** Erreur d'ouverture du fichier4 '%s'** \n", nom);
    			   }



}


/////////////////////////////////////////////////////////////////////////////

void close_files()
{

fclose(Profiling_traces_traces_file);
fclose(Profiling_traces_labels_file);
fclose(Attack_traces_traces_file);
fclose(Attack_traces_labels_file);

}
/////////////////////////////////////////////////////////////////////////////
void capture_template(int message_index)
{
    int i,j;

int label;




for(i=0;i<NB_SAMPLES;i++)bufferA[i]=0.0;

    // initialisation des donn�es de l'apdu avec les valeurs de x introduites





fscanf(Profiling_traces_labels_file, "%d\n", &label);

//printf("trace=%d label=%d \n", message_index, label);
//scanf("%*c");

for(j=0;j<NB_SAMPLES;j++)
       		{
			fscanf(Profiling_traces_traces_file, "%f,", &bufferA[j]);


			//printf("bufferA[%d]=%f \n", j,bufferA[j]);
			//scanf("%*c");
			}
//printf("je suis la \n");




	//t[q]=label;// ii;
	n[label^(TRUE_KEY)]++;// ii^la clé
//printf("je suis la la \n");
    //printf("t[%d]=%d",q,t[q]);
    for (j=0;j<NBMESURES;j++)
    {
	//printf("\n\nje suis à capture_template message_index=%d, j=%d\n\n",message_index,j);
			   //[i0-22256-(int)(NBMESURES/2)+j]
    //x[j][q]
	x[j][message_index]=bufferA[350-(int)(NBMESURES/2)+j];//39999-17744// 350=700/2 700 représente D qui est la longeur de la trace

	//ymoyen[j][message_index^(TRUE_KEY)]
        ymoyen[j][label]+=bufferA[350-(int)(NBMESURES/2)+j];
 //^la clé
    }




}


////////////////////////////////////////////////////////////////////////////

void template_estimation()
{
printf("\n\nje suis à template_estimation\n\n");

//float xtrans[Q][NBMESURES], ymoyentrans[256][NBMESURES];

float x_xtrans[NBMESURES][NBMESURES], ymoyen_ymoyentrans[NBMESURES][NBMESURES];
//float segma[NBMESURES][NBMESURES];
//float segmainverse[NBMESURES][NBMESURES];

int i, p, j;

for (i=0;i<256;i++){
 				 n[i]=0;
				for (j=0;j<NBMESURES;j++) // NBMESURES=D
				{
				 xmoyentran[j][i]=0.0;
				ymoyentrans[i][j]=0.0;


                ymoyen[j][i]=0.0;
				}
			  }

//for (i=0;i<Q;i++){ //Q nombre de traces du profiling (N), avant la coalescence
			 	//t[i]=0;
//				 for (j=0;j<NBMESURES;j++)
//				 {
//			 	x[j][i]=0.0;
			 	//xtrans[i][j]=0.0;
//			 	}
//			}

for (i=0;i<NBMESURES;i++){for (j=0;j<NBMESURES;j++)
						{x_xtrans[i][j]=0.0;
						 ymoyen_ymoyentrans[i][j]=0.0;
						}
				 }

/*
for (rr=0;rr<1;rr++)
 { printf("\n\nrr=%d\n\n",rr);
     for (r=0;r<10;r++)
    { printf("\n\nr=%d\n\n",r);


    }
 }
*/


      for (p=0;p<Q;p++)
      {
        //printf("."); fflush(stdout);
        capture_template(p);
      }
printf("fin capture_template(s) \n");

for (i=0;i<NBM;i++){

for (j=0;j<NBMESURES;j++)
						{
 if(n[i]) ymoyen[j][i]/=n[i];}
			  }

//calcule de transposé de x[][]
/*for (i=0;i<Q;i++){

for (j=0;j<NBMESURES;j++)
						{
 xtrans[i][j]=x[j][i];}
			  }
*/
//calcule de transposé de model
for (i=0;i<256;i++){

for (j=0;j<NBMESURES;j++)
						{
 ymoyentrans[i][j]=ymoyen[j][i];}
			  }


//mutiplication_matrice(x[0], NBMESURES, Q, NBMESURES, xtrans[0], x_xtrans[0])
;

mutiplication_matrice_par_sa_transposee( x[0], NBMESURES, Q, x_xtrans[0]);

reelfoismatrice(x_xtrans[0], (float)1/Q, NBMESURES, NBMESURES);

//mutiplication_matrice(ymoyen[0], NBMESURES, 256, NBMESURES, ymoyentrans[0], ymoyen_ymoyentrans[0])
;

mutiplication_matrice_par_sa_transposee( ymoyen[0], NBMESURES, 256,  ymoyen_ymoyentrans[0])
;

reelfoismatrice(ymoyen_ymoyentrans[0], 1/256, NBMESURES, NBMESURES);


//soustraction_matrice(x_xtrans[0], NBMESURES, NBMESURES, ymoyen_ymoyentrans[0], segma[0]);
soustraction_matrice(x_xtrans[0], NBMESURES, NBMESURES, ymoyen_ymoyentrans[0], x_xtrans[0]);// x_xtrans devient segma

///inverse_matrice_pivot(x_xtrans[0], NBMESURES, x_xtrans[0]);
InversionCholesky(x_xtrans[0],x_xtrans[0]);
							 //x_xtrans devient segmainverse
//inverse_matrice_pivot(segma[0], NBMESURES, segmainverse[0]);



//mutiplication_matrice(segmainverse[0], NBMESURES,  NBMESURES,256, ymoyen[0], segmainverse_fois_ymoyen[0])
;

mutiplication_matrice(x_xtrans[0], NBMESURES,  NBMESURES,256, ymoyen[0], segmainverse_fois_ymoyen[0]); //ytilde
mutiplication_matrice(ymoyentrans[0], 256,  NBMESURES,256, segmainverse_fois_ymoyen[0], ymoyen_fois_segmainverse_fois_ymoyen[0]); //ycal
for(j=0;j<NBMESURES;j++){
				//Fourier_float(xmoyentran[j],8);
				Fourier_float(segmainverse_fois_ymoyen[j],8);
				}
for(j=0;j<256;j++){
				//Fourier_float(xmoyentran[j],8);
				Fourier_float(ymoyen_fois_segmainverse_fois_ymoyen[j],8);
				}
}
////////////////////////////////////////////////////////
int All_Attack_labels[ATTACKQ];
float All_Attack_traces[ATTACKQ][NB_SAMPLES];

///////////////////////////////////////////////////////

void capture(int label, int message_index)//ii cest le label
{

    int i,j;







for(i=0;i<NB_SAMPLES;i++)bufferA[i]=0.0;

    // initialisation des donn�es de l'apdu avec les valeurs de x introduites





 //t[q]=label;
 n[label]++;
 // printf("t[%d]=%d",q,t[q]);
 for (j=0;j<NBMESURES;j++)
 {

  xmoyentran[j][label]+=All_Attack_traces[message_index][350-(int)(NBMESURES/2)+j];
  //printf("All_Attack_traces[%d][%d]=%f\n",message_index, 350-(int)(NBMESURES/2)+j ,All_Attack_traces[message_index][350-(int)(NBMESURES/2)+j]);
  //scanf("%*c");
 }




}

////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
float trace(float *mat, int dim,int k)
{
int i;
float somme=0.0;
for (i=0;i<dim;i++)somme+=*(mat + (i)*dim + (i^k));//mat[i*dim+(i^k)]
return somme;
}


/////////////////////////////////////////////////////////////////////////////
int AnalyserXor()
{

    	int i,z,j;
//

    	//printf("\nAnalyse\n");

//calcul de la transposée
/*
    for (j=0;j<NBMESURES;j++)
    {

for (i=0;i<256;i++)
      {
        xmoyentran[j][i]=xmoyen[j][i]
;
      }
    }
*/





//float xmoyentran_fois_segmainverse_fois_ymoyen[256][256]; Gain grace a fourier
float xmoyentran_fourier[NBMESURES][256];
float somme_des_produits[256];

for(i=0;i<256;i++){
	for(j=0;j<NBMESURES;j++)
	{
	  xmoyentran_fourier[j][i]=xmoyentran[j][i];
	}
			}

float temps;
clock_t t1, t2;
t1 = clock();
// Ton programme


//mutiplication_matrice(xmoyentran[0], 256,NBMESURES,  256, segmainverse_fois_ymoyen[0], xmoyentran_fois_segmainverse_fois_ymoyen[0])
;

//rappelant que  xmoyentran[NBMESURES][256];
//rappelant que  segmainverse_fois_ymoyen[NBMESURES][256];

for(j=0;j<NBMESURES;j++){
				Fourier_float(xmoyentran_fourier[j],8);
				//Fourier_float(segmainverse_fois_ymoyen[j],8);
				}

for(i=0;i<256;i++)somme_des_produits[i]=0.0;

for(i=0;i<256;i++){
	for(j=0;j<NBMESURES;j++)
	{
	  somme_des_produits[i]+=xmoyentran_fourier[j][i]*segmainverse_fois_ymoyen[j][i];
	}
			}

Fourier_float(somme_des_produits,8);
//t2 = clock();

 z=0;
//float m=trace(xmoyentran_fois_segmainverse_fois_ymoyen[0], 256,0);

float m=somme_des_produits[0];
for(i=1;i<256;i++)
		{
		//if (trace(xmoyentran_fois_segmainverse_fois_ymoyen[0], 256,i)>m)
		if (somme_des_produits[i]>m)
			{
			//m=trace(xmoyentran_fois_segmainverse_fois_ymoyen[0], 256,i);
			m=somme_des_produits[i];
			z=i;
			}
		}
t2 = clock();
temps = (float)(t2-t1)/CLOCKS_PER_SEC;
//printf("temps = %f\n", temps);
//printf(" z = %d\n", z);
return z;

/*printf("x =% 4d (0x%02x %c),  m = %f\n",
 z,z,(z>=32)&&(z<=127) ?z:'.',m);
*/


}


////////////////////////////////////////////////////////////////////////////
void travailler()
{

 int NB_ITE_ATTACK=NB_ITE_ATTACKK; //6 nombre maximal de message (traces) lors de l'attaque;
 int success_rate[NB_ITE_ATTACK];
 int p,r,i,j;
 char nom[128];

FILE  *template_succes_rate_file; //

//sprintf(nom,CHEMIN "template_walsh_hadamar_succes_rate_autour_point_corrige.txt");
sprintf(nom,CHEMIN "template_walsh_hadamar_succes_rate_autour_point_coalescance_CENTRAL_POINT%d.txt", CENTRAL_POINT);
template_succes_rate_file=fopen(nom,"a");
if (template_succes_rate_file==NULL)

{
printf("**  erreur lors de l'ouverture du fichier template_succes_rate.txt'%s' **\n", nom);
}





for(i=0;i<ATTACKQ;i++) //ATTACKQ cest le N de la phase de matching
{
	fscanf(Attack_traces_labels_file, "%d\n", &All_Attack_labels[i]);
	//printf("trace=%d label=%d \n", i, All_Attack_labels[i]);
	//scanf("%*c");

	for(j=0;j<NB_SAMPLES;j++)
       	{
		fscanf(Attack_traces_traces_file, "%f,", &All_Attack_traces[i][j]);

		//printf("All_Attack_traces[%d][%d]=%f \n", i,j,All_Attack_traces[i][j]);
		//scanf("%*c");
		}
}//fin for i




srand(time(NULL)); //reinitialisation de randomisation

for (i=0;i<NB_ITE_ATTACK;i++) success_rate[i]=0;



printf("\n\n\nNBMESURES=%d\n",NBMESURES);

 for (r=0;r<1000;r++)
//100// 500
 { printf("\nr=%d\n",r);
   //fprintf(template_succes_rate_file,"round:%d\n\n\n\n",r);

	for (i=0;i<256;i++)
	{
 		n[i]=0;
		for (j=0;j<NBMESURES;j++)
			{
			  xmoyentran[j][i]=0.0;
			}
	}



for (p=0;p<NB_ITE_ATTACK;p++) // 7
   	 {
//printf("p=%d\n",p);
        //printf("."); fflush(stdout);
	  //scanf("%*c");

	  int message_index=(int)(rand()/(double)RAND_MAX*(ATTACKQ-1));
	  //printf("message_index=%d\n",message_index);
	  //scanf("%*c");

        capture(All_Attack_labels[message_index],message_index);
//p,r*256+p



	  /*for (i=0;i<256;i++){for (j=0;j<NBMESURES;j++)
					 if(n[i]) xmoyentran[j][i]/=n[i];  //*********** xcumul

			 	    }
      */
	   // printf("Analyse...\n");

   	   if(AnalyserXor()==TRUE_KEY)success_rate[p]++; //***********

		//if(AnalyserXor()==200)success_rate[p]++;

	   //printf("guessed_key=%d\n",AnalyserXor());

	  /*for (i=0;i<256;i++){for (j=0;j<NBMESURES;j++)
						xmoyentran[j][i]*=n[i]; //***********

			 	    }*/
	  }
//fin for (p

/*
 	for (i=0;i<NB_ITE_ATTACK;i++)
			     {
				printf("template_NBM=%d %d  \n",i, success_rate[i]);
			     }
	scanf("%*c");
*/

 }
//fin for(r

//printf("\n\n\n\n");

 fprintf(template_succes_rate_file,"\n\n\n# NBMESURES=%d\n",NBMESURES);
 fprintf(template_succes_rate_file,"# round:%d\n\n\n\n",r-1);

 for(i=0; i<NB_ITE_ATTACK; i++)
 {
   printf("template_NBM=%d %d      \n",i, success_rate[i]);
//   fprintf(template_succes_rate_file,"m=%d:%f \n",i,(float)success_rate[i]/(r));
   fprintf(template_succes_rate_file,"%d\t%f\n",i,(float)success_rate[i]/(r));
 }

 fflush(template_succes_rate_file);
}

///////////////////////////////////////////////////////////////
//commencer
// termine le programme


int main( int argc, char ** argv )
{
	if( argc ==2 )
	{
		CENTRAL_POINT = atoi( argv[1] );
	}
  printf("Hello analyse de consommation!\n");

    // initialisation

  open_files();

        printf("Mesures...");

		template_estimation(); //profiling

		printf("\n\nje suis à travailler()\n\n");

        travailler(); //Matching

  printf("\a\a\a\a");

  close_files();

  return 0;
}


