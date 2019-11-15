ISTRUZIONI COMPILAZIONE E UTILIZZO

PREREQUISITI DA INSTALLARE
make	(per compilare)
gcc o altro compilatore c compatibile
gmp (libreria aritmetica multiprecisione)

PASSI PER COMPILARE E ESEGUIRE (I COMANDI VANNO LANCIATI DALLA DIRECTORY PRINCIPALE DEL PROGETTO)
- compilare con make install
- eseguire con ./start.o N B1 [B2] [MAX_ITER]
						B2 può essere omesso in tal caso si userà		B2 = B1 * 100
						MAX_ITER può essere omesso in tal caso si userà		MAX_ITER = 2000

CONFIGURAZIONE DELLE MODALITA' DI ESECUZIONE
- aprire con un editor di testo il file m_ellcurv_fact.h
- all'inizio del file dovrebbo essere presenti (potrebbero essere commentati con //)
#define NO_LINUX =>	disattiva funzioni native su linux 
#define NO_THREAD =>	disattiva il multithread 
- queste modalità possono essere utilizzate qualora il programma abbia qualche problema altrimenti è preferibile lasciarle commentate (precedute da //)

#define THREAD_NUM => numero thread (deve essere > 0) (in caso di attivazione modalità NO_THREAD non c'è bisogno di modificarlo verrà ignorato dal programma)
