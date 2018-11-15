
/**
 * @file gradienteconjugado.c
 * @author GRR20176143 Cláudio Torres Júnior
 * @author GRR20171607 Gabriela Stein
 * @date 16 Sep 2018
 * @brief Arquivo de implementação das principais funções.
 */

#include "gradienteconjugado.h"
#include "utils.h"

/**
 * @brief Função que multiplica duas matrizes
 * @param pri Primeiro parametro: Matriz
 * @param sec Segundo parametro: Matriz
 * @param mult Resultado da multiplicação
 * @param tam Ordem da Matriz
*/

inline void multMatMat(double *pri, double *sec, long int dgn, long int tam, double *mult){
	long int i, j, k;
	double soma = 0.0;

	long int numZerosJ;	
	long int col;
	long int colJ;
	long int tamJ;
	long int iniJ;
	long int cJ;

	long int cK;
	long int numZerosK;	
	long int salvaZero;	
	long int salvaColK;	
	long int colK;

	col = dgn;
	colK = dgn - dgn/2;
	for(i = 0; i < dgn/2; ++i){
		numZerosK = dgn/2;
		cK = 0;

		for(j = 0; j < col; ++j){
			for(k = cK; (k < colK) && k < (dgn/2 + 1 + j); ++k){
				soma = soma + pri[i*dgn + (dgn/2 - i) + k]*sec[k*dgn + (dgn/2 - k) + j];
				//printf("(%ld, %ld)*(%ld, %ld)\n", i, k, k, j);
			}
			mult[i*(dgn*2 - 1) + ((dgn*2 - 1)/2 - i) + j] = soma;
			//printf("%lf\n", soma);
			soma = 0.0;
			numZerosK--;
			if(numZerosK < 0){
				cK++;
			}
		}
		colK++;
		col++;
	}
	
	col = dgn + dgn/2;
	colJ = dgn;
	numZerosJ = 0;
	tamJ = dgn - dgn/2;
	salvaZero = dgn/2;
	iniJ = 0;
	cK = 0;
	cJ = 0;
	for(i = dgn/2; i < (tam - dgn/2); ++i){
		numZerosK = salvaZero;
		colK = dgn - dgn/2;
		salvaColK = dgn - dgn/2 + iniJ;

		for(j = iniJ; (j < col) && (j < tam); ++j){
			
			if(salvaColK > colJ){
				colK = colJ;
			}else{
				colK = salvaColK;
			}

			for(k = cK; k < colK && k < tam; ++k){
				soma = soma + pri[i*dgn + (dgn/2 - i) + k]*sec[k*dgn + (dgn/2 - k) + j];
				//printf("(%ld, %ld)*(%ld, %ld)-->cK %ld, numZerosK %ld \n", i, k, k, j, cK, numZerosK);
			}
			mult[i*(dgn*2 - 1) + ((dgn*2 - 1)/2 - i) + j] = soma;
			//printf("%lf\n", soma);
			soma = 0.0;
			numZerosK--;
			colK++;
			tamJ++;
			salvaColK++;
			if(numZerosK < 0){
				if(abs(numZerosK) > abs(numZerosJ)){
					cK = abs(numZerosK);
				}
			}

		}
		col++;
		colJ++;
		numZerosJ--;
		tamJ = dgn - dgn/2;
		if(numZerosJ < 0){
			if(abs(numZerosJ) > dgn/2){
				iniJ++;
				tamJ += iniJ;
				salvaZero--;
			}	
			cJ++;
			if(abs(numZerosJ) >= abs(numZerosK)){
				cK = abs(numZerosJ);
			}else{
				cK = cJ;
			}
		}else{
			cK = 0;
		}
	}

	colK = 1;
	for(i = tam - dgn/2; i < tam; ++i){
		numZerosK = dgn - 1;
		cK = tam - dgn + 1;
		colK = cK;
		colJ = tam - (dgn + 2);

		for(j = colJ; j < tam; ++j){
			for(k = cK; (k < colK + 1) && (k < tam); ++k){
				soma = soma + pri[i*dgn + (dgn/2 - i) + k]*sec[k*dgn + (dgn/2 - k) + j];
				//printf("(%ld, %ld)*(%ld, %ld)\n", i, k, k, j);
			}
			mult[i*(dgn*2 - 1) + ((dgn*2 - 1)/2 - i) + j] = soma;
			//printf("%lf\n", soma);
			soma = 0.0;
			numZerosK--;
			colK++;
			if(numZerosK < 0){
				cK++;
			}
		}
		colJ++;
	}
}

/**
 * @brief Função que multiplica uma matriz por um vetor
 * @param pri Primeiro parametro: Matriz
 * @param sec Segundo parametro: Vetor
 * @param mult Resultado da multiplicação
 * @param tam Ordem da Matriz
*/

inline void multMatVet(double *pri, double *sec, long int inicio, long int dgn, long int tam, double *mult){
	long int i, j, c = 0;
	long int numZeros = dgn/2;
	long int col = dgn - dgn/2;	
	double soma = 0.0;

	for(i = 0; i < tam; i++){
		for(j = c; (j < col) && (j < tam); j++){
			soma = soma + pri[i*dgn + (dgn/2 - i) + j]*sec[inicio + j];
			//printf("(%ld, %ld)*(%ld + %ld)\n", i, j, inicio, j);
		}
		mult[inicio + i] = soma;
		soma = 0.0;
		col++;
		numZeros--;
		if(numZeros < 0){
			c++;
		}
	} 
}

/**
 * @brief Função que multiplica dois vetores
 * @param pri Primeiro parametro: Vetor
 * @param sec Segundo parametro: Vetor
 * @param tam Ordem da Matriz
*/

inline double multVetVet(double *pri, double *sec, long int dgn, long int tam){
	long int i;
	double soma1, soma2, soma3, soma4, soma5;
	long int numZeros = dgn/2;

	soma1 = soma2 = soma3 = soma4 = soma5 = 0.0;
	for(i = 0; i < tam; i+=4){
		soma1 = soma1 + pri[i]*sec[i];
		soma2 = soma2 + pri[i+1]*sec[i+1]; 
		soma3 = soma3 + pri[i+2]*sec[i+2]; 
		soma4 = soma4 + pri[i+3]*sec[i+3]; 
	}
	for(i = tam; i < (tam + numZeros); ++i){
		soma5 = soma5 + pri[i]*sec[i]; 
	}
	return (soma1 + soma2 + soma3 + soma4 + soma5);
}

/**
 * @brief Função que transforma um sistema em um sistema simetrico e positivo definido
 * @param A matriz original
 * @param B vetor de termos independentes
 * @param tam Ordem da Matriz
*/

inline void trasformaSistema(double *A, double *B, double *Atf, double *Btf, parametro par){
	long int tam = par.n*par.k;

	double *T = (double*)malloc(tam*sizeof(double));

	memset(T, 0, sizeof(*T));

	transposta(A, T, par);
	multMatMat(T, A, par.k, par.n, Atf);
	multMatVet(T, B, par.k/2, par.k, par.n, Btf);

	free(T);
}

/**
 * @brief Função que calcula a transposta de uma matriz
 * @param A matriz original
 * @param T matriz transposta
 * @param tam Ordem da Matriz
*/

inline void transposta(double *A, double *T, parametro par){
	long int i, j;
	long int numZeros = par.k/2;
	long int col = par.k - par.k/2;
	long int c = 0;

	for(i = 0; i < par.n; i++){
		for(j = c; (j < col) && (j < par.n); j++){
			T[j*par.k + (par.k/2 - j) + i] = A[i*par.k + (par.k/2 - i) + j];
		}
		col++;
		numZeros--;
		if(numZeros < 0){
			c++;
		}
	}
}

/**
 * @brief Função que condiciona uma matriz
 * @param A matriz original
 * @param M matriz resultante do uso do pré condicionador
 * @param p Indica o pré-condicionador a ser utilizado
 * @param tam Ordem da Matriz
*/

inline void preCondicionador(double p, double *M, double *A, parametro par){
	long int i, j;

	if(p == 0.0){//sem pre-condicionador
		for(i = 0; i < par.n; i++){
			for(j = 0; j < par.n; j++){
				if(i == j){
					M[i*par.n + j] = 1.0;
				}else{
					M[i*par.n + j] = 0.0;
				}
			}
		}
	}else if(p > 0.0 && p < 1.0){//pre-condicionador de jacobi
		long int numZeros = par.k/2;
		long int passo = (par.k*2 - 1)/2;

		for(i = par.k/2; i < (par.n + numZeros); ++i){
			M[i] = A[passo];
			passo += par.k*2 - 1;
		}

	}else{//pré-condicionador de Gauss-Seidel p/ x = 1 e pré-condicionador SSOR p/ 1 < x < 2 
		double *L = (double*)malloc(par.n*par.n*sizeof(double)); //matriz lower
		double *U = (double*)malloc(par.n*par.n*sizeof(double)); //matriz upper
		double *D = (double*)malloc(par.n*par.n*sizeof(double)); //matriz diagonal
		double *R1 = (double*)malloc(par.n*par.n*sizeof(double)); //matriz resposta1
		double *R2 = (double*)malloc(par.n*par.n*sizeof(double)); //matriz resposta2

		criaMatrizes(A, L, U, D, par.n);
		for(i = 0; i < par.n; i++){
			for(j = 0; j < par.n; j++){
					R1[i*par.n + j] = D[i*par.n + j] + p*L[i*par.n + j];
			}
		}
		for(i = 0; i < par.n; i++){
			for(j = 0; j < par.n; j++){
					R2[i*par.n + j] = D[i*par.n + j] + p*U[i*par.n + j];
			}
		}
		for(i = 0; i < par.n; i++){
			for(j = 0; j < par.n; j++){
				if(i == j){
					D[i*par.n + j] = 1 / D[i*par.n + j];
				}
			}
		}
		multMatMat(R1, D, par.k, par.n, U); //matriz U sendo usada como intermediaria da operação
		multMatMat(U, R2,par.k, par.n, M);
		Cholesky(M, par.n);
		free(L);
		free(U);
		free(D);
		free(R1);
		free(R2);
	}
}

/**
 * @brief Função que encontra o maior valor de um vetor
 * @param V Vetor qualquer
 * @param tam Ordem da Matriz
 * @return Retorna o maior valor
*/

inline double maxVetor(double *V, parametro par)
{	
	long int i;
	long int numZeros = par.k/2;
	double maximo = V[0];
	for (i = 0; i < (par.n + numZeros); ++i)
	{
		if (V[i] > maximo)
			maximo = V[i];
	}
	return maximo;
}

/**
 * @brief Função que imprime a saída de dados
 * @param erroIT Vetor que contem o erro em todas iterações
 * @param X Vetor solução
 * @param norma Norma do resíduo
 * @param pc Tempo para cálculo do pré-condicionador
 * @param it Tempo para resolver uma iteração
 * @param r Tempo usado para calcular o residuo do SL
 * @param par Struct contendo os parametros passados na execução do programa
 * @param iter Quantidade de iterações usadas pelo método
*/

inline void imprime_dados(double *erroIt, double *X, double norma, double pc, double it, double r, parametro par, long int iter)
{	
	long int numZeros = par.k/2;
	FILE *arqOut;
	arqOut = fopen(par.o, "w");
	if(arqOut == NULL){
		fprintf(stderr, "erro ao criar arquivo\n");
		exit(1);
	}

	fprintf(arqOut, "# ctj17 Cláudio Torres Júnior\n"); //# login1 Nome1
	fprintf(arqOut, "# gs17 Gabriela Stein\n#\n"); //# login2 Nome2
	for(long int i = 0; i < iter; i++){
		fprintf(arqOut, "# iter %ld: <||%lf||>\n", i+1, erroIt[i]); //# iter k: <||x||>
	}
	fprintf(arqOut, "# residuo: <||%lf||>\n", norma); //# residuo: <||r||>
	fprintf(arqOut, "# Tempo PC: <%.15g>\n", pc); //# Tempo PC: <tempo para cálculo do pré-condicionador>
	fprintf(arqOut, "# Tempo iter: <%.15g>\n", it); //# Tempo iter: <tempo para resolver uma iteração do método>
	fprintf(arqOut, "# Tempo residuo: <%.15g>\n#\n%ld\n", r, par.n); //# Tempo residuo: <tempo para calcular o residuo do SL> 
	for(long int i = par.k/2; i < (par.n + numZeros); i++){
		fprintf(arqOut, "%.15g ", X[i]); //x_1 x_12 ... x_n
	}
	fprintf(arqOut, "\n");
} 

/**
 * @brief Função que calcula a Fatoração Incompleta de Cholesky
 * @param M Matriz condionada usando Gauss Seidel ou SSOR
 * @param tam Ordem da Matriz
*/

inline void Cholesky(double *M, long int tam)
{	
	long int i, j, k;
	
	for(k = 0; k < tam; k++){
		M[k*tam + k] = sqrtf(M[k*tam + k]);
		for(i = k + 1; i < tam; i++){
			if(M[i*tam + k] != 0.0){
				M[i*tam + k] = M[i*tam + k] / M[k*tam + k];
			}
		}
		for(j = k + 1; j < tam; j++){
			for(i = j; i < tam; i++){
				if(M[i*tam + j] != 0.0){
					M[i*tam + j] = M[i*tam + j] - M[i*tam + k]*M[j*tam + k];
				}
			}
		}
	}
	for(i = 0; i < tam; i++){
		for(j = i + 1; j < tam; j++){
			M[i*tam + j] = 0.0;
		}
	}	
}

/**
 * @brief Função que aloca e atribui as matrizes criadas para os pré condionantes Gauss Seidel e SSOR
 * @param A Matriz original
 * @param L Matriz inferior (lower)
 * @param U Matriz superior (upper)
 * @param D Matriz diagonal
 * @param tam Ordem da Matriz
*/

inline void criaMatrizes(double *A, double *L, double *U, double *D, long int tam)
{
	long int i, j;
	//cria a matriz upper
		for(i = 0; i < tam; i++){
			for(j = 0; j < tam; j++){
				if(j > i){
					U[i*tam + j] = A[i*tam + j];
				}else{
					U[i*tam + j] = 0.0;
				}

				if(j < i){
					L[i*tam + j] = A[i*tam + j];
				}else{
					L[i*tam + j] = 0.0;
				}

				if(j == i){
					D[i*tam + j] = A[i*tam + j];
				}else{
					D[i*tam + j] = 0.0;
				}
			}
		}
}

/**
 * @brief Função que libera os vetores criados durante a execução do programa
 * @param M Matriz condicionada
 * @param X Vetor solução
 * @param Xant Vetor solução da iteração anterior
 * @param r Resíduo
 * @param v M^(-1)*B
 * @param z Resultado para a matriz verdadeira * vetor resposta de M
 * @param y vetor soluçao para o sistema M*y = r
 * @param T matriz transposta
 * @param erroAproximadoR vetor para o erro aproximado relativo
 * @param erroAproximadoA vetor para o erro aproximado absoluto
 * @param erroIt vetor para o erro maximo de cada iteração
*/

inline void liberaVet(double *M, double *X, double *r, double *v, double *z, 
	double *y, double *Xant, double *erroAproximadoA, double *erroIt, double *Atf, double *Btf){

	free(M);
	free(X);
	free(r);
	free(v);
	free(z);
	free(y);
	free(Xant);
	free(erroAproximadoA);
	free(erroIt);
	free(Atf);
	free(Btf);
}

/**
 * @brief Função que calcula a solução do sistema linear pelo método do Gradiente Conjugado
 * @param A Matriz do sistema 
 * @param B Vetor solução
 * @param par Struct contendo os parametros passados durante a chamada do programa
*/

int gradienteConjugado(double *A, double *B, parametro par){

	LIKWID_MARKER_INIT;

	int convergiu = 0;
	long int numZeros = par.k/2;

	double *M = (double*)malloc((par.n + numZeros)*sizeof(double)); //pre-condicionador
	double *X = (double*)malloc((par.n + numZeros)*sizeof(double)); 		//vetor de 'chutes' iniciais
	double *Xant = (double*)malloc((par.n + numZeros)*sizeof(double));
	double *r = (double*)malloc((par.n + numZeros)*sizeof(double));		//residuo
	double *v = (double*)malloc((par.n + numZeros)*sizeof(double));		
	double *z = (double*)malloc((par.n + numZeros)*sizeof(double));		
	double *y = (double*)malloc((par.n + numZeros)*sizeof(double));	
	////double *erroAproximadoR = malloc(par.n*sizeof(double));           //vetor com o erro relativo
	double *erroAproximadoA = malloc((par.n + numZeros)*sizeof(double));	//vetor com o erro absoluto
	double *erroIt = malloc(par.i*sizeof(double));	//vetor com o erro max absoluto em cada iteraçao
	double *Atf = (double*)malloc((par.k*2 - 1)*par.n*sizeof(double)); 
	double *Btf = (double*)malloc((par.n + numZeros)*sizeof(double)); 		

	double s, aux, aux1, m, norma; //Xprox;

	long int i, j, it;

	memset(Atf, 0, sizeof(*Atf));
	memset(M, 0, sizeof(*M));

	tempo t_pc; //struct para usar o timestamp precodicionador
	tempo t_it; //struct para usar o timestamp cada iteraçao
	tempo t_r; //struct para usar o timestamp residuo

	t_it.ini = timestamp();
	t_pc.ini = timestamp();

	LIKWID_MARKER_START("OP1"); //iteraçao

	//transforma A em uma matriz positiva simetrica
	trasformaSistema(A, B, Atf, Btf, par);
	//acha pre-condicionador
	preCondicionador(par.p, M, Atf, par);

	t_pc.fim = timestamp();
	t_pc.dif = t_pc.fim - t_pc.ini;

	
	for(i = 0; i < (par.n + numZeros); ++i){
		X[i] = 0; //X0 = 0
		r[i] = B[i]; //r = B
	}
		

	if(par.p < 1){
//		//v = M-¹ * B e y = M-¹ * r pois r == B
		for(i = par.k/2; i < (par.n + numZeros); ++i){
			v[i] = Btf[i] / M[i];
			y[i] = v[i];
		}
	}else{
		double soma = 0.0;
		v[0] = Btf[0] / M[0*par.n + 0];
		y[0] = v[0];
		for(i = 1; i < par.n; i++){
			soma = Btf[i];
			for(j = i - 1; j >= 0; j--){
				soma = soma - M[i*par.n + j]*v[j];
			}
			v[i] = soma / M[i*par.n + i];
			y[i] = v[i];
		}
	}

	//aux = y^t * r
	aux = multVetVet(y, r, par.k, par.n);

	t_it.fim = timestamp();
	t_it.dif = t_it.fim - t_it.ini;
	//it == 1 pois a it 0 foi feito fora do for
	for(it = 0; it < par.i; it++){
		t_it.ini = timestamp();

		//z = A*v
		multMatVet(Atf, v, par.k/2, (par.k*2 - 1), par.n, z);

		s = aux/multVetVet(v, z, par.k, par.n);

		//salva o vetor
		//for(i = 0; i < (par.n + numZeros); ++i)
		//	Xant[i] = X[i];

		//x = x + s*v
		//erro aproximado absoluto
		for(i = 0; i < par.n; i+=4){
			erroAproximadoA[i] = fabs(X[i] - (X[i] + s*v[i])); 
			erroAproximadoA[i+1] = fabs(X[i+1] - (X[i+1] + s*v[i+1])); 
			erroAproximadoA[i+2] = fabs(X[i+2] - (X[i+2] + s*v[i+2])); 
			erroAproximadoA[i+3] = fabs(X[i+3] - (X[i+3] + s*v[i+3])); 
			X[i] = X[i] + s*v[i];
			X[i+1] = X[i+1] + s*v[i+1];
			X[i+2] = X[i+2] + s*v[i+2];
			X[i+3] = X[i+3] + s*v[i+3];
		}
		for(i = par.n; i < (par.n + numZeros); ++i){
			erroAproximadoA[i] = fabs(X[i] - (X[i] + s*v[i])); 
			X[i] = X[i] + s*v[i];
		}


		erroIt[it] = maxVetor(erroAproximadoA, par);

		//r = r - s*z
		LIKWID_MARKER_START("OP2");
		t_r.ini = timestamp();
		//calculo do residuo
		for(i = 0; i < par.n; i+=4){
			r[i] = r[i] - s*z[i]; 
			r[i+1] = r[i+1] - s*z[i+1]; 
			r[i+2] = r[i+2] - s*z[i+2]; 
			r[i+3] = r[i+3] - s*z[i+3]; 
		}
		for(i = par.n; i < (par.n + numZeros); ++i){
			r[i] = r[i] - s*z[i]; 
		}
		t_r.fim = timestamp();
		t_r.dif = t_r.fim - t_r.ini;

		LIKWID_MARKER_STOP("OP2");


		//y = M-¹ * r
		if(par.p < 1){
			for(i = par.k/2; i < (par.n + numZeros); ++i){
				y[i] = r[i] / M[i];
			}
		}else{
			double soma = 0.0;
			y[0] = r[0] / M[0*par.n + 0];
			for(i = 1; i < par.n; i++){
				soma = r[i];
				for(j = i - 1; j >= 0; j--){
					soma = soma - M[i*par.n + j]*y[j];
				}
				y[i] = soma / M[i*par.n + i];
			}
		}

		if(maxVetor(erroAproximadoA, par) < par.e && !convergiu){
			//achou resultado
			//t_r.ini = timestamp();
			//calcula a norma
			norma = sqrtf(multVetVet(r, r, par.k, par.n)); //norma euclidiana do residuo
			//t_r.fim = timestamp();
			//t_r.dif = t_r.fim - t_r.ini;
			//imprime dados no arquivo
			imprime_dados(erroIt, X, norma, t_pc.dif/1000, t_it.dif/(it + 1)/1000, t_r.dif/1000, par, it);
			/*for(i = 0; i < par.n; i++)
				printf("%lf ", X[i]);
			printf("\n ");*/
			if (par.op == 0){ //se falso, para as iterações
				liberaVet(M, X, r, v, z, y, Xant, erroAproximadoA, erroIt, Atf, Btf); 
				return 0;
			}else{
				convergiu = 1;
			}
		}

		//aux1 = y^t * r;
		aux1 = multVetVet(y, r, par.k, par.n);

		//m = aux1/aux
		m = aux1/aux;

		//aux = aux1
		aux = aux1;

		//v = y + m*v
		for(i = 0; i < par.n; i+=4){
			v[i] = y[i] - m*v[i]; 
			v[i+1] = y[i+1] - m*v[i+1]; 
			v[i+2] = y[i+2] - m*v[i+2]; 
			v[i+3] = y[i+3] - m*v[i+3]; 
		}
		for(i = par.n; i < (par.n + numZeros); ++i){
			v[i] = y[i] - m*v[i]; 
		}

	t_it.fim = timestamp();
	t_it.dif = t_it.dif + (t_it.fim - t_it.ini);
	}
	if (!convergiu){		
		norma = sqrtf(multVetVet(r, r, par.k, par.n)); //norma euclidiana do residuo

		imprime_dados(erroIt, X, norma, t_pc.dif/1000, t_it.dif/(it + 1)/1000, t_r.dif/1000, par, it);
		//fprintf(stderr, "O método não convergiu!\n");
	}

	liberaVet(M, X, r, v, z, y, Xant, erroAproximadoA, erroIt, Atf, Btf); 
	LIKWID_MARKER_STOP("OP1");

  	LIKWID_MARKER_CLOSE;
	return -1; 
}
