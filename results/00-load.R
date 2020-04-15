## 2.1 Algoritmo One-sided Jacobi para aproximar descomposicion SVD

svd_jacobi_aprox <- function(A,TOL,maxsweep){
  # Solucion del sistema A11 y = b1, usando algoritmo One-Sided Jacobi
  svd <-svd_jacobi_aprox(bloq$A11,TOL,maxsweep)
  U <-svd$U
  S <- svd$S 
  V<-svd$V
  y = solver(U,S,V,bloq$b1)
  
  # Solucion del sistema A11 Y = A12 
  Y= bloq$A12
  
  Y = solver(U,S,V,bloq$A12)
  
  # Construye matrix S = A22 A11^-1 A12
  S = bloq$A22 - bloq$A21%*%Y
  
  # Construye b_hat = b2 - A21 A11^-1 b1
  b_hat = bloq$b2-bloq$A21%*%y
  
  # Solucion del sistema S x2 = b_hat algoritmo One-Sided Jacobi
  x_2 = sel_solver(S,b_hat,TOL, maxsweep)
  
  # Construye b_hat2 S x2 = b1 - A12 b_hat
  b_hat2 = bloq$b1 -bloq$A12%*%x_2
  
  # Solucion del sistema A11 x2 = b_hat algoritmo One-Sided Jacobi
  x_1 = solver(U,S,V,b_hat2)
  
  # Solucion del sistema original
  x <- c(x_1,x_2)
  x
}

## 2.2 Solver de sistemas lineales usando aproximacion SVD
sel_solver<-function(A,b,TOL=10**-8,maxsweep=20){
  # Aproxima la solución de un sistema de ecuaciones lineales (SEL) de la forma Ax=b 
  # usando la descomposición SVD de A obtenido por medio del método de One-sided Jacobi.
  # Args: 
  #    A (matriz): matriz de coeficientes del SEL, de números reales y de dimension m x n.
  #    b (vector): vector de lado derecho del sistema, de números reales y de dimension m.
  #    TOL (numeric): real positivo, que sirve para controlar la convergencia del metodo iterativo. 
  #             Nota: se usa un valor default de 10^-8 como tolerancia
  #    maxsweep (int): numero entero posiitivo, que se usar para controlar el número máximo de sweeps
  #                    que se ejecutan como parte del metodo iterativo.
  # Returns: 
  #    x (vector): vector solución del SEL, con dimension n.
  
  svd<-svd_jacobi_aprox(A,TOL,maxsweep)
  x<-solver(svd$U,svd$S,svd$V,b)
  x
}

##  2.3 Eliminación por bloques basada en solver linear que usa SVD

bloques<-function(A,b,corte) {
  # Dadas matriz A de n x ny un vector b de dimension n, así como un parametro corte (de dimension), 
  # genera una descomposicion en 4 y 2 bloques de estos.
  # En concreto, el parametro corte sirve para dividir A en los bloques de las siguientes dimensiones
  #    A11: bloque dado por A[1:corte,1:corte] (cuadrante superior izquierdo de A)
  #    A12: bloque dado por A[1:corte, corte+1:n] (cuadrante superior derecho de A)
  #    A21: bloque dado por A[corte+1:n,1:corte] (cuadrante inferior izquierdo de A)
  #    A22: bloque dado por A[corte+1:n,corte+1:n] (cuadrante inferior derecho de A)
  # Asimismo, b se divide como:
  #    b1: primeras corte-entradas de b
  #    b2: entradas restantes de b
  
  # Args: 
  #    A (matriz): matriz de dimensiones n x n sobre la que se hara la division en bloques, 
  #    b (vector): vector de dimensiones n sobre el cual que se hara la division en bloques,
  #    corte (int): entero positivo, que se usa como parametro de dimension con respecto al cual se hace la division de A (valor entre 1 y n)
  #Returns: 
  #     lista con la matriz A dividida en 4 bloques (lista): A11,A12, A21, A22 y el vector b dividido en 2 bloques b1, b2
  #     Notas: dimensiones de los bloques 
  #       A11 (corte x corte), A12 (corte x n-corte), A11 (n-corte x corte), A22 (n-corte x n-corte)
  #       b1 (corte) y b2 (n-corte)
  
  
  # Obtiene la descomposicion en bloques de A, segun parametro corte
  a11 <- A[c(1:corte),c(1:corte)]
  a12 <- A[c(1:corte),c((corte+1):dim(A)[2])]
  a21 <- A[c((corte+1):dim(A)[1]),c(1:corte)]
  a22 <- A[c((corte+1):dim(A)[1]),c((corte+1):dim(A)[2])]
  
  # Obtiene la descomposicion en bloques de b, segun parametro corte
  b1 <- b[1:(corte)]
  b2 <- b[((corte)+1):length(b)]
  
  list(A11 = a11,A12 =a12, A21=a21, A22=a22, b1 = b1, b2 = b2)
}

## 2.4 Solver de sistema lineal por eliminacion por bloques, basado en SVD

eliminacion_bloques <- function(A,b,corte, TOL, maxsweep){
  # Aproxima la solución de un sistema de ecuaciones lineales (SEL) de la forma Ax=b, usando el método
  # de eliminación por bloques. La solución de los sistemas lineales inducidos por el complemento de Schur
  # se aproxima usando la descomposicion SVD aproximada con el algoritmo One-Sided Jacobi
  #Args: 
  #    A (matriz):  matriz de dimensiones n x n, asociada al sistema del sistema Ax=b que se dividira en bloques para su solucion
  #          Nota: para que exista solucion A y el el bloque A11 deben ser no singulares
  #    b (vector): vector de dimension n, asociado al sistema del sistema Ax=b que se dividira en bloques para su solucion
  #    corte (int): entero positivo, que se usara como parametro para divir en bloques a A y B, (valor entre 1 y n)
  #    TOL (numeric): real positivo, que sirve para controlar la convergencia del algoritmo One-Sided Jacobi para estimar SVD
  #    maxsweep (int): numero entero positivo, que sirve como parametro del máximo de sweeps del algoritmo One-Sided Jacobi
  #                    para estimar SVD (valor entre 1 y (n-1)*n/2)
  #Returns: 
  #    x: vector de dimensiones n x n
  
  # Descomposicion en bloques A11, A12, A21, A22, b1 y b2, segun parametro de dimension
  bloq = bloques(A,b,corte)
  
  # Evalua la singularidad de A y el bloque A11, en caso de que una de tales sea singular
  # arroja mensaje al usuario, en caso contrario intenta aproximar sistema inducidos por 
  # complemento de Schur
  if(is.singular.matrix(A,TOL) ==FALSE & is.singular.matrix(bloq$A11,TOL) == FALSE){
    
    # Solucion del sistema A11 y = b1, usando algoritmo One-Sided Jacobi
    y = sel_solver(bloq$A11,bloq$b1,TOL, maxsweep)
    
    # Solucion del sistema A11 Y = A12 
    Y= bloq$A12
    for(i in 1:dim(bloq$A12)[2]){
      Y[,i] = sel_solver(bloq$A11,bloq$A12[,i],TOL, maxsweep)
    }
    
    # Construye matrix S = A22 A11^-1 A12
    S = bloq$A22 - bloq$A21%*%Y
    
    # Construye b_hat = b2 - A21 A11^-1 b1
    b_hat = bloq$b2-bloq$A21%*%y
    
    # Solucion del sistema S x2 = b_hat algoritmo One-Sided Jacobi
    x_2 = sel_solver(S,b_hat,TOL, maxsweep)
    
    # Construye b_hat2 S x2 = b1 - A12 b_hat
    b_hat2 = bloq$b1 -bloq$A12%*%x_2
    
    # Solucion del sistema A11 x2 = b_hat algoritmo One-Sided Jacobi
    x_1 = sel_solver(bloq$A11,b_hat2,TOL, maxsweep)
    
    # Solucion del sistema original
    x <- c(x_1,x_2)
    x
  }else{
    print("Singularidad de las matrices involucradas; no hay solucion")
  }
  
}