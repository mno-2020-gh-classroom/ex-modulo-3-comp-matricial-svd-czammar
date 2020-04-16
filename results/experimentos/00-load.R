## 2.1 Algoritmo One-sided Jacobi para aproximar descomposicion SVD

svd_jacobi_aprox <- function(A,TOL,maxsweep){
  # Calcula la descomposición de una matriz A en sus componentes U, S V, 
  # utilizando el método de Jacobi para calcular la factorización SVD.De esta forma 
  # la matriz A queda descompuesta de la siguiente forma: A = U*S*t(V).
  # Args: 
  #    A (matriz): Matriz de entrada (nxm) de números reales a la que se le calculará la descomposición SVD.
  #    TOL (numeric): controla la convergencia del método, siendo un valor real de 10^-8 (sugerido en la nota 3.3.d.SVD)
  #    Nota: Se sugiere una TOL mayor a 10^-32.
  #    maxsweep (numeric): número máximo de sweeps,donde cada sweep consiste de un número máximo(nmax)
  #    de rotaciones; y en cada sweep se ortogonalizan 2 columnas.
  # Returns: 
  #   Lista con 3 elementos, donde el primer elemento representa a la matriz S(mxm) matriz diagonal,el segundo a la matriz U(nxm)
  #   y el tercero y último a la matriz V (mxm).En conjunto estas tres matrices componen la factorización SVD de la matriz de entrada A.
  # Nota: Esta función estima la SVD thin,la cual calcula unicamente las m columnas de U correspondientes a los m renglones de V. De esta
  # manera las columnas restantes de U no son calculadas, provocando una mejora significativa en velocidad de ejecución comparada con la 
  # la Full SVD. Referencia: https://en.wikipedia.org/wiki/Singular_value_decomposition#Thin_SVD.
  
  #dimensiones
  n<-dim(A)[2] #numero de columnas
  m<-dim(A)[1] #numero de filas
  nmax<-n*(n-1)/2
  
  #inicializa valores del ciclo
  ak<-A
  vk<-diag(n)
  sig <- NULL
  uk <- ak
  num_col_ortogonal<-0
  k<-0
  stop<-FALSE
  
  while(k<=maxsweep & num_col_ortogonal<nmax){
    num_col_ortogonal<-0
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        col_j<-ak[,j]
        col_i<-ak[,i]
        
        #comprueba ortogonalidad  
        if(ortogonal(col_i,col_j,TOL)==1){
          num_col_ortogonal<-num_col_ortogonal+1
        }
        else{
          #calcula coeficientes de la matriz
          a<-col_i%*%col_i
          b<-col_j%*%col_j
          c<-col_i%*%col_j
          
          #si c es cercano a cero no actualiza
          if(c<TOL){
            stop<-TRUE
            break}
          
          #calcula la rotacion givens que diagonaliza
          epsilon<-(b-a)/(2*c)
          t<-signo(epsilon)/(abs(epsilon)+sqrt(1+epsilon**2))
          cs<-1/sqrt(1+t**2)
          sn<-cs*t
          
          #actualiza las columnas de la matriz ak
          temp<-ak[,i] 
          ak[,i]<-c(cs)*temp-c(sn)*ak[,j]
          ak[,j]<-c(sn)*temp+c(cs)*ak[,j]
          
          
          #actualiza las columnas de la matriz vk
          temp<-vk[,i] #cambio
          vk[,i]<-c(cs)*temp-c(sn)*vk[,j]
          vk[,j]<-c(sn)*temp+c(cs)*vk[,j]             
        }#cierra else
      }#cierra for j
      if(stop==TRUE){
        stop<-FALSE
        break
      }
    }#cierra for i
    k<-k+1
  }#cierra while
  
  #Obtener sigma
  sig<-apply(ak, 2, function(x){norm(x,"2")})
  
  #Obtener U
  for(i in 1:n){
    if (sig[i]<TOL){
      uk[,i]<-0  
    } else{
      uk[,i] <- ak[,i]/sig[i]
    }
  }
  
  # Indices de sigma ordenada en forma decreciente para ordenar V,S,U
  index <- order(sig,decreasing = TRUE)
  vk <- vk[,index]
  S <- diag(sig[index])
  uk <- uk[,index]
  
  list(S = S, U = uk, V= vk)
}

## 2.2 Linear solver aproximating SVD decomposition using One-sided Jacobi algorithm

sel_solver<-function(A,b,TOL=10**-8,maxsweep){
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