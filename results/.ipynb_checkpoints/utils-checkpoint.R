### 1.1 Generación de índices
indices <- function(n) {
  # Crea una lista de tamaño (n-1)n/2 con pares de índices de la siguiente
  #  manera: (1,2),..,(1,n),(2,3),..,(2,n),...,(n-1,n)
  # Args: 
  #    n (int): número entero postivo 
  #       se refiere al número de columnas de una matriz
  #Returns:
  #    lista con pares de índices
  a <- NULL
  b <- NULL
  indices <- NULL
  for (i in 1:(n-1)){
    a <- append(a,rep(i,n-i))
    b <- append(b,seq(i+1,n))    
  }
  for(i in 1:round(n*(n-1)/2))
    indices[[i]] <- list(c(a[i], b[i]))
  indices
}

### 1.2 Verificación de ortogonalidad entre vectores
ortogonal <- function(u,v,TOL=10^-8){
  # Verifica si dos vectores son ortogonales, de acuerdo a cierto nivel de tolerancia, 
  # arrojando un 1 si lo es, y un 0 si no lo es.
  # Args: 
  #   u (vector): vector de dimension n,
  #   v (vector): vector de dimension n, 
  #   TOL (numeric): real positivo, que sirve como parametro de tolerancia para evaluar ortogonalidad de u y v. 
  #   Notas: 
  #   1) Valor por default TOL es 10^-8
  #   2) Se sugiere una TOL mayor a 10^-32.
  # Returns: 
  #   Valor booleano 0 (no son ortongoales), 1 (son ortogonales)
  if ( norm(u,type ="2") < TOL | norm(v,type ="2") < TOL){ret<-0} 
  else{ 
    if( (u%*%v)/(norm(u,type ="2")*norm(v,type ="2")) < TOL){ret<-1}
    else{ret<-0}  
  }
  ret
}

### 1.3 Función signo
signo<-function(x) {
  # Indica el signo de un número x
  # Args: 
  #    x (numeric): número a revisar
  # Returns:
  #    1 si el número es positivo o cero
  #    -1 si el número es negativo
  
  ifelse(x<0,-1,1)
}

### 1.4 Solver dada descomposición SVD
solver <- function(U,S,V,b){
  # Construye la solución de un sistema de ecuaciones a partir de matrices 
  # U, S, V, y vector b. Se asume que S es diagonal. 
  # Para ello resuelve S d = U^Tb, y construye x=Vd.
  # Notas:
  # 1) Se utilizó la función backsolve para resolver el sistema triangular.
  # 2) Al ser S diagonal, es indistinto si usar un solver para matrices traingulares inferiores o superiores.
  # Args: 
  #     U (matriz): matriz para lado derecho de sistema S d = U^Tb, con entrada reales y dimension m x n,
  #     S (matriz): matriz diagonal, que define sistema sistema S d = U^Tb, con entrada reales y dimension n x n,
  #     V (matriz): para construir x, con entrada reales y dimension n x n, 
  #     b (vector): vector con el que se forma lado derecho de primer sistema, de dimension m.
  # Returns: 
  #     x (vector): vector formado por la solucion de S d = U^Tb, multiplicado por V, con dimension n
  d = backsolve(S, t(U)%*%b)
  x = V%*%d
  return(x)
}







