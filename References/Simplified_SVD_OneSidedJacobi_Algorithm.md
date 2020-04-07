# Versión simplificada Algoritmo One-Sided Jacobi para SVD

### 0. Introducción

Este documento presenta una versión simplificada a nivel de pseudo-código del algoritmo que los equipos de programación y revisión deben de consolidar a través de código de R.

Para facilitar la implementación y documentación de la funciones que compondra estas labores, se ha divido el proyecto a realizarse en una serie de etapas que se detallan a continuación.

## 1. Funciones auxiliares 

### 1.1 Generación de índices

**Propósito:**

Dado un número entero positivo $n$, generar una función *indices* que devuelva una lista de índices o sublistas  siguiente:

$$\{(1,2), (1,3), \ldots (1,n)\} \cup  \{(2,3), (2,4), \ldots (2,n)\} \cup  \{(3,4),\ldots (2,n)\} \cup \{(n-1,n)\} $$

**Notas:** 

* La cardinalidad del conjunto de índices es 

$$\frac{n(n+1)}{2}$$

* La función debe regresar elementos un elemento que se pueda recorrer y que a su vez se puedan acceder a los elementos de las tuplas. Por ejemplo, si contiene a la tupla (1,3), queremos poder acceder a 1 y 3 posteriormente en la salida de la función.

Por favor, incluir la documentación correspondiente:

```R
indices <- function(...){
    # Breve explicacion del objetivo de la funcion.
    #Args: # Argumentos involucrados en la definicion de la funcion
    #    indexes (function): funcion de ...
    #    argumento (tipo de argumento, por ejemplo int): explicacion 
    #    ...
    #Returns:
    #    ...
  
  # Corpus de la funcion
  
}
```



### 1.2 Verificación de ortogonalidad entre vectores

**Propósito:**

Dados dos vectores $u,v\in \mathbb{R}^n$ y un parámetro $TOL$ de tolerancia, generar una función *ortogonal* que verifique si el producto punto de la versión normalizada de ambos vectores se considera ortogonal; es decir, que evalue la verdad o falsedad de la afirmación:

$$ \frac{u \cdot v^T}{|| u ||_2 || v ||_2 } < TOL$$



**Notas:** 

* $TOL$ debe ser un parámetro (número real positivo) de este funcion que el usuario pueda definiar a voluntad,

* La norma de los vectores a involucrarse en la definición es la norma-2,
* La función debe regresar un valor booleano, que nos diga si dada la tolerancia la versión normalizada de ambos vectores se considera ortogonal.

Por favor, incluir la documentación correspondiente (ver chuck en el numeral 1.1)

```R
ortogonal <- function(...){
    # Breve explicacion del objetivo de la funcion.
    #Args: # Argumentos involucrados en la definicion de la funcion
    #    ...
    #Returns:
    #    ...
  
  # Corpus de la funcion
  
}
```

### 1.3 Función signo

**Propósito:**

Generar una función *signo* que verifique el signo de un número real; es decir:

$$signo(a) = \begin{cases} 1 & \mbox{ si } a \geq 0 \\ -1 & \mbox{ si } a < 0 \end{cases}$$

Por favor, incluir la documentación correspondiente (ver chuck en el numeral 1.1)

```R
signo <- function(...){
    # Breve explicacion del objetivo de la funcion.
    #Args: # Argumentos involucrados en la definicion de la funcion
    #    ...
    #Returns:
    #    ...
  
  # Corpus de la funcion
  
}
```

### 1.4 Solver dada descomposición

**Propósito:**

Dada una matriz $ A \in \mathbb{R}^{m \times n}$ y su descomposición SVD, dado por las consabidas matrices $U \in \mathbb{R}^{m \times m}$, $\Sigma\in \mathbb{R}^{m \times n}$ y $V \in \mathbb{R}^{n \times n}$, así como un vector $b \in \mathbb{R}^m$. Formar una función *solver* que devuelva la soluci´on el sistema lineal $Ax=b$ a través del siguiente procediemiento:

* Resolver el sistema de ecuaciones $$ \Sigma d = U^T b$$, empleando alguna de las librerías [backsolve/forwardsolve](https://stat.ethz.ch/R-manual/R-devel/library/base/html/backsolve.html) de R.
* Formar al vector $x:=Vd$ la cual será solución del sistema $Ax=b$.

**Notas:** 

* La matrices $U$, $\Sigma$, $V$ deben ser parte de las entradas la función.

Por favor, incluir la documentación correspondiente (ver chuck en el numeral 1.1).

```R
solver <- function(...){
    # Breve explicacion del objetivo de la funcion.
    #Args: # Argumentos involucrados en la definicion de la funcion
    #    ...
    #Returns:
    #    ...
  
  # Corpus de la funcion
  
}
```

## 2. Algoritmo SVD y solución de sistema lineal

### 2.1 One-sided Jacobi numerical aproximación

**Propósito:**

Dada una matriz $$A \in \mathbb{R}^{ m \times n}$$, crear una función *SVD_jacobi* que permita obtener aproximaciones numéricas asociadas a la descomposición SVD de esta,. En concreto, si la descomposición de $A$ está dada por las consabidas matrices $U \in \mathbb{R}^{m \times m}$, $\Sigma\in \mathbb{R}^{m \times n}$ y $V \in \mathbb{R}^{n \times n}$,  queremos obtener versiones numéricas de 1) $V$ y 2) $W:=U \Sigma$.

Ello se debe realizar con base en las siguientes especificaciones, que circunscriben el Algoritmo One-sided Jacobi para la descomposición SVD:

**Entradas:**

* $$A \in \mathbb{R}^{ m \times n}$$,
* TOL (numero real, parámetro de tolerancia prefijado en $10^{-8}$),
* *maxsweeps* (entero mayor a cero). **Nota:** Las rotaciones que se aplicarán sucesivamente para ortogonalizar las matrices a partir de A se realizan en una secuencia con nombre *sweep*. Cada *sweep* consiste de como máximo $\frac{n(n+1)}{2}$ rotaciones (pues depende de cuántas columnas son o no ortogonales) y en cada *sweep* se ortogonalizan 2 columnas. El número de *sweeps* a realizar se controla con esta variable.

**Salidas:**

* $V \in \mathbb{R}^{n \times n}$ (dada por $V\leftarrow V^k$, al final del proceso iterativo),
* $W \in \mathbb{R}^{m \times n}$ (dada por $W\leftarrow A^k$, al final del proceso iterativo),

**Procedure (pseudo-codigo)**

`````````````````python
# Inicializar parametros
k <- 0
A^(k) <- A
V^(k) <- I_n # Identidad de R^n
num_columnas_ortogonales <- 0 # variable auxiliar para contar cantidad de columnas or

# Proceso iterativo

Mientras k <= maxsweeps o bien num_columnas_ortogonales<= n(n+1)/2:
  for (i,j) \in \{(1,2), (1,3), \ldots (1,n)\} \cup  \{(2,3), (2,4), \ldots (2,n)\} \cup  \{(3,4),\ldots (2,n)\} \cup \{(n-1,n)\}:
    num_columnas_ortogonales <- 0 # reinicia variable
    
    #Revisar si columnas i y j de A^k son ortogonales segun tolerancia
    Si a_i^k y a_i^k son ortogonales:
      num_columnas_ortogonales aumenta en 1
    En otro caso:
      # Calcula coeficientes de submatriz (i,j) de A^k
      a <- || a_i^k||_2^2 # norma 2 al cuadrado de la columna i de A^k
      b <- || a_j^k||_2^2 # norma 2 al cuadrado de la columna j de A^k
      c <-  transpose(a_i^k) * a_j^k # Producto punto de columnas i y j de A^k

      # Calcula coeficientes de rotaciones de Givens para diagonalizar submatriz previa
      xi <- (b-a)/(2c)
      t <- signo(xi)/[|xi|+\sqrt(1+xi^2)]
      cs <- 1/[sqrt(1+t^2)]
      sn <- cs * t
      
      # Actualiza columnas i,j de A^k tras aplicar rotacion de Givens
      Para l in {1,...,n}:
        temp <- A_{li}^k # valor l-esimo de columna i de A^k
        A_{li}^k <- cs*temp - sn*A_{lj}^k #A_{lj}^k es valor l-esimo de columna j de A^k
        A_{lj}^k <- sn*temp + cs*A_{lj}^k #A_{lj}^k es valor l-esimo de columna j de A^k

      # Actualiza matriz V^k tras aplicar rotacion de Givens
      Para l in {1,...,n}:
        temp <- V_{li}^k # valor l-esimo de columna i de V^k
        V_{li}^k <- cs*temp - sn*V_{lj}^k #V_{lj}^k es valor l-esimo de columna j de A^k
        V_{lj}^k <- sn*temp + cs*A_{lj}^k #V_{lj}^k es valor l-esimo de columna j de A^k
      
      k <- k+1 # Actualiza contador de sweeps
      {\sigma_1, \sigma_2, ..., \sigma_r} <- normas de columnas de A^k # Estimar valores singulares con normas de columnas de A^k
      U <- A^k con columnas normalizadas # Estimar U con normalizando en norma 2 a las columnas de A^k
`````````````````

**Notas:**

* Por favor, incluir la documentación correspondiente (ver chuck en el numeral 1.1),
* Para la implementación, echar mano de las funciones implementadas en las tareas de los numeras 1.1, 1,2 y 1.3

### 2.2 Linear solver aproximating SVD decomposition using One-sided Jacobi algorithm

**Propósito:**

Dada una matriz $ A \in \mathbb{R}^{m \times n}$ así como un vector $b \in \mathbb{R}^m$ y su descomposición SVD, crear una función *SVD_solver* que permita aproximar la soluci´on el sistema lineal $Ax=b$.

* Estimar  la descomposición SVD de $A$ está dada por matrices $U \in \mathbb{R}^{m \times m}$, $\Sigma\in \mathbb{R}^{m \times n}$ y $V \in \mathbb{R}^{n \times n}$,  obteniendo versiones numéricas de $V$, $U$ y $\Sigma$, con procedimiento descrito en 2.1, mediante la función *SVD_jacobi*.
* A partir de este punto y la implementación realizada en 1.4:
  * Resolver el sistema de ecuaciones $$ \Sigma  d = U^T b$$, empleando alguna de las librerías [backsolve/forwardsolve](https://stat.ethz.ch/R-manual/R-devel/library/base/html/backsolve.html) de R.
  * Formar al vector $x:=Vd$ la cual será solución del sistema $Ax=b$.

**Notas:** 

* La matriz $A$ y el vector $b$ deben ser parte de las entradas la función.

Por favor, incluir la documentación correspondiente (ver chuck en el numeral 1.1).

```R
SVD_solver <- function(...){
    # Breve explicacion del objetivo de la funcion.
    #Args: # Argumentos involucrados en la definicion de la funcion
    #    ...
    #Returns:
    #    ...
  
  # Corpus de la funcion
  
}
```
