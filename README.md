Método de eliminación por bloques  y descomposición SVD para sistemas lineales
==============================

El propósito de este proyecto es consolidar una **implementación del método de eliminación por bloques para la solución de sistemas lineales, empleando el lenguaje de programación R**. Para ello, también se aborda la **implementación del algoritmo One-Sided Jacobi** para aproximar numéricamente la descomposición en valores singulares de una matriz (SVD, por sus siglas en ingles) el cual se emplea como una herramienta que permite resolver numéricamente sistemas lineales más pequeños que surjen al considerar bloques de una matriz.

## Tabla de contenido

1. [Introduction](https://github.com/mno-2020-gh-classroom/ex-modulo-3-comp-matricial-svd-czammar/blob/master/README.md#introducción)
2. [Overview](https://github.com/mno-2020-gh-classroom/ex-modulo-3-comp-matricial-svd-czammar/blob/master/README.md#overview)

## Introducción 

Al tratar de encontrar la solución de sistemas lineales de tamaño considerable, una estrategia interesante es recurrir a métodos que intercambien el problema original por hallar la respuesta de una serie de sistemas de ecuaciones de menor dimensión, que puedan ser resueltos con otras técnicas y cuyas soluciones permitan armar la solución del sistema que originalmente era de ínterés. 

Por otro lado, la descomposición SVD de una matriz $A$ consiste hallar una factorización de ésta en términos de dos matrices ortonormales U y V asociadas a su dominio y co-dominio, junto con una matriz diagonal $\Sigma$ cuyas entradas se encuentra determinadas por los **valores singulares** de $A$. En el contexto de un problema de resolver un sistema $Ax=b$, dicha estructura se puede aprovechar para usar la estructura de las matrices $U$, $V$ y $\Sigma$, para replantear el sistema original, en otro más sencillo que se puede resolver con sustitución atrás/adelante, para posteriorme recuperar la solución de interés con una multiplicación *matriz-vector*.

Así pues, para implementar la solución de sistemas lineales la implementación a través del método de eliminación por bloques, se puede  método que echar mano de la factorización SVD para dar solución a los sistemas lineales de menor dimensión que resultar de re-expresar el problema en términos de sus correspondientes bloques.

Estas ideas fueron exploraras en el presente proyecto, para la implementación llevada a cabo del método de eliminación por bloques.

## Overview

**Método de eliminación por bloques**

El método de eliminación por bloques consiste en obtener solución a un sistema lineal $Ax=b$ pensando en que $A$, $b$ y $x$ se puede dividir en términos bloques, $A_{11}$,$A_{12}$,$A_{21}$ y $A_{22}$, $b_1$ y $b_2$ junto con $x_1$ y $x_2$, tal como se aprecia en la siguiente imagen:

![bloques](./images/bloques.png)

En el caso de que $A_{11}$ sea invertible (no singular), las ecuaciones inducidas por los bloques se pueden resolver como sigue:

1. Si conocemos $x_2$, entonces de la primera ecuación de bloques![A_{11}x_1 + A_{12}x_2 = b_1](https://render.githubusercontent.com/render/math?math=A_%7B11%7Dx_1%20%2B%20A_%7B12%7Dx_2%20%3D%20b_1)podemos obtener $x_1$ al resolver el sistema $A_{11} x_1 = b_1-A_{12}x_2$ o de modo equivalente $x_1 = A_{11}^{-1}(b_1-A_{12}x_2)$, que es un sistema lineal de menores dimensiones.
2. Por otro lado, haciendo los cálculos, de la segunda ecuación por bloques se sigue que $S x_2 = b_2 - A_{21}A_{11}^{-1}b_1$, en donde a) $S := A_{22}-A_{21}A_{11}^{-1}A_{12}$ se conoce como el **complemento de Schur** del bloque $A_{11}$ en $A$, y b) $S$ es no singular si y sólo si $A$ es no singular. Así para estimar $x_2$, basta resolver $S x_2 = b_2 - A_{21}A_{11}^{-1}b_1$ que es un sistema más pequeño.
3. La solución del sistema $Ax=b$ se obtiene en términos de $x_1$ y $x_2$.

Con este procedimiento reducimos el problema de hallar la solución de un sistema $Ax = b$, pero debemos contar con un método que permita resolve los sistemas lineales más pequeños que se describen en los numerales 1. y 2.

**Algoritmo de One-Sided Jacobi para descomposición SVD**

Por resultados de álgebra lineal y análisis matemático, se sabe que para una matriz ![$A \in \mathbb{R}^{mxn}$](https://render.githubusercontent.com/render/math?math=A%20%5Cin%20%5Cmathbb%7BR%7D%5E%7Bmxn%7D&mode=inline) es posible encontrar una factorización de ésta en términos de matrices ![$U \in \mathbb{R}^{mxm}, V \in \mathbb{R}^{nxn}$](https://render.githubusercontent.com/render/math?math=U%20%5Cin%20%5Cmathbb%7BR%7D%5E%7Bmxm%7D%2C%20V%20%5Cin%20%5Cmathbb%7BR%7D%5E%7Bnxn%7D&mode=inline) ortogonales tales que: ![$A = U\Sigma V^T$](https://render.githubusercontent.com/render/math?math=A%20%3D%20U%5CSigma%20V%5ET&mode=inline) con ![$\Sigma = diag(\sigma_1, \sigma_2, \dots, \sigma_p) \in \mathbb{R}^{mxn}$](https://render.githubusercontent.com/render/math?math=%5CSigma%20%3D%20diag%28%5Csigma_1%2C%20%5Csigma_2%2C%20%5Cdots%2C%20%5Csigma_p%29%20%5Cin%20%5Cmathbb%7BR%7D%5E%7Bmxn%7D&mode=inline), ![$p = \min\{m,n\}$](https://render.githubusercontent.com/render/math?math=p%20%3D%20%5Cmin%5C%7Bm%2Cn%5C%7D&mode=inline) y ![$\sigma_1 \geq \sigma_2 \geq \dots \geq \sigma_p \geq 0$](https://render.githubusercontent.com/render/math?math=%5Csigma_1%20%5Cgeq%20%5Csigma_2%20%5Cgeq%20%5Cdots%20%5Cgeq%20%5Csigma_p%20%5Cgeq%200&mode=inline).

**Nota:** Existen representaciones alternas dicha factorización, pero similarmente tales se basan en encontrar matrices ortogonales asociadas al dominio y condominio de $A$, junto con sus valores singulares. Véase https://en.wikipedia.org/wiki/Singular_value_decomposition#Reduced_SVDs

Para obtener versiones numéricas de la descomposición SVD de una matriz existe un método iterativo conocido como el **algoritmo One-Sided Jacobi** que se basa en la aplicación sucesiva de rotaciones de Givens, para ortogonalizar las columnas de las matrices involucradas, así como estimar los respectivos valores singulares. En concreto, la idea en que se basa tal método es multiplicar a la matriz ![$A \in \mathbb{R}^{m \times n}$](https://render.githubusercontent.com/render/math?math=A%20%5Cin%20%5Cmathbb%7BR%7D%5E%7Bm%20%5Ctimes%20n%7D&mode=inline) por la derecha de forma repetida por matrices ortogonales de nombre **rotaciones Givens** hasta que se converja a ![$U \Sigma$](https://render.githubusercontent.com/render/math?math=U%20%5CSigma&mode=inline).

**Solución de un sistema lineal usando descomposición SVD**

Si conocemos la descomposición SVD de una matriz $A$, podemos resolver un sistema lineal $Ax=b$ como sigue

a) **Descomposición SVD:** Estimar factores ![$U, \Sigma, V$](https://render.githubusercontent.com/render/math?math=U%2C%20%5CSigma%2C%20V&mode=inline) tales que ![$A=U \Sigma V^T$](https://render.githubusercontent.com/render/math?math=A%3DU%20%5CSigma%20V%5ET&mode=inline).

b) **Solución de sistema intermedio:** resolver el sistema diagonal ![$\Sigma d = U^Tb$](https://render.githubusercontent.com/render/math?math=%5CSigma%20d%20%3D%20U%5ETb&mode=inline). Este paso se puede hacer con solución hacia delante o atrás, puesto que es un sistema diagonal.

c) **Armado de solución de sistema original:** hacer la multiplicación matriz vector para obtener ![$x$](https://render.githubusercontent.com/render/math?math=x&mode=inline): ![$x=Vd$](https://render.githubusercontent.com/render/math?math=x%3DVd&mode=inline).

**Roadmap**

En este sentido, el proyecto gira en torno al desarrollo de siguientes ejes:

1. **Funciones auxiliares:** se trata de una serie de funciones atómicas en R que serán base del resto de los algoritmos (por ejemplo, la función signo, función que determina si dos vectores son ortogonales, entre otras).
2. **Algoritmo One-sided Jacobi y descomposición SVD:** implementación de este método iterativo para obtener la descomposición en valores singulares de una matriz arbitraria.
3. **Solver de un sistema lineal usando descomposición SVD**
4. **Método de eliminación por bloques para resolver sistemas lineales:** usando el solver que emplea la descomposición SVD de los sistemas lineales inducidos por los bloques.

## Organización

Para el desarrollo del proyecto, los integrantes se agruparon en torno a dos equipos; uno encargado de la implementación de los diferentes métodos y algoritmos (**Equipo de Programación, o simplemente E-Prog**) y otro más encargado probar los métodos del E-Prog con diferentes parámetros y generar reportes de resultados con las variaciones de los parámetros (**Equipo de Revisión, o simplemente E-Rev**). Ambos equipos fueron coordinados por un project manager (**PM**).

La división anterior se puede resumir mediante la siguiente tabla:

| #    | Rol                   | Persona      | Github    |
| ---- | --------------------- | ------------ | --------- |
| 1    | Grupo de programación | Danahi       | Danahirmt |
| 2    | Grupo de programación | Miguel       | Millan13  |
| 3    | Grupo de programación | Juan Pablo   | Pilo1961  |
| 4    | Grupo de revisión     | Dorely       | DorelyMS  |
| 5    | Grupo de revisión     | Javier       | valencig  |
| 6    | Grupo de revisión     | Leon         | lgarayva  |
| 7    | Project Manager       | César Zamora | czammar   |

Esencialmente, para realizar el trabajo correspondiente se empleó el [project board de Github](https://github.com/mno-2020-gh-classroom/ex-modulo-3-comp-matricial-svd-czammar/projects), en donde se crearon dos proyectos, juntos con sus milestones que sirvieron para agrupar los issues que cada integrante debía resolver. 

En este sentido, el [primer proyecto](https://github.com/mno-2020-gh-classroom/ex-modulo-3-comp-matricial-svd-czammar/projects/1) es relativo a los numerales 1 a 4 del **Roadmap** mencionado en el apartado de [Overview](https://github.com/mno-2020-gh-classroom/ex-modulo-3-comp-matricial-svd-czammar/blob/master/README.md#overview). En tanto el [segundo proyecto](https://github.com/mno-2020-gh-classroom/ex-modulo-3-comp-matricial-svd-czammar/projects/2) al desarrollo del reporte de resultados.

## Flujo de trabajo en Github

Para facilitar el sido desarrollado de forma colaborativa entre **E-Prog** , **E-Rev** y **PM**, se siguió un *Github flow*, consistente, en líneas generales, en la creación de ramas para resolver un issues específico, para solicitar la revisión del PM a través de un *Pull request*, y su posterior aprobación para unir los cambios hacia la rama *master*.

![gitflow](./images/gitflow.png)

**Fuente:** Notas del curso *Programación para Ciencia de Datos* de la Maestría en Ciencia de Datos del ITAM (2019). Véase https://github.com/ITAM-DS/programming-for-data-science-2019/blob/master/handbook.pdf

Cabe destacar que una ves solucionado el issue correspondiente, se borró la rama asociada para facilitar el entendimineto y administración del proyecto.

## Requerimientos de infraestructura
A efecto de que el **E-Prog** , **E-Rev** y **PM** tuvieran un entorno común de trabajo para el desarrollo del proyecto, se empleó la imagen de docker basada en R del curso MNO 2020 (palmoreck/jupyterlab_r_kernel:1.1.0)

```bash
docker run --rm -v ($pwd):/datos --name jupyterlab_r_kernel_local\
-p 8888:8888 -d palmoreck/jupyterlab_r_kernel:1.1.0
```

Con ello se habilitó la posibilidad de realizar el trabajo mediante sucesivos *Jupyter Notebooks*.

## Organización del proyecto

La organización del proyecto se realizó a través una serie de carpetas, entre las cuales destacan:

+ hola

En complemento, se presenta una version esquemática de la organización de repositorio del proyecto:

```bash
├── README.md            <- Archivo readme del proyecto
├── References           <- Carpeta de materiales usados para desarrollo del proyecto
│   ├── 3.3.d.SVD.ipynb  <- Makefile with commands like `make data` or `make train`
│   ├── Images
│   ├── Minutas
│   ├── Readme.md
│   └── Simplified_SVD_OneSidedJacobi_Algorithm.md
├── jupyter
│   └── Ex_CM_SVD.ipynb
├── results               <- Contiene el reporte ejecutivo de resultados
├── test                  <- Reportes derivados de issues del Equipo de Revisión
│   ├── Rev_FuncionSigno.ipynb
│   ├── Rev_GeneracionIndices.ipynb
│   ├── Rev_Generacion_de_Indices.ipynb
│   ├── Rev_One-sided_Jacobi.ipynb
│   ├── Rev_Solver_SVD.ipynb
│   ├── Rev_VerifOrtogonalidad.ipynb
│   └── Rev_parte2_One-sided_Jacobi.ipynb
```

