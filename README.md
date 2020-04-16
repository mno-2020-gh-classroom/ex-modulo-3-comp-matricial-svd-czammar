Método de eliminación por bloques  y descomposición SVD para sistemas lineales
==============================

El propósito de este proyecto es consolidar una **implementación del método de eliminación por bloques para la solución de sistemas lineales, empleando el lenguaje de programación R**. Para ello, también se aborda la **implementación del algoritmo One-Sided Jacobi** para aproximar numéricamente la descomposición en valores singulares de una matriz (SVD, por sus siglas en ingles) el cual se emplea como una herramienta que permite resolver numéricamente sistemas lineales más pequeños que surjen al considerar bloques de una matriz.

## Tabla de contenido

1. [Introduction](https://github.com/dssg/usal_echo#introduction)
2. [Overview](https://github.com/dssg/usal_echo#overview)

## Introducción 

Al tratar de encontrar la solución de sistemas lineales de tamaño considerable, una estrategia interesante es recurrir a métodos que intercambien el problema original por hallar la respuesta de una serie de sistemas de ecuaciones de menor dimensión, en los que  se cuente con otras técnicas que nos permitan encontrar su solución para posteriormente armar la solución del sistema que originalmente era de ínterés. Tal idea es en la que se apoya el método de eliminación por bloques para la solución de un sistema lineal, echando mano de que la inversa de una matriz se puede escribir en términos de sus bloques a través de la fórmula del complemento de Schur.

Sin embargo, ello implica tener que invertir bloques más pequeños de la matriz, o equivalentemente resolver sistemas lineales asociados a estos.

Para ello, es útil la descomposición SVD de una matriz, pues una vez conocida se puede explorar con el propósito de hallar la solución de cualquier sistema lineal asociado a esta, pues al encontrarse en términos de matrices ortogonales y diagonales, se puede aprovechar su estructura para replantear un sistema en uno que puede resolver con sustitución hacia adelante o atrás. Ello se puede lograr, por ejemplo, con la del **algoritmo One-Sided Jacobi** para aproximar numéricamente la descomposición SVD.

## Overview

Para dar solución a un sistema lineal $Ax=b$ con el método de eliminación por bloques se debe:

* 1) dividir una matriz en bloques, 
* 2) re-expresar el 

Por tales motivos, este proyecto también aborda la implementación del **algoritmo One-Sided Jacobi** para aproximar numéricamente la descomposición en valores singulares de una matriz (SVD, por sus siglas en ingles), aprovechando sus propiedades para resolver un sistema lineal específico, junto con su aplicación inmediata para llevar a cabo el método de eliminación por bloques.

En este sentido, el proyecto gira en torno a los siguientes ejes:

1. **Algoritmo One-sided Jacobi y descomposición SVD** [pendiente: desarrollo]
2. **Complemento de Schur** [pendiente: desarrollo].
3. **Método de eliminación por bloques para resolver sistemas lineales** [pendiente: desarrollo].





## Flujo de trabajo en Github

grupo de programación, grupo de revisión y una persona project manager



Este proyecto ha sido desarrollado de forma colaborativa, en donde cada unos de los equipos se ha encargo de implementar módulos de código que permiten [pendiente: desarrollo]

Para re



![](/Users/cesar/github/ex-modulo-3-comp-matricial-svd-czammar/images/gitflow.png)

**Fuente:** Notas del curso *Programación para Ciencia de Datos* de la Maestría en Ciencia de Datos del ITAM (2019). Véase https://github.com/ITAM-DS/programming-for-data-science-2019/blob/master/handbook.pdf

## Requerimientos de infraestructura
A efecto de que el equipo de programación y revisión tuviera un entorno común de trabajo para el desarrollo del proyecto, se empleó la imagen de docker basada en R del curso MNO 2020 (palmoreck/jupyterlab_r_kernel:1.1.0)

```bash
docker run --rm -v ($pwd):/datos --name jupyterlab_r_kernel_local -p 8888:8888 -d palmoreck/jupyterlab_r_kernel:1.1.0
```

Con ello se habilitó la posibilidad de realizar el trabajo mediante *Jupyter Notebooks*.



Para probar el código implementado sobre una matriz de dimensiones $10^4 \times 10^4$, usamos una instancia de AWS EC2, en donde se a su vez se guardaron los resultados obtenidos.

```bash
Infrastructure: AWS

+ AMI: RStudio-1.2.1335_R-3.6.0_CUDA-10.0_cuDNN-7.5.1_ubuntu-18.04-LTS-64bit - ami-0226a8af83fcecb43
+ EC2 instance: t2.2xlarge
    + vCPU: 8
    + RAM: 32 GB
+ OS: ubuntu 18.04 LTS
+ Volumes: 1
    + Type: gp2
    + Size: 20 GB
```

*Nota:* Instancia considerada bajo las indicaciones de la referencia https://towardsdatascience.com/how-to-run-rstudio-on-aws-in-under-3-minutes-for-free-65f8d0b6ccda

## Project Organization

```bash
├── README.md            <- Some text
├── References           <- Carpeta de materiales usados para desarrollo del proyecto
│   ├── 3.3.d.SVD.ipynb  <- Makefile with commands like `make data` or `make train`
│   ├── Images
│   ├── Minutas
│   ├── Readme.md
│   └── Simplified_SVD_OneSidedJacobi_Algorithm.md
├── jupyter
│   └── Ex_CM_SVD.ipynb
├── results               <- Carpeta de resultados
├── test                  <- Reportes derivados de issues del Equipo de Revisión
│   ├── Rev_FuncionSigno.ipynb
│   ├── Rev_GeneracionIndices.ipynb
│   ├── Rev_Generacion_de_Indices.ipynb
│   ├── Rev_One-sided_Jacobi.ipynb
│   ├── Rev_Solver_SVD.ipynb
│   ├── Rev_VerifOrtogonalidad.ipynb
│   └── Rev_parte2_One-sided_Jacobi.ipynb
```

