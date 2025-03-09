# <center>**TAREA #2: Análisis de sensibilidad lineal**</center>


# <center>*Introducción al Análisis de Sensibilidad Lineal en Sistemas de Potencia*</center>
<div style="text-align: justify;">
El análisis de sensibilidad lineal es una herramienta clave en la evaluación de la seguridad y confiabilidad de los sistemas eléctricos de potencia. Su objetivo es determinar cómo cambios en la generación o en la topología del sistema afectan los flujos de potencia en las líneas de transmisión. Estos cálculos se realizan bajo el modelo de flujo de potencia DC, lo que permite obtener resultados de manera rápida y eficiente.
</div><br>
Entre los factores más utilizados en este análisis se encuentran:


### Generation Shift Factor (GSF): 
<div style="text-align: justify;">
El GSF indica la sensibilidad del flujo de una línea ante un cambio en la generación de un nodo específico, asumiendo que la variación se compensa con el nodo de referencia.
</div><br>


### Line Outage Distribution Factor (LODF):
<div style="text-align: justify;">
El LODF mide el impacto de la desconexión de una línea sobre el flujo de potencia en otras líneas del sistema, ayudando a predecir sobrecargas en eventos de contingencia.
</div><br>

### Node-Based Line Outage Distribution Factor (NBLODF):
<div style="text-align: justify;">
Este factor extiende el LODF, considerando la distribución de carga en los nodos afectados por la salida de una línea.
</div><br>

### Matriz Bbus:
<div style="text-align: justify;">
Esencial en el cálculo de estos factores, representa la relación entre los voltajes de los nodos y las inyecciones de potencia, excluyendo el nodo slack.
Estos factores permiten evaluar el impacto de posibles fallas y tomar medidas preventivas para mejorar la estabilidad operativa del sistema eléctrico.
</div><br>

# <center>**Marco teórico**</center>


### 11.3.2 Linear Sensitivity Factors
<div style="text-align: justify;">
Los factores de sensibilidad lineal se utilizan en el análisis de seguridad de sistemas de potencia para evaluar el impacto de contingencias, como la pérdida de una línea o cambios en la generación. Estos factores permiten una evaluación rápida del flujo de potencia utilizando aproximaciones lineales basadas en el flujo de carga DC.
</div><br>

### Matriz Bbus

<div style="text-align: justify;">
La **matriz Bbus** representa la relación entre los flujos de potencia activa y los ángulos de voltaje en un sistema de potencia. Se deriva del modelo de flujo de carga DC y se usa para calcular los factores de sensibilidad. Su expresión general es:
</div><br>

$$
W = B_{bus}^{-1}
$$

donde:
- $ W $ es la matriz inversa de reactancias de la red eléctrica.  
- $ B_{bus} $ es la matriz de susceptancias de la red.  

Esta matriz es fundamental en el cálculo de los factores de sensibilidad lineal.

### Node-Based Line Outage Distribution Factor (NBLODF)

El **NBLODF** mide el impacto en el flujo de una línea específica cuando otra línea en el sistema se desconecta. Se calcula como:

$$
NBLODF_{i,k} = \frac{\Delta P_i}{P_k}
$$

donde:
- $ \Delta P_i $ es el cambio en el flujo de potencia de la línea monitoreada $i$.
- $ P_k $ es el flujo previo a la contingencia en la línea desconectada $ k $.

Este factor ayuda a predecir posibles sobrecargas en líneas después de una falla en la red.

### Line Outage Distribution Factor (LODF)

El **LODF** determina cómo cambia el flujo de potencia en una línea cuando otra línea es desconectada. Se expresa como:

$$
LODF_{i,k} = \frac{\Delta P_i}{P_k^{pre}}
$$

donde:
- $ \Delta P_i $ es el cambio en el flujo de la línea $ i $.
- $ P_k^{pre} $ es el flujo de la línea $ k $ antes de la contingencia.

Este factor es crucial en el análisis de contingencias porque permite estimar los efectos de fallas en la red de transmisión.

### Generation Shift Factor (GSF)

El **GSF** mide la sensibilidad del flujo de potencia en una línea debido a un cambio en la generación en un nodo específico. Se define como:

$$
GSF_{i,j} = \frac{\Delta P_i}{\Delta P_j}
$$

donde:
- $ \Delta P_i $ es el cambio en el flujo de la línea $i$.
- $ \Delta P_j $ es el cambio en la generación en el bus $j$.

Otra expresión utilizada para este factor es:

$$
\Delta f_e = \alpha_{l,i} \cdot \Delta P_i
$$

donde:
- $ \Delta f_e $ es el cambio en el flujo de potencia en la línea $l$,
- $ \alpha_{l,i} $ es el factor de desplazamiento de generación,
- $ \Delta P_i $ es el cambio en la generación en el bus $i$.

El **GSF** es útil para analizar cómo el desplazamiento de generación entre unidades afecta el flujo de potencia en la red y puede ayudar en la redistribución de generación para evitar sobrecargas.


# <center>**Funciones**</center>

La primera función que se realiza es para calcular la matriz Bbus.

*Requiere*

    """
    matriz_Bbus(lines, nodes)

    Calcula la matriz de susceptancia Bbus para una red eléctrica utilizando la información de las líneas y nodos.

    ### Entradas:
    - lines: DataFrame que contiene la información de las líneas de transmisión. Debe tener las siguientes columnas:
        - FROM: Nodo de envío de la línea.
        - TO: Nodo de recibo de la línea.
        - X: Reactancia de la línea (en ohmios).

    - Nodes: DataFrame que contiene la información de los nodos de la red.

    ### Salida:
    - Bbus: Matriz de susceptancia Bbus de la red, la matriz contiene los valores de susceptancia entre los nodos del sistema, que se utilizarán en el análisis de flujo de potencia y en la simulación de contingencias.

    ### Excepciones:
    - Si la reactancia (`X`) de una línea es cero o inválida, se lanza un error.
    - Si algún nodo de las líneas está fuera del rango de nodos disponibles, se lanza un error.

    ### Requiere:
    - `using LinearAlgebra`
    - `using DataFrames`
    - `using CSV`
    - `using Plots`
    """

La segunda función que se utiliza es la del Generation Shift Factor

*Requiere*

    """
    Generation_shift_α(lines, nodes, Bbus)

    Calcula el Generation Shift Factor (GSF) para una red eléctrica utilizando la matriz de susceptancia nodal.

    ### Entradas:
        - `lines`: DataFrame que contiene la información de las líneas de transmisión. Debe tener las siguientes columnas:
        - `FROM`: Nodo de envío de la línea.
        - `TO`: Nodo de recibo de la línea.
        - `X`: Reactancia de la línea (en ohmios).

    - `nodes`: DataFrame que contiene la información de los nodos de la red. Debe incluir:
        - `NUMBER`: Identificación del nodo.
        - `TYPE`: Tipo de nodo (ejemplo: 3 para nodo slack).

    - `Bbus`: Matriz de susceptancia nodal sin incluir el nodo slack.

    ### Salida:
    - `alpha_α`: Matriz con los factores de sensibilidad por cambio de generación, indicando cómo una variación en la generación en un nodo afecta el flujo en cada línea.

    ### Excepciones:
    - Si no hay exactamente un nodo slack (`TYPE == 3`), se lanza un error.
    - Si la matriz `Bbus` es singular, se utilizará la pseudo-inversa (`pinv(Bbus)`) en lugar de la inversa regular.
    - Si la reactancia (`X`) de una línea es cero o inválida, se lanza un error.

    ### Notas:
    - Se expande la matriz inversa de susceptancia para incluir el nodo slack en la matriz `W`.
    - Se calcula el **GSF** usando la diferencia de términos en la matriz expandida `W` dividida por la reactancia de la línea correspondiente.

    ### Requiere:
    - `using LinearAlgebra`
    - `using DataFrames`
    """

La tercera función que se utiliza es la del del Factor de Distribución de Corte de Línea por Nodo (NBLODF)

*Requiere*

    """
    Node_Based_Line_Outage_Distribution_Factor(lines, nodes, Bbus)

    Calcula el Factor de Distribución de Corte de Línea por Nodo (**NBLODF**), que mide la sensibilidad del flujo de potencia en los nodos cuando una línea es desconectada.

    ### Entradas:
        - `lines`: DataFrame que contiene la información de las líneas de transmisión. Debe tener las siguientes columnas:
        - `FROM`: Nodo de envío de la línea.
        - `TO`: Nodo de recibo de la línea.
        - `X`: Reactancia de la línea (en ohmios).

    - `nodes`: DataFrame que contiene la información de los nodos de la red. Debe incluir:
        - `NUMBER`: Identificación del nodo.
        - `TYPE`: Tipo de nodo (ejemplo: `3` para nodo slack).

    - `Bbus`: Matriz de susceptancia nodal sin incluir el nodo slack.

    ### Salida:
    - `NBLODF`: Matriz de tamaño `(total_nodes, total_lines)`, donde cada elemento representa el impacto en el flujo de un nodo específico ante la pérdida de una línea.

    ### Excepciones:
    - Si no hay exactamente un nodo slack (`TYPE == 3`), se lanza un error.
    - Si la matriz `Bbus` es singular, se utilizará la pseudo-inversa (`pinv(Bbus)`) en lugar de la inversa regular.
    - Si la reactancia de la línea (`X`) es inválida o cero, se lanza un error.

    ### Requiere:
    - `using LinearAlgebra`
    - `using DataFrames`
    """


La cuarta función que se utiliza es la del Factor de Distribución de Corte de Línea (LODF).

*Requiere*

    """
    Line_Outage_Distribution_Factor(lines, nodes, Bbus)

    Calcula el Factor de Distribución de Corte de Línea (**LODF**), que mide la sensibilidad del flujo en una línea cuando otra línea es desconectada.

    ### Entradas:
    - `lines`: DataFrame que contiene la información de las líneas de transmisión. Debe incluir las siguientes columnas:
        - `FROM`: Nodo de envío de la línea.
        - `TO`: Nodo de recibo de la línea.
        - `X`: Reactancia de la línea (en ohmios).

    - `nodes`: DataFrame que contiene la información de los nodos de la red. Debe incluir:
        - `NUMBER`: Identificación del nodo.
        - `TYPE`: Tipo de nodo (`3` para el nodo slack).

    - `Bbus`: Matriz de susceptancia nodal sin incluir el nodo slack.

    ### Salida:
    - `beta`: Matriz de tamaño `(total_lines × total_lines)`, donde cada elemento `beta[l, h]` indica la sensibilidad del flujo en la línea `l` ante la desconexión de la línea `h`.

    ### Excepciones:
    - Si no hay exactamente un nodo slack (`TYPE == 3`), se lanza un error.
    - Si la matriz `Bbus` es singular, se utilizará la pseudo-inversa (`pinv(Bbus)`) en lugar de la inversa regular.
    - Si la reactancia de una línea (`X`) es inválida o cero, se lanza un error.

    ### Requiere:
    - `using LinearAlgebra`
    - `using DataFrames`
    """


**Licencia**

Programa realizado por: Juan Manuel Agudelo Ocampo

Correo: m.agudelo@utp.edu.co

