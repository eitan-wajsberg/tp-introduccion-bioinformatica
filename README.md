# TP Introducción a la Bioinformática — Parte 1

La enfermedad hereditaria elegida es la **Enfermedad de Huntington (HD)**, catalogada en OMIM bajo el código [#143100](https://www.omim.org/entry/143100?search=Huntington&highlight=huntington). El gen asociado es **HTT** (Huntingtin), ubicado en el cromosoma 4p16.3. La secuencia de referencia utilizada es el transcripto [NM_001388492.1](https://www.ncbi.nlm.nih.gov/nuccore/NM_001388492.1) (huntingtin isoforma 1, Homo sapiens), obtenida de la base de datos [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/3064) en formato GenBank.

## La Enfermedad de Huntington y el gen HTT

### ¿Qué es la Enfermedad de Huntington?

La Enfermedad de Huntington es una enfermedad neurodegenerativa progresiva e incurable que afecta al sistema nervioso central. Sus síntomas principales son la pérdida del control motor (movimientos involuntarios llamados corea), deterioro cognitivo y alteraciones psiquiátricas. Los síntomas suelen aparecer entre los 30 y 50 años, aunque existe una forma juvenil que puede manifestarse antes de los 20. Una vez que los síntomas aparecen, la enfermedad avanza inevitablemente hasta requerir cuidado total, con un promedio de supervivencia de 15 a 20 años desde el diagnóstico.

### Información por niveles biológicos

**Nivel genómico (DNA)**
El gen HTT se encuentra en el cromosoma 4, en la posición 4p16.3, y tiene una longitud de aproximadamente 170.000 pares de bases. Contiene 67 exones. La mutación causante de la enfermedad se encuentra en el exón 1 y consiste en una expansión del triplete CAG (citosina-adenina-guanina) más allá de un umbral crítico. En personas sanas, este triplete se repite entre 10 y 35 veces. A partir de 40 repeticiones, la enfermedad es prácticamente inevitable.

**Nivel transcripcional (mRNA)**
El gen HTT se transcribe en un mRNA que, luego del procesamiento (splicing), produce el transcripto maduro que usamos en este TP: NM_001388492.1, de 13.472 pares de bases. La región CAG expandida en el mRNA produce un tramo de codones CAG en tándem que codifica una cadena de glutaminas (Q) anormalmente larga en la proteína resultante.

**Nivel proteico**
La proteína huntingtina (HTT) tiene 3.144 aminoácidos y un peso molecular de aproximadamente 348 kDa. Es una proteína muy grande que se expresa en prácticamente todos los tejidos del cuerpo, pero en mayor concentración en el cerebro, especialmente en las neuronas. Su función exacta no está del todo clara, pero se sabe que participa en el transporte intracelular, en la supervivencia neuronal y en la regulación de la expresión de otros genes. Interactúa con más de 100 proteínas distintas. La versión mutante (mHTT) contiene una cadena de glutaminas expandida en su extremo N-terminal que tiende a formar agregados tóxicos dentro de las células, especialmente en las neuronas espinosas medianas del estriado, una región cerebral clave para el control del movimiento.

**Nivel celular y tisular**
La acumulación de la proteína mutante causa disfunción y muerte neuronal progresiva, principalmente en el estriado y la corteza cerebral. Las neuronas espinosas medianas del estriado son las más vulnerables. A medida que estas neuronas degeneran, el cerebro literalmente se encoge, lo que puede verse claramente en imágenes de resonancia magnética de pacientes con HD avanzada.

### Mutaciones conocidas

La única mutación causante de HD es la **expansión del triplete CAG en el exón 1 del gen HTT**. El número de repeticiones determina el estado clínico:

| Repeticiones CAG | Estado |
|-----------------|--------|
| ≤ 26 | Normal |
| 27 – 35 | Normal pero inestable (puede expandirse en la siguiente generación) |
| 36 – 39 | Penetrancia reducida (puede o no desarrollar HD) |
| ≥ 40 | Penetrancia completa (desarrollará HD inevitablemente) |
| ≥ 60 | Forma juvenil (síntomas antes de los 20 años) |

### Patrón de herencia

La Enfermedad de Huntington se hereda de forma **autosómica dominante**. Esto significa que el gen HTT está en uno de los 22 pares de cromosomas no sexuales (autosomas), y que basta con heredar **una sola copia mutada** del gen para desarrollar la enfermedad. En otras palabras, si uno de tus padres tiene HD, tenés un 50% de probabilidad de heredar la mutación.

A diferencia de enfermedades ligadas al cromosoma X (como el daltonismo, donde los hombres son mucho más afectados porque solo tienen un cromosoma X), la Enfermedad de Huntington afecta por igual a hombres y mujeres, ya que el gen HTT está en el cromosoma 4 y todos tenemos dos copias de ese cromosoma.

### Curiosidad: anticipación genética

Un fenómeno particular de la HD es la **anticipación genética**: a medida que el gen mutado se transmite de generación en generación, el número de repeticiones CAG tiende a aumentar, lo que hace que la enfermedad aparezca más temprano y de forma más severa en las generaciones siguientes. Este fenómeno ocurre especialmente cuando la mutación se hereda del padre, ya que los repetidos CAG son más inestables durante la espermatogénesis (producción de espermatozoides) que durante la ovogénesis. Es decir, un padre con 42 repeticiones CAG puede transmitirle a su hijo una versión con 50 o más repeticiones, acortando el tiempo hasta la aparición de los síntomas.

---

## Estructura del repositorio

```
tp-introduccion-bioinformatica/
├── README.md
├── .gitignore
├── data/
│   ├── sequence.gb          # Archivo GenBank del mRNA de referencia del gen HTT (NM_001388492.1)
│   ├── ORFs_HTT.fas         # Los 6 marcos de lectura traducidos a aminoácidos
│   └── HTT_correcto.fas     # Solo el marco de lectura correcto (+2)
├── scripts/
│   ├── verificar_entorno.sh # Script para verificar que el entorno está correctamente instalado
│   ├── Ex1.pl               # Ejercicio 1: procesamiento de secuencias y traducción
│   ├── Ex2_local.pl         # Ejercicio 2: BLAST local contra SwissProt
│   └── Ex2_remoto.pl        # Ejercicio 2: BLAST remoto contra servidor NCBI
└── results/
    ├── blast_local.out      # Resultado del BLAST local
    └── blast_remoto.out     # Resultado del BLAST remoto
```

> **Nota:** La carpeta `blast_db/` con la base de datos SwissProt no se incluye en el repositorio por su tamaño.

---

## Modo de uso

> **Importante:** todos los scripts deben ejecutarse siempre desde la raíz del repositorio, no desde dentro de la carpeta `scripts/`. Esto garantiza que las rutas de input y output funcionen correctamente.

```bash
# Correcto
cd tp-introduccion-bioinformatica
perl scripts/Ex1.pl data/sequence.gb

# Incorrecto
cd tp-introduccion-bioinformatica/scripts
perl Ex1.pl ../data/sequence.gb
```

---

## Ejercicio 1 — Procesamiento de secuencias

### ¿Qué hace el script?

El script `Ex1.pl` lee un archivo GenBank de un mRNA de referencia, genera los 6 marcos de lectura posibles (3 en la hebra directa y 3 en la hebra complementaria inversa), traduce cada uno a su secuencia de aminoácidos, e identifica cuál es el marco de lectura correcto usando la anotación CDS que viene incluida en el archivo GenBank. Los resultados se escriben en un archivo FASTA.

### ¿Qué es un marco de lectura?

El mRNA es una cadena de nucleótidos que se lee de a 3 (codones) para producir aminoácidos. Dependiendo desde qué posición se empiece a leer, se obtienen proteínas completamente distintas. Hay 3 posiciones posibles en la hebra directa (+1, +2, +3) y 3 en la hebra complementaria inversa (-1, -2, -3), dando un total de 6 marcos de lectura posibles. Solo uno de ellos es el correcto y produce la proteína real.

### Ejecución

```bash
perl scripts/Ex1.pl data/sequence.gb
```

### Output generado

```
data/ORFs_HTT.fas
```

### Resultado

El marco de lectura correcto identificado es el **+2**, lo que significa que la secuencia codificante del gen HTT comienza en la segunda posición del mRNA. Esto coincide con la anotación CDS del archivo GenBank de referencia.

Para trabajar en los ejercicios siguientes se extrajo solo la secuencia del marco correcto:

```bash
awk '/^>.*frame_\+2_MARCO_CORRECTO/{p=1} /^>/ && !/frame_\+2_MARCO_CORRECTO/{p=0} p' data/ORFs_HTT.fas > data/HTT_correcto.fas
```

---

## Ejercicio 2.a — BLAST

### ¿Qué hacen los scripts?

Los scripts `Ex2_local.pl` y `Ex2_remoto.pl` toman como input la secuencia de aminoácidos del marco de lectura correcto y realizan una búsqueda BLAST contra la base de datos SwissProt para encontrar secuencias similares en otros organismos.

Se implementaron dos variantes:
- **BLAST local:** corre directamente en la máquina usando la base de datos SwissProt descargada localmente. Es más rápido y no depende de la conexión.
- **BLAST remoto:** envía la secuencia al servidor del NCBI y espera la respuesta. No requiere tener la base de datos instalada pero puede tardar varios minutos dependiendo del tamaño de la secuencia y el tráfico del servidor.

### Ejecución

BLAST local:
```bash
perl scripts/Ex2_local.pl data/HTT_correcto.fas blast_db/swissprot
```

BLAST remoto:
```bash
perl scripts/Ex2_remoto.pl data/HTT_correcto.fas
```

### Outputs generados

```
results/blast_local.out
results/blast_remoto.out
```

### Resultados obtenidos

| Hit | Organismo | Identidad | E-value | Score |
|-----|-----------|-----------|---------|-------|
| P42858.2 | Homo sapiens | 100.0% | 0.0 | 16594 |
| P42859.2 | Mus musculus (ratón) | 91.2% | 0.0 | 14571 |
| P51111.1 | Rattus norvegicus (rata) | 90.8% | 0.0 | 14494 |
| P51112.1 | Takifugu rubripes (pez globo) | 71.5% | 0.0 | 11468 |
| Q76P24.1 | Dictyostelium discoideum | 30.9% | 3e-18 | 241 |

Los resultados obtenidos por ambas variantes son consistentes: los mismos 5 hits en el mismo orden de relevancia. Las pequeñas diferencias en los scores e identidades (por ejemplo, 91.2% vs 91.1% para el ratón, o E-value de 3e-18 vs 5e-18 para Dictyostelium) son normales y esperables. El servidor remoto del NCBI puede utilizar una versión ligeramente distinta de la base de datos SwissProt o parámetros internos levemente diferentes a los de la instalación local. Lo importante es que las conclusiones biológicas son las mismas en ambos casos.

---

## Ejercicio 2.b — Interpretación del resultado del BLAST

### Las secuencias encontradas

El primer hit es la huntingtina humana (P42858.2) con 100% de identidad, lo cual confirma que la secuencia que estamos analizando es correcta y corresponde exactamente a la proteína HTT humana almacenada en SwissProt.

Los hits 2 y 3 son la huntingtina de ratón y rata, con identidades del 91.2% y 90.8% respectivamente. Esto refleja la alta conservación evolutiva de esta proteína entre mamíferos: el gen HTT es esencial para el desarrollo neurológico y ha sido muy conservado a lo largo de la evolución.

El hit 4, el pez globo (Takifugu rubripes), muestra una identidad del 71.5%. A pesar de ser un vertebrado mucho más distante evolutivamente que los mamíferos, la proteína sigue siendo reconociblemente similar, lo que indica que HTT cumple funciones fundamentales conservadas en todos los vertebrados.

El hit 5, Dictyostelium discoideum (un moho mucilaginoso unicelular), es el más interesante y sorprendente: con solo 30.9% de identidad pero un E-value de 3e-18, la similitud sigue siendo estadísticamente significativa. Esto sugiere que algunas regiones funcionales de la huntingtina ya existían en organismos muy primitivos, antes de la aparición de los animales multicelulares.

### Significado de los valores estadísticos

**Score (puntuación):** Es un número que refleja qué tan bien se alinean dos secuencias. Se calcula sumando puntos por cada posición donde los aminoácidos coinciden o son similares, y restando puntos por los gaps (espacios que se introducen para alinear mejor). Cuanto más alto el score, mejor el alineamiento.

**E-value (valor esperado):** Es el valor estadístico más importante del BLAST. Representa cuántos alineamientos con ese score o mejor se esperaría encontrar por pura casualidad en una base de datos del tamaño de SwissProt. Un E-value de 0.0 significa que la probabilidad de que ese alineamiento sea producto del azar es prácticamente nula: la similitud es real y tiene significado biológico. Un E-value de 3e-18 (equivalente a 0.000000000000000003) sigue siendo extremadamente significativo aunque la similitud sea menor. Como regla general, se considera significativo cualquier E-value menor a 0.001, y altamente significativo cualquier valor menor a 1e-10.

**Identidad:** Es el porcentaje de posiciones en el alineamiento donde los dos aminoácidos son exactamente iguales. Un 100% indica que las secuencias son idénticas. Un 30.9% puede parecer bajo, pero en proteínas largas con E-values muy pequeños, ese nivel de identidad es suficiente para concluir que las proteínas comparten un ancestro común y probablemente funciones similares.

La combinación de un score alto, un E-value cercano a cero y una identidad elevada es la señal más confiable de que dos proteínas están genuinamente relacionadas evolutivamente.

