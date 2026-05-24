# TP Introducción a la Bioinformática: Parte 1

Elegimos para realizar el trabajo la **Enfermedad de Huntington (HD)**, catalogada en OMIM con el código [#143100](https://www.omim.org/entry/143100?search=Huntington&highlight=huntington). El gen asociado es **HTT** (Huntingtin), ubicado en el cromosoma 4p16.3. La secuencia de referencia utilizada es el transcripto [NM_001388492.1](https://www.ncbi.nlm.nih.gov/nuccore/NM_001388492.1) (huntingtina isoforma 1, Homo sapiens), obtenida de la base de datos [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/3064) en formato GenBank.

&nbsp;

## La Enfermedad de Huntington y el gen HTT

### ¿Qué es la Enfermedad de Huntington?
La Enfermedad de Huntington es una enfermedad hereditaria que afecta al cerebro y por el momento no tiene cura. Las personas que la tienen van perdiendo el control de sus movimientos (hacen movimientos involuntarios llamados corea), se les va deteriorando la memoria y el pensamiento, y también tienen cambios de conducta y estado de ánimo. Los síntomas suelen aparecer entre los 30 y 50 años, y una vez que aparecen la enfermedad avanza sin parar. En promedio, las personas viven entre 15 y 20 años más después del diagnóstico.

### Información por niveles
Para entender esta enfermedad desde el punto de vista bioinformático, hay que mirarla en distintos niveles, desde el ADN hasta la función de la proteína que produce.

1. **SNP ([dbSNP - NCBI](https://www.ncbi.nlm.nih.gov/snp/?term=htt)):**
Los SNPs son variaciones de una sola letra en el ADN entre distintas personas. En el gen HTT hay muchos SNPs documentados, pero la causa de la enfermedad no es un SNP: es una expansión de un triplete (ver sección de mutación más abajo).
 
2. **Secuencia nucleotídica ([NCBI Nucleotide](https://www.ncbi.nlm.nih.gov/nuccore/NG_009378.1)):**
El gen HTT está en el cromosoma 4, en la posición 4p16.3. Ocupa aproximadamente 180.000 pares de bases y tiene 67 exones. Es un gen muy grande.
 
3. **Gen ([NCBI Gene - HTT](https://www.ncbi.nlm.nih.gov/gene/3064)):**
HTT es el símbolo oficial del gen. También se lo conoce como HD o IT15. A partir de este gen se produce la proteína huntingtina y es esencial para el desarrollo normal del organismo.
 
4. **ARNm ([NCBI Nucleotide - NM_001388492.1](https://www.ncbi.nlm.nih.gov/nuccore/NM_001388492.1)):**
El transcripto que usamos en este TP es el NM_001388492.1 (huntingtina isoforma 1, humana), que tiene 13.472 bases. Este es el ARNm sin intrones, que la célula usa para fabricar la proteína.
 
5. **Secuencia primaria de proteína ([NCBI Protein](https://www.ncbi.nlm.nih.gov/protein/NP_001375421.1)):**
La proteína huntingtina (HTT) tiene 3.142 aminoácidos. Es una de las proteínas más grandes del organismo humano.
 
6. **Estructura de la proteína ([RCSB PDB](https://www.rcsb.org/structure/8R2O)):**
La huntingtina es una proteína muy grande y su estructura tridimensional completa todavía no está del todo resuelta. Se conocen algunas partes, pero entender cómo se dobla en su totalidad sigue investigándose.
 
7. **Función de la proteína ([UniProt](https://www.uniprot.org/uniprot/P42858)):**
La huntingtina está presente en casi todos los tejidos del cuerpo, pero principalmente en el cerebro. Su función exacta todavía no se entiende del todo, pero se sabe que es esencial para que el organismo se desarrolle con normalidad y para que las neuronas funcionen bien.

### La mutación que causa la enfermedad
La Enfermedad de Huntington tiene una sola mutación conocida como causa: **una expansión del triplete CAG en el exón 1 del gen HTT**.

> ¿Qué significa esto? En el ADN hay una región donde la secuencia CAG se repite varias veces seguidas. En personas sanas, este triplete se repite entre 9 y 35 veces. Cuando el número de repeticiones supera cierto umbral, la persona desarrolla la enfermedad. Cuantas más repeticiones hay, más temprano aparecen los síntomas.
 
| Repeticiones CAG | Qué significa |
|---|---|
| 9 – 27 | Normal, no hay enfermedad |
| 27 – 35 | Normal, pero puede expandirse al pasarlo a los hijos |
| 36 – 39 | Zona gris: puede o no desarrollar HD dependiendo de otros factores |
| 40 o más | La persona va a desarrollar HD inevitablemente |
| 60 o más | Forma juvenil: los síntomas aparecen antes de los 20 años |
 
Cada triplete CAG codifica para el aminoácido glutamina (Q). Entonces, una expansión de CAG produce una cadena de glutaminas anormalmente larga en la proteína. Esa cadena larga hace que la proteína se doble mal y forme grumos tóxicos dentro de las neuronas, especialmente en una zona del cerebro llamada estriado, que es clave para controlar el movimiento.

### Cómo se hereda
La Enfermedad de Huntington se hereda de forma **autosómica dominante**. Esto significa:
- El gen HTT está en un autosoma (cromosoma 4), no en los cromosomas sexuales X o Y.
- Con una sola copia mutada alcanza para desarrollar la enfermedad. No hace falta que las dos copias estén afectadas.
- Si uno de tus padres tiene HD, tenés un 50% de chance de heredar la mutación.
- Afecta por igual a hombres y mujeres.

Esto la diferencia de enfermedades ligadas al cromosoma X, como el daltonismo, donde los hombres son mucho más afectados porque solo tienen un cromosoma X. En el caso de HD, como el gen está en el cromosoma 4 y todos tenemos dos copias de ese cromosoma, el riesgo es el mismo para todos.

### Curiosidad: anticipación genética
La enfermedad puede aparecer más temprano y de forma más grave en cada generación. Esto pasa porque el número de repeticiones CAG tiende a aumentar cuando el gen se transmite de padres a hijos, especialmente cuando lo pasa el padre. Por ejemplo, un padre con 42 repeticiones puede transmitirle a su hijo una versión con 50 o más, lo que acorta el tiempo hasta que aparecen los síntomas. Esto se llama **anticipación genética**.

&nbsp;

## Implementación del Trabajo

### Estructura del repositorio
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
> La carpeta `blast_db/` con la base de datos SwissProt no se incluye en el repositorio por su tamaño, pero tener en cuenta que es necesario para correrlo en local.

### Modo de uso de los scripts
Todos los scripts deben ejecutarse siempre desde la raíz del repositorio, no desde dentro de la carpeta `scripts/`. Esto garantiza que las rutas de input y output funcionen correctamente. Ejemplo de uso correcto:
```bash
cd tp-introduccion-bioinformatica
perl scripts/Ex1.pl data/sequence.gb
```

&nbsp;

## Ejercicio 1: Procesamiento de secuencias

### ¿Qué hace el script?
El script `Ex1.pl` lee un archivo GenBank de un mRNA de referencia, genera los 6 marcos de lectura posibles, traduce cada uno a su secuencia de aminoácidos, e identifica cuál es el marco de lectura correcto usando la anotación CDS que viene incluida en el archivo GenBank. Los resultados se escriben en un archivo FASTA.

### ¿Qué es un marco de lectura?
El mRNA es una cadena de nucleótidos que se lee de a 3 (codones) para producir aminoácidos. Dependiendo desde qué posición se empiece a leer, se obtienen proteínas completamente distintas. Hay 3 posiciones posibles en la hebra directa (+1, +2, +3) y 3 en la hebra complementaria inversa (-1, -2, -3), dando un total de 6 marcos de lectura posibles. Solo uno de ellos es el correcto y produce la proteína real.

En los marcos incorrectos aparecen codones de stop prematuros (representados como `*`) que cortarían la traducción antes de tiempo. El marco correcto es el que produce una secuencia de aminoácidos continua y coherente, sin esos stops prematuros. Lo que buscamos dentro de ese marco se llama **ORF (Open Reading Frame)**: una región que empieza con un codón de inicio (ATG) y termina con un codón de stop. Encontrar el ORF correcto es equivalente a encontrar el marco de lectura que usa la naturaleza para fabricar esa proteína.

### Ejecución
```bash
perl scripts/Ex1.pl data/sequence.gb
```

### Output generado
```
data/ORFs_HTT.fas
```

### Resultado
El marco de lectura correcto identificado es el +2, lo que significa que la secuencia codificante del gen HTT comienza en la segunda posición del mRNA. Esto coincide con la anotación CDS del archivo GenBank de referencia.

Para los ejercicios siguientes se extrajo solo la secuencia del marco correcto con este comando:
```bash
awk '/^>.*frame_\+2_MARCO_CORRECTO/{p=1} /^>/ && !/frame_\+2_MARCO_CORRECTO/{p=0} p' data/ORFs_HTT.fas > data/HTT_correcto.fas
```

Lo que hace es buscar en el archivo FASTA la secuencia que tiene `frame_+2_MARCO_CORRECTO` en su nombre y guardarla en un archivo nuevo llamado `HTT_correcto.fas`.

&nbsp;

## Ejercicio 2.a: BLAST

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

## Ejercicio 2.b: Interpretación del resultado del BLAST

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

