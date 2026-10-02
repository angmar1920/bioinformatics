# 📄 Práctica 06: Docking Molecular con AutoDock4

> 📖 **Módulo 7** — [Volver al README del módulo](../README.md) | Sección relacionada: [§9 Docking molecular](../README.md#9-docking-molecular)

## Docking Molecular

El docking molecular es una técnica computacional que predice la afinidad de unión de los ligandos a las proteínas receptoras.  Existen varias herramientas computacionales y algoritmos disponibles para las técnicas de acoplamiento molecular, tanto comerciales como gratuitos. Estos programas y herramientas se han desarrollado y se utilizan actualmente en la investigación de fármacos y en el ámbito académico. Según Sahoo et al., algunos de los programas de acoplamiento más utilizados son AutoDock Vina, Discovery Studio, Surflex, AutoDock GOLD, Glide, MCDock, MOE-Dock, FlexX, DOCK, LeDock, rDock, ICM, Cdcker, LigandFit, FRED y UCSF Dock. Entre estos programas, AutoDock Vina, Glide y AutoDock GOLD han sido identificados como las opciones mejor valoradas con las mejores puntuaciones. Además, algunos de estos programas han sido eficaces a la hora de predecir desviaciones cuadráticas medias (RMSD) que oscilan entre 1,5 y 2 Å, dependiendo de las poses experimentales. Sin embargo, el acoplamiento de receptores flexibles, en concreto la flexibilidad de la espina dorsal del receptor, sigue siendo un reto para los programas de acoplamiento actuales.

La electrostática computacional del complejo ligando-receptor puede evaluarse, analizarse y predecirse mediante el estudio de acoplamiento (docking), como afirman Sahoo et al. Según Mohapatra et al, este estudio suele seguir dos pasos distintos. En segundo lugar, las conformaciones se clasifican según una función de puntuación. Los algoritmos de muestreo deben reproducir teóricamente los modos de unión experimentales, y las confirmaciones obtenidas deben clasificarse según una función de puntuación, según Dash et al. El enfoque de laboratorio seco ofrece una ventaja significativa sobre los estudios de laboratorio in vivo en términos de inversión de recursos y tiempo, como señalan Sahoo et al. y Pramanik et al. Nanda et al. explicaron que este enfoque predice la orientación del ligando en un complejo formado por el propio ligando con proteínas o enzimas. Además, la forma del complejo acoplado y la interacción electrostática cuantifican la interacción.

La utilidad del Molecular Docking en el descubrimiento y diseño de fármacos está bien establecida. Sin embargo, Tao et al. han informado recientemente de un aumento del interés por la aplicación de este método en la ciencia de los alimentos. En concreto, el acoplamiento molecular se está utilizando para autentificar las dianas moleculares de los nutracéuticos en el tratamiento de enfermedades.
[leer mas](https://doi.org/10.1038/s41598-023-40160-2)


### Ejercicio 1. Docking con AutoDock4.

Luego de dos años y numerosos intentos de modelación por homología y diseño de fármacos _in silico_, se ha logrado diseñar un fármaco, pero necesita evaluar la afinidad del ligando con el receptor. 

#### Setup

1. Antes de empezar, necesita instalar los siguientes programas en su ordenador:
    - **ChimeraX**. Puede descargarse [aquí](https://www.cgl.ucsf.edu/chimerax/download.html). ChimeraX es la nueva versión de Chimera, desarrollada por el mismo laboratorio (Resource for Biocomputing, Visualization, and Informatics, RBVI, UCSF). Reemplaza gradualmente a Chimera, ofreciendo mejor rendimiento, soporte para estructuras más grandes (incluidos mapas de crio-EM), una interfaz renovada y un sistema de "bundles" (complementos) instalables desde el propio programa, entre ellos uno para ejecutar AutoDock Vina directamente desde su interfaz. Usaremos ChimeraX en lugar de Chimera (ya discontinuado) para visualizar y preparar las estructuras en este ejercicio.
    - **Autodock**. Puede descargarse [aqui](https://autodock.scripps.edu/download-autodock4/). AutoDock 4.2.6 (la versión estable actual) ya cuenta con binarios nativos para Windows, Mac OS X y Linux, por lo que **ya no es necesario instalar Cygwin ni copiar los ejecutables a ninguna carpeta especial**: basta con instalarlo y, opcionalmente, añadir su carpeta de instalación al PATH del sistema para poder llamarlo desde cualquier terminal. Véanse las [instrucciones de instalación](http://autodock.scripps.edu/downloads/autodock-4-2-x-installation-on-windows) y el [manual](https://autodock.scripps.edu/wp-content/uploads/sites/56/2021/10/AutoDock4.2.6_UserGuide.pdf) para más información. Los archivos principales son Autodock y AutoGrid, necesarios para ejecutar el predocking (mapas de energía), el docking y el cálculo de la puntuación.
    - ~~**AutoDock tools (MGLtools)**~~. **Ya no la utilizaremos.** MGLTools/AutoDockTools era la interfaz gráfica clásica para preparar PDBQT y archivos de parámetros, pero está basada en Python 2, sin mantenimiento activo y su instalación suele dar problemas en sistemas operativos recientes (especialmente en macOS con Apple Silicon). En su lugar usaremos herramientas modernas, mantenidas activamente y 100% multiplataforma (Windows/Mac/Linux), instalables con `pip`/`conda`:
        - **[Meeko](https://github.com/forlilab/Meeko)** (`pip install meeko`): desarrollado por el propio laboratorio de AutoDock (Forli Lab, Scripps Research). Reemplaza a `prepare_ligand4.py` para convertir el ligando a `.pdbqt`.
        - **[Open Babel](https://openbabel.org/)** (`conda install -c conda-forge openbabel`, o `brew install open-babel` en Mac): herramienta universal de conversión de formatos químicos, útil aquí para preparar el receptor (`obabel ... -xr`) y como alternativa rápida para el ligando.
        - **[AutoDockTools_py3](https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3)** (`pip install git+https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3`): un port a Python 3 de los scripts de línea de comandos de AutoDockTools (`prepare_receptor4.py`, `prepare_gpf4.py`, `prepare_dpf4.py`), sin necesidad de instalar la suite gráfica completa de MGLTools. Lo usaremos para generar los archivos `.gpf`/`.dpf` que AutoDock4 necesita.

2. Ahora que ya tenemos los programas necesarios, vamos a descargar la proteina y el ligando que vamos a utilizar. La estructura cristalizada de Imipenem-hydrolyzing beta-lactamase SME-1 (1DY6) se puede descargar del [Protein Data Bank](http://www.rcsb.org). La estructura de imipenem-hydrolyzing beta-lactamase SME-1 descargada tiene dos cadenas, puede eliminar una pues el imipenem puede interactuar con cualquiera de las dos cadenas, hágalo en ChimeraX (`open 1dy6` y luego `delete /B`). Esta tiene 267 residuos de aminoácidos y el monómero tiene una resolución de 2,13 Angstroms. Puede obtener el ligando (Imipenem) de las dos mayores bases de datos, [drugbank](http://www.drugbank.ca/) o [pubchem](http://pubchem.ncbi.nlm.nih.gov/), descarguelo y luego abralo en ChimeraX y guárdelo en formato .pdb (`save ligand.pdb`).

    **Nota**: la referencia clásica de esta práctica sugería regenerar target.pdb y ligand.pdb con el Swiss-PdbViewer; hoy en día ese programa está descontinuado (sin actualizaciones desde ~2011) y ChimeraX cubre esa misma función (limpieza de estructura, adición de hidrógenos, exportación a PDB) de forma más moderna, por lo que puede omitir el Swiss-PdbViewer y hacer todo desde ChimeraX.

#### Preparación de archivos

Los siguientes pasos son críticos porque dictan el procedimiento para ejecutar AutoGrid y AutoDock y proporcionan parámetros de acoplamiento precisos. Los archivos coordinados y la información correspondiente deben crearse en un formato específico denominado PDBQT, que contiene tipos de átomos/vínculos, cargas atómicas parciales, etc. En esta versión actualizada, la limpieza visual se hace en ChimeraX y la conversión a PDBQT y los archivos de parámetros se generan por línea de comandos con Open Babel, Meeko y AutoDockTools_py3 (ninguno requiere MGLTools). En esta sección, limitaremos nuestro experimento de docking a la configuración por defecto.

3. Preparación del receptor (target):
   1. En ChimeraX, abra `target.pdb`, elimine H<sub>2</sub>O y heteroátomos ajenos a la proteína, y añada hidrógenos:
      ```
      open target.pdb
      delete solvent
      delete ligand
      addh
      save target_clean.pdb
      ```
   2. Convierta el receptor limpio a PDBQT. Puede usar Open Babel (rápido, trata la estructura como rígida con `-xr`):
      ```
      obabel target_clean.pdb -O target.pdbqt -xr
      ```
      o, si prefiere el flujo más apegado al original de AutoDock4 (asigna cargas de Gasteiger y tipos de átomo AD4 explícitamente), use el script `prepare_receptor4.py` de AutoDockTools_py3:
      ```
      prepare_receptor4.py -r target_clean.pdb -o target.pdbqt -A checkhydrogens
      ```

4. Preparación del ligando:
   1. Asegúrese de que el ligando esté en 3D con hidrógenos (puede exportarlo con conformación 3D desde PubChem/DrugBank en formato `.sdf`, o generarla con ChimeraX/Avogadro/Open Babel).
   2. Convierta el ligando a PDBQT con Meeko, que detecta automáticamente el árbol de torsiones y los átomos aromáticos (sustituye los pasos manuales de "Detect Root", "Set Number of Torsions" y "Aromaticity Criterion" que antes se hacían a mano en MGLTools):
      ```
      mk_prepare_ligand.py -i ligand.sdf -o ligand.pdbqt
      ```

Al finalizar estos pasos tendrá `target.pdbqt` y `ligand.pdbqt`, listos para calcular la rejilla de energía.

#### Preparación de los parameters de Grid

5. Determine el centro de la caja de búsqueda alrededor del ligando. En ChimeraX puede obtener el centroide del ligando cargado con:
   ```
   open ligand.pdbqt
   info centroid #1
   ```
   Anote las coordenadas X, Y, Z que devuelve.
6. Genere el archivo `dock.gpf` con el script `prepare_gpf4.py` de AutoDockTools_py3, indicando el centro obtenido en el paso anterior y un tamaño de caja (npts) adecuado para cubrir el sitio de unión (por defecto 60×60×60 puntos a 0.375 Å):
   ```
   prepare_gpf4.py -l ligand.pdbqt -r target.pdbqt -o dock.gpf -p gridcenter="X,Y,Z" -p npts="60,60,60"
   ```
   Esto genera el mismo archivo `dock.gpf` que antes se creaba manualmente desde el menú Grid de AutoDockTools (mapa de energía por tipo de átomo, mapa electrostático y mapa de desolvatación), pero de forma reproducible y sin GUI.

#### Preparación de los parameters para Docking

7. Genere el archivo `dock.dpf` con `prepare_dpf4.py`, indicando el receptor y ligando en formato PDBQT (por defecto usa el algoritmo genético Lamarckiano, LGA, igual que en el flujo original de AutoDockTools):
   ```
   prepare_dpf4.py -l ligand.pdbqt -r target.pdbqt -o dock.dpf
   ```
   Si desea ajustar parámetros de búsqueda (número de corridas, tamaño de población, etc.) puede editar `dock.dpf` directamente con un editor de texto, ya que es un archivo plano.

Ahora tienes todos los archivos necesarios para el docking (target.pdbqt, ligand.pdbqt, dock.gpf, dock.dpf).

#### Ejecutando autodock

Despues de descargar e instalar AutoDock (ya no se necesita Cygwin ni copiar ejecutables a ninguna carpeta: los binarios de AutoDock4/AutoGrid4 se ejecutan de forma nativa).

    Para usuarios de Windows, Start > Run y escriba "cmd.exe" luego escriba el commando: "C:\Program Files (x86)\The Scripps Research Institute\Autodock\4.2.6\autodock4.exe"
    
    Luego debería ver un mensaje como este:

```
C:\Users\mgl>"C:\Program Files\The Scripps
Research Institute\Autodock"\autodock4.exe
usage: AutoDock -p parameter_filename
-l log_filename
-k (Keep original residue numbers)
-i (Ignore header-checking)
-t (Parse the PDBQT file to check torsions,
then stop.)
-d (Increment debug level)
-C (Print copyright notice)
--version (Print autodock version)
--help (Display this message)
C:\Users\mgl>
```
    Para entornos operativos tipo Unix (Mac/Linux), basta con que el ejecutable esté en el PATH (por ejemplo en `/usr/local/bin`) o invocarlo con su ruta completa desde la terminal; no se requiere Cygwin.

Start > Run y escriba "cmd.exe", cambie su directorio de trabajo a ~Desktop\autodock (usando el comando cd ).

Escriba en la consola: autogrid4.exe -p dock.gpf -l dock.glg
Escriba en la consola: autodock4.exe -p dock.dpf -l dock.dlg

Esto llevará algún tiempo, dependiendo de la capacidad de tu CPU y de tu memoria.

El archivo dlg contiene toda la información sobre las ejecuciones de acoplamiento, la energía de enlace estimada en Kcal/mol y otra información como la RMSD frente a la pose de enlace del cristal.

#### Analizando los resultados

Como ya no usamos MGLTools, analizaremos el archivo `.dlg` sin su menú Analyze, con un script corto y multiplataforma (funciona igual en Windows/Mac/Linux porque solo requiere Python 3, ya instalado junto con Meeko/AutoDockTools_py3):

1. Los resultados de docking se encuentran en el archivo de registro `dock.dlg`. Revise la tabla de puntuaciones (energía de unión estimada y RMSD por corrida) buscando la sección `RANKING` dentro del archivo, por ejemplo abriéndolo con un editor de texto o con:
   ```
   python -c "print(open('dock.dlg').read()[open('dock.dlg').read().find('RANKING'):])"
   ```
2. Extraiga las poses acopladas (líneas `DOCKED`) a un archivo PDBQT independiente para visualizarlas:
   ```python
   with open("dock.dlg") as f, open("docked_poses.pdbqt", "w") as out:
       for line in f:
           if line.startswith("DOCKED"):
               out.write(line[8:])
   ```
3. Convierta las poses a PDB si lo prefiere (`obabel docked_poses.pdbqt -O docked_poses.pdb`) y ábralas en ChimeraX para recorrer cada conformación del ligando unido a la betalactamasa, tal como antes se hacía con Analyze > Conformations > Play. La mejor conformación tiene una energía de unión (ΔG) de -5,75 kcal/mol y una constante de inhibición (Ki) de 60,87 µM y una RMSD (desviación cuadrática media de las posiciones atómicas) de la estructura de referencia de 1,22 Å. Esto demuestra que los resultados de Autodock son fiables y precisos (en el rango nanomolar para un inhibidor conocido). El docking y el cribado virtual serían una baza importante para identificar nuevos inhibidores de BACE1.

### Ejercicio 2. Docking con AutoDock vina.

Se han desarrollado dos métodos de acoplamiento en paralelo, para responder a dos necesidades diferentes. El desarrollo comenzó con AutoDock, y sigue siendo la plataforma de experimentación en métodos de acoplamiento. AutoDock Vina se desarrolló más recientemente para satisfacer la necesidad de un método de docking llave en mano que no requiriera amplios conocimientos expertos por parte de los usuarios1. Está altamente optimizado para realizar experimentos de docking utilizando métodos por defecto bien probados. Ambos métodos están actualmente disponibles de forma gratuita. AutoDock Vina es rápido y eficaz para la mayoría de los sistemas, mientras que AutoDock está disponible para los sistemas que requieren mejoras metodológicas adicionales.

Ambos métodos están diseñados como herramientas genéricas de acoplamiento computacional, aceptan archivos de coordenadas para el receptor y el ligando y predicen conformaciones acopladas óptimas. Normalmente, los usuarios parten de coordenadas del receptor obtenidas por cristalografía o espectroscopia de RMN, y de coordenadas del ligando generadas a partir de cadenas SMILES u otros métodos.[leer mas](https://doi.org/10.1038/nprot.2016.051)

Para este ejercicio usaremos AutoDock Vina integrado directamente en ChimeraX, lo que simplifica enormemente el flujo de trabajo respecto al Ejercicio 1 (ya no es necesario usar MGLTools, Cygwin ni la línea de comandos por separado).

#### Setup

1. **AutoDock Vina**. Puede descargarse [aquí](https://vina.scripps.edu/downloads/). Instálelo y anote la ruta del ejecutable (`vina.exe` en Windows o `vina` en Mac/Linux), ya que ChimeraX se la pedirá la primera vez que lo use.
2. **ChimeraX** (ver sección de instalación arriba). Verifique que tiene una versión reciente (1.6 o superior), ya que el bundle de AutoDock Vina viene incluido de forma nativa a partir de esa versión.

#### Preparación del receptor y el ligando

3. Abra ChimeraX y cargue la estructura del receptor directamente desde el PDB con el comando:
   ```
   open 1dy6
   ```
4. Elimine la cadena duplicada, las moléculas de agua y los heteroátomos que no correspondan al ligando de interés:
   ```
   delete /B
   delete solvent
   delete ligand
   ```
5. Prepare el receptor con la herramienta **Dock Prep** (menú Tools > Structure Editing > Dock Prep, o el comando `dockprep`). Esto añade hidrógenos, asigna cargas parciales y completa cadenas laterales incompletas, de forma equivalente a los pasos manuales de AutoDockTools del Ejercicio 1.
6. Cargue el ligando (por ejemplo desde un archivo `.pdb`, `.mol2` o `.sdf` descargado de PubChem/DrugBank):
   ```
   open ligand.mol2
   ```
   y aplique también `dockprep` sobre el ligando para añadir hidrógenos y cargas.

#### Ejecutando AutoDock Vina desde ChimeraX

7. Abra la herramienta de docking desde el menú **Tools > Structure Analysis > AutoDock Vina**, o ejecute el comando:
   ```
   vina receptor #1 ligand #2
   ```
   donde `#1` y `#2` son los números de modelo del receptor y el ligando cargados.
8. En el panel de la herramienta (o mediante la opción `resMargin`/`binding site`), defina la región de búsqueda: puede seleccionar los residuos del sitio de unión conocido en el visor 3D y usar la opción "usar la selección actual" para centrar la caja de búsqueda automáticamente, en lugar de definir manualmente las coordenadas X, Y, Z como en el Ejercicio 1.
9. Indique la ruta del ejecutable de `vina` instalado en el paso 1 la primera vez que se le solicite (ChimeraX la recordará en usos posteriores).
10. Ejecute el docking. ChimeraX llama a AutoDock Vina en segundo plano y, al finalizar, carga automáticamente las poses generadas como modelos adicionales en el visor, ordenadas por afinidad de unión (kcal/mol), sin necesidad de manipular archivos `.pdbqt`, `.gpf` o `.dpf` manualmente.

#### Analizando los resultados

11. Use el panel **Log** o el comando `vina list` para ver la tabla de puntuaciones (afinidad estimada en kcal/mol) de cada pose generada.
12. Recorra las poses con el selector de modelos en el panel lateral (o el comando `modelinspection`) para visualizar cada conformación del ligando acoplada al receptor.
13. Compare la mejor pose obtenida con AutoDock Vina frente a la obtenida con AutoDock4 en el Ejercicio 1: calcule el RMSD entre ambas con el comando `rmsd` y discuta las diferencias en energía de unión y en la pose predicha.

