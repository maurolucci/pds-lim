# pds-lim

Requisitos

* Tener instalado gurobi. La instalación asume que se encuentra en /opt/gurobi1103
* Tener instalada la libreria program options de boost

Instrucciones de instalación

1. Crear directorio _deps
2. Clonar los repositorios fmt (https://github.com/fmtlib/fmt) y tinyxml2 (https://github.com/leethomason/tinyxml2)
3. Hacer cmake . en ambos directorios
4. Volver al directorio raiz y hacer make

Intrucciones de ejecución
cd experiments/
nohup ./script1.sh &