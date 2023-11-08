# Código Picket Fence (PF)
Autores: M. en C. Alvaro Daniel Cruz Cortes (FM, Médica Sur) y  M. en C. Luis Guitierrez Malgarejo (FM, Instituto Nacional de Pediatría)
         
Información del repositorio: Este repositorio tiene como finalidad albergar los códigos que surjan durante el concurso de la SMFM de PF.

Consiste en un par de códigos para evaluar cuantitativamente la prueba de Picket-Fance propuesta por el TG-142 de la AAPM. Se pueden evaluar imágenes de dosimetría portal DICOM y dosimetría con película radiográfica.  El código obtiene las desviaciones de las posiciones de las hojas respecto al eje central de cada franja. 

# Documentación
Los códigos fueron realizados con lenguaje MATLAB. Utilice la versión MATLAB R2022b en adelante para correr los programas correspondientes. 

# Archivos de entrada
Para el programa "PF_Portal", es necesario ingresar archivos DICOM de imagenes de dosimetría portal. En las siguientes líneas ingrese el nombre del archivo. No olvide agregar la terminación ".dcm". 

   
     %% Abrimos la imagen
     im = dicomread("PF_Varian.dcm"); %Aquí leemos la imagen DICOM 
     info=dicominfo('PF_Varian.dcm'); %En caso que se use, extraer datos del DICOM para generar un PDF
     I = dicomread(info); 
     %imshow(im)

Posteriormente, corra el programa con "RUN".


Para el programa "PF_Peli", es necesario ingresar imagenes JPG. Las condiciones de irradiación que recomendamos es: al menos 500 UM por cada franja, 600 UM/min, haz de 6X, una geometría de SSD=100 cm, colocando placas de 5cm por debajo y 1.5cm por encima de la película. Colocar 4 marcas de referencia en la película con ayuda del crosshair.  Escanear la pelicula en modo reflexión con los siguientes parámetros: archivo .jpg, 300 dpi, 48 bits de profunidad, En las siguientes líneas ingrese el nombre del archivo. No olvide agregar la terminación ".jpg". 

   
       % Apertura de la imagen y ajuste de constraste
        im=rgb2gray(imread('bk_300006.jpg'));
        imc=imread('bk_300006.jpg');

Posteriormente, corra el programa con "RUN".

En este programa, el usuario deberá seleccionar las marcas de plumon puestas en la película con el fin de centrar la imagen para obtener los perfiles. Siga las instrucciones que menciona en cada ventana que se desplegue. 

 # Resultados
En ambos programas, al terminar las iteraciones, le generará el conjunto de resultados que pide el concurso y un mensaje para ingresar el nombre del PDF con los resultados adecuados, no olvide agregar la terminación ".pdf". 


     La hoja 23 del banco izquierdo en la franja 7 tiene la maxima desviación siendo de 1.1627 mm respecto a la 
     tolerancia de 1 mm.
     En la prueba que consto de 10 franjas, hubo un total de 58 hojas desviadas, siendo la desviación media de: 
     0.30482 mm respecto a la tolerancia de 1 mm
     Hay 19 pares de hojas desviadas y el par de hojas con maxima desviación es el 23 en la franja 7 con desviación 
     de 1.5036
     Ingrese el nombre con que se exportara el reporte (no olvide colocar al final <.pdf>): prueba.pdf


Los programas generá 3 archivos pdf:

-"prueba" son los resultados necesarios que pide la convocatoria

-"extra1" son las hojas que presentan las posiciones y las hojas desviaciones mayores a la tolerancia

-"extra2" son las  hojas que presentan las hojas con desviaciones mayores a la tolerancia

Además, el código genera una imagen JPG llamado "Setup", donde se pueden observar los perfiles que fueron tomados para su evaluación (franjas y hojas).

# Notas
Para probar el código "PF_Portal", utilice el archivo "PF_Varian".
Para probar el código "PF_Peli", utilice el archivo "bk_300006".

# Discusión
Cualquier duda, anexamos nuestros correos:

-<alvarodan.crco@ciencias.unam.mx>

-<luisgtz1796@gmail.com>
