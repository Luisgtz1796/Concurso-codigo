close all %Para eliminar todas las ventanas que abrimos 
clear all %Para cerrar todas las variables 
clc %Para limpiar la ventana de comandos
%% Abrimos la imagen
im = dicomread("PF_Varian.dcm"); %Aquí leemos la imagen DICOM 
info=dicominfo('PF_Varian.dcm'); %En caso que se use, extraer datos del DICOM para generar un PDF
I = dicomread(info); 
%imshow(im)

%% Parametros de la imagen. Tamaño de pixel, columnas, filas
tamx=size(im); %[filas, columnas, canal]. En caso de una imagen DICOM, no tiene canal. 
tamx= tamx(2)/2; %numero de columnas
tmx=tamx*2; %para el numero de columnas de la matriz
tpix=1/2.56; %resolucion de portal en mm, tamaño de pixel de la imagen
lv=size(im); %longitud vertical en pixeles
lv=lv(1)/2;  %mitad de la longitud vertical
lh=size(im); %longitud horizontal en pixeles
lh=lh(2)/2; %mitad de la longitud horizontal
tmy=lv*2; % para el número de renglones 
%% Vectores de distancia (Tamaño de imagen)
mitad=floor(tamx);
zz=tpix:tpix:tpix*tmx; %Creamos un vector donde cada elemento de pixel en mm en pasos de tpix
zzt= zz-zz(1,mitad); %Los elementos de ese vector se nombran de -mitad, ..., 0, .... mitad  OJO ESTE LO USAREMOS DESPUES
ww=tpix:tpix:tpix*lv*2; %analogamente para l
%% Conversión de imagenes
imrt35=cast(im,'double'); %convertir los elementos de la matriz a formato double

%% Crear perfiles para identificar las franjas de mayor transmisión
pf1=imrt35(:,30); %(filas,columnas)
pf2=imrt35(:,lh-30); %Tomamos 3 perfiles, uno central y dos laterales, 30 pixeles a los laterales de los bordes del lado horizontal
pf3=imrt35(:,lh/2);
pff=smoothdata(((pf1+pf2+pf3)),'gaussian',20); %Agregamos un filtro de 20 gaussiano para eliminar el ruido, haciendo un promedio de los 3 perfiles
%plot(ww,pff,'k') %Esta linea es para visualizar el perfil obtenido
%% Localización de los bordes entre hojas vecinas continuas en cada franja
pfft=pff; %perfil de prueba
[zpkst,zlocst]=[findpeaks(pfft,ww/tpix)]; %[zpks se refiere a las alturas de los picos, zlocs las posiciones en el vector zz]
zpksmean= mean(zpkst); % para elegir los picos que sean más altos que un umbral (en este caso el promedio)
[zpks,zlocs]=[findpeaks(pff,ww/tpix, 'MinPeakHeight',zpksmean)]; %nuevamente que encuentre los picos pero a una altura más allá de zpksmean
% figure();
% plot(ww/tpix,pff,'k') %para encontrar las posiciones y picos
% hold on;
fance_pos=zlocs; %Definimos esta variable para usarla despues
 fance_pm= ceil((zlocs(2)-zlocs(1))/2); %punto medio entre cada franja
% plot(zlocs, zpks,'o')
% findpeaks(pff,ww/tpix, 'MinPeakHeight',zpksmean)
%% Perfiles sobre las franjas
i=0; %contador para el if 
zmatr3= zeros(length(zlocs),tmx); %creamos una matriz donde guadaremos todos los valores de las franjas encontradas, de dimensiones (longitud(zlocs), tmx)
zmatr2s= zeros(length(zlocs),tmx); %otra matriz para suavizar
c=0;
for i=1:1:length(zlocs) %va a leer desde el primer perfil hasta la long de zlocss
    pos=zlocs(i);  %guadar la posición de la i-franja
    c=c+1;
    zmatr3(c,:)=imrt35(pos,:); %rellena los valores en la franja c todos los valores de posicion y valor de pixel
    zmatr2s(c,:)= imrt35(pos-30,:); %analogo pero se realiza un suavizado 
    %zmatr2s(c,:)= smoothdata((100*imrt35(pos,:)),"gaussian",13); %analogo pero se realiza un suavizado 
    %ESTA ULTIMA LINEA PARA QUE ERA?
end 
%zmatr3, tiene los valores de las columnas en las filas que estan los picos
%de las franjas. 
%Al final obtendremos una matriz donde estan guardadas los valores de pix
%en cada uno de los perfiles (franjas)
%zmatr3=zmatr2s;
%% Promedio de la cada uno de los renglones y obtener un perfil final
zmatrizsum=(sum(zmatr3)+sum(zmatr2s));
zmatrizprom= mean2(zmatrizsum);
%zmatrizsum1=smoothdata((zmatrizsum-zmatrizprom),"gaussian",5);
zmatrizsum1=smoothdata(((zmatrizsum-zmatrizprom)/zmatrizprom),"gaussian",1);
zmatrizabs=abs(zmatrizsum1);
zmatrizsum1= zmatrizsum1+zmatrizabs;
%zmatrizsum2=smoothdata((zmatrizsum)/zmatrizprom*100,"movmean",100);
%zmatrizsum=zmatrizsum1+zmatrizsum2; 
%zmatrizprom= mean(zmatrizsum);
zdmin= 10; %es la distancia minima para encontrar el siguiente pico, porque es el grosor de la hoja mas chica =5mm
% plot(zz,zmatrizsum1,'k'); %para encontrar las posiciones y picos
[xpkst,xlocst,w,p]=[findpeaks(zmatrizsum1,zz/tpix)]; %nuevamente que encuentre los picos 
zaltur=mean(p);
%alturprom= median(w);
[xpks,xlocs]=[findpeaks(zmatrizsum1,zz/tpix, 'MinPeakDistance',zdmin, "MinPeakHeight",0)]; %nuevamente que encuentre los picos 
%findpeaks(zmatrizsum1,zz/tpix,'MinPeakDistance',zdmin, "MinPeakHeight",0);
% text(xlocs+.02,xpks,num2str((1:numel(xpks))'))
%findpeaks(zmatrizsum,zz/tpix,'MaxPeakWidth', alturprom)
%% Obtengamos los picos verdaderos y encontrar la mitad de las hojas
% Define la diferencia mínima requerida
delta = 0.04; %Este valor es mi referencia para elegir picos con altura delta
picos_true = []; % Aqui vamos a guardar los picos y sus posiciones que entraran al for
pos_true = [];
for i = 1:numel(xpks) %Va a leer todos los elementos de xpks encontrados anteriormente
    pico = xpks(i);
    posicion = xlocs(i);
    % Calcula la diferencia entre el pico y elementos a 3 posiciones a la izquierda y derecha
    if posicion > 3 && posicion < numel(zmatrizsum1) - 3
        dif_izq = pico - zmatrizsum1(posicion - 3);
        dif_der = pico - zmatrizsum1(posicion + 3);
        % Verifica si ambas diferencias son mayores o iguales a 0.4
        if abs(dif_izq) >= delta && abs(dif_der) >= delta
            % Guarda el pico y su posición original
            picos_true = [picos_true, pico];
            pos_true = [pos_true, zz(posicion)/tpix];
        end
    end
end
% figure();
% plot(zz/tpix,zmatrizsum1,'k'); %para encontrar las posiciones y picos
% hold on;
% plot(pos_true,picos_true,'o'); %para verificar que son los correctos
%% Aquí vamos a tomar los perfiles horizontales a la mitad de las hojas
leaf_pos = zeros(1, numel(pos_true) - 1); %Creamos un vector de para que guarde las posiciones, desde 1 hasta n-dim
for i = 1:numel(pos_true) - 1
    % Calcula el punto medio entre los elementos consecutivos
    pm = (pos_true(i) + pos_true(i + 1)) / 2;
    leaf_pos(i) = pm;
end
%A continuación, agregaremos la 1era y ultima hoja tomando en cuenta la
%distancia de 2da y 3era hoja 
delta2= leaf_pos(2)-leaf_pos(1);
leaf1= leaf_pos(1)-delta2;
leafn= leaf_pos(numel(leaf_pos))+delta2;
leaf_pos= [leaf1, leaf_pos, leafn];

figure;
imshow(im)
%imcontrast()
hold on 
xline(leaf_pos, 'y')
yline(fance_pos, 'r' )
[X, Y] = meshgrid(leaf_pos, fance_pos);
scatter(X(:), Y(:), 20, 'b', 'filled');
xline(pos_true, 'r')
xline(xlocs, 'w')

ax = gca;
exportgraphics(ax,'Setup.jpg','Resolution',300) %AQUI CAMBIE EL NOMBRE DE LA IMAGEN
close
%% Aqui crearemos los perfiles para obtener la distancia del gap entre hojas

[nr,nc]=size(X); 
c=zeros(nr,(fance_pm*2)+1); %Matriz para guardar los perfiles
za=cell(nr,2); %Celda que contiene la intensidad del perfil, za{i} es para el perfil izquierdo, mientras que za{i,2} es para el lado derecho
b=0; %contador
for i=1:nr %For para obtener la intensidad de pixel en el gap de las hojas en los lados izquierdo y derecho
    dizq=Y(i)-fance_pm:Y(i);
    dder=Y(i):Y(i)+fance_pm;
    zdelta=Y(i)-20:Y(i)+20;
    for j=1:length(X)
        for k=1:length(dizq)
            zbizq(k)=imrt35(dizq(k),round(X(1,j)));
            zbizq=smoothdata(zbizq,"gaussian",3);
        end
        for k=1:length(dder)
            zbder(k)=imrt35(dder(k),round(X(1,j)));
            zbder=smoothdata(zbder,"gaussian",3);
        end
        for k=1:length(zbizq)
            cizq(j,k)=(zbizq(k)*100)/max(zbizq);
        end
        for k=1:length(zbder)
            cder(j,k)=(zbder(k)*100)/max(zbder);
        end
    end
    za{i}=cizq;
    za{i,2}=cder;
end
%% Aqui obtenemos el valor donde el perfil obtiene un valor del 50% en el lado izq y derecho
pos_izq=zeros(60,1); %Vector que contiene las posiciones al 50% de la intensidad en el lado izquierdo
pos_der=zeros(60,1);  %Vector que contiene las posiciones al 50% de la intensidad en el lado derecho
vector_izq=-7.8:0.39:0; %Vector de distancia para hacer la interpolación en el lado izquierdo
vector_der=0:0.39:7.8; %Vector de distancia para hacer la interpolación en el lado derecho
zab=cell(nr,2); %contendra la distancia a la cual se encuentra la caida del 50% de intensidad. zab{i} es para el perfil izquierdo, mientras que zab{i,2} es para el lado derecho
%Aqui se supone que es un EPID modelo AS1000 con resolución de 0.39 mm/pp
for i=1:nr %For para obtener las posiciones al 50%
    for j=1:length(X)
        pos_izq(j)=interp1(za{i}(j,:),vector_izq,50);
        pos_der(j)=interp1(za{i,2}(j,:),vector_der,50);
    end
    zab{i}=pos_izq;
    zab{i,2}=pos_der;

end
%% Promedios de posición al 50% para cada  hoja en las n franjas
promf=cell(nr,2);
for i=1:nr
    promf{i}=mean(zab{i});
    promf{i,2}=mean(zab{i,2});
end
%% Promedio de gap obtenido a partir de las n franjas
zxizq=0;
zxder=0;
for i=1:nr
    zxizq=promf{i}+zxizq;
    zxder=promf{i,2}+zxder;
end
pfizq=zxizq/nr; %Promedio de apertura del gap en el lado izq
pfder=zxder/nr; %Promedio de apertura del gap en el lado der
%% Desviaciones del banco izquierdo
maxdesvizqi=0;
hmaxdesvizqi=0;
frmaxdesvizqi=0;
maxdesvderi=0;
hmaxdesvderi=0;
frmaxdesvderi=0;
cbi=0; %Cuenta el número de hojas desviadas en el banco izquierdo
vhdi=0; %Cuenta la desviacion de todas las hojas en el banco izquierdo

for i=1:nr
    for j=1:length(X)
        if zab{i}(j) > pfizq+1
            desv=abs(pfizq-zab{i}(j))-1;
            if desv>maxdesvizqi
                maxdesvizqi=desv;
                hmaxdesvizqi=j;
                frmaxdesvizqi=i;
            end
            % disp(['La hoja ', num2str(j), ' del banco izquierdo en la franja ', num2str(i), ' esta desviada ', num2str(desv), ' mm a la derecha respecto a la tolerancia']);
            cbi=cbi+1;
            vhdi=vhdi+desv;

        end
    end
    for j=1:length(X)
        if zab{i}(j) < pfizq-1
            desv=abs(zab{i}(j)-pfizq)-1;
            if desv>maxdesvderi
                maxdesvderi=desv;
                hmaxdesvderi=j;
                frmaxdesvderi=i;
            end
            % disp(['La hoja ', num2str(j), ' del banco izquierdo en la franja ', num2str(i), ' esta desviada ', num2str(desv), ' mm a la izquierda respecto a la tolerancia']);
            cbi=cbi+1;
            vhdi=vhdi+desv;
        end
    end
end
%disp(['*La hoja ', num2str(hmaxdesvizqi), ' del banco izquierdo en la franja ', num2str(frmaxdesvizqi), ' tiene la maxima desviación siendo de ', num2str(maxdesvizqi), ' mm']);
%disp(['*La hoja ', num2str(hmaxdesvderi), ' del banco derecho en la franja ', num2str(frmaxdesvderi), ' tiene la maxima desviación siendo de ', num2str(maxdesvderi), ' mm respecto a la tolerancia']);
%% Desviaciones del banco derecho
%Los for son para desplegar la hoja que tiene una desviación mayor a +-1mm
%de cada lado en el banco derecho
maxdesvizqd=0;
hmaxdesvizqd=0;
frmaxdesvizqd=0;
maxdesvderd=0;
hmaxdesvderd=0;
frmaxdesvderd=0;
cbd=0; %Cuenta el número de hojas desviadas en el banco izquierdo
vhdd=0; %Cuenta la desviacion de todas las hojas en el banco izquierdo
for i=1:nr
    for j=1:length(X)
        if zab{i,2}(j) > pfder+1
            desv=abs(zab{i,2}(j)-pfder)-1;
            if desv>maxdesvizqd
                maxdesvizqd=desv;
                hmaxdesvizqd=j;
                frmaxdesvizqd=i;
            end
            % disp(['La hoja ', num2str(j), ' del banco derecho en la franja ', num2str(i), ' esta desviada ', num2str(desv), ' mm a la derecha respecto a la tolerancia']);
            cbd=cbd+1;
            vhdd=vhdd+desv;
        end
    end
    for j=1:length(X)
        if zab{i,2}(j) < pfder-1
            desv=abs(pfder-zab{i,2}(j))-1;
            if desv>maxdesvderd
                maxdesvderd=desv;
                hmaxdesvderd=j;
                frmaxdesvderd=i;
            end
            % disp(['La hoja ', num2str(j), ' del banco derecho en la franja ', num2str(i), ' esta desviada ', num2str(desv), ' mm a la izquierda respecto a la tolerancia']);
            cbd=cbd+1;
            vhdd=vhdd+desv;
        end
    end
end
%% Obtenemos la hoja con maxima desviación
maxdes=0; %Valor de de la hoja con maxima desviación
hmd=0; %Numero de hoja con la maxima desviación
fmd=0; %Numero de franja con la hoja de maxima desviación
if maxdesvizqi>maxdes
    maxdes=maxdesvizqi;
    hmd=hmaxdesvizqi;
    fmd=frmaxdesvizqi;
    banco='izquierdo';

end
if maxdesvizqd>maxdes
    maxdes=maxdesvizqd;
    hmd=hmaxdesvizqd;
    fmd=frmaxdesvizqd;
    banco='derecho';
end
if maxdesvderd>maxdes
    maxdes=maxdesvderd;
    hmd=hmaxdesvderd;
    fmd=frmaxdesvderd;
    banco='derecho';
end
if maxdesvderi>maxdes
    maxdes=maxdesvderi;
    hmd=hmaxdesvderi;
    fmd=frmaxdesvderi;
    banco='izquierdo';
end
%% Aqui vamos a encontrar el par de hojas con maxima desviación
pmaxdes=0; %En esta variable se guardara el valor del par de hojas con maxima desviacion
phmd=0; %Aqui guardamos el par de hojas con maxima desviacipin
pfmd=0; %Aqui guardamos la franja en la que esta el par mas desviado
chd=0; %Contador de hojas desviadas
for i=1:nr
    for j=1:length(X)
        if zab{i}(j) > pfizq+1 && zab{i,2}(j) < pfder-1
            desvi1=abs(pfizq-zab{i}(j))-1;
            desvd1=abs(pfder-zab{i,2}(j))-1;
            desv1=desvi1+desvd1;
            chd=chd+1;
            if desv1>pmaxdes
                pmaxdes=desv1;
                phmd=j;
                pfmd=i;
            end
            %disp(['El par de hojas ', num2str(j), ' en la franja ', num2str(i), ' estan desviadas']);         
        end

        if zab{i}(j) > pfizq+1 && zab{i,2}(j) > pfder+1
            desvi1=abs(pfizq-zab{i}(j))-1;
            desvd2=abs(zab{i,2}(j)-pfder)-1;
            desv2=desvi1+desvd2;
            chd=chd+1;
            if desv2>pmaxdes
                pmaxdes=desv2;
                phmd=j;
                pfmd=i;
            end
            %disp(['El par de hojas ', num2str(j), ' en la franja ', num2str(i), ' estan desviadas']);
        end

        if zab{i}(j) < pfizq-1 && zab{i,2}(j) < pfder-1
            desvi2=abs(zab{i}(j)-pfizq)-1;
            desvd1=abs(pfder-zab{i,2}(j))-1;
            desv3=desvi2+desvd1;
            chd=chd+1;
            if desv3>pmaxdes
                pmaxdes=desv3;
                phmd=j;
                pfmd=i;
            end
            %disp(['El par de hojas ', num2str(j), ' en la franja ', num2str(i), ' estan desviadas']);
        end

        if zab{i}(j) < pfizq-1 && zab{i,2}(j) > pfder+1
            desvi2=abs(zab{i}(j)-pfizq)-1;
            desvd2=abs(zab{i,2}(j)-pfder)-1;
            desv4=desvi2+desvd2;
            chd=chd+1;
            if desv4>pmaxdes
                pmaxdes=desv4;
                phmd=j;
                pfmd=i;
            end
            %disp(['El par de hojas ', num2str(j), ' en la franja ', num2str(i), ' estan desviadas']);
        end
    end
end
%% Mapa de color para desviaciones fuera de la tolerancia en el banco der e izq
zmftdti=zeros(nr,length(X)); %Matriz que contendra las desviaciones izquierdas del banco derecho
zmftdtd=zeros(nr,length(X)); %Matriz que contendra las desviaciones derechas del banco derecho 

zmftiti=zeros(nr,length(X)); %Matriz que contendra las desviaciones izquierdas del banco izquierdo
zmftitd=zeros(nr,length(X)); %Matriz que contendra las desviaciones derechas del banco izquierdo

%% for para llenar la matriz del banco derecho

for i=1:nr
    for j=1:length(X)
        if zab{i,2}(j) < pfder-1 %desviación a la izquierda
            desv=zab{i,2}(j)-pfder; %Valor siempre negativo
            zmftdti(i,j)=desv; %guardamos el valor en la matriz
        else
            zmftdti(i,j)=0; %si no se cumple, es porque esta en tolerancia
        end
    end
    for j=1:length(X)
        if zab{i,2}(j) > pfder+1 % Desviación a la derecha
            desv=zab{i,2}(j)-pfder; %valor siempre positivo
            zmftdtd(i,j)=desv;
        else
            zmftdtd(i,j)=0;
        end
    end
end
zmftd=zmftdtd+zmftdti; % MAtriz con todas las deviaciones del banco derecho

%% for para llenar la matriz del banco izquierdo

for i=1:nr
    for j=1:length(X)
        if zab{i}(j) < pfizq-1 %desviación a la izquierda
            desv=zab{i}(j)-pfizq; %Valor siempre negativo
            zmftiti(i,j)=desv; %guardamos el valor en la matriz
        else
            zmftiti(i,j)=0; %si no se cumple, es porque esta en tolerancia
        end
    end
    for j=1:length(X)
        if zab{i}(j) > pfizq+1 % Desviación a la derecha
            desv=zab{i}(j)-pfizq; %valor siempre positivo
            zmftitd(i,j)=desv;
        else
            zmftitd(i,j)=0;
        end
    end
end
zmfti=zmftitd+zmftiti; % MAtriz con todas las deviaciones del banco izquierdo
%% limites del plot
lsdft=max(zmftd); %limite superior de todas las columnas
lsdft=max(lsdft); %Limite superior global
lidft=min(zmftd); %Limite infierior de todas las columnas
lidft=min(lidft); %Limite inferior global
lsift=max(zmfti); %limite superior de todas las columnas
lsift=max(lsift); %Limite superior global
liift=min(zmfti); %Limite infierior de todas las columnas
liift=min(liift); %Limite inferior global
%% Mapa de color para desviaciones del banco izquierdo
zmii=zeros(nr,length(X));
zmid=zeros(nr,length(X));

zmdi=zeros(nr,length(X));
zmdd=zeros(nr,length(X));
%% Banco izquierdo (ABAJO)
zmii=zeros(nr,length(X));
zmid=zeros(nr,length(X));

for i=1:nr
    for j=1:length(X)
        if zab{i}(j) < pfizq %Desvio a la izquierda
            desv=zab{i}(j)-pfizq; %valor negativo
            zmii(i,j)=desv;
        else
            zmii(i,j)=0;
        end
    end
    for j=1:length(X)
        if zab{i}(j) > pfizq %Desvio a la derecha
            desv=zab{i}(j)-pfizq; %Valor positivo
            zmid(i,j)=desv;
        else
            zmid(i,j)=0;
        end
    end
end
zmi=zmid+zmii;

%% Banco derecho (BANCO B ARRIBA)
for i=1:nr
    for j=1:length(X)
        if zab{i,2}(j) < pfder %desviación a la izquierda
            desv=zab{i,2}(j)-pfder; %Valor negativo
            zmdi(i,j)=desv;
        else
            zmdi(i,j)=0;
        end
    end
    for j=1:length(X)
        if zab{i,2}(j) > pfder % Desviación a la derecha
            desv=zab{i,2}(j)-pfder; %Valor positivo
            zmdd(i,j)=desv;
        else
            zmdd(i,j)=0;
        end
    end
end
zmd=zmdd+zmdi;
%% limites del plot
lsd=max(zmd); %limite superior de todas las columnas
lsd=max(lsd); %Limite superior global
lid=min(zmd); %Limite infierior de todas las columnas
lid=min(lid); %Limite inferior global

lsi=max(zmi); %limite superior de todas las columnas
lsi=max(lsi); %Limite superior global
lii=min(zmi); %Limite infierior de todas las columnas
lii=min(lii); %Limite inferior global

%% Resultados
disp(['La hoja ', num2str(hmd), ' del banco ', banco, ' en la franja ', num2str(fmd), ' tiene la maxima desviación siendo de ', num2str(maxdes), ' mm respecto a la tolerancia de 1 mm.']);
disp(['En la prueba que consto de ', num2str(nr), ' franjas, hubo un total de ', num2str((cbi+cbd)), ' hojas desviadas, siendo la desviación media de: ', num2str(((vhdi+vhdd)/(cbi+cbd))),' mm respecto a la tolerancia de 1 mm'])
disp(['Hay ',num2str(chd), ' pares de hojas desviadas y el par de hojas con maxima desviación es el ', num2str(phmd), ' en la franja ', num2str(pfmd), ' con desviación de ', num2str(pmaxdes)]);

%% Plot desviaciones de cada hoja
% figure();
% subplot(2,1,1)
% imagesc(zmi);
% ax=gca;
% ax.CLim=[lii lsi];
% colormap(jet(8))
% title('Desviaciones para el banco de hojas izquierdo')
% xticks(0:5:60)
% colorbar
% 
% subplot(2,1,2)
% imagesc(zmd);
% ax=gca;
% ax.CLim=[lid lsd];
% colormap(jet(8))
% title('Desviaciones para el banco de hojas derecho')
% xticks(0:5:60)
% colorbar
%% Plot de hojas fuera de la tolerancia
% figure();
% subplot(2,1,1)
% imagesc(zmfti);
% ax=gca;
% ax.CLim=[liift lsift];
% colormap(jet(8))
% title('Hojas del banco izquierdo fuera de la tolerancia')
% xticks(0:5:60)
% colorbar
% 
% subplot(2,1,2)
% imagesc(zmftd);
% ax=gca;
% ax.CLim=[lidft lsdft];
% colormap(jet(8))
% title('Hojas del banco derecho fuera de la tolerancia')
% xticks(0:5:60)
% colorbar

%%
import mlreportgen.report.* 
import mlreportgen.dom.* 
nomrep = input('Ingrese el nombre con que se exportara el reporte (no olvide colocar al final <.pdf>): ', 's');
rpt = Report(nomrep,'pdf');
tp = TitlePage; 
tp.Title = 'Prueba picket fence'; 
%tp.Subtitle = 'Instituto Nacional de Pédiatria, Médica Sur'; 
tp.Author = 'Luis Gutiérrez Melgarejo, Alvaro Daniel Cruz Cortes'; 
append(rpt,tp); 


ch1 = Section; 
ch1.Title = 'Resultados'; 

sec0 = Section;
sec0.Title = 'Procesamiento de la imagen';
para = Paragraph(['En la carpeta que a ejecutado este programa se creo un archivo llamado Setup.jpg, en el cual se muestra una imagen de como fue procesada la imagen' ...
    'Primero se tomaron perfiles verticales para poder saber cuantas franjas tiene la prueba picket fence y posteriormente sobre el punto de maxima intensidad de cada franja' ...
    'se tomaron los perfiles para poder encontrar las posiciones de las imagenes y en los puntos azules que se observan se tomaron los perfiles para poder encontrar' ...
    'la caida al 50% de los perfiles que corresponden al tamaño del gap que tienen los pares de hojas']); 
append(sec0,para) 
append(ch1,sec0) 


sec1 = Section; 
sec1.Title = 'Desviación maxima encontrada'; 
para = Paragraph(['La hoja ', num2str(hmd), ...
    ' del banco ', banco, ' en la franja ', num2str(fmd), ' tiene la maxima desviación siendo de ', ...
    num2str(maxdes), ' mm respecto a la tolerancia de 1 mm.']); 
append(sec1,para) 
append(ch1,sec1) 

sec2= Section;
sec2.Title = 'Par de laminas con la máxima desviación'
para2 = Paragraph(['Hay ',num2str(chd), ' pares de hojas desviadas y el par de hojas con maxima desviación es el ', ...
    num2str(phmd), ' en la franja ', num2str(pfmd), ' con desviación de ', num2str(pmaxdes), ' mm.']); 
append(sec2,para2) 
append(ch1,sec2) 

sec3= Section;
sec3.Title = 'Desviación media'
para3 = Paragraph(['En la prueba que consto de ', num2str(nr), ' franjas, hubo un total de ', num2str((cbi+cbd)), ...
    ' hojas desviadas, siendo la desviación media de: ', num2str(((vhdi+vhdd)/(cbi+cbd))),' mm respecto a la tolerancia de 1 mm.']); 
append(sec3,para3) 
append(ch1,sec3) 

append(rpt,ch1)

close(rpt);
% % Invokes rptview method
% rptview(rpt);
% % Invokes rptview function  
% rptview(nomrep);

%%
import mlreportgen.report.* 
import mlreportgen.dom.* 
nomrep = 'extra1';
rpt = Report(nomrep,'pdf');


ch1 = Section; 
ch1.Title = 'Extra 1'; 



sec1 = Section; 
sec1.Title = 'Hojas que presentan desviaciones mayores a la tolerancia'; 
para = Paragraph(['En la figura observamos las hojas de cada banco que superan la tolerancia de desviación de 1 mm o -1 mm' ...
    'Las desviaciones que son positivas indican que las hojas se contrayendo (cierran el gap), mientras que las desviaciones' ...
    'negativas indican que las hojas se estan contrayendo (se abre el gap). Los valores igual a 0 en el mapa de colores, indica que ' ...
    'esas hojas estan dentro de la tolerancia']); 
fig=Figure();
figure();
subplot(2,1,1)
imagesc(zmfti);
ax=gca;
ax.CLim=[liift lsift];
colormap(jet(8))
title('Hojas del banco izquierdo fuera de la tolerancia')
xticks(0:5:60)
colorbar

subplot(2,1,2)
imagesc(zmftd);
ax=gca;
ax.CLim=[lidft lsdft];
colormap(jet(8))
title('Hojas del banco derecho fuera de la tolerancia')
xticks(0:5:60)
colorbar
fig.Snapshot.Caption = 'Mapas de color que indican las hojas que superan la tolerancia permitida';
fig.Snapshot.Height = '3in';
fig.Snapshot.Height = '5in';
% fig.Snapshot.ScaleToFit = true;



append(sec1,para) 
append(ch1,sec1) 
append(ch1,fig);
append(rpt,ch1)

close(rpt);
% % Invokes rptview method
% rptview(rpt);
% % Invokes rptview function  
% rptview(nomrep);

%%
import mlreportgen.report.* 
import mlreportgen.dom.* 
nomrep = 'extra2';
rpt = Report(nomrep,'pdf');


ch1 = Section; 
ch1.Title = 'Extra 2'; 



sec1 = Section; 
sec1.Title = 'Hojas que presentan desviaciones mayores a la tolerancia'; 
para = Paragraph(['En la figura observamos los valores de desviación que poseen las hojas en cada uno de los bancos.' ...
    'Los valores positivos en la desviación de las hojas indican que estas se estan contrayendo, es decir, que el gap se cierra, ' ...
    'mientras que cuando la desviación es negativa indica que las hojas se contraen y por lo tanto el gap es más grande de lo que deberia.']); 
fig=Figure();
subplot(2,1,1)
imagesc(zmi);
ax=gca;
ax.CLim=[lii lsi];
colormap(jet(8))
title('Desviaciones para el banco de hojas izquierdo')
xticks(0:5:60)
colorbar

subplot(2,1,2)
imagesc(zmd);
ax=gca;
ax.CLim=[lid lsd];
colormap(jet(8))
title('Desviaciones para el banco de hojas derecho')
xticks(0:5:60)
colorbar
fig.Snapshot.Caption = 'Mapas de color que indican las desviaciones de cada hojas para cada banco';
fig.Snapshot.Height = '3in';
fig.Snapshot.Height = '5in';
% fig.Snapshot.ScaleToFit = true;



append(sec1,para) 
append(ch1,sec1) 
append(ch1,fig);
append(rpt,ch1)

close(rpt);
% % Invokes rptview method
% rptview(rpt);
% % Invokes rptview function  
% rptview(nomrep);