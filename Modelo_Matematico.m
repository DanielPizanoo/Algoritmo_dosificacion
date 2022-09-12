%Algoritmo Metodos de Minimos Cuadrados Doc Santiago

close all
clear all
clc
%Programa para el ajuste polinomial de datos
%Datos de entrada
xd=[0 7 21 35 49 63 77 91 115 129 143 157] % dias de muestreo
yd=[20 43.67 68.13 83.83 129.90 137.34 216.21 228.70 315.19 363.68 484 550] % peso medido

%U=length(xd);
%&orden del ajuste
orden=input('Defina el orden de la aproximacion: ')
%#datos en x en y
n=length(xd); %mide la longitud del vector de dias
%filas
for i=[1:orden+1]
    %columnas
    for j=[1:orden+1] %valor inicial
        a(i,j)=0;
        for r=[1:n] %suma de terminos por elemento
            a(i,j)=a(i,j)+xd(r).^(j+i-2); % matriz que contiene los coeficientes de la matriz A de mínimos cuadrados
        end
    end
end

for j=[1:orden+1]
    b(j)=0; % valor inicial para los coeficientes de mínimos cuadrados
    for r=[1:n]
        b(j)=b(j)+yd(r).*xd(r)^(j-1); % coeficientes del vector b de mínimos cuadrados.
    end
end

mw=a\b' % resuelve para calculas los coeficientes del polinomio de ajuste
%ws=@(x) mw(1)+mw(2)*x
%ws=@(x) mw(1)+mw(2)*x+mw(3)*x^2
ws=@(x) mw(1)+mw(2)*x+mw(3)*x^2+mw(4)*x.^3 % formula del polinomio de ajuste de orden 3
%ws=@(x) mw(1)+mw(2)*x +mw(3)*x^2+mw(4)*x^3+mw(5)*x^4;

h=1 %tamaño del intervalo

x = xd(1):h:xd(end); % vector que considera los dias de muestreo

U=length(x);
ys=[]; %inicializa el vector ys que servira para calcular los valores teoricos para el peso de los peces calculado por el polinomio de ajuste
for i=1:U
    ys(i)=ws(x(i)); % vector que contiene los valores teoricos para el peso, calulados por el polinomio de ajuste
end

%Calculo del error en el polinomio

z=[];
for i=1:length(xd)
    z(end+1)=ws(xd(i)); % calcula los valores del peso de los peces deacuerdo al polinomio de ajuste
end

Error=(yd'-z');
abError=abs(Error); % calcula el valor absoluto para el error
disp('El promedio del error es: ')
mean(abError)

plot(xd,z)

%---------------- r constantes ---------------------%
Wi=zeros; %inicializa la matriz Wi
p=1;
for r=0.1:0.1:3.0
    %PARÁMETROS DEL SISTEMA        
    n=0.81; m=0.67;  b=0.62; a=0.53;   %se asume b de otro articulo
    hache=0.8; kmin=0.00133; j=0.0132; s=21.38;
    %CONDICIONES AMBIENTALES DEL ESTANQUE
    %Temp=33; OD=4; A=0.83;  
    Temp=40; OD=0.3; A=0.83;  
    %CONDICIONES CRITICAS
    Tcri=41; ODcrit=1; Acri=0.001;
    %CONDICIONES MINIMAS
    Tmin= 15; ODmin=0.3; Amax=1.4;
    %DADO QUE EL ESTANQUE ES MUY GRANDE SE ASUME QUE LA PRESENCIA DE LA TILAPIA
    %NO AFECTA AL MEDIO ACUATICO
    tau=0.9948; delta=0.9429; niu=.4254; betha=0.7750; 
   
    %CONDICIONES DE OPERACION 
    t0=0; %condicion inicial de tiempo (eje x)
    tf=157; %tiempo final
    h=1;             % step size
    x = t0:h:tf;        % vector de dias 
    y = zeros(1,length(x)) ; % vector de pesos 
    y(1) = 20;   % valor inicial del peso
    
    %ECUACIONES ALGEBRAICAS
    k=kmin.*exp(j.*(Temp-Tmin));
    
    %EDO que describe tasa de crecimiento donde x es el peso o variable independiente   
    F_xy=@(t,y) ((1-a).*b.*tau.*delta.*niu*betha*(r).*y.^m)-k*y.^n;

    %Metodo de Runge-Kutta para resolver la EDO la resuelve 1 vez para cada valor de r
    for i=1:(length(x)-1)                     
        Wi(p,i)=y(i);
        k_1 = F_xy(x(i),y(i));
        k_2 = F_xy(x(i)+0.5*h,y(i)+0.5*h*k_1);
        k_3 = F_xy((x(i)+0.5*h),(y(i)+0.5*h*k_2));
        k_4 = F_xy((x(i)+h),(y(i)+k_3*h));
        y(i+1) = y(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;
        Wi(p,i+1)=y(i+1);
    end
    p=p+1; % aumenta en 1 el valor de p
    
    ew=length(Wi(:,1));
    aw=length(Wi(1,:));

    error=zeros;
    for z1=1:ew
        for z2=1:aw
            error(z1,z2)=abs(Wi(z1,z2)-ys(1,z2)); % calcula 1 matriz del error
        end
    end
end

%-------------- Modelo San Buenaventura -------------------%

dd=[];
for j=1:aw
    [M,I]=min(error(:,j)); % almacena el error para cada dia y cada una de las curvas verdes. 
    dd(end+1)=I;
end
dd=dd/10;

%PARÁMETROS DEL SISTEMA      % repeticion del bloque de datos de la linea 75  
n=0.81; m=0.67;  b=0.62; a=0.53;   %se asume b de otro articulo
hache=0.8; kmin=0.00133; j=0.0132; s=21.38;
%Condiciones ambientales
A=0.83;
%Condiciones críticas
Tcri=41; ODcrit=1; Acrit=0.001;
%Condiciones minimas, maximas y optimas
Tmin=15; ODmin=0.3; Amax=1.4; Topt=33;
betha=0.7750;

%CONDICIONES DE OPERACION 
t0=0; %condicion inicial de tiempo (eje x)
tf=157; %tiempo final
h=1; % step size
xs = t0:h:tf;         
yd = zeros(1,length(xs)); 
yd(1) = 20; %yb calcula el peso de los peces pero usando los valores de r que se calcularon para zz

if xs>=7 & xs<21
    T=32.6; OD=2.985;
    
    if T<Topt
        tau = exp(-4.6*(Topt-T)/(Topt-Tmin).^(4));
    elseif T>=Topt
        tau = exp(-4.6*(T-Topt)/(Tmax-Topt).^(4));
    end
    
    if OD>ODcrit
        delta=1;
    elseif ODmin <= OD <= ODcrit
        delta=(OD-ODmin)/(ODcrit-ODmin);
    elseif OD<ODmin
        delta=0.0;
    end
    
    if A<Acrit
        niu=1.0;
    elseif Acrit<=A<=Amax
        niu=(Amax-A)/(Amax-Acrit);
    elseif A>Amax
        niu=0.0;
    end
    
elseif xs>=21 & xs<35
    T=32.4; OD=2.025;
    
    if T<Topt
        tau = exp(-4.6*(Topt-T)/(Topt-Tmin).^(4));
    elseif T>=Topt
        tau = exp(-4.6*(T-Topt)/(Tmax-Topt).^(4));
    end
    
    if OD>ODcrit
        delta=1;
    elseif ODmin <= OD <= ODcrit
        delta=(OD-ODmin)/(ODcrit-ODmin);
    elseif OD<ODmin
        delta=0.0;
    end
    
    if A<Acrit
        niu=1.0;
    elseif Acrit<=A<=Amax
        niu=(Amax-A)/(Amax-Acrit);
    elseif A>Amax
        niu=0.0;
    end
    
elseif xs>=35 & xs<49
    T=31.25; OD=3.485;
    
    if T<Topt
        tau = exp(-4.6*(Topt-T)/(Topt-Tmin).^(4));
    elseif T>=Topt
        tau = exp(-4.6*(T-Topt)/(Tmax-Topt).^(4));
    end
    
    if OD>ODcrit
        delta=1;
    elseif ODmin <= OD <= ODcrit
        delta=(OD-ODmin)/(ODcrit-ODmin);
    elseif OD<ODmin
        delta=0.0;
    end
    
    if A<Acrit
        niu=1.0;
    elseif Acrit<=A<=Amax
        niu=(Amax-A)/(Amax-Acrit);
    elseif A>Amax
        niu=0.0;
    end
    
elseif xs>=49 & xs<63
    T=31.1; OD=4.2;
    
    if T<Topt
        tau = exp(-4.6*(Topt-T)/(Topt-Tmin).^(4));
    elseif T>=Topt
        tau = exp(-4.6*(T-Topt)/(Tmax-Topt).^(4));
    end
    
    if OD>ODcrit
        delta=1;
    elseif ODmin <= OD <= ODcrit
        delta=(OD-ODmin)/(ODcrit-ODmin);
    elseif OD<ODmin
        delta=0.0;
    end
    
    if A<Acrit
        niu=1.0;
    elseif Acrit<=A<=Amax
        niu=(Amax-A)/(Amax-Acrit);
    elseif A>Amax
        niu=0.0;
    end
    
elseif xs>=63 & xs<77
    T=29.85; OD=3.515;
    
    if T<Topt
        tau = exp(-4.6*(Topt-T)/(Topt-Tmin).^(4));
    elseif T>=Topt
        tau = exp(-4.6*(T-Topt)/(Tmax-Topt).^(4));
    end
    
    if OD>ODcrit
        delta=1;
    elseif ODmin <= OD <= ODcrit
        delta=(OD-ODmin)/(ODcrit-ODmin);
    elseif OD<ODmin
        delta=0.0;
    end
    
    if A<Acrit
        niu=1.0;
    elseif Acrit<=A<=Amax
        niu=(Amax-A)/(Amax-Acrit);
    elseif A>Amax
        niu=0.0;
    end
    
elseif xs>=77 & xs<91
    T=30.1; OD=2.18;
    
    if T<Topt
        tau = exp(-4.6*(Topt-T)/(Topt-Tmin).^(4));
    elseif T>=Topt
        tau = exp(-4.6*(T-Topt)/(Tmax-Topt).^(4));
    end
    
    if OD>ODcrit
        delta=1;
    elseif ODmin <= OD <= ODcrit
        delta=(OD-ODmin)/(ODcrit-ODmin);
    elseif OD<ODmin
        delta=0.0;
    end
    
    if A<Acrit
        niu=1.0;
    elseif Acrit<=A<=Amax
        niu=(Amax-A)/(Amax-Acrit);
    elseif A>Amax
        niu=0.0;
    end
    
elseif xs>=91 & xs<105
    T=28.45; OD=2.82;
    
    if T<Topt
        tau = exp(-4.6*(Topt-T)/(Topt-Tmin).^(4));
    elseif T>=Topt
        tau = exp(-4.6*(T-Topt)/(Tmax-Topt).^(4));
    end
    
    if OD>ODcrit
        delta=1;
    elseif ODmin <= OD <= ODcrit
        delta=(OD-ODmin)/(ODcrit-ODmin);
    elseif OD<ODmin
        delta=0.0;
    end
    
    if A<Acrit
        niu=1.0;
    elseif Acrit<=A<=Amax
        niu=(Amax-A)/(Amax-Acrit);
    elseif A>Amax
        niu=0.0;
    end
    
elseif xs>=115 & xs<129
    T=28.95; OD=1.45;
    
    if T<Topt
        tau = exp(-4.6*(Topt-T)/(Topt-Tmin).^(4));
    elseif T>=Topt
        tau = exp(-4.6*(T-Topt)/(Tmax-Topt).^(4));
    end
    
    if OD>ODcrit
        delta=1;
    elseif ODmin <= OD <= ODcrit
        delta=(OD-ODmin)/(ODcrit-ODmin);
    elseif OD<ODmin
        delta=0.0;
    end
    
    if A<Acrit
        niu=1.0;
    elseif Acrit<=A<=Amax
        niu=(Amax-A)/(Amax-Acrit);
    elseif A>Amax
        niu=0.0;
    end
    
elseif xs>=129 & xs<143
    T=29.15; OD=2.08;
    
    if T<Topt
        tau = exp(-4.6*(Topt-T)/(Topt-Tmin).^(4));
    elseif T>=Topt
        tau = exp(-4.6*(T-Topt)/(Tmax-Topt).^(4));
    end
    
    if OD>ODcrit
        delta=1;
    elseif ODmin <= OD <= ODcrit
        delta=(OD-ODmin)/(ODcrit-ODmin);
    elseif OD<ODmin
        delta=0.0;
    end
    
    if A<Acrit
        niu=1.0;
    elseif Acrit<=A<=Amax
        niu=(Amax-A)/(Amax-Acrit);
    elseif A>Amax
        niu=0.0;
    end
    
elseif xs>=143 & xs<157
    T=28.45; OD=2.16;
    
    if T<Topt
        tau = exp(-4.6*(Topt-T)/(Topt-Tmin).^(4));
    elseif T>=Topt
        tau = exp(-4.6*(T-Topt)/(Tmax-Topt).^(4));
    end
    
    if OD>ODcrit
        delta=1;
    elseif ODmin <= OD <= ODcrit
        delta=(OD-ODmin)/(ODcrit-ODmin);
    elseif OD<ODmin
        delta=0.0;
    end
    
    if A<Acrit
        niu=1.0;
    elseif Acrit<=A<=Amax
        niu=(Amax-A)/(Amax-Acrit);
    elseif A>Amax
        niu=0.0;
    end
    
else
    T=27.26; OD=3.83;
    
    if T<Topt
        tau = exp(-4.6*(Topt-T)/(Topt-Tmin).^(4));
    elseif T>=Topt
        tau = exp(-4.6*(T-Topt)/(Tmax-Topt).^(4));
    end
    
    if OD>ODcrit
        delta=1;
    elseif ODmin <= OD <= ODcrit
        delta=(OD-ODmin)/(ODcrit-ODmin);
    elseif OD<ODmin
        delta=0.0;
    end
    
    if A<Acrit
        niu=1.0;
    elseif Acrit<=A<=Amax
        niu=(Amax-A)/(Amax-Acrit);
    elseif A>Amax
        niu=0.0;
    end
end

%ECUACIONES ALGEBRAICAS
k=kmin.*exp(j.*(T-Tmin));

%Dosificación de alimento (%peso de alimento respecto al peso del animal) con base al tiempo en el interior del estanque

ii=[dd(1)];
%EDO que describe tasa de crecimiento donde x es el peso o variable independiente   
aa=1;
ee=[];
%Metodo de Runge-Kutta para resolver la EDO

for jj=1:(length(xs)-1)
    rb=dd(jj);
    gr=(aa*rb);
    ee(1)=abs(yd(1)-ys(1));
    Fb_xy=@(t,yd) ((1-a).*b.*tau.*delta.*niu*betha*(gr).*yd.^m)-k*yd.^n;
    ks_1 = Fb_xy(xs(jj),yd(jj));
    ks_2 = Fb_xy(xs(jj)+0.5*h,yd(jj)+0.5*h*ks_1);
    ks_3 = Fb_xy((xs(jj)+0.5*h),(yd(jj)+0.5*h*ks_2));
    ks_4 = Fb_xy((xs(jj)+h),(yd(jj)+ks_3*h));
    yd(jj+1) = yd(jj) + (1/6)*(ks_1+2*ks_2+2*ks_3+ks_4)*h;
    ii(end+1)=gr;
    ee(end+1) = abs(yd(jj+1)-ys(jj+1));
    if ee(jj) > 0.1*ys(jj)
        aa = 0.70;
    else
        aa = 1;
    end
end

%-------------- Modelo Doc Santiago -------------------%

zz=[];
for j=1:aw
    [M,I]=min(error(:,j)); % almacena el error para cada dia y cada una de las curvas verdes. 
    zz(end+1)=I;
end
zz=zz/10;

%PARÁMETROS DEL SISTEMA      % repeticion del bloque de datos de la linea 75  
n=0.81; m=0.67;  b=0.62; a=0.53;   %se asume b de otro articulo
hache=0.8; kmin=0.00133; j=0.0132; s=21.38;
%CONDICIONES AMBIENTALES DEL ESTANQUE
Temp=33; OD=4; A=0.83;  
%CONDICIONES CRITICAS
Tcri=41; ODcrit=1; Acri=0.001;
%CONDICIONES MINIMAS
Tmin= 15; ODmin=0.3; Amax=1.4;
%DADO QUE EL ESTANQUE ES MUY GRANDE SE ASUME QUE LA PRESENCIA DE LA TILAPIA
%NO AFECTA AL MEDIO ACUATICO
tau=0.9948; delta=0.9429; niu=.4254; betha=0.7750; 
%ECUACIONES ALGEBRAICAS
k=kmin.*exp(j.*(Temp-Tmin));

%CONDICIONES DE OPERACION 
t0=0; %condicion inicial de tiempo (eje x)
tf=157; %tiempo final
h=1; % step size
xb = t0:h:tf;         
yb = zeros(1,length(xb)); 
yb(1) = 20; %yb calcula el peso de los peces pero usando los valores de r que se calcularon para zz

%Dosificación de alimento (%peso de alimento respecto al peso del animal) con base al tiempo en el interior del estanque

si=[zz(1)];
%EDO que describe tasa de crecimiento donde x es el peso o variable independiente   
aa=1;
ee=[];
%Metodo de Runge-Kutta para resolver la EDO

for jj=1:(length(xb)-1)
    rb=zz(jj);
    gg=(aa*rb);
    ee(1)=abs(yb(1)-ys(1));
    Fb_xy=@(t,yb) ((1-a).*b.*tau.*delta.*niu*betha*(gg).*yb.^m)-k*yb.^n;
    kb_1 = Fb_xy(xb(jj),yb(jj));
    kb_2 = Fb_xy(xb(jj)+0.5*h,yb(jj)+0.5*h*kb_1);
    kb_3 = Fb_xy((xb(jj)+0.5*h),(yb(jj)+0.5*h*kb_2));
    kb_4 = Fb_xy((xb(jj)+h),(yb(jj)+kb_3*h));
    yb(jj+1) = yb(jj) + (1/6)*(kb_1+2*kb_2+2*kb_3+kb_4)*h;
    si(end+1)=gg;
    ee(end+1) = abs(yb(jj+1)-ys(jj+1));
    if ee(jj) > 0.1*ys(jj)
        aa = 0.70;
        %e;
        %ii(end+1)=aa*rb
        %elseif abs(yb(i+1)-ys(i+1)) < 0.8*ys(i+1)
        %aa = 0.8;
        %yb(i+1)-yb(i)
        %0
        %ii(end+1)=aa*rb
    else
        aa = 1;
        %1
        %ii(end+1)=aa*rb
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xw = t0:h:tf;         
yw = zeros(1,length(xw)); 
yw(1) = 20;   

for i=1:(length(xw)-1)                     
    %1  2    3   4   5   6   7   8   9  10   11  12  13   14  15  16  17  18  19  20  21  22  23
    %0 0.9 2.1 2.1 2.0 1.9 1.7 1.6 1.5  1.5 1.4 1.4 1.3  1.5 1.5 1.5 1.5 1.4 1.2 1.2 1.2 1.2 1.2

    if i>0 & i<3
        rw=1.0;
    elseif i>=3 & i<6
        rw=1.8;
    elseif i>=6 & i<= 7
        rw=1.5;
    elseif i>7 & i<=18 
        rw=1.2;
    else
        rw=0.8;
    end
    
    F_xy=@(t,yw) ((1-a).*b.*tau.*delta.*niu*betha*(rw).*yw.^m)-k*yw.^n;

    kw_1 = F_xy(xw(i),yw(i));
    kw_2 = F_xy(xw(i)+0.5*h,yw(i)+0.5*h*kw_1);
    kw_3 = F_xy((xw(i)+0.5*h),(yw(i)+0.5*h*kw_2));
    kw_4 = F_xy((xw(i)+h),(yw(i)+kw_3*h));
    yw(i+1) = yw(i) + (1/6)*(kw_1+2*kw_2+2*kw_3+kw_4)*h;
end

%-------------- Experimental ------------------------
%PARÁMETROS DEL SISTEMA        
n=0.81; m=0.67;  b=0.62; a=0.53;   %se asume b de otro articulo
hache=0.8; kmin=0.00133; j=0.0132; s=21.38;
%CONDICIONES AMBIENTALES DEL ESTANQUE
Temp=33; OD=4; A=0.83;  
%CONDICIONES CRITICAS
Tcri=41; ODcrit=1; Acri=0.001;
%CONDICIONES MINIMAS
Tmin= 15; ODmin=0.3; Amax=1.4;
%DADO QUE EL ESTANQUE ES MUY GRANDE SE ASUME QUE LA PRESENCIA DE LA TILAPIA
%NO AFECTA AL MEDIO ACUATICO
tau=0.9948; delta=0.9429; niu=.4254; betha=0.7750; 
%ECUACIONES ALGEBRAICAS
k=kmin.*exp(j.*(Temp-Tmin));

%CONDICIONES DE OPERACION 
t0=0; %condicion inicial de tiempo (eje x)
tf=157; %tiempo final
h=1;             % step size
xq = t0:h:tf;         
yq = zeros(1,length(xq)); 
yq(1) = 20;   

%Dosificación de alimento (%peso de alimento respecto al peso del animal) con base al tiempo en el interior del estanque

%RESTRICCIONES DE ALIMENTACIÓN DE TILAPIA EN FUNCION DEL PESO DE LA TILAPIA

%Metodo de Runge-Kutta para resolver la EDO
for i=1:(length(xq)-1)
    if xq>=0 & xq<=21
        rq=0.84;
    elseif xq>21 & xq<=49
        rq=0.9;
    elseif xq>49 & xq<=70 
        rq=1.03;
    elseif xq>71 & xq<=91
        rq=0.92;
    elseif xq>92 & xq<=105
        rq=1.15;
    elseif xq>105 & xq<=119
        rq=1.25;
    elseif xq>119 & xq<=133
        rq=1.26;
    elseif xq>133 & xq<=140
        rq=1.28;
    else
        rq=1.30;
    end
    
    %EDO que describe tasa de crecimiento donde x es el peso o variable independiente   
    Fq_xy=@(t,yq) ((1-a).*b.*tau.*delta.*niu*betha*(rq).*yq.^m)-k*yq.^n;
    
    kq_1 = Fq_xy(xq(i),yq(i));
    kq_2 = Fq_xy(xq(i)+0.5*h,yq(i)+0.5*h*kq_1);
    kq_3 = Fq_xy((xq(i)+0.5*h),(yq(i)+0.5*h*kq_2));
    kq_4 = Fq_xy((xq(i)+h),(yq(i)+kq_3*h));
    yq(i+1) = yq(i) + (1/6)*(kq_1+2*kq_2+2*kq_3+kq_4)*h; 
end

%-------------- Levenberg-Marquardt --------------------
%PARÁMETROS DEL SISTEMA        
n=0.81; m=0.67;  b=0.62; a=0.53;   %se asume b de otro articulo
hache=0.8; kmin=0.00133; j=0.0132; s=21.38;
%CONDICIONES AMBIENTALES DEL ESTANQUE
Temp=33; OD=4; A=0.83;  
%CONDICIONES CRITICAS
Tcri=41; ODcrit=1; Acri=0.001;
%CONDICIONES MINIMAS
Tmin= 15; ODmin=0.3; Amax=1.4;
%DADO QUE EL ESTANQUE ES MUY GRANDE SE ASUME QUE LA PRESENCIA DE LA TILAPIA
%NO AFECTA AL MEDIO ACUATICO
tau=0.9948; delta=0.9429; niu=.4254; betha=0.7750; 
%ECUACIONES ALGEBRAICAS
k=kmin.*exp(j.*(Temp-Tmin));

%CONDICIONES DE OPERACION 
t0=0; %condicion inicial de tiempo (eje x)
tf=157; %tiempo final
h=1;             % step size
xa = t0:h:tf;         
ya = zeros(1,length(xq)); 
ya(1) = 20;

%Dosificación de alimento (%peso de alimento respecto al peso del animal) con base al tiempo en el interior del estanque
%RESTRICCIONES DE ALIMENTACIÓN DE TILAPIA EN FUNCION DEL PESO DE LA TILAPIA

%Metodo de Runge-Kutta para resolver la EDO
for i=1:(length(xa)-1)
    
    if i>=0 & i<=7
        ra=0.1; %0.1
    elseif i>7 & i<=14
        ra=3.0; %ra=2.4;
    elseif i>14 & i<=21 
        ra=2.5; %ra=2.0
    elseif i>21 & i<=28
        ra=2.1; %ra=1.68
    elseif i>28 & i<=35
        ra=1.90; %ra=1.52
    elseif i>35 & i<=42
        ra=1.80; %ra=1.44
    elseif i>42 & i<=49
        ra=1.70; %ra=1.36
    elseif i>49 & i<=56
        ra=1.60;  %ra=1.28
    elseif i>56 & i<=63
        ra=1.50;  %ra=1.20
    elseif i>63 & i<=77
        ra=1.40;  %ra=1.12
    elseif i>77 & i<=91
        ra=1.30; %ra=1.12
    else
        ra=1.20;  %ra=1.04
    end
    
    %EDO que describe tasa de crecimiento donde x es el peso o variable independiente   
    Fa_xy=@(t,ya) ((1-a).*b.*tau.*delta.*niu*betha*(ra).*ya.^m)-k*ya.^n;
    
    ka_1 = Fa_xy(xa(i),ya(i));
    ka_2 = Fa_xy(xa(i)+0.5*h,ya(i)+0.5*h*ka_1);
    ka_3 = Fa_xy((xa(i)+0.5*h),(ya(i)+0.5*h*ka_2));
    ka_4 = Fa_xy((xa(i)+h),(ya(i)+ka_3*h));
    ya(i+1) = ya(i) + (1/6)*(ka_1+2*ka_2+2*ka_3+ka_4)*h;
end

%Comienza bloque de graficacion de las diferentes curvas de peso
hold on
for i=1:ew    
    plot(x,Wi(i,:),'g') % curvas verdes calculadas con valores fijos de "r"
end

plot(x,ys,'r','linewidth', 4)  % curva calculada con el polinomio de ajuste
%plot(xw,yw,'dk','linewidth',4)
plot(xq,yq,'o')    % curva calculada con los valores de r que encontraron en los artículos
plot(xs,yd,'*r')   % curva calculada con datos de san buenaventura
plot(xb,yb,'^b')   % curva calculada con el r variable (metodo Doc Santiago)
%scatter(xd,ym,'+k','linewidth',8)  % datos de peso medidos por mena
plot(xa,ya,'+k') % curva calculada con los r de lebenger -macuard
xlabel('Tiempo (días)');
ylabel('Crecimiento (gr)');
plot(xq,rq,'+g')

%Bloque en el que se guradan datos de salida en archivos .txt
for i=1:ew
    for j=1:aw
        %filename = 'PesosVar.txt';
        filename=sprintf('ModeloRK.txt');
        fid = fopen(filename, "a");
        fprintf(fid, "%f ",[xq(j), Wi(i,j)]);
        fprintf(fid, "\n");
        fclose(fid);
    end
end

for j=1:length(xd)
    %filename = 'PesosVar.txt';
    filename=sprintf('DatosTilapia.txt');
    fid = fopen(filename, "a");
    fprintf(fid, "%f ",[xd(j), yd(j)]);
    fprintf(fid, "\n");
    fclose(fid);
end

for j=1:7:length(xs)
    %filename = 'PesosVar.txt';
    filename=sprintf('ModeloSBuenaventura.txt');
    fid = fopen(filename, "a");
    fprintf(fid, "%f ",[xs(j), yb(j)]);
    fprintf(fid, "\n");
    fclose(fid);
end

for j=1:7:length(xb)
    %filename = 'PesosVar.txt';
    filename=sprintf('ModeloDSantiago.txt');
    fid = fopen(filename, "a");
    fprintf(fid, "%f ",[xb(j), yb(j)]);
    fprintf(fid, "\n");
    fclose(fid);
end

for j=1:7:length(ys)
    filename=sprintf('DatosTilapiaAjusteO3.txt');
    fid = fopen(filename, "a");
    fprintf(fid, "%f ",[x(j), ys(j)]);
    fprintf(fid, "\n");
    fclose(fid);
end

for j=1:7:length(ys)
    filename=sprintf('DatosTilapiaLM.txt');
    fid = fopen(filename, "a");
    fprintf(fid, "%f ",[xq(j), yq(j)]);
    fprintf(fid, "\n");
    fclose(fid);
end