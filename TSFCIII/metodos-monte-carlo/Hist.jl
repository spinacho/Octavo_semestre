function histograma(X::Array,divisiones::Int)
    #esta función toma como entrada un arreglo y un número de divisiones, y a la salida escupe el número de elementos que
    #había en cacda uno de los intervalos (las divisiones son homogéneas)
    #además escupe un arreglo que tiene un elemento en la mitad de cada división para poder graficar

    Xmin=minimum(X)
    Xmax=maximum(X)
    intervalo=Xmax-Xmin
    INTERVALOS=[Xmin:intervalo/divisiones:Xmax]
    contadora=zeros(Int,divisiones)

    for i in 1:divisiones-1
        contadora[i]=count(x-> INTERVALOS[i]<=x<INTERVALOS[i+1],X)
    end

    contadora[end]=count(x-> INTERVALOS[end-1]<=x<=INTERVALOS[end],X)
    Plotting_itvls=zeros(divisiones)

    for i in 1:divisiones
        Plotting_itvls[i]=(INTERVALOS[i+1]+INTERVALOS[i])/2
    end

    return contadora,Plotting_itvls
end

function normal_histo(X::Array,divisiones::Int)
    #hace lo mismo que el histograma pero de forma normalizada

    H,INTERVALOS=histograma(X,divisiones)
    TOTAL=sum(H)
    Hnorm=H/TOTAL
    return Hnorm,INTERVALOS
end;
