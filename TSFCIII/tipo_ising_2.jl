#Eduardo Ibarra García Padilla

#ELIMINAMOS EL PROBABILITY CHANGE DEL TIPO ISING Y VEMOS SI ESTO AYUDA A QUE CORRA MÁS RÁPIDO

using PyPlot

"""
Conjunto de tipos para trabajar con el algoritmo
de Metropolis Ising, específicamente el modelo
de Ising. Se trata de hacer lo más general los
tipos para que éstos puedan ser utilizados después
bajo otras circunstancias.
"""

#TIPO CONFIGURACIÓN
"""
Este tipo permite trabajar con un arreglo,
cuyas dimensiones vienen especificadas en el vector
L, el arreglo viene dado por M (que es sobre el cual
se opera). Asimismo trae la energía de la
configuración y su magnetización.
"""

type Configuracion
  L::(Int,Int)
  M::Array{Int,2}
  E::Int
  Mag::Int
end


#TIPO METROPOLIS (ISING)
"""
Este tipo trae las ideas esenciales del método
de Metroplis-Hastings, para cualquier configuración.
La idea es que uno pueda suministrar cualquier
configuración y nuestro código opere correctamente.
"""

type Metropolis_Ising
  spin_choose::Function
  spin_flipper!::Function
  acceptance::Function
end

function range_rand(l::Int)
    int(floor(l*rand()))+1
end

function spin_choose_2D(S::Configuracion)
  i=range_rand(S.L[1])
  j=range_rand(S.L[2])
  return i,j
end

function spin_flipper_2D!(S::Configuracion,i,j)
  S.M[i,j]=-S.M[i,j]
end

function acceptance(ΔE::Function,T,S::Configuracion,i,j)
  deltaE=ΔE(S,i,j)
  α=min(exp(-(1/T)*(deltaE)),1)

  if rand()<α
    true,deltaE
  else
    false,0
  end
end


#TIPO ISING 2D
"""
Este tipo trae las ideas esenciales del problema
del modelo de ising, como cómo calcular la energía,
la magnetización del sistema, sus ΔE, su
distribución de probabilidad de cambio, etc...
"""

type Ising_model
  ΔE::Function
  energia_nueva!::Function
  magnetizacion_nueva!::Function
end

function ΔE_2D_ising(S::Configuracion,i,j)
  DeltaE = 2*S.M[i,j]*(S.M[mod1(i-1,S.L[1]),j]+S.M[mod1(i+1,S.L[1]),j]+S.M[i,mod1(j-1,S.L[2])]+S.M[i,mod1(j+1,S.L[2])])
end

function energia_periodica_2D(S::Configuracion)
  E=0
  for i in 1:S.L[1], j in 1:S.L[2]
    E += S.M[i,j]*(S.M[mod1(i-1,S.L[1]),j]+S.M[mod1(i+1,S.L[1]),j]+S.M[i,mod1(j-1,S.L[2])]+S.M[i,mod1(j+1,S.L[2])])
  end
    int(-E/2)
end

function magnetizacion(S::Configuracion)
    sum(S.M)
end

function energia_nueva!(S::Configuracion,deltaE)
  S.E += deltaE
end

function magnetizacion_nueva!(S::Configuracion,aceptor,i,j)
  if aceptor
    if S.M[i,j]==1
      S.Mag += 2
    else
      S.Mag -= 2
    end
  end
end


#CONSTRUCTORES DE TIPOS
"""
Crear los tipos de Configuracion, Metropolis_Ising
e Ising_model. En este caso en 2D.
"""

function spin()
  if rand()<0.5
    -1
  else
    1
  end
end

function config_maker_2D(L)
  M=zeros(Int,L[1],L[2])
  for i in 1:length(M)
    M[i]=spin()
  end
  S=Configuracion(L,M,0,0)
  S.E=energia_periodica_2D(S)
  S.Mag=magnetizacion(S)
  S
end

function metropolis_ising_maker_2D()
  Metropolis=Metropolis_Ising(spin_choose_2D,spin_flipper_2D!,acceptance)
end

function ising_model_maker_2D()
  Ising=Ising_model(ΔE_2D_ising,energia_nueva!,magnetizacion_nueva!)
end


#FUNCIONES PARA HACER REPETICIONES Y ANALIZAR LA FÍSICA DEL PROBLEMA

function one_step_flip(S::Configuracion,Metropolis::Metropolis_Ising,Ising::Ising_model,T)
  #Se basa en tomar un spin al azar y determinar
  #después si se cambia o no a esa nueva distribución
  #Regresa el valor del cambio de la energía,
  #pues la función spin flipper cambia ya el sistema.

  i,j=Metropolis.spin_choose(S)
  aceptacion,deltaE=Metropolis.acceptance(Ising.ΔE,T,S,i,j)

  if aceptacion
    Metropolis.spin_flipper!(S,i,j)
    deltaE,true,i,j
  else
    0,false,i,j
  end
end



function corrida_2D(Metropolis::Metropolis_Ising,Ising::Ising_model,L,T,N)
  #Construye una configuración de spines
  #Calcula las energías de las configuraciones
  #para cada paso, de un total de N, a β dada.

  S=config_maker_2D(L)

  Energias=Int[S.E]
  Magnetizaciones=Int[S.Mag]

  sizehint(Energias,N)
  sizehint(Magnetizaciones,N)

  for n in 1:N
    deltaE,aceptor,i,j=one_step_flip(S,Metropolis,Ising,T)
    Ising.energia_nueva!(S,deltaE)
    Ising.magnetizacion_nueva!(S,aceptor,i,j)
    push!(Energias,S.E)
    push!(Magnetizaciones,S.Mag)
  end
  Energias,Magnetizaciones
end



function graficar_E_M_vs_t(Metropolis::Metropolis_Ising,Ising::Ising_model,L,T,N)
  #Permite graficar Energía y Magnetización vs tiempo.
  Energias,Magnetizaciones=corrida_2D(Metropolis,Ising,L,T,N)
  tiempos=collect(0:N)

  #Gráfica de energía vs tiempo.
  figure(figsize=(6,4))
  plot(tiempos,Energias)
  title(string(L"$E(t)$"," a",L" $T= $"," $T"))
  xlabel(L"t")
  ylabel(L"E(t)")
  axis([0,N,minimum(Energias)-100,maximum(Energias)+100])

  #Gráfica de magnetización vs tiempo.
  figure(figsize=(6,4))
  plot(tiempos,Magnetizaciones)
  title(string(L"$M(t)$"," a",L" $T= $"," $T"))
  xlabel(L"t")
  ylabel(L"M(t)")
  axis([0,N,minimum(Magnetizaciones)-100,maximum(Magnetizaciones)+100])

end


#FUNCIONES PARA GRAFICAR EN FUNCIÓN DE LA TEMPERATURA

function E_M_promedio_vs_T(Metropolis::Metropolis_Ising,Ising::Ising_model,L,T_i,T_f,N,step=0.1)

  Temperaturas=[T for T in T_i:step:T_f]

  Energias=Float64[]
  Magnetizaciones=Float64[]

  sizehint(Energias,length(Temperaturas))
  sizehint(Magnetizaciones,length(Temperaturas))

  for T in Temperaturas
    Energy,Magnetization=corrida_2D(Metropolis,Ising,L,T,N)
    push!(Energias,mean(Energy[int(0.85*end):end]))
    push!(Magnetizaciones,mean(abs(Magnetization[int(0.85*end):end])))
  end
  Temperaturas,Energias,Magnetizaciones
end



function graficar_E_M_vs_T(Metropolis::Metropolis_Ising,Ising::Ising_model,L,T_i,T_f,N,step=0.1)

  Temperaturas,Energias,Magnetizaciones=E_M_promedio_vs_T(Metropolis,Ising,L,T_i,T_f,N,step)

  #Gráfica de la energía
  figure(figsize=(6,4))
  plot(Temperaturas,Energias,".-")
  title(L"$\langle E(T) \rangle$")
  xlabel(L"$T$")
  ylabel("Energía")

  #Gráfica de la magnetización
  figure(figsize=(6,4))
  plot(Temperaturas,Magnetizaciones,".-")
  title(L"$\langle \vert M(T) \vert  \rangle$")
  xlabel(L"$T$")
  ylabel("Magnetización")

end

"""
#Ahora vamos a correrlo
Metropolis=metropolis_ising_maker_2D()
Ising=ising_model_maker_2D()
graficar_E_M_vs_t(Metropolis,Ising,(2,2),0.1,10)
graficar_E_M_vs_T(Metropolis,Ising,(2,2),0.05,10,10,true,0.1)

graficar_E_M_vs_t(Metropolis,Ising,(20,20),0.1,10^4)
graficar_E_M_vs_T(Metropolis,Ising,(20,20),0.05,10,10^4,0.1)
"""