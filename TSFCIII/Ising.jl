using PyPlot

module Ising

#Creamos el tipo de variable Spin, el cual es un Int8 (sólo usaremos arriba=+1 y abajo=-1.

type Spin
  spin::Int8
end

function spin_flip!(s::Spin)
  s.spin=-s.spin
end

type Spin_Configuration
  L::Int
  N::Int
  T::Float64
  Mat::Array
  β::Float64
  J
  E
  Mag
end

function periodic_energy(S::Array,L::Int)

    E=0

    for i in 1:L, j in 1:L
        E+=S[i,j]*(S[mod1(i-1,L),j]+S[mod1(i+1,L),j]+S[i,mod1(j-1,L)]+S[i,mod1(j+1,L)])
    end

    return S.J*int(-E/2)

end;

function periodic_energy(S::Spin_Configuration)

    E=0

    for i in 1:S.L, j in 1:S.L
        E+=S.Mat[i,j]*(S.Mat[mod1(i-1,L),j]+S.Mat[mod1(i+1,L),j]+S.Mat[i,mod1(j-1,L)]+S.Mat[i,mod1(j+1,L)])
    end

    return S.J*int(-E/2)

end;

function magnetization(S::Array)
  return sum(S)
end

function magnetization(S::Spin_Configuration)
  return sum(S.Mat)
end

function Spin_Configuration(L::Int,T::Float64,Mat::Array,J)

  Spin_Configuration(L,L*L,T,Mat,inv(T),J,periodic_energy(Mat,L),magnetization(Mat))

end


function config_maker(L::Int,T::Float64,J=1)
  SM=Array(Spin,(L,L))

  for i in 1:L, j in 1:L
          if rand()<0.5
                SM[i,j].spin=-1
          else
              SM[i,j].spin=1
          end
      end

      return Spin_Configuration(L,T,SM,J)

end;

function spin_choose(L::Int)

    i=rand(1:L)
    j=rand(1:L)

    return i,j
end;

function spin_flip!(S::Spin_Configuration,i,j)

    spinflip!(S.Mat[i,j])

end;

type Changeor
  i::Int
  j::Int
  ΔE
  α::Bool
end

function Changeor(L::Int)
  i,j=spin_choose(L)
  return Changeor(i,j,0,false)
end


function ΔE(S,i,j)

    return S.J*2*S.Mat[i,j]*(S.Mat[mod1(i+1,S.L),j]+S.Mat[mod1(i-1,S.L),j]+S.Mat[i,mod1(j-1,S.L)]+S.Mat[i,mod1(j+1,S.L)])

end;

function ΔE!(S::Spin_Configuration,C::Changeor)

    C.ΔE = S.J*2*S.Mat[C.i,C.j]*(S.Mat[mod1(C.i+1,S.L),C.j]+S.Mat[mod1(C.i-1,S.L),C.j]+S.Mat[C.i,mod1(C.j-1,S.L)]+S.Mat[C.i,mod1(C.j+1,S.L)])

end;

function acceptance!(S::Spin_Configuration,C::Changeor)
  ΔE!(S,C)
  C.α=exp(-(S.β)*(C.ΔE))

  if rand()<α
      C.α=true
  else
      C.α=false
  end
end;

function new_energy!(S::Spin_Configuration,C::Changeor)
   S.E+=C.ΔE
end;

function new_mag!(S::Spin_Configuration,C::Changeor)
  if S.Mat[C.i,C.j]==1
    S.Mag+=-2
  else
    S.Mag+=+2
  end
end

function one_step_flip!(S::Spin_Configuration,energy==false,mag=false)
  C=Changeor(S.L)
  acceptance!(S,C)

    if C.α==true
        spin_flip!(S,C.i,C.j)
    end

  if energy
    new_energy!(S,C)
  end

  if mag
    new_mag!(S,C)
  end

end;

function montecarlo_config_run(L::Int,steps::Int,T)

    S=config_maker(L,T)

    for i in 1:steps
        one_step_flip!(S)
    end

    return S.Mat
end;


function montecarlo_energy_run(L::Int,steps::Int,T)
    S=config_maker(L,T)
    Energy=[S.E]
    sizehint(Energy,steps)

    for i in 1:steps
        one_step_flip!(S,true)
        push!(Energy,S.E)
    end

    return E
end;

function montecarlo_mag_run(L::Int,steps::Int,T)
    S=config_maker(L,T)
    Mags=[S.Mag]
    sizehint(Mags,steps)

    for i in 1:steps
        one_step_flip!(S,false,true)
        push!(Mags,S.Mag)
    end

    return Mag
end;

  function Energies(L::Int,steps::Int,Temp::Array)
    for T in Temp
      figure()
      E=montecarlo_energy_run(L,steps,T)
      plot(E)
      title("Energía a temperatura $T con respecto al tiempo")
      xlabel("Tiempo")
      ylabel("Energía")
    end
  end

  function Magnetizations(L::Int,steps,Temp::Array)
   for T in Temp
      figure()
      M=montecarlo_mag_run(L,steps,T)
      plot(abs(M))
      title("Valores absolutos de la magnetización a temperatura $T con respecto al tiempo")
      xlabel("Tiempo")
      ylabel("Valor absoluto de la magnetización")
    end
  end

  function Mean_Energies_T(L::Int,steps::Int,Temp::Array)
    ET=Float64[]
    T=Temp
    sizehint(ET,length(T))
    for temp in T
        E=montecarlo_energy_run(L,steps,temp)
        push!(ET,mean(E))
    end
    plot(T,ET)
    xlabel("Temperatura")
    ylabel(L"Energia promedio $<E>_{T}$")
  title("Energía promedio con respecto a la temperatura");
  end

  function Mean_Energies_β(L::Int,steps::Int,Temp::Array)
    ET=Float64[]
    T=Temp
    sizehint(ET,length(T))
    @time for temp in T
        E=montecarlo_energy_run(L,steps,temp)
        push!(ET,mean(E))
    end
    plot([1/t for t in T],ET)
    xlabel(L"$\beta$")
    ylabel(L"Energia promedio $<E>_{\beta}$")
    title(L"Energia promedio con respecto a $\beta$");
  end

  function Mean_Magnetizations_T(L::Int,steps::Int,Temp::Array)
    MT=Float64[]
    T=Temp
    sizehint(MT,length(T))
    @time for temp in T
        M=montecarlo_mag_run(L,steps,temp)
        push!(MT,mean(M))
    end
    plot(T,MT,"o-")
    xlabel("Temperatura")
    ylabel(L"Magnetizacion promedio $<M>_{T}$")
    title("Magnetización promedio con respecto a la temperatura");
  end

  function Mean_Abs_Magnetizations(L::Int,steps::Int,Temp::Array)
    MT=Float64[]
    T=Temp
    sizehint(MT,length(T))
    @time for temp in T
        M=montecarlo_mag_run(L,steps,temp)
        push!(MT,mean(M))
    end
    plot(T,abs(MT),"o-")
    xlabel("Temperatura")
    ylabel(L"Valor absoluto de la magnetizacion promedio $|<M>_{T}|$")
  title("Valor absoluto de la magnetización promedio con respecto a la temperatura \n");
  end

  function Configs(L::Int,steps::Int,Temp::Array)
    @time for T in Temp
        figure()
        S=montecarlo_config_run(L,steps,T)
        imshow(S)
    end
  end
end



