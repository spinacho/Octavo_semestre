using PyPlot

module XY_2D

  using PyPlot

  type Spin_Configuration
    L::Int
    T::Float64
    Mat::Array{Float64,2}
    β::Float64
    J::Float64
    E::Float64
    Mag::Complex64
  end

  type Changeor
    i::Int
    j::Int
    θ::Float64
    ΔE::Float64
    α::Bool
    θ_old::Float64
  end

  function Changeor(L::Int)
    i,j=spin_choose(L)
    θ=2*π*rand()
    return Changeor(i,j,θ,0,false,θ)
  end

  function periodic_energy(S::Spin_Configuration)
      E=0

      for i in 1:S.L, j in 1:S.L
          E+=cos(S.Mat[i,j]-S.Mat[mod1(i-1,S.L),j])+cos(S.Mat[i,j]-S.Mat[mod1(i+1,S.L),j])+cos(S.Mat[i,j]-S.Mat[i,mod1(j-1,S.L)])+cos(S.Mat[i,j]-S.Mat[i,mod1(j+1,S.L)])
      end

      return S.J*int(-E/2)
  end

  function magnetization(S::Spin_Configuration)
    return sum(cos(S.Mat))+im*(sum(sin(S.Mat)))
  end

  function config_maker(L::Int,T::Real,J=1)
    SM=Array(Float64,(L,L))

    for i in 1:L, j in 1:L
            SM[i,j]=2*π*rand()
    end

    S=Spin_Configuration(L,T,SM,inv(T),J,0,0)
    S.E=periodic_energy(S)
    S.Mag=magnetization(S)
    return S
  end

  function spin_choose(L::Int)
      i=int(floor(L*rand()+1))
      j=int(floor(L*rand()+1))

      return i,j
  end

  function spin_flip!(S::Spin_Configuration,C::Changeor)
      S.Mat[C.i,C.j]=C.θ
  end

  function ΔE!(S::Spin_Configuration,C::Changeor)
    E_i=cos(S.Mat[C.i,C.j]-S.Mat[mod1(C.i-1,S.L),C.j])+cos(S.Mat[C.i,C.j]-S.Mat[mod1(C.i+1,S.L),C.j])+cos(S.Mat[C.i,C.j]-S.Mat[C.i,mod1(C.j-1,S.L)])+cos(S.Mat[C.i,C.j]-S.Mat[C.i,mod1(C.j+1,S.L)])
    E_f=cos(C.θ-S.Mat[mod1(C.i-1,S.L),C.j])+cos(C.θ-S.Mat[mod1(C.i+1,S.L),C.j])+cos(C.θ-S.Mat[C.i,mod1(C.j-1,S.L)])+cos(C.θ-S.Mat[C.i,mod1(C.j+1,S.L)])
    C.ΔE=S.J*(E_i-E_f)
  end

  function acceptance!(S::Spin_Configuration,C::Changeor)
    ΔE!(S,C)

    if rand()<exp(-(S.β)*(C.ΔE))
        C.α=true
        C.θ_old=S.Mat[C.i,C.j]
    end
  end

  function new_energy!(S::Spin_Configuration,C::Changeor)
     S.E+=C.ΔE
  end

  function new_mag!(S::Spin_Configuration,C::Changeor)
    S.Mag+=-cos(C.θ_old)+cos(C.θ)+im*(-sin(C.θ_old)+sin(C.θ))
  end

  function one_step_flip!(S::Spin_Configuration,C::Changeor)
    acceptance!(S,C)

      if C.α==true
          spin_flip!(S,C)
      end
  end

  function config_runner(L::Int,steps::Int,T)
      Accepted=0
      S=config_maker(L,T)
      Energy=[S.E]
      sizehint(Energy,steps)
      Mags=[S.Mag]
      sizehint(Mags,steps)

      for i in 1:steps
          C=Changeor(L)
          one_step_flip!(S,C)
          Accepted+=C.α
            if C.α==true
              push!(Energy,Energy[i]+C.ΔE)
            else
              push!(Energy,Energy[i])
            end
          new_mag!(S,C)
          push!(Mags,S.Mag)
      end

      return S.Mat,Energy,Mags,Accepted
  end

  function plotter(L::Int,steps::Int,T,IMPLOT=true,QUIPLOT=true,EPLOT=true,MAGPLOT=true)
    S,E,M,Accepted=config_runner(L,steps,T)

    if IMPLOT==true
      figure()
      imshow(cos(S))
      colorbar()
    end

    if QUIPLOT==true
      figure()
      x=[1:L]
      y=[L:-1:1]
      quiver(x,y,cos(S),sin(S),scale=L)
    end

    if EPLOT==true
      figure()
      plot(E)
      title("Energía a temperatura $T con respecto al tiempo")
      xlabel("Tiempo")
      ylabel("Energía")
    end

    if MAGPLOT==true
      figure()
      plot(abs(M))
      title("Magnetización a temperatura $T con respecto al tiempo")
      xlabel("Tiempo")
      ylabel("Magnetización")
    end

  end

  function mean_runner(L::Int,steps::Int,T,BC=0.85)
    #BC (balance constant) es la constante que vamos a utilizar para suponer que el sistema ya está en estado de equilibrio
    #por lo general vamos a utilizarla como 0.85, sea que sólo tomamos en cuenta el último 15% de la simulación
    ET=Float64[]
    MT=Float64[]
    AT=Float64[]
    sizehint(ET,length(T))
    sizehint(MT,length(T))
    sizehint(AT,length(T))

    for temp in T
      M,E,T,Accepted=config_runner(L,steps,temp)
      push!(ET,mean(E[int(BC*end):end]))
      push!(MT,mean(M[int(BC*end):end]))
      push!(AT,Accepted)
    end

    AT=AT/steps

    return ET,MT,AT
  end

  function Mean_Plotter(L::Int,steps::Int,T,EPLOT=true,MAGPLOT=true,APLOT=true,BC=0.85)
    ET,MT,AT=mean_runner(L,steps,T,BC)

    if EPLOT==true
        figure()
        plot(T,ET)
        title("Energía promedio con respecto a la temperatura")
        xlabel("Temperatura")
        ylabel("Energía promedio")
      end

      if MAGPLOT==true
        figure()
        plot(T,abs(MT)/(L^2))
        title("Magnetización promedio con respecto a la temperatura")
        xlabel("Temperatura")
        ylabel("Magnetización promedio")
      end

      if APLOT==true
        figure()
        plot(T,100*AT)
        title("Promedio de aceptación con respecto a la temperatura")
        xlabel("Temperatura")
        ylabel("Porcentaje de aceptados")
      end
  end

  function config_teller_for_anim(L::Int,steps::Int,jump::Int,T)
    S=config_maker(L,T)
    MATS=Array[]
    sizehint(MATS,int(steps/jump))
      for i in 1:steps
          C=Changeor(L)
          one_step_flip!(S,C)
            if i%jump==0
              push!(MATS,S.Mat)
            end
      end

      return MATS
  end
end

import XY_2D.Mean_Plotter
import XY_2D.config_runner
import XY_2D.plotter
