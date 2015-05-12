using PyPlot

module XY

  using PyPlot

  type Spin_Configuration
    L::Int
    N::Int
    T::Float64
    Mat::Array{Float64,2}
    β::Float64
    J::Float64
    E::Float64
    Mag::Complex64
    Abs_Mag::Float64
  end

  function periodic_energy(S::Array,L::Int,J)

      E=0

      for i in 1:L, j in 1:L
          E+=cos(S[i,j]-S[mod1(i-1,L),j])+cos(S[i,j]-S[mod1(i+1,L),j])+cos(S[i,j]-S[i,mod1(j-1,L)])+cos(S[i,j]-S[i,mod1(j+1,L)])
      end

      return J*int(-E/2)

  end;

  function periodic_energy(S::Spin_Configuration)

      E=0

      for i in 1:S.L, j in 1:S.L
          E+=cos(S.Mat[i,j]-S.Mat[mod1(i-1,S.L),j])+cos(S.Mat[i,j]-S.Mat[mod1(i+1,S.L),j])+cos(S.Mat[i,j]-S.Mat[i,mod1(j-1,S.L)])+cos(S.Mat[i,j]-S.Mat[i,mod1(j+1,S.L)])
      end

      return S.J*int(-E/2)

  end;

  function magnetization(S::Array)
    return sum(cos(S))+im*(sum(sin(S)))
  end

  function magnetization(S::Spin_Configuration)
    return sum(cos(S.Mat))+im*(sum(sin(S.Mat)))
  end

  function Spin_Configuration(L::Int,T::Real,Mat::Array,J)

    Spin_Configuration(L,L*L,T,Mat,inv(T),J,periodic_energy(Mat,L,J),magnetization(Mat),abs(magnetization(Mat)))

  end


  function config_maker(L::Int,T::Real,J=1)
    SM=Array(Float64,(L,L))

    for i in 1:L, j in 1:L
            SM[i,j]=2*π*rand()
    end

    return Spin_Configuration(L,T,SM,J)
  end

  function spin_choose(L::Int)

      i=int(floor(L*rand()+1))
      j=int(floor(L*rand()+1))

      return i,j
  end

  type Changeor
    i::Int
    j::Int
    θ::Float64
    ΔE::Float64
    α::Bool
    θ_old::Float64
  end

  function spin_flip!(S::Spin_Configuration,C::Changeor)

      S.Mat[C.i,C.j]=C.θ

  end

  function Changeor(L::Int)
    i,j=spin_choose(L)
    θ=2*π*rand()
    return Changeor(i,j,θ,0,false,θ)
  end

  function ΔE!(S::Spin_Configuration,C::Changeor)
    E_i=cos(S.Mat[C.i,C.j]-S.Mat[mod1(C.i-1,S.L),C.j])+cos(S.Mat[C.i,C.j]-S.Mat[mod1(C.i+1,S.L),C.j])+cos(S.Mat[C.i,C.j]-S.Mat[C.i,mod1(C.j-1,S.L)])+cos(S.Mat[C.i,C.j]-S.Mat[C.i,mod1(C.j+1,S.L)])
    E_f=cos(C.θ-S.Mat[mod1(C.i-1,S.L),C.j])+cos(C.θ-S.Mat[mod1(C.i+1,S.L),C.j])+cos(C.θ-S.Mat[C.i,mod1(C.j-1,S.L)])+cos(C.θ-S.Mat[C.i,mod1(C.j+1,S.L)])
    C.ΔE=E_f-E_i
  #println("la delta de energía es $(C.ΔE)")
  end

  #function ΔE!(S::Spin_Configuration,C::Changeor)
  #    S_new=copy(S.Mat)
  #    S.Mat[C.i,C.j]=C.θ
  #    C.ΔE = periodic_energy(S_new,S.L,S.J)-periodic_energy(S)
  #end

  function acceptance!(S::Spin_Configuration,C::Changeor)
    ΔE!(S,C)
  #println("la prob de aceptación es $(exp(-(S.β)*(C.ΔE)))")

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
    S.Abs_Mag=abs(S.Mag)
  end

  function one_step_flip!(S::Spin_Configuration,C::Changeor)
    acceptance!(S,C)

      if C.α==true
          spin_flip!(S,C)
      end

  end;

  function montecarlo_config_run(L::Int,steps::Int,T)
      Acc=0
      S=config_maker(L,T)
      for i in 1:steps
          C=Changeor(L)
          one_step_flip!(S,C)
          Acc+=C.α
      end

      return S.Mat,Acc
  end;


  function montecarlo_energy_run(L::Int,steps::Int,T)
      S=config_maker(L,T)
      Energy=[S.E]
      sizehint(Energy,steps)
      for i in 1:steps
          C=Changeor(S.L)
          one_step_flip!(S,C)
          if C.α==true
            push!(Energy,Energy[i]+C.ΔE)
          else
            push!(Energy,Energy[i])
          end
      end

      return Energy
  end

  function montecarlo_mag_run(L::Int,steps::Int,T)
      S=config_maker(L,T)
      Mags=[S.Abs_Mag]
      sizehint(Mags,steps)
      for i in 1:steps
          C=Changeor(S.L)
          one_step_flip!(S,C)
          new_mag!(S,C)
          push!(Mags,S.Abs_Mag)
      end

      return Mags
  end

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
    end

    title("Valores absolutos de la magnetización a temperatura $T con respecto al tiempo")
    xlabel("Tiempo")
    ylabel("Valor absoluto de la magnetización")
  end

  function Mean_Energies_T(L::Int,steps::Int,Temp::Array)
    ET=Float64[]
    T=Temp
    sizehint(ET,length(T))

    for temp in T
        E=montecarlo_energy_run(L,steps,temp)
        push!(ET,mean(E[int(0.85*end):end]))
    end

    plot(T,ET)
    xlabel("Temperatura")
    ylabel(L"Energia promedio $<E>_T$")
    title("Energía promedio con respecto a la temperatura");
  end

  function Mean_Energies_β(L::Int,steps::Int,Temp::Array)
    ET=Float64[]
    T=Temp
    sizehint(ET,length(T))
    @time for temp in T
        E=montecarlo_energy_run(L,steps,temp)
        push!(ET,mean(E[int(0.85*end):end]))
    end
    plot([1/t for t in T],ET)
    xlabel(L"$\beta$")
    ylabel(L"Energia promedio $<E>_\beta$")
    title(L"Energia promedio con respecto a $\beta$");
  end

  function Mean_Magnetizations_T(L::Int,steps::Int,Temp::Array)
    MT=Float64[]
    T=Temp
    sizehint(MT,length(T))
    for temp in T
        M=montecarlo_mag_run(L,steps,temp)
        #push!(MT,mean(M[int(0.85*end):end]))
        push!(MT,mean(M))
    end
    plot(T,MT,"o-")
    xlabel("Temperatura")
    ylabel(L"Magnetizacion promedio $<M>_T$")
    title("Magnetización promedio con respecto a la temperatura");
  end

  function Mean_Abs_Magnetizations(L::Int,steps::Int,Temp::Array)
    MT=Complex[]
    T=Temp
    sizehint(MT,length(T))
    for temp in T
        M=montecarlo_mag_run(L,steps,temp)
    #push!(MT,mean(abs(M[int(0.85*end):end])))
    push!(MT,mean((M[int(0.85*end):end])))
    end
    plot(T,abs(MT),"o-")
    xlabel("Temperatura")
    ylabel(L"Valor absoluto de la magnetizacion promedio $|<M>_T|$")
    title("Valor absoluto de la magnetización promedio con respecto a la temperatura \n");
  end

  function Mean_Abs_Magnetizations_β(L::Int,steps::Int,Temp::Array)
    MT=Float64[]
    T=Temp
    sizehint(MT,length(T))
    for temp in T
        M=montecarlo_mag_run(L,steps,temp)
    push!(MT,mean(abs(M[int(0.85*end):end])))
    end
    plot([1/t for t in T],MT)
    xlabel(L"$\beta$")
    ylabel(L"Valor absoluto de la magnetizacion promedio $|<M>_β|$")
    title("Valor absoluto de la magnetización promedio con respecto a la β \n");
  end

  function Configs(L::Int,steps::Int,Temp::Array)
    @time for T in Temp
        figure()
        S=montecarlo_config_run(L,steps,T)
        imshow(cos(S))
    end
  end

end

