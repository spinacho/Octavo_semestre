using PyPlot

module Spin
  function config_maker(L::Int)

      S=zeros(L,L)

      for i in 1:L, j in 1:L
          if rand()<0.5
              S[i,j]=-1
          else
              S[i,j]=1
          end
      end

      return S
  end;

  function spin_choose(L)

      i=rand(1:L)
      j=rand(1:L)

      return i,j
  end;

  function spin_flip!(S::Array,i,j,N=length(S),L=int(sqrt(N)))

      S[i,j]*=-1

      return S
  end;

  function periodic_energy(S,N=length(S),L=int(sqrt(N)),J=1)

      E=0

      for i in 1:L, j in 1:L
          E+=S[i,j]*(S[mod1(i-1,L),j]+S[mod1(i+1,L),j]+S[i,mod1(j-1,L)]+S[i,mod1(j+1,L)])
      end

      return J*int(-E/2)

  end;

  function ΔE(S,i,j,N=length(S),L=int(sqrt(N)),J=1)

      return J*2*S[i,j]*(S[mod1(i+1,L),j]+S[mod1(i-1,L),j]+S[i,mod1(j-1,L)]+S[i,mod1(j+1,L)])

  end;

  function acceptance(S,T,i,j,N=length(S_old),L=int(sqrt(N)),J=1)

      Δe=ΔE(S,i,j,N,L,J)
      α=exp(-(1/T)*(Δe))

      if rand()<α
          return true,Δe
      else
          return false,0
      end
  end;

  function one_step_flip(S::Array,T,N=length(S),L=int(sqrt(N)),J=1)

      i,j=spin_choose(L)
      α,Δe=acceptance(S,T,i,j,N,L,J)
      if α==true
          spin_flip!(S,i,j,N,L)
          return Δe
      else
          return 0
      end
  end;

  function new_energy(E_old,Δe)

      return E_old+Δe

  end;

  function montecarlo_config_run(L::Int,steps::Int,T)
      N=L*L
      S=config_maker(L)

      for i in 1:steps
          one_step_flip(S,T,N,L,1)
      end

      return S
  end;

  for T in [0.5:0.5:3]
      figure()
      S=montecarlo_config_run(10,10^5,T)
      imshow(S)
  end

  function montecarlo_energy_run(L::Int,steps::Int,T)
      N=L*L
      S=config_maker(L)
      E=[periodic_energy(S)]
      sizehint(E,steps)

      for i in 1:steps
          Δe=one_step_flip(S,T,N,L,1)
          push!(E,new_energy(E[i],Δe))
      end

      return E
  end;

  magnetization(S)=sum(S);

  function montecarlo_mag_run(L::Int,steps::Int,T)
      N=L*L
      S=config_maker(L)
      M=[magnetization(S)]
      sizehint(M,steps)

      for i in 1:steps
          one_step_flip(S,T,N,L,1)
          push!(M,magnetization(S))
      end

      return M
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

import Spin.Mean_Abs_Magnetizations() as Mean_Abs_Magnetizations()
import Mean
