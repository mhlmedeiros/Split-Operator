#*******************************************************************************
#
#
#                Aplicação do split-operator em um sistema 1D
#
#
#*******************************************************************************
using PyPlot

# constantes físicas
Ry = 13.6056917253e3; # meV
a0 = 0.0529177208319; # nm
h2m = (Ry*a0^2); # /0.067; massa efetiva GaAs
hbar = 0.65821188926; # meV.ps


# definição da rede espacial
L = 1000;
Nx = 2*L;
x = linspace(-L/2, L/2, Nx+1);
x = x[1:(end-1)];
Δx = x[2] - x[1];

# definição da rede do "espaço-recíproco"
k = fftshift(linspace(-pi/Δx, pi/Δx, Nx+1)[1:end-1]);
Δk = k[2]-k[1];

# cond. inicial
σ = 15 ; # L/100
k0_max = pi/(2*Δx);
k₀ = 0.1 * k0_max;  # o valor de "k0" não pode ultrapassar π/(2dx)
x_0 = -200;
Ψ_0 = exp(-(x-x_0).^2/(2*σ^2) + 1im*k₀*x); # Gaussian packet

# Hamiltoniano separado em "T + V"
α, β, γ = 10, 0.07, 0;
T_k = h2m * k.^2;
V_x = α./(e.^(-β*(x-γ))+1);


function psi_dt(ψ,Δt)
  # espaço de posições 'x'
  psi_spinor = ψ;
  psi_spinor = [exp(-1im*V_x[i]*Δt/(2*hbar))*psi_spinor[i] for i=1:Nx ];
  # espaço de momentos 'k'
  psi_spinor = fft(psi_spinor);
  psi_spinor = [ exp(-1im*T_k[i]*Δt/hbar)*psi_spinor[i] for i=1:Nx ];
  # espaço de posições 'x'
  psi_spinor = ifft(psi_spinor);
  psi_spinor = [exp(-1im*V_x[i]*Δt/(2*hbar))*psi_spinor[i] for i=1:Nx ];
  # global Ψ_0 = psi_spinor;
  return psi_spinor;
end


function psit(Ψ, t_final, N_steps = 20)
  tspan = linspace(0, t_final, N_steps);
  delta_t = tspan[2]-tspan[1];

  for i = 1:N_steps
    psi_t = psi_dt(Ψ,delta_t);
    Ψ = psi_t;
  end
  Ψ
end


ψ_1 = psit(Ψ_0, 10.0,100)
plot(x,abs(ψ_1).^2)
plot(x,10.0^(-1)*V_x, linewidth=2.0);
# ylim(-0.05,1.05);
# xlim(-20,20)
title("Gaussian distribution");
xlabel(L"$x$");
ylabel(L"$|\psi(x)|^2$");
