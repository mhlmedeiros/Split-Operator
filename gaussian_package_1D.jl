#*******************************************************************************
#
#
#                Aplicação do split-operator em um sistema 1D
#
#
#*******************************************************************************

include("save_the_data.jl")

# constantes físicas
Ry = 13.6056917253e3; # meV
a0 = 0.0529177208319; # nm
h2m = (Ry*a0^2)/0.067; # massa efetiva GaAs
hbar = 0.65821188926; # meV.ps


# definição da rede espacial
L = 1000;
Nx = 2*L;
x_grid = linspace(-L/2, L/2, Nx+1)
x_grid = x_grid[1:(end-1)];
Δx = x_grid[2] - x_grid[1];

# definição da rede do "espaço-recíproco"
k = fftshift(linspace(-pi/Δx, pi/Δx, Nx+1)[1:end-1]);
Δk = k[2]-k[1];

# Condição inicial
σ = 15 ;
k₀_max = pi/(2*Δx);
frac_max = 0.01; # fração do valor máximo de "k" que será usado na cond. inicial
k₀ = frac_max * k₀_max;  # o valor de "k₀" não pode ultrapassar π/(2*Δx)
x₀ = -200;
Ψ₀ = exp(-(x_grid-x₀).^2/(2*σ^2) + 1im*k₀*x_grid); # Gaussian packet


# Hamiltoniano separado em "T + V"
T_k = h2m * k.^2;

α, β, γ = 10, 0.07, 0;
V_x = α./(e.^(-β*(x_grid-γ))+1);


function psi_dt(ψ,Δt)
  # espaço de posições 'x'
  psi = ψ;
  psi = [exp(-1im*V_x[i]*Δt/(2*hbar))*psi[i] for i=1:Nx ];
  # espaço de momentos 'k'
  fft!(psi);
  psi = [ exp(-1im*T_k[i]*Δt/hbar)*psi[i] for i=1:Nx ];
  # espaço de posições 'x'
  ifft!(psi);
  psi = [exp(-1im*V_x[i]*Δt/(2*hbar))*psi[i] for i=1:Nx ];
  return psi;
end

function psit(Ψ; t_final = 10.0, N_steps = 100)
  tspan = linspace(0, t_final, N_steps);
  delta_t = tspan[2]-tspan[1];

  for i = 1:N_steps
    psi_t = psi_dt(Ψ,delta_t);
    Ψ = psi_t;
  end
  Ψ
end

t_f = 2.5;
ψ_1 = psit(Ψ_0, t_final = t_f, N_steps = 10);
nome_dados = string("psi_t_",round(t_f,2))

nome_arquivo = "dados_split_operator.h5";
nome_grupo = "split_1D";
nome_subg = "potencial_degrau";

atr_cond_ini = Dict("largura_pacote" => round(σ,1),
                        "x_0" => round(x₀,1),
                        "frac_k_max" => round(frac_max,3));

atr_pot = Dict("altura" => round(α,1),
                        "beta" => round(β,2),
                        "posição" => round(γ,1));


save_them_all( nome_arquivo, 
               nome_grupo,
               nome_subg,
               ψ_1,
               nome_dados,
               atr_cond_ini,
               atr_pot);
