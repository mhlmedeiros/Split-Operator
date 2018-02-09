#*******************************************************************************
#
#
#                Aplicação do split-operator em um sistema 1D
#
#
#*******************************************************************************
using PyPlot;
using HDF5;

PyPlot.matplotlib[:rc]("text", usetex=true) # allow tex rendering
PyPlot.matplotlib[:rc]("font", family="serif")

# constantes físicas
Ry = 13.6056917253e3; # meV
a0 = 0.0529177208319; # nm
h2m = (Ry*a0^2)/0.067; # massa efetiva GaAs
hbar = 0.65821188926; # meV.ps


# definição da rede espacial
L = 1000;
Nx = 2*L;
x = linspace(-L/2, L/2, Nx+1)
x = x[1:(end-1)];
Δx = x[2] - x[1];

# definição da rede do "espaço-recíproco"
k = fftshift(linspace(-pi/Δx, pi/Δx, Nx+1)[1:end-1]);
Δk = k[2]-k[1];

# cond. inicial
σ = 15 ; # L/100
k0_max = pi/(2*Δx);
k₀ = 0.01 * k0_max;  # o valor de "k0" não pode ultrapassar π/(2dx)
x₀ = -200;
Ψ_0 = exp(-(x-x_0).^2/(2*σ^2) + 1im*k₀*x); # Gaussian packet

# Hamiltoniano separado em "T + V"
α, β, γ = 10, 0.07, 0;
T_k = h2m * k.^2;
V_x = α./(e.^(-β*(x-γ))+1);


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


t_f = 1.0;
ψ_1 = psit(Ψ_0, t_final = t_f, N_steps = 10)
# plot(x,abs(ψ_1).^2)
# plot(x,0.5*10.0^(-2)*V_x, linewidth=2.0);
# # ylim(-0.05,1.05);
# # xlim(-20,20)
# title("Gaussian distribution");
# xlabel(L"$x$");
# ylabel(L"$|\psi(x)|^2$");





nome_dados_real = string("psi_t_",round(t_f,2),"_real")
nome_dados_imag = string("psi_t_",round(t_f,2),"_imag")

h5open("dados_split_operator.h5","w") do file
    # Cria os grupos, ou abre aqueles já existentes:
    try
        g = g_create(file, "split_1D");
        h = g_create(g, "potencial_degrau");
    catch
        g = file["split_1D"];
        h = g["potencial_degrau"];
    end

    # Caso não tenha salvo o potencial, salve-o:
    try
        h["Potencial"] = V_x;
        attrs(h["Potencial"])["alpha"] = string(α);
        attrs(h["Potencial"])["beta"] = string(β);
        attrs(h["Potencial"])["gamma"] = string(γ);
    end

    # Caso não tenha salvo a condição inicial, salve-a:
    try
        h["Cond_Ini_real"] = real(Ψ_0);
        h["Cond_Ini_imag"] = imag(Ψ_0);
        h["grid"] = x;
        attrs(g["potencial_degrau"])["x_0"] = string(round(x₀,1));
        attrs(g["potencial_degrau"])["k_0"] = string(round(k₀,1));
        attrs(g["potencial_degrau"])["sigma"] = string(round(σ,1));
    end

    # Caso não tenha salvo a solução para o t_f corrente, salve-a:
    try
        h[nome_dados_real] = real(ψ_1);
        h[nome_dados_imag] = imag(ψ_1);
    end
end





primeira_camada = h5open("dados_split_operator.h5","r") do file
    names(file)
end
#
#
#
segunda_camada = h5open("dados_split_operator.h5","r") do file
    segunda_camada = [];
    for name in primeira_camada
        append!(segunda_camada,names(file[name]))
    end
    return segunda_camada
end
#

terceira_camada = h5open("dados_split_operator.h5","r") do file
    terceira_camada = [];
    for name in primeira_camada
        g = file[name];
        for name_grupo in segunda_camada
            append!(terceira_camada,names(g[name_grupo]))
        end
    end
    return terceira_camada
end

#
Atr_grupo = h5readattr("dados_split_operator.h5", "split_1D/potencial_degrau")
#
#
#
# h5writeattr("dados_split_operator.h5", "split_1D/potencial_degrau",
#         Dict( "x_0" => string(round(x₀,1)),
#               "k_0" => string(round(k₀,1)),
#               "sigma" => string(round(σ,1))));
#
h5readattr("dados_split_operator.h5","split_1D/potencial_degrau/")
#
# name_path_real = string("split_1D/potencial_degrau/",nome_dados_real)
# name_path_imag = string("split_1D/potencial_degrau/",nome_dados_imag)
#
# psi_lido_real = h5read("dados_split_operator.h5",name_path_real)
# psi_lido_imag = h5read("dados_split_operator.h5",name_path_imag)
#
# psi_lido_total = complex(psi_lido_real,psi_lido_imag)
