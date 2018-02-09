#*******************************************************************************
#
#
#                Aplicação do split-operator em um sistema 2D
#               com condições de contorno tipo hardwall.
#
#
#*******************************************************************************

using PyPlot
using PyCall

@pyimport matplotlib.animation as anim # necessário para gerar vídeo da simulação

PyPlot.matplotlib[:rc]("text", usetex=true) # allow tex rendering
PyPlot.matplotlib[:rc]("font", family="serif")

# constantes físicas
Ry = 13.6056917253e3; # meV
a0 = 0.0529177208319; # nm
h2m = (Ry*a0^2); #  massa efetiva GaAs
hbar = 0.65821188926; # meV.ps


################################################################################
                    # Grid Bidimensional (Espaço real)
################################################################################
Lx = 1000; # Tamanho da dimensão "x" em nanometros (nm)
Ly = 200; # Tamanho da dimensão "y" em nanometros (nm)
Nx = div(Lx,3) # número de pontos na dimensão x
Ny = div(Ly,1) # número de pontos na dimensão y
# Nx = 2*Lx;
# Ny = 2*Ly;

x = linspace(-Lx/2, Lx/2, Nx+1)[1:end-1];
y = linspace(-Ly/2, Ly/2, Ny+1)[1:end-1];
dx = x[2] - x[1];
dy = y[2] - y[1];
X = kron(x', ones(Ny))
Y = kron(ones(Nx)', y)

################################################################################
                    # Espaço-k (na dimensão longitudinal à tira)
################################################################################
# A array "kx" deve ser 'invertida' via 'fftshift' para que todas as
# as funções de 'kx' possam ser multiplicadas elementwise sem confusão
# quando estivermos trabalhando com trasformadas de Fourier.

# Se houver interesse em observar como transformadas de Fourier e funções de
# 'kx' é interessante antes do 'plotting' reverter o 'shift'.

kx = fftshift(linspace(-pi/dx, pi/dx, Nx+1)[1:end]); # Valores possíveis de k_x
k0_x = pi/(2*dx);
k0_y = pi/(2*dy);


################################################################################
                    # Condição Inicial: Pacote Gaussiano
################################################################################

# Note que a condição inicial está representada no espaço de posição.

x_0 = -50; #  centro do pacote
y_0 = 0; #  centro do pacote
sigma_x = 30.0; # largura longitudinal
sigma_y = 30.0; # largura transversal
k0_x *= 0.5; # fração do kx_max
k0_y *= 0;   # fração do ky_max
# Norm = ((pi^2)*(sigma_x*sigma_y)^2)^(-1/4); # Fator para normalização
Norm = 1 ; # Fator para pacote NÃO normalizado

# pacote de ondas na representação de posições:
Ψ_0 = Norm * exp(-(X-x_0).^2/(2*sigma_x^2)-(Y-y_0).^2/(2*sigma_y^2)
                    + 1im*(k0_x*X + k0_y*Y))


################################################################################
                        # Hamiltoniano
################################################################################

#  O operador Hamiltoniano aqui é simplesmente aquele para uma partícula
# (Elétron) submetida a um potencial dependente da posição.

# Potencial:
α, β, γ = 10., .1, 50.;
V_degrau = α./(e.^(-β*(X-γ))+1);

# Cinética:
# matrizes Tridiagonais (Ny x Ny)
d2dy2 = 1/(dy^2) * SymTridiagonal(-2*ones(Ny), ones(Ny-1))

T_y  = - h2m * full(d2dy2)
T_kx =   h2m * kx.^2


################################################################################
                        # Funções para evolução temporal
################################################################################

function plot_simples(Ψ)
  pcolormesh(x,y,abs(Ψ).^2)
  ylim(-200,200)
  title("Gaussian distribution");
  xlabel(L"$x$");
  ylabel(L"$y$");
  return 0
end

function psi_dt(Ψ, Δt)
  # espaço de posições 'x'
  Ψ = [ exp(-1im*V_degrau[j,i]*Δt/(2*hbar)) * Ψ[j,i] for j=1:Ny, i=1:Nx ];

  # representação mista (kx,y)
  fft!(Ψ,2);
  for i = 1:Nx
    Ψ[:,i] = expm(-1im*(T_kx[i]*eye(Ny) + T_y)*Δt/hbar) * Ψ[:,i] ;
  end

  # espaço de posições 'x'
  ifft!(Ψ,2);
  Ψ = [exp(-1im*V_degrau[j,i]*Δt/(2*hbar)) * Ψ[j,i] for j=1:Ny, i=1:Nx ];

  return Ψ;
end

function psi_dt_animado(Δt)
  global Ψ_0
  if Δt == 0
    # Essa  condição existe para incluir a
    # cond. inicial no video da simulação
    clf()
    plot_simples(Ψ_0)
    return 0
  else
    # espaço de posições 'x'
    Ψ_0 = [ exp(-1im*V_degrau[j,i]*Δt/(2*hbar)) * Ψ_0[j,i] for j=1:Ny, i=1:Nx ];

    # representação mista (kx,y)
    fft!(Ψ_0,2);
    for i = 1:Nx
      Ψ_0[:,i] = expm(-1im*(T_kx[i]*eye(Ny) + T_y)*Δt/hbar) * Ψ_0[:,i] ;
    end

    # espaço de posições 'x'
    ifft!(Ψ_0,2);
    Ψ_0 = [exp(-1im*V_degrau[j,i]*Δt/(2*hbar)) * Ψ_0[j,i] for j=1:Ny, i=1:Nx ];

    clf()
    plot_simples(Ψ_0)
    return 0;
  end
end

function psit(Ψ, t_final, N_steps = 20)
  tspan = linspace(0, t_final, N_steps-1);
  delta_t = tspan[2]-tspan[1];

  for i = 1:N_steps
    psi_aux = psi_dt(Ψ,delta_t);
    Ψ = psi_aux;
  end
  Ψ
end

function psit_livre(Ψ, t)
  fft!(Ψ,2); # substitui Ψ por sua transformada de Fourier em "x"
  for i = 1:Nx
    Ψ[:,i] = expm(-1im*(T_kx[i]*eye(Ny) + T_y)* t/hbar)*Ψ[:,i] ;
  end
  ifft!(Ψ,2)
end

function simula_potencial(;t_final = 2.0, N_steps = 10, intervalo_frames = 100)

  Δt = t_final/N_steps;
  array_of_deltas = Δt*[0; ones(N_steps)];

  fig = figure();

  video = anim.FuncAnimation(fig,
                            psi_dt_animado,
                            frames = array_of_deltas,
                            interval = intervalo_frames);

  path = "/home/marcos/Dropbox/projetos/julia/split operator/";
  nome = string("simulação_com_potencial_degrau_em_x_50.mp4")

  video[:save](string(path,nome), extra_args=["-vcodec", "libx264"]);
  return 0
end


################################################################################
                        # Apresentação dos resultados
################################################################################

Ψ_1 = psit(Ψ_0, 2.0, 6)
# Ψ_1 = psit_livre(Ψ_0, 2)

figure(1)
pcolormesh(x,y,abs(Ψ_1).^2)
ylim(-200,200)
title("Gaussian distribution");
xlabel(L"$x$");
ylabel(L"$y$");


figure(2)
pcolormesh(X,Y,V_degrau)
ylim(-200,200)
title("Potential function");
xlabel(L"$x$");
ylabel(L"$y$");
colorbar()

# simula_potencial(t_final=10.0,N_steps=100,intervalo_frames=50)
