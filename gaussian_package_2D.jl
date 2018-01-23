#*******************************************************************************
#
#
#                Aplicação do split-operator em um sistema 2D
#               com condições de contorno tipo hardwall.
#
#*******************************************************************************
using PyPlot

# constantes físicas
Ry = 13.6056917253e3; # meV
a0 = 0.0529177208319; # nm
h2m = (Ry*a0^2); #  massa efetiva GaAs
hbar = 0.65821188926; # meV.ps


################################################################################
                    # Grid Bidimensional (Espaço real)
################################################################################
Lx = 1000; # Tamanho da dimensão "x" em nanometros (nm)
Ly = 400; # Tamanho da dimensão "y" em nanometros (nm)
Nx = div(Lx,2) # número de pontos na dimensão x
Ny = div(Ly,2) # número de pontos na dimensão y
# Nx = 2*Lx;
# Ny = 2*Ly;

x = linspace(-Lx/2, Lx/2, Nx+1)[1:end-1];
y = linspace(-Ly/2, Ly/2, Ny+1)[1:end-1];
dx = x[2] - x[1];
dy = y[2] - y[1];
X = kron(x', ones(Ny));
Y = kron(ones(Nx)', y);

################################################################################
                    # Espaço-k (na dimensãolongitudinal à tira)
################################################################################
kx = fftshift(linspace(-pi/dx, pi/dx, Nx+1)[1:end-1]); # Valores possíveis de k_x
k0x_max = pi/(2*dx);
k0y_max = pi/(2*dy);


################################################################################
                    # Condição Inicial: Pacote Gaussiano
################################################################################
x_0 = 0; #  centro do pacote
y_0 = 0; #  centro do pacote
sigma_x = 50.0; # largura longitudinal
sigma_y = 50.0; # largura transversal
k0_x = 0.1 * k0x_max; # momento inicial
k0_y = 0; # momento inicial
# Norm = ((pi^2)*(sigma_x*sigma_y)^2)^(-1/4); # Fator para normalização
Norm = 1 ; # Fator para pacote NÃO normalizado

# pacote de ondas na representação de posições:
Ψ_0 = Norm * exp(-(X-x_0).^2/(2*sigma_x^2)-(Y-y_0).^2/(2*sigma_y^2)
                    + 1im*(k0_x*X + k0_y*Y))
# fft!(Ψ_0,2)
#
# contourf(X,Y,fftshift(abs(Ψ_0).^2,2))
# ylim(-200,200)
################################################################################
                        # Hamiltoniano
################################################################################

# Potencial
α, β, γ = 5, 0.01, 0;
V_degrau = α./(e.^(-β*(X-γ))+1);

# pcolormesh(x,y,V_degrau)

# Cinética:
# matrizes Tridiagonais (Ny x Ny)
d2dy2 = 1/(dy^2) * SymTridiagonal(-2*ones(Ny), ones(Ny-1));
# ddy = 1/(2*dy) * Tridiagonal(-1 * ones(Ny-1), zeros(Ny), ones(Ny-1));

T_y  = - h2m * full(d2dy2)
T_kx =   h2m * kx


function psi_dt(ψ,Δt)
  # espaço de posições 'x'
  psi_spinor = ψ;
  psi_spinor = [exp(-1im*V_degrau[j,i]*Δt/(2*hbar))*psi_spinor[j,i] for j=1:Ny, i=1:Nx ];
  # representação mista (kx,y)
  psi_spinor = fft(psi_spinor,2);

  for i = 1:Nx
    psi_spinor[:,i] = expm(-1im*(T_kx[i] + T_y)*Δt/hbar)*psi_spinor[:,i] ;
  end

  # espaço de posições 'x'
  psi_spinor = ifft(psi_spinor,2);
  psi_spinor = [exp(-1im*V_degrau[j,i]*Δt/(2*hbar))*psi_spinor[j,i] for j=1:Ny, i=1:Nx ];

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

function psit_livre(Ψ,Δt)
  fft!(Ψ,2)
  psi = copy(Ψ);
  for i = 1:Nx
    psi[:,i] = expm(-1im*(T_kx[i] + T_y)*Δt/hbar)*Ψ[:,i] ;
  end
  ifft(psi,2)
end

# ψ_1 = psit(Ψ_0, 0.01,5)
ψ_1 = psit_livre(Ψ_0, 0.01)
pcolormesh(x,y,abs(ψ_1).^2)
ylim(-200,200)
# plot(x,10.0^(-1)*V_x, linewidth=2.0);
# # ylim(-0.05,1.05);
# # xlim(-20,20)
# title("Gaussian distribution");
# xlabel(L"$x$");
# ylabel(L"$|\psi(x)|^2$");
