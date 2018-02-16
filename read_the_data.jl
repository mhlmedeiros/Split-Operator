using HDF5;
using PyPlot;


#*******************************************************************************

primeira_camada = h5open("dados_split_operator.h5","r") do file
    names(file)
end

segunda_camada = h5open("dados_split_operator.h5","r") do file
    segunda_camada = [];
    for name in primeira_camada
        append!(segunda_camada,names(file[name]))
    end
    return segunda_camada
end

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

#*******************************************************************************

Atr_grupo = h5readattr("dados_split_operator.h5", "split_1D/potencial_degrau")

Atr_potencial = h5readattr("dados_split_operator.h5","split_1D/potencial_degrau/Potencial")

#*******************************************************************************

V_lido = h5open("dados_split_operator.h5","r") do file
    group = file[primeira_camada[1]];
    subgroup = group[segunda_camada[1]];
    read(subgroup["Potencial"])
end


x_lido = h5open("dados_split_operator.h5","r") do file
    group = file[primeira_camada[1]];
    subgroup = group[segunda_camada[1]];
    read(subgroup["grid"])
end

# plot(x_lido, V_lido)

#*******************************************************************************

terceira_camada[6]

nome_dados_real = terceira_camada[6]
nome_dados_imag = terceira_camada[5]

name_path_real = string("split_1D/potencial_degrau/",nome_dados_real)
name_path_imag = string("split_1D/potencial_degrau/",nome_dados_imag)

psi_lido_real = h5read("dados_split_operator.h5",name_path_real)
psi_lido_imag = h5read("dados_split_operator.h5",name_path_imag)

psi_lido_total = complex(psi_lido_real,psi_lido_imag)


PyPlot.matplotlib[:rc]("text", usetex=true) # allow tex rendering
PyPlot.matplotlib[:rc]("font", family="serif")

figure(2)
plot(x_lido, abs(psi_lido_total).^2)
plot(x_lido, 0.5*10.0^(-2)*V_lido, linewidth=2.0);
title("Gaussian distribution");
xlabel(L"$x$");
ylabel(L"$|\psi(x)|^2$");
# ylim(-0.05,1.05);
# xlim(-20,20);
